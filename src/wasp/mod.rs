//! WASP allele-specific-mapping filter.
//!
//! Ported from STAR-rs `crates/star-align/src/star_wasp.rs`, which itself ports
//! STAR's `ReadAlign_waspMap.cpp`, `Transcript::variationAdjust` and
//! `Variation::loadVCF`. See `docs-old/dev/porting-from-star-rs.md`.
//!
//! A read overlapping heterozygous SNP(s) is re-mapped with every alternative
//! allele combination; if all re-maps land on the identical locus the read passes
//! (`vW:i:1`), otherwise a failure code is reported. Reads overlapping no variant
//! get no `vW` tag (code `-1`).
//!
//! Covers **single-end** ([`annotate_records_se`]) and **paired-end**
//! ([`annotate_records_pe`]). For a pair, variants overlapping mate1 and mate2
//! are walked in mate1-then-mate2 order (STAR's combined `mate1 ++ rc(mate2)`
//! transcript frame), a single `vW` code is computed by re-mapping the whole
//! pair with each allele combination, and both mate records share `vW`/`vA`/`vG`.
//!
//! Everything works in rustar-aligner's base-code space (`A=0,C=1,G=2,T=3,N=4`),
//! and variant/read overlap is computed by walking the transcript CIGAR in
//! SAM/forward-genomic order (the same convention as `io::sam::build_md_tag`),
//! so it does not depend on the exon read-coordinate frame.

use std::path::Path;

use noodles::sam::alignment::record::cigar;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;

use crate::align::read_align::{
    PairedAlignment, PairedAlignmentResult, align_paired_read, align_read,
};
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::io::fastq::complement_base;
use crate::params::{Parameters, SamAttributes};

/// One heterozygous SNV: absolute 0-based genomic `loci`, and `nt = [ref, allele0,
/// allele1]` as base codes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Snp {
    pub loci: u64,
    pub nt: [u8; 3],
}

fn nt_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Reverse-complement a base-code sequence (`N`=4 preserved).
fn rc_codes(v: &[u8]) -> Vec<u8> {
    v.iter().rev().map(|&b| complement_base(b)).collect()
}

/// STAR `scanVCF`: single-char REF and ALT, genotype not homozygous-reference and
/// heterozygous (`gt0 != gt1`); `nt = [REF, altV[gt0], altV[gt1]]` with
/// `altV = [REF, ALT...]`. Returns SNPs sorted by absolute genomic position (only
/// ACGT variants on known chromosomes). `chr_names`/`chr_starts` come from the
/// loaded genome.
pub fn load_vcf(
    path: &Path,
    chr_names: &[String],
    chr_starts: &[u64],
) -> std::io::Result<Vec<Snp>> {
    let text = std::fs::read_to_string(path)?;
    let mut snps = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let (chr, pos, r#ref, alt, sample) = (f[0], f[1], f[3], f[4], f[9]);
        if r#ref.len() != 1 {
            continue;
        }
        let alt_v: Vec<&str> = std::iter::once(r#ref).chain(alt.split(',')).collect();
        if alt_v.iter().skip(1).any(|a| a.len() != 1) {
            continue;
        }
        let sb = sample.as_bytes();
        if sb.len() < 3 {
            continue;
        }
        let (Some(gt0), Some(gt1)) = ((sb[0] as char).to_digit(10), (sb[2] as char).to_digit(10))
        else {
            continue;
        };
        // Skip homozygous-reference / homozygous genotypes (heteroOnly default).
        if (gt0 == 0 && gt1 == 0) || gt0 == gt1 {
            continue;
        }
        let Some(chr_idx) = chr_names.iter().position(|c| c == chr) else {
            continue;
        };
        let (Some(a0), Some(a1)) = (
            alt_v.get(gt0 as usize).map(|s| s.as_bytes()[0]),
            alt_v.get(gt1 as usize).map(|s| s.as_bytes()[0]),
        ) else {
            continue;
        };
        let nt = [nt_code(r#ref.as_bytes()[0]), nt_code(a0), nt_code(a1)];
        if nt.iter().any(|&c| c >= 4) {
            continue;
        }
        let Ok(pos1) = pos.parse::<u64>() else {
            continue;
        };
        snps.push(Snp {
            loci: pos1 - 1 + chr_starts[chr_idx],
            nt,
        });
    }
    snps.sort_by_key(|s| s.loci);
    Ok(snps)
}

/// The matched-allele code for a read base at a SNP: `1`/`2` for the SNP's two
/// alleles, `3` if neither, `4` if the read base is `N`.
fn classify_allele(snp: &Snp, read_code: u8) -> u8 {
    if read_code > 3 {
        4
    } else if snp.nt[1] == read_code {
        1
    } else if snp.nt[2] == read_code {
        2
    } else {
        3
    }
}

/// Variants overlapping an alignment, by walking the CIGAR in SAM/forward-genomic
/// order. `sam_codes` is the read in that frame (reverse-complemented for a reverse
/// alignment). Returns `(snp_index, read_pos_in_sam_frame, allele)` per overlapping
/// SNP, in genomic order. `snps` must be sorted by `loci`.
fn variation_overlap(snps: &[Snp], sam_codes: &[u8], tr: &Transcript) -> Vec<(usize, usize, u8)> {
    use cigar::op::Kind;
    let mut out = Vec::new();
    let mut genome_pos = tr.genome_start;
    let mut read_pos: usize = 0;
    for op in &tr.cigar {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    if let Ok(isnp) = snps.binary_search_by(|s| s.loci.cmp(&genome_pos)) {
                        let read_code = sam_codes.get(read_pos).copied().unwrap_or(4);
                        out.push((isnp, read_pos, classify_allele(&snps[isnp], read_code)));
                    }
                    genome_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Deletion | Kind::Skip => {
                genome_pos += op.len() as u64;
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += op.len();
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    out
}

/// STAR `Transcript::variationAdjust`: the `vG`/`vA` tag values for an alignment
/// overlapping heterozygous SNP(s). `vg` = each SNP's chromosome-relative 0-based
/// genomic coordinate; `va` = the matched-allele code, in genomic order. Empty when
/// the alignment overlaps no variant. `read_codes` is the forward read (base codes);
/// `chr_start` is the alignment chromosome's genome offset.
pub fn wasp_variants(
    chr_start: u64,
    snps: &[Snp],
    read_codes: &[u8],
    tr: &Transcript,
) -> (Vec<i32>, Vec<u8>) {
    let sam_codes = if tr.is_reverse {
        rc_codes(read_codes)
    } else {
        read_codes.to_vec()
    };
    let mut vg = Vec::new();
    let mut va = Vec::new();
    for (isnp, _rp, allele) in variation_overlap(snps, &sam_codes, tr) {
        vg.push((snps[isnp].loci - chr_start) as i32);
        va.push(allele);
    }
    (vg, va)
}

fn same_exons(a: &Transcript, b: &Transcript) -> bool {
    a.exons.len() == b.exons.len()
        && a.exons.iter().zip(b.exons.iter()).all(|(x, y)| {
            x.read_start == y.read_start
                && x.genome_start == y.genome_start
                && x.genome_end == y.genome_end
        })
}

/// STAR `waspMap` (single-end): the `vW` code for one alignment, or `-1` for no `vW`
/// tag (read overlaps no variant). `1` = passed; `2` multimaps; `3` variant base is
/// `N`; `4` a remap is unmapped; `5` a remap multimaps; `6` a remap lands elsewhere;
/// `7` too many variants. `n_tr` is the read's mapped-locus count; `remap_params`
/// is the WASP-relaxed [`wasp_remap_params`].
pub fn wasp_type(
    index: &GenomeIndex,
    snps: &[Snp],
    read_codes: &[u8],
    read_name: &str,
    tr: &Transcript,
    n_tr: usize,
    remap_params: &Parameters,
) -> Result<i32, Error> {
    let sam_codes = if tr.is_reverse {
        rc_codes(read_codes)
    } else {
        read_codes.to_vec()
    };
    let vars = variation_overlap(snps, &sam_codes, tr);
    if vars.is_empty() {
        return Ok(-1);
    }
    if n_tr > 1 {
        return Ok(2);
    }
    if vars.len() > 10 {
        return Ok(7);
    }
    if vars.iter().any(|&(_, _, allele)| allele > 3) {
        return Ok(3);
    }

    let actual: Vec<u8> = vars.iter().map(|&(_, _, a)| a).collect();
    let n = vars.len();
    // All 2^n allele combinations of {1, 2}.
    for mask in 0..(1u32 << n) {
        let combo: Vec<u8> = (0..n)
            .map(|i| if mask & (1 << i) != 0 { 2 } else { 1 })
            .collect();
        if combo == actual {
            continue;
        }
        // Flip each variant's base to the combination's allele in the SAM frame,
        // then map back to the forward read and re-align.
        let mut modified = sam_codes.clone();
        for (iv, &(isnp, read_pos, _)) in vars.iter().enumerate() {
            modified[read_pos] = snps[isnp].nt[combo[iv] as usize];
        }
        let remap_read = if tr.is_reverse {
            rc_codes(&modified)
        } else {
            modified
        };
        let (trs, _chim, _n_for_mapq, _reason) =
            align_read(&remap_read, read_name, index, remap_params)?;
        if trs.is_empty() {
            return Ok(4); // remap unmapped or too-many-loci
        }
        if trs.len() > 1 {
            return Ok(5);
        }
        if !same_exons(&trs[0], tr) {
            return Ok(6);
        }
    }
    Ok(1)
}

/// The parameter set used for WASP re-mapping: relaxed match filter and a fixed
/// multimap window, matching STAR's `waspMap` (STAR-rs passes `match_nmin=0`,
/// `multimap_score_range=1`, `multimap_nmax=10`). Built once and reused per read.
pub fn wasp_remap_params(params: &Parameters) -> Parameters {
    // STAR re-maps the allele-swapped read with the UNMODIFIED user parameters:
    // `waspRA = new ReadAlign(Pin,…)` shares `P`, then runs the normal
    // mapOneRead() → multMapSelect() → mappedFilter() (ReadAlign_waspMap.cpp:78-92).
    // A read that fails mappedFilter under the user's own thresholds becomes
    // waspType=4 (fails WASP). Relaxing outFilterMatchNminOverLread to 0 or
    // hardcoding outFilterMultimapNmax=10 (as the original port did) is a STAR-rs
    // divergence that lets swapped reads pass WASP that STAR would reject — so we
    // re-map with the caller's params verbatim.
    params.clone()
}

/// Insert `vW`/`vA`/`vG` tags on a record, gated on the requested attributes.
fn insert_wasp_tags(record: &mut RecordBuf, vw: i32, vg: &[i32], va: &[u8], attrs: SamAttributes) {
    let data = record.data_mut();
    if attrs.contains(SamAttributes::VW) {
        data.insert(Tag::new(b'v', b'W'), Value::from(vw));
    }
    if attrs.contains(SamAttributes::VA) && !va.is_empty() {
        data.insert(
            Tag::new(b'v', b'A'),
            Value::Array(Array::Int8(va.iter().map(|&x| x as i8).collect())),
        );
    }
    if attrs.contains(SamAttributes::VG) && !vg.is_empty() {
        data.insert(
            Tag::new(b'v', b'G'),
            Value::Array(Array::Int32(vg.to_vec())),
        );
    }
}

/// Loaded WASP state for a run: the heterozygous SNVs and the relaxed re-map
/// parameter set. Built once and shared read-only across the per-read parallel loop.
pub struct WaspContext {
    pub snps: Vec<Snp>,
    pub remap_params: Parameters,
}

impl WaspContext {
    /// Load the VCF and build the WASP re-map parameters for a run.
    pub fn load(
        vcf_path: &Path,
        chr_names: &[String],
        chr_starts: &[u64],
        params: &Parameters,
    ) -> std::io::Result<Self> {
        Ok(Self {
            snps: load_vcf(vcf_path, chr_names, chr_starts)?,
            remap_params: wasp_remap_params(params),
        })
    }
}

/// Compute WASP tags for a single-end read and stamp them onto its already-built
/// SAM records. `records[i]` corresponds to `transcripts[i]`. A uniquely-mapped read
/// is re-mapped (computing `vW`); a multi-mapped read overlapping variants gets
/// `vW:i:2`. Reads overlapping no variant are left untagged.
pub fn annotate_records_se(
    records: &mut [RecordBuf],
    transcripts: &[Transcript],
    read_codes: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    ctx: &WaspContext,
    attrs: SamAttributes,
) -> Result<(), Error> {
    let snps = &ctx.snps;
    if transcripts.len() == 1 {
        let tr = &transcripts[0];
        let vw = wasp_type(index, snps, read_codes, read_name, tr, 1, &ctx.remap_params)?;
        if vw != -1 {
            let chr_start = index.genome.chr_start[tr.chr_idx];
            let (vg, va) = wasp_variants(chr_start, snps, read_codes, tr);
            if let Some(rec) = records.get_mut(0) {
                insert_wasp_tags(rec, vw, &vg, &va, attrs);
            }
        }
    } else {
        // STAR keys waspType off the PRIMARY alignment (trBest = transcripts[0],
        // ReadAlign_waspMap.cpp): a multimapper whose primary overlaps ≥1 variant
        // gets waspType=2 stamped on ALL of its records; if the primary overlaps no
        // variant, no record gets vW (even if a secondary locus does). vG/vA stay
        // per-record (STAR runs variationAdjust per output alignment).
        let primary = &transcripts[0];
        let p_chr_start = index.genome.chr_start[primary.chr_idx];
        let (p_vg, _) = wasp_variants(p_chr_start, snps, read_codes, primary);
        if !p_vg.is_empty() {
            for (rec, tr) in records.iter_mut().zip(transcripts.iter()) {
                let chr_start = index.genome.chr_start[tr.chr_idx];
                let (vg, va) = wasp_variants(chr_start, snps, read_codes, tr);
                insert_wasp_tags(rec, 2, &vg, &va, attrs);
            }
        }
    }
    Ok(())
}

/// STAR `Transcript::variationAdjust` for a paired-end alignment: the combined
/// `vG`/`vA` tag values from variants overlapping mate1 then mate2 (STAR's
/// `mate1 ++ rc(mate2)` combined-transcript order). Both mate records share this
/// array. `chr_start1`/`chr_start2` are the mates' chromosome genome offsets;
/// `m1_codes`/`m2_codes` are the forward mate reads (base codes).
#[allow(clippy::too_many_arguments)]
pub fn wasp_variants_pe(
    chr_start1: u64,
    chr_start2: u64,
    snps: &[Snp],
    m1_codes: &[u8],
    m2_codes: &[u8],
    tr1: &Transcript,
    tr2: &Transcript,
) -> (Vec<i32>, Vec<u8>) {
    let mut vg = Vec::new();
    let mut va = Vec::new();
    for (chr_start, codes, tr) in [(chr_start1, m1_codes, tr1), (chr_start2, m2_codes, tr2)] {
        let sam_codes = if tr.is_reverse {
            rc_codes(codes)
        } else {
            codes.to_vec()
        };
        for (isnp, _rp, allele) in variation_overlap(snps, &sam_codes, tr) {
            vg.push((snps[isnp].loci - chr_start) as i32);
            va.push(allele);
        }
    }
    (vg, va)
}

/// STAR `waspMap` (paired-end): the `vW` code for a pair, or `-1` for no `vW` tag
/// (neither mate overlaps a variant). Codes match [`wasp_type`]: `1` passed, `2`
/// pair multimaps, `3` a variant base is `N`, `4` a re-map is not both-mapped,
/// `5` a re-map multimaps, `6` a re-map lands elsewhere, `7` too many variants.
/// Each allele combination is applied per mate in its SAM frame, mapped back to
/// forward mates, and re-mapped as a pair with the WASP-relaxed [`remap_params`].
#[allow(clippy::too_many_arguments)]
pub fn wasp_type_pe(
    index: &GenomeIndex,
    snps: &[Snp],
    m1_codes: &[u8],
    m2_codes: &[u8],
    read_name: &str,
    tr1: &Transcript,
    tr2: &Transcript,
    n_tr: usize,
    remap_params: &Parameters,
) -> Result<i32, Error> {
    let sam1 = if tr1.is_reverse {
        rc_codes(m1_codes)
    } else {
        m1_codes.to_vec()
    };
    let sam2 = if tr2.is_reverse {
        rc_codes(m2_codes)
    } else {
        m2_codes.to_vec()
    };
    let vars1 = variation_overlap(snps, &sam1, tr1);
    let vars2 = variation_overlap(snps, &sam2, tr2);
    if vars1.is_empty() && vars2.is_empty() {
        return Ok(-1);
    }
    if n_tr > 1 {
        return Ok(2);
    }
    let n1 = vars1.len();
    let n = n1 + vars2.len();
    if n > 10 {
        return Ok(7);
    }
    if vars1
        .iter()
        .chain(vars2.iter())
        .any(|&(_, _, allele)| allele > 3)
    {
        return Ok(3);
    }

    // Combined actual alleles, mate1 then mate2 (matches STAR's combined read).
    let actual: Vec<u8> = vars1
        .iter()
        .chain(vars2.iter())
        .map(|&(_, _, a)| a)
        .collect();
    for mask in 0..(1u32 << n) {
        let combo: Vec<u8> = (0..n)
            .map(|i| if mask & (1 << i) != 0 { 2 } else { 1 })
            .collect();
        if combo == actual {
            continue;
        }
        // Flip each variant's base to the combination's allele in each mate's SAM
        // frame, then map back to the forward mate reads.
        let mut mod1 = sam1.clone();
        for (iv, &(isnp, read_pos, _)) in vars1.iter().enumerate() {
            mod1[read_pos] = snps[isnp].nt[combo[iv] as usize];
        }
        let mut mod2 = sam2.clone();
        for (iv, &(isnp, read_pos, _)) in vars2.iter().enumerate() {
            mod2[read_pos] = snps[isnp].nt[combo[n1 + iv] as usize];
        }
        let remap1 = if tr1.is_reverse {
            rc_codes(&mod1)
        } else {
            mod1
        };
        let remap2 = if tr2.is_reverse {
            rc_codes(&mod2)
        } else {
            mod2
        };

        let (results, _chim, _n_for_mapq, _reason) =
            align_paired_read(&remap1, &remap2, read_name, index, remap_params)?;
        let both: Vec<&PairedAlignment> = results
            .iter()
            .filter_map(|r| match r {
                PairedAlignmentResult::BothMapped(pa) => Some(pa.as_ref()),
                PairedAlignmentResult::HalfMapped { .. } => None,
            })
            .collect();
        if both.is_empty() {
            return Ok(4); // re-map unmapped / half-mapped / too-many-loci
        }
        if both.len() > 1 {
            return Ok(5);
        }
        let pa = both[0];
        if !same_exons(&pa.mate1_transcript, tr1) || !same_exons(&pa.mate2_transcript, tr2) {
            return Ok(6);
        }
    }
    Ok(1)
}

/// Compute WASP tags for a paired-end read and stamp them onto its already-built
/// SAM records. `records` is `[mate1, mate2, mate1, mate2, ...]` (two per pair, as
/// [`crate::io::sam::SamWriter::build_paired_records`] emits); `paired[i]`
/// corresponds to `records[2i]`/`records[2i+1]`. A uniquely-mapped pair is
/// re-mapped (computing a single `vW` shared by both mate lines); a multi-mapped
/// pair overlapping variants gets `vW:i:2`. Pairs overlapping no variant are left
/// untagged.
#[allow(clippy::too_many_arguments)]
pub fn annotate_records_pe(
    records: &mut [RecordBuf],
    paired: &[PairedAlignment],
    m1_codes: &[u8],
    m2_codes: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    ctx: &WaspContext,
    attrs: SamAttributes,
) -> Result<(), Error> {
    let snps = &ctx.snps;
    if paired.len() == 1 {
        let pa = &paired[0];
        let (tr1, tr2) = (&pa.mate1_transcript, &pa.mate2_transcript);
        let vw = wasp_type_pe(
            index,
            snps,
            m1_codes,
            m2_codes,
            read_name,
            tr1,
            tr2,
            1,
            &ctx.remap_params,
        )?;
        if vw != -1 {
            let cs1 = index.genome.chr_start[tr1.chr_idx];
            let cs2 = index.genome.chr_start[tr2.chr_idx];
            let (vg, va) = wasp_variants_pe(cs1, cs2, snps, m1_codes, m2_codes, tr1, tr2);
            for rec in records.iter_mut().take(2) {
                insert_wasp_tags(rec, vw, &vg, &va, attrs);
            }
        }
    } else {
        // Key waspType off the PRIMARY pair (paired[0], = STAR trBest): a multimapping
        // pair whose primary overlaps ≥1 variant gets vW:2 on ALL records; none if the
        // primary overlaps no variant. vG/vA stay per-pair.
        let primary = &paired[0];
        let (p1, p2) = (&primary.mate1_transcript, &primary.mate2_transcript);
        let (pcs1, pcs2) = (
            index.genome.chr_start[p1.chr_idx],
            index.genome.chr_start[p2.chr_idx],
        );
        let (p_vg, _) = wasp_variants_pe(pcs1, pcs2, snps, m1_codes, m2_codes, p1, p2);
        if !p_vg.is_empty() {
            for (i, pa) in paired.iter().enumerate() {
                let (tr1, tr2) = (&pa.mate1_transcript, &pa.mate2_transcript);
                let cs1 = index.genome.chr_start[tr1.chr_idx];
                let cs2 = index.genome.chr_start[tr2.chr_idx];
                let (vg, va) = wasp_variants_pe(cs1, cs2, snps, m1_codes, m2_codes, tr1, tr2);
                for k in 0..2 {
                    if let Some(rec) = records.get_mut(2 * i + k) {
                        insert_wasp_tags(rec, 2, &vg, &va, attrs);
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{Exon, Transcript};
    use noodles::sam::alignment::record::cigar::{Op, op::Kind};

    fn tr_fwd_50m(genome_start: u64) -> Transcript {
        Transcript {
            chr_idx: 0,
            genome_start,
            genome_end: genome_start + 50,
            is_reverse: false,
            exons: vec![Exon {
                genome_start,
                genome_end: genome_start + 50,
                read_start: 0,
                read_end: 50,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 49,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        }
    }

    #[test]
    fn nt_code_maps_bases() {
        assert_eq!(nt_code(b'A'), 0);
        assert_eq!(nt_code(b'c'), 1);
        assert_eq!(nt_code(b'G'), 2);
        assert_eq!(nt_code(b'T'), 3);
        assert_eq!(nt_code(b'N'), 4);
    }

    #[test]
    fn rc_codes_reverses_and_complements() {
        // A C G T N (0 1 2 3 4) -> reverse-complement -> N A C G T (4 0 1 2 3)
        assert_eq!(rc_codes(&[0, 1, 2, 3, 4]), vec![4, 0, 1, 2, 3]);
    }

    #[test]
    fn classify_allele_codes() {
        let snp = Snp {
            loci: 10,
            nt: [0, 0, 2],
        }; // ref A, alleles A/G
        assert_eq!(classify_allele(&snp, 0), 1); // A -> allele 1
        assert_eq!(classify_allele(&snp, 2), 2); // G -> allele 2
        assert_eq!(classify_allele(&snp, 1), 3); // C -> neither
        assert_eq!(classify_allele(&snp, 4), 4); // N
    }

    #[test]
    fn load_vcf_parses_het_snp() {
        use std::io::Write as _;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
        )
        .unwrap();
        // het 0/1 at chrI:100 (A->G); hom-ref skipped; multi-char REF skipped
        writeln!(f, "chrI\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(f, "chrI\t200\t.\tA\tG\t.\t.\t.\tGT\t0|0").unwrap();
        writeln!(f, "chrI\t300\t.\tAC\tG\t.\t.\t.\tGT\t0|1").unwrap();
        let names = vec!["chrI".to_string()];
        let starts = vec![0u64];
        let snps = load_vcf(f.path(), &names, &starts).unwrap();
        assert_eq!(snps.len(), 1);
        assert_eq!(snps[0].loci, 99); // 0-based, chr_start 0
        assert_eq!(snps[0].nt, [0, 0, 2]); // ref A, allele0 A, allele1 G
    }

    #[test]
    fn variation_overlap_forward() {
        // SNP at genomic 110, transcript maps read[0..50] to genome[100..150] fwd.
        let snps = vec![Snp {
            loci: 110,
            nt: [0, 0, 2],
        }];
        let tr = tr_fwd_50m(100);
        let mut read = vec![0u8; 50];
        read[10] = 2; // read base at the SNP is G (allele 2)
        let out = variation_overlap(&snps, &read, &tr);
        assert_eq!(out, vec![(0, 10, 2)]);
    }

    #[test]
    fn wasp_variants_reverse_strand_frame() {
        // Reverse alignment: sam frame is rc(read). SNP at genome 110.
        let snps = vec![Snp {
            loci: 110,
            nt: [0, 0, 2],
        }];
        let mut tr = tr_fwd_50m(100);
        tr.is_reverse = true;
        // sam_codes[10] must be the read base at the SNP; build read so rc(read)[10]=G(2).
        let mut sam = vec![0u8; 50];
        sam[10] = 2;
        let read = rc_codes(&sam); // forward read whose rc is `sam`
        let (vg, va) = wasp_variants(0, &snps, &read, &tr);
        assert_eq!(vg, vec![110]); // chr-relative coord
        assert_eq!(va, vec![2]);
    }

    #[test]
    fn wasp_variants_empty_when_no_overlap() {
        let snps = vec![Snp {
            loci: 500,
            nt: [0, 0, 2],
        }];
        let tr = tr_fwd_50m(100);
        let (vg, va) = wasp_variants(0, &snps, &[0u8; 50], &tr);
        assert!(vg.is_empty() && va.is_empty());
    }

    #[test]
    fn wasp_variants_pe_combines_mate1_then_mate2() {
        // SNP over mate1 (genome 110) and mate2 (genome 310); both forward.
        let snps = vec![
            Snp {
                loci: 110,
                nt: [0, 0, 2],
            },
            Snp {
                loci: 310,
                nt: [0, 0, 2],
            },
        ];
        let tr1 = tr_fwd_50m(100);
        let tr2 = tr_fwd_50m(300);
        let mut m1 = vec![0u8; 50];
        m1[10] = 2; // allele 2 at mate1 SNP
        let mut m2 = vec![0u8; 50];
        m2[10] = 2; // allele 2 at mate2 SNP
        let (vg, va) = wasp_variants_pe(0, 0, &snps, &m1, &m2, &tr1, &tr2);
        // mate1 variant first, then mate2 (STAR's mate1 ++ rc(mate2) order).
        assert_eq!(vg, vec![110, 310]);
        assert_eq!(va, vec![2, 2]);
    }

    #[test]
    fn wasp_variants_pe_reverse_mate2_frame() {
        // mate2 aligns reverse: its SAM frame is rc(mate2 read).
        let snps = vec![Snp {
            loci: 310,
            nt: [0, 0, 2],
        }];
        let tr1 = tr_fwd_50m(100);
        let mut tr2 = tr_fwd_50m(300);
        tr2.is_reverse = true;
        let m1 = vec![0u8; 50]; // no mate1 overlap
        // Build mate2 forward read whose rc has allele 2 at the SNP offset (10).
        let mut sam2 = vec![0u8; 50];
        sam2[10] = 2;
        let m2 = rc_codes(&sam2);
        let (vg, va) = wasp_variants_pe(0, 0, &snps, &m1, &m2, &tr1, &tr2);
        assert_eq!(vg, vec![310]);
        assert_eq!(va, vec![2]);
    }
}
