//! STAR-faithful splice junction insertion into the genome index.
//!
//! Ports `source/sjdbPrepare.cpp` + `source/sjdbBuildIndex.cpp` from STAR.
//! At `genomeGenerate` time, STAR extracts the flanking `sjdbOverhang`
//! bases on each side of every GTF-derived splice junction, concatenates
//! them into a `Gsj` buffer, appends that buffer to the `Genome` binary,
//! and extends the suffix array to index the new bases. This module
//! provides the same machinery for rustar-aligner so that the generated
//! `Genome` / `SA` / `SAindex` / `sjdbInfo.txt` / `sjdbList.out.tab` files
//! match STAR's byte-for-byte.
//!
//! The module is orchestrated from `index::GenomeIndex::build` after the
//! base suffix array has been built.

use std::fs::File;
use std::io::{BufWriter, Write as _};
use std::path::Path;

use crate::align::score::detect_splice_motif;
use crate::error::Error;
use crate::genome::Genome;
use crate::junction::encode_motif;

/// STAR's inter-SJ spacer byte in the Gsj buffer (same value STAR uses
/// for inter-chromosome padding — `IncludeDefine.h::GENOME_spacingChar`).
const GSJ_SPACING: u8 = 5;

/// Compute STAR's `(sjdbShiftLeft, sjdbShiftRight)` for an intron whose
/// 0-based donor/acceptor positions are `s` and `e`.
///
/// STAR defines these as the number of bases the intron boundary can shift
/// left / right while preserving the donor/acceptor base identity
/// (`sjdbPrepare.cpp:52-73`) — i.e. the repeat length across the junction.
/// Intended to land the junction at its left-most canonical position so
/// identical splice events produce identical SA indices regardless of
/// which exon pair they came from.
///
/// Stops at genome bounds, on any N-base (code ≥ 4), or at the 255 cap.
pub fn compute_shifts(genome: &Genome, s: u64, e: u64, n_genome_real: u64) -> (u8, u8) {
    let forward = &genome.sequence[..n_genome_real as usize];
    let si = s as usize;
    let ei = e as usize;

    let mut jj_l: u8 = 0;
    while jj_l < 255 && (jj_l as usize) < si && ei >= jj_l as usize {
        let a = forward[si - 1 - jj_l as usize];
        let b = forward[ei - jj_l as usize];
        if a != b || a >= 4 {
            break;
        }
        jj_l += 1;
    }

    let mut jj_r: u8 = 0;
    while jj_r < 255 && (ei + 1 + jj_r as usize) < forward.len() {
        let a = forward[si + jj_r as usize];
        let b = forward[ei + 1 + jj_r as usize];
        if a != b || a >= 4 {
            break;
        }
        jj_r += 1;
    }

    (jj_l, jj_r)
}

/// A single junction with all metadata STAR needs to write the sjdb
/// files and extend the genome. Positions are 0-based absolute genome
/// coordinates; `start_pos` / `end_pos` are the FIRST and LAST bases of
/// the intron (inclusive), already shifted left by `shift_left` so they
/// sit at the canonical (left-most) representation of the motif.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedJunction {
    /// Chromosome index the junction belongs to.
    pub chr_idx: usize,
    /// Shift-adjusted 0-based genome position of the first intron base.
    pub start_pos: u64,
    /// Shift-adjusted 0-based genome position of the last intron base.
    pub end_pos: u64,
    /// STAR motif code (0 = non-canonical, 1-6 = canonical variants).
    pub motif: u8,
    /// Repeat length to the left of the (pre-shift) donor.
    pub shift_left: u8,
    /// Repeat length to the right of the (pre-shift) acceptor.
    pub shift_right: u8,
    /// STAR strand code (0 = unknown/dot, 1 = +, 2 = -).
    pub strand: u8,
}

impl PreparedJunction {
    /// Value STAR writes into `mapGen.sjdbStart` after its post-dedup
    /// sort: ORIGINAL pre-shift start for canonical motifs, shifted
    /// start for non-canonical (matches `sjdbPrepare.cpp:127,174`).
    pub fn stored_start(&self) -> u64 {
        if self.motif == 0 {
            self.start_pos
        } else {
            self.start_pos + self.shift_left as u64
        }
    }

    /// Companion to `stored_start` — the `mapGen.sjdbEnd` value.
    pub fn stored_end(&self) -> u64 {
        if self.motif == 0 {
            self.end_pos
        } else {
            self.end_pos + self.shift_left as u64
        }
    }

    /// Original (pre-shift) 0-based position of the first intron base.
    /// Used for Gsj flanking-sequence extraction (STAR uses this regardless
    /// of motif).
    pub fn original_start(&self) -> u64 {
        self.start_pos + self.shift_left as u64
    }

    /// Companion to `original_start`.
    pub fn original_end(&self) -> u64 {
        self.end_pos + self.shift_left as u64
    }
}

/// Convert a splice-junction database entry to a fully-prepared entry
/// carrying motif, shifts, strand, and shift-adjusted coordinates.
///
/// `db_strand` is STAR's 0/1/2 (unknown/+/-); when it's 0 and the motif
/// is canonical, STAR derives the strand from the motif via
/// `2 - motif % 2` (see `sjdbPrepare.cpp:184-188`).
pub fn prepare_junction(
    chr_idx: usize,
    intron_start: u64,
    intron_end: u64,
    db_strand: u8,
    genome: &Genome,
    n_genome_real: u64,
) -> PreparedJunction {
    let intron_len = (intron_end - intron_start + 1) as u32;
    let motif = encode_motif(detect_splice_motif(intron_start, intron_len, genome));
    let (shift_left, shift_right) = compute_shifts(genome, intron_start, intron_end, n_genome_real);

    // sjdbPrepare.cpp:71-72 — land the junction at its left-most canonical
    // representation so identical splice events produce identical indices.
    let shifted_start = intron_start - shift_left as u64;
    let shifted_end = intron_end - shift_left as u64;

    let strand = match db_strand {
        1 | 2 => db_strand,
        _ if motif == 0 => 0,
        _ => 2 - (motif % 2), // 1/3/5 → 1 (+), 2/4/6 → 2 (-)
    };

    PreparedJunction {
        chr_idx,
        start_pos: shifted_start,
        end_pos: shifted_end,
        motif,
        shift_left,
        shift_right,
        strand,
    }
}

/// Sort a prepared junction list into STAR's post-dedup order and apply
/// the cross-strand deduplication that STAR does after its second sort
/// (`sjdbPrepare.cpp:141-192`).
///
/// STAR's first-pass (intra-strand) dedup collapses duplicate sjdb
/// entries from the same source at the same `(start, end, strand)`.
/// rustar-aligner's `SpliceJunctionDb` already deduplicates on that key at the
/// HashMap level, so those first-pass branches never trigger here; the
/// second-pass cross-strand collision dedup does.
///
/// Dedup rules when two surviving junctions share `(stored_start,
/// stored_end)` but have different strand assignments:
///
/// - Undefined strand vs defined strand → keep the defined-strand one.
/// - Both non-canonical → collapse to a single entry with strand = 0
///   (undefined).
/// - One canonical + one not → keep the canonical one.
/// - Both canonical but on correct vs wrong strand relative to motif —
///   keep the one whose strand matches `2 - motif % 2`.
pub fn sort_and_dedup(mut junctions: Vec<PreparedJunction>) -> Vec<PreparedJunction> {
    junctions.sort_by(|a, b| {
        a.stored_start()
            .cmp(&b.stored_start())
            .then_with(|| a.stored_end().cmp(&b.stored_end()))
    });

    let mut out: Vec<PreparedJunction> = Vec::with_capacity(junctions.len());
    for j in junctions {
        match out.last() {
            Some(last)
                if last.stored_start() == j.stored_start()
                    && last.stored_end() == j.stored_end() =>
            {
                if let Some(winner) = merge_cross_strand(last, &j) {
                    *out.last_mut().unwrap() = winner;
                }
                // else: keep `last` unchanged.
            }
            _ => out.push(j),
        }
    }
    out
}

/// Decide what `(stored_start, stored_end)` duplicate to keep.
/// Returns `Some(new)` to replace the stored entry, `None` to keep it.
/// For the "both non-canonical on opposite strands" case we keep the
/// existing entry but force its strand to 0 in-place (handled by the
/// caller via a special-case — represented here as returning a cloned
/// `old` with strand=0).
fn merge_cross_strand(old: &PreparedJunction, new: &PreparedJunction) -> Option<PreparedJunction> {
    // Strand 0 = undefined.
    if old.strand > 0 && new.strand == 0 {
        return None; // keep old
    }
    if old.strand == 0 && new.strand > 0 {
        return Some(new.clone()); // replace
    }
    // Both non-canonical → collapse to undefined strand on the old one.
    if old.motif == 0 && new.motif == 0 {
        let mut merged = old.clone();
        merged.strand = 0;
        return Some(merged);
    }
    // One canonical, one not: prefer canonical.
    if old.motif > 0 && new.motif == 0 {
        return None;
    }
    if old.motif == 0 && new.motif > 0 {
        return Some(new.clone());
    }
    // Both canonical with defined strands. Keep the one on the correct
    // strand for its motif (2 - motif % 2). If the old one is on the
    // correct strand, skip the new one; otherwise replace.
    let old_expected = 2 - (old.motif % 2);
    if old.strand == old_expected {
        None
    } else {
        Some(new.clone())
    }
}

/// Build the Gsj buffer — the concatenated splice-junction flanking
/// sequences STAR appends to the Genome binary at `genomeGenerate`.
///
/// Each junction contributes exactly `2 * sjdb_overhang + 1` bytes:
/// `sjdb_overhang` donor-side bases (from the genome position
/// `original_start - sjdb_overhang`), `sjdb_overhang` acceptor-side
/// bases (from `original_end + 1`), and one `GSJ_SPACING` spacer byte
/// at the end. Matches `sjdbPrepare.cpp:203-215`.
///
/// Callers supply junctions in sorted/deduplicated order (via
/// [`sort_and_dedup`]). `genome` is read through
/// `Genome::sequence[..n_genome_real]`.
pub fn build_gsj(
    junctions: &[PreparedJunction],
    genome: &Genome,
    n_genome_real: u64,
    sjdb_overhang: u32,
) -> Result<Vec<u8>, Error> {
    let overhang = sjdb_overhang as usize;
    let sjdb_length = 2 * overhang + 1;
    let forward = &genome.sequence[..n_genome_real as usize];
    let mut gsj = vec![GSJ_SPACING; junctions.len() * sjdb_length];

    for (i, pj) in junctions.iter().enumerate() {
        let donor_start = pj
            .original_start()
            .checked_sub(overhang as u64)
            .ok_or_else(|| sjdb_bounds_err(pj, overhang, "donor underflows"))?;
        let acceptor_start = pj.original_end() + 1;

        let d0 = donor_start as usize;
        let a0 = acceptor_start as usize;
        if d0 + overhang > forward.len() || a0 + overhang > forward.len() {
            return Err(sjdb_bounds_err(pj, overhang, "flank overruns genome"));
        }

        let base = i * sjdb_length;
        gsj[base..base + overhang].copy_from_slice(&forward[d0..d0 + overhang]);
        gsj[base + overhang..base + 2 * overhang].copy_from_slice(&forward[a0..a0 + overhang]);
        // base + 2*overhang = trailing spacer, already initialized.
    }

    Ok(gsj)
}

fn sjdb_bounds_err(pj: &PreparedJunction, overhang: usize, reason: &str) -> Error {
    Error::Index(format!(
        "sjdb flanking window ({} bp overhang) {} for junction at original coords {}..{} on chr {}",
        overhang,
        reason,
        pj.original_start(),
        pj.original_end(),
        pj.chr_idx,
    ))
}

/// STAR's strand character table (`sjdbPrepare.cpp:198`): strand 0→'.',
/// 1→'+', 2→'-'. Any other value is an internal error.
fn strand_char(strand: u8) -> char {
    match strand {
        0 => '.',
        1 => '+',
        2 => '-',
        _ => '?',
    }
}

/// Write STAR's `sjdbInfo.txt` (`sjdbPrepare.cpp:196-215`).
///
/// Format: one header line `<n>\t<sjdbOverhang>\n` followed by one row per
/// junction with six tab-separated fields:
/// `stored_start  stored_end  motif  shift_left  shift_right  strand\n`.
///
/// `junctions` must already be in STAR's final sorted/deduped order.
pub fn write_sjdb_info_tab(
    path: &Path,
    junctions: &[PreparedJunction],
    sjdb_overhang: u32,
) -> Result<(), Error> {
    let file = File::create(path).map_err(|e| Error::io(e, path))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "{}\t{}", junctions.len(), sjdb_overhang).map_err(|e| Error::io(e, path))?;
    for pj in junctions {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            pj.stored_start(),
            pj.stored_end(),
            pj.motif,
            pj.shift_left,
            pj.shift_right,
            pj.strand,
        )
        .map_err(|e| Error::io(e, path))?;
    }
    w.flush().map_err(|e| Error::io(e, path))?;
    Ok(())
}

/// Write STAR's `sjdbList.out.tab` (`sjdbPrepare.cpp:197-219`).
///
/// Format: one row per junction with four tab-separated fields:
/// `chr  <1-based start>  <1-based end>  <strand_char>\n`. Coordinates
/// are chromosome-local and always reflect the pre-shift ("original")
/// intron boundaries, regardless of motif type.
pub fn write_sjdb_list_out_tab(
    path: &Path,
    junctions: &[PreparedJunction],
    genome: &Genome,
) -> Result<(), Error> {
    let file = File::create(path).map_err(|e| Error::io(e, path))?;
    let mut w = BufWriter::new(file);
    for pj in junctions {
        let chr_start = genome.chr_start[pj.chr_idx];
        let name = &genome.chr_name[pj.chr_idx];
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            name,
            pj.original_start() - chr_start + 1,
            pj.original_end() - chr_start + 1,
            strand_char(pj.strand),
        )
        .map_err(|e| Error::io(e, path))?;
    }
    w.flush().map_err(|e| Error::io(e, path))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a single-chromosome Genome from a forward-strand byte slice.
    // Callers who only exercise compute_shifts / prepare_junction don't
    // need a valid reverse complement; we fill the RC half with padding.
    fn make_test_genome(forward: Vec<u8>) -> Genome {
        let n = forward.len() as u64;
        let mut seq = vec![5u8; (n * 2) as usize];
        seq[..forward.len()].copy_from_slice(&forward);
        Genome {
            sequence: seq,
            n_genome: n,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![n],
            chr_start: vec![0, n],
        }
    }

    #[test]
    fn shifts_no_repeat() {
        let mut f = vec![4u8; 300];
        f[100] = 0;
        f[101] = 1;
        f[102] = 0;
        f[103] = 1;
        f[104] = 0;
        f[105] = 2; // intron start
        f[106] = 3;
        f[203] = 0;
        f[204] = 2; // intron end
        f[205] = 3;
        f[206] = 2;
        f[207] = 3;
        f[208] = 2;
        f[209] = 3;
        let g = make_test_genome(f);
        let (l, r) = compute_shifts(&g, 105, 204, g.n_genome);
        assert_eq!(l, 0);
        assert_eq!(r, 0);
    }

    #[test]
    fn shifts_with_repeat_on_left() {
        let mut f = vec![4u8; 300];
        f[103] = 0;
        f[104] = 1;
        f[105] = 2;
        f[106] = 3;
        f[203] = 0;
        f[204] = 1;
        f[205] = 3;
        let g = make_test_genome(f);
        let (l, _r) = compute_shifts(&g, 105, 204, g.n_genome);
        assert_eq!(l, 2);
    }

    #[test]
    fn shifts_with_repeat_on_right() {
        let mut f = vec![4u8; 300];
        f[105] = 2;
        f[106] = 3;
        f[107] = 0;
        f[204] = 1;
        f[205] = 2;
        f[206] = 3;
        f[207] = 1;
        let g = make_test_genome(f);
        let (_l, r) = compute_shifts(&g, 105, 204, g.n_genome);
        assert_eq!(r, 2);
    }

    #[test]
    fn shifts_cap_at_255() {
        let g = make_test_genome(vec![0u8; 2000]);
        let (l, r) = compute_shifts(&g, 500, 1000, g.n_genome);
        assert_eq!(l, 255);
        assert_eq!(r, 255);
    }

    #[test]
    fn shifts_stop_at_n_base() {
        let mut f = vec![0u8; 300];
        f[100] = 4;
        let g = make_test_genome(f);
        let (l, _) = compute_shifts(&g, 150, 249, g.n_genome);
        assert_eq!(l, 49);
    }

    #[test]
    fn prepare_gt_ag_forward_no_repeat() {
        let mut f = vec![4u8; 400];
        f[100] = 0;
        f[101] = 0;
        f[102] = 2; // G
        f[103] = 3; // T
        f[198] = 0; // A
        f[199] = 2; // G
        f[200] = 1;
        f[201] = 1;
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 102, 199, 1, &g, g.n_genome);
        assert_eq!(pj.motif, 1); // GT/AG
        assert_eq!(pj.shift_left, 0);
        assert_eq!(pj.shift_right, 0);
        assert_eq!(pj.start_pos, 102);
        assert_eq!(pj.end_pos, 199);
        assert_eq!(pj.strand, 1);
    }

    #[test]
    fn prepare_dot_strand_derived_from_motif() {
        // Non-canonical motif with db_strand=0 → strand 0.
        let g = make_test_genome(vec![0u8; 300]);
        let pj = prepare_junction(0, 100, 200, 0, &g, g.n_genome);
        assert_eq!(pj.motif, 0);
        assert_eq!(pj.strand, 0);

        // GT/AG motif (forward) with db_strand=0 → strand 1.
        let mut f = vec![4u8; 300];
        f[100] = 2;
        f[101] = 3;
        f[199] = 0;
        f[200] = 2;
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 100, 200, 0, &g, g.n_genome);
        assert_eq!(pj.motif, 1);
        assert_eq!(pj.strand, 1);

        // CT/AC motif (reverse) with db_strand=0 → strand 2.
        let mut f = vec![4u8; 300];
        f[100] = 1;
        f[101] = 3;
        f[199] = 0;
        f[200] = 1;
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 100, 200, 0, &g, g.n_genome);
        assert_eq!(pj.motif, 2);
        assert_eq!(pj.strand, 2);
    }

    #[test]
    fn prepare_shift_left_applied_to_coords() {
        // One base of repeat to the left → shift_left=1, coords decrement by 1.
        let mut f = vec![4u8; 400];
        f[102] = 2;
        f[103] = 3;
        f[198] = 0;
        f[199] = 2;
        f[101] = 2; // repeat
        f[100] = 3; // break
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 102, 199, 1, &g, g.n_genome);
        assert_eq!(pj.shift_left, 1);
        assert_eq!(pj.start_pos, 101);
        assert_eq!(pj.end_pos, 198);
    }

    fn pj(
        chr_idx: usize,
        start: u64,
        end: u64,
        motif: u8,
        shift_left: u8,
        strand: u8,
    ) -> PreparedJunction {
        PreparedJunction {
            chr_idx,
            start_pos: start,
            end_pos: end,
            motif,
            shift_left,
            shift_right: 0,
            strand,
        }
    }

    #[test]
    fn stored_and_original_coords_differ_for_canonical_only() {
        let canon = pj(0, 100, 200, 1, 3, 1);
        // Canonical stores ORIGINAL (shifted + shift_left).
        assert_eq!(canon.stored_start(), 103);
        assert_eq!(canon.stored_end(), 203);
        assert_eq!(canon.original_start(), 103);
        assert_eq!(canon.original_end(), 203);

        let noncanon = pj(0, 100, 200, 0, 3, 0);
        // Non-canonical stores SHIFTED.
        assert_eq!(noncanon.stored_start(), 100);
        assert_eq!(noncanon.stored_end(), 200);
        // original_* still recovers pre-shift coords.
        assert_eq!(noncanon.original_start(), 103);
        assert_eq!(noncanon.original_end(), 203);
    }

    #[test]
    fn sort_orders_by_stored_coords() {
        let a = pj(0, 100, 200, 1, 0, 1); // stored 100..200
        let b = pj(0, 50, 150, 1, 0, 1); // stored 50..150
        let c = pj(0, 80, 180, 1, 0, 1); // stored 80..180
        let out = sort_and_dedup(vec![a.clone(), b.clone(), c.clone()]);
        assert_eq!(out, vec![b, c, a]);
    }

    #[test]
    fn dedup_prefers_defined_strand_over_undefined() {
        let defined = pj(0, 100, 200, 1, 0, 1);
        let undefined = pj(0, 100, 200, 0, 0, 0);
        let out = sort_and_dedup(vec![defined.clone(), undefined]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], defined);
    }

    #[test]
    fn dedup_collapses_two_non_canonical_to_undefined_strand() {
        // Two non-canonical at same stored coords, opposite strands → one
        // entry with strand = 0.
        let a = pj(0, 100, 200, 0, 0, 1);
        let b = pj(0, 100, 200, 0, 0, 2);
        let out = sort_and_dedup(vec![a, b]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, 0);
    }

    #[test]
    fn dedup_prefers_canonical_over_non_canonical() {
        // Motif=1 (canonical) beats motif=0 (non) at same stored coords.
        let canon = pj(0, 100, 200, 1, 0, 1);
        let noncan = pj(0, 100, 200, 0, 0, 2);
        let out = sort_and_dedup(vec![noncan, canon.clone()]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], canon);
    }

    #[test]
    fn build_gsj_canonical_no_shift() {
        // Single canonical junction at original [100..200] with 3bp overhang.
        // Donor bytes = forward[97..100] (3 bytes), acceptor bytes =
        // forward[201..204] (3 bytes), then 1 spacer.
        let mut f = vec![4u8; 300];
        // Overhang-donor region: 97,98,99 = 0,1,2 (A,C,G)
        f[97] = 0;
        f[98] = 1;
        f[99] = 2;
        // Intron body filled with noise (irrelevant)
        f[100] = 2; // GT
        f[101] = 3;
        f[199] = 0; // AG
        f[200] = 2;
        // Acceptor-flank bytes: 201,202,203 = 3,0,1 (T,A,C)
        f[201] = 3;
        f[202] = 0;
        f[203] = 1;
        let g = make_test_genome(f);
        let junction = pj(0, 100, 200, 1, 0, 1);
        let gsj = build_gsj(&[junction], &g, g.n_genome, 3).unwrap();
        assert_eq!(gsj, vec![0, 1, 2, 3, 0, 1, GSJ_SPACING]);
    }

    #[test]
    fn build_gsj_noncanonical_reverts_shift_for_extraction() {
        // Non-canonical junction whose STORED coords sit 2bp left of
        // the ORIGINAL. Flank extraction must use ORIGINAL (the
        // `original_start()` / `original_end()` helpers) so the same
        // bytes land in Gsj regardless of motif type.
        let mut f = vec![4u8; 300];
        f[98] = 0; // original_start - overhang = 102 - 4 = 98
        f[99] = 1;
        f[100] = 2;
        f[101] = 3;
        // intron garbage at 102..=199
        f[200] = 2;
        f[201] = 3;
        f[202] = 0;
        f[203] = 1;
        let g = make_test_genome(f);
        // start_pos/end_pos are SHIFTED (stored form for motif==0).
        // With shift_left=2, original_start = 102, original_end = 199.
        let junction = pj(0, 100, 197, 0, 2, 0);
        let gsj = build_gsj(&[junction], &g, g.n_genome, 4).unwrap();
        // Donor: forward[98..102] = [0,1,2,3]
        // Acceptor: forward[200..204] = [2,3,0,1]
        assert_eq!(gsj, vec![0, 1, 2, 3, 2, 3, 0, 1, GSJ_SPACING]);
    }

    #[test]
    fn build_gsj_multiple_junctions_concatenate() {
        let mut f = vec![4u8; 400];
        // Junction A: canonical at [100..200], overhang 2 → donor
        // bytes 98,99; acceptor 201,202.
        f[98] = 0;
        f[99] = 0;
        f[201] = 3;
        f[202] = 3;
        // Junction B: canonical at [300..350], overhang 2 → donor
        // bytes 298,299; acceptor 351,352.
        f[298] = 1;
        f[299] = 1;
        f[351] = 2;
        f[352] = 2;
        let g = make_test_genome(f);
        let a = pj(0, 100, 200, 1, 0, 1);
        let b = pj(0, 300, 350, 1, 0, 1);
        let gsj = build_gsj(&[a, b], &g, g.n_genome, 2).unwrap();
        // sjdb_length = 5 (2 donor + 2 acceptor + 1 spacer) per junction.
        assert_eq!(
            gsj,
            vec![
                0,
                0,
                3,
                3,
                GSJ_SPACING, // junction A
                1,
                1,
                2,
                2,
                GSJ_SPACING, // junction B
            ]
        );
    }

    #[test]
    fn build_gsj_errors_when_flank_underflows() {
        // Junction too close to the left edge → donor_start underflows.
        let g = make_test_genome(vec![0u8; 300]);
        let bad = pj(0, 3, 200, 1, 0, 1); // original_start=3; overhang 10 underflows
        let err = build_gsj(&[bad], &g, g.n_genome, 10).unwrap_err();
        assert!(format!("{}", err).contains("underflows"));
    }

    #[test]
    fn build_gsj_errors_when_flank_overruns_genome() {
        // Acceptor overruns the right edge.
        let g = make_test_genome(vec![0u8; 300]);
        // original_end = 295, overhang = 10 → acceptor region [296..306), exceeds n_genome=300.
        let bad = pj(0, 100, 295, 1, 0, 1);
        let err = build_gsj(&[bad], &g, g.n_genome, 10).unwrap_err();
        assert!(format!("{}", err).contains("overruns"));
    }

    #[test]
    fn write_sjdb_info_matches_star_format() {
        // Bytes exactly match STAR's `sjdbPrepare.cpp:200+215` format:
        // header "n\toverhang\n" then six tab-separated fields per row.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let junctions = vec![
            // Canonical: stored = original = 87387..87499, motif 1, shiftR 1, strand 1.
            PreparedJunction {
                chr_idx: 0,
                start_pos: 87387, // shifted — shift_left=0, so stored=original
                end_pos: 87499,
                motif: 1,
                shift_left: 0,
                shift_right: 1,
                strand: 1,
            },
            // Non-canonical: stored = shifted (139187..139217).
            PreparedJunction {
                chr_idx: 0,
                start_pos: 139187,
                end_pos: 139217,
                motif: 0,
                shift_left: 0,
                shift_right: 0,
                strand: 1,
            },
        ];
        write_sjdb_info_tab(tmp.path(), &junctions, 99).unwrap();
        let bytes = std::fs::read(tmp.path()).unwrap();
        assert_eq!(
            bytes,
            b"2\t99\n87387\t87499\t1\t0\t1\t1\n139187\t139217\t0\t0\t0\t1\n"
        );
    }

    #[test]
    fn write_sjdb_list_is_1_based_chr_local() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut genome = make_test_genome(vec![0u8; 200_000]);
        genome.chr_name = vec!["I".to_string()];
        // Canonical junction — stored == original — (motif>0 path).
        let canon = PreparedJunction {
            chr_idx: 0,
            start_pos: 87387,
            end_pos: 87499,
            motif: 1,
            shift_left: 0,
            shift_right: 0,
            strand: 1,
        };
        // Non-canonical with shift_left=3 — STAR writes
        // `stored + shift_left + 1`, which is `original + 1`.
        let noncan = PreparedJunction {
            chr_idx: 0,
            start_pos: 139184, // shifted
            end_pos: 139214,
            motif: 0,
            shift_left: 3,
            shift_right: 0,
            strand: 0,
        };
        write_sjdb_list_out_tab(tmp.path(), &[canon, noncan], &genome).unwrap();
        let bytes = std::fs::read(tmp.path()).unwrap();
        assert_eq!(bytes, b"I\t87388\t87500\t+\nI\t139188\t139218\t.\n");
    }

    #[test]
    fn write_sjdb_list_uses_chr_local_coords_for_second_chromosome() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        // Two chromosomes: chrA length 1000, chrB starts at 1000, length 1000.
        let mut seq = vec![5u8; 4000];
        seq[..2000].copy_from_slice(&vec![0u8; 2000]);
        let genome = Genome {
            sequence: seq,
            n_genome: 2000,
            n_chr_real: 2,
            chr_name: vec!["chrA".to_string(), "chrB".to_string()],
            chr_length: vec![1000, 1000],
            chr_start: vec![0, 1000, 2000],
        };
        let pj_b = PreparedJunction {
            chr_idx: 1,
            start_pos: 1500,
            end_pos: 1900,
            motif: 1,
            shift_left: 0,
            shift_right: 0,
            strand: 2,
        };
        write_sjdb_list_out_tab(tmp.path(), &[pj_b], &genome).unwrap();
        let bytes = std::fs::read(tmp.path()).unwrap();
        // chrB-local: 1500-1000+1 = 501, 1900-1000+1 = 901.
        assert_eq!(bytes, b"chrB\t501\t901\t-\n");
    }

    #[test]
    fn dedup_prefers_strand_matching_motif() {
        // motif=1 (GT/AG +) stored on wrong strand (2) vs a competing
        // motif=2 (CT/AC -) on its correct strand. Same stored coords.
        // STAR keeps the one whose strand matches `2 - motif%2`.
        // For motif=1: expected strand = 2 - 1%2 = 1.
        // For motif=2: expected strand = 2 - 2%2 = 2.
        let old_wrong = pj(0, 100, 200, 1, 0, 2); // motif 1 wants strand 1, has 2
        let new_right = pj(0, 100, 200, 2, 0, 2); // motif 2 wants strand 2, has 2
        let out = sort_and_dedup(vec![old_wrong, new_right.clone()]);
        assert_eq!(out.len(), 1);
        // `old_wrong` is on wrong strand for its motif — STAR replaces.
        assert_eq!(out[0], new_right);
    }
}
