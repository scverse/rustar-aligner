//! `--chimOutType WithinBAM`: chimeric records in the main BAM.
//!
//! A port of `ChimericAlign::chimericBAMoutput` and the mapped branch of
//! `ReadAlign::alignBAM`, which it calls once per chimeric segment. The segments
//! are STAR's whole transcripts ([`ChimBam`]), so a segment that covers both
//! mates of a pair produces two records, exactly as STAR's does.
//!
//! STAR writes these in place of the read's normal alignments
//! (`ReadAlign_oneRead.cpp:99`); the caller is responsible for that.

use crate::chimeric::detect::{ChimExon, WinTr};
use crate::chimeric::{ChimBam, ChimericAlignment};
use crate::error::Error;
use crate::genome::Genome;
use crate::io::fastq::{complement_base, decode_base};
use crate::io::sam::{apply_sam_flag_or_and, fastq_qual_to_phred, maybe_insert_rg_tag};
use crate::params::{Parameters, SamAttributes};
use bstr::BString;
use noodles::sam;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

/// STAR's `canonSJ` for the gap between mates.
const SJ_MATE_GAP: i32 = -3;
/// STAR's `canonSJ` for a deletion.
const SJ_DELETION: i32 = -1;

/// The read a chimera came from, as sequenced.
pub struct ChimReadInput<'a> {
    pub name: &'a str,
    /// Mates as sequenced (encoded bases), before clipping. `[1]` is empty for
    /// single-end.
    pub seq: [&'a [u8]; 2],
    /// FASTQ qualities (Phred+33) for the same bases.
    pub qual: [&'a [u8]; 2],
    /// Bases clipped before alignment, per mate: `[5', 3']`.
    pub clip: [[usize; 2]; 2],
}

/// STAR's `alignType` for a chimeric record.
#[derive(Clone, Copy, PartialEq, Eq)]
enum AlignType {
    /// `-10`: written as a normal alignment.
    Normal,
    /// `-11`: supplementary, hard-clipped on the left (junction on the left).
    HardLeft,
    /// `-12`: supplementary, hard-clipped on the right.
    HardRight,
    /// `-13`: supplementary, soft-clipped.
    Soft,
}

/// Mate information passed to `alignBAM` for a single-mate segment.
#[derive(Clone, Copy)]
struct MateRef {
    chr: usize,
    /// Absolute genome position of the mate's first base.
    start: u64,
    /// The mate aligns to the reverse strand.
    reverse: bool,
}

/// Every chimeric BAM record for one read, in STAR's order.
pub fn build_chimeric_bam_records(
    chims: &[ChimericAlignment],
    read: &ChimReadInput,
    genome: &Genome,
    params: &Parameters,
) -> Result<Vec<RecordBuf>, Error> {
    let mut out = Vec::new();
    for chim in chims {
        if let Some(b) = &chim.bam {
            out.extend(chimeric_bam_output(b, read, genome, params)?);
        }
    }
    Ok(out)
}

fn first_frag(t: &WinTr) -> u8 {
    t.exons[0].frag
}

fn last_frag(t: &WinTr) -> u8 {
    t.exons[t.exons.len() - 1].frag
}

/// `ChimericAlign::chimericBAMoutput` (`ChimericAlign_chimericBAMoutput.cpp`).
fn chimeric_bam_output(
    b: &ChimBam,
    read: &ChimReadInput,
    genome: &Genome,
    params: &Parameters,
) -> Result<Vec<RecordBuf>, Error> {
    let tr = &b.tr;
    // Which segment is written as the representative (normal) alignment.
    // chimType 1: one segment covers both mates; 2: the mates bracket the
    // junction (both written as normal); 3: both segments in one mate (SE).
    let (represent, chim_type): (Option<usize>, u8) = if first_frag(&tr[0]) != last_frag(&tr[0]) {
        (Some(0), 1)
    } else if first_frag(&tr[1]) != last_frag(&tr[1]) {
        (Some(1), 1)
    } else if first_frag(&tr[0]) != first_frag(&tr[1]) {
        (None, 2)
    } else {
        (Some(usize::from(b.max_score[0] <= b.max_score[1])), 3)
    };
    let hard_clip = params.chim_out_bam_hard_clip();

    let mut recs: Vec<RecordBuf> = Vec::new();
    let (mut i_repr, mut i_suppl): (Option<usize>, Option<usize>) = (None, None);
    for itr in 0..2 {
        let other = &tr[1 - itr];
        let (mate, align_type);
        if chim_type == 2 {
            mate = Some(MateRef {
                chr: other.chr,
                start: other.exons[0].g,
                reverse: other.str_ != other.exons[0].frag,
            });
            align_type = AlignType::Normal;
        } else if represent == Some(itr) {
            mate = None;
            align_type = AlignType::Normal;
            // For a two-mate representative, the chimerically split mate is
            // its second record.
            i_repr = Some(recs.len() + usize::from(first_frag(&tr[itr]) != first_frag(other)));
        } else {
            align_type = if !hard_clip {
                AlignType::Soft
            } else if itr % 2 == tr[itr].str_ as usize {
                AlignType::HardRight
            } else {
                AlignType::HardLeft
            };
            i_suppl = Some(recs.len());
            mate = if chim_type == 1 {
                // Mate info for the supplementary: the representative's mate
                // that is not this segment's.
                let r = &tr[represent.expect("chimType 1 has a representative")];
                let mut iex = 0;
                while iex < r.exons.len() - 1 && r.exons[iex].frag == first_frag(&tr[itr]) {
                    iex += 1;
                }
                Some(MateRef {
                    chr: r.chr,
                    start: r.exons[iex].g,
                    reverse: r.str_ != r.exons[iex].frag,
                })
            } else {
                None
            };
        }
        recs.extend(align_bam(
            &tr[itr],
            b.max_score[itr],
            b.n_mm[itr],
            b,
            mate,
            align_type,
            read,
            genome,
            params,
        )?);
    }

    // SA tags: the representative and the supplementary point at each other,
    // using the other record's final CIGAR, MAPQ and NM.
    if let (Some(r), Some(s)) = (i_repr, i_suppl) {
        let sa_for_repr = sa_tag(&recs[s], genome);
        let sa_for_suppl = sa_tag(&recs[r], genome);
        recs[r]
            .data_mut()
            .insert(Tag::new(b'S', b'A'), Value::String(sa_for_repr));
        recs[s]
            .data_mut()
            .insert(Tag::new(b'S', b'A'), Value::String(sa_for_suppl));
    }
    Ok(recs)
}

/// `chr,pos,strand,CIGAR,MAPQ,NM;` from a finished record, as STAR builds it
/// with `bam_cigarString` and the record's `NM` (always present under
/// WithinBAM, which adds the attribute).
fn sa_tag(rec: &RecordBuf, genome: &Genome) -> BString {
    use std::fmt::Write;
    let chr = rec
        .reference_sequence_id()
        .map_or("*", |i| genome.chr_name[i].as_str());
    let pos = rec.alignment_start().map_or(0, usize::from);
    let strand = if rec.flags().is_reverse_complemented() {
        '-'
    } else {
        '+'
    };
    let mut cigar = String::new();
    for op in rec.cigar().as_ref() {
        let c = match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        };
        let _ = write!(cigar, "{}{c}", op.len());
    }
    let mapq = rec.mapping_quality().map_or(255, u8::from);
    let nm = match rec.data().get(&Tag::EDIT_DISTANCE) {
        Some(v) => v.as_int().unwrap_or(0),
        None => 0,
    };
    format!("{chr},{pos},{strand},{cigar},{mapq},{nm};").into()
}

/// The combined read in a transcript's strand frame (`Read1[0]` or `Read1[2]`),
/// over the clipped mates the exons index.
fn strand_frame(read: &ChimReadInput, b: &ChimBam, str_: u8) -> Vec<u8> {
    let clipped = |m: usize| -> &[u8] {
        let s = read.seq[m];
        let [c5, c3] = read.clip[m];
        if s.len() < c5 + c3 {
            &[]
        } else {
            &s[c5..s.len() - c3]
        }
    };
    let mut fwd: Vec<u8> = clipped(0).to_vec();
    if b.paired {
        fwd.push(crate::align::stitch::PE_SPACER_BASE);
        fwd.extend(clipped(1).iter().rev().map(|&x| complement_base(x)));
    }
    if str_ == 0 {
        fwd
    } else {
        fwd.iter().rev().map(|&x| complement_base(x)).collect()
    }
}

/// `ReadAlign::alignBAM` for a mapped chimeric transcript: one record per mate
/// the transcript covers.
#[allow(clippy::too_many_arguments)]
fn align_bam(
    t: &WinTr,
    max_score: i32,
    n_mm: u32,
    b: &ChimBam,
    mate: Option<MateRef>,
    align_type: AlignType,
    read: &ChimReadInput,
    genome: &Genome,
    params: &Parameters,
) -> Result<Vec<RecordBuf>, Error> {
    let n_ex = t.exons.len();
    let (mut n_mates, mut i_ex_mate) = (1usize, n_ex - 1);
    for i in 0..n_ex - 1 {
        if t.canon_sj[i] == SJ_MATE_GAP {
            n_mates = 2;
            i_ex_mate = i;
            break;
        }
    }
    let rl = b.read_length;
    let rl_orig = [read.seq[0].len(), read.seq[1].len()];
    let left_mate = if b.paired { t.str_ as usize } else { 0 };
    let chr_start = genome.chr_start[t.chr];
    let frame = strand_frame(read, b, t.str_);
    let attrs = params.out_sam_attributes;
    let rg_owned = params.primary_rg_id()?;

    let mapq = match b.chim_n {
        n if n >= 5 => 0,
        n if n >= 3 => 1,
        2 => 3,
        _ => params.out_sam_mapq_unique,
    };

    let mut out = Vec::with_capacity(n_mates);
    for imate in 0..n_mates {
        let (i_ex1, i_ex2) = if imate == 0 {
            (0, i_ex_mate)
        } else {
            (i_ex_mate + 1, n_ex - 1)
        };
        let e1 = t.exons[i_ex1];
        let e2 = t.exons[i_ex2];
        let mate_idx = e1.frag as usize;
        let str_ = u16::from(t.str_);

        // FLAG
        let mut flag: u16 = 0;
        if b.paired {
            flag = 0x1;
            if i_ex_mate == n_ex - 1 {
                if mate.is_none() {
                    flag |= 0x8;
                }
            } else {
                flag |= 0x2;
            }
        }
        if align_type != AlignType::Normal {
            flag |= 0x800;
        }
        if !b.is_best {
            flag |= 0x100;
        }
        if mate_idx == 0 {
            flag |= str_ * 0x10;
            if n_mates == 2 {
                flag |= (1 - str_) * 0x20;
            }
        } else {
            flag |= (1 - str_) * 0x10;
            if n_mates == 2 {
                flag |= str_ * 0x20;
            }
        }
        if b.paired {
            flag |= if mate_idx == 0 { 0x40 } else { 0x80 };
            if n_mates == 1 && mate.is_some_and(|m| m.reverse) {
                flag |= 0x20;
            }
        }

        // Clipped bases on this record's left, in genome orientation.
        let trim_l = match (t.str_, mate_idx) {
            (0, 0) => read.clip[0][0],
            (0, _) => read.clip[1][1],
            (_, 0) => read.clip[0][1],
            (_, _) => read.clip[1][0],
        };

        // CIGAR
        let mut cigar: Vec<Op> = Vec::new();
        let trim_l1 = trim_l + e1.r
            - if e1.r < rl[left_mate] {
                0
            } else {
                rl[left_mate] + 1
            };
        if trim_l1 > 0 {
            let k = if align_type == AlignType::HardLeft {
                Kind::HardClip
            } else {
                Kind::SoftClip
            };
            cigar.push(Op::new(k, trim_l1));
        }
        let mut sj_motif: Vec<i8> = Vec::new();
        let mut sj_intron: Vec<i32> = Vec::new();
        for ii in i_ex1..=i_ex2 {
            let x = t.exons[ii];
            if ii > i_ex1 {
                let p = t.exons[ii - 1];
                let gap_g = x.g - (p.g + p.l as u64);
                let gap_r = x.r - p.r - p.l;
                if gap_r > 0 {
                    cigar.push(Op::new(Kind::Insertion, gap_r));
                }
                if t.canon_sj[ii - 1] >= 0 || t.sj_annot[ii - 1] {
                    cigar.push(Op::new(Kind::Skip, gap_g as usize));
                    sj_motif
                        .push((t.canon_sj[ii - 1] + if t.sj_annot[ii - 1] { 20 } else { 0 }) as i8);
                    sj_intron.push((p.g + p.l as u64 + 1 - chr_start) as i32);
                    sj_intron.push((x.g - chr_start) as i32);
                } else if gap_g > 0 {
                    cigar.push(Op::new(Kind::Deletion, gap_g as usize));
                }
            }
            if x.l > 0 {
                cigar.push(Op::new(Kind::Match, x.l));
            }
        }
        if sj_motif.is_empty() {
            sj_motif.push(-1);
            sj_intron.push(-1);
        }
        let trim_r1 = (if e1.r < rl[left_mate] {
            rl_orig[left_mate]
        } else {
            rl[left_mate] + 1 + rl_orig[mate_idx]
        }) as i64
            - (e2.r + e2.l) as i64
            - trim_l as i64;
        if trim_r1 > 0 {
            let k = if align_type == AlignType::HardRight {
                Kind::HardClip
            } else {
                Kind::SoftClip
            };
            cigar.push(Op::new(k, trim_r1 as usize));
        }

        let mut rec = RecordBuf::default();
        rec.name_mut().replace(read.name.into());
        *rec.flags_mut() = sam::alignment::record::Flags::from(flag);
        *rec.reference_sequence_id_mut() = Some(t.chr);
        let pos = (e1.g - chr_start + 1) as usize;
        *rec.alignment_start_mut() = Some(
            pos.try_into()
                .map_err(|e| Error::Alignment(format!("invalid chimeric position {pos}: {e}")))?,
        );
        *rec.mapping_quality_mut() = MappingQuality::new(mapq);
        *rec.cigar_mut() = cigar.into_iter().collect();

        // RNEXT / PNEXT / TLEN
        if n_mates > 1 {
            *rec.mate_reference_sequence_id_mut() = Some(t.chr);
            let other = t.exons[if imate == 0 { i_ex_mate + 1 } else { 0 }];
            let mpos = (other.g - chr_start + 1) as usize;
            *rec.mate_alignment_start_mut() = mpos.try_into().ok();
            let last = t.exons[n_ex - 1];
            let tlen = if params.out_sam_tlen == 2 {
                let mate1_end = t.exons[i_ex_mate].g + t.exons[i_ex_mate].l as u64;
                let end = (last.g + last.l as u64).max(mate1_end);
                let start = t.exons[0].g.min(t.exons[i_ex_mate + 1].g);
                let left_most = usize::from(t.exons[0].g > t.exons[i_ex_mate + 1].g);
                let v = (end - start) as i32;
                if imate == left_most { v } else { -v }
            } else {
                let v = (last.g + last.l as u64 - t.exons[0].g) as i32;
                if imate == 0 { v } else { -v }
            };
            *rec.template_length_mut() = tlen;
        } else if let Some(m) = mate {
            *rec.mate_reference_sequence_id_mut() = Some(m.chr);
            let mpos = (m.start - genome.chr_start[m.chr] + 1) as usize;
            *rec.mate_alignment_start_mut() = mpos.try_into().ok();
        }

        // SEQ / QUAL: the mate as sequenced, reverse-complemented when this
        // record is on the other strand, then trimmed where hard-clipped.
        let seq = read.seq[mate_idx];
        let qual = read.qual[mate_idx];
        let same_strand = mate_idx == t.str_ as usize;
        let mut seq_out: Vec<u8> = if same_strand {
            seq.iter().map(|&x| decode_base(x)).collect()
        } else {
            seq.iter()
                .rev()
                .map(|&x| decode_base(complement_base(x)))
                .collect()
        };
        let mut qual_out = fastq_qual_to_phred(qual);
        if !same_strand {
            qual_out.reverse();
        }
        match align_type {
            AlignType::HardLeft => {
                let k = trim_l1.min(seq_out.len());
                seq_out.drain(..k);
                if qual_out.len() >= k {
                    qual_out.drain(..k);
                }
            }
            AlignType::HardRight => {
                let k = (trim_r1.max(0) as usize).min(seq_out.len());
                seq_out.truncate(seq_out.len() - k);
                if qual_out.len() >= k {
                    qual_out.truncate(qual_out.len() - k);
                }
            }
            _ => {}
        }
        *rec.sequence_mut() = Sequence::from(seq_out);
        *rec.quality_scores_mut() = QualityScores::from(qual_out);

        // Attributes, in the order the main writer uses.
        let (nm, md) = nm_md(t, i_ex1, i_ex2, &frame, genome);
        let data = rec.data_mut();
        if attrs.contains(SamAttributes::NH) {
            data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(b.chim_n as i32));
        }
        if attrs.contains(SamAttributes::HI) {
            data.insert(
                Tag::HIT_INDEX,
                Value::from(b.i_tr as i32 + params.out_sam_attr_ih_start as i32),
            );
        }
        if attrs.contains(SamAttributes::AS) {
            data.insert(Tag::ALIGNMENT_SCORE, Value::from(max_score));
        }
        if attrs.contains(SamAttributes::NMM) {
            data.insert(Tag::new(b'n', b'M'), Value::from(n_mm as i32));
        }
        if attrs.contains(SamAttributes::NM) {
            data.insert(Tag::EDIT_DISTANCE, Value::from(nm as i32));
        }
        if attrs.contains(SamAttributes::XS) {
            let (p, m) = (t.intron_motifs[1], t.intron_motifs[2]);
            if p > 0 && m == 0 {
                data.insert(Tag::new(b'X', b'S'), Value::Character(b'+'));
            } else if m > 0 && p == 0 {
                data.insert(Tag::new(b'X', b'S'), Value::Character(b'-'));
            }
        }
        if attrs.contains(SamAttributes::JM) {
            data.insert(Tag::new(b'j', b'M'), Value::Array(Array::Int8(sj_motif)));
        }
        if attrs.contains(SamAttributes::JI) {
            data.insert(Tag::new(b'j', b'I'), Value::Array(Array::Int32(sj_intron)));
        }
        if attrs.contains(SamAttributes::MD) {
            data.insert(Tag::new(b'M', b'D'), Value::String(md.into()));
        }
        maybe_insert_rg_tag(&mut rec, rg_owned.as_deref());
        apply_sam_flag_or_and(&mut rec, params);
        out.push(rec);
    }
    Ok(out)
}

/// `ReadAlign::samAttrNM_MD`: edit distance and MD string over one mate's
/// exons. Unlike `alignScore`, an `N` on either side counts as a mismatch.
fn nm_md(t: &WinTr, i_ex1: usize, i_ex2: usize, frame: &[u8], genome: &Genome) -> (u32, String) {
    use std::fmt::Write;
    let nt = |b: u8| -> char {
        match b {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        }
    };
    let mut md = String::new();
    let (mut match_n, mut n_mm, mut n_i, mut n_d) = (0u32, 0u32, 0u32, 0u32);
    for iex in i_ex1..=i_ex2 {
        let x: ChimExon = t.exons[iex];
        for ii in 0..x.l {
            let r1 = frame.get(x.r + ii).copied().unwrap_or(4);
            let g1 = genome.get_base(x.g + ii as u64).unwrap_or(4);
            if r1 != g1 || r1 == 4 || g1 == 4 {
                n_mm += 1;
                let _ = write!(md, "{match_n}{}", nt(g1));
                match_n = 0;
            } else {
                match_n += 1;
            }
        }
        if iex < i_ex2 {
            let y = t.exons[iex + 1];
            if t.canon_sj[iex] < 0 {
                n_d += (y.g - (x.g + x.l as u64)) as u32;
            }
            n_i += (y.r - x.r - x.l) as u32;
            if t.canon_sj[iex] == SJ_DELETION {
                let _ = write!(md, "{match_n}^");
                for g in x.g + x.l as u64..y.g {
                    md.push(nt(genome.get_base(g).unwrap_or(4)));
                }
                match_n = 0;
            }
        }
    }
    let _ = write!(md, "{match_n}");
    (n_mm + n_i + n_d, md)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chimeric::detect::{ChimExon, WinTr};

    const CHR_PAD: u64 = 1024;

    /// Two 1000-base chromosomes of pseudo-random sequence.
    fn genome() -> Genome {
        let n_genome = CHR_PAD * 2;
        let mut seq = vec![5u8; 2 * n_genome as usize];
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        for chr in 0..2u64 {
            for i in 0..1000u64 {
                x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                seq[(chr * CHR_PAD + i) as usize] = ((x >> 33) % 4) as u8;
            }
        }
        Genome {
            transform_blocks: None,
            sequence: seq.into(),
            n_genome,
            n_genome_real: n_genome,
            n_chr_real: 2,
            chr_name: vec!["chr0".into(), "chr1".into()],
            chr_length: vec![1000, 1000],
            chr_start: vec![0, CHR_PAD, n_genome],
        }
    }

    fn slice(g: &Genome, start: u64, len: usize) -> Vec<u8> {
        (0..len as u64)
            .map(|i| g.get_base(start + i).unwrap())
            .collect()
    }

    fn params() -> Parameters {
        Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "12",
            "--chimOutType",
            "WithinBAM",
            "--outSAMtype",
            "BAM",
            "Unsorted",
        ])
    }

    fn tr(chr: usize, str_: u8, g: u64, r: usize, l: usize, frag: u8) -> WinTr {
        WinTr {
            chr,
            str_,
            exons: vec![ChimExon { g, r, l, frag }],
            canon_sj: Vec::new(),
            sj_annot: Vec::new(),
            max_score: 0,
            intron_motifs: [0; 3],
            r_length: l,
            ro_start: r,
        }
    }

    fn tag(rec: &RecordBuf, t: [u8; 2]) -> Option<&Value> {
        rec.data().get(&Tag::new(t[0], t[1]))
    }

    fn cigar_string(rec: &RecordBuf) -> String {
        let mut s = String::new();
        for op in rec.cigar().as_ref() {
            let c = match op.kind() {
                Kind::Match => 'M',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                _ => '?',
            };
            let _ = std::fmt::Write::write_fmt(&mut s, format_args!("{}{c}", op.len()));
        }
        s
    }

    /// A single-end chimera, half chr0 and half chr1: the higher-scoring
    /// segment is written as the normal record, the other as a hard-clipped
    /// supplementary with only its own bases in SEQ, and each SA tag carries
    /// the other record's final CIGAR.
    #[test]
    fn single_end_chimera_is_representative_plus_hard_clipped_supplementary() {
        let g = genome();
        let mut read = slice(&g, 100, 60);
        read.extend(slice(&g, CHR_PAD + 500, 60));
        let qual = vec![b'I'; 120];
        let b = ChimBam {
            tr: [tr(0, 0, 100, 0, 60, 0), tr(1, 0, CHR_PAD + 500, 60, 60, 0)],
            max_score: [57, 59],
            n_mm: [0, 0],
            i_tr: 0,
            chim_n: 1,
            is_best: true,
            read_length: [120, 0],
            paired: false,
        };
        let input = ChimReadInput {
            name: "r",
            seq: [&read, &[]],
            qual: [&qual, &[]],
            clip: [[0, 0], [0, 0]],
        };
        let recs = chimeric_bam_output(&b, &input, &g, &params()).unwrap();
        assert_eq!(recs.len(), 2);
        // Segment 0 scores lower, so it is the supplementary; it leads in the
        // read and is forward, so its right end (the junction) is hard-clipped.
        let (suppl, repr) = (&recs[0], &recs[1]);
        assert_eq!(u16::from(suppl.flags()), 0x800);
        assert_eq!(cigar_string(suppl), "60M60H");
        assert_eq!(suppl.sequence().len(), 60);
        assert_eq!(u16::from(repr.flags()), 0);
        assert_eq!(cigar_string(repr), "60S60M");
        assert_eq!(repr.sequence().len(), 120);
        // WithinBAM adds NM; SA points at the other record's final CIGAR.
        assert_eq!(tag(repr, *b"NM"), Some(&Value::from(0)));
        assert_eq!(
            tag(repr, *b"SA"),
            Some(&Value::String("chr0,101,+,60M60H,255,0;".into()))
        );
        assert_eq!(
            tag(suppl, *b"SA"),
            Some(&Value::String("chr1,501,+,60S60M,255,0;".into()))
        );
        assert_eq!(tag(repr, *b"AS"), Some(&Value::from(59)));
    }

    /// Mates on two chromosomes bracket the junction (chimType 2): both are
    /// written as normal paired records pointing at each other, with no SA.
    #[test]
    fn bracketing_mates_are_two_normal_records_without_sa() {
        let g = genome();
        let m1 = slice(&g, 100, 60);
        let m2_rc = slice(&g, CHR_PAD + 500, 60);
        let m2: Vec<u8> = m2_rc.iter().rev().map(|&x| 3 - x).collect();
        let q = vec![b'I'; 60];
        let b = ChimBam {
            tr: [tr(0, 0, 100, 0, 60, 0), tr(1, 0, CHR_PAD + 500, 61, 60, 1)],
            max_score: [59, 59],
            n_mm: [0, 0],
            i_tr: 0,
            chim_n: 1,
            is_best: true,
            read_length: [60, 60],
            paired: true,
        };
        let input = ChimReadInput {
            name: "r",
            seq: [&m1, &m2],
            qual: [&q, &q],
            clip: [[0, 0], [0, 0]],
        };
        let recs = chimeric_bam_output(&b, &input, &g, &params()).unwrap();
        assert_eq!(recs.len(), 2);
        // Mate1 forward, mate reverse: 0x1|0x20|0x40. Mate2 reverse: 0x1|0x10|0x80.
        assert_eq!(u16::from(recs[0].flags()), 0x1 | 0x20 | 0x40);
        assert_eq!(u16::from(recs[1].flags()), 0x1 | 0x10 | 0x80);
        assert_eq!(recs[0].mate_reference_sequence_id(), Some(1));
        assert_eq!(recs[0].mate_alignment_start().map(usize::from), Some(501));
        assert_eq!(recs[1].mate_reference_sequence_id(), Some(0));
        assert_eq!(recs[1].mate_alignment_start().map(usize::from), Some(101));
        assert!(recs.iter().all(|r| tag(r, *b"SA").is_none()));
        assert_eq!(cigar_string(&recs[0]), "60M");
        assert_eq!(cigar_string(&recs[1]), "60M");
    }

    /// `SoftClip` keeps the supplementary's clipped bases in SEQ.
    #[test]
    fn soft_clip_mode_keeps_the_whole_read() {
        let g = genome();
        let mut read = slice(&g, 100, 60);
        read.extend(slice(&g, CHR_PAD + 500, 60));
        let qual = vec![b'I'; 120];
        let b = ChimBam {
            tr: [tr(0, 0, 100, 0, 60, 0), tr(1, 0, CHR_PAD + 500, 60, 60, 0)],
            max_score: [57, 59],
            n_mm: [0, 0],
            i_tr: 0,
            chim_n: 1,
            is_best: true,
            read_length: [120, 0],
            paired: false,
        };
        let input = ChimReadInput {
            name: "r",
            seq: [&read, &[]],
            qual: [&qual, &[]],
            clip: [[0, 0], [0, 0]],
        };
        let p = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "12",
            "--chimOutType",
            "WithinBAM",
            "SoftClip",
            "--outSAMtype",
            "BAM",
            "Unsorted",
        ]);
        let recs = chimeric_bam_output(&b, &input, &g, &p).unwrap();
        assert_eq!(u16::from(recs[0].flags()), 0x800);
        assert_eq!(cigar_string(&recs[0]), "60M60S");
        assert_eq!(recs[0].sequence().len(), 120);
    }
}
