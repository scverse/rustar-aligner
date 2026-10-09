//! Paired-end mate-overlap merging (`--peOverlapNbasesMin`), STAR's `peOverlapMergeMap`.
//!
//! When the two mates of a fragment overlap in genome space (fragment shorter than
//! `mate1_len + mate2_len`), STAR merges them into one contiguous single-end read, aligns that
//! merged read, and converts the resulting single-end transcript back into a two-mate
//! paired-end transcript. This typically yields a cleaner alignment across the overlap region
//! than aligning the mates separately, and STAR unconditionally replaces the separate-mate
//! result with the converted one whenever the merge succeeds and produces at least one
//! convertible transcript.
//!
//! `local_search_n_is_mm` and `pe_merge_mates` are near-verbatim ports of STAR's
//! `localSearchNisMM` (`SequenceFuns.cpp`) and `ReadAlign::peMergeMates`
//! (`ReadAlign_peOverlapMergeMap.cpp`), taken from the sister project STAR-rs
//! (`crates/star-align/src/star_pe_overlap.rs`), which cites the exact upstream line ranges.
//!
//! The SE->PE conversion (`convert_merged_transcript_to_pe`) is NOT a verbatim port: STAR-rs
//! represents a transcript as parallel arrays of relative exon offsets plus a `canon_sj`
//! sentinel array (`-3` marks a mate boundary within one combined transcript), whereas
//! rustar-aligner's [`Transcript`] stores absolute `genome_start`/`genome_end`/`read_start`/
//! `read_end` per [`Exon`] and represents a paired alignment as *two separate* `Transcript`s
//! (`i_frag` on each `Exon` is always 0 within a single mate's own transcript; the PE pairing
//! layer re-tags mate2's exons to `i_frag = 1` only when projecting to a transcriptome BAM).
//! This module reimplements STAR's `Transcript::peOverlapSEtoPE` /
//! `ReadAlign::peOverlapSEtoPE` (`ReadAlign_peOverlapMergeMap.cpp:136-306`) against those
//! actual types instead, and rescores each mate from scratch (mirroring STAR's
//! `Transcript::alignScore`) rather than reusing the WorkingTranscript-based finalization
//! path in `stitch.rs` (that path expects raw SA-space genome coordinates pre-dating
//! `sa_pos_to_forward`, which the already-finalized merged transcript no longer has).

use crate::align::score::{AlignmentScorer, SpliceMotif};
use crate::align::transcript::{Exon, Transcript};
use crate::genome::Genome;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};

/// STAR `localSearchNisMM` (`SequenceFuns.cpp:317-339`): the best start offset of `y` within `x`
/// (slides `y` over `x`), accepting an offset only if the overlap's mismatch/match ratio is
/// within `p_mm`. An `N` (code > 3) in `x` or `y` always counts as a mismatch. Returns `x.len()`
/// if no offset qualifies. Ties (equal matches) break to fewer mismatches.
pub fn local_search_n_is_mm(x: &[u8], y: &[u8], p_mm: f64) -> usize {
    let (nx, ny) = (x.len(), y.len());
    let (mut n_match_best, mut n_mm_best, mut ix_best) = (0usize, 0usize, nx);
    for ix in 0..nx {
        let (mut n_match, mut n_mm) = (0usize, 0usize);
        for iy in 0..ny.min(nx - ix) {
            if x[ix + iy] == y[iy] && y[iy] < 4 {
                n_match += 1;
            } else {
                n_mm += 1;
            }
        }
        // `n_mm / n_match` with n_match == 0 is +inf/NaN, which is never <= p_mm (matches C++ IEEE).
        if (n_match > n_match_best || (n_match == n_match_best && n_mm < n_mm_best))
            && (n_mm as f64) / (n_match as f64) <= p_mm
        {
            ix_best = ix;
            n_match_best = n_match;
            n_mm_best = n_mm;
        }
    }
    ix_best
}

/// The result of merging two overlapping mates into one single-end read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Merge {
    /// The merged single-end read (numericized), length `l1 + l2 - n_ov`.
    pub merged: Vec<u8>,
    /// Where each mate starts in the merged read (`peOv.mateStart`).
    pub mate_start: [u64; 2],
    /// The overlap length (`peOv.nOv`).
    pub n_ov: u64,
}

/// STAR `ReadAlign::peMergeMates` (`ReadAlign_peOverlapMergeMap.cpp:79-134`): given mate1 (`m1`,
/// numericized) and the reverse-complement of mate2 (`m2rc`, numericized, as it sits in the
/// combined read), find the overlap and build the merged single-end read. Returns `None` when
/// the overlap is shorter than `nbases_min` (peOverlap does nothing).
pub fn pe_merge_mates(m1: &[u8], m2rc: &[u8], nbases_min: u64, mmp: f64) -> Option<Merge> {
    let (l1, l2) = (m1.len(), m2rc.len());
    let s1 = local_search_n_is_mm(m1, m2rc, mmp); // rc(mate2) offset within mate1
    let s0 = local_search_n_is_mm(m2rc, m1, mmp); // mate1 offset within rc(mate2)
    let o1 = l2.min(l1.saturating_sub(s1));
    let o0 = l1.min(l2.saturating_sub(s0));
    let n_ov = o0.max(o1);
    if (n_ov as u64) < nbases_min {
        return None;
    }
    // `o1 >= o0`: mate2 sits at/after mate1 -> merged = mate1 ++ tail_of_rc(mate2).
    // else: mate1 sits inside mate2 -> merged = rc(mate2) ++ tail_of_mate1.
    let (merged, mate_start) = if o1 >= o0 {
        let mut merged = m1.to_vec();
        merged.extend_from_slice(&m2rc[o1..]);
        (merged, [0u64, s1 as u64])
    } else {
        let mut merged = m2rc.to_vec();
        merged.extend_from_slice(&m1[o0..]);
        (merged, [s0 as u64, 0u64])
    };
    debug_assert_eq!(merged.len(), l1 + l2 - n_ov);
    Some(Merge {
        merged,
        mate_start,
        n_ov: n_ov as u64,
    })
}

/// STAR's `MAX_N_EXONS` (`IncludeDefine.h`): a transcript with more exons than this cannot be
/// converted (mirrors STAR-rs's `star_pe_overlap::pe_overlap_se_to_pe`, which returns `None`
/// past this cap rather than silently truncating).
const MAX_N_EXONS: usize = 20;

/// Reverse-complement a numericized sequence (codes 0-3 complemented as `3 - b`; codes >= 4
/// such as `N` pass through unchanged).
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| if b < 4 { 3 - b } else { b })
        .collect()
}

/// Build the CIGAR, rescore, and recount mismatches/gaps/junctions for one mate's already-split
/// exon list, mirroring STAR's `Transcript::alignScore` (`Transcript_alignScore.cpp`) but reusing
/// [`AlignmentScorer`]'s own gap-scoring constants/helpers so the numbers stay consistent with
/// every other scored `Transcript` in this crate. `junctions` holds the (motif, annotated) pair
/// for each of *this mate's* interior splice-junction gaps, in left-to-right order (copied
/// straight from the original merged transcript rather than re-detected).
///
/// `exons` must be sorted ascending by `read_start` (equivalently `genome_start`) and use
/// coordinates local to this mate (`read_start`/`read_end` in `0..mate_len`, in the same
/// orientation as `native_read`; `genome_start`/`genome_end` absolute forward-genome).
#[allow(clippy::too_many_arguments)]
fn score_mate_exons(
    exons: &[Exon],
    junctions: &[(SpliceMotif, bool)],
    mate_len: usize,
    native_read: &[u8],
    genome: &Genome,
    scorer: &AlignmentScorer,
) -> (Vec<Op>, i32, u32, u32, u32) {
    let mut cigar_ops: Vec<Op> = Vec::new();
    let mut score = 0i32;
    let mut n_mismatch = 0u32;
    let mut n_gap = 0u32;
    let mut n_junction = 0u32;
    let mut jidx = 0usize;

    let left_clip = exons[0].read_start;
    if left_clip > 0 {
        cigar_ops.push(Op::new(Kind::SoftClip, left_clip));
    }

    for (i, ex) in exons.iter().enumerate() {
        if i > 0 {
            let prev = &exons[i - 1];
            let read_gap = ex.read_start as i64 - prev.read_end as i64;
            let genome_gap = ex.genome_start as i64 - prev.genome_end as i64;
            if read_gap > 0 {
                cigar_ops.push(Op::new(Kind::Insertion, read_gap as usize));
                score += scorer.score_ins_open + scorer.score_ins_base * read_gap as i32;
                n_gap += 1;
            } else if genome_gap > 0 {
                if genome_gap as u32 >= scorer.align_intron_min {
                    cigar_ops.push(Op::new(Kind::Skip, genome_gap as usize));
                    let (motif, annotated) = junctions
                        .get(jidx)
                        .copied()
                        .unwrap_or((SpliceMotif::NonCanonical, false));
                    jidx += 1;
                    score += scorer
                        .score_annotated_junction(scorer.score_splice_junction(motif), annotated);
                    n_junction += 1;
                } else {
                    cigar_ops.push(Op::new(Kind::Deletion, genome_gap as usize));
                    score += scorer.score_del_open + scorer.score_del_base * genome_gap as i32;
                    n_gap += 1;
                }
            }
            // read_gap == 0 && genome_gap == 0 cannot happen for two distinct kept exons.
        }

        let len = ex.read_end - ex.read_start;
        if let Some(last) = cigar_ops.last_mut()
            && last.kind() == Kind::Match
        {
            *last = Op::new(Kind::Match, last.len() + len);
        } else {
            cigar_ops.push(Op::new(Kind::Match, len));
        }

        for ii in 0..len {
            let read_base = native_read.get(ex.read_start + ii).copied().unwrap_or(4);
            let genome_base = genome.get_base(ex.genome_start + ii as u64).unwrap_or(4);
            if read_base > 3 || genome_base > 3 {
                // N (or out-of-bounds/padding): no score impact, matches STAR's convention.
            } else if read_base == genome_base {
                score += 1;
            } else {
                score -= 1;
                n_mismatch += 1;
            }
        }
    }

    let right_clip = mate_len - exons.last().unwrap().read_end;
    if right_clip > 0 {
        cigar_ops.push(Op::new(Kind::SoftClip, right_clip));
    }

    let genomic_span = exons.last().unwrap().genome_end - exons[0].genome_start;
    let final_score = (score + scorer.genomic_length_penalty(genomic_span)).max(0);

    (cigar_ops, final_score, n_mismatch, n_gap, n_junction)
}

/// STAR `Transcript::peOverlapSEtoPE` / `ReadAlign::peOverlapSEtoPE`
/// (`ReadAlign_peOverlapMergeMap.cpp:136-306`): convert one merged-single-end `Transcript`
/// (aligned against `merge.merged`) back into a two-mate paired alignment, rescoring each mate
/// from scratch against its own read frame. Returns `None` when the transcript cannot be validly
/// split (a mate ends up with no exons, or either mate would exceed [`MAX_N_EXONS`]).
///
/// `mate1_seq`/`mate2_seq` are the *original* (non-reverse-complemented) mate sequences, exactly
/// as `align_paired_read` receives them.
pub fn convert_merged_transcript_to_pe(
    merged: &Transcript,
    merge: &Merge,
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    genome: &Genome,
    scorer: &AlignmentScorer,
) -> Option<(Transcript, Transcript)> {
    if merged.exons.is_empty() {
        return None;
    }
    let (l1, l2) = (mate1_seq.len() as u64, mate2_seq.len() as u64);
    let merged_lread = merge.merged.len() as u64;
    let s = usize::from(merged.is_reverse);
    let read_length = [l1, l2];
    let m_len = [read_length[s], read_length[1 - s]];
    let mut m_sta = merge.mate_start;
    if merged.is_reverse {
        m_sta[0] = merged_lread - read_length[0] - merge.mate_start[0];
        m_sta[1] = merged_lread - read_length[1] - merge.mate_start[1];
        m_sta.swap(0, 1);
    }
    let m_end = [m_sta[0] + m_len[0], m_sta[1] + m_len[1]];

    // Precompute, once, which of the ORIGINAL (pre-split) inter-exon gaps are splice junctions
    // and which (motif, annotated) pair each corresponds to -- copied from the merged transcript
    // rather than re-detected (STAR-rs's approach: reuse the classification the SE aligner
    // already made). Interior gaps are unaffected by mate-boundary clipping (only the first/last
    // kept exon of a mate is ever clipped, and clipping shifts read coordinates by a constant
    // without touching genome coordinates), so this classification carries over unchanged to the
    // post-split per-mate exon lists.
    let n = merged.exons.len();
    let mut gap_info: Vec<Option<(SpliceMotif, bool)>> = Vec::with_capacity(n.saturating_sub(1));
    let mut orig_junction_idx = 0usize;
    for iex in 0..n.saturating_sub(1) {
        let cur = &merged.exons[iex];
        let next = &merged.exons[iex + 1];
        let read_gap = next.read_start as i64 - cur.read_end as i64;
        let genome_gap = next.genome_start as i64 - cur.genome_end as i64;
        if read_gap == 0 && genome_gap >= scorer.align_intron_min as i64 {
            let entry = if orig_junction_idx < merged.junction_motifs.len() {
                Some((
                    merged.junction_motifs[orig_junction_idx],
                    merged.junction_annotated[orig_junction_idx],
                ))
            } else {
                None
            };
            gap_info.push(entry);
            orig_junction_idx += 1;
        } else {
            gap_info.push(None);
        }
    }

    // out[0] collects mate1's exons/junctions, out[1] collects mate2's (indexed by final PE
    // mate index i_frag, not by the `imate` scan order below -- see module doc for why these
    // differ when `merged.is_reverse`).
    let mut out_exons: [Vec<Exon>; 2] = [Vec::new(), Vec::new()];
    let mut out_junctions: [Vec<(SpliceMotif, bool)>; 2] = [Vec::new(), Vec::new()];

    for imate in 0..2usize {
        let i_frag = if imate == 0 { s } else { 1 - s };
        let mut prev_kept_iex: Option<usize> = None;
        for (iex, ex) in merged.exons.iter().enumerate() {
            let er = ex.read_start as u64;
            let el = (ex.read_end - ex.read_start) as u64;
            let in_span = er < m_end[imate] && er + el > m_sta[imate];
            if !in_span {
                prev_kept_iex = None;
                continue;
            }
            let (new_r, new_g, mut new_l) = if er >= m_sta[imate] {
                (er - m_sta[imate], ex.genome_start, el)
            } else {
                let delta = m_sta[imate] - er;
                (0u64, ex.genome_start + delta, el - delta)
            };
            if er + el > m_end[imate] {
                new_l -= er + el - m_end[imate];
            }
            if new_l == 0 {
                prev_kept_iex = None;
                continue;
            }
            if let Some(p) = prev_kept_iex
                && p + 1 == iex
                && let Some(j) = gap_info.get(p).copied().flatten()
            {
                out_junctions[i_frag].push(j);
            }
            out_exons[i_frag].push(Exon {
                genome_start: new_g,
                genome_end: new_g + new_l,
                read_start: new_r as usize,
                read_end: (new_r + new_l) as usize,
                i_frag: 0, // local to this mate's own Transcript; see module doc.
            });
            if out_exons[i_frag].len() > MAX_N_EXONS {
                return None;
            }
            prev_kept_iex = Some(iex);
        }
        if out_exons[i_frag].is_empty() {
            return None;
        }
    }

    // mate1 (output index 0) takes on the merged transcript's own strand; mate2 (output index 1)
    // is always the opposite strand (PE mates are always FR/RF -- see module doc derivation).
    let out_is_reverse = [merged.is_reverse, !merged.is_reverse];
    let mate_seqs: [&[u8]; 2] = [mate1_seq, mate2_seq];
    let mate_lens = [l1 as usize, l2 as usize];

    let mut mates: Vec<Transcript> = Vec::with_capacity(2);
    for i in 0..2usize {
        let native_read = if out_is_reverse[i] {
            reverse_complement(mate_seqs[i])
        } else {
            mate_seqs[i].to_vec()
        };
        let (cigar_ops, score, n_mismatch, n_gap, n_junction) = score_mate_exons(
            &out_exons[i],
            &out_junctions[i],
            mate_lens[i],
            &native_read,
            genome,
            scorer,
        );
        let genome_start = out_exons[i].iter().map(|e| e.genome_start).min().unwrap();
        let genome_end = out_exons[i].iter().map(|e| e.genome_end).max().unwrap();
        mates.push(Transcript {
            chr_idx: merged.chr_idx,
            genome_start,
            genome_end,
            is_reverse: out_is_reverse[i],
            exons: std::mem::take(&mut out_exons[i]),
            cigar: cigar_ops,
            score,
            n_mismatch,
            n_gap,
            n_junction,
            junction_motifs: out_junctions[i].iter().map(|(m, _)| *m).collect(),
            junction_annotated: out_junctions[i].iter().map(|(_, a)| *a).collect(),
        });
    }

    let mate2 = mates.pop().unwrap();
    let mate1 = mates.pop().unwrap();
    Some((mate1, mate2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::score::AlignmentScorer;
    use crate::genome::Genome;

    #[test]
    fn local_search_finds_the_overlap_offset() {
        // rc(mate2) = [2,3,0,1,2,3] starts at offset 2 within mate1 = [0,1,2,3,0,1].
        let m1 = [0u8, 1, 2, 3, 0, 1];
        let m2rc = [2u8, 3, 0, 1, 2, 3];
        assert_eq!(local_search_n_is_mm(&m1, &m2rc, 0.01), 2);
        // No qualifying offset -> returns x.len().
        let a = [0u8, 0, 0, 0];
        let b = [1u8, 1, 1, 1];
        assert_eq!(local_search_n_is_mm(&a, &b, 0.01), a.len());
        // An `N` (code 4) in y counts as a mismatch even against an equal x base.
        let x = [0u8, 1, 2, 3];
        let y = [0u8, 1, 4, 3]; // 3 matches, 1 mismatch (the N) at offset 0 -> 1/3 > 0.01, rejected
        assert_eq!(local_search_n_is_mm(&x, &y, 0.01), x.len());
        assert_eq!(local_search_n_is_mm(&x, &y, 0.5), 0); // 1/3 <= 0.5 -> accepted at offset 0
    }

    #[test]
    fn merge_reconstructs_the_fragment() {
        // Fragment ACGTACGT: mate1 = F[0..6], rc(mate2) = F[2..8], 4 bp overlap.
        let m1 = [0u8, 1, 2, 3, 0, 1];
        let m2rc = [2u8, 3, 0, 1, 2, 3];
        let m = pe_merge_mates(&m1, &m2rc, 4, 0.01).expect("should merge");
        assert_eq!(m.n_ov, 4);
        assert_eq!(m.mate_start, [0, 2]);
        assert_eq!(m.merged, vec![0, 1, 2, 3, 0, 1, 2, 3]); // == the original fragment
        // Below the min-overlap threshold -> no merge.
        assert!(pe_merge_mates(&m1, &m2rc, 5, 0.01).is_none());
    }

    #[test]
    fn merge_when_one_mate_is_contained() {
        // rc(mate2) fully contained at the end of mate1 (o1 == l2). Non-repetitive bases so the
        // overlap offset is unambiguous (STAR's strict `>` tie-break keeps the first-found offset).
        let m1 = [0u8, 1, 2, 3, 1, 0, 3, 2];
        let m2rc = [1u8, 0, 3, 2]; // == m1[4..8], unique within m1
        let m = pe_merge_mates(&m1, &m2rc, 3, 0.01).expect("should merge");
        assert_eq!(m.n_ov, 4);
        assert_eq!(m.mate_start, [0, 4]);
        assert_eq!(m.merged, m1.to_vec()); // merged == mate1 (mate2 adds nothing)
    }

    /// Build a tiny genome for conversion tests: one chromosome holding exactly `forward`.
    /// The reverse-complement half of `sequence` is left as padding since
    /// [`convert_merged_transcript_to_pe`]'s rescoring only ever reads forward coordinates
    /// (mate strand is handled by reverse-complementing the *read*, not the genome lookup).
    fn tiny_genome(forward: &[u8]) -> Genome {
        let n = forward.len() as u64;
        let mut sequence = vec![5u8; (n * 2) as usize];
        sequence[..forward.len()].copy_from_slice(forward);
        Genome {
            transform_blocks: None,
            sequence: sequence.into(),
            n_genome: n,
            n_genome_real: n,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![n],
            chr_start: vec![0, n],
        }
    }

    #[test]
    fn convert_ungapped_forward_overlap_splits_cleanly() {
        // Reuse the exact (non-repetitive, unambiguous-overlap) sequences from
        // `merge_reconstructs_the_fragment`: mate1 = [0,1,2,3,0,1] (genome[0..6]),
        // rc(mate2) = [2,3,0,1,2,3] (genome[2..8]), 4bp overlap, merged == genome[0..8].
        let mate1_seq = vec![0u8, 1, 2, 3, 0, 1];
        let rc_mate2 = vec![2u8, 3, 0, 1, 2, 3];
        let mate2_seq = reverse_complement(&rc_mate2); // as actually sequenced (mate2 is reverse)

        let mut genome_full = mate1_seq.clone();
        genome_full.extend_from_slice(&rc_mate2[4..]); // genome[0..8] = [0,1,2,3,0,1,2,3]
        let genome = tiny_genome(&genome_full);

        let merge = pe_merge_mates(&mate1_seq, &rc_mate2, 4, 0.01).expect("should merge");
        assert_eq!(merge.merged, genome_full);
        assert_eq!(merge.mate_start, [0, 2]);

        // The merged read aligns perfectly, forward strand, one ungapped exon.
        let merged_len = merge.merged.len();
        let merged_tr = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: merged_len as u64,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: merged_len as u64,
                read_start: 0,
                read_end: merged_len,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, merged_len)],
            score: merged_len as i32,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let scorer = AlignmentScorer::from_params_minimal();
        let (m1, m2) = convert_merged_transcript_to_pe(
            &merged_tr, &merge, &mate1_seq, &mate2_seq, &genome, &scorer,
        )
        .expect("should convert");

        // mate1: forward, covers genome[0..6], perfect match.
        assert!(!m1.is_reverse);
        assert_eq!(m1.genome_start, 0);
        assert_eq!(m1.genome_end, 6);
        assert_eq!(m1.n_mismatch, 0);
        assert_eq!(
            m1.exons
                .iter()
                .map(|e| e.read_end - e.read_start)
                .sum::<usize>(),
            6
        );

        // mate2: reverse, covers genome[2..8], perfect match.
        assert!(m2.is_reverse);
        assert_eq!(m2.genome_start, 2);
        assert_eq!(m2.genome_end, 8);
        assert_eq!(m2.n_mismatch, 0);
        assert_eq!(
            m2.exons
                .iter()
                .map(|e| e.read_end - e.read_start)
                .sum::<usize>(),
            6
        );
    }

    /// Reverse-strand merged transcript: exercises the `m_sta = Lread - readLen - m_sta`
    /// coordinate swap + the mate-strand flip (mate1 takes the fragment/merged strand,
    /// mate2 the opposite) — the path the review flagged as untested.
    #[test]
    fn convert_reverse_strand_flips_mates_and_swaps_coords() {
        let mate1_seq = vec![0u8, 1, 2, 3, 0, 1];
        let rc_mate2 = vec![2u8, 3, 0, 1, 2, 3];
        let mate2_seq = reverse_complement(&rc_mate2);

        let mut genome_full = mate1_seq.clone();
        genome_full.extend_from_slice(&rc_mate2[4..]); // genome[0..8]
        let genome = tiny_genome(&genome_full);

        let merge = pe_merge_mates(&mate1_seq, &rc_mate2, 4, 0.01).expect("should merge");
        let merged_len = merge.merged.len();

        // Same perfect single-exon alignment, but the merged read maps to the REVERSE
        // strand (is_reverse = true).
        let merged_tr = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: merged_len as u64,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: merged_len as u64,
                read_start: 0,
                read_end: merged_len,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, merged_len)],
            score: merged_len as i32,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let scorer = AlignmentScorer::from_params_minimal();
        let (m1, m2) = convert_merged_transcript_to_pe(
            &merged_tr, &merge, &mate1_seq, &mate2_seq, &genome, &scorer,
        )
        .expect("should convert");

        // Strand flip: mate1 takes the merged (reverse) strand, mate2 the opposite.
        assert!(m1.is_reverse, "mate1 should take the merged reverse strand");
        assert!(!m2.is_reverse, "mate2 should be the opposite (forward)");

        // Each mate still covers its full length with no mismatch, within genome [0,8].
        assert_eq!(m1.n_mismatch, 0);
        assert_eq!(m2.n_mismatch, 0);
        for m in [&m1, &m2] {
            let covered: usize = m.exons.iter().map(|e| e.read_end - e.read_start).sum();
            assert_eq!(covered, 6, "each mate covers its 6 bp");
            assert!(m.genome_start < m.genome_end && m.genome_end <= merged_len as u64);
            assert_eq!(m.genome_end - m.genome_start, 6, "6 bp genomic span");
        }
        // The two mates occupy the mirrored offsets of the forward case: [2,8) and [0,6).
        let spans: std::collections::HashSet<(u64, u64)> = [
            (m1.genome_start, m1.genome_end),
            (m2.genome_start, m2.genome_end),
        ]
        .into_iter()
        .collect();
        assert!(spans.contains(&(0, 6)) && spans.contains(&(2, 8)));
    }
}
