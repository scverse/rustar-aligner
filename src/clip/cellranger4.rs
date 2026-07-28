//! `--clipAdapterType CellRanger4`: the 10x Chromium v4 clipping rules.
//!
//! Two independent trims, ported from STAR's `ClipCR4.cpp` and
//! `ClipMate_clipChunk.cpp`:
//!
//! - a 3' poly-A tail trim, scored base by base from the 3' end;
//! - a 5' template-switch-oligo trim, which is an overlap alignment of the TSO
//!   against the first 91 bases of the read.
//!
//! STAR does the 5' alignment with the Opal SIMD library. Here it goes through
//! [`crate::swalign`], which is required to be bit-identical across
//! instruction sets, so the clip length cannot depend on the machine.

use crate::swalign::{self, Mode, Scoring};

/// Number of 3' bases to trim as a CellRanger4 poly-A tail.
///
/// STAR `ClipCR4::polyTail3p`. Walks in from the 3' end scoring `+1` per `A`
/// and `-2` per non-`A`, and remembers the longest prefix of that walk whose
/// running score still clears a 70% density threshold (`score * 10 >= ib * 7`).
/// It gives up once the score has fallen more than 27 behind the position, and
/// returns nothing unless the remembered score reached 20.
///
/// `seq` is numeric base codes, so `A == 0`.
pub fn poly_tail_3p(seq: &[u8]) -> usize {
    let seq_len = seq.len();
    if seq_len < 20 {
        return 0;
    }
    let mut best_len: i64 = seq_len as i64 - 1;
    let mut score: i64 = 0;
    let mut best_score: i64 = 0;
    for ib in 1..=seq_len as i64 {
        if seq[seq_len - ib as usize] == 0 {
            score += 1;
            if score * 10 >= ib * 7 {
                best_len = ib;
                best_score = score;
            }
        } else {
            score -= 2;
            if ib - score > 27 {
                break;
            }
        }
    }
    if best_score < 20 {
        0
    } else {
        best_len as usize
    }
}

/// How much of the read STAR aligns the TSO against (`ClipCR4::opalFillOneSeq`).
const CR4_TARGET_LEN: usize = 91;

/// Number of 5' bases to trim as the 10x TSO.
///
/// STAR aligns the TSO against the first 91 bases of the read in overlap mode,
/// asking for the score and the position in the target where it ends, then
/// applies an acceptance gate: a score below 20 is rejected outright, and
/// scores of exactly 20 or 21 are rejected if they took too long to reach
/// (more than 26 and 30 bases respectively). A weak alignment that happens to
/// run a long way is what that gate is there to catch.
///
/// The read is padded to 91 bases with `N` when it is shorter, which is why
/// the scoring scheme has to treat `N` against `N` as neutral rather than as a
/// mismatch: otherwise the padding would drag every score down.
///
/// Both arguments are numeric base codes; the return value is a count of 5'
/// bases to clip, `0` when the alignment is rejected.
pub fn tso_clip(read: &[u8], tso: &[u8]) -> usize {
    if tso.is_empty() {
        return 0;
    }
    let take = read.len().min(CR4_TARGET_LEN);
    let mut target = Vec::with_capacity(CR4_TARGET_LEN);
    target.extend_from_slice(&read[..take]);
    target.resize(CR4_TARGET_LEN, 4); // N padding

    let a = swalign::align(tso, &target, Mode::Ov, &Scoring::CLIP_CR4);
    let clip = a.target_end as i64 + 1; // 1-based end == number of bases covered

    let reject = a.score < 20 || (a.score == 20 && clip > 26) || (a.score == 21 && clip > 30);
    if reject { 0 } else { clip as usize }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ACGT text to base codes.
    fn code(s: &str) -> Vec<u8> {
        s.bytes()
            .map(|b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect()
    }

    /// The 10x template switch oligo.
    const TSO: &str = "AAGCAGTGGTATCAACGCAGAGTACATGGG";

    #[test]
    fn cr4_tso_clip_matches_opal() {
        // Frozen vector shared with STAR-rs's `cr4_tso_clip_matches_opal`,
        // which validates the same numbers against Opal itself. A read that
        // starts with the TSO is clipped by exactly its length.
        let read = code(&format!("{TSO}ACGTACGTACGTACGTACGTACGTACGTAC"));
        assert_eq!(tso_clip(&read, &code(TSO)), 30);

        // A read with no TSO is not clipped at all.
        let read = code("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        assert_eq!(tso_clip(&read, &code(TSO)), 0);
    }

    #[test]
    fn tso_clip_is_inert_without_an_adapter() {
        let read = code(&format!("{TSO}ACGTACGT"));
        assert_eq!(tso_clip(&read, &[]), 0);
    }

    #[test]
    fn cr4_polya_trim_matches_star() {
        // A clean 30-base poly-A tail is trimmed whole. The prefix is
        // deliberately A-free: the scan does not stop at the tail boundary, so
        // an `A` just upstream of it would legitimately extend the trim.
        let read = code(&format!("CGTCGTCGTCGTCGTCGTCG{}", "A".repeat(30)));
        assert_eq!(poly_tail_3p(&read), 30);

        // No tail: nothing to trim.
        let read = code("ACGTACGTACGTACGTACGTACGTACGTACGT");
        assert_eq!(poly_tail_3p(&read), 0);

        // A tail shorter than the score-20 floor is not trimmed, however clean.
        let read = code(&format!("CGTCGTCGTCGTCGTCGTCG{}", "A".repeat(10)));
        assert_eq!(poly_tail_3p(&read), 0);
    }

    #[test]
    fn poly_tail_needs_twenty_bases_of_read() {
        assert_eq!(poly_tail_3p(&code("AAAAAAAAAAAAAAAAAAA")), 0); // 19
    }

    #[test]
    fn poly_tail_scan_does_not_stop_at_the_tail_boundary() {
        // STAR keeps scoring past the run of A's, so A-rich sequence just
        // upstream extends the trim. Worth pinning: it looks like an off-by-one
        // otherwise.
        let read = code(&format!("ACGTACGTACGTACGTACGT{}", "A".repeat(30)));
        assert_eq!(poly_tail_3p(&read), 34);
    }

    #[test]
    fn poly_tail_tolerates_a_single_interruption() {
        // 70% density is the threshold, so one non-A inside a long tail is
        // survivable.
        let tail = format!("{}C{}", "A".repeat(15), "A".repeat(15));
        let read = code(&format!("CGTCGTCGTCGTCGTCGTCG{tail}"));
        assert!(
            poly_tail_3p(&read) >= 15,
            "one mismatch should not abandon the tail, got {}",
            poly_tail_3p(&read)
        );
    }
}
