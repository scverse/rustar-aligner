//! `--clipAdapterType CellRanger4`: the 10x Chromium v4 clipping rules.
//!
//! Two independent trims, ported from STAR's `ClipCR4.cpp` and
//! `ClipMate_clipChunk.cpp`:
//!
//! - a 3' poly-A tail trim, scored base by base from the 3' end;
//! - a 5' template-switch-oligo trim, which is an overlap alignment of the TSO
//!   against the first 91 bases of the read.
//!
//! Only the poly-A trim is implemented here. The 5' TSO trim is an overlap
//! alignment, for which STAR links the Opal SIMD library; rustar will take
//! that from a dedicated deterministic-SIMD crate rather than carrying a
//! second aligner in-tree. Until then `--clipAdapterType CellRanger4` is
//! rejected rather than silently doing half the job.

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
