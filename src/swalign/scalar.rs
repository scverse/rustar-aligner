//! Portable reference implementation.
//!
//! This is the definition of the result. It is written for clarity rather than
//! speed: the SIMD backends must agree with it bit-for-bit, so it needs to be
//! obviously correct more than it needs to be fast.

use super::{Alignment, Mode, Scoring};

/// Sentinel for "unreachable". Far below any real score, and far enough from
/// `i32::MIN` that subtracting a gap penalty cannot wrap.
const NEG: i32 = i32::MIN / 4;

/// Align `query` against `target`, filling the DP column by column.
///
/// Only two columns are ever live, so the working set is `O(|query|)` rather
/// than the full matrix. That is also what makes the striped SIMD form a
/// drop-in replacement later: it computes the same columns in the same order.
pub fn align(query: &[u8], target: &[u8], mode: Mode, scoring: &Scoring) -> Alignment {
    let ql = query.len();
    let tl = target.len();
    if ql == 0 || tl == 0 {
        return Alignment {
            score: 0,
            target_end: 0,
        };
    }

    // Whether a gap before the start of each sequence is free.
    let (free_query_start, free_target_start) = match mode {
        Mode::Nw => (false, false),
        Mode::Hw => (false, true),
        Mode::Ov | Mode::Sw => (true, true),
    };

    // Column 0 of the previous iteration: H is the score of aligning the first
    // `r+1` query bases against nothing.
    let mut prev_h = vec![0i32; ql];
    let mut prev_e = vec![NEG; ql];
    if !free_query_start {
        for (r, h) in prev_h.iter_mut().enumerate() {
            *h = -(scoring.gap_open + scoring.gap_extend * r as i32);
        }
    }

    // Best score seen on the final query row, and the column that first
    // achieved it. Ties go to the earlier column.
    let mut best_last_row = NEG;
    let mut best_last_row_col = 0usize;
    // Best score anywhere in the most recent column.
    let mut last_col_max = NEG;
    // Best score anywhere in the matrix, for local mode.
    let mut best_anywhere = 0i32;
    let mut best_anywhere_col = 0usize;

    for (c, &tc) in target.iter().enumerate() {
        // Top of the column: aligning the first `c+1` target bases against
        // nothing.
        let mut up_h = if free_target_start {
            0
        } else {
            -(scoring.gap_open + scoring.gap_extend * c as i32)
        };
        // The diagonal predecessor, i.e. the cell up and to the left. Column 0
        // has no predecessor and scores 0 either way, so it falls out of the
        // free-start case.
        let mut diag_h = if free_target_start || c == 0 {
            0
        } else {
            -(scoring.gap_open + scoring.gap_extend * (c as i32 - 1))
        };
        let mut up_f = NEG;
        let mut col_max = NEG;
        let mut h = NEG;

        for r in 0..ql {
            // Gap in the query (moving right): open from H, or extend E.
            let e = (prev_h[r] - scoring.gap_open).max(prev_e[r] - scoring.gap_extend);
            // Gap in the target (moving down): open from H, or extend F.
            let f = (up_h - scoring.gap_open).max(up_f - scoring.gap_extend);
            // Match or mismatch on the diagonal.
            let d = diag_h + scoring.pair(query[r], tc);

            h = e.max(f).max(d);
            if mode == Mode::Sw {
                // Local alignment never carries a negative prefix forward.
                h = h.max(0);
                if h > best_anywhere {
                    best_anywhere = h;
                    best_anywhere_col = c;
                }
            }
            if h > col_max {
                col_max = h;
            }

            up_f = f;
            up_h = h;
            diag_h = prev_h[r];
            prev_e[r] = e;
            prev_h[r] = h;
        }

        // `h` now holds the last query row for this column.
        if h > best_last_row {
            best_last_row = h;
            best_last_row_col = c;
        }
        last_col_max = col_max;
    }

    match mode {
        Mode::Sw => Alignment {
            score: best_anywhere,
            target_end: best_anywhere_col,
        },
        Mode::Nw => Alignment {
            // Global: the corner cell, which is the last row of the last
            // column.
            score: prev_h[ql - 1],
            target_end: tl - 1,
        },
        Mode::Hw => Alignment {
            // Semi-global on the target: the query must be consumed, so the
            // answer lives on the last query row.
            score: best_last_row,
            target_end: best_last_row_col,
        },
        Mode::Ov => {
            // Overlap: either the query ran out (best on the last row) or the
            // target ran out (best in the last column). A score reached in the
            // final column wins the tie, which is what `OPAL_SEARCH_SCORE_END`
            // specifies and what STAR's clip length depends on.
            let score = last_col_max.max(best_last_row);
            let target_end = if last_col_max >= best_last_row {
                tl - 1
            } else {
                best_last_row_col
            };
            Alignment { score, target_end }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::swalign::Scoring;

    const S: Scoring = Scoring::CLIP_CR4;

    #[test]
    fn global_mode_pays_for_every_unaligned_base() {
        // Query is a prefix of the target; NW must pay to skip the tail.
        let q = [0u8, 1, 2];
        let t = [0u8, 1, 2, 3, 3];
        let nw = align(&q, &t, Mode::Nw, &S);
        let ov = align(&q, &t, Mode::Ov, &S);
        assert_eq!(ov.score, 3, "overlap: the tail is free");
        assert!(
            nw.score < ov.score,
            "global should pay for the tail, got {} vs {}",
            nw.score,
            ov.score
        );
    }

    #[test]
    fn local_mode_ignores_flanking_mismatch() {
        // A clean 4-base core buried in mismatching flanks.
        let q = [3u8, 3, 0, 1, 2, 3, 3, 3];
        let t = [1u8, 1, 0, 1, 2, 3, 1, 1];
        let sw = align(&q, &t, Mode::Sw, &S);
        assert!(
            sw.score >= 4,
            "local should find the shared core, got {}",
            sw.score
        );
    }

    #[test]
    fn local_score_is_never_negative() {
        let q = [0u8, 0, 0, 0];
        let t = [3u8, 3, 3, 3];
        assert_eq!(align(&q, &t, Mode::Sw, &S).score, 0);
    }

    #[test]
    fn overlap_prefers_the_final_column_on_a_tie() {
        // A query that matches equally well at two positions: the one that
        // runs to the end of the target must win, because that is the tie-break
        // STAR's clip length is built on.
        let q = [0u8, 1];
        let t = [0u8, 1, 4, 4, 0, 1];
        let a = align(&q, &t, Mode::Ov, &S);
        assert_eq!(a.target_end, t.len() - 1);
    }

    #[test]
    fn affine_gaps_cost_less_than_repeated_opens() {
        // One 3-base gap should beat three separate 1-base gaps.
        let one_long = S.gap_open + S.gap_extend * 3;
        let three_short = 3 * (S.gap_open + S.gap_extend);
        assert!(one_long < three_short);
    }
}
