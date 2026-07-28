//! aarch64 NEON backend.
//!
//! # Why anti-diagonals rather than Farrar stripes
//!
//! The striped layout is faster, but its lazy-F correction and its
//! stripe-relative indexing make the *end position* awkward to extract, and the
//! end position is exactly what STAR's clip length is built on. A backend that
//! is quicker but disagrees with [`super::scalar`] about where an alignment
//! ends is worthless here.
//!
//! Cells on one anti-diagonal `d = r + c` are mutually independent: `H(r,c)`
//! depends on `E(r,c-1)` and `F(r-1,c)`, both on `d-1`, and on `H(r-1,c-1)` on
//! `d-2`. So a whole anti-diagonal can be computed in parallel with no
//! correction pass and no reordering, which makes agreement with the scalar
//! path a property of the layout rather than something to be tested for and
//! hoped about.
//!
//! The cost is strided access. For the 30×91 matrix STAR's TSO clip actually
//! uses, that is irrelevant.

use std::arch::aarch64::{vaddq_s32, vdupq_n_s32, vld1q_s32, vmaxq_s32, vst1q_s32, vsubq_s32};

use super::{Alignment, Mode, Scoring};

/// Matches `scalar::NEG`: far below any real score, far enough from `i32::MIN`
/// that subtracting a gap penalty cannot wrap.
const NEG: i32 = i32::MIN / 4;

/// Lanes per NEON vector at 32-bit width.
const LANES: usize = 4;

/// Align `query` against `target` using NEON.
///
/// # Safety
///
/// The caller must have established that NEON is available. On aarch64 it is
/// architecturally guaranteed, so [`is_available`] is a constant.
pub fn align(query: &[u8], target: &[u8], mode: Mode, scoring: &Scoring) -> Alignment {
    let ql = query.len();
    let tl = target.len();
    if ql == 0 || tl == 0 {
        return Alignment {
            score: 0,
            target_end: 0,
        };
    }

    let (free_query_start, free_target_start) = match mode {
        Mode::Nw => (false, false),
        Mode::Hw => (false, true),
        Mode::Ov | Mode::Sw => (true, true),
    };

    // Score of aligning the first `r+1` query bases against nothing, i.e. the
    // cell one column to the left of column 0.
    let left_h = |r: usize| -> i32 {
        if free_query_start {
            0
        } else {
            -(scoring.gap_open + scoring.gap_extend * r as i32)
        }
    };
    // Score of aligning the first `c+1` target bases against nothing, i.e. the
    // cell one row above row 0.
    let top_h = |c: i64| -> i32 {
        // A column before the first scores 0 either way, so it folds into the
        // free-start case.
        if free_target_start || c < 0 {
            0
        } else {
            -(scoring.gap_open + scoring.gap_extend * c as i32)
        }
    };

    // Rows of the two previous anti-diagonals, indexed by `r`.
    let mut h1 = vec![NEG; ql]; // H on d-1
    let mut e1 = vec![NEG; ql]; // E on d-1
    let mut f1 = vec![NEG; ql]; // F on d-1
    let mut h2 = vec![NEG; ql]; // H on d-2
    let mut h0 = vec![NEG; ql];
    let mut e0 = vec![NEG; ql];
    let mut f0 = vec![NEG; ql];

    let mut best_last_row = NEG;
    let mut best_last_row_col = 0usize;
    let mut last_col_max = NEG;
    let mut best_anywhere = 0i32;
    let mut best_anywhere_col = 0usize;
    // The scalar sweeps column-major and keeps the first cell to reach a new
    // maximum, so among equal scores the smallest `(c, r)` wins. Visiting
    // anti-diagonals changes the order, so the key has to be compared
    // explicitly rather than inferred from arrival.
    let mut best_anywhere_row = 0usize;

    let go = scoring.gap_open;
    let ge = scoring.gap_extend;

    for d in 0..(ql + tl - 1) {
        let r_lo = d.saturating_sub(tl - 1);
        let r_hi = d.min(ql - 1);

        // SAFETY: every load and store below is a plain lane-wise op on stack
        // scalars gathered by index; no pointer arithmetic escapes the slices,
        // which are all length `ql` and indexed within `r_lo..=r_hi`.
        unsafe {
            let vgo = vdupq_n_s32(go);
            let vge = vdupq_n_s32(ge);

            let mut row = r_lo;
            while row <= r_hi {
                let lanes = LANES.min(r_hi - row + 1);

                // Gather the four predecessor terms for lanes r..r+n.
                let mut e_prev_h = [NEG; LANES]; // H(r, c-1)
                let mut e_prev_e = [NEG; LANES]; // E(r, c-1)
                let mut f_prev_h = [NEG; LANES]; // H(r-1, c)
                let mut f_prev_f = [NEG; LANES]; // F(r-1, c)
                let mut diag = [NEG; LANES]; // H(r-1, c-1)
                let mut sub = [0i32; LANES]; // score(q[r], t[c])

                for k in 0..lanes {
                    let rr = row + k;
                    let col = d - rr;
                    // (rr, c-1) lives on d-1 at row rr.
                    if col == 0 {
                        e_prev_h[k] = left_h(rr);
                        e_prev_e[k] = NEG;
                    } else {
                        e_prev_h[k] = h1[rr];
                        e_prev_e[k] = e1[rr];
                    }
                    // (rr-1, c) lives on d-1 at row rr-1.
                    if rr == 0 {
                        f_prev_h[k] = top_h(col as i64);
                        f_prev_f[k] = NEG;
                    } else {
                        f_prev_h[k] = h1[rr - 1];
                        f_prev_f[k] = f1[rr - 1];
                    }
                    // (rr-1, c-1) lives on d-2 at row rr-1.
                    diag[k] = if rr == 0 {
                        top_h(col as i64 - 1)
                    } else if col == 0 {
                        left_h(rr - 1)
                    } else {
                        h2[rr - 1]
                    };
                    sub[k] = scoring.pair(query[rr], target[col]);
                }

                let ins = vmaxq_s32(
                    vsubq_s32(vld1q_s32(e_prev_h.as_ptr()), vgo),
                    vsubq_s32(vld1q_s32(e_prev_e.as_ptr()), vge),
                );
                let del = vmaxq_s32(
                    vsubq_s32(vld1q_s32(f_prev_h.as_ptr()), vgo),
                    vsubq_s32(vld1q_s32(f_prev_f.as_ptr()), vge),
                );
                let diagv = vaddq_s32(vld1q_s32(diag.as_ptr()), vld1q_s32(sub.as_ptr()));
                let mut cur = vmaxq_s32(vmaxq_s32(ins, del), diagv);
                if mode == Mode::Sw {
                    cur = vmaxq_s32(cur, vdupq_n_s32(0));
                }

                let mut hs = [0i32; LANES];
                let mut es = [0i32; LANES];
                let mut fs = [0i32; LANES];
                vst1q_s32(hs.as_mut_ptr(), cur);
                vst1q_s32(es.as_mut_ptr(), ins);
                vst1q_s32(fs.as_mut_ptr(), del);

                for k in 0..lanes {
                    let rr = row + k;
                    h0[rr] = hs[k];
                    e0[rr] = es[k];
                    f0[rr] = fs[k];
                }
                row += lanes;
            }
        }

        // Bookkeeping, in increasing column order so ties go to the earlier
        // column exactly as the scalar's column loop does.
        for (rr, &cell) in h0.iter().enumerate().take(r_hi + 1).skip(r_lo) {
            let col = d - rr;
            if mode == Mode::Sw
                && (cell > best_anywhere
                    || (cell == best_anywhere
                        && (col, rr) < (best_anywhere_col, best_anywhere_row)))
            {
                best_anywhere = cell;
                best_anywhere_col = col;
                best_anywhere_row = rr;
            }
            if rr == ql - 1 && cell > best_last_row {
                best_last_row = cell;
                best_last_row_col = col;
            }
            if col == tl - 1 {
                last_col_max = last_col_max.max(cell);
            }
        }

        std::mem::swap(&mut h2, &mut h1);
        std::mem::swap(&mut h1, &mut h0);
        std::mem::swap(&mut e1, &mut e0);
        std::mem::swap(&mut f1, &mut f0);
    }

    match mode {
        Mode::Sw => Alignment {
            score: best_anywhere,
            target_end: best_anywhere_col,
        },
        Mode::Nw => Alignment {
            score: h1[ql - 1],
            target_end: tl - 1,
        },
        Mode::Hw => Alignment {
            score: best_last_row,
            target_end: best_last_row_col,
        },
        Mode::Ov => {
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

/// NEON is architecturally guaranteed on aarch64.
pub const fn is_available() -> bool {
    true
}
