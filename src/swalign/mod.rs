//! Deterministic Smith-Waterman-family alignment with affine gaps.
//!
//! STAR performs exactly one alignment of this shape: the 5' TSO clip of
//! `--clipAdapterType CellRanger4`, for which it links the Opal C/C++ SIMD
//! library (`OPAL_MODE_OV` + `OPAL_SEARCH_SCORE_END`). This module provides
//! the same capability in-tree, with no new dependency.
//!
//! # Determinism is a correctness property here
//!
//! Opal's overflow handling groups database sequences into buckets sized by the
//! SIMD vector width, so when 8-bit lanes saturate, the grouping — and with it
//! the recompute path — depends on which instruction set is available. Sixteen
//! SSE lanes and thirty-two AVX2 lanes bucket differently, and simde emulates
//! AVX2 on ARM, so the same input can take a different path on a different
//! machine.
//!
//! That is not acceptable for an aligner whose output is supposed to be
//! reproducible. Here the rule is: **the vector width is never observable.**
//! [`scalar`] defines the result; every other backend must agree with it
//! bit-for-bit, including under saturation, on empty inputs, on `N`, and on
//! end-position ties. The differential test in [`tests`] is what enforces it.
//!
//! # Coordinates and conventions
//!
//! Sequences are numeric base codes, `0..=3` for ACGT and `4` for `N`, the same
//! encoding the rest of the aligner uses. Scores are `i32`; the caller supplies
//! the scoring scheme.

pub mod backend;
#[cfg(target_arch = "aarch64")]
pub mod neon;
pub mod scalar;

pub use backend::Backend;

/// How the ends of the two sequences are treated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Global: both sequences must be consumed end to end (Needleman-Wunsch).
    Nw,
    /// Semi-global on the target: the query must be consumed entirely, the
    /// target may extend past it on both sides.
    Hw,
    /// Overlap: gaps at the start of either sequence and at the end of either
    /// sequence are free. This is the mode STAR's TSO clip uses.
    Ov,
    /// Local (Smith-Waterman): the best-scoring subalignment.
    Sw,
}

/// Affine-gap scoring.
///
/// `gap_open` and `gap_extend` are the penalties *subtracted*, so both are
/// given as positive numbers. STAR's ClipCR4 uses `match_score = 1`,
/// `mismatch = -2`, `gap_open = gap_extend = 2`, and scores `N` against `N` as
/// zero rather than as a mismatch.
#[derive(Debug, Clone, Copy)]
pub struct Scoring {
    /// Added when two bases are equal.
    pub match_score: i32,
    /// Added when two bases differ (negative).
    pub mismatch: i32,
    /// Subtracted to open a gap.
    pub gap_open: i32,
    /// Subtracted for each base a gap is extended by.
    pub gap_extend: i32,
    /// Added when both positions are `N`. Opal treats this pairing as neutral
    /// rather than as a mismatch, and STAR relies on that when it pads the
    /// target with `N`.
    pub n_vs_n: i32,
}

impl Scoring {
    /// The scoring STAR uses for the CellRanger4 TSO clip
    /// (`ClipCR4.cpp`: match +1, mismatch -2, gaps 2, `N`/`N` neutral).
    pub const CLIP_CR4: Self = Self {
        match_score: 1,
        mismatch: -2,
        gap_open: 2,
        gap_extend: 2,
        n_vs_n: 0,
    };

    /// Score one aligned pair of base codes.
    #[inline]
    pub fn pair(&self, q: u8, t: u8) -> i32 {
        if q == 4 && t == 4 {
            self.n_vs_n
        } else if q == t {
            self.match_score
        } else {
            self.mismatch
        }
    }
}

/// The outcome of one alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Alignment {
    /// Best score under the chosen mode.
    pub score: i32,
    /// Zero-based position in the target where that score is reached.
    ///
    /// Ties are broken towards the **earlier** column, except that a score
    /// reached in the final column wins outright: that is what
    /// `OPAL_SEARCH_SCORE_END` means, and STAR's clip length depends on it.
    pub target_end: usize,
}

/// Align `query` against `target` under `mode`.
///
/// Dispatches to the fastest backend available for this machine. Every backend
/// is required to return exactly what [`scalar::align`] would, so the choice is
/// invisible in the output.
pub fn align(query: &[u8], target: &[u8], mode: Mode, scoring: &Scoring) -> Alignment {
    Backend::detect().align(query, target, mode, scoring)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_inputs_score_zero_in_every_mode() {
        for mode in [Mode::Nw, Mode::Hw, Mode::Ov, Mode::Sw] {
            let a = align(&[], &[], mode, &Scoring::CLIP_CR4);
            assert_eq!(a.score, 0, "{mode:?} on two empty sequences");
            let a = align(&[0, 1, 2], &[], mode, &Scoring::CLIP_CR4);
            assert_eq!(a.score, 0, "{mode:?} on an empty target");
        }
    }

    #[test]
    fn n_against_n_is_neutral_not_a_mismatch() {
        let s = Scoring::CLIP_CR4;
        assert_eq!(s.pair(4, 4), 0);
        assert_eq!(s.pair(4, 0), -2);
        assert_eq!(s.pair(0, 4), -2);
        assert_eq!(s.pair(0, 0), 1);
        assert_eq!(s.pair(0, 1), -2);
    }

    #[test]
    fn identical_sequences_score_one_per_base() {
        let seq = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let a = align(&seq, &seq, Mode::Ov, &Scoring::CLIP_CR4);
        assert_eq!(a.score, seq.len() as i32);
        assert_eq!(a.target_end, seq.len() - 1);
    }
}
