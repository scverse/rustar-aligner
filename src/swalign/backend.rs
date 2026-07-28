//! Backend selection, and the differential harness that keeps backends honest.
//!
//! The contract every backend signs: given the same inputs it returns exactly
//! what [`super::scalar::align`] returns. Not "within rounding", not "the same
//! score with a different end position" — the same [`Alignment`].
//!
//! [`differential_check`] is what enforces that. It is written here rather than
//! in a test module so a backend can be checked from a test, a benchmark or a
//! debug session without duplicating the generator.

use super::{Alignment, Mode, Scoring, scalar};

/// Which implementation actually ran.
///
/// Exposed so tests can assert that a machine with the hardware really used it,
/// rather than silently falling back and reporting a pass that proves nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    /// The portable reference. Always available, and the definition of the
    /// result.
    Scalar,
    /// aarch64 NEON, computing one anti-diagonal at a time.
    #[cfg(target_arch = "aarch64")]
    Neon,
}

impl Backend {
    /// The best backend this machine can run.
    pub fn detect() -> Self {
        #[cfg(target_arch = "aarch64")]
        if super::neon::is_available() {
            return Backend::Neon;
        }
        Backend::Scalar
    }

    /// Run this backend.
    pub fn align(self, query: &[u8], target: &[u8], mode: Mode, scoring: &Scoring) -> Alignment {
        match self {
            Backend::Scalar => scalar::align(query, target, mode, scoring),
            #[cfg(target_arch = "aarch64")]
            Backend::Neon => super::neon::align(query, target, mode, scoring),
        }
    }
}

/// A deterministic sequence generator for the differential harness.
///
/// Seeded from the in-tree splitmix64 so the cases are identical on every
/// machine and every run: a backend that fails does so reproducibly, on a case
/// the report can name.
struct SeqGen(crate::rng::SplitMix64);

impl SeqGen {
    fn new(seed: u64) -> Self {
        Self(crate::rng::SplitMix64::seed(seed))
    }

    /// A sequence of `len` base codes. `n_rate` out of 16 bases are `N`, so the
    /// generator covers the `N`-heavy inputs STAR's padding produces as well as
    /// clean ones.
    fn seq(&mut self, len: usize, n_rate: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                let r = self.0.next_u64();
                if r % 16 < n_rate {
                    4
                } else {
                    (r >> 8) as u8 % 4
                }
            })
            .collect()
    }
}

/// Compare `backend` against the scalar reference over a spread of inputs.
///
/// Returns `Err` with a description of the first disagreement, naming the case
/// so it can be reproduced. `Ok(n)` reports how many cases were checked.
///
/// The spread is deliberately awkward: empty and length-1 sequences, queries
/// longer than targets, all-`N` inputs, and long runs that push scores far
/// enough to saturate a narrow lane. Saturation is where a SIMD backend is most
/// likely to diverge, so it is not left to chance.
pub fn differential_check(backend: Backend, seed: u64) -> Result<usize, String> {
    let scoring = Scoring::CLIP_CR4;
    let mut rng = SeqGen::new(seed);
    let mut checked = 0usize;

    // Lengths chosen around the boundaries a vectorised kernel cares about:
    // zero, one, just under and just over a typical lane count, and the 91
    // STAR actually uses.
    const LENS: &[usize] = &[0, 1, 2, 7, 8, 9, 15, 16, 17, 30, 31, 33, 64, 91, 128];

    for &ql in LENS {
        for &tl in LENS {
            for &n_rate in &[0u64, 1, 8, 16] {
                for mode in [Mode::Nw, Mode::Hw, Mode::Ov, Mode::Sw] {
                    let q = rng.seq(ql, n_rate);
                    let t = rng.seq(tl, n_rate);
                    let want = scalar::align(&q, &t, mode, &scoring);
                    let got = backend.align(&q, &t, mode, &scoring);
                    if got != want {
                        return Err(format!(
                            "{backend:?} disagrees with scalar on {mode:?}, \
                             |q|={ql} |t|={tl} n_rate={n_rate}/16: \
                             scalar {want:?}, backend {got:?}\n  query  {q:?}\n  target {t:?}"
                        ));
                    }
                    checked += 1;
                }
            }
        }
    }

    // Saturation: a long exact match scores one per base, which overflows an
    // 8-bit lane well before it overflows the scalar's i32. A backend that
    // escalates lane width incorrectly fails here and nowhere else.
    for &len in &[200usize, 400, 1000] {
        let q = rng.seq(len, 0);
        for mode in [Mode::Nw, Mode::Hw, Mode::Ov, Mode::Sw] {
            let want = scalar::align(&q, &q, mode, &scoring);
            let got = backend.align(&q, &q, mode, &scoring);
            if got != want {
                return Err(format!(
                    "{backend:?} disagrees with scalar under saturation, \
                     {mode:?}, len={len}: scalar {want:?}, backend {got:?}"
                ));
            }
            checked += 1;
        }
    }

    Ok(checked)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every backend this machine can run, not just the one `detect()` picks.
    ///
    /// On x86 that matters: a machine with AVX2 would otherwise never exercise
    /// SSE2, and the baseline is what runs on everything older.
    fn available_backends() -> Vec<Backend> {
        let mut v = vec![Backend::Scalar];
        #[cfg(target_arch = "aarch64")]
        if super::super::neon::is_available() {
            v.push(Backend::Neon);
        }
        v
    }

    #[test]
    fn every_available_backend_agrees_with_scalar() {
        // Not just the detected one: on a machine with AVX2, SSE2 would
        // otherwise go unchecked despite being what older hardware runs.
        for backend in available_backends() {
            if let Err(e) = differential_check(backend, 0x0BAD_5EED) {
                panic!("{e}");
            }
        }
    }

    #[test]
    fn scalar_is_consistent_with_itself() {
        // Tautological for the scalar backend, but it proves the harness runs,
        // covers every mode and length pair, and reaches the saturation cases.
        // A backend added later inherits a harness that is known to work.
        let n = differential_check(Backend::Scalar, 0x5EED).expect("scalar vs scalar");
        assert!(n > 900, "harness covered only {n} cases");
    }

    #[test]
    fn the_detected_backend_agrees_with_scalar_on_this_machine() {
        // The test that matters once SIMD backends exist: whatever this
        // machine selects must match the reference on this machine.
        let backend = Backend::detect();

        // A pass proves nothing if the machine quietly fell back to the
        // reference, so assert that the hardware backend really was selected.
        #[cfg(target_arch = "aarch64")]
        assert_eq!(
            backend,
            Backend::Neon,
            "aarch64 must select NEON; a silent fallback would make this test vacuous"
        );

        if let Err(e) = differential_check(backend, 0x00C0_FFEE) {
            panic!("{e}");
        }
    }

    #[test]
    fn generator_is_reproducible() {
        // The harness is only useful if a failure can be reproduced, which
        // needs the sequences to be identical run to run.
        let a = SeqGen::new(7).seq(64, 4);
        let b = SeqGen::new(7).seq(64, 4);
        assert_eq!(a, b);
        assert!(a.contains(&4), "n_rate 4/16 should produce Ns");
        assert!(a.iter().any(|&b| b < 4), "and also real bases");
    }
}
