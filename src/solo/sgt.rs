//! Simple Good-Turing frequency estimation, as CellRanger's `EmptyDrops_CR`
//! and STAR's port of it use for the ambient profile.
//!
//! Gadsby & Sampson's method by way of David Elworthy's C implementation
//! (`SimpleGoodTuring/sgt.h`), which is what STAR vendors.
//!
//! The problem it solves: the ambient profile is estimated from the counts in
//! empty droplets, and a gene seen zero times there is not a gene with zero
//! probability — it is a gene whose probability the sample was too small to
//! show. Good-Turing reserves mass for those unseen events from the frequency
//! of the events seen exactly once, and smooths the rest along a fitted
//! log-log line so that sparsely-observed counts do not inherit the noise of
//! their raw frequencies.
//!
//! # D17
//!
//! STAR leaves `PZero` uninitialised until `analyse()` runs, and `analyse()`
//! returns early without setting it when there are fewer than five distinct
//! frequencies. A caller that then asks for the probability of an unseen gene
//! reads whatever was on the stack. Here it is 0.0 from construction: with too
//! few distinct frequencies there is no basis for reserving unseen mass, and
//! zero is the answer that says so. Any input large enough to reach the
//! significance test has more than five, so this is a divergence on a path
//! STAR's own output is undefined on.

use std::collections::BTreeMap;

/// One observed frequency and its smoothed probability.
#[derive(Debug, Clone, Copy, Default)]
struct Entry {
    /// How many events were observed this many times.
    freq: u32,
    /// The smoothed probability, filled in by [`Sgt::analyse`].
    estimate: f64,
}

/// A Simple Good-Turing estimator over a frequency-of-frequencies table.
#[derive(Debug, Default)]
pub struct Sgt {
    data: BTreeMap<u32, Entry>,
    /// Total probability reserved for events never observed.
    pzero: f64,
}

fn sq(x: f64) -> f64 {
    x * x
}

impl Sgt {
    pub fn new() -> Self {
        Self::default()
    }

    /// Record that `frequency` distinct events were each seen `observation`
    /// times.
    pub fn add(&mut self, observation: u32, frequency: u32) {
        self.data
            .entry(observation)
            .and_modify(|e| e.freq = e.freq.wrapping_add(frequency))
            .or_insert(Entry {
                freq: frequency,
                estimate: 0.0,
            });
    }

    /// Fit the estimator. Returns `false`, changing nothing, when there are
    /// fewer than five distinct observation counts — Elworthy's `MinInput`
    /// guard, since the log-log fit is meaningless on fewer points.
    pub fn analyse(&mut self) -> bool {
        let rows = self.data.len();
        if rows < 5 {
            return false;
        }
        let obs: Vec<u32> = self.data.keys().copied().collect();
        let freq: Vec<u32> = self.data.values().map(|e| e.freq).collect();

        // The total number of events observed. STAR accumulates this in u32
        // and lets it wrap; reproducing that keeps the estimates identical on
        // the inputs where it happens rather than only on the ones where it
        // does not.
        let mut big_n: u32 = 0;
        for r in 0..rows {
            big_n = big_n.wrapping_add(obs[r].wrapping_mul(freq[r]));
        }

        // The Good-Turing estimate of unseen mass: the share of events seen
        // exactly once.
        self.pzero = match self.data.get(&1) {
            Some(e) => f64::from(e.freq) / f64::from(big_n),
            None => 0.0,
        };

        // Z-transform: each frequency is averaged over the gap to its
        // neighbours, which is what makes the log-log fit stable at the sparse
        // high-count end where most counts are 0 or 1.
        let mut log_obs = vec![0.0f64; rows];
        let mut log_z = vec![0.0f64; rows];
        let (mut mean_x, mut mean_y) = (0.0f64, 0.0f64);
        let mut prev: u32 = 0;
        for r in 0..rows {
            let k = if r + 1 == rows {
                f64::from(2u32.wrapping_mul(obs[r]).wrapping_sub(prev))
            } else {
                f64::from(obs[r + 1])
            };
            let z = f64::from(2u32.wrapping_mul(freq[r])) / (k - f64::from(prev));
            log_obs[r] = f64::from(obs[r]).ln();
            log_z[r] = z.ln();
            mean_x += log_obs[r];
            mean_y += log_z[r];
            prev = obs[r];
        }
        mean_x /= rows as f64;
        mean_y /= rows as f64;

        let (mut xy, mut xx) = (0.0f64, 0.0f64);
        for r in 0..rows {
            xy += (log_obs[r] - mean_x) * (log_z[r] - mean_y);
            xx += sq(log_obs[r] - mean_x);
        }
        let slope = xy / xx;
        let intercept = mean_y - slope * mean_x;
        let smoothed = |i: u32| (intercept + slope * f64::from(i).ln()).exp();

        // For each observation count, the Turing estimate while it still
        // differs from the fitted line by more than 1.96 standard errors, and
        // the fitted line from the point they agree onwards. Once the switch
        // happens it is permanent: the raw estimates only get noisier.
        let mut r_star = vec![0.0f64; rows];
        let mut indifferent = false;
        for r in 0..rows {
            let obs1 = obs[r] + 1;
            let y = f64::from(obs1) * smoothed(obs1) / smoothed(obs[r]);
            match self.data.get(&obs1) {
                None => indifferent = true,
                Some(next) if !indifferent => {
                    let n_next = next.freq;
                    let x = f64::from(obs1.wrapping_mul(n_next)) / f64::from(freq[r]);
                    let threshold = 1.96
                        * (sq(f64::from(obs1)) * f64::from(n_next) / sq(f64::from(freq[r]))
                            * (1.0 + f64::from(n_next) / f64::from(freq[r])))
                        .sqrt();
                    if (x - y).abs() <= threshold {
                        indifferent = true;
                    } else {
                        r_star[r] = x;
                    }
                }
                Some(_) => {}
            }
            if indifferent {
                r_star[r] = y;
            }
        }

        let mut big_n_prime = 0.0f64;
        for r in 0..rows {
            big_n_prime += f64::from(freq[r]) * r_star[r];
        }
        for (r, e) in self.data.values_mut().enumerate() {
            e.estimate = (1.0 - self.pzero) * r_star[r] / big_n_prime;
        }
        true
    }

    /// The estimated probability of an event observed `observation` times.
    ///
    /// `0` gives the reserved unseen mass. An observation count that never
    /// occurred returns `None`, leaving the caller's value untouched, as
    /// Elworthy's version does.
    pub fn estimate(&self, observation: u32) -> Option<f64> {
        if observation == 0 {
            return Some(self.pzero);
        }
        self.data.get(&observation).map(|e| e.estimate)
    }

    /// The reserved unseen mass. Zero until [`analyse`](Self::analyse)
    /// succeeds — see D17 in the module docs.
    pub fn pzero(&self) -> f64 {
        self.pzero
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A well-formed spectrum: mass is reserved for the unseen, every seen
    /// count gets a finite estimate, and the total stays a probability.
    #[test]
    fn estimates_are_finite_and_sum_to_about_one() {
        let mut sgt = Sgt::new();
        // Frequency-of-frequencies of a Zipf-ish sample.
        for (obs, freq) in [(1u32, 120u32), (2, 40), (3, 24), (4, 8), (5, 4), (7, 2)] {
            sgt.add(obs, freq);
        }
        assert!(sgt.analyse());

        let p0 = sgt.estimate(0).unwrap();
        assert!(p0 > 0.0 && p0 < 1.0, "unseen mass must be a probability");

        let mut total = p0;
        for (obs, freq) in [(1u32, 120u32), (2, 40), (3, 24), (4, 8), (5, 4), (7, 2)] {
            let e = sgt.estimate(obs).unwrap();
            assert!(e.is_finite() && e > 0.0, "estimate for {obs} was {e}");
            total += e * f64::from(freq);
        }
        assert!(
            (total - 1.0).abs() < 0.05,
            "estimates should sum to roughly 1, got {total}"
        );
    }

    /// The smoothing is monotone in the observation count: an event seen more
    /// often cannot be estimated as less likely.
    #[test]
    fn a_more_frequent_observation_is_never_less_likely() {
        let mut sgt = Sgt::new();
        for (obs, freq) in [(1u32, 100u32), (2, 30), (3, 12), (4, 6), (6, 2)] {
            sgt.add(obs, freq);
        }
        assert!(sgt.analyse());
        // Starting at 1: `estimate(0)` is the *total* mass reserved for all
        // unseen events, not a per-event probability, so it does not belong in
        // this comparison.
        let mut last = 0.0f64;
        for obs in [1u32, 2, 3, 4, 6] {
            let e = sgt.estimate(obs).unwrap();
            assert!(e >= last, "estimate dropped at {obs}: {e} < {last}");
            last = e;
        }
    }

    /// D17. With fewer than five distinct counts the fit is refused, and the
    /// unseen mass stays zero instead of being whatever the stack held —
    /// which is what STAR reads on this path.
    #[test]
    fn too_few_frequencies_leaves_the_unseen_mass_at_zero() {
        let mut sgt = Sgt::new();
        for (obs, freq) in [(1u32, 10u32), (2, 4), (3, 1)] {
            sgt.add(obs, freq);
        }
        assert!(!sgt.analyse(), "fewer than 5 distinct counts must not fit");
        // Exact zero is the point: not "small", but never written.
        assert_eq!(sgt.pzero().to_bits(), 0.0f64.to_bits());
        assert_eq!(sgt.estimate(0).map(f64::to_bits), Some(0.0f64.to_bits()));
    }

    /// An observation count that never occurred has no estimate, rather than a
    /// zero that would be indistinguishable from a real one.
    #[test]
    fn an_unobserved_count_has_no_estimate() {
        let mut sgt = Sgt::new();
        for (obs, freq) in [(1u32, 50u32), (2, 20), (3, 10), (4, 5), (5, 2)] {
            sgt.add(obs, freq);
        }
        assert!(sgt.analyse());
        assert_eq!(sgt.estimate(9), None);
    }

    /// No observation of count 1 means nothing was seen exactly once, so there
    /// is no evidence of unseen events and no mass is reserved.
    #[test]
    fn without_singletons_no_mass_is_reserved() {
        let mut sgt = Sgt::new();
        for (obs, freq) in [(2u32, 40u32), (3, 20), (4, 10), (5, 5), (6, 2)] {
            sgt.add(obs, freq);
        }
        assert!(sgt.analyse());
        assert_eq!(sgt.estimate(0).map(f64::to_bits), Some(0.0f64.to_bits()));
    }
}
