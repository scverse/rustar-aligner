//! Bit-exact ports of the three libc++ random facilities STARsolo's
//! `EmptyDrops_CR` depends on.
//!
//! STAR's Monte-Carlo rescue draws from `std::mt19937`, converts to doubles
//! with `std::generate_canonical<double, 53>`, and samples categories with
//! `std::discrete_distribution`. All three are implementation-defined in the
//! parts that matter: the standard fixes `mt19937`'s output but not how
//! `generate_canonical` consumes it, and says nothing about how
//! `discrete_distribution` maps a uniform draw onto categories. So "port the
//! algorithm" is not enough — it has to be *libc++'s* algorithm, because that
//! is what STAR was built against and what its numbers come out of.
//!
//! Every value asserted in the tests below was produced by compiling a C++
//! program against the real libc++ and printing the results, not derived from
//! reading the source. That is the only way to be sure.
//!
//! `solo::count` currently samples with a `SplitMix64` stream, under a comment
//! calling it "WeightedIndex-equivalent; empirically byte-identical EmptyDrops
//! cell calls". That claim cannot hold in general: two unrelated generators
//! cannot agree on an arbitrary number of draws, so it is true of the cases
//! that happened to be checked and unknown everywhere else. These types remove
//! the guesswork. Wiring them into the EmptyDrops path is the next step and is
//! deliberately separate, since it moves cell calls.

/// libc++'s `std::mt19937`.
///
/// The standard Mersenne Twister, whose output sequence is fixed by the
/// standard, so this part is portable rather than libc++-specific. It is here
/// because the two facilities that follow are not.
#[derive(Debug, Clone)]
pub struct Mt19937 {
    state: [u32; Self::N],
    index: usize,
}

impl Mt19937 {
    const N: usize = 624;
    const M: usize = 397;
    const MATRIX_A: u32 = 0x9908_b0df;
    const UPPER_MASK: u32 = 0x8000_0000;
    const LOWER_MASK: u32 = 0x7fff_ffff;

    /// Seed exactly as `std::mt19937(seed)` does.
    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; Self::N];
        state[0] = seed;
        for i in 1..Self::N {
            state[i] = 1_812_433_253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self {
            state,
            index: Self::N,
        }
    }

    /// One draw, equivalent to `operator()`.
    pub fn next_u32(&mut self) -> u32 {
        if self.index >= Self::N {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    fn twist(&mut self) {
        for i in 0..Self::N {
            let y = (self.state[i] & Self::UPPER_MASK)
                | (self.state[(i + 1) % Self::N] & Self::LOWER_MASK);
            let mut next = self.state[(i + Self::M) % Self::N] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= Self::MATRIX_A;
            }
            self.state[i] = next;
        }
        self.index = 0;
    }

    /// libc++'s `std::generate_canonical<double, 53>`.
    ///
    /// This is where implementations diverge. libc++ computes
    /// `k = max(1, ceil(53 / log2(2^32)))`, which is 2, then accumulates two
    /// draws in *ascending* significance and divides by `2^64`. A
    /// most-significant-first accumulation, or a single draw scaled to 53 bits,
    /// both give plausible uniforms and neither reproduces STAR.
    pub fn canonical_f64(&mut self) -> f64 {
        // r = 2^32, k = 2, so the base is r^k = 2^64.
        let base = 2f64.powi(64);
        let mut sum = 0f64;
        let mut factor = 1f64;
        for _ in 0..2 {
            sum += f64::from(self.next_u32()) * factor;
            factor *= 2f64.powi(32);
        }
        sum / base
    }
}

/// libc++'s `std::discrete_distribution`.
///
/// Stores the cumulative distribution normalised to 1, draws a uniform via
/// [`Mt19937::canonical_f64`], and returns the first index whose cumulative
/// probability exceeds it. The final category is the fallback, so a draw of
/// exactly 1.0 cannot fall off the end.
#[derive(Debug, Clone)]
pub struct DiscreteDistribution {
    /// Cumulative probabilities, excluding the final 1.0.
    cumulative: Vec<f64>,
}

impl DiscreteDistribution {
    /// Build from unnormalised weights, as `discrete_distribution(w)` does.
    pub fn new(weights: &[f64]) -> Self {
        let total: f64 = weights.iter().sum();
        let mut cumulative = Vec::with_capacity(weights.len().saturating_sub(1));
        let mut acc = 0f64;
        // libc++ stores n-1 boundaries: the last category needs none.
        for w in weights.iter().take(weights.len().saturating_sub(1)) {
            acc += w;
            cumulative.push(if total > 0.0 { acc / total } else { 0.0 });
        }
        Self { cumulative }
    }

    /// One sample.
    pub fn sample(&self, rng: &mut Mt19937) -> usize {
        let u = rng.canonical_f64();
        self.cumulative.partition_point(|&c| c <= u)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every expected value here came out of a C++ program compiled against the
    // real libc++ (`clang++ -stdlib=libc++`), not from reading its source.

    #[test]
    fn mt19937_matches_libcxx_stream() {
        let mut g = Mt19937::new(19_760_110);
        let got: Vec<u32> = (0..10).map(|_| g.next_u32()).collect();
        assert_eq!(
            got,
            vec![
                2_116_612_583,
                978_492_435,
                3_413_959_089,
                853_152_524,
                3_057_288_333,
                81_846_811,
                724_235_003,
                450_930_519,
                3_920_508_463,
                4_192_617_403,
            ]
        );
    }

    #[test]
    fn generate_canonical_matches_libcxx_bit_for_bit() {
        let mut g = Mt19937::new(19_760_110);
        // Bit patterns rather than decimal literals: the decimal forms carry
        // more digits than an f64 holds, and the point is that the double is
        // identical down to the last place. A difference there changes which
        // category a sample lands in.
        let expected: [u64; 5] = [
            0x3fcd_294e_09bf_1479, // 0.22782302356623399
            0x3fc9_6d09_8665_be71, // 0.19864005148291455
            0x3f93_8388_6ed8_ea12, // 0.019056445851882549
            0x3fba_e0a7_572b_2af3, // 0.10499044302120435
            0x3fef_3cc8_777d_35c7, // 0.97616980874743653
        ];
        for (i, &want) in expected.iter().enumerate() {
            let got = g.canonical_f64();
            assert_eq!(
                got.to_bits(),
                want,
                "draw {i}: got {got:.17}, libc++ gives {:.17}",
                f64::from_bits(want)
            );
        }
    }

    #[test]
    fn discrete_distribution_matches_libcxx_integer_weights() {
        let mut g = Mt19937::new(19_760_110);
        let d = DiscreteDistribution::new(&[1.0, 2.0, 3.0, 4.0]);
        let got: Vec<usize> = (0..20).map(|_| d.sample(&mut g)).collect();
        assert_eq!(
            got,
            vec![1, 1, 0, 1, 3, 3, 0, 3, 3, 3, 2, 1, 1, 3, 3, 1, 3, 0, 3, 2]
        );
    }

    #[test]
    fn a_single_category_always_wins() {
        let mut g = Mt19937::new(1);
        let d = DiscreteDistribution::new(&[5.0]);
        for _ in 0..8 {
            assert_eq!(d.sample(&mut g), 0);
        }
    }

    #[test]
    fn zero_weight_categories_are_never_drawn() {
        let mut g = Mt19937::new(42);
        let d = DiscreteDistribution::new(&[0.0, 1.0, 0.0]);
        for _ in 0..64 {
            assert_eq!(d.sample(&mut g), 1);
        }
    }
}
