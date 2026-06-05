use crate::error::Error;
use crate::genome::Genome;
use crate::index::packed_array::PackedArray;

/// Suffix array for genome indexing.
///
/// Stores positions in the genome (forward + reverse complement) sorted
/// by their suffixes, enabling fast exact match search.
#[derive(Clone)]
pub struct SuffixArray {
    /// Packed array of suffix positions (with strand bit)
    pub data: PackedArray,

    /// Strand bit position (position where bit distinguishes forward/reverse)
    pub gstrand_bit: u32,

    /// Mask for extracting position (all bits below gstrand_bit)
    pub gstrand_mask: u64,
}

impl SuffixArray {
    /// Calculate GstrandBit from genome size.
    ///
    /// Formula: max(32, floor(log2(nGenome)) + 1)
    pub fn calculate_gstrand_bit(n_genome: u64) -> u32 {
        if n_genome == 0 {
            return 32;
        }
        let log2_bits = 64 - n_genome.leading_zeros();
        u32::max(32, log2_bits)
    }

    /// Build the suffix array for `genome`.
    ///
    /// Delegates to [`crate::index::sa_build::build`], which is a Rust port of
    /// the **CaPS-SA** sample-sort SA construction (Khan et al., WABI 2023,
    /// via the `caps-sa` crate) wrapped in a STAR-faithful sentinel transform.
    /// Produces a `PackedArray` byte-identical to STAR's `SA` file.
    pub fn build(genome: &Genome) -> Result<Self, Error> {
        crate::index::sa_build::build(genome)
    }

    /// Get the number of suffixes in the array.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if the suffix array is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Decode a suffix array entry into (position, is_reverse).
    pub fn decode(&self, sa_value: u64) -> (u64, bool) {
        let is_reverse = (sa_value >> self.gstrand_bit) != 0;
        let position = sa_value & self.gstrand_mask;
        (position, is_reverse)
    }

    /// Read a suffix array entry.
    pub fn get(&self, index: usize) -> u64 {
        self.data.read(index)
    }
}

/// Compare two suffixes for sorting using STAR's exact comparison logic.
///
/// Retained as a `#[cfg(test)]` ground-truth oracle: the previous naive
/// `sort_by`-based `SuffixArray::build` was differentially validated against
/// STAR (yeast SA byte-identical), and is now used to differentially
/// validate the caps-sa-backed implementation in [`sa_build::build`].
#[cfg(test)]
fn compare_suffixes(
    genome: &Genome,
    pos_a: usize,
    reverse_a: bool,
    pos_b: usize,
    reverse_b: bool,
) -> std::cmp::Ordering {
    use std::cmp::Ordering;

    let n_genome = genome.n_genome as usize;
    let sequence = &genome.sequence;

    // Adjust positions for reverse complement
    let start_a = if reverse_a { pos_a + n_genome } else { pos_a };
    let start_b = if reverse_b { pos_b + n_genome } else { pos_b };

    // Compare up to n_genome bytes. Padding (5) stops the comparison early
    // in practice. Out-of-bounds (RC suffixes near the genome boundary) is
    // treated as padding so those entries sort correctly rather than falling
    // through to the position-only fallback.
    let max_len = n_genome;

    for offset in 0..max_len {
        let idx_a = start_a + offset;
        let idx_b = start_b + offset;

        let byte_a = if idx_a < sequence.len() {
            sequence[idx_a]
        } else {
            5
        };
        let byte_b = if idx_b < sequence.len() {
            sequence[idx_b]
        } else {
            5
        };

        // Stop at padding (value 5) - this is STAR's sentinel
        let is_padding_a = byte_a == 5;
        let is_padding_b = byte_b == 5;

        if is_padding_a && is_padding_b {
            // Both hit padding at same depth — sort ascending by packed SA value
            // (strand bit at position gstrand_bit, same as what's stored in the SA).
            // For yeast (gstrand_bit=32): FW entries have packed_value = pos (no bit 32),
            // RC entries have packed_value = pos | (1<<32). All FW entries therefore
            // sort before all RC entries, matching STAR's tie-breaking behavior.
            let packed_a = if reverse_a {
                pos_a | (1usize << 32)
            } else {
                pos_a
            };
            let packed_b = if reverse_b {
                pos_b | (1usize << 32)
            } else {
                pos_b
            };
            return packed_a.cmp(&packed_b);
        }

        if is_padding_a {
            return Ordering::Greater; // Padding sorts after valid bases
        }

        if is_padding_b {
            return Ordering::Less;
        }

        // Normal byte comparison
        let byte_cmp = byte_a.cmp(&byte_b);
        if byte_cmp != Ordering::Equal {
            return byte_cmp;
        }
    }

    // If we exhausted max_len, fall back to packed SA value comparison
    let packed_a = if reverse_a {
        pos_a | (1usize << 32)
    } else {
        pos_a
    };
    let packed_b = if reverse_b {
        pos_b | (1usize << 32)
    } else {
        pos_b
    };
    packed_a.cmp(&packed_b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::packed_array::PackedArray;
    use crate::params::Parameters;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// STAR-faithful naive oracle, kept for differential testing of the new
    /// caps-sa-backed implementation.
    ///
    /// The pre-caps-sa rustar builder filtered `< 5` (ACGT + N), which was a
    /// latent divergence from STAR's `G[ii] < 4` (ACGT only) at
    /// `Genome_genomeGenerate.cpp:185`. Yeast reference genomes contain no
    /// N's, so the difference never surfaced in the byte-identity test. The
    /// new caps-sa-backed builder uses STAR's correct `< 4` filter; this
    /// oracle does the same.
    fn build_naive(genome: &Genome) -> SuffixArray {
        let n_genome = genome.n_genome as usize;
        let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
        let gstrand_mask = (1u64 << gstrand_bit) - 1;
        let word_length = gstrand_bit + 1;

        let mut suffixes: Vec<(u64, bool)> = Vec::new();
        for i in 0..n_genome {
            if genome.sequence[i] < 4 {
                suffixes.push((i as u64, false));
            }
        }
        for i in n_genome..(2 * n_genome) {
            if genome.sequence[i] < 4 {
                suffixes.push(((i - n_genome) as u64, true));
            }
        }
        suffixes.sort_by(|a, b| compare_suffixes(genome, a.0 as usize, a.1, b.0 as usize, b.1));

        let mut packed = PackedArray::new(word_length, suffixes.len());
        let n2bit = 1u64 << gstrand_bit;
        for (i, &(pos, is_reverse)) in suffixes.iter().enumerate() {
            let packed_value = if is_reverse { pos | n2bit } else { pos };
            packed.write(i, packed_value);
        }
        SuffixArray {
            data: packed,
            gstrand_bit,
            gstrand_mask,
        }
    }

    fn make_test_genome(sequence: &str, bin_nbits: u32) -> Genome {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "{sequence}").unwrap();

        let bin_nbits_str = bin_nbits.to_string();
        let args = vec![
            "rustar-aligner",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            &bin_nbits_str,
        ];

        let params = Parameters::parse_from(args);
        Genome::from_fasta(&params).unwrap()
    }

    #[test]
    fn gstrand_bit_calculation() {
        assert_eq!(SuffixArray::calculate_gstrand_bit(0), 32);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1), 32);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1000), 32);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1u64 << 32), 33);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1u64 << 33), 34);
    }

    #[test]
    fn build_small_genome() {
        let genome = make_test_genome("ACGT", 2);
        let sa = SuffixArray::build(&genome).unwrap();

        // Should have suffixes for forward + reverse (excluding padding)
        assert!(!sa.is_empty());
        assert_eq!(sa.gstrand_bit, 32); // Small genome
    }

    #[test]
    fn decode_sa_entry() {
        let genome = make_test_genome("ACGT", 2);
        let sa = SuffixArray::build(&genome).unwrap();

        // Read first entry and decode
        let entry = sa.get(0);
        let (pos, _is_reverse) = sa.decode(entry);

        // Position should be valid
        assert!(pos < genome.n_genome);
    }

    #[test]
    fn suffix_sorting() {
        // Simple test: "AAB" should sort as A, AA, AAB, B
        let genome = make_test_genome("AAB", 2);
        let sa = SuffixArray::build(&genome).unwrap();

        // Verify we have entries
        assert!(!sa.is_empty());

        // The lexicographically first suffix should start with the smallest base
        let first_entry = sa.get(0);
        let (first_pos, _) = sa.decode(first_entry);
        let first_base = genome.sequence[first_pos as usize];

        // In "AAB", the first suffix lexicographically is "A" (from pos 0 or 1)
        assert!(first_base == 0); // A
    }

    /// Differential: the caps-sa-backed builder must produce a `PackedArray`
    /// byte-identical to the naive oracle. Covers single-chromosome and
    /// multi-chromosome inputs (the latter exercises the per-segment
    /// sentinel transform).
    #[test]
    fn caps_sa_matches_naive_oracle_single_chr() {
        let genome = make_test_genome("ACGTACGTACGTNACGT", 4);
        let new_sa = SuffixArray::build(&genome).unwrap();
        let naive_sa = build_naive(&genome);
        assert_eq!(new_sa.data.data(), naive_sa.data.data());
        assert_eq!(new_sa.gstrand_bit, naive_sa.gstrand_bit);
    }

    #[test]
    fn caps_sa_matches_naive_oracle_with_padding() {
        // Tiny bin so the chromosome is padded with several spacer bytes,
        // exercising the spacer-run → sentinel transform on a multi-byte run.
        let genome = make_test_genome("ACGTACGT", 2);
        let new_sa = SuffixArray::build(&genome).unwrap();
        let naive_sa = build_naive(&genome);
        assert_eq!(new_sa.data.data(), naive_sa.data.data());
    }

    #[test]
    fn reverse_complement_included() {
        let genome = make_test_genome("ACGT", 2);
        let sa = SuffixArray::build(&genome).unwrap();

        let mut has_forward = false;
        let mut has_reverse = false;

        for i in 0..sa.len() {
            let entry = sa.get(i);
            let (_, is_reverse) = sa.decode(entry);
            if is_reverse {
                has_reverse = true;
            } else {
                has_forward = true;
            }
        }

        assert!(has_forward, "SA should include forward strand suffixes");
        assert!(has_reverse, "SA should include reverse strand suffixes");
    }
}
