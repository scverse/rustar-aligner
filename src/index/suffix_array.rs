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

    /// Build suffix array from genome.
    ///
    /// This is a simplified implementation that works for small genomes.
    /// For production, STAR uses prefix bucketing and parallel sorting.
    pub fn build(genome: &Genome) -> Result<Self, Error> {
        let n_genome = genome.n_genome as usize;
        let gstrand_bit = Self::calculate_gstrand_bit(genome.n_genome);
        let gstrand_mask = (1u64 << gstrand_bit) - 1;
        let word_length = gstrand_bit + 1;

        // Create array of (position, is_reverse) tuples for all valid suffixes
        let mut suffixes = Vec::new();

        // Add forward strand suffixes
        for i in 0..n_genome {
            // Only include positions that start with a valid base (not padding)
            if genome.sequence[i] < 5 {
                suffixes.push((i as u64, false));
            }
        }

        // Add reverse strand suffixes
        for i in n_genome..(2 * n_genome) {
            if genome.sequence[i] < 5 {
                suffixes.push(((i - n_genome) as u64, true));
            }
        }

        let sa_length = suffixes.len();
        log::info!(
            "Building suffix array with {} entries (gstrand_bit={}, word_length={})",
            sa_length,
            gstrand_bit,
            word_length
        );

        // Sort suffixes using custom comparator
        suffixes.sort_by(|a, b| compare_suffixes(genome, a.0 as usize, a.1, b.0 as usize, b.1));

        // Pack into PackedArray
        let mut packed = PackedArray::new(word_length, sa_length);
        let n2bit = 1u64 << gstrand_bit;

        for (i, &(pos, is_reverse)) in suffixes.iter().enumerate() {
            let packed_value = if is_reverse {
                pos | n2bit // Set strand bit for reverse
            } else {
                pos
            };
            packed.write(i, packed_value);
        }

        Ok(SuffixArray {
            data: packed,
            gstrand_bit,
            gstrand_mask,
        })
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

/// Compare two suffixes for sorting.
///
/// Implements STAR's comparison logic:
/// - Compares up to 8-byte words at a time
/// - Stops at padding (value 5)
/// - Uses anti-stable sort when both hit padding at same depth
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
        match byte_a.cmp(&byte_b) {
            Ordering::Equal => continue,
            other => return other,
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
    use crate::params::Parameters;
    use clap::Parser;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_test_genome(sequence: &str, bin_nbits: u32) -> Genome {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "{}", sequence).unwrap();

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
