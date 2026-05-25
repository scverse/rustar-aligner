use std::fs::File;
use std::os::unix::fs::FileExt;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use crate::error::Error;
use crate::genome::Genome;
use crate::index::packed_array::PackedArray;
use crate::index::suffix_array::SuffixArray;


/// SA index for fast k-mer lookup during binary search.
///
/// For each k-mer prefix (up to genomeSAindexNbases length), stores the
/// SA range where suffixes starting with that prefix can be found.
#[derive(Clone)]
pub struct SaIndex {
    /// Length of indexed k-mers (typically 14)
    pub nbases: u32,

    /// Cumulative counts for each k-mer level
    /// Length = nbases + 1
    /// genomeSAindexStart[k] = total number of k-mers of length <= k
    pub genome_sa_index_start: Vec<u64>,

    /// Packed array storing SA index entries
    /// Each entry is (gstrand_bit + 3) bits:
    /// - Bits 0..gstrand_bit: SA position (or max value if prefix absent)
    /// - Bit gstrand_bit+1: "contains N" flag
    /// - Bit gstrand_bit+2: "prefix absent" flag
    pub data: PackedArray,

    /// Word length for SAindex entries (gstrand_bit + 3)
    pub word_length: u32,

    /// Strand bit position (from SA)
    pub gstrand_bit: u32,
}

impl SaIndex {
    /// Calculate the total number of k-mer prefixes to index.
    ///
    /// Returns sum of 4^1 + 4^2 + ... + 4^nbases = (4^(nbases+1) - 4) / 3
    pub fn calculate_num_indices(nbases: u32) -> u64 {
        if nbases == 0 {
            return 0;
        }
        // Sum of geometric series: (4^(n+1) - 4) / 3
        let power = 4u64.pow(nbases + 1);
        (power - 4) / 3
    }

    /// Parallel SAindex build from a memory-mapped SA byte buffer.
    ///
    /// The streaming `genomeGenerate` path writes the SA to disk via
    /// [`PackedStreamWriter`][crate::index::packed_stream::PackedStreamWriter]
    /// and then calls this entry. caps-sa's emit no longer touches
    /// the SAindex, so the ~16 min of per-emit serial k-mer work
    /// that previously sat on top of caps-sa's parallel SA build is
    /// replaced by a single parallel pass here.
    ///
    /// ## Algorithm
    ///
    /// For each k-mer slot `(k, kmer_idx)` we want the **smallest**
    /// `sa_idx` whose suffix starts with that k-mer (which is the
    /// "first occurrence" in lex order). The sequential algorithm
    /// processes SA entries in order and only writes if the slot is
    /// still "absent"; we generalise to parallel by using
    /// `AtomicU64::fetch_min` over a temporary `Vec<AtomicU64>` of
    /// length `num_indices` (~358 M entries × 8 B = ~2.86 GB on
    /// the human genome). Each rayon worker takes a contiguous
    /// chunk of SA indices; per-entry it decodes the packed value
    /// from `sa_bytes` (read via [`read_packed_entry`] — same bit
    /// layout as [`PackedArray::read`]), computes the genome offset,
    /// extracts k-mers incrementally (one base read per k, break on
    /// N) and atomic-min's `sa_idx` into the slot.
    ///
    /// After all workers finish, a final sequential pass walks the
    /// `firsts[]` array and packs each entry into the
    /// [`PackedArray`] SAindex output (`u64::MAX` sentinel →
    /// `absent_marker`, else → `sa_idx`).
    ///
    /// ## Memory
    ///
    /// - `firsts: Vec<AtomicU64>`: ~2.86 GB transient.
    /// - `PackedArray` output: ~1.5 GB (kept).
    /// - SA bytes: not in process RSS (kernel page cache, mmap'd).
    /// - `genome.sequence`: shared with the caller, ~6.3 GB.
    pub fn build_parallel(
        genome: &Genome,
        sa_file: &File,
        sa_word_length: u32,
        gstrand_bit: u32,
        gstrand_mask: u64,
        n_entries: usize,
        nbases: u32,
    ) -> Result<Self, Error> {
        let sai_word_length = gstrand_bit + 3;
        let mut genome_sa_index_start = vec![0u64; (nbases + 1) as usize];
        for k in 1..=nbases {
            genome_sa_index_start[k as usize] =
                genome_sa_index_start[(k - 1) as usize] + 4u64.pow(k);
        }
        let num_indices = Self::calculate_num_indices(nbases) as usize;
        log::info!(
            "Building SA index in parallel: nbases={nbases}, num_indices={num_indices}, \
             n_entries={n_entries}, threads={}",
            rayon::current_num_threads()
        );

        // Sentinel `u64::MAX` means "no SA index has reached this
        // k-mer yet". `fetch_min` makes any real `sa_idx` smaller
        // than the sentinel, so the first thread to touch a slot
        // wins; subsequent threads only narrow the value to the
        // smallest sa_idx seen across all workers. Total order is
        // monotonic, so the final value is exactly the
        // first-in-lex-order `sa_idx`.
        let firsts: Vec<AtomicU64> = (0..num_indices)
            .into_par_iter()
            .map(|_| AtomicU64::new(u64::MAX))
            .collect();

        let sa_mask: u64 = if sa_word_length == 64 {
            u64::MAX
        } else {
            (1u64 << sa_word_length) - 1
        };
        let n_genome = genome.n_genome as usize;
        let genome_seq: &[u8] = &genome.sequence;

        // Chunk by entries — each worker reads a contiguous byte
        // range from `sa_file` via `read_at` (= pread; thread-safe
        // on a shared `&File`). This avoids `memmap2::Mmap`, whose
        // touched pages get counted in process RSS (`RssFile`); a
        // 24 GB SA would push the peak up by 24 GB. `read_at` just
        // hits the kernel page cache, which is kernel-side memory
        // and **not** counted in the process's resident set.
        //
        // Chunk size is 1 M entries (~4 MB at 33-bit `word_length`),
        // giving 32 × 4 MB ≈ 128 MB of in-flight per-worker buffer
        // peak. Picking too small a chunk amplifies the per-task
        // overhead (the SAindex inner loop only does ~14 work units
        // per entry); too large blows up the buffer per worker.
        const ENTRIES_PER_CHUNK: usize = 1 << 20;
        let n_chunks = n_entries.div_ceil(ENTRIES_PER_CHUNK);

        (0..n_chunks)
            .into_par_iter()
            .try_for_each(|chunk_idx| -> std::io::Result<()> {
                let chunk_start = chunk_idx * ENTRIES_PER_CHUNK;
                let chunk_end = (chunk_start + ENTRIES_PER_CHUNK).min(n_entries);
                let chunk_n = chunk_end - chunk_start;

                // Bit and byte ranges in the SA file covered by this
                // chunk's entries. The first entry may start
                // mid-byte (`start_bit_shift > 0`); we include the
                // partial start byte and pad the read by 8 bytes
                // so the 8-byte LE load at the last entry never
                // goes out of bounds.
                let start_bit = chunk_start as u64 * sa_word_length as u64;
                let end_bit = chunk_end as u64 * sa_word_length as u64;
                let start_byte = start_bit / 8;
                let start_bit_shift = (start_bit % 8) as u32;
                let end_byte_excl = end_bit.div_ceil(8);
                let bytes_to_read = (end_byte_excl - start_byte) as usize + 8;
                let mut buf = vec![0u8; bytes_to_read];
                // read_at may return short reads near EOF — that's
                // fine because we pre-zeroed `buf`. The +8 padding
                // bytes don't exist on disk and stay zero, matching
                // what `PackedArray::data_byte_len_for`'s `+ 8`
                // would have given an in-RAM PackedArray.
                let _ = sa_file.read_at(&mut buf, start_byte)?;

                for i in 0..chunk_n {
                    let local_bit =
                        i as u64 * sa_word_length as u64 + start_bit_shift as u64;
                    let local_byte = (local_bit / 8) as usize;
                    let local_shift = (local_bit % 8) as u32;
                    let bytes: &[u8; 8] = (&buf[local_byte..local_byte + 8])
                        .try_into()
                        .expect("padded buf must have 8 bytes at every entry's start");
                    let word = u64::from_le_bytes(*bytes);
                    let packed = (word >> local_shift) & sa_mask;

                    let pos = packed & gstrand_mask;
                    let is_reverse = (packed >> gstrand_bit) != 0;
                    let genome_pos = if is_reverse {
                        pos as usize + n_genome
                    } else {
                        pos as usize
                    };
                    let sa_idx = chunk_start + i;

                    let mut kmer_idx: u64 = 0;
                    for k in 1..=nbases {
                        if genome_pos + (k as usize) > genome_seq.len() {
                            break;
                        }
                        let next_base = genome_seq[genome_pos + (k - 1) as usize];
                        if next_base >= 4 {
                            break;
                        }
                        kmer_idx = (kmer_idx << 2) | (next_base as u64);
                        let sai_pos = genome_sa_index_start[(k - 1) as usize] + kmer_idx;
                        firsts[sai_pos as usize].fetch_min(sa_idx as u64, Ordering::Relaxed);
                    }
                }
                Ok(())
            })
            .map_err(|e| Error::Index(format!("pread during SAindex build: {e}")))?;

        // Final sequential pack into the output `PackedArray`. The
        // `+8`-byte pad at the end of the buffer is reserved by
        // `PackedArray::new`. Per-slot work is ~10 ns × 358 M ≈
        // 3.6 s on the human genome — negligible vs. the parallel
        // scan.
        let mut data = PackedArray::new(sai_word_length, num_indices);
        let absent_marker = (1u64 << (gstrand_bit + 2)) | ((1u64 << gstrand_bit) - 1);
        for (i, slot) in firsts.iter().enumerate() {
            let v = slot.load(Ordering::Relaxed);
            if v == u64::MAX {
                data.write(i, absent_marker);
            } else {
                data.write(i, v);
            }
        }
        drop(firsts);

        Ok(SaIndex {
            nbases,
            genome_sa_index_start,
            data,
            word_length: sai_word_length,
            gstrand_bit,
        })
    }

    /// Streaming builder: accepts (sa_idx, packed_value) pairs in SA
    /// order and updates the per-k-mer first-occurrence entries
    /// online. Used by the streaming `genomeGenerate` path so the
    /// 25 GB SA `PackedArray` never has to be materialised in RAM —
    /// caps-sa emits each entry, the SA file writer takes one copy,
    /// and this builder takes another (a few bytes per emit).
    ///
    /// Holds a shared reference to `genome` for the k-mer extraction
    /// `genome.sequence[genome_pos..]` reads.
    pub fn streaming_builder<'a>(
        genome: &'a Genome,
        gstrand_bit: u32,
        gstrand_mask: u64,
        nbases: u32,
    ) -> SaIndexBuilder<'a> {
        let word_length = gstrand_bit + 3;
        let mut genome_sa_index_start = vec![0u64; (nbases + 1) as usize];
        for k in 1..=nbases {
            genome_sa_index_start[k as usize] =
                genome_sa_index_start[(k - 1) as usize] + 4u64.pow(k);
        }
        let num_indices = Self::calculate_num_indices(nbases);
        let mut data = PackedArray::new(word_length, num_indices as usize);
        let absent_marker = (1u64 << (gstrand_bit + 2)) | ((1u64 << gstrand_bit) - 1);
        for i in 0..num_indices as usize {
            data.write(i, absent_marker);
        }
        SaIndexBuilder {
            nbases,
            word_length,
            gstrand_bit,
            gstrand_mask,
            genome_sa_index_start,
            data,
            genome,
            sa_idx: 0,
        }
    }

    /// Build SA index from genome and sorted suffix array.
    ///
    /// # Arguments
    /// * `genome` - The genome sequence
    /// * `sa` - The sorted suffix array
    /// * `nbases` - K-mer length to index (genomeSAindexNbases, typically 14)
    pub fn build(genome: &Genome, sa: &SuffixArray, nbases: u32) -> Result<Self, Error> {
        let gstrand_bit = sa.gstrand_bit;
        let word_length = gstrand_bit + 3; // +2 for flags, +1 for STAR's formula

        // Calculate genomeSAindexStart array
        let mut genome_sa_index_start = vec![0u64; (nbases + 1) as usize];
        genome_sa_index_start[0] = 0;
        for k in 1..=nbases {
            genome_sa_index_start[k as usize] =
                genome_sa_index_start[(k - 1) as usize] + 4u64.pow(k);
        }

        let num_indices = Self::calculate_num_indices(nbases);

        log::info!(
            "Building SA index: nbases={nbases}, num_indices={num_indices}, word_length={word_length}"
        );

        // Initialize packed array with "absent" markers
        let mut data = PackedArray::new(word_length, num_indices as usize);
        let absent_marker = (1u64 << (gstrand_bit + 2)) | ((1u64 << gstrand_bit) - 1);

        for i in 0..num_indices as usize {
            data.write(i, absent_marker);
        }

        // Iterate through SA and record first occurrence of each k-mer.
        // Inner k-loop maintains `kmer_idx` **incrementally** — one
        // base read per k iteration (vs. the original `O(k²)` read
        // pattern that re-scanned the prefix for every k). When an N
        // is encountered at position `genome_pos + (k - 1)` we
        // `break` rather than `continue`: every longer k-mer at this
        // same `genome_pos` necessarily includes that N too.
        for sa_idx in 0..sa.len() {
            let sa_entry = sa.get(sa_idx);
            let (pos, is_reverse) = sa.decode(sa_entry);
            let genome_pos = if is_reverse {
                pos as usize + genome.n_genome as usize
            } else {
                pos as usize
            };

            let mut kmer_idx: u64 = 0;
            for k in 1..=nbases {
                if genome_pos + (k as usize) > genome.sequence.len() {
                    break;
                }
                let next_base = genome.sequence[genome_pos + (k - 1) as usize];
                if next_base >= 4 {
                    break;
                }
                kmer_idx = (kmer_idx << 2) | (next_base as u64);

                let sai_pos = genome_sa_index_start[(k - 1) as usize] + kmer_idx;
                let current_entry = data.read(sai_pos as usize);
                let is_absent = (current_entry >> (gstrand_bit + 2)) & 1 != 0;
                if is_absent {
                    data.write(sai_pos as usize, sa_idx as u64);
                }
            }
        }

        Ok(SaIndex {
            nbases,
            genome_sa_index_start,
            data,
            word_length,
            gstrand_bit,
        })
    }
}

/// Online accumulator for [`SaIndex`]'s per-k-mer first-occurrence
/// entries. See [`SaIndex::streaming_builder`] for the constructor;
/// feed each emitted SA entry through [`SaIndexBuilder::add`] in SA
/// order, then call [`SaIndexBuilder::finish`] to obtain the
/// [`SaIndex`].
///
/// The streaming `genomeGenerate` path interleaves this with the
/// on-disk SA writer: each caps-sa emit goes to both, and the
/// 25 GB-class SA `PackedArray` is never materialised in RAM.
pub struct SaIndexBuilder<'a> {
    nbases: u32,
    word_length: u32,
    gstrand_bit: u32,
    gstrand_mask: u64,
    genome_sa_index_start: Vec<u64>,
    data: PackedArray,
    genome: &'a Genome,
    /// 0-based index of the next entry to be added. Updated by
    /// [`add`][Self::add] after each call.
    sa_idx: usize,
}

impl<'a> SaIndexBuilder<'a> {
    /// Feed the next SA entry (the strand-bit-encoded packed value,
    /// exactly what would be at `SuffixArray::get(self.sa_idx())`).
    /// Updates the first-occurrence entry of every k-mer of length
    /// 1..=nbases that starts at the SA entry's genome position,
    /// matching the body of [`SaIndex::build`]'s main loop.
    pub fn add(&mut self, packed_value: u64) {
        let pos = packed_value & self.gstrand_mask;
        let is_reverse = (packed_value >> self.gstrand_bit) != 0;
        let genome_pos = if is_reverse {
            pos as usize + self.genome.n_genome as usize
        } else {
            pos as usize
        };

        // Incremental k-mer index: one base read per k, break on N.
        // See the matching loop in `SaIndex::build` for the
        // correctness argument.
        let mut kmer_idx: u64 = 0;
        for k in 1..=self.nbases {
            if genome_pos + (k as usize) > self.genome.sequence.len() {
                break;
            }
            let next_base = self.genome.sequence[genome_pos + (k - 1) as usize];
            if next_base >= 4 {
                break;
            }
            kmer_idx = (kmer_idx << 2) | (next_base as u64);

            let sai_pos = self.genome_sa_index_start[(k - 1) as usize] + kmer_idx;
            let current_entry = self.data.read(sai_pos as usize);
            let is_absent = (current_entry >> (self.gstrand_bit + 2)) & 1 != 0;
            if is_absent {
                self.data.write(sai_pos as usize, self.sa_idx as u64);
            }
        }
        self.sa_idx += 1;
    }

    /// Number of entries fed so far.
    pub fn len(&self) -> usize {
        self.sa_idx
    }

    /// Finalise the builder into a [`SaIndex`].
    pub fn finish(self) -> SaIndex {
        SaIndex {
            nbases: self.nbases,
            genome_sa_index_start: self.genome_sa_index_start,
            data: self.data,
            word_length: self.word_length,
            gstrand_bit: self.gstrand_bit,
        }
    }
}

impl SaIndex {
    /// Hierarchical SAindex lookup (STAR's maxMappableLength2strands approach).
    ///
    /// Starts with full k-mer at level min(k, nbases), progressively shortens
    /// until a present entry is found. Returns narrowed SA range + matched level.
    ///
    /// Returns: Some((sa_start, sa_end_exclusive, matched_level, bounds_tight))
    ///   - sa_start: first SA index in range
    ///   - sa_end_exclusive: past-the-end SA index
    ///   - matched_level: how many bases the SAindex resolved
    ///   - bounds_tight: both bounds came from present SAindex entries
    ///     (safe to skip first matched_level bases in binary search)
    ///
    /// Returns None if no prefix exists in the index (all levels absent).
    pub fn hierarchical_lookup(
        &self,
        kmer_idx: u64,
        k: u32,
        n_sa: usize,
    ) -> Option<(usize, usize, usize, bool)> {
        let mut lind = k.min(self.nbases);
        let mut ind = kmer_idx;

        // If k > nbases, truncate to nbases by removing trailing bases
        if k > self.nbases {
            ind >>= 2 * (k - self.nbases);
        }

        // Walk down levels until we find a present entry
        while lind > 0 {
            let sai_pos = self.genome_sa_index_start[(lind - 1) as usize] + ind;
            if sai_pos >= self.data.len() as u64 {
                lind -= 1;
                ind >>= 2;
                continue;
            }

            let entry = self.data.read(sai_pos as usize);
            let is_absent = (entry >> (self.gstrand_bit + 2)) & 1 != 0;

            if !is_absent {
                // Found present entry — extract SA position
                let sa_pos_mask = (1u64 << self.gstrand_bit) - 1;
                let sa_start = (entry & sa_pos_mask) as usize;

                // Get upper bound from next k-mer at same level
                let level_end = self.genome_sa_index_start[lind as usize];
                let next_pos = self.genome_sa_index_start[(lind - 1) as usize] + ind + 1;

                let (sa_end, bounds_tight) = if next_pos < level_end {
                    let next_entry = self.data.read(next_pos as usize);
                    let next_absent = (next_entry >> (self.gstrand_bit + 2)) & 1 != 0;
                    if !next_absent {
                        ((next_entry & sa_pos_mask) as usize, true)
                    } else {
                        (n_sa, false)
                    }
                } else {
                    (n_sa, false)
                };

                return Some((sa_start, sa_end, lind as usize, bounds_tight));
            }

            lind -= 1;
            ind >>= 2;
        }

        None
    }

    /// Get SA range for a k-mer lookup.
    ///
    /// Returns (start_sa_index, is_present) where is_present indicates
    /// whether this k-mer exists in the genome.
    pub fn lookup(&self, kmer_idx: u64, k: u32) -> (u64, bool) {
        if k == 0 || k > self.nbases {
            return (0, false);
        }

        let sai_pos = self.genome_sa_index_start[(k - 1) as usize] + kmer_idx;
        if sai_pos >= self.data.len() as u64 {
            return (0, false);
        }

        let entry = self.data.read(sai_pos as usize);

        // Check "prefix absent" flag (bit gstrand_bit + 2)
        let is_absent = (entry >> (self.gstrand_bit + 2)) & 1 != 0;

        if is_absent {
            return (0, false);
        }

        // Extract SA position
        let sa_pos = entry & ((1u64 << self.gstrand_bit) - 1);
        (sa_pos, true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_test_index(sequence: &str, bin_nbits: u32, sa_nbases: u32) -> SaIndex {
        let (sai, _) = make_test_index_with_sa(sequence, bin_nbits, sa_nbases);
        sai
    }

    fn make_test_index_with_sa(sequence: &str, bin_nbits: u32, sa_nbases: u32) -> (SaIndex, usize) {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "{sequence}").unwrap();

        let bin_nbits_str = bin_nbits.to_string();
        let sa_nbases_str = sa_nbases.to_string();
        let args = vec![
            "rustar-aligner",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            &bin_nbits_str,
            "--genomeSAindexNbases",
            &sa_nbases_str,
        ];

        let params = Parameters::parse_from(args);
        let genome = Genome::from_fasta(&params).unwrap();
        let sa = SuffixArray::build(&genome).unwrap();
        let sa_len = sa.len();
        let sai = SaIndex::build(&genome, &sa, sa_nbases).unwrap();
        (sai, sa_len)
    }

    #[test]
    fn calculate_num_indices() {
        // For nbases=1: 4^1 = 4
        assert_eq!(SaIndex::calculate_num_indices(1), 4);

        // For nbases=2: 4^1 + 4^2 = 4 + 16 = 20
        assert_eq!(SaIndex::calculate_num_indices(2), 20);

        // For nbases=3: 4 + 16 + 64 = 84
        assert_eq!(SaIndex::calculate_num_indices(3), 84);

        // For nbases=14 (STAR default): ~357 million
        let n14 = SaIndex::calculate_num_indices(14);
        assert_eq!(n14, 357_913_940);
    }

    #[test]
    fn build_simple_index() {
        let sai = make_test_index("ACGT", 2, 2);

        assert_eq!(sai.nbases, 2);
        assert_eq!(sai.genome_sa_index_start, vec![0, 4, 20]);
        assert!(!sai.data.is_empty());
    }

    #[test]
    fn lookup_present_kmer() {
        let sai = make_test_index("AAAA", 2, 2);

        // K-mer "AA" (00 in 2-bit) should be present
        let kmer_idx = 0b00_00; // AA
        let (sa_pos, is_present) = sai.lookup(kmer_idx, 2);

        assert!(is_present, "K-mer 'AA' should be present in 'AAAA'");
        assert!(sa_pos < sai.data.len() as u64);
    }

    #[test]
    fn lookup_absent_kmer() {
        let sai = make_test_index("ACAC", 2, 2);

        // K-mer "GG" (22 in 2-bit) should be absent from "ACAC" (and its revcomp "GTGT")
        let kmer_idx = 0b10_10; // GG
        let (_sa_pos, is_present) = sai.lookup(kmer_idx, 2);

        assert!(!is_present, "K-mer 'GG' should be absent from 'ACAC'");
    }

    #[test]
    fn genome_sa_index_start_progression() {
        let sai = make_test_index("ACGT", 2, 4);

        // Should be [0, 4, 20, 84, 340]
        assert_eq!(sai.genome_sa_index_start[0], 0);
        assert_eq!(sai.genome_sa_index_start[1], 4); // 4^1
        assert_eq!(sai.genome_sa_index_start[2], 20); // 4 + 4^2
        assert_eq!(sai.genome_sa_index_start[3], 84); // 4 + 16 + 4^3
        assert_eq!(sai.genome_sa_index_start[4], 340); // 4 + 16 + 64 + 4^4
    }

    #[test]
    fn hierarchical_lookup_full_kmer_present() {
        // Genome "AAAA" with sa_nbases=2: 2-mer "AA" (idx=0) should be present
        let (sai, n_sa) = make_test_index_with_sa("AAAA", 2, 2);
        let kmer_idx = 0b00_00; // AA

        let result = sai.hierarchical_lookup(kmer_idx, 2, n_sa);
        assert!(result.is_some(), "AA should be found in AAAA");

        let (sa_start, sa_end, matched_level, _bounds_tight) = result.unwrap();
        assert_eq!(matched_level, 2, "Should match at full level 2");
        assert!(sa_start < sa_end, "SA range should be non-empty");
    }

    #[test]
    fn hierarchical_lookup_fallback_to_shorter() {
        // Genome "ACAC" with sa_nbases=2.
        // Forward: A,C,A,C. RC: GTGT → G,T,G,T.
        // 2-mers present: AC(01), CA(10) from forward; GT(23), TG(32) from RC.
        // 2-mer "AT" (idx=0b00_11=3) is absent.
        // But 1-mer "A" (idx=0) is present.
        let (sai, n_sa) = make_test_index_with_sa("ACAC", 2, 2);
        let kmer_idx = 0b00_11; // AT

        let result = sai.hierarchical_lookup(kmer_idx, 2, n_sa);
        assert!(result.is_some(), "Should fall back to 1-mer 'A'");

        let (sa_start, sa_end, matched_level, _bounds_tight) = result.unwrap();
        assert_eq!(matched_level, 1, "Should match at level 1 (1-mer 'A')");
        assert!(sa_start < sa_end, "SA range should be non-empty");
    }

    #[test]
    fn hierarchical_lookup_no_prefix_exists() {
        // Genome "AAAA" → forward: AAAA, RC: TTTT
        // Only 1-mers A(0) and T(3) present. C(1) and G(2) absent at all levels.
        let (sai, n_sa) = make_test_index_with_sa("AAAA", 2, 2);
        let kmer_idx = 0b10; // G (1-mer)

        let result = sai.hierarchical_lookup(kmer_idx, 1, n_sa);
        assert!(result.is_none(), "G should not be found in AAAA genome");
    }

    #[test]
    fn hierarchical_lookup_tight_vs_nontight() {
        // Genome "ACGTACGT" with sa_nbases=2: many 2-mers present
        // AC(01), CG(12), GT(23) and their RC: AC(01), CG(12), GT(23)
        // Also from RC read: CA? TG? etc.
        // "AC" (idx=0b00_01=1) is present. Next 2-mer "AG" (idx=0b00_10=2)?
        // If "AG" is present → tight bounds. If absent → not tight.
        let (sai, n_sa) = make_test_index_with_sa("ACGTACGT", 2, 2);
        let kmer_idx_ac = 0b00_01; // AC

        let result = sai.hierarchical_lookup(kmer_idx_ac, 2, n_sa);
        assert!(result.is_some());

        let (sa_start, sa_end, matched_level, bounds_tight) = result.unwrap();
        assert_eq!(matched_level, 2);
        assert!(sa_start < sa_end);

        // Verify consistency: tight bounds means sa_end came from a present entry
        if bounds_tight {
            // sa_end should be a valid SA index (not n_sa)
            assert!(sa_end <= n_sa);
        }
    }

    #[test]
    fn hierarchical_lookup_matches_lookup_when_present() {
        // When full k-mer is present, hierarchical_lookup should give the same
        // sa_start as lookup()
        let (sai, n_sa) = make_test_index_with_sa("ACGTACGT", 2, 2);

        for kmer_idx in 0..16u64 {
            let (old_sa_pos, old_present) = sai.lookup(kmer_idx, 2);
            let new_result = sai.hierarchical_lookup(kmer_idx, 2, n_sa);

            if old_present {
                let (new_sa_start, _, matched_level, _) = new_result.unwrap();
                assert_eq!(
                    new_sa_start, old_sa_pos as usize,
                    "SA start should match for k-mer {kmer_idx}"
                );
                assert_eq!(
                    matched_level, 2,
                    "Should match at full level for k-mer {kmer_idx}"
                );
            }
        }
    }
}
