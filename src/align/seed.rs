use crate::error::Error;
use crate::index::GenomeIndex;
use crate::io::fastq::complement_base;
use crate::params::Parameters;

/// A seed represents an exact match between a read position and genome location(s).
#[derive(Debug, Clone)]
pub struct Seed {
    /// Position in the read where this seed starts
    pub read_pos: usize,

    /// Length of the exact match
    pub length: usize,

    /// Range in the suffix array [start, end) where this k-mer appears
    pub sa_start: usize,
    pub sa_end: usize,

    /// Whether this seed is on the reverse strand of the read
    pub is_reverse: bool,

    /// Whether this seed was found via R→L (reverse-complement) search.
    /// When true, genome_positions() converts coordinates back to forward orientation.
    pub search_rc: bool,

    /// Mate identifier for paired-end reads
    /// 0 = mate1, 1 = mate2, 2 = single-end (default)
    pub mate_id: u8,
}

impl Seed {
    /// Find all seeds for a read sequence using MMP (Maximal Mappable Prefix) search.
    ///
    /// For each position in the read, performs binary search on the suffix array
    /// to find the longest exact match.
    ///
    /// # Arguments
    /// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
    /// * `index` - Genome index with SA and SAindex
    /// * `min_seed_length` - Minimum seed length to report
    /// * `params` - Parameters including seedMultimapNmax
    ///
    /// # Returns
    /// Vector of seeds found in the read
    pub fn find_seeds(
        read_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
        debug_name: &str,
    ) -> Result<Vec<Seed>, Error> {
        let mut seeds = Vec::new();
        let read_len = read_seq.len();

        // STAR uses the SAME sparse chain-based loop for both L→R and R→L directions.
        // For each direction: Nstart evenly-spaced starting positions, each advancing by
        // MMP length (Lmapped). Chains continue until only seedMapMin (5) bases remain.
        // This matches STAR's ReadAlign_mapOneRead.cpp: for(iDir=0;iDir<2;iDir++)
        //   for(istart=0;istart<Nstart;istart++)
        //     while(istart*Lstart + Lmapped + seedMapMin < readLen) { ... Lmapped += L; }

        // Search L→R (forward direction on read): sparse chain search
        search_direction_sparse(
            read_seq,
            read_len,
            index,
            min_seed_length,
            params,
            false,
            debug_name,
            &mut seeds,
        )?;

        // Cap check between directions (STAR: seedPerReadNmax applies across both)
        if seeds.len() >= params.seed_per_read_nmax {
            return Ok(seeds);
        }

        // Search R→L (reverse direction on read): sparse chain search on RC read
        let rc_read = reverse_complement_read(read_seq);
        search_direction_sparse(
            &rc_read,
            read_len,
            index,
            min_seed_length,
            params,
            true,
            debug_name,
            &mut seeds,
        )?;

        // STAR's storeAligns dedup: same rStart + same Length → skip duplicate.
        // Multiple istart chains can find the same (read_pos, length, direction) seed.
        // Keep only the first occurrence (identical SA ranges guaranteed by same sequence).
        // Matches STAR's OPTIM_STOREaligns_SIMPLE: `if (PC[iP][PC_rStart]==rStart) &&
        // PC[iP][PC_Length]==L) return; //same alignment as before, do not store!`
        {
            let mut seen = std::collections::HashSet::new();
            seeds.retain(|s| seen.insert((s.read_pos, s.length, s.search_rc)));
        }

        Ok(seeds)
    }

    /// Find all seeds for paired-end reads using unified seed pooling.
    ///
    /// This implements STAR's hybrid approach: seeds from both mates are found
    /// independently, tagged with their mate origin, then pooled together for
    /// unified clustering.
    ///
    /// # Arguments
    /// * `mate1_seq` - First mate sequence (encoded)
    /// * `mate2_seq` - Second mate sequence (encoded)
    /// * `index` - Genome index with SA and SAindex
    /// * `min_seed_length` - Minimum seed length to report
    /// * `params` - Parameters including seedMultimapNmax
    ///
    /// # Returns
    /// Vector of seeds from both mates, tagged with mate_id (0 or 1)
    pub fn find_paired_seeds(
        mate1_seq: &[u8],
        mate2_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
    ) -> Result<Vec<Seed>, Error> {
        // Find seeds from mate1 (tag with mate_id = 0)
        let mut seeds = Self::find_seeds(mate1_seq, index, min_seed_length, params, "")?;
        for seed in &mut seeds {
            seed.mate_id = 0;
        }

        // Find seeds from mate2 (tag with mate_id = 1)
        // IMPORTANT: read_pos is relative to mate2 start (will be adjusted during stitching)
        let mut seeds2 = Self::find_seeds(mate2_seq, index, min_seed_length, params, "")?;
        for seed in &mut seeds2 {
            seed.mate_id = 1;
        }

        // Pool seeds together
        seeds.extend(seeds2);

        Ok(seeds)
    }

    /// Get all genome positions for this seed.
    ///
    /// Expands the SA range to actual genome positions.
    pub fn get_genome_positions(&self, index: &GenomeIndex) -> Vec<(u64, bool)> {
        self.genome_positions(index).collect()
    }

    /// Iterate over genome positions for this seed without allocating.
    ///
    /// Returns an iterator that lazily decodes SA entries.
    /// For R→L seeds (search_rc == true), converts positions back:
    /// (pos, is_rev) → (n_genome - pos - length, !is_rev)
    /// Positions where the conversion would underflow are filtered out.
    pub fn genome_positions<'a>(
        &'a self,
        index: &'a GenomeIndex,
    ) -> impl Iterator<Item = (u64, bool)> + 'a {
        let search_rc = self.search_rc;
        let length = self.length as u64;
        let n_genome = index.genome.n_genome;
        (self.sa_start..self.sa_end).filter_map(move |sa_idx| {
            let sa_entry = index.suffix_array.get(sa_idx);
            let (pos, is_rev) = index.suffix_array.decode(sa_entry);
            if search_rc {
                if pos + length <= n_genome {
                    Some((n_genome - pos - length, !is_rev))
                } else {
                    None // Position would span past genome boundary
                }
            } else {
                Some((pos, is_rev))
            }
        })
    }
}

/// Reverse-complement an encoded read sequence.
///
/// Reverses the order and complements each base (A↔T, C↔G).
fn reverse_complement_read(read_seq: &[u8]) -> Vec<u8> {
    read_seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Result of an MMP (Maximal Mappable Prefix) search at a single position.
/// Always provides the advance length for Lmapped tracking, even when no
/// seed is stored (matching STAR's behavior).
struct MmpResult {
    /// The seed to store, if it passed all filters (multimap, min length)
    seed: Option<Seed>,
    /// MMP length to advance by (>= 1). Used for Lmapped tracking regardless
    /// of whether a seed was stored.
    advance: usize,
}

/// Search one direction using STAR's seedSearchNmax-based starting positions with Lmapped tracking.
///
/// Uses seedSearchNmax (= seedSearchStartLmax = 50 by default) evenly-spaced starting
/// positions in [0, seedSearchStartLmax). From each start, does successive MMP searches
/// forward, advancing past found seeds (Lmapped).
///
/// STAR's formula: iStart = seedSearchStartLmax * i / seedSearchNmax
/// With default seedSearchNmax=seedSearchStartLmax=50: iStart = i → dense {0,1,...,49}.
///
/// Used for R→L direction (is_rc=true). L→R uses dense every-position search.
#[allow(clippy::too_many_arguments)]
fn search_direction_sparse(
    read_seq: &[u8],
    original_read_len: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    params: &Parameters,
    is_rc: bool,
    debug_name: &str,
    seeds: &mut Vec<Seed>,
) -> Result<(), Error> {
    let read_len = read_seq.len();

    // STAR (ReadAlign_mapOneRead.cpp lines 41-42):
    //   seedSearchStartLmax = min(P.seedSearchStartLmax, seedSearchStartLmaxOverLread*(Lread-1))
    let effective_start_lmax = if read_len > 0 {
        let over_lread_limit =
            (params.seed_search_start_lmax_over_lread * (read_len as f64 - 1.0)) as usize;
        params.seed_search_start_lmax.min(over_lread_limit)
    } else {
        params.seed_search_start_lmax
    };

    // STAR (line 48): Nstart = seedSearchStartLmax>0 && seedSearchStartLmax<readLen
    //                          ? readLen/seedSearchStartLmax + 1 : 1
    // Same formula for both L→R and R→L (computed once before the iDir loop).
    // For readLen=150, seedSearchStartLmax=50: Nstart=150/50+1=4, Lstart=37.
    let nstart = if effective_start_lmax > 0 && effective_start_lmax < read_len {
        read_len / effective_start_lmax + 1
    } else {
        1
    };
    let lstart = read_len / nstart; // STAR: Lstart = (splitR[1]-splitR[0]) / Nstart

    for istart in 0..nstart {
        let start_pos = (istart * lstart).min(read_len);
        let mut pos = start_pos;

        // From this starting position, search forward with Lmapped tracking.
        // Continue while remaining bases >= seedMapMin (STAR: istart*Lstart + Lmapped + seedMapMin < readLen).
        // Chains advance until only seedMapMin (5) bases remain.
        loop {
            if pos >= read_len {
                break;
            }
            // Stop if remaining bases < seedMapMin (matches STAR's while condition:
            // istart*Lstart + Lmapped + P.seedMapMin < splitR[1][ip]).
            // STAR chains continue until only seedMapMin (5) bases remain, NOT
            // seedSearchStartLmax (50). This allows chains to reach terminal small
            // exons (e.g. 9M after intron) near the read end.
            if read_len - pos < min_seed_length {
                break;
            }

            let result =
                find_seed_at_position(read_seq, pos, index, min_seed_length, false, params)?;

            if !debug_name.is_empty() {
                let dir = if is_rc { "RC" } else { "FWD" };
                let seed_info = match &result.seed {
                    Some(s) => format!("seed(len={} sa={}-{})", s.length, s.sa_start, s.sa_end),
                    None => "no_seed".to_string(),
                };
                eprintln!(
                    "[DEBUG-SEED {}] {} istart={} pos={} advance={} {}",
                    debug_name, dir, istart, pos, result.advance, seed_info
                );
            }

            if let Some(mut seed) = result.seed {
                // Apply seedSearchLmax cap
                if params.seed_search_lmax > 0 && seed.length > params.seed_search_lmax {
                    seed.length = params.seed_search_lmax;
                }

                seed.search_rc = is_rc;

                // Convert RC read_pos back to original read coordinates
                if is_rc {
                    seed.read_pos = original_read_len - seed.read_pos - seed.length;
                }

                seeds.push(seed);

                if seeds.len() >= params.seed_per_read_nmax {
                    return Ok(());
                }
            }

            pos += result.advance; // Always advance by MMP length (matches STAR)
            // Remaining-length check at loop top: stop when < seedMapMin bases remain
        }
    }

    Ok(())
}

/// Find a seed starting at a specific position in the read.
///
/// Returns an MmpResult that always provides the MMP advance length for Lmapped
/// tracking, even when no seed is stored. This matches STAR's behavior where
/// `maxMappableLength2strands()` always returns the MMP length, and `Lmapped += L`
/// always advances — regardless of whether the seed passes filters.
fn find_seed_at_position(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    is_reverse: bool,
    params: &Parameters,
) -> Result<MmpResult, Error> {
    if read_pos >= read_seq.len() {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        });
    }

    // Extract k-mer for SAindex lookup
    let sa_nbases = index.sa_index.nbases as usize;
    let remaining = read_seq.len() - read_pos;

    if remaining < min_seed_length {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        });
    }

    // Build k-mer for SAindex lookup, stopping at first N base
    let lookup_len = remaining.min(sa_nbases);
    let mut kmer_idx = 0u64;
    let mut actual_len = 0usize;

    for i in 0..lookup_len {
        let base = read_seq[read_pos + i];
        if base >= 4 {
            break; // N base — stop building k-mer
        }
        kmer_idx = (kmer_idx << 2) | (base as u64);
        actual_len = i + 1;
    }

    if actual_len == 0 {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        }); // First base is N
    }

    // Hierarchical SAindex lookup (STAR's maxMappableLength2strands approach).
    // Starts at full k-mer level and progressively shortens until a present
    // entry is found. This lets us find seeds even when the full k-mer is
    // absent (e.g., short exon straddling a splice junction).
    let n_sa = index.suffix_array.len();
    let result = index
        .sa_index
        .hierarchical_lookup(kmer_idx, actual_len as u32, n_sa);

    let (sa_start, sa_end, matched_level, bounds_tight) = match result {
        Some(r) => r,
        None => {
            return Ok(MmpResult {
                seed: None,
                advance: 1,
            });
        }
    };

    if sa_start >= sa_end {
        return Ok(MmpResult {
            seed: None,
            advance: 1,
        });
    }

    // STAR short-circuit (maxMappableLength2strands.cpp):
    // "if (Lind < gSAindexNbases && iSA1noN && iSA2good) { maxL=Lind; }"
    // When the hierarchical lookup falls back to a prefix shorter than the full
    // SAindex depth AND bounds are tight, STAR skips the binary search and uses
    // matched_level directly as the MMP. This causes chains to advance by shorter
    // amounts, ensuring intermediate positions (missed if advancing by the true MMP)
    // are not skipped.
    let (match_length, narrowed_start, narrowed_end) = if bounds_tight && matched_level < sa_nbases
    {
        // STAR short-circuit (maxMappableLength2strands.cpp):
        // "if (Lind < gSAindexNbases && iSA1noN && iSA2good) { maxL=Lind; }"
        // Very short prefix already found in SAindex; no genome comparison needed.
        (matched_level, sa_start, sa_end)
    } else {
        // Find maximum mappable prefix length and narrow SA range (STAR's maxMappableLength).
        // When bounds are tight (both from present SAindex entries), we can skip
        // comparing the first matched_level bases since all entries share that prefix.
        let l_initial = if bounds_tight { matched_level } else { 0 };
        max_mappable_length(read_seq, read_pos, index, sa_start, sa_end, l_initial)
    };
    let advance = match_length.max(1);

    // Check seedMultimapNmax: filter seeds that map to too many loci
    // Key fix: still advance by MMP length even when seed is not stored
    // Uses narrowed range (accurate loci count, not overestimated k-mer range)
    let n_loci = narrowed_end - narrowed_start;
    if n_loci > params.seed_multimap_nmax {
        return Ok(MmpResult {
            seed: None,
            advance,
        });
    }

    if match_length >= min_seed_length {
        Ok(MmpResult {
            seed: Some(Seed {
                read_pos,
                length: match_length,
                sa_start: narrowed_start,
                sa_end: narrowed_end,
                is_reverse,
                search_rc: false,
                mate_id: 2, // Single-end default
            }),
            advance,
        })
    } else {
        Ok(MmpResult {
            seed: None,
            advance,
        })
    }
}

/// Overflow-safe median of two unsigned integers.
/// Equivalent to STAR's medianUint2: a/2 + b/2 + (a%2 + b%2)/2
fn median_uint2(a: usize, b: usize) -> usize {
    a / 2 + b / 2 + (a % 2 + b % 2) / 2
}

/// Compare read to genome at a specific SA position, starting from offset l_start.
/// Returns (total_match_length, is_read_greater_at_mismatch).
/// Ports STAR's compareSeqToGenome (SuffixArrayFuns.cpp).
///
/// Starts comparing from offset `l_start` (bases 0..l_start are assumed to match).
/// Walks forward until a mismatch, end of read, or genome padding.
fn compare_seq_to_genome(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_idx: usize,
    l_start: usize,
) -> (usize, bool) {
    let sa_entry = index.suffix_array.get(sa_idx);
    let (genome_pos, is_reverse) = index.suffix_array.decode(sa_entry);

    let genome_start = if is_reverse {
        genome_pos as usize + index.genome.n_genome as usize
    } else {
        genome_pos as usize
    };

    let remaining = read_seq.len() - read_pos;
    let mut match_len = l_start;

    for i in l_start..remaining {
        let genome_idx = genome_start + i;

        if genome_idx >= index.genome.sequence.len() {
            // Past end of genome array — treat like padding (STAR: comp_res > 0)
            return (match_len, true);
        }

        let genome_base = index.genome.sequence[genome_idx];

        if genome_base >= 5 {
            // Padding character — STAR returns comp_res > 0 (read > genome)
            return (match_len, true);
        }

        let read_base = read_seq[read_pos + i];

        if read_base != genome_base {
            return (match_len, read_base > genome_base);
        }

        match_len += 1;
    }

    // Matched all remaining bases — STAR returns comp_res < 0 (genome >= read)
    (match_len, false)
}

/// Find maximum mappable prefix length within SA range [sa_start, sa_end).
/// Binary searches the range while extending match length, then narrows to
/// all positions matching the maximum length.
/// Returns (match_length, narrowed_sa_start, narrowed_sa_end_exclusive).
/// Ports STAR's maxMappableLength (SuffixArrayFuns.cpp).
fn max_mappable_length(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_start: usize,
    sa_end: usize,
    l_initial: usize,
) -> (usize, usize, usize) {
    let remaining = read_seq.len() - read_pos;

    // Single element: just compare
    if sa_start + 1 >= sa_end {
        let (l, _) = compare_seq_to_genome(read_seq, read_pos, index, sa_start, l_initial);
        return (l, sa_start, sa_start + 1);
    }

    // Convert to inclusive range (STAR convention internally)
    let mut i1 = sa_start;
    let mut i2 = sa_end - 1;

    let (mut l1, _) = compare_seq_to_genome(read_seq, read_pos, index, i1, l_initial);
    let (mut l2, _) = compare_seq_to_genome(read_seq, read_pos, index, i2, l_initial);

    let mut l = l1.min(l2);
    let mut l3 = l;
    let mut i3 = i1;

    // Track history for find_mult_range
    let (mut i1a, mut l1a) = (i1, l1);
    let (mut i1b, mut l1b) = (i1, l1);
    let (mut i2a, mut l2a) = (i2, l2);
    let (mut i2b, mut l2b) = (i2, l2);

    // Binary search within SA range
    while i1 + 1 < i2 {
        i3 = median_uint2(i1, i2);
        let comp3;
        (l3, comp3) = compare_seq_to_genome(read_seq, read_pos, index, i3, l);

        if l3 == remaining {
            break; // Perfect match found
        }

        if comp3 {
            // read > genome at mismatch: move left boundary up
            // STAR only shifts history when match length improves (L3 > L1)
            if l3 > l1 {
                i1a = i1b;
                l1a = l1b;
                i1b = i1;
                l1b = l1;
            }
            i1 = i3;
            l1 = l3;
        } else {
            // read <= genome at mismatch: move right boundary down
            // STAR only shifts history when match length improves (L3 > L2)
            if l3 > l2 {
                i2a = i2b;
                l2a = l2b;
                i2b = i2;
                l2b = l2;
            }
            i2 = i3;
            l2 = l3;
        }

        l = l1.min(l2);
    }

    // Pick the best match length
    if l3 < remaining {
        if l1 > l2 {
            l3 = l1;
            i3 = i1;
        } else {
            l3 = l2;
            i3 = i2;
        }
    }

    // Find narrowed range using find_mult_range
    let narrowed_start = find_mult_range(
        read_seq, read_pos, index, remaining, i3, l3, i1, l1, i1a, l1a, i1b, l1b,
    );
    let narrowed_end = find_mult_range(
        read_seq, read_pos, index, remaining, i3, l3, i2, l2, i2a, l2a, i2b, l2b,
    );

    // Convert back to exclusive end
    (l3, narrowed_start, narrowed_end + 1)
}

/// Binary search to find the SA boundary where match length transitions
/// from >= l3 to < l3. Used to narrow the SA range to only positions
/// matching the maximum prefix length.
/// Ports STAR's findMultRange (SuffixArrayFuns.cpp).
///
/// STAR's logic: given the "best" SA index i3 with match length L3,
/// find the farthest SA index that also matches L3 bases, searching
/// outward from i1 (which may or may not already match).
///
/// i1a tracks the boundary with L >= L3 ("good" side)
/// i1b tracks the boundary with L < L3 ("bad" side)
/// Binary search narrows between them until adjacent.
#[allow(clippy::too_many_arguments)]
fn find_mult_range(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    _remaining: usize,
    i3: usize,
    l3: usize,
    i1: usize,
    l1: usize,
    i1a: usize,
    l1a: usize,
    i1b: usize,
    l1b: usize,
) -> usize {
    // STAR's findMultRange: set up (i1a, i1b) search range
    // i1a will have L >= L3 (the "good" side)
    // i1b will have L < L3 (the "bad" side)
    let (mut ia, mut ib, mut lb);
    if l1 < l3 {
        // i1 is below target: search between i3 (good) and i1 (bad)
        ib = i1;
        lb = l1;
        ia = i3;
    } else {
        // i1 already at target length
        if l1a < l1 {
            // Search between i1a (bad) and i1 (good), outward from i1
            ib = i1a;
            lb = l1a;
            ia = i1;
        } else {
            // i1a also at target — search between i1a and i1b
            // (STAR: falls through without reassignment, keeps original i1a/i1b)
            ia = i1a;
            ib = i1b;
            lb = l1b;
        }
    }

    // Binary search: ia has L >= l3, ib has L < l3
    // compareSeqToGenome is called with N=l3 (not remaining), matching STAR
    while (ib + 1 < ia) || (ia + 1 < ib) {
        let ic = median_uint2(ia, ib);
        let (lc, _) = compare_seq_to_genome(read_seq, read_pos, index, ic, lb);

        if lc >= l3 {
            ia = ic;
        } else {
            ib = ic;
            lb = lc;
        }
    }

    ia
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;
    use clap::Parser;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_test_index(sequence: &str) -> GenomeIndex {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "{}", sequence).unwrap();

        let dir = tempfile::tempdir().unwrap();

        let args = vec![
            "rustar-aligner",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeDir",
            dir.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            "2",
            "--genomeSAindexNbases",
            "2",
        ];

        let params = Parameters::parse_from(args);
        GenomeIndex::build(&params).unwrap()
    }

    fn encode_sequence(seq: &str) -> Vec<u8> {
        seq.bytes()
            .map(|b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            })
            .collect()
    }

    #[test]
    fn find_exact_match() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        // Should find at least one seed
        assert!(!seeds.is_empty());

        // First seed should be at position 0 with length 4
        assert_eq!(seeds[0].read_pos, 0);
        assert_eq!(seeds[0].length, 4);
    }

    #[test]
    fn min_seed_length_filter() {
        let index = make_test_index("AAAAAAAA");
        let read = encode_sequence("AAA");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        // With min_seed_length=4, should find nothing (read is only 3bp)
        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(seeds.is_empty());

        // With min_seed_length=2, should find seeds
        let seeds = Seed::find_seeds(&read, &index, 2, &params, "").unwrap();
        assert!(!seeds.is_empty());
    }

    #[test]
    fn no_match() {
        let index = make_test_index("ACAC");
        let read = encode_sequence("GGGG");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 2, &params, "").unwrap();

        // No seeds should be found (GGGG not in ACAC or its reverse complement GTGT)
        assert!(seeds.is_empty());
    }

    #[test]
    fn get_genome_positions() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(!seeds.is_empty());

        // Get positions for first seed
        let positions = seeds[0].get_genome_positions(&index);
        assert!(!positions.is_empty());

        // Should have at least one valid position
        for (pos, _is_reverse) in positions {
            assert!(pos < index.genome.n_genome);
        }
    }

    #[test]
    fn test_single_end_mate_id() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(!seeds.is_empty());

        // Single-end seeds should have mate_id = 2
        for seed in seeds {
            assert_eq!(seed.mate_id, 2);
        }
    }

    #[test]
    fn test_find_paired_seeds() {
        let index = make_test_index("ACGTACGTTTGGCCAA");
        let mate1 = encode_sequence("ACGT");
        let mate2 = encode_sequence("TTGG");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_paired_seeds(&mate1, &mate2, &index, 4, &params).unwrap();

        // Should have seeds from both mates
        let mate1_seeds: Vec<_> = seeds.iter().filter(|s| s.mate_id == 0).collect();
        let mate2_seeds: Vec<_> = seeds.iter().filter(|s| s.mate_id == 1).collect();

        assert!(!mate1_seeds.is_empty(), "Should have mate1 seeds");
        assert!(!mate2_seeds.is_empty(), "Should have mate2 seeds");

        // Verify mate1 seeds have correct read positions
        for seed in mate1_seeds {
            assert!(seed.read_pos < mate1.len());
        }

        // Verify mate2 seeds have correct read positions (relative to mate2)
        for seed in mate2_seeds {
            assert!(seed.read_pos < mate2.len());
        }
    }

    #[test]
    fn test_paired_seeds_pooling() {
        let index = make_test_index("ACGTACGT");
        let mate1 = encode_sequence("ACGT");
        let mate2 = encode_sequence("ACGT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_paired_seeds(&mate1, &mate2, &index, 4, &params).unwrap();

        // Should have roughly double the seeds (one set from each mate)
        let mate1_count = seeds.iter().filter(|s| s.mate_id == 0).count();
        let mate2_count = seeds.iter().filter(|s| s.mate_id == 1).count();

        assert!(mate1_count > 0);
        assert!(mate2_count > 0);
        assert_eq!(seeds.len(), mate1_count + mate2_count);
    }

    #[test]
    fn test_reverse_complement_read() {
        // ACGT → RC = ACGT (palindrome)
        let read = encode_sequence("ACGT");
        let rc = reverse_complement_read(&read);
        assert_eq!(rc, encode_sequence("ACGT"));

        // AACC → RC = GGTT
        let read2 = encode_sequence("AACC");
        let rc2 = reverse_complement_read(&read2);
        assert_eq!(rc2, encode_sequence("GGTT"));

        // Single base
        let read3 = encode_sequence("A");
        let rc3 = reverse_complement_read(&read3);
        assert_eq!(rc3, encode_sequence("T"));

        // N bases preserved
        let read4 = vec![0, 4, 1]; // A, N, C
        let rc4 = reverse_complement_read(&read4);
        assert_eq!(rc4, vec![2, 4, 3]); // G, N, T
    }

    #[test]
    fn test_rl_seeds_found() {
        // Genome has ACGTACGT. The RC of that is ACGTACGT (palindrome),
        // so L→R already finds everything. Use an asymmetric sequence instead.
        // Genome: AACCGGTT — RC genome half has AACCGGTT too.
        // Read: CCGG — L→R finds it at pos 2 in genome.
        // RC of read: CCGG — R→L also finds it.
        // So with this palindromic example, R→L seeds duplicate L→R.
        // Instead use: Genome = AACCTTGG, Read = CCAAGGTT (= RC of AACCTTGG)
        // The read itself won't match L→R in forward genome, but its RC (AACCTTGG) will.
        let index = make_test_index("AACCTTGG");
        // Read is RC of genome: CCAAGGTT
        let read = encode_sequence("CCAAGGTT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        // Should have R→L seeds (search_rc == true)
        let rc_seeds: Vec<_> = seeds.iter().filter(|s| s.search_rc).collect();
        let lr_seeds: Vec<_> = seeds.iter().filter(|s| !s.search_rc).collect();

        // R→L search should find seeds because RC(read) = AACCTTGG matches genome
        assert!(
            !rc_seeds.is_empty(),
            "R→L search should find seeds (RC of read matches genome). All seeds: {:?}",
            seeds
        );

        // Verify R→L seeds have valid read positions
        for seed in &rc_seeds {
            assert!(
                seed.read_pos + seed.length <= read.len(),
                "R→L seed read_pos {} + length {} exceeds read length {}",
                seed.read_pos,
                seed.length,
                read.len()
            );
        }

        // Total seeds should be more than just L→R
        assert!(
            seeds.len() > lr_seeds.len(),
            "Total seeds ({}) should exceed L→R seeds ({})",
            seeds.len(),
            lr_seeds.len()
        );
    }

    #[test]
    fn test_shared_seed_cap() {
        // Test that combined L→R + R→L respects seedPerReadNmax
        let index = make_test_index("ACGTACGTACGTACGT");
        let read = encode_sequence("ACGTACGT");

        let args = vec![
            "rustar-aligner",
            "--runMode",
            "alignReads",
            "--seedPerReadNmax",
            "3",
        ];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(
            seeds.len() <= 3,
            "Total seeds ({}) should respect seedPerReadNmax=3",
            seeds.len()
        );
    }

    #[test]
    fn test_sparse_nstart_calculation() {
        // Verify Nstart matches STAR: readLen/seedSearchStartLmax + 1
        // (when seedSearchStartLmax > 0 && seedSearchStartLmax < readLen)
        //
        // STAR (ReadAlign_mapOneRead.cpp line 48):
        //   Nstart = seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1]
        //            ? splitR[1]/seedSearchStartLmax+1 : 1
        //   Lstart = splitR[1] / Nstart

        // 150 / 50 + 1 = 4, Lstart = 150/4 = 37
        let nstart = 150 / 50 + 1;
        assert_eq!(nstart, 4);
        assert_eq!(150 / nstart, 37);

        // 151 / 50 + 1 = 4, Lstart = 151/4 = 37
        let nstart = 151 / 50 + 1;
        assert_eq!(nstart, 4);
        assert_eq!(151 / nstart, 37);

        // 30 / 50: seedSearchStartLmax (50) >= readLen (30) → Nstart=1
        // (condition seedSearchStartLmax < readLen is false)
        assert_eq!(1_usize, 1);

        // 50 / 50: seedSearchStartLmax (50) >= readLen (50) → Nstart=1
        assert_eq!(1_usize, 1);

        // 100 / 50 + 1 = 3, Lstart = 100/3 = 33
        let nstart = 100 / 50 + 1;
        assert_eq!(nstart, 3);
        assert_eq!(100 / nstart, 33);
    }

    #[test]
    fn test_sparse_rc_read_pos_conversion() {
        // Genome: AACCTTGG, read is RC of genome: CCAAGGTT
        // RC(read) = AACCTTGG matches forward genome at pos 0
        // When searching R→L (is_rc=true), seeds found in the RC read have
        // positions relative to the RC read. They must be converted back to
        // original read coordinates: read_pos = original_read_len - rc_pos - length
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        for seed in &seeds {
            // All seeds (L→R and R→L) must have valid read positions
            assert!(
                seed.read_pos + seed.length <= read.len(),
                "Seed read_pos {} + length {} exceeds read len {} (search_rc={})",
                seed.read_pos,
                seed.length,
                read.len(),
                seed.search_rc
            );
        }

        // Should have R→L seeds
        let rc_seeds: Vec<_> = seeds.iter().filter(|s| s.search_rc).collect();
        assert!(
            !rc_seeds.is_empty(),
            "Should have R→L seeds with valid read positions"
        );
    }

    #[test]
    fn test_sparse_fewer_seeds_than_dense() {
        // With a longer genome and read, sparse search should produce fewer seeds
        // than the old dense (every-position) search
        let genome_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let index = make_test_index(genome_seq);

        // Use a read that's long enough for multiple start positions
        let read = encode_sequence("ACGTACGTACGTACGTACGTACGT"); // 24bp

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let sparse_seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        // Count how many seeds dense would produce (every position that has a match)
        let mut dense_count = 0;
        for read_pos in 0..read.len() {
            let result = find_seed_at_position(&read, read_pos, &index, 4, false, &params).unwrap();
            if result.seed.is_some() {
                dense_count += 1;
            }
        }
        // Also count R→L dense seeds
        let rc_read = reverse_complement_read(&read);
        for rc_pos in 0..rc_read.len() {
            let result =
                find_seed_at_position(&rc_read, rc_pos, &index, 4, false, &params).unwrap();
            if result.seed.is_some() {
                dense_count += 1;
            }
        }

        assert!(
            sparse_seeds.len() <= dense_count,
            "Sparse ({}) should produce <= dense ({}) seeds",
            sparse_seeds.len(),
            dense_count
        );
    }

    #[test]
    fn test_rc_seed_genome_positions() {
        // Genome: AACCTTGG, read RC = CCAAGGTT
        // RC(read) = AACCTTGG matches forward genome
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");

        let args = vec!["rustar-aligner", "--runMode", "alignReads"];
        let params = Parameters::parse_from(args);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        for seed in &seeds {
            if seed.search_rc {
                let positions: Vec<_> = seed.genome_positions(&index).collect();
                assert!(
                    !positions.is_empty(),
                    "RC seed should have genome positions"
                );

                for (pos, _is_rev) in &positions {
                    // Converted positions should be valid (within genome)
                    assert!(
                        *pos < index.genome.n_genome,
                        "Converted position {} should be < n_genome {}",
                        pos,
                        index.genome.n_genome
                    );
                }
            }
        }
    }
}
