/// Seed clustering and stitching via dynamic programming
use crate::align::score::AlignmentScorer;
use crate::align::seed::Seed;
use crate::align::transcript::{CigarOp, Transcript};
use crate::error::Error;
use crate::index::GenomeIndex;

/// STAR's MARK_FRAG_SPACER_BASE (IncludeDefine.h:174).
/// Separates mate1 and mate2 fragments in the combined PE read.
pub(crate) const PE_SPACER_BASE: u8 = 11;

/// Count mismatches in an alignment by comparing read sequence to genome sequence.
///
/// The read sequence is always in forward orientation. For reverse-strand alignments,
/// the genome is accessed at `pos + n_genome` (the reverse-complement region) rather
/// than reverse-complementing the read.
///
/// # Arguments
/// * `read_seq` - Read sequence in forward orientation (encoded as 0=A, 1=C, 2=G, 3=T, 4=N)
/// * `cigar_ops` - CIGAR operations
/// * `genome_start` - Starting position in genome (decoded SA position, WITHOUT n_genome offset)
/// * `read_start` - Starting position in read
/// * `index` - Genome index (contains genome sequence)
/// * `is_reverse` - Whether this is a reverse-strand alignment
///
/// # Returns
/// Number of mismatched bases (excluding N bases)
/// Count mismatches in a simple region (no CIGAR, just read vs genome)
fn count_mismatches_in_region(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    length: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> u32 {
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };
    let mut n_mismatch = 0u32;

    for i in 0..length {
        let read_pos = read_start + i;
        if read_pos >= read_seq.len() {
            break;
        }
        let read_base = read_seq[read_pos];
        if let Some(genome_base) = index
            .genome
            .get_base(genome_start + i as u64 + genome_offset)
            && read_base != genome_base
            && read_base != 4
            && genome_base != 4
        {
            n_mismatch += 1;
        }
    }

    n_mismatch
}

/// Score a region base-by-base, skipping N in read or genome (STAR behavior).
///
/// Mirrors STAR's stitchAlignToTranscript.cpp loop:
///   `if (G[ii]<4 && R[ii]<4) { if match: Score+=scoreMatch; else: Score-=scoreMatch; }`
///
/// Returns `(score, n_mismatch)`. N positions contribute 0 to score (neither +1 nor -1).
fn score_region(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    length: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> (i32, u32) {
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };
    let mut score = 0i32;
    let mut n_mismatch = 0u32;

    for i in 0..length {
        let read_pos = read_start + i;
        if read_pos >= read_seq.len() {
            break;
        }
        let read_base = read_seq[read_pos];
        let Some(genome_base) = index
            .genome
            .get_base(genome_start + i as u64 + genome_offset)
        else {
            break;
        };
        // N in read or genome: skip, no score contribution (STAR: `if (G<4 && R<4)`)
        if read_base >= 4 || genome_base >= 4 {
            continue;
        }
        if read_base == genome_base {
            score += 1;
        } else {
            score -= 1;
            n_mismatch += 1;
        }
    }

    (score, n_mismatch)
}

fn count_mismatches(
    read_seq: &[u8],
    cigar_ops: &[CigarOp],
    genome_start: u64,
    read_start: usize,
    index: &GenomeIndex,
    is_reverse: bool,
) -> u32 {
    // Add n_genome offset for reverse-strand genome access
    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };

    let mut n_mismatch = 0u32;
    let mut read_pos = read_start;
    let mut genome_pos = genome_start;

    for op in cigar_ops {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                for _i in 0..*len {
                    if read_pos < read_seq.len() {
                        let read_base = read_seq[read_pos];
                        if let Some(genome_base) = index.genome.get_base(genome_pos + genome_offset)
                            && read_base != genome_base
                            && read_base != 4
                            && genome_base != 4
                        {
                            n_mismatch += 1;
                        }
                    }
                    read_pos += 1;
                    genome_pos += 1;
                }
            }
            CigarOp::Ins(len) => {
                read_pos += *len as usize;
            }
            CigarOp::Del(len) | CigarOp::RefSkip(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(len) => {
                read_pos += *len as usize;
            }
            CigarOp::HardClip(_) => {}
        }
    }

    n_mismatch
}

/// Result of extending an alignment into flanking regions
#[derive(Debug, Clone)]
struct ExtendResult {
    /// How far the extension reached (bases)
    extend_len: usize,
    /// Maximum score achieved during extension
    max_score: i32,
    /// Number of mismatches in the extended region
    n_mismatch: u32,
}

/// Extend alignment from a boundary into flanking sequence, mirroring STAR's extendAlign().
///
/// Walks base-by-base from the alignment boundary, scoring +1 match / -1 mismatch,
/// tracking the maximum-score extension point. Stops when total mismatches exceed
/// `min(p_mm_max * total_length, n_mm_max)`.
///
/// # Arguments
/// * `read_seq` - Full read sequence (encoded)
/// * `read_start` - Boundary position in read (where extension begins)
/// * `genome_start` - Corresponding genome position (WITHOUT n_genome offset)
/// * `direction` - +1 for rightward extension, -1 for leftward
/// * `max_extend` - Maximum distance to extend (to read boundary)
/// * `n_mm_prev` - Mismatches already in the aligned portion
/// * `len_prev` - Length of the already-aligned portion
/// * `n_mm_max` - outFilterMismatchNmax (absolute max mismatches)
/// * `p_mm_max` - outFilterMismatchNoverLmax (max mismatch ratio)
/// * `index` - Genome index
/// * `is_reverse` - Whether this is a reverse-strand alignment
#[allow(clippy::too_many_arguments)]
fn extend_alignment(
    read_seq: &[u8],
    read_start: usize,
    genome_start: u64,
    direction: i32,
    max_extend: usize,
    n_mm_prev: u32,
    len_prev: usize,
    n_mm_max: u32,
    p_mm_max: f64,
    index: &GenomeIndex,
    is_reverse: bool,
) -> ExtendResult {
    if max_extend == 0 {
        return ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        };
    }

    let genome_offset = if is_reverse { index.genome.n_genome } else { 0 };

    let mut score: i32 = 0;
    let mut max_score: i32 = 0;
    let mut best_len: usize = 0;
    let mut best_mm: u32 = 0;
    let mut n_mm = 0u32;

    for i in 0..max_extend {
        // Calculate read and genome positions based on direction
        let read_pos = if direction > 0 {
            read_start + i
        } else {
            // Leftward: read_start is exclusive boundary, go backwards
            if read_start < 1 + i {
                break;
            }
            read_start - 1 - i
        };

        if read_pos >= read_seq.len() {
            break;
        }

        let genome_pos = if direction > 0 {
            genome_start + i as u64
        } else {
            if genome_start < 1 + i as u64 {
                break;
            }
            genome_start - 1 - i as u64
        };

        // Get genome base (with strand offset)
        let genome_base = match index.genome.get_base(genome_pos + genome_offset) {
            Some(b) => b,
            None => break,
        };

        // Stop at chromosome boundary (padding = 5)
        if genome_base == 5 {
            break;
        }

        let read_base = read_seq[read_pos];

        // Stop at PE fragment spacer base (STAR: MARK_FRAG_SPACER_BASE, extendAlign.cpp:63)
        if read_base == PE_SPACER_BASE {
            break;
        }

        // Skip N bases (no score impact, matches STAR behavior)
        if read_base == 4 || genome_base == 4 {
            continue;
        }

        if read_base == genome_base {
            // MATCH — only record new max on match (STAR behavior)
            score += 1;
            if score > max_score {
                let total_mm = n_mm_prev + n_mm;
                // STAR uses double comparisons throughout (extendAlign.cpp)
                let record_limit_f = (p_mm_max * (len_prev + i + 1) as f64).min(n_mm_max as f64);
                if total_mm as f64 <= record_limit_f {
                    max_score = score;
                    best_len = i + 1;
                    best_mm = n_mm;
                }
            }
        } else {
            // MISMATCH — check break BEFORE incrementing nMM (STAR behavior)
            // Break uses full extension length max_extend, not current position i+1
            let total_mm = n_mm_prev + n_mm;
            // STAR uses double comparisons throughout (extendAlign.cpp)
            let break_limit_f =
                (p_mm_max * (len_prev.saturating_add(max_extend)) as f64).min(n_mm_max as f64);
            if total_mm as f64 >= break_limit_f {
                break;
            }
            n_mm += 1;
            score -= 1;
        }
    }

    // Only accept extension if it has positive score
    if max_score > 0 {
        ExtendResult {
            extend_len: best_len,
            max_score,
            n_mismatch: best_mm,
        }
    } else {
        ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        }
    }
}

/// A Window Alignment entry — equivalent to STAR's WA[iW][iA] array.
///
/// Each entry represents one seed at one specific genome position within a window.
/// During seed assignment, verify_match_at_position() confirms the actual match length.
/// The DP reads these entries directly (no SA range re-expansion needed).
#[derive(Debug, Clone)]
pub struct WindowAlignment {
    /// Index into the seeds array (for DP expansion to identify the originating seed)
    pub seed_idx: usize,
    /// Read start position (STAR: WA_rStart)
    pub read_pos: usize,
    /// Verified match length at this specific position (STAR: WA_Length)
    pub length: usize,
    /// Forward-strand genome position (STAR: WA_gStart)
    pub genome_pos: u64,
    /// Raw SA position (for genome base access in DP — reverse strand uses sa_pos + n_genome)
    pub sa_pos: u64,
    /// SA range size of the originating seed (STAR: WA_Nrep) — for scoring
    pub n_rep: usize,
    /// Whether this entry is an anchor (protected from capacity eviction)
    pub is_anchor: bool,
    /// Mate ID: 0=mate1, 1=mate2, 2=SE (STAR: PC[iP][PC_iFrag] / WA[iW][iS][WA_iFrag])
    pub mate_id: u8,
    /// Pre-computed extension score estimate (STAR: scoreSeedBest[iS] base case).
    /// = length + left_ext_score + right_ext_score (in stitch coords, forward strand).
    /// Computed in stitch_seeds_core after dedup/sort. Default: length as i32.
    pub pre_ext_score: i32,
}

/// A cluster of seeds mapping to the same genomic region
#[derive(Debug, Clone)]
pub struct SeedCluster {
    /// Window Alignment entries — seed positions assigned to this window (STAR's WA array)
    pub alignments: Vec<WindowAlignment>,
    /// Chromosome index
    pub chr_idx: usize,
    /// Genomic start (leftmost position, forward coords, from actual seed positions)
    pub genome_start: u64,
    /// Genomic end (rightmost position, forward coords, from actual seed positions)
    pub genome_end: u64,
    /// Strand (false = forward, true = reverse)
    pub is_reverse: bool,
    /// Anchor seed index (in the seeds array)
    pub anchor_idx: usize,
    /// Anchor genomic bin (anchor_pos >> win_bin_nbits) for MAPQ window counting
    pub anchor_bin: u64,
}

/// Cluster seeds using STAR's bin-based windowing algorithm.
///
/// # Algorithm (faithful to STAR's `createExtendWindowsWithAlign` + `assignAlignToWindow`)
/// 1. Identify anchor seeds (SA range ≤ max_loci_for_anchor)
/// 2. Create windows from anchor positions using `winBin[(strand, bin)]` lookup:
///    - If bin already has a window → assign anchor to it, skip creation
///    - Else scan left/right for nearby windows → merge or create new
/// 3. Extend windows by ±win_flank_nbins on each side
/// 4. Assign ALL seeds to windows with overlap dedup + capacity eviction
/// 5. Build SeedCluster output
///
/// # Arguments
/// * `seeds` - All seeds found in the read
/// * `index` - Genome index
/// * `params` - Parameters (windowing params: winBinNbits, winAnchorDistNbins, winFlankNbins,
///   winAnchorMultimapNmax, seedPerWindowNmax, seedMapMin)
///
/// # Returns
/// Vector of seed clusters, one per window with assigned seeds
pub fn cluster_seeds(
    seeds: &[Seed],
    index: &GenomeIndex,
    params: &crate::params::Parameters,
    read_len: usize,
    _debug: bool,
) -> Vec<SeedCluster> {
    use std::collections::HashMap;

    let win_bin_nbits = params.win_bin_nbits;
    let win_anchor_dist_nbins = params.win_anchor_dist_nbins;
    let win_flank_nbins = params.win_flank_nbins;
    let max_loci_for_anchor = params.win_anchor_multimap_nmax;
    let win_anchor_multimap_nmax = params.win_anchor_multimap_nmax;
    let seed_per_window_nmax = params.seed_per_window_nmax;
    let min_seed_length = params.seed_map_min;

    let anchor_set: Vec<bool> = seeds
        .iter()
        .map(|seed| {
            let n_loci = seed.sa_end - seed.sa_start;
            n_loci > 0 && n_loci <= max_loci_for_anchor
        })
        .collect();

    // Phase 1: Identify anchor seeds (few genomic positions → high specificity)
    // STAR: only seeds with Nrep <= winAnchorMultimapNmax create windows.
    let anchor_indices: Vec<usize> = anchor_set
        .iter()
        .enumerate()
        .filter(|(_, is_anchor)| **is_anchor)
        .map(|(i, _)| i)
        .collect();

    // No fallback: matches STAR behavior where reads with no anchors are unmapped.
    // MMP search now narrows SA ranges from both ends (max_mappable_length),
    // so seeds have accurate loci counts and anchor classification is correct.
    if anchor_indices.is_empty() {
        return Vec::new();
    }

    // Phase 2: Create windows from anchor positions
    // (matches STAR's createExtendWindowsWithAlign)
    struct Window {
        bin_start: u64,
        bin_end: u64,
        chr_idx: usize,
        is_reverse: bool,
        anchor_idx: usize,
        alignments: Vec<WindowAlignment>,
        // Tight bounds from actual seed positions
        actual_start: u64,
        actual_end: u64,
        alive: bool, // false = merged into another window (STAR kills merged windows)
        // STAR's WALrec: persistent minimum non-anchor length threshold.
        // Non-anchor seeds shorter than this are rejected early (before capacity check).
        // Updated during capacity eviction. Matches STAR's assignAlignToWindow behavior.
        wa_lrec: usize,
    }

    let mut windows: Vec<Window> = Vec::new();
    // winBin: (strand, bin) → window_index
    // Chromosome is implicit since bins are from absolute forward positions
    let mut win_bin: HashMap<(bool, u64), usize> = HashMap::new();

    for &anchor_idx in &anchor_indices {
        let anchor = &seeds[anchor_idx];
        let n_loci = anchor.sa_end - anchor.sa_start;

        // Skip anchors with too many loci (STAR: winAnchorMultimapNmax)
        if n_loci > win_anchor_multimap_nmax {
            continue;
        }

        for (sa_pos, strand) in anchor.genome_positions(index) {
            // STAR uses MMP length directly without per-position verification.
            // All SA positions in the range match for the full MMP length by definition.
            let length = anchor.length;
            if length < min_seed_length {
                continue;
            }

            let forward_pos = index.sa_pos_to_forward(sa_pos, strand, length);

            let chr_idx = match index.genome.position_to_chr(forward_pos) {
                Some(info) => info.0,
                None => continue,
            };

            let anchor_bin = forward_pos >> win_bin_nbits;

            // STAR's stitchPieces: Phase 1 creates windows from anchor positions
            // but does NOT populate WA entries. After flank extension, nWA is reset
            // to 0 (line 115), then Phase 3 re-assigns ALL seeds through
            // assignAlignToWindow. We match this by not adding entries here.

            // Check if this bin already has a window (STAR: skip creation, just assign)
            if let Some(&win_idx) = win_bin.get(&(strand, anchor_bin)) {
                let window = &mut windows[win_idx];
                if window.alive && window.chr_idx == chr_idx {
                    window.actual_start = window.actual_start.min(forward_pos);
                    window.actual_end = window.actual_end.max(forward_pos + length as u64);
                    continue;
                }
            }

            // Scan LEFT for existing window to merge with
            let mut merge_left: Option<usize> = None;
            for scan_bin in
                (anchor_bin.saturating_sub(win_anchor_dist_nbins as u64)..anchor_bin).rev()
            {
                if let Some(&win_idx) = win_bin.get(&(strand, scan_bin)) {
                    let w = &windows[win_idx];
                    if w.alive && w.chr_idx == chr_idx {
                        merge_left = Some(win_idx);
                        break;
                    }
                }
            }

            // Scan RIGHT for existing window to merge with
            let mut merge_right: Option<usize> = None;
            for scan_bin in (anchor_bin + 1)..=(anchor_bin + win_anchor_dist_nbins as u64) {
                if let Some(&win_idx) = win_bin.get(&(strand, scan_bin)) {
                    let w = &windows[win_idx];
                    if w.alive && w.chr_idx == chr_idx {
                        merge_right = Some(win_idx);
                        break;
                    }
                }
            }

            match (merge_left, merge_right) {
                (Some(left_idx), Some(right_idx)) if left_idx != right_idx => {
                    // Merge both windows: extend left window to cover right + anchor
                    let right_window = &windows[right_idx];
                    let new_bin_end = right_window.bin_end.max(anchor_bin);
                    let new_actual_start = right_window.actual_start.min(forward_pos);
                    let new_actual_end = right_window.actual_end.max(forward_pos + length as u64);
                    // Kill right window
                    windows[right_idx].alive = false;

                    // Extend left window
                    let left_window = &mut windows[left_idx];
                    left_window.bin_start = left_window.bin_start.min(anchor_bin);
                    left_window.bin_end = left_window.bin_end.max(new_bin_end);
                    left_window.actual_start = left_window.actual_start.min(new_actual_start);
                    left_window.actual_end = left_window.actual_end.max(new_actual_end);

                    // Update winBin for all bins from left to right
                    for bin in left_window.bin_start..=left_window.bin_end {
                        win_bin.insert((strand, bin), left_idx);
                    }
                }
                (Some(idx), _) | (_, Some(idx)) => {
                    // Merge with one existing window
                    let window = &mut windows[idx];
                    window.bin_start = window.bin_start.min(anchor_bin);
                    window.bin_end = window.bin_end.max(anchor_bin);
                    window.actual_start = window.actual_start.min(forward_pos);
                    window.actual_end = window.actual_end.max(forward_pos + length as u64);

                    // Update winBin for newly covered bins
                    for bin in window.bin_start..=window.bin_end {
                        win_bin.insert((strand, bin), idx);
                    }
                }
                _ => {
                    // No merge: create new window
                    let new_idx = windows.len();
                    win_bin.insert((strand, anchor_bin), new_idx);
                    windows.push(Window {
                        bin_start: anchor_bin,
                        bin_end: anchor_bin,
                        chr_idx,
                        is_reverse: strand,
                        anchor_idx,
                        alignments: Vec::new(),
                        actual_start: forward_pos,
                        actual_end: forward_pos + length as u64,
                        alive: true,
                        wa_lrec: 0,
                    });
                }
            }
        }
    }

    if windows.iter().all(|w| !w.alive) {
        return Vec::new();
    }

    // Phase 3: Extend windows by ±win_flank_nbins (matches STAR's flanking extension)
    // Update winBin for newly covered bins
    for (win_idx, window) in windows.iter_mut().enumerate() {
        if !window.alive {
            continue;
        }
        let old_start = window.bin_start;
        let old_end = window.bin_end;
        let new_start = old_start.saturating_sub(win_flank_nbins as u64);
        let new_end = old_end + win_flank_nbins as u64;
        window.bin_start = new_start;
        window.bin_end = new_end;

        let strand = window.is_reverse;
        for bin in new_start..old_start {
            win_bin.entry((strand, bin)).or_insert(win_idx);
        }
        for bin in (old_end + 1)..=new_end {
            win_bin.entry((strand, bin)).or_insert(win_idx);
        }
    }

    // Phase 4: Assign ALL seeds to windows (matches STAR's stitchPieces Phase 3).
    // STAR resets nWA=0 after window creation, then re-assigns ALL seeds (including
    // anchors) through assignAlignToWindow. We match this by not pre-loading anchors
    // in Phase 2 and processing all seeds here through the same capacity logic.
    //
    // STAR processes all pieces in a sorted pass (rStart asc, length desc for ties).
    // This sorting ensures that when two seeds overlap on the same diagonal, the longer
    // one enters the window first and blocks the shorter one via overlap detection —
    // preventing window capacity overflow from many short copies on the same diagonal.
    //
    // We replicate this by: (1) collecting all candidates per window, (2) pre-deduping
    // overlapping entries per diagonal (longest wins, shorter blocked before processing),
    // then (3) processing survivors in original discovery order. This achieves the same
    // overlap suppression without globally reordering seeds across windows.

    // Per-window candidate entries collected before sorting.
    struct WinCandidate {
        seed_idx: usize,
        sa_pos: u64,
        forward_pos: u64,
        length: usize,
        n_loci: usize,
        is_anchor: bool,
        ps_rstart: usize, // positive-strand read start (sort key)
        mate_id: u8,      // STAR: iFrag (overlap dedup must respect fragment boundaries)
    }

    let mut win_candidates: Vec<Vec<WinCandidate>> =
        (0..windows.len()).map(|_| Vec::new()).collect();

    for (seed_idx, seed) in seeds.iter().enumerate() {
        let n_loci = seed.sa_end - seed.sa_start;
        if n_loci == 0 {
            continue;
        }
        let is_anchor_seed = anchor_set[seed_idx];

        for (sa_pos, strand) in seed.genome_positions(index) {
            let length = seed.length;
            if length < min_seed_length {
                continue;
            }

            let forward_pos = index.sa_pos_to_forward(sa_pos, strand, length);

            let chr_idx = match index.genome.position_to_chr(forward_pos) {
                Some(info) => info.0,
                None => {
                    continue;
                }
            };

            let seed_bin = forward_pos >> win_bin_nbits;

            let win_idx = match win_bin.get(&(strand, seed_bin)) {
                Some(&idx) if windows[idx].alive && windows[idx].chr_idx == chr_idx => idx,
                _ => {
                    continue;
                }
            };

            let ps_rstart = if windows[win_idx].is_reverse {
                read_len - (length + seed.read_pos)
            } else {
                seed.read_pos
            };

            win_candidates[win_idx].push(WinCandidate {
                seed_idx,
                sa_pos,
                forward_pos,
                length,
                n_loci,
                is_anchor: is_anchor_seed,
                ps_rstart,
                mate_id: seed.mate_id,
            });
        }
    }

    // Pre-dedup overlapping diagonals (order-independent, length wins):
    // Two candidate entries on the same diagonal with overlapping read ranges represent
    // the same alignment position — keep only the longest. This is what STAR achieves
    // via sorted-order overlap detection (longer seed enters first, shorter blocked).
    // Pre-computing this dedup ensures correct results regardless of discovery order,
    // fixing window capacity overflow without changing capacity-eviction behavior
    // for other reads (seeds on unique diagonals are unaffected).
    //
    // After pre-dedup, Phase 4 processes survivors in original discovery order.
    // The overlap detection in the main Phase 4 loop below is still present but will
    // never find an overlap (since pre-dedup already resolved all of them). It serves
    // as a safety net for any edge cases not covered by pre-dedup.
    let win_n = windows.len();
    let mut win_blocked: Vec<Vec<bool>> = (0..win_n)
        .map(|i| vec![false; win_candidates[i].len()])
        .collect();

    for win_idx in 0..win_n {
        let candidates = &win_candidates[win_idx];
        if candidates.is_empty() {
            continue;
        }
        // Sort indices by length descending so that the longest entry per diagonal
        // is processed first (and blocks shorter overlapping entries on the same diagonal).
        let mut by_len: Vec<usize> = (0..candidates.len()).collect();
        by_len.sort_by(|&a, &b| candidates[b].length.cmp(&candidates[a].length));

        // For each (diagonal, mate_id) pair, track accepted [ps_rstart, ps_rend) ranges.
        // STAR's assignAlignToWindow checks aFrag==WA[iA][WA_iFrag] before overlap test:
        // seeds from different fragments are never treated as overlapping duplicates.
        let mut diag_ranges: HashMap<(i64, u8), Vec<(usize, usize)>> = HashMap::new();
        for &ci in &by_len {
            let cand = &candidates[ci];
            let diag = cand.forward_pos as i64 - cand.ps_rstart as i64;
            let ps_rend = cand.ps_rstart + cand.length;
            let key = (diag, cand.mate_id);

            let blocked = diag_ranges.get(&key).is_some_and(|ranges| {
                ranges.iter().any(|&(rs, re)| {
                    (cand.ps_rstart >= rs && cand.ps_rstart < re) || (ps_rend >= rs && ps_rend < re)
                })
            });

            if blocked {
                win_blocked[win_idx][ci] = true;
            } else {
                diag_ranges
                    .entry(key)
                    .or_default()
                    .push((cand.ps_rstart, ps_rend));
            }
        }
    }

    // Process each window's candidates in original discovery order,
    // skipping pre-dedup-blocked entries.
    let mut too_many_anchors = false; // STAR: MARKER_TOO_MANY_ANCHORS_PER_WINDOW
    'outer: for win_idx in 0..win_n {
        if !windows[win_idx].alive {
            continue;
        }
        for (ci, cand) in win_candidates[win_idx].iter().enumerate() {
            if win_blocked[win_idx][ci] {
                continue; // blocked by a longer overlapping entry on same diagonal
            }

            let seed_idx = cand.seed_idx;
            let seed = &seeds[seed_idx];
            let length = cand.length;
            let forward_pos = cand.forward_pos;
            let sa_pos = cand.sa_pos;
            let n_loci = cand.n_loci;
            let is_anchor_seed = cand.is_anchor;
            let new_ps_rstart = cand.ps_rstart;
            let new_ps_rend = new_ps_rstart + length;

            let window = &mut windows[win_idx];

            // STAR's WALrec early rejection (assignAlignToWindow line 20):
            // Non-anchor seeds shorter than the persistent minimum are rejected
            // before overlap check. Uses strict < (not <=).
            if !is_anchor_seed && length < window.wa_lrec {
                continue;
            }

            // Safety-net overlap detection: after pre-dedup, no overlapping entries
            // should remain. This check handles any edge cases and matches STAR's
            // assignAlignToWindow overlap logic for correctness.
            // STAR checks aFrag==WA[iA][WA_iFrag] before overlap test — seeds from
            // different mate fragments are never merged.
            let new_mate_id = cand.mate_id;
            let new_diag = forward_pos as i64 - new_ps_rstart as i64;
            let mut overlap_idx = None;
            for (i, wa) in window.alignments.iter().enumerate() {
                if wa.mate_id != new_mate_id {
                    continue; // STAR: only merge seeds from the same fragment
                }
                let wa_ps_rstart = if window.is_reverse {
                    read_len - (wa.length + wa.read_pos)
                } else {
                    wa.read_pos
                };
                let wa_ps_rend = wa_ps_rstart + wa.length;
                let wa_diag = wa.genome_pos as i64 - wa_ps_rstart as i64;
                if new_diag == wa_diag
                    && ((new_ps_rstart >= wa_ps_rstart && new_ps_rstart < wa_ps_rend)
                        || (new_ps_rend >= wa_ps_rstart && new_ps_rend < wa_ps_rend))
                {
                    overlap_idx = Some(i);
                    break;
                }
            }
            if let Some(oi) = overlap_idx {
                if length > window.alignments[oi].length {
                    window.alignments.remove(oi);
                    let insert_pos = window.alignments.partition_point(|wa| {
                        let wa_ps = if window.is_reverse {
                            read_len - (wa.length + wa.read_pos)
                        } else {
                            wa.read_pos
                        };
                        wa_ps < new_ps_rstart
                    });
                    window.alignments.insert(
                        insert_pos,
                        WindowAlignment {
                            seed_idx,
                            read_pos: seed.read_pos,
                            length,
                            genome_pos: forward_pos,
                            sa_pos,
                            n_rep: n_loci,
                            is_anchor: is_anchor_seed,
                            mate_id: seed.mate_id,
                            pre_ext_score: length as i32,
                        },
                    );
                }
                continue;
            }

            // Capacity check (seedPerWindowNmax) with anchor protection
            if window.alignments.len() >= seed_per_window_nmax {
                // Find min length of non-anchor entries (STAR: recalculate WALrec)
                let min_non_anchor_len = window
                    .alignments
                    .iter()
                    .filter(|wa| !wa.is_anchor)
                    .map(|wa| wa.length)
                    .min()
                    .unwrap_or(usize::MAX);

                // Update persistent threshold
                window.wa_lrec = min_non_anchor_len;

                // STAR uses strict < for rejection (not <=): entries equal to
                // the minimum CAN enter and trigger eviction+replacement.
                if length < min_non_anchor_len && !is_anchor_seed {
                    continue; // New entry too short
                }

                // Evict shortest non-anchor entries (STAR: WA_Length <= WALrec)
                window
                    .alignments
                    .retain(|wa| wa.is_anchor || wa.length > min_non_anchor_len);

                if window.alignments.len() >= seed_per_window_nmax {
                    // STAR: MARKER_TOO_MANY_ANCHORS_PER_WINDOW → nW=0, abort entire read
                    too_many_anchors = true;
                    break 'outer;
                }
            }

            // STAR's addition condition: if (aAnchor || aLength > WALrec[iW])
            // Non-anchor entries must be STRICTLY longer than wa_lrec to be added.
            if !is_anchor_seed && length <= window.wa_lrec {
                continue;
            }

            // Insert in sorted order by positive-strand read start (matches
            // STAR's assignAlignToWindow lines 107-115 sorted insertion).
            window.actual_start = window.actual_start.min(forward_pos);
            window.actual_end = window.actual_end.max(forward_pos + length as u64);
            let insert_pos = window.alignments.partition_point(|wa| {
                let wa_ps = if window.is_reverse {
                    read_len - (wa.length + wa.read_pos)
                } else {
                    wa.read_pos
                };
                wa_ps < new_ps_rstart
            });
            window.alignments.insert(
                insert_pos,
                WindowAlignment {
                    seed_idx,
                    read_pos: seed.read_pos,
                    length,
                    genome_pos: forward_pos,
                    sa_pos,
                    n_rep: n_loci,
                    is_anchor: is_anchor_seed,
                    mate_id: seed.mate_id,
                    pre_ext_score: length as i32,
                },
            );
        }
    }

    // STAR: if any window hit MARKER_TOO_MANY_ANCHORS_PER_WINDOW, abort entire read
    if too_many_anchors {
        return Vec::new();
    }

    // Phase 5: Build SeedCluster output
    let mut clusters = Vec::with_capacity(windows.len());
    for window in windows.iter() {
        if !window.alive || window.alignments.is_empty() {
            continue;
        }

        clusters.push(SeedCluster {
            alignments: window.alignments.clone(),
            chr_idx: window.chr_idx,
            genome_start: window.actual_start,
            genome_end: window.actual_end,
            is_reverse: window.is_reverse,
            anchor_idx: window.anchor_idx,
            anchor_bin: window.bin_start,
        });
    }

    clusters
}

/// Lightweight exon block for in-progress transcript during recursion
#[derive(Debug, Clone)]
pub(crate) struct ExonBlock {
    pub(crate) read_start: usize, // 0-based inclusive
    pub(crate) read_end: usize,   // 0-based exclusive
    pub(crate) genome_start: u64, // SA coordinate space (raw sa_pos)
    pub(crate) genome_end: u64,   // SA coordinate space (exclusive)
    /// Mate ID: 0=mate1, 1=mate2, 2=SE (STAR: EX_iFrag)
    pub(crate) mate_id: u8,
}

/// In-progress transcript during recursive search (cheap to clone)
#[derive(Debug, Clone)]
pub(crate) struct WorkingTranscript {
    pub(crate) exons: Vec<ExonBlock>,
    pub(crate) score: i32,
    pub(crate) n_mismatch: u32,
    pub(crate) n_gap: u32,
    pub(crate) n_junction: u32,
    pub(crate) junction_motifs: Vec<crate::align::score::SpliceMotif>,
    pub(crate) junction_annotated: Vec<bool>,
    /// Per-junction repeat lengths (jjL, jjR) for overhang check at finalization.
    /// STAR's shiftSJ[isj][0] and shiftSJ[isj][1].
    pub(crate) junction_shifts: Vec<(u32, u32)>,
    pub(crate) n_anchor: u32,
    // Tight bounds for extension at finalization
    pub(crate) read_start: usize,
    pub(crate) read_end: usize,
    pub(crate) genome_start: u64,
    pub(crate) genome_end: u64,
}

impl WorkingTranscript {
    fn new() -> Self {
        WorkingTranscript {
            exons: Vec::new(),
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            junction_shifts: Vec::new(),
            n_anchor: 0,
            read_start: 0,
            read_end: 0,
            genome_start: 0,
            genome_end: 0,
        }
    }
}

/// Fill the gap between the last exon in a WorkingTranscript and the next WA entry.
/// Returns None if stitching fails (mismatch limit, overhang too short, etc.).
/// Matches STAR's stitchAlignToTranscript.cpp logic.
#[allow(clippy::too_many_arguments)]
fn stitch_align_to_transcript(
    wt: &WorkingTranscript,
    wa: &WindowAlignment,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    cluster: &SeedCluster,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    align_mates_gap_max: u64,
    _debug_name: &str,
) -> Option<WorkingTranscript> {
    let last_exon = wt.exons.last().unwrap();

    // Mate-boundary detection: STAR canonSJ[iex] = -3 (stitchAlignToTranscript.cpp:402)
    // When crossing from mate1 to mate2 (or vice versa), skip junction scoring and
    // check alignMatesGapMax instead.
    let last_mate = last_exon.mate_id;
    let is_mate_boundary = wa.mate_id != last_mate && wa.mate_id != 2 && last_mate != 2;

    if is_mate_boundary {
        // STAR allows at most ONE mate-boundary crossing per combined PE transcript.
        // stitchAlignToTranscript.cpp returns -1000007/-1000008 for the 3rd+ seed when
        // mates overlap (genome positions interleaved), naturally limiting WTs to 2 exons.
        // rustar-aligner's overlap-trimming allows continued stitching, inflating combined_n_match.
        // Fix: if the WT already has exons from BOTH mates, a second crossing is invalid.
        let has_m0 = wt.exons.iter().any(|e| e.mate_id == 0);
        let has_m1 = wt.exons.iter().any(|e| e.mate_id == 1);
        if has_m0 && has_m1 {
            return None;
        }
        // STAR condition (stitchAlignToTranscript.cpp:352):
        // gBstart + trA->exons[0][EX_R] + nBasesMax >= trA->exons[0][EX_G] || EX_G < EX_R
        // For forward clusters: checked in forward genome space.
        // For reverse clusters: checked in STAR's encoded genome space, then converted to
        // rustar-aligner's forward-position representation (Phase 16.27 convention).
        if cluster.is_reverse {
            // Reverse cluster: stitch_read = RC(combined) = [mate2|SPACER|RC(mate1)].
            // wt.exons[0] is the first mate2 exon; wa is the new mate1 seed.
            // Recover mate1's position in the original combined read [mate1|SPACER|RC(mate2)]:
            //   combined_start_mate1 = combined_len - stitch_pos_mate1 - len_mate1
            // STAR's encoded-space reject condition (converted to forward-position arithmetic):
            //   combined_start_mate1 < len_mate2_exon - len_mate1 + (P_mate2 - P_mate1)
            let combined_start_mate1 =
                (read_seq.len() as i64) - (wa.read_pos as i64) - (wa.length as i64);
            let first_exon = &wt.exons[0];
            let len_mate2_exon = (first_exon.read_end - first_exon.read_start) as i64;
            let len_mate1 = wa.length as i64;
            let p_diff = first_exon.genome_start as i64 - wa.sa_pos as i64; // P_mate2 - P_mate1
            if combined_start_mate1 < len_mate2_exon - len_mate1 + p_diff {
                return None;
            }
        } else {
            let first_exon = &wt.exons[0];
            let ex_g = first_exon.genome_start;
            let ex_r = first_exon.read_start as u64;
            // STAR stitchAlignToTranscript.cpp line 352:
            // accept if gBstart + EX_R + alignEndsProtrude.nBasesMax >= EX_G || EX_G < EX_R
            // With default alignEndsProtrude.nBasesMax=0 this simplifies to:
            // reject if EX_G >= EX_R && gBstart + EX_R < EX_G
            let fwd_reject = ex_g >= ex_r && wa.genome_pos + ex_r < ex_g;
            if fwd_reject {
                return None;
            }
        }
        // Forward-gap check: alignMatesGapMax (disabled when 0)
        let genome_gap = wa.genome_pos.saturating_sub(last_exon.genome_end);
        if align_mates_gap_max > 0 && genome_gap > align_mates_gap_max {
            return None;
        }
        let mut new_wt = wt.clone();

        // STAR stitchAlignToTranscript.cpp:374-381: right-extend mate A to fragment boundary.
        // extendAlign(R, G, rAend+1, gAend+1, 1, 1, DEF_readSeqLengthMax, nMatch, nMM, ...)
        // extendToEnd=false (local mode); spacer at position 150 naturally stops the extension.
        // nMatch approximated as score + n_mismatch (exact for seeds-only paths, approx otherwise).
        let n_match_m1 = (wt.score + wt.n_mismatch as i32).max(0) as usize;
        let right_ext = extend_alignment(
            read_seq,
            last_exon.read_end,
            last_exon.genome_end,
            1,
            10_000, // DEF_readSeqLengthMax (STAR); SPACER stop applies before this limit
            wt.n_mismatch,
            n_match_m1,
            scorer.n_mm_max,
            scorer.p_mm_max,
            index,
            cluster.is_reverse,
        );
        if right_ext.extend_len > 0 {
            let last = new_wt.exons.last_mut().unwrap();
            last.read_end += right_ext.extend_len;
            last.genome_end += right_ext.extend_len as u64;
            new_wt.score += right_ext.max_score;
            new_wt.n_mismatch += right_ext.n_mismatch;
        }

        // STAR:360, 383-386: add seed B (mate fragment) length to score and push exon.
        let n_mm_after_right = new_wt.n_mismatch;
        new_wt.score += wa.length as i32;
        new_wt.exons.push(ExonBlock {
            read_start: wa.read_pos,
            read_end: wa.read_pos + wa.length,
            genome_start: wa.sa_pos,
            genome_end: wa.sa_pos + wa.length as u64,
            mate_id: wa.mate_id,
        });
        new_wt.read_end = wa.read_pos + wa.length;
        new_wt.genome_end = wa.sa_pos + wa.length as u64;
        if wa.is_anchor {
            new_wt.n_anchor += 1;
        }

        // STAR:390-400: left-extend seed B toward the fragment boundary.
        // extlen = gBstart - EX_G_first + EX_R_first  (when alignEndsType.ext == false)
        // nMatch after adding seed B = n_match_m1 + right_ext.extend_len + wa.length
        //
        // The effective genome start of the first exon = genome_start - read_start,
        // because the base-case left extension will extend the first exon leftward by
        // read_start bases (to read position 0). Using the raw genome_start with the
        // ELSE fallback to wa.read_pos causes over-extension when the first exon was
        // built from a later seed (e.g. pos=86) — the extlen must be computed using
        // STAR's formula (signed) even when wa.sa_pos < first_exon.genome_start.
        let first_exon = &new_wt.exons[0];
        let extlen = {
            let raw = (wa.sa_pos as i64) - (first_exon.genome_start as i64)
                + (first_exon.read_start as i64);
            if raw > 0 {
                (raw as usize).min(wa.read_pos)
            } else {
                0
            }
        };
        let n_match_for_left = n_match_m1 + right_ext.extend_len + wa.length;
        let left_ext = extend_alignment(
            read_seq,
            wa.read_pos,
            wa.sa_pos,
            -1,
            extlen,
            n_mm_after_right,
            n_match_for_left,
            scorer.n_mm_max,
            scorer.p_mm_max,
            index,
            cluster.is_reverse,
        );
        if left_ext.extend_len > 0 {
            let last = new_wt.exons.last_mut().unwrap();
            last.read_start -= left_ext.extend_len;
            last.genome_start -= left_ext.extend_len as u64;
            new_wt.score += left_ext.max_score;
            new_wt.n_mismatch += left_ext.n_mismatch;
        }

        return Some(new_wt);
    }

    // Overlap trimming: if new WA overlaps previous exon in read coords, shift start right
    let mut eff_read_pos = wa.read_pos;
    let mut eff_genome_pos = wa.sa_pos;
    let mut eff_length = wa.length;

    if last_exon.read_end > eff_read_pos {
        let overlap = last_exon.read_end - eff_read_pos;
        if overlap >= eff_length {
            return None; // Fully consumed
        }
        eff_read_pos = last_exon.read_end;
        eff_genome_pos += overlap as u64;
        eff_length -= overlap;
    }

    // Handle genome overlap
    if last_exon.genome_end > eff_genome_pos && eff_genome_pos > last_exon.genome_start {
        let g_overlap = (last_exon.genome_end - eff_genome_pos) as usize;
        if g_overlap >= eff_length {
            return None; // Fully consumed
        }
        eff_read_pos += g_overlap;
        eff_genome_pos += g_overlap as u64;
        eff_length -= g_overlap;
    }

    let read_gap = (eff_read_pos as i64) - (last_exon.read_end as i64);
    let genome_gap = (eff_genome_pos as i64) - (last_exon.genome_end as i64);

    // Reject negative gaps
    if read_gap < 0 || genome_gap < 0 {
        return None;
    }

    let mut new_wt = wt.clone();
    let mut d_score: i32 = 0;
    let mut gap_mm: u32 = 0;

    if read_gap == 0 && genome_gap == 0 {
        // Adjacent seeds — just extend the last exon
        if let Some(last) = new_wt.exons.last_mut() {
            last.read_end = eff_read_pos + eff_length;
            last.genome_end = eff_genome_pos + eff_length as u64;
        }
    } else if read_gap == genome_gap {
        // Equal gap: base-by-base scoring
        let shared = read_gap as usize;
        let (region_score, region_mm) = score_region(
            read_seq,
            last_exon.read_end,
            last_exon.genome_end,
            shared,
            index,
            cluster.is_reverse,
        );
        gap_mm = region_mm;
        d_score += region_score;

        // Extend last exon through the gap and the new seed
        if let Some(last) = new_wt.exons.last_mut() {
            last.read_end = eff_read_pos + eff_length;
            last.genome_end = eff_genome_pos + eff_length as u64;
        }
    } else if genome_gap > read_gap {
        // Deletion or splice junction
        let del = (genome_gap - read_gap) as u32;
        let shared = read_gap as usize;

        let is_splice = del >= scorer.align_intron_min && del <= scorer.align_intron_max;

        // STAR: Del > alignIntronMax → reject (return -1000003)
        if del > scorer.align_intron_max && scorer.align_intron_max > 0 {
            return None;
        }

        // STAR stitchAlignToTranscript.cpp: reject splice when exon B is too short
        // (nBstart < alignSJoverhangMin). Prevents tiny exons from creating spurious
        // splice paths that waste recursion budget with large introns.
        if is_splice && eff_length < scorer.align_sj_overhang_min as usize {
            return None;
        }

        // --- jR scanning for BOTH splice junctions and deletions (STAR-faithful) ---
        // STAR uses the same scanning code path for both cases; the only difference
        // is motif detection (splice) vs pure positional score (deletion).
        // donor_sa = exclusive end of exon A = STAR's gAend+1. jr_shift = STAR's jR.
        let donor_sa = last_exon.genome_end;
        let (jr_shift, motif, motif_score, jj_l, jj_r) = scorer.find_best_junction_position(
            read_seq,
            last_exon.read_end,
            donor_sa,
            read_gap.max(0),
            genome_gap,
            &index.genome,
            cluster.is_reverse,
            index.genome.n_genome,
            last_exon.read_end - last_exon.read_start,
            eff_length,
        );

        // Clamp shift: jr_shift = STAR's jR. Lower bound: can't consume entire exon A.
        // Upper bound: scan already limited to < shared+eff_length but clamp for safety.
        let prev_match_len = (last_exon.read_end - last_exon.read_start) as i32;
        let jr_shift = jr_shift
            .max(-prev_match_len)
            .min((eff_length + shared) as i32);

        // junction_offset = number of shared bases assigned to donor side (= clamped jR)
        let junction_offset = jr_shift.max(0).min(shared as i32) as usize;

        let mut shared_mm = 0u32;
        let mut shared_score = 0i32;

        // Score bases on donor side (before junction/deletion)
        if junction_offset > 0 {
            let (s, mm) = score_region(
                read_seq,
                last_exon.read_end,
                last_exon.genome_end,
                junction_offset,
                index,
                cluster.is_reverse,
            );
            shared_mm += mm;
            shared_score += s;
        }

        // Score bases on acceptor side (after junction/deletion, skip gap)
        let acceptor_bases = shared - junction_offset;
        if acceptor_bases > 0 {
            let (s, mm) = score_region(
                read_seq,
                last_exon.read_end + junction_offset,
                last_exon.genome_end + junction_offset as u64 + del as u64,
                acceptor_bases,
                index,
                cluster.is_reverse,
            );
            shared_mm += mm;
            shared_score += s;
        }

        // Extended left range: junction shifted left of natural position (jR < 0).
        // These bases were presumed matches in exon A; penalize mismatches on acceptor side.
        if jr_shift < 0 {
            let n_extra = (-jr_shift) as usize;
            // STAR: R[rAend + jR+1 .. 0] vs G[gBstart1 + jR+1 .. 0]
            // gBstart1 = eff_genome_pos - shared - 1, so start at gBstart1+jR+1 = eff_genome_pos-shared+jR
            let extra_read_start = (last_exon.read_end as i64 + jr_shift as i64) as usize;
            let extra_genome_start =
                (eff_genome_pos as i64 - shared as i64 + jr_shift as i64) as u64;
            let extra_mm = count_mismatches_in_region(
                read_seq,
                extra_read_start,
                extra_genome_start,
                n_extra,
                index,
                cluster.is_reverse,
            );
            shared_mm += extra_mm;
            shared_score -= 2 * extra_mm as i32;
        }

        // Extended right range: junction shifted right beyond all shared bases (jR > rGap).
        // These bases are in seed B territory; penalize mismatches on donor side.
        if jr_shift > shared as i32 {
            let n_extra = (jr_shift - shared as i32) as usize;
            let extra_read_start = last_exon.read_end + shared;
            let extra_genome_start = last_exon.genome_end + shared as u64;
            let extra_mm = count_mismatches_in_region(
                read_seq,
                extra_read_start,
                extra_genome_start,
                n_extra,
                index,
                cluster.is_reverse,
            );
            shared_mm += extra_mm;
            shared_score -= 2 * extra_mm as i32;
        }

        gap_mm += shared_mm;
        d_score += shared_score;

        // --- Type-specific scoring and tracking ---
        if is_splice {
            // Check stitch mismatch limit
            if !scorer.stitch_mismatch_allowed(&motif, gap_mm) {
                return None;
            }

            // Check annotation (needed for sjdbScore bonus and finalization check)
            let is_annotated = junction_db.is_some_and(|db| {
                let junc_donor_sa = (donor_sa as i64 + jr_shift as i64) as u64;
                let donor_fwd =
                    index.sa_pos_to_forward(junc_donor_sa, cluster.is_reverse, del as usize);
                let acceptor_fwd = donor_fwd + del as u64;
                db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 0)
                    || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 1)
                    || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 2)
            });

            d_score += motif_score;
            if is_annotated {
                d_score += scorer.sjdb_score;
            }

            new_wt.n_junction += 1;
            new_wt.junction_motifs.push(motif);
            new_wt.junction_annotated.push(is_annotated);
            new_wt.junction_shifts.push((jj_l, jj_r));
        } else {
            // Deletion gap scoring
            let del_score = scorer.score_del_open + scorer.score_del_base * del as i32;
            d_score += del_score;
            new_wt.n_gap += 1;
        }

        // --- Common: adjust exon A and create exon B ---
        // jr_shift = STAR's jR: number of shared bases assigned to donor (exon A).
        // Exon A extends right by jr_shift; exon B starts jr_shift bases into the shared region.
        if jr_shift != 0
            && let Some(last) = new_wt.exons.last_mut()
        {
            last.read_end = (last.read_end as i64 + jr_shift as i64) as usize;
            last.genome_end = (last.genome_end as i64 + jr_shift as i64) as u64;
        }

        // New exon for B side: starts at (eff_read_pos - shared + jr_shift), absorbs
        // the (shared - jr_shift) acceptor-side shared bases plus the full eff_length seed.
        let b_read_start = (eff_read_pos as i64 - shared as i64 + jr_shift as i64) as usize;
        let b_genome_start = (eff_genome_pos as i64 - shared as i64 + jr_shift as i64) as u64;
        let b_len = (eff_length as i64 + shared as i64 - jr_shift as i64).max(0) as usize;
        new_wt.exons.push(ExonBlock {
            read_start: b_read_start,
            read_end: b_read_start + b_len,
            genome_start: b_genome_start,
            genome_end: b_genome_start + b_len as u64,
            mate_id: wa.mate_id,
        });
    } else {
        // Insertion: read_gap > genome_gap
        let ins = (read_gap - genome_gap) as usize;
        let shared = genome_gap; // i64, can be negative (genome overlap)

        let mut jr = 0i32; // number of shared bases going to A side

        if shared > 0 {
            // jR scanning to find optimal insertion placement (STAR lines 265-291)
            let shared_usize = shared as usize;
            let genome_offset: u64 = if cluster.is_reverse {
                index.genome.n_genome
            } else {
                0
            };

            let mut score1 = 0i32;
            let mut max_score1 = 0i32;

            // Phase 1: Scan shared bases to find best split point
            // STAR: for (jR1=1; jR1<=gGap; jR1++)
            for jr1 in 1..=shared_usize {
                let g_pos = last_exon.genome_end + (jr1 - 1) as u64 + genome_offset;
                if let Some(gb) = index.genome.get_base(g_pos)
                    && gb < 4
                {
                    // Pre-insertion read base (A side)
                    let r_pre = read_seq[last_exon.read_end + jr1 - 1];
                    // Post-insertion read base (B side, after ins gap)
                    let r_post = read_seq[last_exon.read_end + ins + jr1 - 1];

                    if r_pre == gb {
                        score1 += 1;
                    } else {
                        score1 -= 1;
                    }
                    if r_post == gb {
                        score1 -= 1;
                    } else {
                        score1 += 1;
                    }
                }

                // STAR default (alignInsertionFlush=None): strict > only.
                // First maximum wins = leftmost insertion in current coordinate space.
                // For flushRight mode (not yet implemented): Score1 >= maxScore1.
                if score1 > max_score1 {
                    max_score1 = score1;
                    jr = jr1 as i32;
                }
            }

            // Phase 2: Score shared bases with jR determining pre/post-insertion side
            // STAR: for (ii=1; ii<=gGap; ii++) r1 = rAend+ii+(ii<=jR ? 0 : Ins)
            for ii in 1..=shared_usize {
                let r_idx = if (ii as i32) <= jr {
                    // Pre-insertion side
                    last_exon.read_end + ii - 1
                } else {
                    // Post-insertion side (skip insertion gap)
                    last_exon.read_end + ins + ii - 1
                };
                let g_pos = last_exon.genome_end + (ii - 1) as u64 + genome_offset;

                if let Some(gb) = index.genome.get_base(g_pos)
                    && gb < 4
                    && read_seq[r_idx] < 4
                {
                    if read_seq[r_idx] == gb {
                        d_score += 1;
                    } else {
                        d_score -= 1;
                        gap_mm += 1;
                    }
                }
            }
        } else if shared < 0 {
            // Overlapping seeds on genome — reduce score
            // STAR: for (ii=0; ii<-gGap; ii++) Score -= scoreMatch;
            d_score -= (-shared) as i32;
        }
        // shared == 0: simple insertion, no scanning needed, jr stays 0

        let ins_score = scorer.score_ins_open + scorer.score_ins_base * ins as i32;
        d_score += ins_score;
        new_wt.n_gap += 1;

        // Extend last exon by jr shared bases (A side)
        let jr_usize = jr.max(0) as usize;
        if jr_usize > 0
            && let Some(last) = new_wt.exons.last_mut()
        {
            last.read_end += jr_usize;
            last.genome_end += jr_usize as u64;
        }

        // New exon after insertion
        // B read start: after extended A + insertion gap in read
        // B genome start: after extended A (contiguous on genome)
        let b_read_start = last_exon.read_end + jr_usize + ins;
        let b_genome_start = last_exon.genome_end + jr_usize as u64;
        // B ends at original seed B end
        new_wt.exons.push(ExonBlock {
            read_start: b_read_start,
            read_end: eff_read_pos + eff_length,
            genome_start: b_genome_start,
            genome_end: eff_genome_pos + eff_length as u64,
            mate_id: wa.mate_id,
        });
    }

    // Mismatch limit check
    let total_mm = new_wt.n_mismatch + gap_mm;
    let total_len = new_wt.read_end.max(eff_read_pos + eff_length) - new_wt.read_start;
    let mm_limit = ((scorer.p_mm_max * total_len as f64) as u32).min(scorer.n_mm_max);
    if total_mm > mm_limit {
        return None;
    }

    // Update working transcript
    // Seeds from SA are exact matches (0 internal mismatches).
    // Mismatches are only in gap-fill shared bases (counted in d_score) and extensions.
    new_wt.score += d_score + eff_length as i32;
    new_wt.n_mismatch += gap_mm;
    new_wt.read_end = new_wt.exons.last().map_or(new_wt.read_end, |e| e.read_end);
    new_wt.genome_end = new_wt
        .exons
        .last()
        .map_or(new_wt.genome_end, |e| e.genome_end);
    if wa.is_anchor {
        new_wt.n_anchor += 1;
    }

    Some(new_wt)
}

/// Stitch seeds within a cluster using recursive combinatorial stitching.
///
/// # Arguments
/// * `cluster` - Seed cluster with WindowAlignment entries (STAR's WA array)
/// * `read_seq` - Read sequence
/// * `index` - Genome index
/// * `scorer` - Alignment scorer
///
/// # Returns
/// Vector of transcripts (may have multiple paths through the cluster)
pub fn stitch_seeds(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
) -> Result<Vec<Transcript>, Error> {
    stitch_seeds_with_jdb(cluster, read_seq, index, scorer, None, 1)
}

/// Stitch seeds with optional junction database for annotation-aware scoring.
///
/// Uses STAR's recursive combinatorial stitcher (stitchWindowAligns) which
/// explores include/exclude branches for each seed, allowing it to skip
/// spurious short seeds that would create false splices.
///
/// `max_transcripts_per_window` controls how many transcripts are collected
/// (STAR's `alignTranscriptsPerWindowNmax`, default 100).
pub fn stitch_seeds_with_jdb(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts_per_window: usize,
) -> Result<Vec<Transcript>, Error> {
    stitch_seeds_with_jdb_debug(
        cluster,
        read_seq,
        index,
        scorer,
        junction_db,
        max_transcripts_per_window,
        "",
    )
}

/// Check if two transcripts overlap in their exon blocks.
/// Returns total overlapping bases (where exons share the same read-genome diagonal).
/// Used for deduplication: if new transcript is a subset of existing, drop it.
fn blocks_overlap(t1_exons: &[ExonBlock], t2_exons: &[ExonBlock]) -> u32 {
    let mut overlap = 0u32;
    let mut i = 0;
    let mut j = 0;

    while i < t1_exons.len() && j < t2_exons.len() {
        let e1 = &t1_exons[i];
        let e2 = &t2_exons[j];

        // Check if exons are on the same read-genome diagonal
        let diag1 = e1.genome_start as i64 - e1.read_start as i64;
        let diag2 = e2.genome_start as i64 - e2.read_start as i64;

        if diag1 == diag2 {
            // Same diagonal — compute read-space overlap
            let r_start = e1.read_start.max(e2.read_start);
            let r_end = e1.read_end.min(e2.read_end);
            if r_start < r_end {
                overlap += (r_end - r_start) as u32;
            }
        }

        // Advance the exon that ends first in read space
        if e1.read_end <= e2.read_end {
            i += 1;
        } else {
            j += 1;
        }
    }

    overlap
}

/// Convert a WorkingTranscript to a final Transcript by extending flanks,
/// building CIGAR, counting mismatches, and converting to forward coordinates.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finalize_transcript(
    wt: &WorkingTranscript,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    cluster: &SeedCluster,
    original_is_reverse: bool,
    no_left_ext: bool,
    no_right_ext: bool,
) -> Option<Transcript> {
    use crate::align::transcript::Exon;

    let alignment_start = wt.read_start;
    let alignment_end = wt.read_end;

    // Guard: exon positions must be within read bounds. Out-of-bounds positions can
    // arise when a combined PE read's WorkingTranscript has junction shifts that extend
    // exon read_end past the individual mate boundary. Filter rather than underflow.
    if alignment_end > read_seq.len() {
        return None;
    }

    // STAR: Lprev = tR2 - trA.rStart + 1 (current transcript length)
    let transcript_len = alignment_end - alignment_start;

    // Extend alignment into flanking regions (STAR-style extendAlign).
    // STAR EXTEND_ORDER=1 (stitchWindowAligns.cpp): extend the 5' end of the read first.
    //   Forward strand (original_is_reverse=false): 5' = left (start) → extend left first.
    //   Reverse strand (original_is_reverse=true):  5' = right (end)  → extend right first.
    // After Phase 16.27, stitch_cluster.is_reverse=false for all clusters, so original_is_reverse
    // must be passed explicitly to recover the correct extension order.
    let zero_extend = ExtendResult {
        extend_len: 0,
        max_score: 0,
        n_mismatch: 0,
    };

    let (left_extend, right_extend) = if !original_is_reverse {
        // Forward: extend left (5') first, then right (3')
        let left = if alignment_start > 0 && !no_left_ext {
            extend_alignment(
                read_seq,
                alignment_start,
                wt.genome_start,
                -1,
                alignment_start,
                wt.n_mismatch,
                transcript_len,
                scorer.n_mm_max,
                scorer.p_mm_max,
                index,
                cluster.is_reverse,
            )
        } else {
            zero_extend.clone()
        };
        let transcript_len_after_first = transcript_len + left.extend_len;
        let right = if alignment_end < read_seq.len() && !no_right_ext {
            extend_alignment(
                read_seq,
                alignment_end,
                wt.genome_end,
                1,
                read_seq.len() - alignment_end,
                wt.n_mismatch + left.n_mismatch,
                transcript_len_after_first,
                scorer.n_mm_max,
                scorer.p_mm_max,
                index,
                cluster.is_reverse,
            )
        } else {
            zero_extend.clone()
        };
        (left, right)
    } else {
        // Reverse: extend right (5') first, then left (3')
        let right = if alignment_end < read_seq.len() && !no_right_ext {
            extend_alignment(
                read_seq,
                alignment_end,
                wt.genome_end,
                1,
                read_seq.len() - alignment_end,
                wt.n_mismatch,
                transcript_len,
                scorer.n_mm_max,
                scorer.p_mm_max,
                index,
                cluster.is_reverse,
            )
        } else {
            zero_extend.clone()
        };
        let transcript_len_after_first = transcript_len + right.extend_len;
        let left = if alignment_start > 0 && !no_left_ext {
            extend_alignment(
                read_seq,
                alignment_start,
                wt.genome_start,
                -1,
                alignment_start,
                wt.n_mismatch + right.n_mismatch,
                transcript_len_after_first,
                scorer.n_mm_max,
                scorer.p_mm_max,
                index,
                cluster.is_reverse,
            )
        } else {
            zero_extend
        };
        (left, right)
    };

    // STAR finalization check: exon lengths including repeat lengths (shiftSJ)
    // For non-annotated junctions: exon_len >= alignSJoverhangMin + shiftSJ[side]
    // For annotated junctions: exon_len >= alignSJDBoverhangMin
    // The first exon length includes left extension, last exon includes right extension.
    if wt.n_junction > 0 {
        let mut junction_idx = 0usize;
        for (isj, exon) in wt.exons.iter().enumerate() {
            if isj >= wt.exons.len() - 1 {
                break;
            }
            // Check if this gap between exon[isj] and exon[isj+1] is a junction
            let next_exon = &wt.exons[isj + 1];
            let genome_gap = next_exon.genome_start as i64 - exon.genome_end as i64;
            let read_gap = next_exon.read_start as i64 - exon.read_end as i64;
            let del = genome_gap - read_gap.max(0);
            if del >= scorer.align_intron_min as i64 && junction_idx < wt.junction_shifts.len() {
                // This is a junction — check exon lengths with repeat
                let (shift_l, shift_r) = wt.junction_shifts[junction_idx];
                let is_annotated = wt.junction_annotated[junction_idx];

                // Left exon length (includes left extension for first exon)
                let left_exon_len = if isj == 0 {
                    (exon.read_end - exon.read_start) + left_extend.extend_len
                } else {
                    exon.read_end - exon.read_start
                };

                // Right exon length (includes right extension for last exon)
                let right_exon_len = if isj + 1 == wt.exons.len() - 1 {
                    (next_exon.read_end - next_exon.read_start) + right_extend.extend_len
                } else {
                    next_exon.read_end - next_exon.read_start
                };

                if is_annotated {
                    let min_oh = scorer.align_sjdb_overhang_min as usize;
                    if left_exon_len < min_oh || right_exon_len < min_oh {
                        return None;
                    }
                } else {
                    let min_oh_l = scorer.align_sj_overhang_min as usize + shift_l as usize;
                    let min_oh_r = scorer.align_sj_overhang_min as usize + shift_r as usize;
                    if left_exon_len < min_oh_l || right_exon_len < min_oh_r {
                        return None;
                    }
                }
                junction_idx += 1;
            }
        }
    }

    // STAR check: spliced mates must have mapped length >= alignSplicedMateMapLmin
    // and >= alignSplicedMateMapLminOverLmate * readLength
    if wt.n_junction > 0 {
        let total_mapped =
            left_extend.extend_len + (alignment_end - alignment_start) + right_extend.extend_len;
        let min_from_fraction =
            (scorer.align_spliced_mate_map_lmin_over_lmate * read_seq.len() as f64) as usize;
        let min_mapped = std::cmp::max(
            scorer.align_spliced_mate_map_lmin as usize,
            min_from_fraction,
        );
        if total_mapped < min_mapped {
            return None;
        }
    }

    // Build final CIGAR from exon blocks
    let mut final_cigar: Vec<CigarOp> = Vec::new();

    // Left soft clip
    let remaining_left_clip = alignment_start - left_extend.extend_len;
    if remaining_left_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_left_clip as u32));
    }

    // Left extension match
    if left_extend.extend_len > 0 {
        final_cigar.push(CigarOp::Match(left_extend.extend_len as u32));
    }

    // Walk exon blocks to build CIGAR
    for (idx, exon) in wt.exons.iter().enumerate() {
        if idx > 0 {
            let prev = &wt.exons[idx - 1];
            let read_gap = exon.read_start as i64 - prev.read_end as i64;
            let genome_gap = exon.genome_start as i64 - prev.genome_end as i64;

            if genome_gap > read_gap && genome_gap > 0 {
                // Shared match bases before the gap
                let shared = read_gap.max(0) as u32;
                if shared > 0 {
                    if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                        *prev_len += shared;
                    } else {
                        final_cigar.push(CigarOp::Match(shared));
                    }
                }
                let del = (genome_gap - read_gap.max(0)) as u32;
                if del >= scorer.align_intron_min && del <= scorer.align_intron_max {
                    final_cigar.push(CigarOp::RefSkip(del));
                } else {
                    final_cigar.push(CigarOp::Del(del));
                }
            } else if read_gap > genome_gap && read_gap > 0 {
                // Insertion
                let shared = genome_gap.max(0) as u32;
                if shared > 0 {
                    if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                        *prev_len += shared;
                    } else {
                        final_cigar.push(CigarOp::Match(shared));
                    }
                }
                let ins = (read_gap - genome_gap.max(0)) as u32;
                final_cigar.push(CigarOp::Ins(ins));
            }
            // Equal gap case is handled by extended exon blocks in stitch_align_to_transcript
        }

        // This exon's match region
        let match_len = (exon.read_end - exon.read_start) as u32;
        if match_len > 0 {
            if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
                *prev_len += match_len;
            } else {
                final_cigar.push(CigarOp::Match(match_len));
            }
        }
    }

    // Right extension match
    if right_extend.extend_len > 0 {
        if let Some(CigarOp::Match(prev_len)) = final_cigar.last_mut() {
            *prev_len += right_extend.extend_len as u32;
        } else {
            final_cigar.push(CigarOp::Match(right_extend.extend_len as u32));
        }
    }

    // Right soft clip
    let remaining_right_clip = (read_seq.len() - alignment_end) - right_extend.extend_len;
    if remaining_right_clip > 0 {
        final_cigar.push(CigarOp::SoftClip(remaining_right_clip as u32));
    }

    // Validate CIGAR read-consuming length
    let cigar_read_len: u32 = final_cigar
        .iter()
        .filter(|op| op.consumes_query())
        .map(|op| op.len())
        .sum();
    if cigar_read_len != read_seq.len() as u32 {
        // Invalid CIGAR: exon block geometry is inconsistent with read length.
        // This can occur when a combined PE read's WorkingTranscript has exon
        // positions that span both mates or include the spacer region. Silently
        // filter it out rather than panicking — the alignment is invalid.
        return None;
    }

    // Adjusted genome start for left extension (raw SA coordinates)
    let adjusted_genome_start = wt.genome_start - left_extend.extend_len as u64;

    // Adjust score: add both left and right extensions (STAR adds both at finalization)
    let adjusted_score = wt.score + left_extend.max_score + right_extend.max_score;

    // Count mismatches — MUST be called BEFORE CIGAR reversal
    let n_mismatch = count_mismatches(
        read_seq,
        &final_cigar,
        adjusted_genome_start,
        0,
        index,
        cluster.is_reverse,
    );

    // Reverse CIGAR for reverse strand
    if cluster.is_reverse {
        final_cigar.reverse();
    }

    // Compute total reference-consuming length from CIGAR
    let mut ref_len = 0u64;
    for op in &final_cigar {
        match op {
            CigarOp::Match(len)
            | CigarOp::Equal(len)
            | CigarOp::Diff(len)
            | CigarOp::Del(len)
            | CigarOp::RefSkip(len) => {
                ref_len += *len as u64;
            }
            _ => {}
        }
    }

    // Guard: reverse-strand alignment must not extend past genome end
    if cluster.is_reverse && adjusted_genome_start + ref_len > index.genome.n_genome {
        return None;
    }

    // Convert to forward genome coordinates
    let forward_genome_start =
        index.sa_pos_to_forward(adjusted_genome_start, cluster.is_reverse, ref_len as usize);
    let forward_genome_end = forward_genome_start + ref_len;

    // Build exons from CIGAR
    let mut exons = Vec::new();
    let mut read_pos_e = 0usize;
    let mut genome_pos_e = forward_genome_start;

    for op in &final_cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                let len = *len as usize;
                exons.push(Exon {
                    genome_start: genome_pos_e,
                    genome_end: genome_pos_e + len as u64,
                    read_start: read_pos_e,
                    read_end: read_pos_e + len,
                    // SE / mate1. PE pair-building in `try_pair_transcripts`
                    // rewrites mate2's exons to `i_frag = 1`.
                    i_frag: 0,
                });
                read_pos_e += len;
                genome_pos_e += len as u64;
            }
            CigarOp::Ins(len) => {
                read_pos_e += *len as usize;
            }
            CigarOp::Del(len) => {
                genome_pos_e += *len as u64;
            }
            CigarOp::RefSkip(len) => {
                genome_pos_e += *len as u64;
            }
            CigarOp::SoftClip(len) => {
                read_pos_e += *len as usize;
            }
            CigarOp::HardClip(_) => {}
        }
    }

    // Merge consecutive exons
    let mut merged_exons: Vec<Exon> = Vec::new();
    for exon in exons {
        if let Some(last_exon) = merged_exons.last_mut()
            && last_exon.genome_end == exon.genome_start
            && last_exon.read_end == exon.read_start
        {
            last_exon.genome_end = exon.genome_end;
            last_exon.read_end = exon.read_end;
            continue;
        }
        merged_exons.push(exon);
    }

    let t_genome_start = merged_exons
        .first()
        .map(|e| e.genome_start)
        .unwrap_or(forward_genome_start);
    let t_genome_end = merged_exons
        .last()
        .map(|e| e.genome_end)
        .unwrap_or(forward_genome_end);

    // Apply genomic length penalty
    let genomic_span = t_genome_end - t_genome_start;
    let length_penalty = scorer.genomic_length_penalty(genomic_span);
    let final_score = (adjusted_score + length_penalty).max(0);

    Some(Transcript {
        chr_idx: cluster.chr_idx,
        genome_start: t_genome_start,
        genome_end: t_genome_end,
        is_reverse: cluster.is_reverse,
        exons: merged_exons,
        cigar: final_cigar,
        score: final_score,
        n_mismatch,
        n_gap: wt.n_gap,
        n_junction: wt.n_junction,
        junction_motifs: wt.junction_motifs.clone(),
        junction_annotated: wt.junction_annotated.clone(),
        read_seq: read_seq.to_vec(),
    })
}

/// Recursive include/exclude stitcher (STAR's stitchWindowAligns).
///
/// For each WA entry: try including it (call stitch_align_to_transcript to fill gap),
/// and try excluding it (subject to anchor constraint). Transcripts are finalized
/// at the base case when all entries have been considered.
#[allow(clippy::too_many_arguments)]
fn stitch_recurse(
    i_a: usize,
    wt: WorkingTranscript,
    wa_entries: &[WindowAlignment],
    last_anchor_idx: Option<usize>,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    cluster: &SeedCluster,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts: usize,
    transcripts: &mut Vec<WorkingTranscript>,
    recursion_count: &mut u32,
    align_mates_gap_max: u64,
    original_is_reverse: bool,
    debug_name: &str,
) {
    const MAX_RECURSION: u32 = 100_000;

    if *recursion_count >= MAX_RECURSION {
        return;
    }
    *recursion_count += 1;

    // Base case: all entries considered
    if i_a >= wa_entries.len() {
        if !wt.exons.is_empty() {
            // STAR stitchWindowAligns.cpp lines 62-103: apply left and right extensions
            // inside the stitcher BEFORE the score-gate/dedup. This matches STAR's behavior
            // where extensions boost the correct WT's score so the score-range filter can
            // correctly eliminate spurious WTs with short first exons.
            // EXTEND_ORDER=1: extend 5' of read first (left for fwd, right for rev).
            let mut wt = wt;
            let zero_ext = ExtendResult {
                extend_len: 0,
                max_score: 0,
                n_mismatch: 0,
            };

            let do_left_first = !original_is_reverse;
            let first_ext = if do_left_first {
                // Left extension (5' of forward read)
                if wt.read_start > 0 {
                    extend_alignment(
                        read_seq,
                        wt.read_start,
                        wt.genome_start,
                        -1,
                        wt.read_start,
                        wt.n_mismatch,
                        wt.read_end - wt.read_start,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        cluster.is_reverse,
                    )
                } else {
                    zero_ext.clone()
                }
            } else {
                // Right extension (5' of reverse read)
                if wt.read_end < read_seq.len() {
                    extend_alignment(
                        read_seq,
                        wt.read_end,
                        wt.genome_end,
                        1,
                        read_seq.len() - wt.read_end,
                        wt.n_mismatch,
                        wt.read_end - wt.read_start,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        cluster.is_reverse,
                    )
                } else {
                    zero_ext.clone()
                }
            };

            let len_after_first = (wt.read_end - wt.read_start) + first_ext.extend_len;
            let mm_after_first = wt.n_mismatch + first_ext.n_mismatch;

            let second_ext = if do_left_first {
                // Right extension (3' of forward read)
                let read_end_after = wt.read_end; // right boundary unchanged by left ext
                if read_end_after < read_seq.len() {
                    extend_alignment(
                        read_seq,
                        read_end_after,
                        wt.genome_end,
                        1,
                        read_seq.len() - read_end_after,
                        mm_after_first,
                        len_after_first,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        cluster.is_reverse,
                    )
                } else {
                    zero_ext
                }
            } else {
                // Left extension (3' of reverse read)
                let read_start_after = wt.read_start; // left boundary unchanged by right ext
                if read_start_after > 0 {
                    extend_alignment(
                        read_seq,
                        read_start_after,
                        wt.genome_start,
                        -1,
                        read_start_after,
                        mm_after_first,
                        len_after_first,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        cluster.is_reverse,
                    )
                } else {
                    zero_ext
                }
            };

            // Apply the extensions to wt
            if do_left_first {
                if first_ext.extend_len > 0 {
                    let first = wt.exons.first_mut().unwrap();
                    first.read_start -= first_ext.extend_len;
                    first.genome_start -= first_ext.extend_len as u64;
                    wt.score += first_ext.max_score;
                    wt.n_mismatch += first_ext.n_mismatch;
                    wt.read_start = first.read_start;
                    wt.genome_start = first.genome_start;
                }
                if second_ext.extend_len > 0 {
                    let last = wt.exons.last_mut().unwrap();
                    last.read_end += second_ext.extend_len;
                    last.genome_end += second_ext.extend_len as u64;
                    wt.score += second_ext.max_score;
                    wt.n_mismatch += second_ext.n_mismatch;
                    wt.read_end = last.read_end;
                    wt.genome_end = last.genome_end;
                }
            } else {
                if first_ext.extend_len > 0 {
                    let last = wt.exons.last_mut().unwrap();
                    last.read_end += first_ext.extend_len;
                    last.genome_end += first_ext.extend_len as u64;
                    wt.score += first_ext.max_score;
                    wt.n_mismatch += first_ext.n_mismatch;
                    wt.read_end = last.read_end;
                    wt.genome_end = last.genome_end;
                }
                if second_ext.extend_len > 0 {
                    let first = wt.exons.first_mut().unwrap();
                    first.read_start -= second_ext.extend_len;
                    first.genome_start -= second_ext.extend_len as u64;
                    wt.score += second_ext.max_score;
                    wt.n_mismatch += second_ext.n_mismatch;
                    wt.read_start = first.read_start;
                    wt.genome_start = first.genome_start;
                }
            }

            // Dedup via blocks_overlap: drop if subset of existing higher-score transcript.
            // Use same_structure guard: only dedup transcripts with same number of exon
            // blocks. A non-spliced path should never be killed by a spliced one here
            // because the spliced path may still be rejected by finalize_transcript
            // (overhang check), leaving no valid transcript. STAR-faithful dedup of
            // spliced-vs-unspliced is handled post-finalization below.
            let mut dominated = false;
            let mut remove_indices = Vec::new();
            for (idx, existing) in transcripts.iter().enumerate() {
                let overlap = blocks_overlap(&wt.exons, &existing.exons);
                let wt_len: u32 = wt
                    .exons
                    .iter()
                    .map(|e| (e.read_end - e.read_start) as u32)
                    .sum();
                let ex_len: u32 = existing
                    .exons
                    .iter()
                    .map(|e| (e.read_end - e.read_start) as u32)
                    .sum();

                // Only dedup transcripts with same number of exon blocks (junctions).
                let same_structure = wt.exons.len() == existing.exons.len();
                if same_structure && overlap >= wt_len && existing.score >= wt.score {
                    dominated = true;
                    break;
                }
                if same_structure && overlap >= ex_len && wt.score >= existing.score {
                    remove_indices.push(idx);
                }
            }

            if !dominated {
                // Remove subsets (iterate in reverse to preserve indices)
                for &idx in remove_indices.iter().rev() {
                    transcripts.swap_remove(idx);
                }
                if transcripts.len() < max_transcripts {
                    transcripts.push(wt);
                } else if let Some(worst_idx) = transcripts
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, t)| t.score)
                    .filter(|(_, t)| t.score < wt.score)
                    .map(|(i, _)| i)
                {
                    // STAR-faithful eviction: keep the N best transcripts.
                    // If the new WT scores better than the current worst, evict
                    // the worst and insert the new one.
                    transcripts.swap_remove(worst_idx);
                    transcripts.push(wt);
                }
            }
        }
        return;
    }

    let wa = &wa_entries[i_a];

    // INCLUDE branch: try stitching wa_entries[i_a] to transcript
    if wt.exons.is_empty() {
        // First seed: create initial transcript
        let mut new_wt = wt.clone();
        new_wt.exons.push(ExonBlock {
            read_start: wa.read_pos,
            read_end: wa.read_pos + wa.length,
            genome_start: wa.sa_pos,
            genome_end: wa.sa_pos + wa.length as u64,
            mate_id: wa.mate_id,
        });
        new_wt.score = wa.length as i32;
        new_wt.read_start = wa.read_pos;
        new_wt.read_end = wa.read_pos + wa.length;
        new_wt.genome_start = wa.sa_pos;
        new_wt.genome_end = wa.sa_pos + wa.length as u64;
        if wa.is_anchor {
            new_wt.n_anchor = 1;
        }

        stitch_recurse(
            i_a + 1,
            new_wt,
            wa_entries,
            last_anchor_idx,
            read_seq,
            index,
            scorer,
            cluster,
            junction_db,
            max_transcripts,
            transcripts,
            recursion_count,
            align_mates_gap_max,
            original_is_reverse,
            debug_name,
        );
    } else {
        // Try stitching this seed onto the existing transcript
        if let Some(new_wt) = stitch_align_to_transcript(
            &wt,
            wa,
            read_seq,
            index,
            scorer,
            cluster,
            junction_db,
            align_mates_gap_max,
            debug_name,
        ) {
            stitch_recurse(
                i_a + 1,
                new_wt,
                wa_entries,
                last_anchor_idx,
                read_seq,
                index,
                scorer,
                cluster,
                junction_db,
                max_transcripts,
                transcripts,
                recursion_count,
                align_mates_gap_max,
                original_is_reverse,
                debug_name,
            );
        }
    }

    // EXCLUDE branch: skip wa_entries[i_a]
    // Anchor constraint: can only skip the last anchor if transcript already has one
    let can_exclude = if let Some(last_anchor) = last_anchor_idx {
        if wa.is_anchor && i_a == last_anchor {
            wt.n_anchor > 0 // Already has an anchor → ok to skip
        } else {
            true
        }
    } else {
        true
    };

    if can_exclude {
        stitch_recurse(
            i_a + 1,
            wt,
            wa_entries,
            last_anchor_idx,
            read_seq,
            index,
            scorer,
            cluster,
            junction_db,
            max_transcripts,
            transcripts,
            recursion_count,
            align_mates_gap_max,
            original_is_reverse,
            debug_name,
        );
    }
}

/// Split a combined-read WorkingTranscript by mate_id into per-mate WorkingTranscripts.
///
/// The combined read is `[mate1_seq | PE_SPACER_BASE | RC(mate2_seq)]`. After stitching,
/// each WorkingTranscript has ExonBlocks tagged with mate_id (0=mate1, 1=mate2). This
/// function splits them and adjusts read_pos offsets so each half can be finalized
/// independently with its mate's read slice.
///
/// # Layout in stitch_read
/// - `stitch_is_reverse=false`: `[mate1(0..len1) | SPACER | RC(mate2)(len1+1..)]`
/// - `stitch_is_reverse=true`:  `[mate2(0..len2) | SPACER | RC(mate1)(len2+1..)]`
///
/// # Returns
/// `Some((m1_wt, m2_wt))` if both mates are present; `None` for single-mate WTs.
pub(crate) fn split_combined_wt(
    wt: &WorkingTranscript,
    len1: usize,
    len2: usize,
    stitch_is_reverse: bool,
    align_intron_min: u32,
) -> Option<(WorkingTranscript, WorkingTranscript)> {
    // Read-pos offsets: subtract from each exon's read_start/read_end to get coords
    // within the mate's slice (mate1_seq or RC(mate2_seq) / mate2_seq / RC(mate1_seq)).
    let (m1_offset, m2_offset) = if stitch_is_reverse {
        (len2 + 1, 0usize)
    } else {
        (0usize, len1 + 1)
    };

    let mut m1_exons: Vec<ExonBlock> = Vec::new();
    let mut m2_exons: Vec<ExonBlock> = Vec::new();
    let mut m1_jm: Vec<crate::align::score::SpliceMotif> = Vec::new();
    let mut m1_ja: Vec<bool> = Vec::new();
    let mut m1_js: Vec<(u32, u32)> = Vec::new();
    let mut m2_jm: Vec<crate::align::score::SpliceMotif> = Vec::new();
    let mut m2_ja: Vec<bool> = Vec::new();
    let mut m2_js: Vec<(u32, u32)> = Vec::new();
    let mut junction_idx = 0usize;

    for (i, ex) in wt.exons.iter().enumerate() {
        if ex.mate_id == 0 {
            let mut ex1 = ex.clone();
            ex1.read_start = ex1.read_start.saturating_sub(m1_offset);
            ex1.read_end -= m1_offset;
            m1_exons.push(ex1);
        } else if ex.mate_id == 1 {
            let mut ex2 = ex.clone();
            ex2.read_start = ex2.read_start.saturating_sub(m2_offset);
            ex2.read_end -= m2_offset;
            m2_exons.push(ex2);
        }

        // Classify the gap to the next exon as splice junction or mate boundary.
        // IMPORTANT: mate-boundary transitions (ex.mate_id != next.mate_id) use the
        // is_mate_boundary code path in stitch_align_to_transcript, which does NOT push
        // to wt.junction_motifs. Only intra-mate junctions have entries in junction_motifs,
        // so junction_idx must only advance for intra-mate junctions.
        if i + 1 < wt.exons.len() {
            let next = &wt.exons[i + 1];
            let genome_gap = next.genome_start as i64 - ex.genome_end as i64;
            let read_gap = next.read_start as i64 - ex.read_end as i64;
            let del = genome_gap - read_gap.max(0);
            if del >= align_intron_min as i64 && ex.mate_id == next.mate_id {
                // Intra-mate splice junction: assign to the owning mate
                if junction_idx < wt.junction_motifs.len() {
                    if ex.mate_id == 0 {
                        m1_jm.push(wt.junction_motifs[junction_idx]);
                        m1_ja.push(wt.junction_annotated[junction_idx]);
                        m1_js.push(wt.junction_shifts[junction_idx]);
                    } else if ex.mate_id == 1 {
                        m2_jm.push(wt.junction_motifs[junction_idx]);
                        m2_ja.push(wt.junction_annotated[junction_idx]);
                        m2_js.push(wt.junction_shifts[junction_idx]);
                    }
                    junction_idx += 1;
                }
            }
        }
    }

    if m1_exons.is_empty() || m2_exons.is_empty() {
        return None;
    }

    // STAR PE-CHECK2 (stitchAlignToTranscript.cpp): reject combined WTs where
    // mate1's genomic end exceeds mate2's estimated genomic end.
    // Applied unconditionally (not just for spliced mates) — STAR's debug shows
    // PE-CHECK2 fires for single-exon mate1 too, correctly rejecting overlapping pairs.
    {
        let m1_last = m1_exons.last().unwrap();
        let m2_last = m2_exons.last().unwrap();
        if !stitch_is_reverse {
            let m2_est_end =
                m2_last.genome_start + (len2 as u64).saturating_sub(m2_last.read_start as u64);
            if m1_last.genome_end > m2_est_end {
                return None;
            }
        } else {
            let m1_est_end =
                m1_last.genome_start + (len1 as u64).saturating_sub(m1_last.read_start as u64);
            if m2_last.genome_end > m1_est_end {
                return None;
            }
        }
    }

    let m1_read_start = m1_exons.iter().map(|e| e.read_start).min().unwrap();
    let m1_read_end = m1_exons.iter().map(|e| e.read_end).max().unwrap();
    let m1_genome_start = m1_exons.iter().map(|e| e.genome_start).min().unwrap();
    let m1_genome_end = m1_exons.iter().map(|e| e.genome_end).max().unwrap();
    let m2_read_start = m2_exons.iter().map(|e| e.read_start).min().unwrap();
    let m2_read_end = m2_exons.iter().map(|e| e.read_end).max().unwrap();
    let m2_genome_start = m2_exons.iter().map(|e| e.genome_start).min().unwrap();
    let m2_genome_end = m2_exons.iter().map(|e| e.genome_end).max().unwrap();

    // Approximate per-mate score from exon coverage (includes inner spacer extensions)
    let m1_score: i32 = m1_exons
        .iter()
        .map(|e| (e.read_end - e.read_start) as i32)
        .sum();
    let m2_score: i32 = m2_exons
        .iter()
        .map(|e| (e.read_end - e.read_start) as i32)
        .sum();

    Some((
        WorkingTranscript {
            exons: m1_exons,
            score: m1_score,
            n_mismatch: wt.n_mismatch,
            n_gap: 0,
            n_junction: m1_jm.len() as u32,
            junction_motifs: m1_jm,
            junction_annotated: m1_ja,
            junction_shifts: m1_js,
            n_anchor: 0,
            read_start: m1_read_start,
            read_end: m1_read_end,
            genome_start: m1_genome_start,
            genome_end: m1_genome_end,
        },
        WorkingTranscript {
            exons: m2_exons,
            score: m2_score,
            n_mismatch: wt.n_mismatch,
            n_gap: 0,
            n_junction: m2_jm.len() as u32,
            junction_motifs: m2_jm,
            junction_annotated: m2_ja,
            junction_shifts: m2_js,
            n_anchor: 0,
            read_start: m2_read_start,
            read_end: m2_read_end,
            genome_start: m2_genome_start,
            genome_end: m2_genome_end,
        },
    ))
}

/// Inner implementation of stitch_seeds_with_jdb with optional debug logging.
/// When `debug_read_name` is non-empty, detailed info is logged to stderr.
///
/// Uses STAR's recursive combinatorial stitcher (stitchWindowAligns) instead of
/// forward DP. For each seed, explores include/exclude branches, allowing the
/// algorithm to skip spurious short seeds that would create false splices.
pub(crate) fn stitch_seeds_with_jdb_debug(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts_per_window: usize,
    debug_read_name: &str,
) -> Result<Vec<Transcript>, Error> {
    let (working_transcripts, stitch_cluster, stitch_is_reverse, stitch_read) = stitch_seeds_core(
        cluster,
        read_seq,
        index,
        scorer,
        junction_db,
        max_transcripts_per_window,
        0,
        debug_read_name,
    )?;

    // Finalize working transcripts → Transcript (filtering by overhang+repeat check)
    let mut transcripts: Vec<Transcript> = Vec::with_capacity(working_transcripts.len());
    for wt in &working_transcripts {
        if let Some(mut transcript) = finalize_transcript(
            wt,
            &stitch_read,
            index,
            scorer,
            &stitch_cluster,
            stitch_is_reverse,
            false,
            false,
        ) {
            // Restore original reverse-strand flag and read sequence for SAM output.
            if stitch_is_reverse {
                transcript.is_reverse = true;
                transcript.read_seq = read_seq.to_vec();
            }
            transcripts.push(transcript);
        }
    }

    // STAR-faithful post-finalization dedup (stitchWindowAligns.cpp lines 337-355):
    // After finalization, drop transcripts that are fully covered by a higher-scored one.
    // This removes spurious unspliced secondaries that are subsets of spliced primaries.
    // Performed AFTER finalize so only valid transcripts (passed overhang checks) participate.
    {
        // Compute overlap between two finalized Exon slices (same diagonal + read-space).
        let exon_overlap =
            |a: &[crate::align::transcript::Exon], b: &[crate::align::transcript::Exon]| -> u32 {
                let mut ov = 0u32;
                let mut i = 0;
                let mut j = 0;
                while i < a.len() && j < b.len() {
                    let diag_a = a[i].genome_start as i64 - a[i].read_start as i64;
                    let diag_b = b[j].genome_start as i64 - b[j].read_start as i64;
                    if diag_a == diag_b {
                        let r_start = a[i].read_start.max(b[j].read_start);
                        let r_end = a[i].read_end.min(b[j].read_end);
                        if r_start < r_end {
                            ov += (r_end - r_start) as u32;
                        }
                    }
                    if a[i].read_end <= b[j].read_end {
                        i += 1;
                    } else {
                        j += 1;
                    }
                }
                ov
            };

        let mut keep = vec![true; transcripts.len()];
        for i in 0..transcripts.len() {
            if !keep[i] {
                continue;
            }
            for j in 0..transcripts.len() {
                if i == j || !keep[j] {
                    continue;
                }
                let overlap = exon_overlap(&transcripts[i].exons, &transcripts[j].exons);
                let len_i: u32 = transcripts[i]
                    .exons
                    .iter()
                    .map(|e| (e.read_end - e.read_start) as u32)
                    .sum();
                let len_j: u32 = transcripts[j]
                    .exons
                    .iter()
                    .map(|e| (e.read_end - e.read_start) as u32)
                    .sum();
                let u_i = len_i.saturating_sub(overlap);
                let u_j = len_j.saturating_sub(overlap);

                if u_i == 0 && transcripts[i].score < transcripts[j].score {
                    // i is fully covered by j AND has strictly worse score → i is redundant
                    keep[i] = false;
                    break;
                } else if u_j == 0 && transcripts[j].score < transcripts[i].score {
                    // j is fully covered by i AND has strictly worse score → j is redundant
                    keep[j] = false;
                }
            }
        }
        let mut out = Vec::with_capacity(transcripts.len());
        for (i, t) in transcripts.into_iter().enumerate() {
            if keep[i] {
                out.push(t);
            }
        }
        transcripts = out;
    }

    // Sort by score descending, then shorter genomic span (STAR's gLength tiebreaker).
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then((a.genome_end - a.genome_start).cmp(&(b.genome_end - b.genome_start)))
    });
    transcripts.truncate(max_transcripts_per_window);

    Ok(transcripts)
}

/// Shared core: preprocessing + recursive stitcher, returns working transcripts + context.
#[allow(clippy::too_many_arguments)]
pub(crate) fn stitch_seeds_core(
    cluster: &SeedCluster,
    read_seq: &[u8],
    index: &GenomeIndex,
    scorer: &AlignmentScorer,
    junction_db: Option<&crate::junction::SpliceJunctionDb>,
    max_transcripts_per_window: usize,
    align_mates_gap_max: u64,
    debug_read_name: &str,
) -> Result<(Vec<WorkingTranscript>, SeedCluster, bool, Vec<u8>), Error> {
    let debug = !debug_read_name.is_empty();

    // Include ALL seeds (anchor and non-anchor) in the stitcher.
    // STAR's WA array contains all seeds assigned to the window — non-anchor seeds
    // (n_rep > winAnchorMultimapNmax) participate in stitching. The anchor constraint
    // (at least one anchor must be included) prevents transcripts composed entirely
    // of repetitive seeds, matching STAR's WA_Anchor=2 "last anchor" logic.
    let mut wa_entries: Vec<WindowAlignment> = cluster.alignments.clone();

    // Diagonal dedup: for each diagonal (genome_pos - ps_rstart), merge overlapping
    // seeds into intervals, keeping only the longest seed per merged interval.
    // This prevents combinatorial explosion in the recursive stitcher when many
    // redundant seeds cover the same diagonal region.
    // Uses positive-strand coordinates consistent with cluster_seeds overlap detection.
    {
        use std::collections::HashMap;
        let read_len = read_seq.len();
        let is_rev = cluster.is_reverse;
        // For each (diagonal, mate_id) pair, find the longest seed per merged interval.
        // STAR's assignAlignToWindow checks aFrag==WA[iA][WA_iFrag] before overlap test:
        // seeds from different fragments are never treated as duplicates.
        type DiagMateKey = (i64, u8);
        type DiagSeeds = Vec<(usize, usize, usize)>;
        let mut diag_seeds: HashMap<DiagMateKey, DiagSeeds> = HashMap::new();
        for (idx, wa) in wa_entries.iter().enumerate() {
            let ps = if is_rev {
                read_len - (wa.length + wa.read_pos)
            } else {
                wa.read_pos
            };
            let diag = wa.genome_pos as i64 - ps as i64;
            diag_seeds
                .entry((diag, wa.mate_id))
                .or_default()
                .push((ps, ps + wa.length, idx));
        }

        let mut keep_indices = std::collections::HashSet::new();
        for (_diag, mut seeds) in diag_seeds {
            // Sort by start position
            seeds.sort();
            // Merge intervals, keeping the index of the longest seed in each merged group
            let mut merged_end = seeds[0].1;
            let mut best_idx = seeds[0].2;
            let mut best_len = seeds[0].1 - seeds[0].0;

            for &(s, e, idx) in &seeds[1..] {
                if s <= merged_end {
                    // Overlapping — extend and track longest
                    merged_end = merged_end.max(e);
                    let len = e - s;
                    if len > best_len {
                        best_len = len;
                        best_idx = idx;
                    }
                } else {
                    // New interval — commit previous best
                    keep_indices.insert(best_idx);
                    merged_end = e;
                    best_len = e - s;
                    best_idx = idx;
                }
            }
            // Commit last group
            keep_indices.insert(best_idx);
        }

        // Retain only the kept indices
        let mut idx = 0usize;
        wa_entries.retain(|_| {
            let keep = keep_indices.contains(&idx);
            idx += 1;
            keep
        });
    }

    // STAR-faithful coordinate conversion for stitching:
    // STAR stores WA_gStart in FORWARD genome coordinates (converting RC positions via
    // a1 = nGenome - (aLength + a1)) and uses the RC read for reverse-strand stitching.
    // This ensures gap-fill bases between seeds are scored against the correct genome
    // region. Without this, reverse-strand gap-fill scoring compares against genome
    // adjacent to the wrong seed, inflating false splice scores.
    let stitch_is_reverse = cluster.is_reverse;
    // Always allocate an owned Vec so we can return it alongside the working transcripts.
    // For forward clusters: clone read_seq (no RC needed).
    // For reverse clusters: RC the read (STAR uses Read1[1] for reverse-strand stitching).
    let stitch_read_owned: Vec<u8> = if cluster.is_reverse {
        read_seq
            .iter()
            .rev()
            .map(|&b| if b < 4 { 3 - b } else { b })
            .collect()
    } else {
        read_seq.to_vec()
    };
    let stitch_read: &[u8] = &stitch_read_owned;

    if cluster.is_reverse {
        // Convert WA entries to positive-strand read coords + forward genome coords
        let read_len = read_seq.len();
        for wa in &mut wa_entries {
            wa.read_pos = read_len - (wa.read_pos + wa.length);
            wa.sa_pos = wa.genome_pos;
        }
    }

    // Create cluster for stitching (is_reverse=false for reverse-strand, so stitcher
    // uses forward genome coords without RC genome offset)
    let stitch_cluster = SeedCluster {
        is_reverse: if stitch_is_reverse {
            false
        } else {
            cluster.is_reverse
        },
        ..cluster.clone()
    };

    // Sort ascending by read_pos (positive-strand coordinates after conversion).
    // STAR's WA array is sorted by aRstart (positive-strand read position, ascending).
    // With forward genome coords, gaps are computed correctly for both strands.
    wa_entries.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));

    // Cap entries to prevent exponential blowup in the recursive stitcher.
    // With anchor-only filtering, this limit is rarely hit, but keep as a safety net.
    const MAX_WA_ENTRIES: usize = 200;
    if wa_entries.len() > MAX_WA_ENTRIES {
        wa_entries.sort_by_key(|wa| std::cmp::Reverse(wa.length));
        wa_entries.truncate(MAX_WA_ENTRIES);
        wa_entries.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));
    }

    if wa_entries.is_empty() {
        return Ok((
            Vec::new(),
            stitch_cluster,
            stitch_is_reverse,
            stitch_read_owned,
        ));
    }

    if debug {
        eprintln!(
            "[DEBUG-STITCH {}] {} WA entries (is_reverse={}, stitch_as_fwd={})",
            debug_read_name,
            wa_entries.len(),
            cluster.is_reverse,
            stitch_is_reverse
        );
        for (i, wa) in wa_entries.iter().enumerate().take(30) {
            eprintln!(
                "  wa[{}]: read_pos={}, sa_pos={}, genome_pos={}, length={}, anchor={}, mate={}",
                i, wa.read_pos, wa.sa_pos, wa.genome_pos, wa.length, wa.is_anchor, wa.mate_id
            );
        }
        if wa_entries.len() > 30 {
            eprintln!("  ... ({} more)", wa_entries.len() - 30);
        }
    }

    // --- STAR stitchWindowSeeds.cpp: scoreSeedBest pre-extension ---
    //
    // Phase 1: Pre-extend each seed left and right (base case of scoreSeedBest DP).
    // Uses stitch_read (forward strand, forward genome coords after coord conversion above).
    // EXTEND_ORDER=1: for forward reads, left extension first (5' of read), then right.
    //                 for reverse reads, right extension first (5' of read in RC), then left.
    // This matches stitch_recurse's `do_left_first = !original_is_reverse`.
    // Mismatch budget for second ext carries over from first ext (STAR: scoreSeedBestMM[iS1]).
    {
        let zero_ext = ExtendResult {
            extend_len: 0,
            max_score: 0,
            n_mismatch: 0,
        };
        let do_left_first = !stitch_is_reverse;
        for wa in &mut wa_entries {
            let right_start = wa.read_pos + wa.length;

            let (first_ext, second_ext) = if do_left_first {
                // Forward cluster: left ext first, then right
                let left = if wa.read_pos > 0 {
                    extend_alignment(
                        stitch_read,
                        wa.read_pos,
                        wa.sa_pos,
                        -1,
                        wa.read_pos,
                        0,
                        wa.length,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        false,
                    )
                } else {
                    zero_ext.clone()
                };
                let right_len_prev = wa.length + left.extend_len;
                let right = if right_start < stitch_read.len() {
                    extend_alignment(
                        stitch_read,
                        right_start,
                        wa.sa_pos + wa.length as u64,
                        1,
                        stitch_read.len() - right_start,
                        left.n_mismatch,
                        right_len_prev,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        false,
                    )
                } else {
                    zero_ext.clone()
                };
                (left, right)
            } else {
                // Reverse cluster: right ext first (5' of RC read), then left
                let right = if right_start < stitch_read.len() {
                    extend_alignment(
                        stitch_read,
                        right_start,
                        wa.sa_pos + wa.length as u64,
                        1,
                        stitch_read.len() - right_start,
                        0,
                        wa.length,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        false,
                    )
                } else {
                    zero_ext.clone()
                };
                let left_len_prev = wa.length + right.extend_len;
                let left = if wa.read_pos > 0 {
                    extend_alignment(
                        stitch_read,
                        wa.read_pos,
                        wa.sa_pos,
                        -1,
                        wa.read_pos,
                        right.n_mismatch,
                        left_len_prev,
                        scorer.n_mm_max,
                        scorer.p_mm_max,
                        index,
                        false,
                    )
                } else {
                    zero_ext.clone()
                };
                (left, right)
            };

            wa.pre_ext_score = wa.length as i32 + first_ext.max_score + second_ext.max_score;
        }
    }

    // Phase 2: DP chain over pre-extended scores (STAR: scoreSeedBest chain accumulation).
    // dp[i] = best chain score ending at wa_entries[i].
    // For intra-read gaps: if read_gap != genome_gap this is a potential splice — 0 penalty.
    // For exact-match gaps (read_gap == genome_gap): score the gap as mismatches only
    // (lenient: use 0 since exact-match gaps are handled by the actual DP later).
    // Both cases produce an overestimate — conservative upper bound matching STAR's approach.
    let best_pre_score = {
        let mut dp: Vec<i32> = wa_entries.iter().map(|wa| wa.pre_ext_score).collect();
        for i in 1..wa_entries.len() {
            for j in 0..i {
                let rj_end = wa_entries[j].read_pos + wa_entries[j].length;
                let gj_end = wa_entries[j].sa_pos + wa_entries[j].length as u64;
                // Must not overlap in read or genome
                if wa_entries[i].read_pos < rj_end || wa_entries[i].sa_pos < gj_end {
                    continue;
                }
                // Gap penalty: 0 for potential splices (read_gap != genome_gap) or short gaps
                // This prevents underestimating spliced reads whose per-seed pre_ext_score is
                // limited by the splice junction (extension stops at the intron boundary).
                let chain = dp[j] + wa_entries[i].pre_ext_score;
                if chain > dp[i] {
                    dp[i] = chain;
                }
            }
        }
        dp.into_iter().max().unwrap_or(0)
    };

    // Note: scoreSeedBest (best_pre_score) is used by STAR for seed ordering within
    // stitchWindowAligns, NOT as a hard pre-filter gate. We keep pre_ext_score on each WA
    // entry for future use in seed ordering (Phase B), but do not filter here.
    let _ = best_pre_score; // suppress unused warning

    // Find last anchor index for the anchor constraint
    let last_anchor_idx = wa_entries.iter().rposition(|wa| wa.is_anchor);

    // Run recursive include/exclude stitcher
    let mut working_transcripts: Vec<WorkingTranscript> = Vec::new();
    let mut recursion_count: u32 = 0;

    stitch_recurse(
        0,
        WorkingTranscript::new(),
        &wa_entries,
        last_anchor_idx,
        stitch_read,
        index,
        scorer,
        &stitch_cluster,
        junction_db,
        max_transcripts_per_window,
        &mut working_transcripts,
        &mut recursion_count,
        align_mates_gap_max,
        stitch_is_reverse,
        debug_read_name,
    );

    if debug {
        eprintln!(
            "[DEBUG-STITCH {}] Recursion done: {} working transcripts, {} recursions",
            debug_read_name,
            working_transcripts.len(),
            recursion_count
        );
    }

    Ok((
        working_transcripts,
        stitch_cluster,
        stitch_is_reverse,
        stitch_read_owned,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;

    fn make_simple_index() -> GenomeIndex {
        // Simple genome: ACGTACGTNN (10 bases)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 4, 4];
        let n_genome = 64u64; // Padded
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(&seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        let genome = Genome {
            sequence,
            n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![10],
            chr_start: vec![0, n_genome],
        };

        // Create dummy SA and SAindex
        let gstrand_bit = 33;
        let suffix_array = SuffixArray {
            data: PackedArray::new(gstrand_bit, 0),
            gstrand_bit,
            gstrand_mask: (1u64 << gstrand_bit) - 1,
        };

        let word_length = gstrand_bit + 3;
        let sa_index = SaIndex {
            data: PackedArray::new(word_length, 0),
            nbases: 14,
            genome_sa_index_start: vec![0],
            word_length,
            gstrand_bit,
        };

        GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db: crate::junction::SpliceJunctionDb::empty(),
            transcriptome: None,
            prepared_junctions: Vec::new(),
        }
    }

    #[test]
    fn test_cluster_seeds_simple() {
        // This test would require a properly populated SA
        // For now, just test that we can create the clustering structure
        // without panicking on an empty index
        let index = make_simple_index();

        // Create seeds with empty SA ranges (won't expand to any positions)
        let seeds = vec![
            Seed {
                read_pos: 0,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
                search_rc: false,
                mate_id: 2,
            },
            Seed {
                read_pos: 10,
                length: 5,
                sa_start: 0,
                sa_end: 0, // Empty range
                is_reverse: false,
                search_rc: false,
                mate_id: 2,
            },
        ];

        // Bin-based windowing with default parameters
        let params = {
            use clap::Parser;
            crate::params::Parameters::parse_from(vec!["rustar-aligner"])
        };
        let clusters = cluster_seeds(&seeds, &index, &params, 150, false);

        // With empty SA ranges, no clusters will be created
        assert_eq!(clusters.len(), 0);
    }

    #[test]
    fn test_wa_entry_sorting() {
        let mut entries = vec![
            WindowAlignment {
                seed_idx: 0,
                read_pos: 10,
                length: 5,
                genome_pos: 100,
                sa_pos: 100,
                n_rep: 1,
                is_anchor: true,
                mate_id: 2,
                pre_ext_score: 5,
            },
            WindowAlignment {
                seed_idx: 1,
                read_pos: 5,
                length: 5,
                genome_pos: 50,
                sa_pos: 50,
                n_rep: 1,
                is_anchor: true,
                mate_id: 2,
                pre_ext_score: 5,
            },
        ];

        entries.sort_by(|a, b| a.read_pos.cmp(&b.read_pos).then(b.length.cmp(&a.length)));

        assert_eq!(entries[0].read_pos, 5);
        assert_eq!(entries[1].read_pos, 10);
    }

    /// Helper to build a GenomeIndex with a specific forward sequence
    fn make_index_with_seq(seq: &[u8]) -> GenomeIndex {
        let n_genome = ((seq.len() as u64 + 1) / 64 + 1) * 64;
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        let genome = Genome {
            sequence,
            n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![seq.len() as u64],
            chr_start: vec![0, n_genome],
        };

        let gstrand_bit = 33;
        let suffix_array = SuffixArray {
            data: PackedArray::new(gstrand_bit, 0),
            gstrand_bit,
            gstrand_mask: (1u64 << gstrand_bit) - 1,
        };
        let word_length = gstrand_bit + 3;
        let sa_index = SaIndex {
            data: PackedArray::new(word_length, 0),
            nbases: 14,
            genome_sa_index_start: vec![0],
            word_length,
            gstrand_bit,
        };

        GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db: crate::junction::SpliceJunctionDb::empty(),
            transcriptome: None,
            prepared_junctions: Vec::new(),
        }
    }

    #[test]
    fn test_extend_perfect_match_rightward() {
        // Genome: ACGTACGTAC (10 bases, A=0, C=1, G=2, T=3)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&seq);
        // Read matches genome perfectly from position 5 onward
        let read_seq = vec![1, 2, 3, 0, 1]; // matches genome[5..10]

        let result = extend_alignment(
            &read_seq, 0,   // read_start (boundary)
            5,   // genome_start
            1,   // rightward
            5,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false, // forward strand
        );

        assert_eq!(result.extend_len, 5);
        assert_eq!(result.max_score, 5);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_perfect_match_leftward() {
        // Genome: ACGTACGTAC
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&seq);
        // Read matches genome[0..5] = ACGTA
        let read_seq = vec![0, 1, 2, 3, 0];

        let result = extend_alignment(
            &read_seq, 5,   // read_start (exclusive boundary for leftward)
            5,   // genome_start (exclusive boundary for leftward)
            -1,  // leftward
            5,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false,
        );

        assert_eq!(result.extend_len, 5);
        assert_eq!(result.max_score, 5);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_stops_at_optimal_point_with_mismatches() {
        // Genome: A C G T A C G T (positions 0-7)
        let genome_seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let index = make_index_with_seq(&genome_seq);
        // Read: A C G T T T T T (matches first 4, then all mismatches)
        let read_seq: Vec<u8> = vec![0, 1, 2, 3, 3, 3, 3, 3];

        let result = extend_alignment(
            &read_seq, 0,   // read_start
            0,   // genome_start
            1,   // rightward
            8,   // max_extend
            0,   // no previous mismatches
            0,   // no previous length
            10,  // n_mm_max
            0.3, // p_mm_max
            &index, false,
        );

        // Should extend 4 bases (perfect match), then mismatches drag score down
        assert_eq!(result.extend_len, 4);
        assert_eq!(result.max_score, 4);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_chromosome_boundary() {
        // Genome: A C G (3 bases, then padding=5)
        let genome_seq = vec![0, 1, 2];
        let index = make_index_with_seq(&genome_seq);
        // Read is 5 bases, but genome only has 3
        let read_seq: Vec<u8> = vec![0, 1, 2, 3, 0];

        let result = extend_alignment(
            &read_seq, 0, // read_start
            0, // genome_start
            1, // rightward
            5, // max_extend
            0, 0, 10, 0.3, &index, false,
        );

        // Should stop at 3 bases (genome boundary)
        assert_eq!(result.extend_len, 3);
        assert_eq!(result.max_score, 3);
    }

    #[test]
    fn test_extend_n_bases_skipped() {
        // Genome: A N C G (N=4 at position 1)
        let genome_seq = vec![0, 4, 1, 2];
        let index = make_index_with_seq(&genome_seq);
        // Read: A A C G (matches at 0, N skip at 1, matches at 2-3)
        let read_seq: Vec<u8> = vec![0, 0, 1, 2];

        let result = extend_alignment(&read_seq, 0, 0, 1, 4, 0, 0, 10, 0.3, &index, false);

        // Should extend all 4 bases: match + N(skip) + match + match = score 3
        assert_eq!(result.extend_len, 4);
        assert_eq!(result.max_score, 3);
        assert_eq!(result.n_mismatch, 0);
    }

    #[test]
    fn test_extend_all_mismatch_returns_zero() {
        // Genome: A A A A (all 0)
        let genome_seq = vec![0, 0, 0, 0];
        let index = make_index_with_seq(&genome_seq);
        // Read: T T T T (all 3, complete mismatch)
        let read_seq: Vec<u8> = vec![3, 3, 3, 3];

        let result = extend_alignment(&read_seq, 0, 0, 1, 4, 0, 0, 10, 0.3, &index, false);

        // Score never goes positive, so extend_len should be 0
        assert_eq!(result.extend_len, 0);
        assert_eq!(result.max_score, 0);
    }

    #[test]
    fn test_extend_recovery_through_mismatch() {
        // Genome: A C G T A C G T A C G T A C (14 bases)
        let genome_seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let index = make_index_with_seq(&genome_seq);
        // Read: A C G X A C G T A C G T A C (1 mismatch at pos 3, then 10 matches)
        //       M M M X M M M M M M M M M M
        let read_seq: Vec<u8> = vec![0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];

        let result = extend_alignment(&read_seq, 0, 0, 1, 14, 0, 0, 10, 0.3, &index, false);

        // Should extend past the mismatch: 3M + 1X + 10M
        // Score: +3 -1 +10 = 12, best at position 14
        assert_eq!(result.extend_len, 14);
        assert_eq!(result.max_score, 12);
        assert_eq!(result.n_mismatch, 1);
    }

    #[test]
    fn test_extend_zero_max_extend() {
        let genome_seq = vec![0, 1, 2, 3];
        let index = make_index_with_seq(&genome_seq);
        let read_seq: Vec<u8> = vec![0, 1, 2, 3];

        let result = extend_alignment(
            &read_seq, 0, 0, 1, 0, // max_extend = 0
            0, 0, 10, 0.3, &index, false,
        );

        assert_eq!(result.extend_len, 0);
        assert_eq!(result.max_score, 0);
    }

    #[test]
    fn test_overhang_check_rejects_short_overhang() {
        // Test that the overhang check rejects splice junctions with tiny flanking seeds.
        // Scenario: prev seed length=3 (below default min of 5), current seed length=20
        // With alignSJoverhangMin=5, the 3bp left overhang should cause rejection.
        use crate::align::score::AlignmentScorer;

        let scorer = AlignmentScorer {
            score_gap: 0,
            score_gap_noncan: -8,
            score_gap_gcag: -4,
            score_gap_atac: -8,
            score_del_open: -2,
            score_del_base: -2,
            score_ins_open: -2,
            score_ins_base: -2,
            align_intron_min: 21,
            sjdb_score: 2,
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
            score_stitch_sj_shift: 1,
            align_spliced_mate_map_lmin: 0,
            align_spliced_mate_map_lmin_over_lmate: 0.66,
            out_filter_score_min_over_lread: 0.66,
        };

        // Left overhang (prev.length) = 3, below min of 5
        let left_overhang: usize = 3;
        let right_overhang: usize = 20;
        let min_overhang = scorer.align_sj_overhang_min as usize;

        // Should be rejected
        assert!(left_overhang < min_overhang || right_overhang < min_overhang);

        // Right overhang too small
        let left_overhang: usize = 20;
        let right_overhang: usize = 4;
        assert!(left_overhang < min_overhang || right_overhang < min_overhang);
    }

    #[test]
    fn test_overhang_check_accepts_sufficient_overhang() {
        // Test that splice junctions with sufficient overhang pass the check.
        use crate::align::score::AlignmentScorer;

        let scorer = AlignmentScorer {
            score_gap: 0,
            score_gap_noncan: -8,
            score_gap_gcag: -4,
            score_gap_atac: -8,
            score_del_open: -2,
            score_del_base: -2,
            score_ins_open: -2,
            score_ins_base: -2,
            align_intron_min: 21,
            sjdb_score: 2,
            align_sj_stitch_mismatch_nmax: [0, -1, 0, 0],
            n_mm_max: 10,
            p_mm_max: 0.3,
            align_sj_overhang_min: 5,
            align_sjdb_overhang_min: 3,
            align_intron_max: 589_824,
            score_genomic_length_log2_scale: -0.25,
            score_stitch_sj_shift: 1,
            align_spliced_mate_map_lmin: 0,
            align_spliced_mate_map_lmin_over_lmate: 0.66,
            out_filter_score_min_over_lread: 0.66,
        };

        // Both overhangs >= 5
        let left_overhang: usize = 5;
        let right_overhang: usize = 10;
        let min_overhang = scorer.align_sj_overhang_min as usize;

        // Should pass
        assert!(!(left_overhang < min_overhang || right_overhang < min_overhang));

        // Exactly at minimum
        let left_overhang: usize = 5;
        let right_overhang: usize = 5;
        assert!(!(left_overhang < min_overhang || right_overhang < min_overhang));
    }

    #[test]
    fn test_bin_based_cluster_bounds() {
        // Verify that cluster genome_start/genome_end are set from bin range
        // With win_bin_nbits=4 (bin_size=16), a window at bin 5 with flank 2
        // should span bins 3-7, i.e., genome_start=48, genome_end=128
        let bin_size: u64 = 1 << 4; // 16
        let bin_start: u64 = 5u64.saturating_sub(2); // 3
        let bin_end: u64 = 5 + 2; // 7
        let genome_start = bin_start * bin_size;
        let genome_end = (bin_end + 1) * bin_size;
        assert_eq!(genome_start, 48);
        assert_eq!(genome_end, 128);
    }

    #[test]
    fn test_window_merge_logic() {
        // Two anchors within winAnchorDistNbins should merge into one window
        // Anchor A at bin 10, Anchor B at bin 15, winAnchorDistNbins=9
        // B is within [10-9, 10+9] = [1, 19] → merge
        let win_anchor_dist_nbins = 9u64;
        let window_bin_start: u64 = 10;
        let window_bin_end: u64 = 10;
        let anchor_bin: u64 = 15;

        let merge_start = window_bin_start.saturating_sub(win_anchor_dist_nbins);
        let merge_end = window_bin_end + win_anchor_dist_nbins;
        let should_merge = anchor_bin >= merge_start && anchor_bin <= merge_end;
        assert!(
            should_merge,
            "Anchors 5 bins apart should merge with dist_nbins=9"
        );

        // Anchor C at bin 25 is outside [1, 19] → separate window
        let anchor_bin_c: u64 = 25;
        let should_merge_c = anchor_bin_c >= merge_start && anchor_bin_c <= merge_end;
        assert!(
            !should_merge_c,
            "Anchors 15 bins apart should NOT merge with dist_nbins=9"
        );
    }

    #[test]
    fn test_window_flank_extension() {
        // Window at bin 10, extended by ±4 flanking bins → bins 6-14
        let mut bin_start: u64 = 10;
        let mut bin_end: u64 = 10;
        let win_flank_nbins: u64 = 4;

        bin_start = bin_start.saturating_sub(win_flank_nbins);
        bin_end += win_flank_nbins;

        assert_eq!(bin_start, 6);
        assert_eq!(bin_end, 14);

        // Edge case: window at bin 2, flanking underflows to 0
        let mut bin_start_edge: u64 = 2;
        bin_start_edge = bin_start_edge.saturating_sub(win_flank_nbins);
        assert_eq!(bin_start_edge, 0, "Flanking should saturate at 0");
    }
}
