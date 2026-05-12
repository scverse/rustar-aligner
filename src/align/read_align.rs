/// Read alignment driver function
use crate::align::score::{AlignmentScorer, SpliceMotif};
use crate::align::seed::Seed;
use crate::align::stitch::{
    PE_SPACER_BASE, cluster_seeds, finalize_transcript, split_combined_wt, stitch_seeds_core,
    stitch_seeds_with_jdb_debug,
};
use crate::align::transcript::{Exon, Transcript};
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::{IntronMotifFilter, IntronStrandFilter, Parameters};
use crate::stats::UnmappedReason;
use rand::{SeedableRng, rngs::StdRng, seq::SliceRandom};
use std::hash::{DefaultHasher, Hash, Hasher};

/// Derive a deterministic per-read RNG seed from `run_rng_seed` + the read name.
///
/// STAR seeds `std::mt19937` once per chunk/thread (`runRNGseed*(iChunk+1)`),
/// then advances the state sequentially per read. rustar-aligner parallelises per-read
/// via rayon, so we instead fold the read name into the seed — this keeps tie
/// breaks reproducible regardless of thread count while still honoring the
/// user's `--runRNGseed` value.
pub(crate) fn per_read_seed(run_rng_seed: u64, read_name: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    read_name.hash(&mut hasher);
    run_rng_seed.wrapping_mul(hasher.finish().wrapping_add(1))
}

/// Shuffle the prefix of `items` whose `score_fn` equals the first element's score.
///
/// Mirrors STAR's `ReadAlign_multMapSelect` / `funPrimaryAlignMark`: best-scoring
/// alignments are randomized so primary selection (index 0) is not biased by
/// upstream sort order. Non-tied elements are left alone.
fn shuffle_tied_prefix<T>(items: &mut [T], score_fn: impl Fn(&T) -> i32, seed: u64) {
    let Some(first) = items.first() else {
        return;
    };
    let best = score_fn(first);
    let tied = items.iter().take_while(|t| score_fn(t) == best).count();
    if tied < 2 {
        return;
    }
    items[..tied].shuffle(&mut StdRng::seed_from_u64(seed));
}

/// Result of aligning a single read: (transcripts, chimeric_alignments, n_for_mapq, unmapped_reason)
pub type AlignReadResult = (
    Vec<Transcript>,
    Vec<crate::chimeric::ChimericAlignment>,
    usize,
    Option<UnmappedReason>,
);

/// Paired-end alignment result
#[derive(Debug, Clone)]
pub struct PairedAlignment {
    /// Transcript for mate1
    pub mate1_transcript: Transcript,
    /// Transcript for mate2
    pub mate2_transcript: Transcript,
    /// Read positions for mate1 in transcript (start, end)
    pub mate1_region: (usize, usize),
    /// Read positions for mate2 in transcript (start, end)
    pub mate2_region: (usize, usize),
    /// Whether this is a proper pair (same chr, concordant orientation, distance)
    pub is_proper_pair: bool,
    /// Signed insert size (TLEN) - genomic distance between mate starts
    pub insert_size: i32,
    /// Combined pair score: sum of per-mate finalized scores (each includes genomic length penalty).
    /// Used for multi-mapper score-range ranking and mappedFilter quality check.
    pub combined_wt_score: i32,
    /// Combined coverage: sum of exon read spans from both mates.
    /// Mirrors STAR's nMatch check in mappedFilter.
    pub combined_n_match: u32,
}

impl PairedAlignment {
    /// Build a STAR-style combined two-mate `Transcript` for transcriptome
    /// projection.
    ///
    /// Matches STAR's single-`Transcript`-per-pair model: mate1 exons with
    /// `i_frag = 0`, then mate2 exons rewritten to `i_frag = 1`. Only
    /// meaningful for pairs on the same chromosome and strand — both are
    /// invariants of a `PairedAlignment` (checked in `try_pair_transcripts`).
    ///
    /// The returned transcript's `cigar` is empty: transcriptome BAM
    /// emission generates per-mate CIGARs from the split exon list rather
    /// than consuming a combined one.
    pub fn combined_transcript_for_projection(&self) -> Transcript {
        let m1 = &self.mate1_transcript;
        let m2 = &self.mate2_transcript;

        let mut exons: Vec<Exon> = Vec::with_capacity(m1.exons.len() + m2.exons.len());
        for e in &m1.exons {
            let mut ee = e.clone();
            ee.i_frag = 0;
            exons.push(ee);
        }
        for e in &m2.exons {
            let mut ee = e.clone();
            ee.i_frag = 1;
            exons.push(ee);
        }

        Transcript {
            chr_idx: m1.chr_idx,
            genome_start: m1.genome_start.min(m2.genome_start),
            genome_end: m1.genome_end.max(m2.genome_end),
            is_reverse: m1.is_reverse,
            exons,
            cigar: Vec::new(),
            score: m1.score + m2.score,
            n_mismatch: m1.n_mismatch + m2.n_mismatch,
            n_gap: m1.n_gap + m2.n_gap,
            n_junction: m1.n_junction + m2.n_junction,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            read_seq: Vec::new(),
        }
    }
}

/// Result of paired-end alignment, covering all mapping outcomes.
#[derive(Debug, Clone)]
pub enum PairedAlignmentResult {
    /// Both mates mapped and paired successfully
    BothMapped(Box<PairedAlignment>),
    /// Only one mate mapped; rescue failed or was not attempted for the other
    HalfMapped {
        /// Transcript of the mapped mate
        mapped_transcript: Transcript,
        /// true = mate1 is the mapped mate, false = mate2 is mapped
        mate1_is_mapped: bool,
    },
}

/// Align a read to the genome.
///
/// # Algorithm
/// 1. Find seeds (exact matches) using MMP search
/// 2. Cluster seeds by genomic proximity
/// 3. Stitch seeds within each cluster using DP
/// 4. Filter transcripts by quality thresholds
/// 5. Sort by score and limit to top N
/// 6. Detect chimeric alignments if enabled
///
/// # Arguments
/// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
/// * `read_name` - Read name (needed for chimeric output)
/// * `index` - Genome index
/// * `params` - User parameters
///
/// # Returns
/// Tuple of (transcripts, chimeric alignments, n_for_mapq, unmapped_reason):
/// - transcripts: sorted by score (best first)
/// - chimeric alignments: sorted by score (best first)
/// - n_for_mapq: effective alignment count for MAPQ calculation (max of transcript count
///   and valid cluster count, to avoid undercounting from coordinate dedup on tandem repeats)
/// - unmapped_reason: `Some(reason)` if no alignments produced, `None` if mapped
pub fn align_read(
    read_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<AlignReadResult, Error> {
    let debug_read = !params.read_name_filter.is_empty() && read_name == params.read_name_filter;

    // Step 1: Find seeds (seedMapMin from params)
    let min_seed_length = params.seed_map_min;
    let seeds = Seed::find_seeds(
        read_seq,
        index,
        min_seed_length,
        params,
        if debug_read { read_name } else { "" },
    )?;

    if debug_read {
        let total_positions: usize = seeds.iter().map(|s| s.sa_end - s.sa_start).sum();
        eprintln!(
            "[DEBUG {}] Seeds: {} seeds, {} total SA positions, read_len={}",
            read_name,
            seeds.len(),
            total_positions,
            read_seq.len()
        );
        let n_lr = seeds.iter().filter(|s| !s.search_rc).count();
        let n_rl = seeds.iter().filter(|s| s.search_rc).count();
        eprintln!(
            "  {} seeds total: {} L→R (sparse), {} R→L (sparse)",
            seeds.len(),
            n_lr,
            n_rl
        );
        for (i, s) in seeds.iter().enumerate() {
            let n_loci = s.sa_end - s.sa_start;
            eprintln!(
                "  seed[{}]: read_pos={}, len={}, n_loci={}, search_rc={}, sa=[{},{})",
                i, s.read_pos, s.length, n_loci, s.search_rc, s.sa_start, s.sa_end
            );
        }
    }

    if seeds.is_empty() {
        if debug_read {
            eprintln!("[DEBUG {}] No seeds found — unmapped", read_name);
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Step 2: Cluster seeds (STAR's bin-based windowing)
    // seed_per_window_nmax capacity eviction is handled inside cluster_seeds()
    let clusters = cluster_seeds(&seeds, index, params, read_seq.len(), debug_read);

    if debug_read {
        eprintln!(
            "[DEBUG {}] Clusters: {} clusters",
            read_name,
            clusters.len()
        );
        for (i, cluster) in clusters.iter().enumerate() {
            let chr_name = if cluster.chr_idx < index.genome.chr_name.len() {
                &index.genome.chr_name[cluster.chr_idx]
            } else {
                "unknown"
            };
            eprintln!(
                "  cluster[{}]: chr={}, is_reverse={}, seeds={}, anchor_bin={}",
                i,
                chr_name,
                cluster.is_reverse,
                cluster.alignments.len(),
                cluster.anchor_bin,
            );
            for (j, wa) in cluster.alignments.iter().enumerate() {
                let chr_pos = wa.genome_pos.saturating_sub(
                    if cluster.chr_idx < index.genome.chr_start.len() {
                        index.genome.chr_start[cluster.chr_idx]
                    } else {
                        0
                    },
                ) + 1; // 1-based
                eprintln!(
                    "    wa[{}]: read_pos={}, len={}, genome_pos={} ({}:{}), n_rep={}, is_anchor={}",
                    j,
                    wa.read_pos,
                    wa.length,
                    wa.genome_pos,
                    chr_name,
                    chr_pos,
                    wa.n_rep,
                    wa.is_anchor
                );
                if j >= 5 {
                    eprintln!(
                        "    ... ({} more WA entries)",
                        cluster.alignments.len() - j - 1
                    );
                    break;
                }
            }
        }
    }

    if clusters.is_empty() {
        if debug_read {
            eprintln!("[DEBUG {}] No clusters — unmapped", read_name);
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Cap total clusters (alignWindowsPerReadNmax)
    let mut clusters = clusters;
    clusters.truncate(params.align_windows_per_read_nmax);

    // NOTE: STAR's winReadCoverageRelativeMin filter is long-reads-only
    // (#ifdef COMPILE_FOR_LONG_READS in stitchPieces.cpp). Standard STAR
    // does NOT filter clusters by seed coverage. Removed to match STAR.

    // Step 2b: Detect chimeric alignments from multi-cluster seeds (Tier 2)
    let mut chimeric_alignments = Vec::new();
    if params.chim_segment_min > 0 && clusters.len() > 1 {
        use crate::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        chimeric_alignments
            .extend(detector.detect_from_multi_clusters(&clusters, read_seq, read_name, index)?);
    }

    // Step 3: Stitch seeds within each cluster
    let scorer = AlignmentScorer::from_params(params);
    let mut transcripts = Vec::new();
    // Collect all raw (pre-dedup) transcripts for chimericDetectionOld (Tier 1).
    let mut all_raw_transcripts: Vec<crate::align::transcript::Transcript> = Vec::new();

    // Use junction DB for annotation-aware scoring if available
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };

    for (ci, cluster) in clusters.iter().enumerate() {
        let debug_name = if debug_read { read_name } else { "" };
        let cluster_transcripts = stitch_seeds_with_jdb_debug(
            cluster,
            read_seq,
            index,
            &scorer,
            junction_db,
            params.align_transcripts_per_window_nmax,
            debug_name,
        )?;
        if debug_read {
            eprintln!(
                "[DEBUG {}] Cluster[{}]: {} transcripts from DP",
                read_name,
                ci,
                cluster_transcripts.len()
            );
            for (ti, t) in cluster_transcripts.iter().enumerate().take(5) {
                let chr_name = if t.chr_idx < index.genome.chr_name.len() {
                    &index.genome.chr_name[t.chr_idx]
                } else {
                    "unknown"
                };
                let cigar_str: String = t.cigar.iter().map(|op| format!("{}", op)).collect();
                eprintln!(
                    "  transcript[{}]: chr={}:{}-{} ({}) score={} mm={} junctions={} cigar={}",
                    ti,
                    chr_name,
                    t.genome_start,
                    t.genome_end,
                    if t.is_reverse { "-" } else { "+" },
                    t.score,
                    t.n_mismatch,
                    t.n_junction,
                    cigar_str
                );
            }
        }
        if params.chim_segment_min > 0 {
            all_raw_transcripts.extend(cluster_transcripts.iter().cloned());
        }
        transcripts.extend(cluster_transcripts);
    }

    // Step 4a: Deduplicate and score-range filter — BEFORE quality filters.
    // STAR order: multMapSelect (score-range) → mappedFilter (quality gates).
    // Doing quality filters first is wrong: it can remove the high-scoring primary,
    // leaving a lower-scoring secondary that then passes as the "best" alignment.

    // Deduplicate transcripts with identical genomic coordinates AND CIGAR.
    transcripts.sort_by(|a, b| {
        (
            a.chr_idx,
            a.genome_start,
            a.genome_end,
            a.is_reverse,
            &a.cigar,
        )
            .cmp(&(
                b.chr_idx,
                b.genome_start,
                b.genome_end,
                b.is_reverse,
                &b.cigar,
            ))
            .then_with(|| b.score.cmp(&a.score))
    });
    transcripts.dedup_by(|a, b| {
        a.chr_idx == b.chr_idx
            && a.genome_start == b.genome_start
            && a.genome_end == b.genome_end
            && a.is_reverse == b.is_reverse
            && a.cigar == b.cigar
    });

    // Sort by score descending (deterministic tie-breaking).
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.n_junction.cmp(&b.n_junction))
            .then_with(|| a.chr_idx.cmp(&b.chr_idx))
            .then_with(|| a.genome_start.cmp(&b.genome_start))
            .then_with(|| a.is_reverse.cmp(&b.is_reverse))
    });

    // Randomize primary among best-scoring ties (ReadAlign_multMapSelect.cpp:71-79).
    shuffle_tied_prefix(
        &mut transcripts,
        |t| t.score,
        per_read_seed(params.run_rng_seed, read_name),
    );

    // Score-range filter: keep only alignments within outFilterMultimapScoreRange of the best.
    // (STAR's multMapSelect step — must run before quality filters.)
    if !transcripts.is_empty() {
        let max_score = transcripts[0].score;
        let score_threshold = max_score - params.out_filter_multimap_score_range;
        transcripts.retain(|t| t.score >= score_threshold);
    }

    // Multimap count check: too many loci → unmapped.
    if transcripts.len() > params.out_filter_multimap_nmax as usize {
        let n_loci = transcripts.len();
        transcripts.clear();
        return Ok((
            transcripts,
            chimeric_alignments,
            n_loci,
            Some(UnmappedReason::TooManyLoci),
        ));
    }

    // Step 4: Quality filters (STAR's mappedFilter — runs after score-range selection).
    // STAR uses (Lread-1) for relative thresholds and casts to integer
    // (ReadAlign_mappedFilter.cpp lines 8-9)
    let read_length = read_seq.len() as f64;
    let lread_m1 = (read_seq.len() as f64) - 1.0;

    // Log filtering statistics
    let pre_filter_count = transcripts.len();
    let mut filter_reasons = std::collections::HashMap::new();

    transcripts.retain(|t| {
        // Absolute score threshold
        if t.score < params.out_filter_score_min {
            *filter_reasons.entry("score_min").or_insert(0) += 1;
            return false;
        }

        // Relative score threshold: STAR casts to intScore (i32)
        if t.score < (params.out_filter_score_min_over_lread * lread_m1) as i32 {
            *filter_reasons.entry("score_min_relative").or_insert(0) += 1;
            return false;
        }

        // Absolute mismatch count
        if t.n_mismatch > params.out_filter_mismatch_nmax {
            *filter_reasons.entry("mismatch_max").or_insert(0) += 1;
            log::debug!(
                "Filtered {}: {} mismatches > {} max (read_len={}, score={})",
                read_name,
                t.n_mismatch,
                params.out_filter_mismatch_nmax,
                read_length,
                t.score
            );
            return false;
        }

        // Relative mismatch count (mismatches / read_length)
        let mismatch_rate = t.n_mismatch as f64 / read_length;
        if mismatch_rate > params.out_filter_mismatch_nover_lmax {
            *filter_reasons.entry("mismatch_rate").or_insert(0) += 1;
            log::debug!(
                "Filtered {}: {:.1}% mismatch rate > {:.1}% max ({}/{} bases, score={})",
                read_name,
                mismatch_rate * 100.0,
                params.out_filter_mismatch_nover_lmax * 100.0,
                t.n_mismatch,
                read_length,
                t.score
            );
            return false;
        }

        // Absolute matched bases
        let n_matched = t.n_matched();
        if n_matched < params.out_filter_match_nmin {
            *filter_reasons.entry("match_min").or_insert(0) += 1;
            return false;
        }

        // Relative matched bases: STAR casts to uint (u32)
        if n_matched < (params.out_filter_match_nmin_over_lread * lread_m1) as u32 {
            *filter_reasons.entry("match_min_relative").or_insert(0) += 1;
            return false;
        }

        // Junction motif filtering
        match params.out_filter_intron_motifs {
            IntronMotifFilter::None => {
                // Accept all motifs
            }
            IntronMotifFilter::RemoveNoncanonical => {
                // Reject if any junction is non-canonical
                if t.junction_motifs.contains(&SpliceMotif::NonCanonical) {
                    *filter_reasons.entry("noncanonical_junction").or_insert(0) += 1;
                    return false;
                }
            }
            IntronMotifFilter::RemoveNoncanonicalUnannotated => {
                // Only reject if a non-canonical junction is NOT annotated in GTF
                if t.junction_motifs
                    .iter()
                    .zip(t.junction_annotated.iter())
                    .any(|(m, annotated)| *m == SpliceMotif::NonCanonical && !annotated)
                {
                    *filter_reasons
                        .entry("noncanonical_unannotated_junction")
                        .or_insert(0) += 1;
                    return false;
                }
            }
        }

        // Intron strand consistency filtering (outFilterIntronStrands)
        // STAR's RemoveInconsistentStrands removes transcripts that have junctions
        // with mixed intron strand (some imply + strand, some imply - strand).
        // This handles chimeric/impossible transcripts spanning both strands.
        // Note: a reverse-strand read CAN have + strand motifs (antisense reads
        // from + strand genes in unstranded RNA-seq) — this is valid and STAR
        // keeps such reads. Only mixed-strand within one transcript is filtered.
        if params.out_filter_intron_strands == IntronStrandFilter::RemoveInconsistentStrands {
            let mut has_plus = false;
            let mut has_minus = false;
            for motif in &t.junction_motifs {
                match motif.implied_strand() {
                    Some('+') => has_plus = true,
                    Some('-') => has_minus = true,
                    None => {}
                    _ => {}
                }
            }
            if has_plus && has_minus {
                *filter_reasons.entry("inconsistent_strand").or_insert(0) += 1;
                return false;
            }
        }

        true
    });

    // Log filtering summary if anything was filtered
    if pre_filter_count > transcripts.len() {
        let filtered = pre_filter_count - transcripts.len();
        log::debug!(
            "Read {}: Filtered {}/{} transcripts: {:?}",
            read_name,
            filtered,
            pre_filter_count,
            filter_reasons
        );
    }

    if debug_read {
        eprintln!(
            "[DEBUG {}] After quality filters: {}/{} transcripts remain (reasons: {:?})",
            read_name,
            transcripts.len(),
            pre_filter_count,
            filter_reasons
        );
    }

    // Step 3b: chimericDetectionOld (Tier 1) — STAR-faithful post-stitching transcript-pair search.
    // Uses the best post-dedup transcript as the primary segment and searches all raw transcripts
    // (pre-dedup, from all clusters) for the best complementary partner.
    if params.chim_segment_min > 0
        && !all_raw_transcripts.is_empty()
        && let Some(tr_best) = transcripts.first()
    {
        use crate::chimeric::detect_chimeric_old;
        let chims = detect_chimeric_old(
            &all_raw_transcripts,
            tr_best,
            read_seq,
            read_name,
            params,
            index,
        )?;
        chimeric_alignments.extend(chims);
    }

    // Step 3c: Soft-clip re-mapping (Phase 12.2) — try to align soft-clipped bases when
    // detect_chimeric_old found no chimeric partner in the existing transcript pool.
    if params.chim_segment_min > 0
        && chimeric_alignments.is_empty()
        && let Some(tr_best) = transcripts.first()
    {
        use crate::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        if let Some(chim) = detector.detect_from_soft_clips(tr_best, read_seq, read_name, index)? {
            chimeric_alignments.push(chim);
        }
    }

    // Step 3d: Tier 3 — re-seed outer uncovered read regions of each chimeric pair (Phase 17.10).
    // Extends 2-segment chimeras toward multi-junction fusions by seeding the read bases
    // that lie outside both chimeric segments.
    if params.chim_segment_min > 0 && !chimeric_alignments.is_empty() {
        use crate::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        let mut tier3 = Vec::new();
        for chim in &chimeric_alignments {
            let extras =
                detector.detect_from_chimeric_residuals(chim, read_seq, read_name, index)?;
            tier3.extend(extras);
        }
        chimeric_alignments.extend(tier3);
    }

    // Note: STAR sometimes finds 2 equivalent indel placements in homopolymer runs
    // via its recursive stitcher's seed exploration (NH=2 instead of NH=1 for ~5 reads).
    // Generating equivalents post-hoc causes more harm than good (41 false NH=2 vs 5 fixed).
    // The root cause is jR scanning placing insertions at different positions — fixing that
    // would be a better approach than post-hoc enumeration.

    // Step 6: Filter chimeric alignments
    if params.chim_segment_min > 0 {
        chimeric_alignments.retain(|chim| {
            chim.meets_min_segment_length(params.chim_segment_min)
                && chim.meets_min_score(params.chim_score_min)
        });
    }

    // n_for_mapq = transcripts.len() after dedup and filtering.
    // Multi-transcript DP (Phase 16.10) produces multiple transcripts per window
    // for tandem repeats (e.g. rDNA), yielding correct NH → correct MAPQ.
    let n_for_mapq = transcripts.len();

    if debug_read {
        eprintln!(
            "[DEBUG {}] Final: {} transcripts, n_for_mapq={}",
            read_name,
            transcripts.len(),
            n_for_mapq
        );
        for (i, t) in transcripts.iter().enumerate() {
            let chr_name = if t.chr_idx < index.genome.chr_name.len() {
                &index.genome.chr_name[t.chr_idx]
            } else {
                "unknown"
            };
            let cigar_str: String = t.cigar.iter().map(|op| format!("{}", op)).collect();
            eprintln!(
                "  FINAL[{}]: chr={}:{}-{} ({}) score={} mm={} junctions={} cigar={}",
                i,
                chr_name,
                t.genome_start,
                t.genome_end,
                if t.is_reverse { "-" } else { "+" },
                t.score,
                t.n_mismatch,
                t.n_junction,
                cigar_str
            );
        }
    }

    let unmapped_reason = if transcripts.is_empty() {
        // Transcripts were generated by DP but all filtered out
        Some(UnmappedReason::TooShort)
    } else {
        None
    };

    Ok((
        transcripts,
        chimeric_alignments,
        n_for_mapq,
        unmapped_reason,
    ))
}

type PairedAlignResult = (
    Vec<PairedAlignmentResult>,
    Vec<crate::chimeric::ChimericAlignment>,
    usize,
    Option<UnmappedReason>,
);

/// Align paired-end reads using STAR's combined-read approach.
///
/// # Algorithm
/// 1. Build combined read: [mate1_seq | PE_SPACER_BASE | RC(mate2_seq)]
/// 2. Seed each fragment (mate1_seq and RC(mate2_seq)) independently with per-fragment Nstart
/// 3. Tag seeds with mate_id (0=mate1, 1=mate2); adjust read_pos to combined-read coords
/// 4. Cluster and stitch combined seeds → WorkingTranscripts spanning both mates
/// 5. Split each WT by mate_id → finalize each half → pair
/// 6. Decision tree: dedup → score-range → TooManyLoci → quality filter
/// 7. Half-mapped fallback from single-mate WTs
///
/// # Arguments
/// * `mate1_seq` - First mate sequence (encoded)
/// * `mate2_seq` - Second mate sequence (encoded)
/// * `index` - Genome index
/// * `params` - Parameters (includes alignMatesGapMax)
///
/// # Returns
/// Tuple of (paired alignment results, n_for_mapq, unmapped_reason)
pub fn align_paired_read(
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<PairedAlignResult, Error> {
    let len1 = mate1_seq.len();
    let len2 = mate2_seq.len();

    let debug_pe = !params.read_name_filter.is_empty() && read_name == params.read_name_filter;
    let scorer = AlignmentScorer::from_params(params);
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };
    let debug_name: &str = if debug_pe { read_name } else { "" };

    // Build combined read: [mate1_seq | PE_SPACER_BASE | RC(mate2_seq)]
    // STAR ReadAlign_oneRead.cpp: Read1[0][readLength[0]] = MARK_FRAG_SPACER_BASE
    let rc_mate2: Vec<u8> = mate2_seq
        .iter()
        .rev()
        .map(|&b| if b < 4 { 3 - b } else { b })
        .collect();
    let mut combined_read = Vec::with_capacity(len1 + 1 + len2);
    combined_read.extend_from_slice(mate1_seq);
    combined_read.push(PE_SPACER_BASE);
    combined_read.extend_from_slice(&rc_mate2);
    let combined_len = combined_read.len();

    // STAR-faithful per-fragment seeding: seed each mate fragment separately using
    // the fragment length for Nstart/Lstart, then merge into combined_read coords.
    // STAR uses qualitySplit() starting positions based on fragment length (e.g. 150bp
    // → Nstart=4, Lstart=37, starts={0,37,74,111}), NOT the combined length (301bp
    // → Nstart=7, starts={0,43,...,129,...}). Using combined length creates a spurious
    // start at position 129 (between mates) that can produce anchors widening windows
    // beyond STAR's range, causing window overflow and eviction of valid 7M exon seeds.
    let mut combined_seeds = Seed::find_seeds(
        &combined_read[..len1],
        index,
        params.seed_map_min,
        params,
        debug_name,
    )?;
    let mut m2_seeds = Seed::find_seeds(
        &combined_read[len1 + 1..],
        index,
        params.seed_map_min,
        params,
        if debug_pe { debug_name } else { "" },
    )?;
    for s in &mut m2_seeds {
        s.read_pos += len1 + 1;
    }
    combined_seeds.extend(m2_seeds);
    // mate_id: positions 0..len1 → mate1(0); positions len1+1.. → RC(mate2)(1).
    for s in &mut combined_seeds {
        s.mate_id = if s.read_pos < len1 { 0 } else { 1 };
    }

    // Cluster combined seeds using the combined read length
    let clusters = cluster_seeds(&combined_seeds, index, params, combined_len, debug_pe);

    // PE chimeric pre-pass: intra-mate multi-cluster detection (Tier 2).
    // Split clusters by mate_id and run per-mate chimeric detection, mirroring SE behavior.
    let mut pe_chimeric: Vec<crate::chimeric::ChimericAlignment> = Vec::new();
    if params.chim_segment_min > 0 && clusters.len() >= 2 {
        use crate::chimeric::ChimericDetector;

        let mate1_clusters: Vec<_> = clusters
            .iter()
            .filter(|c| c.alignments.iter().all(|wa| wa.mate_id == 0))
            .cloned()
            .collect();

        // Mate2 clusters: adjust read_pos to be relative to mate2_seq (subtract len1+1)
        let mate2_clusters: Vec<_> = clusters
            .iter()
            .filter(|c| c.alignments.iter().all(|wa| wa.mate_id == 1))
            .map(|c| {
                let mut c2 = c.clone();
                for wa in &mut c2.alignments {
                    wa.read_pos -= len1 + 1;
                }
                c2
            })
            .collect();

        let detector = ChimericDetector::new(params);
        if mate1_clusters.len() >= 2 {
            pe_chimeric.extend(detector.detect_from_multi_clusters(
                &mate1_clusters,
                mate1_seq,
                read_name,
                index,
            )?);
        }
        if mate2_clusters.len() >= 2 {
            pe_chimeric.extend(detector.detect_from_multi_clusters(
                &mate2_clusters,
                mate2_seq,
                read_name,
                index,
            )?);
        }
        pe_chimeric.retain(|c| {
            c.meets_min_segment_length(params.chim_segment_min)
                && c.meets_min_score(params.chim_score_min)
        });
    }

    // Combined score threshold: use len1+len2 as denominator
    let combined_score_threshold =
        (params.out_filter_score_min_over_lread * (len1 + len2) as f64) as i32;

    let mut joint_pairs: Vec<PairedAlignment> = Vec::new();
    let mut single_mate1_transcripts: Vec<Transcript> = Vec::new();
    let mut single_mate2_transcripts: Vec<Transcript> = Vec::new();
    // All finalized mate transcripts (from both joint pairs and single-mate WTs) used
    // as the search pool for chimericDetectionOld (Tier 1) on each mate independently.
    let mut all_m1_transcripts: Vec<Transcript> = Vec::new();
    let mut all_m2_transcripts: Vec<Transcript> = Vec::new();

    // Stitch combined clusters, split WTs by mate_id, finalize each half
    for cluster in clusters.iter().take(params.align_windows_per_read_nmax) {
        let (wts, stitch_cluster, stitch_is_reverse, stitch_read) = stitch_seeds_core(
            cluster,
            &combined_read,
            index,
            &scorer,
            junction_db,
            params.align_transcripts_per_window_nmax,
            params.align_mates_gap_max.into(),
            debug_name,
        )?;

        for wt in &wts {
            let split_result =
                split_combined_wt(wt, len1, len2, stitch_is_reverse, scorer.align_intron_min);
            match split_result {
                Some((m1_wt, m2_wt)) => {
                    let (m1_read_slice, m1_orig_rev, m2_read_slice, m2_orig_rev) =
                        if stitch_is_reverse {
                            // stitch_read = [mate2(0..len2) | SPACER | RC(mate1)(len2+1..)]
                            (
                                &stitch_read[len2 + 1..], // RC(mate1_seq)
                                true,                     // mate1 5' at right in RC
                                &stitch_read[..len2],     // mate2_seq
                                false,                    // mate2 5' at left
                            )
                        } else {
                            // stitch_read = [mate1(0..len1) | SPACER | RC(mate2)(len1+1..)]
                            (
                                &stitch_read[..len1],     // mate1_seq
                                false,                    // mate1 5' at left
                                &stitch_read[len1 + 1..], // RC(mate2_seq)
                                true,                     // mate2 5' at right in RC
                            )
                        };

                    // Suppress inner-side extensions for each mate.
                    // Inner = 3' end: right for forward (orig_is_rev=false), left for reverse.
                    let Some(mut t1) = finalize_transcript(
                        &m1_wt,
                        m1_read_slice,
                        index,
                        &scorer,
                        &stitch_cluster,
                        m1_orig_rev,
                        m1_orig_rev, // no_left_ext = inner for reverse (orig_is_rev=true)
                        !m1_orig_rev, // no_right_ext = inner for forward (orig_is_rev=false)
                    ) else {
                        continue;
                    };
                    let Some(mut t2) = finalize_transcript(
                        &m2_wt,
                        m2_read_slice,
                        index,
                        &scorer,
                        &stitch_cluster,
                        m2_orig_rev,
                        m2_orig_rev, // no_left_ext = inner for reverse (orig_is_rev=true)
                        !m2_orig_rev, // no_right_ext = inner for forward (orig_is_rev=false)
                    ) else {
                        continue;
                    };

                    if stitch_is_reverse {
                        t1.is_reverse = true;
                        t2.is_reverse = false;
                    } else {
                        t1.is_reverse = false;
                        t2.is_reverse = true;
                    }
                    t1.read_seq = mate1_seq.to_vec();
                    t2.read_seq = mate2_seq.to_vec();

                    if params.chim_segment_min > 0 {
                        all_m1_transcripts.push(t1.clone());
                        all_m2_transcripts.push(t2.clone());
                    }

                    let combined_span =
                        t1.genome_end.max(t2.genome_end) - t1.genome_start.min(t2.genome_start);
                    let combined_wt_score = wt.score + scorer.genomic_length_penalty(combined_span);

                    if let Some(pair) = try_pair_transcripts(
                        &t1,
                        &t2,
                        len1,
                        len2,
                        params,
                        combined_score_threshold,
                        combined_wt_score,
                    ) {
                        joint_pairs.push(pair);
                    }
                }
                None => {
                    // Single-mate WT: save for half-mapped fallback
                    let all_m1 = wt.exons.iter().all(|e| e.mate_id == 0);
                    let all_m2 = wt.exons.iter().all(|e| e.mate_id == 1);
                    if all_m1 {
                        let (read_slice, orig_rev) = if stitch_is_reverse {
                            (&stitch_read[len2 + 1..], true)
                        } else {
                            (&stitch_read[..len1], false)
                        };
                        if let Some(mut t) = finalize_transcript(
                            wt,
                            read_slice,
                            index,
                            &scorer,
                            &stitch_cluster,
                            orig_rev,
                            false,
                            false,
                        ) {
                            t.is_reverse = stitch_is_reverse;
                            t.read_seq = mate1_seq.to_vec();
                            if params.chim_segment_min > 0 {
                                all_m1_transcripts.push(t.clone());
                            }
                            single_mate1_transcripts.push(t);
                        }
                    } else if all_m2 {
                        let (read_slice, orig_rev) = if stitch_is_reverse {
                            (&stitch_read[..len2], false)
                        } else {
                            (&stitch_read[len1 + 1..], true)
                        };
                        if let Some(mut t) = finalize_transcript(
                            wt,
                            read_slice,
                            index,
                            &scorer,
                            &stitch_cluster,
                            orig_rev,
                            false,
                            false,
                        ) {
                            t.is_reverse = !stitch_is_reverse;
                            t.read_seq = mate2_seq.to_vec();
                            if params.chim_segment_min > 0 {
                                all_m2_transcripts.push(t.clone());
                            }
                            single_mate2_transcripts.push(t);
                        }
                    }
                }
            }
        }
    }

    // --- Decision tree: dedup, score-filter, quality-filter, then half-mapped fallback ---

    // Step 1: position dedup — remove exact (chr, mate1_pos, mate2_pos, strand, CIGAR) duplicates.
    // Run dedup BEFORE score-range filter so the backup pool is already deduplicated.
    // (STAR's ordering is multMapSelect → dedup, but dedup before multMapSelect is equivalent
    // since removing exact duplicates doesn't change the best score.)
    joint_pairs.sort_by(|a, b| {
        let pos_cmp = (
            a.mate1_transcript.chr_idx,
            a.mate1_transcript.genome_start,
            a.mate1_transcript.is_reverse,
            a.mate2_transcript.genome_start,
            a.mate2_transcript.is_reverse,
        )
            .cmp(&(
                b.mate1_transcript.chr_idx,
                b.mate1_transcript.genome_start,
                b.mate1_transcript.is_reverse,
                b.mate2_transcript.genome_start,
                b.mate2_transcript.is_reverse,
            ));
        if pos_cmp != std::cmp::Ordering::Equal {
            return pos_cmp;
        }
        b.combined_wt_score
            .cmp(&a.combined_wt_score)
            .then_with(|| a.mate1_transcript.cigar.cmp(&b.mate1_transcript.cigar))
            .then_with(|| a.mate2_transcript.cigar.cmp(&b.mate2_transcript.cigar))
    });
    joint_pairs.dedup_by(|a, b| {
        a.mate1_transcript.chr_idx == b.mate1_transcript.chr_idx
            && a.mate1_transcript.genome_start == b.mate1_transcript.genome_start
            && a.mate1_transcript.is_reverse == b.mate1_transcript.is_reverse
            && a.mate1_transcript.cigar == b.mate1_transcript.cigar
            && a.mate2_transcript.genome_start == b.mate2_transcript.genome_start
            && a.mate2_transcript.is_reverse == b.mate2_transcript.is_reverse
            && a.mate2_transcript.cigar == b.mate2_transcript.cigar
    });

    // Post-finalization mate2-exon-subset dedup.
    {
        use crate::align::transcript::Exon;
        let exons_subset = |b: &[Exon], a: &[Exon]| -> bool {
            let total_b: u32 = b.iter().map(|e| (e.read_end - e.read_start) as u32).sum();
            if total_b == 0 {
                return false;
            }
            let mut covered = 0u32;
            for be in b {
                let b_diag = be.genome_start as i64 - be.read_start as i64;
                for ae in a {
                    let a_diag = ae.genome_start as i64 - ae.read_start as i64;
                    if a_diag == b_diag {
                        let r_start = be.read_start.max(ae.read_start);
                        let r_end = be.read_end.min(ae.read_end);
                        if r_start < r_end {
                            covered += (r_end - r_start) as u32;
                        }
                    }
                }
            }
            covered == total_b
        };
        let n = joint_pairs.len();
        let mut keep = vec![true; n];
        for i in 0..n {
            if !keep[i] {
                continue;
            }
            for j in 0..n {
                if i == j || !keep[j] {
                    continue;
                }
                let same_pos = joint_pairs[i].mate1_transcript.chr_idx
                    == joint_pairs[j].mate1_transcript.chr_idx
                    && joint_pairs[i].mate1_transcript.genome_start
                        == joint_pairs[j].mate1_transcript.genome_start
                    && joint_pairs[i].mate1_transcript.genome_end
                        == joint_pairs[j].mate1_transcript.genome_end
                    && joint_pairs[i].mate1_transcript.is_reverse
                        == joint_pairs[j].mate1_transcript.is_reverse
                    && joint_pairs[i].mate2_transcript.genome_start
                        == joint_pairs[j].mate2_transcript.genome_start
                    && joint_pairs[i].mate2_transcript.is_reverse
                        == joint_pairs[j].mate2_transcript.is_reverse;
                if same_pos
                    && joint_pairs[i].combined_wt_score < joint_pairs[j].combined_wt_score
                    && exons_subset(
                        &joint_pairs[i].mate2_transcript.exons,
                        &joint_pairs[j].mate2_transcript.exons,
                    )
                {
                    keep[i] = false;
                    break;
                }
            }
        }
        joint_pairs = joint_pairs
            .into_iter()
            .enumerate()
            .filter(|(i, _)| keep[*i])
            .map(|(_, p)| p)
            .collect();
    }

    // Step 2: score-range filter (STAR's multMapSelect).
    if !joint_pairs.is_empty() {
        let best_score = joint_pairs
            .iter()
            .map(|pa| pa.combined_wt_score)
            .max()
            .unwrap_or(0);
        let score_threshold = best_score - params.out_filter_multimap_score_range;
        joint_pairs.retain(|pa| pa.combined_wt_score >= score_threshold);
    }

    // Step 3: TooManyLoci check (post-dedup, matching STAR's ordering: multMapSelect → dedup → TooManyLoci).
    if joint_pairs.len() > params.out_filter_multimap_nmax as usize {
        return Ok((
            Vec::new(),
            pe_chimeric,
            0,
            Some(UnmappedReason::TooManyLoci),
        ));
    }

    joint_pairs.sort_by(|a, b| {
        b.combined_wt_score
            .cmp(&a.combined_wt_score)
            .then_with(|| a.mate1_transcript.chr_idx.cmp(&b.mate1_transcript.chr_idx))
            .then_with(|| {
                a.mate1_transcript
                    .genome_start
                    .cmp(&b.mate1_transcript.genome_start)
            })
            .then_with(|| {
                a.mate1_transcript
                    .is_reverse
                    .cmp(&b.mate1_transcript.is_reverse)
            })
    });

    // Randomize primary among best-scoring pairs (STAR's funPrimaryAlignMark).
    shuffle_tied_prefix(
        &mut joint_pairs,
        |pa| pa.combined_wt_score,
        per_read_seed(params.run_rng_seed, read_name),
    );

    // Step 4: quality filter (mappedFilter).
    filter_paired_transcripts(&mut joint_pairs, params);

    // PE Tier 1: chimericDetectionOld per-mate — mirrors SE behavior but run independently
    // on each mate's transcript pool (joint-pair halves + single-mate WTs combined).
    // Runs before the BothMapped early return so chimeras are reported for all pair outcomes.
    if params.chim_segment_min > 0 {
        use crate::chimeric::detect_chimeric_old;
        if let Some(tr_best_m1) = all_m1_transcripts.iter().max_by_key(|t| t.score) {
            let chims = detect_chimeric_old(
                &all_m1_transcripts,
                tr_best_m1,
                mate1_seq,
                read_name,
                params,
                index,
            )?;
            pe_chimeric.extend(chims);
        }
        if let Some(tr_best_m2) = all_m2_transcripts.iter().max_by_key(|t| t.score) {
            let chims = detect_chimeric_old(
                &all_m2_transcripts,
                tr_best_m2,
                mate2_seq,
                read_name,
                params,
                index,
            )?;
            pe_chimeric.extend(chims);
        }
        pe_chimeric.retain(|chim| {
            chim.meets_min_segment_length(params.chim_segment_min)
                && chim.meets_min_score(params.chim_score_min)
        });
    }

    if !joint_pairs.is_empty() {
        let pe_mapq_n = joint_pairs.len().max(1);
        let results = joint_pairs
            .into_iter()
            .map(|pa| PairedAlignmentResult::BothMapped(Box::new(pa)))
            .collect();
        return Ok((results, pe_chimeric, pe_mapq_n, None));
    }

    // Inter-mate chimeric detection: fires when the best single-mate transcripts are discordant
    // (different chr, same strand, or >1Mb apart). Runs before half-mapped fallback consumes
    // the transcript vecs.
    if params.chim_segment_min > 0 {
        use crate::chimeric::detect_inter_mate_chimeric;
        let best_m1_chim = single_mate1_transcripts.iter().max_by_key(|t| t.score);
        let best_m2_chim = single_mate2_transcripts.iter().max_by_key(|t| t.score);
        if let (Some(t1), Some(t2)) = (best_m1_chim, best_m2_chim)
            && let Some(chim) =
                detect_inter_mate_chimeric(t1, t2, mate1_seq, read_name, params, index)
        {
            pe_chimeric.push(chim);
        }
    }

    // Half-mapped fallback: report the best-scoring single-mate transcript.
    // STAR applies the quality filter to the COMBINED read (Lread-1 = len1+len2), so we
    // use the same threshold for each mate here.
    let single_mate_threshold = combined_score_threshold.max(params.out_filter_score_min);

    let best_m1 = single_mate1_transcripts
        .into_iter()
        .filter(|t| t.score >= single_mate_threshold)
        .max_by_key(|t| t.score);
    let best_m2 = single_mate2_transcripts
        .into_iter()
        .filter(|t| t.score >= single_mate_threshold)
        .max_by_key(|t| t.score);

    match (best_m1, best_m2) {
        (Some(t1), None) => Ok((
            vec![PairedAlignmentResult::HalfMapped {
                mapped_transcript: t1,
                mate1_is_mapped: true,
            }],
            pe_chimeric,
            1,
            None,
        )),
        (None, Some(t2)) => Ok((
            vec![PairedAlignmentResult::HalfMapped {
                mapped_transcript: t2,
                mate1_is_mapped: false,
            }],
            pe_chimeric,
            1,
            None,
        )),
        (Some(t1), Some(t2)) => {
            // Both have single-mate alignments but couldn't form a valid pair.
            // Report the higher-scoring mate as half-mapped.
            if t1.score >= t2.score {
                Ok((
                    vec![PairedAlignmentResult::HalfMapped {
                        mapped_transcript: t1,
                        mate1_is_mapped: true,
                    }],
                    pe_chimeric,
                    1,
                    None,
                ))
            } else {
                Ok((
                    vec![PairedAlignmentResult::HalfMapped {
                        mapped_transcript: t2,
                        mate1_is_mapped: false,
                    }],
                    pe_chimeric,
                    1,
                    None,
                ))
            }
        }
        (None, None) => Ok((Vec::new(), pe_chimeric, 0, Some(UnmappedReason::TooShort))),
    }
}

/// Attempt to pair two per-mate transcripts into a PairedAlignment.
///
/// Returns `None` if the mates are incompatible (same strand, different chr, too far, etc.).
#[allow(clippy::too_many_arguments)]
fn try_pair_transcripts(
    t1: &Transcript,
    t2: &Transcript,
    len1: usize,
    len2: usize,
    params: &Parameters,
    combined_score_threshold: i32,
    combined_wt_score: i32,
) -> Option<PairedAlignment> {
    // Must be same chromosome
    if t1.chr_idx != t2.chr_idx {
        return None;
    }
    // Must be opposite strands (FR or RF)
    if t1.is_reverse == t2.is_reverse {
        return None;
    }

    // Determine left mate (smaller genome_start) and right mate for distance/consistency checks
    let (left, right) = if t1.genome_start <= t2.genome_start {
        (t1, t2)
    } else {
        (t2, t1)
    };

    // Reject degenerate pairs (right ends before left starts → negative insert)
    if right.genome_end <= left.genome_start {
        return None;
    }

    // Genomic span check: use alignMatesGapMax if set, else fall back to win_bin_window_dist
    // (STAR's effective limit when alignMatesGapMax=0 is the window distance ~589kb)
    let span = right.genome_end - left.genome_start;
    let max_span = if params.align_mates_gap_max > 0 {
        params.align_mates_gap_max as u64
    } else {
        params.win_bin_window_dist()
    };
    if span > max_span {
        return None;
    }

    // SCORE-GATE: reject pairs where score is below the absolute floor
    if combined_wt_score + params.out_filter_multimap_score_range < combined_score_threshold {
        return None;
    }

    // Junction consistency in overlap region
    if !pe_junctions_consistent(left, right) {
        return None;
    }

    // Combined coverage: sum of exon read spans from both mates
    let combined_n_match: u32 = t1
        .exons
        .iter()
        .map(|e| (e.read_end - e.read_start) as u32)
        .sum::<u32>()
        + t2.exons
            .iter()
            .map(|e| (e.read_end - e.read_start) as u32)
            .sum::<u32>();

    let is_proper_pair = check_proper_pair(t1, t2, params);
    let insert_size = calculate_insert_size(t1, t2);

    Some(PairedAlignment {
        mate1_transcript: t1.clone(),
        mate2_transcript: t2.clone(),
        mate1_region: (0, len1),
        mate2_region: (0, len2),
        is_proper_pair,
        insert_size,
        combined_wt_score,
        combined_n_match,
    })
}

/// Check if paired alignment is a proper pair
fn check_proper_pair(
    mate1_trans: &Transcript,
    mate2_trans: &Transcript,
    params: &Parameters,
) -> bool {
    // Proper pair criteria:
    // 1. Both mates mapped (checked by caller)
    // 2. Same chromosome (checked by caller)
    // 3. Distance within alignMatesGapMax

    if params.align_mates_gap_max == 0 {
        return true; // Auto mode = unlimited
    }

    // Calculate genomic distance
    let start = mate1_trans.genome_start.min(mate2_trans.genome_start);
    let end = mate1_trans.genome_end.max(mate2_trans.genome_end);
    let genomic_span = end - start;

    genomic_span <= params.align_mates_gap_max as u64
}

/// Calculate signed insert size (TLEN)
fn calculate_insert_size(mate1_trans: &Transcript, mate2_trans: &Transcript) -> i32 {
    // STAR outSAMtlen=1 (default): tlen is computed from the combined PE transcript span,
    // not from max/min of individual mate endpoints.
    //
    // Forward cluster (mate1 forward, mate2 reverse — FR pair):
    //   tlen = mate2.genome_end - mate1.genome_start
    //   mate1 (imate=0 in combined transcript) gets +tlen
    //
    // Reverse cluster (mate1 reverse, mate2 forward — RF pair):
    //   tlen = mate1.genome_end - mate2.genome_start
    //   mate2 (imate=0 in combined transcript) gets +tlen, so mate1 gets -tlen
    //
    // This correctly handles:
    //   - Same-start overlapping pairs (uses trailing mate's end, not max of both)
    //   - RF pairs at same position (mate2 gets +tlen, not mate1)
    if !mate1_trans.is_reverse {
        // Forward cluster: mate1 on left
        (mate2_trans.genome_end - mate1_trans.genome_start) as i32
    } else {
        // Reverse cluster: mate2 on left, so mate1 gets negative tlen
        -((mate1_trans.genome_end - mate2_trans.genome_start) as i32)
    }
}

/// Filter paired transcripts by quality thresholds.
/// STAR's mappedFilter applies ALL quality checks to trBest (the highest-scoring transcript).
fn filter_paired_transcripts(paired_alns: &mut Vec<PairedAlignment>, params: &Parameters) {
    // STAR's mappedFilter (ReadAlign_mappedFilter.cpp) applies ALL quality thresholds to trBest
    // (the highest-scoring transcript), NOT to each individual transcript. If trBest passes,
    // all transcripts in the score window are included (they affect NH/MAPQ). If trBest fails,
    // the read is unmapped.
    //
    // Step 1: find the best pair and check quality thresholds on it.
    let best_pa = paired_alns.iter().max_by_key(|pa| pa.combined_wt_score);
    if let Some(best) = best_pa {
        let mate1_len = (best.mate1_region.1 - best.mate1_region.0) as f64;
        let mate2_len = (best.mate2_region.1 - best.mate2_region.0) as f64;
        let combined_lread_m1 = mate1_len + mate2_len;
        let combined_nm = best.mate1_transcript.n_mismatch + best.mate2_transcript.n_mismatch;
        let combined_score = best.combined_wt_score;

        if combined_score < params.out_filter_score_min
            || combined_score < (params.out_filter_score_min_over_lread * combined_lread_m1) as i32
        {
            paired_alns.clear();
            return;
        }

        let combined_match = best.combined_n_match;
        if combined_match < params.out_filter_match_nmin
            || combined_match < (params.out_filter_match_nmin_over_lread * combined_lread_m1) as u32
        {
            paired_alns.clear();
            return;
        }

        if combined_nm > params.out_filter_mismatch_nmax
            || (combined_nm as f64)
                > params.out_filter_mismatch_nover_lmax * (mate1_len + mate2_len)
        {
            paired_alns.clear();
            return;
        }
    }

    if paired_alns.len() > params.out_filter_multimap_nmax as usize {
        paired_alns.clear();
    }
}

/// Extract splice junctions from a Transcript's CIGAR as (donor, acceptor) pairs.
/// Junction coords are in genomic space (0-based). Junctions only exist where CigarOp::RefSkip is.
fn extract_junctions_from_cigar(t: &Transcript) -> Vec<(u64, u64)> {
    let mut junctions = Vec::new();
    let mut genome_pos = t.genome_start;
    for op in &t.cigar {
        use crate::align::transcript::CigarOp;
        match op {
            CigarOp::Match(n) | CigarOp::Equal(n) | CigarOp::Diff(n) | CigarOp::Del(n) => {
                genome_pos += *n as u64;
            }
            CigarOp::RefSkip(n) => {
                let donor = genome_pos;
                let acceptor = genome_pos + *n as u64;
                junctions.push((donor, acceptor));
                genome_pos = acceptor;
            }
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
        }
    }
    junctions
}

/// D5: Check junction consistency in the overlapping region of paired-end mates.
/// When mates overlap in the genome, every splice junction in the overlap from the
/// left mate must appear in the right mate, and vice versa.
/// Implements STAR stitchWindowAligns.cpp check after overlap detection.
///
/// `left` is the mate with the lower genome_start, `right` is the other mate.
pub(crate) fn pe_junctions_consistent(left: &Transcript, right: &Transcript) -> bool {
    // Overlapping region: [overlap_start, overlap_end)
    let overlap_start = left.genome_start.max(right.genome_start);
    let overlap_end = left.genome_end.min(right.genome_end);
    if overlap_start >= overlap_end {
        return true; // No overlap — nothing to check
    }

    let left_juncs = extract_junctions_from_cigar(left);
    let right_juncs = extract_junctions_from_cigar(right);

    // Every junction from left that falls within the overlap must be in right too
    for (donor, acceptor) in &left_juncs {
        if *donor >= overlap_start
            && *acceptor <= overlap_end
            && !right_juncs.iter().any(|(d, a)| d == donor && a == acceptor)
        {
            return false;
        }
    }
    // Every junction from right that falls within the overlap must be in left too
    for (donor, acceptor) in &right_juncs {
        if *donor >= overlap_start
            && *acceptor <= overlap_end
            && !left_juncs.iter().any(|(d, a)| d == donor && a == acceptor)
        {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;
    use clap::Parser;

    fn make_test_params() -> Parameters {
        // Parse empty args to get default parameters
        Parameters::try_parse_from(vec!["rustar-aligner"]).unwrap()
    }

    fn make_test_index() -> GenomeIndex {
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

        // Create dummy SA and SAindex (would need real index for actual alignment)
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
    fn combined_transcript_for_projection_rewrites_mate2_ifrag() {
        use crate::align::transcript::CigarOp;

        let make_tr = |gs: u64, ge: u64, rs: usize, re: usize| Transcript {
            chr_idx: 0,
            genome_start: gs,
            genome_end: ge,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: gs,
                genome_end: ge,
                read_start: rs,
                read_end: re,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match((ge - gs) as u32)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };
        let pair = PairedAlignment {
            mate1_transcript: make_tr(1000, 1100, 0, 100),
            mate2_transcript: make_tr(1300, 1400, 0, 100),
            mate1_region: (0, 100),
            mate2_region: (0, 100),
            is_proper_pair: true,
            insert_size: 400,
            combined_wt_score: 200,
            combined_n_match: 200,
        };
        let combined = pair.combined_transcript_for_projection();
        assert_eq!(combined.exons.len(), 2);
        assert_eq!(combined.exons[0].i_frag, 0);
        assert_eq!(combined.exons[1].i_frag, 1);
        assert_eq!(combined.genome_start, 1000);
        assert_eq!(combined.genome_end, 1400);
        assert_eq!(combined.score, 200);
    }

    #[test]
    fn test_align_read_no_seeds() {
        let index = make_test_index();
        let params = make_test_params();

        // Read with all N's (no seeds possible)
        let read_seq = vec![4, 4, 4, 4, 4, 4, 4, 4, 4, 4];

        let result = align_read(&read_seq, "READ_001", &index, &params);
        assert!(result.is_ok());

        let (transcripts, chimeras, n_for_mapq, unmapped_reason) = result.unwrap();
        assert_eq!(transcripts.len(), 0); // No alignment
        assert_eq!(chimeras.len(), 0); // No chimeric alignments
        assert_eq!(n_for_mapq, 0);
        assert_eq!(unmapped_reason, Some(UnmappedReason::Other));
    }

    #[test]
    fn test_transcript_filtering_score() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_score_min = 50;

        // Would need actual seeds and alignment to test this properly
        // This test just verifies the function doesn't crash
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_002", &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_filtering_mismatch() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_mismatch_nmax = 2;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_003", &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_multimap_limit() {
        let index = make_test_index();
        let mut params = make_test_params();
        params.out_filter_multimap_nmax = 5;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_004", &index, &params);
        assert!(result.is_ok());

        let (transcripts, _chimeras, _n_for_mapq, _reason) = result.unwrap();
        assert!(transcripts.len() <= 5);
    }

    #[test]
    fn test_align_paired_read_no_seeds() {
        let index = make_test_index();
        let params = make_test_params();

        // Both mates with all N's
        let mate1 = vec![4, 4, 4, 4, 4, 4, 4, 4];
        let mate2 = vec![4, 4, 4, 4, 4, 4, 4, 4];

        let result = align_paired_read(&mate1, &mate2, "test", &index, &params);
        assert!(result.is_ok());
        let (paired_alns, _chimeric, n_for_mapq, unmapped_reason) = result.unwrap();
        assert_eq!(paired_alns.len(), 0);
        assert_eq!(n_for_mapq, 0);
        assert!(unmapped_reason.is_some());
    }

    #[test]
    fn test_check_proper_pair_distance() {
        use crate::align::transcript::{CigarOp, Exon};

        let params = make_test_params();

        // Create two transcripts on same chromosome
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1200,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Distance = 300bp, within default limit (auto mode = unlimited)
        assert!(check_proper_pair(&t1, &t2, &params));
    }

    #[test]
    fn test_check_proper_pair_too_far() {
        use crate::align::transcript::{CigarOp, Exon};

        let mut params = make_test_params();
        params.align_mates_gap_max = 100;

        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1300,
            genome_end: 1400,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1300,
                genome_end: 1400,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Distance = 400bp, exceeds limit of 100bp
        assert!(!check_proper_pair(&t1, &t2, &params));
    }

    #[test]
    fn test_calculate_insert_size_positive() {
        use crate::align::transcript::{CigarOp, Exon};

        // Mate1 is leftmost
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 1200,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, 300); // Positive because mate1 is leftmost
    }

    #[test]
    fn test_strand_consistency_filter() {
        use crate::align::transcript::{CigarOp, Exon, Transcript};
        use crate::params::IntronStrandFilter;

        // Create a transcript with conflicting strand motifs (mixed + and - within one transcript)
        let t_inconsistent = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc], // +strand and -strand
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Create a transcript with consistent strand motifs (all + strand)
        let t_consistent = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::GcAg], // both + strand
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Note: STAR's RemoveInconsistentStrands filters transcripts where
        // junctions have MIXED implied strand (both + and - within one transcript).
        // It does NOT compare junction strand vs alignment strand — a reverse-strand
        // read at a GT/AG junction (antisense of + strand gene) is valid and kept.

        // Verify mixed-strand is detected
        let mut has_plus = false;
        let mut has_minus = false;
        for motif in &t_inconsistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && has_minus); // Inconsistent (mixed strands)

        // Verify consistent transcript has no conflict
        has_plus = false;
        has_minus = false;
        for motif in &t_consistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && !has_minus); // Consistent (all +)

        // Also verify: single CT/AC on forward-strand or single GT/AG on reverse-strand
        // are NOT filtered (STAR keeps these — antisense reads are valid)
        let ctac_only = vec![SpliceMotif::CtAc];
        let (mut hp, mut hm) = (false, false);
        for m in &ctac_only {
            match m.implied_strand() {
                Some('+') => hp = true,
                Some('-') => hm = true,
                _ => {}
            }
        }
        assert!(!hp && hm); // Only minus → NOT mixed → NOT filtered

        // Verify the filter enum
        assert_ne!(
            IntronStrandFilter::None,
            IntronStrandFilter::RemoveInconsistentStrands
        );
    }

    #[test]
    fn test_calculate_insert_size_negative() {
        use crate::align::transcript::{CigarOp, Exon};

        // RF pair: mate1 reverse (right side), mate2 forward (left side).
        // STAR reverse cluster: tlen = mate1.genome_end - mate2.genome_start = 1300 - 1000 = 300.
        // mate1 gets -tlen (negative) because mate2 (imate=0 in reverse cluster) is the left mate.
        let t1 = Transcript {
            chr_idx: 0,
            genome_start: 1200,
            genome_end: 1300,
            is_reverse: true, // mate1 is reverse strand → reverse cluster
            exons: vec![Exon {
                genome_start: 1200,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let t2 = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, -300); // Negative for mate1 in RF pair (mate2 is the left mate)
    }

    #[test]
    fn test_noncanonical_unannotated_filter() {
        use crate::align::score::SpliceMotif;
        use crate::align::transcript::{CigarOp, Exon, Transcript};

        // Helper: check if a transcript would be filtered by RemoveNoncanonicalUnannotated
        // (mirrors the logic in the retain closure)
        let would_filter = |t: &Transcript| -> bool {
            t.junction_motifs
                .iter()
                .zip(t.junction_annotated.iter())
                .any(|(m, annotated)| *m == SpliceMotif::NonCanonical && !annotated)
        };

        let base_transcript = || Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1200,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1200,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Case 1: NonCanonical + unannotated → should be filtered
        let mut t1 = base_transcript();
        t1.junction_motifs = vec![SpliceMotif::NonCanonical];
        t1.junction_annotated = vec![false];
        assert!(
            would_filter(&t1),
            "NonCanonical + unannotated should be filtered"
        );

        // Case 2: NonCanonical + annotated → should be KEPT
        let mut t2 = base_transcript();
        t2.junction_motifs = vec![SpliceMotif::NonCanonical];
        t2.junction_annotated = vec![true];
        assert!(
            !would_filter(&t2),
            "NonCanonical + annotated should be kept"
        );

        // Case 3: Canonical + unannotated → should be KEPT
        let mut t3 = base_transcript();
        t3.junction_motifs = vec![SpliceMotif::GtAg];
        t3.junction_annotated = vec![false];
        assert!(!would_filter(&t3), "Canonical + unannotated should be kept");

        // Case 4: Mixed — one canonical + one non-canonical unannotated → filtered
        let mut t4 = base_transcript();
        t4.junction_motifs = vec![SpliceMotif::GtAg, SpliceMotif::NonCanonical];
        t4.junction_annotated = vec![true, false];
        assert!(
            would_filter(&t4),
            "Mixed with unannotated non-canonical should be filtered"
        );

        // Case 5: Mixed — one canonical + one non-canonical annotated → kept
        let mut t5 = base_transcript();
        t5.junction_motifs = vec![SpliceMotif::GtAg, SpliceMotif::NonCanonical];
        t5.junction_annotated = vec![false, true];
        assert!(
            !would_filter(&t5),
            "Mixed with annotated non-canonical should be kept"
        );
    }

    #[test]
    fn test_align_paired_both_unmapped() {
        // Both mates are all N's → both unmapped → empty Vec
        let index = make_test_index();
        let params = make_test_params();

        let mate1 = vec![4, 4, 4, 4, 4, 4, 4, 4];
        let mate2 = vec![4, 4, 4, 4, 4, 4, 4, 4];

        let (results, _chimeric, n_for_mapq, unmapped_reason) =
            align_paired_read(&mate1, &mate2, "test", &index, &params).unwrap();
        assert!(results.is_empty(), "Both unmapped should return empty Vec");
        assert_eq!(n_for_mapq, 0);
        assert!(unmapped_reason.is_some(), "Should have unmapped reason");
    }

    #[test]
    fn test_paired_alignment_result_enum_variants() {
        use crate::align::transcript::{CigarOp, Exon};

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1100,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1100,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Test BothMapped variant
        let both = PairedAlignmentResult::BothMapped(Box::new(PairedAlignment {
            mate1_transcript: transcript.clone(),
            mate2_transcript: transcript.clone(),
            mate1_region: (0, 100),
            mate2_region: (0, 100),
            is_proper_pair: true,
            insert_size: 200,
            combined_wt_score: 0,
            combined_n_match: 200,
        }));
        assert!(matches!(both, PairedAlignmentResult::BothMapped(_)));

        // Test HalfMapped variant
        let half = PairedAlignmentResult::HalfMapped {
            mapped_transcript: transcript,
            mate1_is_mapped: true,
        };
        assert!(matches!(half, PairedAlignmentResult::HalfMapped { .. }));

        // Verify mate1_is_mapped
        if let PairedAlignmentResult::HalfMapped {
            mate1_is_mapped, ..
        } = half
        {
            assert!(mate1_is_mapped);
        }
    }

    #[test]
    fn shuffle_tied_prefix_is_deterministic() {
        // Same seed + same input → same permutation on reruns.
        let items: Vec<(i32, u32)> = (0..8).map(|i| (100, i)).collect();
        let mut a = items.clone();
        let mut b = items.clone();
        shuffle_tied_prefix(&mut a, |t| t.0, 12345);
        shuffle_tied_prefix(&mut b, |t| t.0, 12345);
        assert_eq!(a, b);
    }

    #[test]
    fn shuffle_tied_prefix_respects_ties() {
        // Only the top-score prefix gets shuffled; lower-scored tail is left alone.
        let mut items = vec![(100, 0u32), (100, 1), (100, 2), (50, 3), (40, 4)];
        shuffle_tied_prefix(&mut items, |t| t.0, 777);
        // Last two elements (non-tied) stay in place.
        assert_eq!(items[3], (50, 3));
        assert_eq!(items[4], (40, 4));
        // Tied prefix contains the original three items in some order.
        let mut top: Vec<u32> = items[..3].iter().map(|t| t.1).collect();
        top.sort();
        assert_eq!(top, vec![0, 1, 2]);
    }

    #[test]
    fn shuffle_tied_prefix_different_seeds_can_diverge() {
        // Probabilistic: for a tied set of 8, at least two seeds should disagree
        // on the chosen primary. (Exhaustive over a small seed range is fine.)
        let base: Vec<(i32, u32)> = (0..8).map(|i| (100, i)).collect();
        let mut firsts = std::collections::HashSet::new();
        for seed in 0..32u64 {
            let mut v = base.clone();
            shuffle_tied_prefix(&mut v, |t| t.0, seed);
            firsts.insert(v[0].1);
        }
        assert!(
            firsts.len() >= 2,
            "expected different seeds to pick different primaries, got {:?}",
            firsts
        );
    }

    #[test]
    fn shuffle_tied_prefix_noop_when_no_ties() {
        let mut items = vec![(100, 0u32), (90, 1), (80, 2)];
        let before = items.clone();
        shuffle_tied_prefix(&mut items, |t| t.0, 42);
        assert_eq!(items, before);
    }
}
