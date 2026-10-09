/// Read alignment driver function
use crate::align::score::AlignmentScorer;
use crate::align::seed::Seed;
use crate::align::stitch::{
    PE_SPACER_BASE, cluster_seeds, finalize_transcript, split_combined_wt, stitch_seeds_core,
    stitch_seeds_with_jdb_debug,
};
use crate::align::transcript::{Exon, Transcript};
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::{MultimapperOrder, Parameters};
use crate::stats::UnmappedReason;
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
    crate::rng::shuffle_deterministic(&mut items[..tied], seed);
}

/// STAR's `nMatch`: aligned bases where the read equals the genome. Mismatches
/// and positions with an `N` on either side are not counted
/// (`extendAlign.cpp:35-41`, `stitchAlignToTranscript.cpp:68-87`), so this is
/// smaller than the aligned length whenever the alignment has a mismatch.
///
/// `read` is the read as it was sequenced; a reverse-strand transcript's CIGAR
/// walks its reverse complement.
fn star_n_match(t: &Transcript, read: &[u8], genome: &crate::genome::Genome) -> u32 {
    use noodles::sam::alignment::record::cigar::op::Kind;
    // Base `i` of the frame the CIGAR walks, without building the reverse
    // complement: this runs once per read.
    let frame = |i: usize| -> u8 {
        let b = if t.is_reverse {
            read.len().checked_sub(i + 1).and_then(|j| read.get(j))
        } else {
            read.get(i)
        };
        match b {
            Some(&b) if t.is_reverse && b < 4 => 3 - b,
            Some(&b) => b,
            None => 4,
        }
    };
    let (mut r, mut g) = (0usize, t.genome_start);
    let mut n = 0u32;
    for op in &t.cigar {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for k in 0..len {
                    let rb = frame(r + k);
                    let gb = genome.get_base(g + k as u64).unwrap_or(4);
                    if rb < 4 && gb < 4 && rb == gb {
                        n += 1;
                    }
                }
                r += len;
                g += len as u64;
            }
            Kind::Insertion | Kind::SoftClip => r += len,
            Kind::Deletion | Kind::Skip => g += len as u64,
            _ => {}
        }
    }
    n
}

/// STAR's `rLength`: the sum of exon lengths, i.e. the aligned read bases
/// including mismatches. The denominator of `--outFilterMismatchNoverLmax`.
fn star_r_length(t: &Transcript) -> u32 {
    t.n_matched() as u32
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
    align_read_inner(read_seq, read_name, index, params, true)
}

/// `align_read` with or without STAR's `mappedFilter`. The `--peOverlapNbasesMin`
/// merged read is mapped without it: STAR runs only `mapOneRead` on it
/// (`ReadAlign_peOverlapMergeMap.cpp`) and filters once the alignments are back
/// in the paired frame.
fn align_read_inner(
    read_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
    mapped_filter: bool,
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
            eprintln!("[DEBUG {read_name}] No seeds found — unmapped");
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
            eprintln!("[DEBUG {read_name}] No clusters — unmapped");
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Cap total clusters (alignWindowsPerReadNmax)
    let mut clusters = clusters;
    clusters.truncate(params.align_windows_per_read_nmax);

    // NOTE: STAR's winReadCoverageRelativeMin filter is long-reads-only
    // (#ifdef COMPILE_FOR_LONG_READS in stitchPieces.cpp). Standard STAR
    // does NOT filter clusters by seed coverage. Removed to match STAR.

    // Step 3: Stitch seeds within each cluster
    let scorer = AlignmentScorer::from_params(params);
    let mut transcripts = Vec::new();
    // STAR's `trAll`: every window's transcripts, before any cross-window
    // dedup or filtering, which is what its chimeric detection reads.
    let mut chim_windows: Vec<Vec<crate::chimeric::WinTr>> = Vec::new();

    // Use junction DB for annotation-aware scoring if available
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };

    for (ci, cluster) in clusters.iter().enumerate() {
        // STAR's `alignTranscriptsPerReadNmax` headroom break
        // (`ReadAlign_stitchPieces.cpp`): stop collecting windows once one more
        // window's worth of transcripts could overrun the per-read cap. The
        // test is on the headroom, not on the running total, so the cap is
        // never exceeded rather than merely noticed afterwards.
        if transcripts.len() + params.align_transcripts_per_window_nmax
            >= params.align_transcripts_per_read_nmax
        {
            break;
        }
        let debug_name = if debug_read { read_name } else { "" };
        let cluster_transcripts = stitch_seeds_with_jdb_debug(
            cluster,
            read_seq,
            index,
            &scorer,
            junction_db,
            params.align_transcripts_per_window_nmax,
            debug_name,
        );
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
                    t.cigar_string()
                );
            }
        }
        if params.chim_segment_min > 0 {
            chim_windows.push(
                cluster_transcripts
                    .iter()
                    .filter_map(|t| {
                        crate::chimeric::WinTr::single(
                            t,
                            u8::from(t.is_reverse),
                            0,
                            0,
                            read_seq.len(),
                        )
                    })
                    .collect(),
            );
        }
        transcripts.extend(cluster_transcripts);
    }

    // Step 4a: Deduplicate and score-range filter — BEFORE quality filters.
    // STAR order: multMapSelect (score-range) → mappedFilter (quality gates).
    // Doing quality filters first is wrong: it can remove the high-scoring primary,
    // leaving a lower-scoring secondary that then passes as the "best" alignment.

    // Deduplicate transcripts with identical genomic coordinates AND CIGAR,
    // PRESERVING discovery (window) order. STAR's primary tie-break falls back
    // to the earliest-discovered alignment, so we must not reorder here (the
    // old code sorted by coordinate, destroying that order).
    {
        let mut seen = std::collections::HashSet::new();
        transcripts.retain(|t| {
            seen.insert((
                t.chr_idx,
                t.genome_start,
                t.genome_end,
                t.is_reverse,
                t.cigar_string(),
            ))
        });
    }

    // Deterministic primary tie-break. STAR compares `maxScore`, then
    // `gLength` — the alignment's genomic span — and leaves anything still
    // tied to whichever window it reached first
    // (`ReadAlign_stitchPieces.cpp:340`). The span is reproducible here; the
    // window order is not, since windows are built per read in parallel, so
    // the remaining keys are positional and fixed (see DIVERGENCE.md §1.1).
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| (a.genome_end - a.genome_start).cmp(&(b.genome_end - b.genome_start)))
            .then_with(|| a.chr_idx.cmp(&b.chr_idx))
            .then_with(|| a.genome_start.cmp(&b.genome_start))
            .then_with(|| a.is_reverse.cmp(&b.is_reverse))
    });

    // Primary selection — STAR's multMapSelect (ReadAlign_multMapSelect.cpp).
    //
    // STAR's DEFAULT (`--outMultimapperOrder Old_2.4`) does NOT consult the RNG
    // for primary selection; it marks the deterministic best alignment primary.
    // Only `--outMultimapperOrder Random` shuffles. Previously rustar-aligner
    // shuffled unconditionally, which randomised the primary among equal-score
    // loci and diverged from STAR's deterministic choice. Gate the shuffle on
    // the Random mode so the default is deterministic and STAR-faithful.
    //
    // Under Random, we shuffle the tied top-score prefix with a per-read seed
    // (deterministic per read → thread-count invariant).
    if params.out_multimapper_order == MultimapperOrder::Random {
        shuffle_tied_prefix(
            &mut transcripts,
            |t| t.score,
            per_read_seed(params.run_rng_seed, read_name),
        );
    }

    // Score-range filter: keep only alignments within outFilterMultimapScoreRange of the best.
    // (STAR's multMapSelect step — must run before quality filters.)
    if !transcripts.is_empty() {
        let max_score = transcripts[0].score;
        let score_threshold = max_score - params.out_filter_multimap_score_range;
        transcripts.retain(|t| t.score >= score_threshold);
    }

    // Step 4: STAR's mappedFilter (`ReadAlign_mappedFilter.cpp`), which runs after
    // multMapSelect. It looks at the best transcript only and decides for the whole
    // read: too short (score or matched bases), then too many mismatches, then too
    // many loci. The intron filters are not here: STAR applies them per transcript
    // at stitch time, which is where they now run.
    let lread_m1 = (read_seq.len() as f64) - 1.0;
    let unmapped_by_filter = transcripts
        .first()
        .filter(|_| mapped_filter)
        .and_then(|best| {
            let n_match = star_n_match(best, read_seq, &index.genome);
            let r_length = star_r_length(best);
            if best.score < params.out_filter_score_min
                || best.score < (params.out_filter_score_min_over_lread * lread_m1) as i32
                || n_match < params.out_filter_match_nmin
                || n_match < (params.out_filter_match_nmin_over_lread * lread_m1) as u32
            {
                Some(UnmappedReason::TooShort)
            } else if best.n_mismatch > params.out_filter_mismatch_nmax
                || f64::from(best.n_mismatch) / f64::from(r_length.max(1))
                    > params.out_filter_mismatch_nover_lmax
            {
                Some(UnmappedReason::TooManyMismatches)
            } else {
                None
            }
        });
    if debug_read {
        eprintln!("[DEBUG {read_name}] mappedFilter on trBest: {unmapped_by_filter:?}");
    }

    // Note: STAR sometimes finds 2 equivalent indel placements in homopolymer runs
    // via its recursive stitcher's seed exploration (NH=2 instead of NH=1 for ~5 reads).
    // Generating equivalents post-hoc causes more harm than good (41 false NH=2 vs 5 fixed).
    // The root cause is jR scanning placing insertions at different positions — fixing that
    // would be a better approach than post-hoc enumeration.

    // STAR's chimericDetection, which runs after multMapSelect/mappedFilter
    // whatever their outcome, over the window transcripts.
    let chimeric_alignments = if params.chim_segment_min > 0 {
        let read = crate::chimeric::ChimRead::new(
            read_seq,
            [read_seq.len(), 0],
            read_name,
            [read_seq, &[]],
        );
        crate::chimeric::chimeric_detection(chim_windows, &read, &index.genome, &scorer, params)
    } else {
        Vec::new()
    };

    // n_for_mapq = transcripts.len() after dedup and filtering.
    // Multi-transcript DP (Phase 16.10) produces multiple transcripts per window
    // for tandem repeats (e.g. rDNA), yielding correct NH → correct MAPQ.
    let mut n_for_mapq = transcripts.len();

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
                t.cigar_string()
            );
        }
    }

    // Too short and too many mismatches drop every alignment; too many loci is
    // reached only when the best alignment passed both (STAR's unmapType 1, 2, 3).
    let unmapped_reason = if transcripts.is_empty() {
        // No window kept a transcript: STAR's `nW==0`, counted as "other".
        Some(UnmappedReason::Other)
    } else if let Some(reason) = unmapped_by_filter {
        // Only too-many-loci carries a loci count to the stats; a read that
        // failed the quality gates is not counted as mapped anywhere.
        transcripts.clear();
        n_for_mapq = 0;
        Some(reason)
    } else if mapped_filter && transcripts.len() > params.out_filter_multimap_nmax as usize {
        // `n_for_mapq` already holds the loci count; drop the alignments.
        transcripts.clear();
        Some(UnmappedReason::TooManyLoci)
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
    // The piece offsets are STAR's `splitR[0][ip]`: mate1 starts the concatenated
    // read, mate2 starts one base past the spacer. Only the `flagDirMap`
    // shortcut consults them.
    let mut combined_seeds = Seed::find_seeds_at(
        &combined_read[..len1],
        0,
        index,
        params.seed_map_min,
        params,
        debug_name,
    )?;
    let mut m2_seeds = Seed::find_seeds_at(
        &combined_read[len1 + 1..],
        len1 + 1,
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
        s.mate_id = u8::from(s.read_pos >= len1);
    }

    // Cluster combined seeds using the combined read length
    let clusters = cluster_seeds(&combined_seeds, index, params, combined_len, debug_pe);

    // Combined score threshold: use len1+len2 as denominator
    let combined_score_threshold =
        (params.out_filter_score_min_over_lread * (len1 + len2) as f64) as i32;

    let mut joint_pairs: Vec<PairedAlignment> = Vec::new();
    let mut single_mate1_transcripts: Vec<Transcript> = Vec::new();
    let mut single_mate2_transcripts: Vec<Transcript> = Vec::new();
    // STAR's `trAll` for chimeric detection: per window, the combined-read
    // transcripts — whole pairs and single mates alike.
    let mut chim_windows: Vec<Vec<crate::chimeric::WinTr>> = Vec::new();
    let chim_on = params.chim_segment_min > 0;
    // Whether any window kept a transcript. STAR counts a read with none as
    // "other" (`nW==0`), not "too short"; a pair our pairing step later
    // rejects still counts, since STAR's `trAll` has it.
    let mut any_transcript = false;

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
        );
        let mut chim_window: Vec<crate::chimeric::WinTr> = Vec::new();
        let str_ = u8::from(stitch_is_reverse);

        for wt in &wts {
            let split_result =
                split_combined_wt(wt, len1, len2, stitch_is_reverse, scorer.align_intron_min);
            if let Some((m1_wt, m2_wt)) = split_result {
                let (m1_read_slice, m1_orig_rev, m2_read_slice, m2_orig_rev) = if stitch_is_reverse
                {
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
                    m1_orig_rev,  // no_left_ext = inner for reverse (orig_is_rev=true)
                    !m1_orig_rev, // no_right_ext = inner for forward (orig_is_rev=false)
                    0,            // mate1
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
                    m2_orig_rev,  // no_left_ext = inner for reverse (orig_is_rev=true)
                    !m2_orig_rev, // no_right_ext = inner for forward (orig_is_rev=false)
                    1,            // mate2
                ) else {
                    continue;
                };

                any_transcript = true;
                if stitch_is_reverse {
                    t1.is_reverse = true;
                    t2.is_reverse = false;
                } else {
                    t1.is_reverse = false;
                    t2.is_reverse = true;
                }

                let combined_span =
                    t1.genome_end.max(t2.genome_end) - t1.genome_start.min(t2.genome_start);
                let combined_wt_score = wt.score + scorer.genomic_length_penalty(combined_span);

                let pair = try_pair_transcripts(
                    &t1,
                    &t2,
                    len1,
                    len2,
                    params,
                    combined_score_threshold,
                    combined_wt_score,
                );
                // STAR's `trAll` holds every pair stitching kept. The absolute
                // score floor in `try_pair_transcripts` is an early `mappedFilter`,
                // not a stitching rule, so chimeric detection must still see a
                // pair that only fails that floor — it is often `trBest`.
                if chim_on
                    && (pair.is_some()
                        || try_pair_transcripts(
                            &t1,
                            &t2,
                            len1,
                            len2,
                            params,
                            i32::MIN,
                            combined_wt_score,
                        )
                        .is_some())
                {
                    // The strand frame starts with mate1 forward, mate2
                    // reversed: `[mate1 | spacer | RC(mate2)]` and its reverse
                    // complement.
                    let (first, ff, second, sf, off) = if stitch_is_reverse {
                        (&t2, 1, &t1, 0, len2 + 1)
                    } else {
                        (&t1, 0, &t2, 1, len1 + 1)
                    };
                    chim_window.extend(crate::chimeric::WinTr::pair(
                        first,
                        ff,
                        second,
                        sf,
                        off,
                        str_,
                        combined_wt_score,
                        combined_len,
                    ));
                }
                if let Some(pair) = pair {
                    joint_pairs.push(pair);
                }
            } else {
                // Single-mate WT: save for half-mapped fallback
                let all_m1 = wt.exons.iter().all(|e| e.mate_id == 0);
                let all_m2 = wt.exons.iter().all(|e| e.mate_id == 1);
                if all_m1 {
                    let (read_slice, orig_rev) = if stitch_is_reverse {
                        (&stitch_read[len2 + 1..], true)
                    } else {
                        (&stitch_read[..len1], false)
                    };
                    let m1_wt = crate::align::stitch::rebase_single_mate_wt(
                        wt,
                        crate::align::stitch::mate_read_offset(0, len1, len2, stitch_is_reverse),
                    );
                    if let Some(mut t) = finalize_transcript(
                        &m1_wt,
                        read_slice,
                        index,
                        &scorer,
                        &stitch_cluster,
                        orig_rev,
                        false,
                        false,
                        0, // mate1
                    ) {
                        t.is_reverse = stitch_is_reverse;
                        any_transcript = true;
                        if chim_on {
                            chim_window.extend(crate::chimeric::WinTr::single(
                                &t,
                                str_,
                                crate::align::stitch::mate_read_offset(
                                    0,
                                    len1,
                                    len2,
                                    stitch_is_reverse,
                                ),
                                0,
                                combined_len,
                            ));
                        }
                        single_mate1_transcripts.push(t);
                    }
                } else if all_m2 {
                    let (read_slice, orig_rev) = if stitch_is_reverse {
                        (&stitch_read[..len2], false)
                    } else {
                        (&stitch_read[len1 + 1..], true)
                    };
                    let m2_wt = crate::align::stitch::rebase_single_mate_wt(
                        wt,
                        crate::align::stitch::mate_read_offset(1, len1, len2, stitch_is_reverse),
                    );
                    if let Some(mut t) = finalize_transcript(
                        &m2_wt,
                        read_slice,
                        index,
                        &scorer,
                        &stitch_cluster,
                        orig_rev,
                        false,
                        false,
                        1, // mate2
                    ) {
                        t.is_reverse = !stitch_is_reverse;
                        any_transcript = true;
                        if chim_on {
                            chim_window.extend(crate::chimeric::WinTr::single(
                                &t,
                                str_,
                                crate::align::stitch::mate_read_offset(
                                    1,
                                    len1,
                                    len2,
                                    stitch_is_reverse,
                                ),
                                1,
                                combined_len,
                            ));
                        }
                        single_mate2_transcripts.push(t);
                    }
                }
            }
        }
        if chim_on {
            chim_windows.push(chim_window);
        }
    }

    // STAR's chimericDetection over the combined-read window transcripts. It
    // runs after multMapSelect/mappedFilter whatever their outcome, and reads
    // nothing they decide, so it is done here, before the decision tree
    // consumes the pairs.
    let pe_chimeric = if chim_on {
        let read = crate::chimeric::ChimRead::new(
            &combined_read,
            [len1, len2],
            read_name,
            [mate1_seq, mate2_seq],
        );
        crate::chimeric::chimeric_detection(chim_windows, &read, &index.genome, &scorer, params)
    } else {
        Vec::new()
    };

    // --peOverlapNbasesMin: if the mates overlap in genome space, merge them into one
    // single-end read, align that merged read through the same SE pipeline (`align_read`,
    // this module), and convert each resulting SE transcript back into a two-mate pair,
    // rescoring from scratch in the PE frame (`crate::align::pe_overlap`). STAR's semantics
    // (mirrored here, per `star_pe_overlap.rs` in the sister project STAR-rs) are an
    // unconditional overwrite: if the merge succeeds and at least one transcript converts,
    // the separate-mate `joint_pairs` computed above is replaced outright — no additional
    // score comparison, and the multimap-range filter below is not skipped, it just runs once
    // on the replaced set instead of twice.
    if params.pe_overlap_nbases_min > 0 {
        let m1_slice = &combined_read[..len1];
        let m2rc_slice = &combined_read[len1 + 1..];
        if let Some(merge) = crate::align::pe_overlap::pe_merge_mates(
            m1_slice,
            m2rc_slice,
            params.pe_overlap_nbases_min,
            params.pe_overlap_mmp,
        ) {
            // STAR's peMergeRA->mapOneRead() is a FIND step: it locates windows on the
            // merged read WITHOUT the read-length quality gates (mappedFilter runs later,
            // at the PE level, on the reconstructed pair — ReadAlign_oneRead.cpp:87-91).
            // Applying the SE mappedFilter here (on the longer merged length) would drop
            // valid merges the PE stage would keep. So find on the merged read with the
            // length-relative gates disabled; the PE decision tree below does the real
            // filtering after convert_merged_transcript_to_pe.
            let mut merge_params = params.clone();
            merge_params.out_filter_match_nmin = 0;
            merge_params.out_filter_match_nmin_over_lread = 0.0;
            merge_params.out_filter_score_min_over_lread = 0.0;
            let (merged_transcripts, _merged_chim, _n_mapq, _unmapped) =
                align_read_inner(&merge.merged, read_name, index, &merge_params, false)?;
            let mut converted: Vec<PairedAlignment> = Vec::new();
            for t in &merged_transcripts {
                let Some((m1, m2)) = crate::align::pe_overlap::convert_merged_transcript_to_pe(
                    t,
                    &merge,
                    mate1_seq,
                    mate2_seq,
                    &index.genome,
                    &scorer,
                ) else {
                    continue;
                };
                let m1_span = m1.genome_end - m1.genome_start;
                let m2_span = m2.genome_end - m2.genome_start;
                let combined_span =
                    m1.genome_end.max(m2.genome_end) - m1.genome_start.min(m2.genome_start);
                let combined_wt_score = (m1.score - scorer.genomic_length_penalty(m1_span))
                    + (m2.score - scorer.genomic_length_penalty(m2_span))
                    + scorer.genomic_length_penalty(combined_span);
                let is_proper_pair = check_proper_pair(&m1, &m2, params);
                let insert_size = calculate_insert_size(&m1, &m2);
                converted.push(PairedAlignment {
                    mate1_transcript: m1,
                    mate2_transcript: m2,
                    mate1_region: (0, len1),
                    mate2_region: (0, len2),
                    is_proper_pair,
                    insert_size,
                    combined_wt_score,
                });
            }
            if !converted.is_empty() {
                joint_pairs = converted;
            }
        }
    }

    // --- Decision tree: dedup, score-filter, quality-filter, then half-mapped fallback ---

    // Step 1: position dedup — remove exact (chr, mate1_pos, mate2_pos, strand, CIGAR) duplicates.
    // Run dedup BEFORE score-range filter so the backup pool is already deduplicated.
    // (STAR's ordering is multMapSelect → dedup, but dedup before multMapSelect is equivalent
    // since removing exact duplicates doesn't change the best score.)
    // Order-preserving exact dedup (keep first in discovery order). We must not
    // sort by position here: STAR's primary tie-break falls back to the
    // earliest-discovered pair, so discovery order has to survive to the
    // primary-selection sort below.
    {
        let mut seen = std::collections::HashSet::new();
        joint_pairs.retain(|p| {
            seen.insert((
                p.mate1_transcript.chr_idx,
                p.mate1_transcript.genome_start,
                p.mate1_transcript.is_reverse,
                p.mate1_transcript.cigar_string(),
                p.mate2_transcript.genome_start,
                p.mate2_transcript.is_reverse,
                p.mate2_transcript.cigar_string(),
            ))
        });
    }

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

    // Deterministic primary tie-break (combined score, then a fixed positional
    // order on mate1).
    joint_pairs.sort_by(|a, b| {
        b.combined_wt_score.cmp(&a.combined_wt_score).then_with(|| {
            (
                a.mate1_transcript.chr_idx,
                a.mate1_transcript.genome_start,
                a.mate1_transcript.is_reverse,
            )
                .cmp(&(
                    b.mate1_transcript.chr_idx,
                    b.mate1_transcript.genome_start,
                    b.mate1_transcript.is_reverse,
                ))
        })
    });

    // Primary selection — STAR's multMapSelect. STAR's default does not use the
    // RNG for the primary; only `--outMultimapperOrder Random` shuffles. Gate
    // the (previously unconditional) shuffle so the default is deterministic
    // and STAR-faithful; under Random, shuffle the tied top-score prefix with a
    // per-read seed (deterministic per read → thread-count invariant).
    if params.out_multimapper_order == MultimapperOrder::Random {
        shuffle_tied_prefix(
            &mut joint_pairs,
            |pa| pa.combined_wt_score,
            per_read_seed(params.run_rng_seed, read_name),
        );
    }

    // Step 4: quality filter (mappedFilter). The best pair is STAR's `trBest`
    // whenever a pair exists (a single mate cannot clear the combined-length
    // gates), so a failure here unmaps the read with that reason rather than
    // falling through to the half-mapped fallback.
    if let Some(reason) = filter_paired_transcripts(
        &mut joint_pairs,
        mate1_seq,
        mate2_seq,
        &index.genome,
        params,
    ) {
        return Ok((Vec::new(), pe_chimeric, 0, Some(reason)));
    }

    // Step 5: too-many-loci — STAR checks `multi` AFTER mappedFilter, only when the
    // best pair passes the score/match/mismatch gates (ReadAlign_mappedFilter.cpp).
    // Pairs that survive quality but exceed outFilterMultimapNmax → too many loci;
    // a read whose best pair fails quality was already emptied above (→ unmapped).
    if joint_pairs.len() > params.out_filter_multimap_nmax as usize {
        let n_loci = joint_pairs.len();
        return Ok((
            Vec::new(),
            pe_chimeric,
            n_loci,
            Some(UnmappedReason::TooManyLoci),
        ));
    }

    if !joint_pairs.is_empty() {
        let pe_mapq_n = joint_pairs.len().max(1);
        let results = joint_pairs
            .into_iter()
            .map(|pa| PairedAlignmentResult::BothMapped(Box::new(pa)))
            .collect();
        return Ok((results, pe_chimeric, pe_mapq_n, None));
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
        // STAR splits this case in two (`ReadAlign_mappedFilter.cpp`): a read
        // with no good window at all is `unmappedOther` ("other"), and only a
        // read that *had* a window whose best transcript failed the score or
        // length thresholds is `unmappedShort` ("too short"). Reporting every
        // unmapped pair as too short empties the "other" bucket, which is what
        // #48 measured on paired-end data.
        (None, None) => {
            let reason = if clusters.is_empty() || !any_transcript {
                UnmappedReason::Other
            } else {
                UnmappedReason::TooShort
            };
            Ok((Vec::new(), pe_chimeric, 0, Some(reason)))
        }
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
    })
}

/// Check if paired alignment is a proper pair
pub(crate) fn check_proper_pair(
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
pub(crate) fn calculate_insert_size(mate1_trans: &Transcript, mate2_trans: &Transcript) -> i32 {
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
fn filter_paired_transcripts(
    paired_alns: &mut Vec<PairedAlignment>,
    mate1_seq: &[u8],
    mate2_seq: &[u8],
    genome: &crate::genome::Genome,
    params: &Parameters,
) -> Option<UnmappedReason> {
    // STAR's mappedFilter (ReadAlign_mappedFilter.cpp) applies ALL quality thresholds to trBest
    // (the highest-scoring transcript), NOT to each individual transcript. If trBest passes,
    // all transcripts in the score window are included (they affect NH/MAPQ). If trBest fails,
    // the read is unmapped.
    //
    // Returns why the read is unmapped when the best pair fails, in STAR's
    // order: too short (score or matched bases), then too many mismatches.
    // The pairs arrive sorted by score (and shuffled among ties under
    // `--outMultimapperOrder Random`), so the first is the primary, as in SE.
    // `max_by_key` would return the *last* of tied pairs.
    let best_pa = paired_alns.first();
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
            return Some(UnmappedReason::TooShort);
        }

        // STAR's `nMatch` (read == genome) and `rLength` (aligned bases) of the
        // combined transcript: the mismatch ratio is over `rLength`, not the
        // read length (`ReadAlign_mappedFilter.cpp`).
        let combined_match = star_n_match(&best.mate1_transcript, mate1_seq, genome)
            + star_n_match(&best.mate2_transcript, mate2_seq, genome);
        let combined_r_length =
            star_r_length(&best.mate1_transcript) + star_r_length(&best.mate2_transcript);
        if combined_match < params.out_filter_match_nmin
            || combined_match < (params.out_filter_match_nmin_over_lread * combined_lread_m1) as u32
        {
            paired_alns.clear();
            return Some(UnmappedReason::TooShort);
        }

        if combined_nm > params.out_filter_mismatch_nmax
            || f64::from(combined_nm) / f64::from(combined_r_length.max(1))
                > params.out_filter_mismatch_nover_lmax
        {
            paired_alns.clear();
            return Some(UnmappedReason::TooManyMismatches);
        }
    }
    // Too many loci is the caller's to decide, after these gates and with the
    // loci count intact; clearing here made it unreachable and turned those
    // reads into "too short" (STAR: unmapType 3 vs 1).
    None
}

/// Extract splice junctions from a Transcript's CIGAR as (donor, acceptor) pairs.
/// Junction coords are in genomic space (0-based). Junctions only exist where CigarOp::RefSkip is.
fn extract_junctions_from_cigar(t: &Transcript) -> Vec<(u64, u64)> {
    let mut junctions = Vec::new();
    let mut genome_pos = t.genome_start;
    for op in &t.cigar {
        use noodles::sam::alignment::record::cigar::op::Kind;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                genome_pos += op.len() as u64;
            }
            Kind::Skip => {
                let donor = genome_pos;
                let acceptor = genome_pos + op.len() as u64;
                junctions.push((donor, acceptor));
                genome_pos = acceptor;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
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
    use crate::align::score::SpliceMotif;
    use crate::params::IntronMotifFilter;

    use crate::genome::Genome;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;
    use noodles::sam::alignment::record::cigar;

    fn default_params() -> Parameters {
        // Parse empty args to get default parameters
        Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"])
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
            transform_blocks: None,
            sequence: sequence.into(),
            n_genome,
            n_genome_real: n_genome,
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
            sjdb_overhang: 0,
        }
    }

    #[test]
    fn combined_transcript_for_projection_rewrites_mate2_ifrag() {
        use cigar::op::{Kind, Op};
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
            cigar: vec![Op::new(Kind::Match, (ge - gs) as usize)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };
        let pair = PairedAlignment {
            mate1_transcript: make_tr(1000, 1100, 0, 100),
            mate2_transcript: make_tr(1300, 1400, 0, 100),
            mate1_region: (0, 100),
            mate2_region: (0, 100),
            is_proper_pair: true,
            insert_size: 400,
            combined_wt_score: 200,
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
        let params = default_params();

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
        let mut params = default_params();
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
        let mut params = default_params();
        params.out_filter_mismatch_nmax = 2;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_003", &index, &params);
        assert!(result.is_ok());
    }

    /// A gapless forward or reverse transcript over the test genome.
    fn tx_at(start: u64, len: usize, is_reverse: bool) -> Transcript {
        Transcript {
            chr_idx: 0,
            genome_start: start,
            genome_end: start + len as u64,
            is_reverse,
            exons: vec![Exon {
                genome_start: start,
                genome_end: start + len as u64,
                read_start: 0,
                read_end: len,
                i_frag: 0,
            }],
            cigar: vec![cigar::Op::new(cigar::op::Kind::Match, len)],
            score: len as i32,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
        }
    }

    /// STAR's `nMatch` counts read == genome only; `rLength` is the aligned
    /// length. They differ by the mismatches and by `N` on either side, which
    /// is what makes `--outFilterMatchNminOverLread` stricter than counting
    /// the `M` operations.
    #[test]
    fn star_n_match_counts_only_bases_that_match() {
        let index = make_test_index(); // genome ACGTACGTNN
        let exact = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let t = tx_at(0, 8, false);
        assert_eq!(star_n_match(&t, &exact, &index.genome), 8);
        assert_eq!(star_r_length(&t), 8);

        let mut one_mm = exact.clone();
        one_mm[2] = 0; // G -> A
        assert_eq!(star_n_match(&t, &one_mm, &index.genome), 7);
        assert_eq!(star_r_length(&t), 8, "rLength still counts the mismatch");

        let mut read_n = exact.clone();
        read_n[5] = 4;
        assert_eq!(star_n_match(&t, &read_n, &index.genome), 7, "read N");

        // The genome's trailing NN: aligned, but neither match nor mismatch.
        let over_n = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 0];
        assert_eq!(
            star_n_match(&tx_at(0, 10, false), &over_n, &index.genome),
            8
        );
    }

    /// A reverse transcript's CIGAR walks the reverse complement of the read
    /// as sequenced.
    #[test]
    fn star_n_match_reverse_strand_uses_the_reverse_complement() {
        let index = make_test_index();
        let sequenced = vec![0, 1, 2, 3, 0, 1, 2, 3]; // RC of ACGTACGT is ACGTACGT
        assert_eq!(
            star_n_match(&tx_at(0, 8, true), &sequenced, &index.genome),
            8
        );
        let sequenced = vec![3, 3, 3, 3, 3, 3, 3, 3]; // RC is all A: matches A at 0 and 4
        assert_eq!(
            star_n_match(&tx_at(0, 8, true), &sequenced, &index.genome),
            2
        );
    }

    #[test]
    fn intron_filters_follow_star() {
        use crate::align::score::IntronFilter;
        let check =
            |m: &[SpliceMotif], a: &[bool], f: &IntronFilter| f.passes(m.iter().zip(a.iter()));
        let f = AlignmentScorer::from_params(&default_params()).intron_filter;
        assert!(check(
            &[SpliceMotif::GtAg, SpliceMotif::GcAg],
            &[false, false],
            &f
        ));
        assert!(
            !check(&[SpliceMotif::GtAg, SpliceMotif::CtAc], &[false, false], &f),
            "+ and - junctions in one transcript"
        );
        // Non-canonical junctions carry no strand.
        assert!(check(
            &[SpliceMotif::GtAg, SpliceMotif::NonCanonical],
            &[false, false],
            &f
        ));

        let mut p = default_params();
        p.out_filter_intron_motifs = IntronMotifFilter::RemoveNoncanonicalUnannotated;
        let f = AlignmentScorer::from_params(&p).intron_filter;
        assert!(!check(&[SpliceMotif::NonCanonical], &[false], &f));
        assert!(check(&[SpliceMotif::NonCanonical], &[true], &f));
        p.out_filter_intron_motifs = IntronMotifFilter::RemoveNoncanonical;
        let f = AlignmentScorer::from_params(&p).intron_filter;
        assert!(!check(&[SpliceMotif::NonCanonical], &[true], &f));

        // `--outSAMstrandField intronMotif` also drops a spliced transcript
        // whose strand is undefined (`sjN>0 && sjMotifStrand==0`), but not an
        // unspliced one.
        let mut p = default_params();
        p.out_sam_strand_field = "intronMotif".into();
        let f = AlignmentScorer::from_params(&p).intron_filter;
        assert!(!check(&[SpliceMotif::NonCanonical], &[false], &f));
        assert!(check(
            &[SpliceMotif::NonCanonical, SpliceMotif::GtAg],
            &[false, false],
            &f
        ));
        assert!(check(&[], &[], &f));
    }

    #[test]
    fn test_transcript_multimap_limit() {
        let index = make_test_index();
        let mut params = default_params();
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
        let params = default_params();

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
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let params = default_params();

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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Distance = 300bp, within default limit (auto mode = unlimited)
        assert!(check_proper_pair(&t1, &t2, &params));
    }

    #[test]
    fn test_check_proper_pair_too_far() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let mut params = default_params();
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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Distance = 400bp, exceeds limit of 100bp
        assert!(!check_proper_pair(&t1, &t2, &params));
    }

    #[test]
    fn test_calculate_insert_size_positive() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, 300); // Positive because mate1 is leftmost
    }

    #[test]
    fn test_strand_consistency_filter() {
        use crate::align::transcript::{Exon, Transcript};
        use crate::params::IntronStrandFilter;
        use cigar::op::{Kind, Op};

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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc], // +strand and -strand
            junction_annotated: vec![],
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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::GcAg], // both + strand
            junction_annotated: vec![],
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
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let tlen = calculate_insert_size(&t1, &t2);
        assert_eq!(tlen, -300); // Negative for mate1 in RF pair (mate2 is the left mate)
    }

    #[test]
    fn test_noncanonical_unannotated_filter() {
        use crate::align::score::SpliceMotif;
        use crate::align::transcript::{Exon, Transcript};
        use cigar::op::{Kind, Op};

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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
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
        let params = default_params();

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
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

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
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
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
        top.sort_unstable();
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
            "expected different seeds to pick different primaries, got {firsts:?}"
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
