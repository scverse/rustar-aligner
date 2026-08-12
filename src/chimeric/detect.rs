// Chimeric alignment detection algorithms

use crate::align::SeedCluster;
use crate::align::score::AlignmentScorer;
use crate::align::seed::Seed;
use crate::align::stitch::{cluster_seeds, stitch_seeds, stitch_seeds_with_jdb};
use crate::align::transcript::Transcript;
use crate::chimeric::score::{calculate_repeat_length, classify_junction_type};
use crate::chimeric::segment::{ChimericAlignment, ChimericSegment, ExonSpan};
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::Parameters;

/// Chimeric alignment detector
pub struct ChimericDetector<'a> {
    params: &'a Parameters,
}

impl<'a> ChimericDetector<'a> {
    /// Create a new chimeric detector
    pub fn new(params: &'a Parameters) -> Self {
        Self { params }
    }

    /// Detect chimeric alignments by re-seeding soft-clipped bases (Tier 1 soft-clip re-mapping).
    ///
    /// When the primary alignment has a large soft-clip (>= chimSegmentMin), extract that
    /// clipped sequence and run a new seed search.  If a valid alignment is found it is paired
    /// with the primary transcript to form a chimeric alignment.  Right clips are tried first,
    /// then left clips.  This complements `detect_chimeric_old`, which only searches transcripts
    /// already found during normal seeding.
    pub fn detect_from_soft_clips(
        &self,
        transcript: &Transcript,
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Option<ChimericAlignment>, Error> {
        let params = self.params;
        let min_seg = params.chim_segment_min as usize;
        if min_seg == 0 || transcript.exons.is_empty() {
            return Ok(None);
        }

        let read_len = read_seq.len();
        let [left_clip, right_clip] = transcript.count_soft_clips();
        let score_min = params.chim_score_min;
        let score_drop_max = params.chim_score_drop_max;
        let non_gtag_penalty = params.chim_score_junction_non_gtag;
        let intron_max = params.align_intron_max as u64;
        let overhang_min = params.chim_junction_overhang_min as usize;

        // Try right clip first, then left (match STAR's ordering)
        let candidates = [(right_clip, true), (left_clip, false)];

        for (clip_len, is_right) in candidates {
            if clip_len < min_seg {
                continue;
            }

            let clip_start = if is_right { read_len - clip_len } else { 0 };
            let clip_seq = if is_right {
                &read_seq[clip_start..]
            } else {
                &read_seq[..clip_len]
            };

            // Re-seed the soft-clipped sub-sequence
            let seeds = Seed::find_seeds(clip_seq, index, params.seed_map_min, params, "")?;
            if seeds.is_empty() {
                continue;
            }

            let clusters = cluster_seeds(&seeds, index, params, clip_seq.len(), false);
            if clusters.is_empty() {
                continue;
            }

            let scorer = AlignmentScorer::from_params(params);
            let jdb = if index.junction_db.is_empty() {
                None
            } else {
                Some(&index.junction_db)
            };
            let clip_trs = stitch_seeds_with_jdb(&clusters[0], clip_seq, index, &scorer, jdb, 1);

            let Some(clip_tr_raw) = clip_trs.into_iter().next() else {
                continue;
            };

            if clip_tr_raw.exons.is_empty() {
                continue;
            }

            let clip_aligned =
                clip_tr_raw.exons.last().unwrap().read_end - clip_tr_raw.exons[0].read_start;
            if clip_aligned < min_seg {
                continue;
            }

            // Shift sub-seq read coords into full-read space for right clips
            let clip_tr = if is_right {
                adjust_read_positions(clip_tr_raw, clip_start)
            } else {
                clip_tr_raw
            };

            // Determine donor / acceptor by read order
            let primary_rs = transcript.exons[0].read_start;
            let clip_rs = clip_tr.exons[0].read_start;
            let (tr_donor, tr_acceptor): (&Transcript, &Transcript) = if primary_rs <= clip_rs {
                (transcript, &clip_tr)
            } else {
                (&clip_tr, transcript)
            };

            // Overhang: each segment must cover >= chimJunctionOverhangMin at junction boundary
            let donor_overhang =
                tr_donor.exons.last().unwrap().read_end - tr_donor.exons[0].read_start;
            let acceptor_overhang =
                tr_acceptor.exons.last().unwrap().read_end - tr_acceptor.exons[0].read_start;
            if donor_overhang < overhang_min || acceptor_overhang < overhang_min {
                continue;
            }

            // Classify junction for score adjustment
            let junction_type = classify_junction_type(
                &index.genome,
                tr_donor.chr_idx,
                tr_donor.genome_end,
                tr_donor.is_reverse,
                tr_acceptor.chr_idx,
                tr_acceptor.genome_start,
                tr_acceptor.is_reverse,
            );

            let combined_score = tr_donor.score + tr_acceptor.score;
            let effective_score = if junction_type == 0 {
                combined_score + non_gtag_penalty
            } else {
                combined_score
            };

            if effective_score < score_min {
                continue;
            }
            if effective_score + score_drop_max < read_len as i32 {
                continue;
            }

            // Must be geometrically chimeric (different chr/strand, or span > alignIntronMax)
            let is_chimeric = tr_donor.chr_idx != tr_acceptor.chr_idx
                || tr_donor.is_reverse != tr_acceptor.is_reverse
                || {
                    let span = if tr_donor.genome_end <= tr_acceptor.genome_start {
                        tr_acceptor.genome_start - tr_donor.genome_end
                    } else {
                        tr_donor.genome_start.saturating_sub(tr_acceptor.genome_end)
                    };
                    intron_max > 0 && span > intron_max
                };

            if !is_chimeric {
                continue;
            }

            let donor_seg = transcript_to_segment(tr_donor)
                .map_err(|e| Error::Chimeric(format!("soft-clip donor: {e}")))?;
            let acceptor_seg = transcript_to_segment(tr_acceptor)
                .map_err(|e| Error::Chimeric(format!("soft-clip acceptor: {e}")))?;

            let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
                &index.genome,
                donor_seg.chr_idx,
                donor_seg.genome_end,
                acceptor_seg.chr_idx,
                acceptor_seg.genome_start,
                20,
            );

            let chim = ChimericAlignment::new(
                donor_seg,
                acceptor_seg,
                junction_type,
                repeat_len_donor,
                repeat_len_acceptor,
                read_seq.to_vec(),
                read_name.to_string(),
            );

            return Ok(Some(chim));
        }

        Ok(None)
    }

    /// Re-seed the outer uncovered read regions of an existing chimeric pair (Tier 3).
    ///
    /// After Tiers 1 and 2 find a donor+acceptor chimeric alignment, the read may still
    /// have uncovered bases at the left of the donor or the right of the acceptor.  If
    /// either uncovered span is >= chimSegmentMin, re-seed it and attempt to form an
    /// additional chimeric alignment with the adjacent segment.  This enables detection
    /// of multi-junction chimeric reads (3-way gene fusions).
    pub fn detect_from_chimeric_residuals(
        &self,
        chim: &ChimericAlignment,
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Vec<ChimericAlignment>, Error> {
        let params = self.params;
        let min_seg = params.chim_segment_min as usize;
        if min_seg == 0 {
            return Ok(vec![]);
        }

        let read_len = read_seq.len();
        let score_min = params.chim_score_min;
        let score_drop_max = params.chim_score_drop_max;
        let non_gtag_penalty = params.chim_score_junction_non_gtag;
        let intron_max = params.align_intron_max as u64;
        let overhang_min = params.chim_junction_overhang_min as usize;

        // Outer boundaries of the existing chimeric pair in read space
        let left_covered = chim.donor.read_start.min(chim.acceptor.read_start);
        let right_covered = chim.donor.read_end.max(chim.acceptor.read_end);

        // Which segment is at the left / right boundary
        let left_partner = if chim.donor.read_start <= chim.acceptor.read_start {
            &chim.donor
        } else {
            &chim.acceptor
        };
        let right_partner = if chim.donor.read_end >= chim.acceptor.read_end {
            &chim.donor
        } else {
            &chim.acceptor
        };

        // [clip_start, clip_end) → partner it would be paired with
        let candidates = [
            (0usize, left_covered, left_partner),
            (right_covered, read_len, right_partner),
        ];

        let mut results = Vec::new();

        for (clip_start, clip_end, partner_seg) in candidates {
            let clip_len = clip_end - clip_start;
            if clip_len < min_seg {
                continue;
            }

            let clip_seq = &read_seq[clip_start..clip_end];

            let seeds = Seed::find_seeds(clip_seq, index, params.seed_map_min, params, "")?;
            if seeds.is_empty() {
                continue;
            }

            let clusters = cluster_seeds(&seeds, index, params, clip_seq.len(), false);
            if clusters.is_empty() {
                continue;
            }

            let scorer = AlignmentScorer::from_params(params);
            let jdb = if index.junction_db.is_empty() {
                None
            } else {
                Some(&index.junction_db)
            };
            let clip_trs = stitch_seeds_with_jdb(&clusters[0], clip_seq, index, &scorer, jdb, 1);

            let Some(clip_tr_raw) = clip_trs.into_iter().next() else {
                continue;
            };
            if clip_tr_raw.exons.is_empty() {
                continue;
            }

            let clip_aligned =
                clip_tr_raw.exons.last().unwrap().read_end - clip_tr_raw.exons[0].read_start;
            if clip_aligned < min_seg {
                continue;
            }

            // Shift sub-seq read coords into full-read space
            let clip_tr = if clip_start > 0 {
                adjust_read_positions(clip_tr_raw, clip_start)
            } else {
                clip_tr_raw
            };

            let new_seg = transcript_to_segment(&clip_tr)
                .map_err(|e| Error::Chimeric(format!("tier3 segment: {e}")))?;

            // Donor / acceptor ordered by read position
            let (donor_seg, acceptor_seg): (&ChimericSegment, &ChimericSegment) =
                if new_seg.read_start <= partner_seg.read_start {
                    (&new_seg, partner_seg)
                } else {
                    (partner_seg, &new_seg)
                };

            // Overhang check
            let donor_overhang = donor_seg.read_end - donor_seg.read_start;
            let acceptor_overhang = acceptor_seg.read_end - acceptor_seg.read_start;
            if donor_overhang < overhang_min || acceptor_overhang < overhang_min {
                continue;
            }

            // Junction classification and score
            let junction_type = classify_junction_type(
                &index.genome,
                donor_seg.chr_idx,
                donor_seg.genome_end,
                donor_seg.is_reverse,
                acceptor_seg.chr_idx,
                acceptor_seg.genome_start,
                acceptor_seg.is_reverse,
            );

            let combined_score = donor_seg.score + acceptor_seg.score;
            let effective_score = if junction_type == 0 {
                combined_score + non_gtag_penalty
            } else {
                combined_score
            };

            if effective_score < score_min || effective_score + score_drop_max < read_len as i32 {
                continue;
            }

            // Geometry check: must be genuinely chimeric
            let is_chimeric = donor_seg.chr_idx != acceptor_seg.chr_idx
                || donor_seg.is_reverse != acceptor_seg.is_reverse
                || {
                    let span = if donor_seg.genome_end <= acceptor_seg.genome_start {
                        acceptor_seg.genome_start - donor_seg.genome_end
                    } else {
                        donor_seg
                            .genome_start
                            .saturating_sub(acceptor_seg.genome_end)
                    };
                    intron_max > 0 && span > intron_max
                };

            if !is_chimeric {
                continue;
            }

            let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
                &index.genome,
                donor_seg.chr_idx,
                donor_seg.genome_end,
                acceptor_seg.chr_idx,
                acceptor_seg.genome_start,
                20,
            );

            results.push(ChimericAlignment::new(
                donor_seg.clone(),
                acceptor_seg.clone(),
                junction_type,
                repeat_len_donor,
                repeat_len_acceptor,
                read_seq.to_vec(),
                read_name.to_string(),
            ));
        }

        Ok(results)
    }

    /// Detect chimeric alignments from multi-cluster seeds (Tier 2)
    ///
    /// Triggers when:
    /// - Seeds cluster on different chromosomes
    /// - Seeds cluster on different strands (same chromosome)
    /// - Seeds cluster with large genomic distance (>1Mb, same chr/strand)
    pub fn detect_from_multi_clusters(
        &self,
        clusters: &[SeedCluster],
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Vec<ChimericAlignment>, Error> {
        let mut chimeras = Vec::new();

        // Find cluster pairs with chimeric signatures
        for i in 0..clusters.len() {
            for j in (i + 1)..clusters.len() {
                if is_chimeric_signature(&clusters[i], &clusters[j]) {
                    // Try to build chimeric alignment from these clusters
                    if let Some(chim) = self.build_chimeric_from_clusters(
                        &clusters[i],
                        &clusters[j],
                        read_seq,
                        read_name,
                        index,
                    )? {
                        chimeras.push(chim);
                    }
                }
            }
        }

        Ok(chimeras)
    }

    /// Build chimeric alignment from two clusters
    fn build_chimeric_from_clusters(
        &self,
        cluster1: &SeedCluster,
        cluster2: &SeedCluster,
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Option<ChimericAlignment>, Error> {
        if cluster1.alignments.is_empty() || cluster2.alignments.is_empty() {
            return Ok(None);
        }

        // Stitch each cluster independently using existing stitch_seeds
        use crate::align::score::AlignmentScorer;
        let scorer = AlignmentScorer::from_params(self.params);

        let transcripts1 = stitch_seeds(cluster1, read_seq, index, &scorer);
        let transcripts2 = stitch_seeds(cluster2, read_seq, index, &scorer);

        if transcripts1.is_empty() || transcripts2.is_empty() {
            return Ok(None);
        }

        // Take best transcript from each cluster
        let t1 = &transcripts1[0];
        let t2 = &transcripts2[0];

        // Check that both transcripts have exons
        if t1.exons.is_empty() || t2.exons.is_empty() {
            return Ok(None);
        }

        // Determine donor/acceptor based on read position
        let (donor_t, acceptor_t) = if t1.exons[0].read_start < t2.exons[0].read_start {
            (t1, t2)
        } else {
            (t2, t1)
        };

        // Convert transcripts to chimeric segments
        let donor = transcript_to_segment(donor_t)?;
        let acceptor = transcript_to_segment(acceptor_t)?;

        // Check minimum segment lengths
        if !donor.meets_min_length(self.params.chim_segment_min)
            || !acceptor.meets_min_length(self.params.chim_segment_min)
        {
            return Ok(None);
        }

        // Classify junction type
        let junction_type = classify_junction_type(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            donor.is_reverse,
            acceptor.chr_idx,
            acceptor.genome_start,
            acceptor.is_reverse,
        );

        // Calculate repeat lengths
        let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            acceptor.chr_idx,
            acceptor.genome_start,
            20, // max check distance
        );

        // Create chimeric alignment
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            junction_type,
            repeat_len_donor,
            repeat_len_acceptor,
            read_seq.to_vec(),
            read_name.to_string(),
        );

        Ok(Some(chim))
    }
}

/// Calculate genomic distance between two clusters
fn genomic_distance(c1: &SeedCluster, c2: &SeedCluster) -> u64 {
    if c1.chr_idx != c2.chr_idx {
        return u64::MAX;
    }

    if c1.genome_end < c2.genome_start {
        c2.genome_start - c1.genome_end
    } else {
        c1.genome_start.saturating_sub(c2.genome_end)
    }
}

/// Check if two clusters represent a chimeric signature
fn is_chimeric_signature(c1: &SeedCluster, c2: &SeedCluster) -> bool {
    // Different chromosomes
    if c1.chr_idx != c2.chr_idx {
        return true;
    }

    // Different strands (same chromosome)
    if c1.is_reverse != c2.is_reverse {
        return true;
    }

    // Large genomic distance (same chr/strand)
    let distance = genomic_distance(c1, c2);
    if distance > 1_000_000 {
        return true;
    }

    false
}

/// Detect inter-mate chimeric alignment from two single-mate transcripts.
///
/// Fires when mate1 and mate2 map to different chromosomes, opposite-orientation
/// same-chromosome positions, or positions too far apart to be a normal PE pair.
/// This is the primary PE-specific chimeric case (gene-fusion detection).
pub fn detect_inter_mate_chimeric(
    t1: &Transcript,
    t2: &Transcript,
    mate1_seq: &[u8],
    read_name: &str,
    params: &Parameters,
    index: &GenomeIndex,
) -> Option<ChimericAlignment> {
    // Only fire if the pair is discordant (different chr, same strand (both FW or both RC =
    // not FR orientation), or too far apart).
    let is_inter_chr = t1.chr_idx != t2.chr_idx;
    // FR pair expects t1.is_reverse=false (mate1 FW) and t2.is_reverse=true (mate2 RC).
    // Chimeric if both same strand.
    let same_strand = t1.is_reverse == t2.is_reverse;
    let too_far = if t1.chr_idx == t2.chr_idx {
        let left_end = t1.genome_end.min(t2.genome_end);
        let right_start = t1.genome_start.max(t2.genome_start);
        right_start > left_end && right_start - left_end > 1_000_000
    } else {
        false
    };

    if !is_inter_chr && !same_strand && !too_far {
        return None;
    }

    if t1.exons.is_empty() || t2.exons.is_empty() {
        return None;
    }

    // Convert transcripts to chimeric segments
    let donor = transcript_to_segment(t1).ok()?;
    let acceptor = transcript_to_segment(t2).ok()?;

    if !donor.meets_min_length(params.chim_segment_min)
        || !acceptor.meets_min_length(params.chim_segment_min)
    {
        return None;
    }

    // Junction type: non-canonical (0) for inter-chromosomal; try motif for same-chr
    let junction_type = if is_inter_chr {
        0
    } else {
        classify_junction_type(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            donor.is_reverse,
            acceptor.chr_idx,
            acceptor.genome_start,
            acceptor.is_reverse,
        )
    };

    let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
        &index.genome,
        donor.chr_idx,
        donor.genome_end,
        acceptor.chr_idx,
        acceptor.genome_start,
        20,
    );

    let chim = ChimericAlignment::new(
        donor,
        acceptor,
        junction_type,
        repeat_len_donor,
        repeat_len_acceptor,
        mate1_seq.to_vec(),
        read_name.to_string(),
    );

    // Same `--chimFilter` treatment as the intra-mate path: no detection route
    // gets to skip it.
    // No `N` check here, deliberately. STAR runs the junction scan only in the
    // "chimeric junction is within one of the mates" branch
    // (`ReadAlign_chimericDetectionOld.cpp:144` onward). The branch above it,
    // for mates that *bracket* the junction
    // (`exons[e0][EX_iFrag] < exons[e1][EX_iFrag]`, `:130-143`), sets
    // `chimMotif=-1` and returns without ever reading `b0`, `b1` or `bR`;
    // `ChimericAlign_chimericStitching.cpp:22` short-circuits the same way.
    // There is no junction position within a mate to scan, so filtering here
    // would drop pairs STAR reports — and would run the scan over a
    // fragment-wide geometry it was never derived for.
    let _ = index;
    Some(chim)
}

/// Shift all exon read_start/read_end values in a transcript by `offset`.
///
/// Used when a transcript was stitched against a sub-slice of the read (e.g. a right soft-clip
/// at position `offset`) so that its read coordinates become relative to the full read.
fn adjust_read_positions(mut tr: Transcript, offset: usize) -> Transcript {
    for exon in &mut tr.exons {
        exon.read_start += offset;
        exon.read_end += offset;
    }
    tr
}

/// Compute read-orientation (5'→3' of original read) start/end for a transcript.
///
/// STAR uses "ro" coords so that clipping amounts are always measured from the 5' end of the
/// original read regardless of mapping strand.  For the SAM CIGAR convention used internally:
/// - Forward: ro_start = exons.first().read_start, ro_end = exons.last().read_end − 1
/// - Reverse:  ro_start = Lread − exons.last().read_end, ro_end = Lread − exons.first().read_start − 1
fn ro_coords(transcript: &Transcript, read_len: usize) -> (usize, usize) {
    if transcript.exons.is_empty() {
        return (0, 0);
    }
    let first = transcript.exons.first().unwrap();
    let last = transcript.exons.last().unwrap();
    if !transcript.is_reverse {
        (first.read_start, last.read_end.saturating_sub(1))
    } else {
        (
            read_len.saturating_sub(last.read_end),
            read_len.saturating_sub(first.read_start + 1),
        )
    }
}

/// Overlap, in read-orientation coordinates, of two segments' covered spans.
fn ro_overlap(a: (usize, usize), b: (usize, usize)) -> usize {
    let ((a_start, a_end), (b_start, b_end)) = (a, b);
    if b_start > a_start {
        if b_start > a_end {
            0
        } else {
            a_end - b_start + 1
        }
    } else if b_end < a_start {
        0
    } else {
        b_end - a_start + 1
    }
}

/// STAR's chimeric strand code for one segment: 0 undefined, 1 same as the RNA,
/// 2 opposite (`ChimericDetection::chimericDetectionMult`, `chimStr`).
///
/// A transcript with no annotated junction motif is undefined and pairs with
/// either strand. STAR derives the code from the presence of a `+` motif alone,
/// so a transcript mixing `+` and `-` motifs is not rejected here — it takes the
/// `+` answer. That is STAR's behaviour and this reproduces it.
fn motif_strand(tr: &Transcript) -> u8 {
    use crate::align::score::SpliceMotif;
    let has_plus = tr
        .junction_motifs
        .iter()
        .any(|m| matches!(m, SpliceMotif::GtAg | SpliceMotif::GcAg | SpliceMotif::AtAc));
    let has_minus = tr
        .junction_motifs
        .iter()
        .any(|m| matches!(m, SpliceMotif::CtAc | SpliceMotif::CtGc | SpliceMotif::GtAt));
    if !has_plus && !has_minus {
        0
    } else if tr.is_reverse != has_plus {
        1
    } else {
        2
    }
}

/// STAR's `ChimericDetection::chimericDetectionMult` (`--chimMultimapNmax > 0`):
/// enumerate *every* chimeric alignment of a read, not just the best one.
///
/// The old path ([`detect_chimeric_old`]) pins the best transcript as one segment
/// and looks for a partner, so it can only ever emit one chimera. This one runs a
/// triangular loop over all transcript pairs, keeps each pair that clears the
/// score floor, and reports those within `--chimMultimapScoreRange` of the best.
///
/// Two STAR behaviours the loop depends on:
///
/// * The floor **ratchets**. It starts at `--chimScoreMin`, is raised above the
///   best linear alignment score, and again to `readLength - chimScoreDropMax`;
///   then every new best chimera raises it to `best - multimapScoreRange`. Pairs
///   found early can therefore be admitted and later dropped, which is why the
///   final `retain` is not redundant with the in-loop test.
/// * More than `--chimMultimapNmax` survivors means *no* output, not a truncated
///   list. The read is chimerically multimapping beyond the cap, and STAR reports
///   nothing rather than an arbitrary subset.
///
/// STAR loops over (window, alignment) pairs and skips `iA2 = iA1 + 1` within a
/// window so each unordered pair is visited once. `all_transcripts` here is the
/// flat, window-ordered pool, so the same de-duplication is the plain `j > i`
/// triangular loop.
#[allow(clippy::too_many_arguments)]
pub fn detect_chimeric_mult(
    all_transcripts: &[Transcript],
    max_nonchim_score: i32,
    read_seq: &[u8],
    read_name: &str,
    params: &Parameters,
    index: &GenomeIndex,
    mate_boundary: Option<usize>,
) -> Result<Vec<ChimericAlignment>, Error> {
    use crate::align::score::SpliceMotif;

    if params.chim_segment_min == 0 || params.chim_multimap_nmax == 0 {
        return Ok(Vec::new());
    }
    let read_len = read_seq.len();
    let min_seg = params.chim_segment_min as usize;
    let gap_max = params.chim_segment_read_gap_max as usize;
    let ban_genomic_n = params.chim_filter.iter().any(|f| f == "banGenomicN");

    // STAR only looks for a chimera when the best linear alignment leaves enough
    // of the read unexplained (`--chimNonchimScoreDropMin`). A read that already
    // aligns end to end is not a fusion candidate.
    if max_nonchim_score > read_len as i32 - params.chim_nonchim_score_drop_min {
        return Ok(Vec::new());
    }

    // A usable segment is long enough and free of non-canonical junctions.
    let seg_ok = |tr: &Transcript| {
        !tr.exons.is_empty()
            && !tr.junction_motifs.contains(&SpliceMotif::NonCanonical)
            && ro_coords(tr, read_len).1 + 1 - ro_coords(tr, read_len).0 >= min_seg
    };

    let mut min_score = params.chim_score_min;
    if max_nonchim_score >= min_score {
        min_score = max_nonchim_score + 1;
    }
    if read_len as i32 - params.chim_score_drop_max > min_score {
        min_score = read_len as i32 - params.chim_score_drop_max;
    }

    let mut chims: Vec<(ChimericAlignment, i32)> = Vec::new();
    let mut chim_score_best = 0i32;

    for (i, tr1) in all_transcripts.iter().enumerate() {
        if !seg_ok(tr1) {
            continue;
        }
        let str1 = motif_strand(tr1);
        let ro1 = ro_coords(tr1, read_len);
        for tr2 in all_transcripts.iter().skip(i + 1) {
            if !seg_ok(tr2) {
                continue;
            }
            let str2 = motif_strand(tr2);
            if str1 != 0 && str2 != 0 && str1 != str2 {
                continue; // chimeric segments must agree on strand
            }
            let ro2 = ro_coords(tr2, read_len);
            let overlap = ro_overlap(ro1, ro2);
            let (len1, len2) = (ro1.1 + 1 - ro1.0, ro2.1 + 1 - ro2.0);
            // STAR writes this as `roE > segmentMin + roS + overlap`, which on
            // inclusive coordinates means a segment must be *longer* than
            // `segmentMin + overlap + 1`, not merely longer than
            // `segmentMin + overlap`.
            if len1 <= min_seg + overlap + 1 || len2 <= min_seg + overlap + 1 {
                continue;
            }
            // Same waiver as the old path: segments in different mates are
            // expected to be far apart in read space.
            let diff_mates = mate_boundary
                .is_some_and(|b| (ro1.1 < b && ro2.0 >= b) || (ro2.1 < b && ro1.0 >= b));
            let gap_ok =
                diff_mates || ((ro1.1 + gap_max + 1 >= ro2.0) && (ro2.1 + gap_max + 1 >= ro1.0));
            if !gap_ok {
                continue;
            }

            let chim_score = tr1.score + tr2.score - overlap as i32;
            if chim_score < min_score {
                continue;
            }

            let Some((chim, score)) = finalize_chimera(
                tr1, tr2, overlap, chim_score, read_len, read_seq, read_name, params, index,
            )?
            else {
                continue;
            };
            // STAR bans an `N`-bearing candidate *inside* stitching, by zeroing
            // its chimScore (`ChimericAlign_chimericStitching.cpp:63-66`). The
            // caller then only pushes candidates that still clear
            // `minScoreToConsider`, so a banned one never enters the list, never
            // raises `chimScoreBest`, never narrows the score floor, and is not
            // counted against `--chimMultimapNmax`
            // (`ChimericDetection_chimericDetectionMult.cpp:79-95`). Filtering
            // after selection instead would let a banned locus evict legitimate
            // ones, or push the survivor count over the cap so nothing is
            // reported at all.
            if !junction_scan_is_clean(&chim, index, ban_genomic_n) {
                continue;
            }
            if score < min_score {
                continue;
            }
            if score > chim_score_best {
                chim_score_best = score;
                if chim_score_best - params.chim_multimap_score_range > min_score {
                    min_score = chim_score_best - params.chim_multimap_score_range;
                }
            }
            chims.push((chim, score));
        }
    }

    if chim_score_best == 0 {
        return Ok(Vec::new());
    }
    chims.retain(|(_, score)| *score >= min_score);
    if chims.len() > params.chim_multimap_nmax {
        return Ok(Vec::new()); // too many chimeric loci: STAR reports none
    }
    // No filter pass here: candidates were banned as they were scored, above.
    Ok(chims.into_iter().map(|(c, _)| c).collect())
}

/// Implement STAR's `chimericDetectionOld()`: find the best chimeric pair from all
/// post-stitching transcripts.
///
/// Algorithm: use the best transcript (`tr_best`) as the primary segment and search
/// every other transcript for a complementary segment that covers a different part of
/// the read at a different genomic location.  Applies STAR's score-drop, uniqueness,
/// segment-length, and read-gap filters, then emits at most one `ChimericAlignment`.
///
/// Called after stitching + dedup for SE reads, and per-mate after split+finalize for PE.
pub fn detect_chimeric_old(
    all_transcripts: &[Transcript],
    tr_best: &Transcript,
    read_seq: &[u8],
    read_name: &str,
    params: &Parameters,
    index: &GenomeIndex,
) -> Result<Vec<ChimericAlignment>, Error> {
    // SE / per-mate PE pool: no combined-read mate boundary, so diffMates never applies.
    let chims = detect_chimeric_old_impl(
        all_transcripts,
        tr_best,
        read_seq,
        read_name,
        params,
        index,
        None,
    )?;
    Ok(apply_chim_filter(chims, params, index))
}

/// Apply `--chimFilter` to detected chimeras.
///
/// STAR treats this as a post-detection filter rather than a condition inside
/// the search, and so does this: the detection paths all funnel through here,
/// so a filter cannot be forgotten on one of them.
///
/// `banGenomicN` (STAR's default) drops a junction whose flanking genomic bases
/// are not real bases. Sequence around an assembly gap produces junctions that
/// look clean by score and are meaningless. `None` keeps everything.
pub fn apply_chim_filter(
    chims: Vec<ChimericAlignment>,
    params: &Parameters,
    index: &GenomeIndex,
) -> Vec<ChimericAlignment> {
    let ban_genomic_n = params.chim_filter.iter().any(|f| f == "banGenomicN");
    chims
        .into_iter()
        .filter(|c| junction_scan_is_clean(c, index, ban_genomic_n))
        .collect()
}

/// STAR's `N` check across the chimeric junction scan.
///
/// This is not a test of the two bases beside the breakpoint. STAR walks the
/// candidate junction positions — the same scan that picks the canonical motif —
/// and abandons the whole chimera on the first `N` it meets
/// (`ReadAlign_chimericDetectionOld.cpp:145-180`,
/// `ChimericAlign_chimericStitching.cpp`):
///
/// ```cpp
/// uint jRmax = roStart1+trChim[1].exons[e1][EX_L];
/// jRmax = jRmax>roStart0 ? jRmax-roStart0-1 : 0;
/// for (jR=0; jR<jRmax; jR++) {
///     char bR=Read1[0][roStart0+jR];
///     ...
///     if ( ( P.pCh.filter.genomicN && (b0>3 || b1>3) ) || bR>3 ) { chimN=0; break; };
/// }
/// ```
///
/// Two details that a two-base check cannot express:
///
/// - **Strand.** For a reverse segment STAR walks the offset from the other end
///   and complements the base (`b0 = G[EX_G + EX_L - 1 - jR]; if (b0<4) b0=3-b0`).
///   Complementing does not change whether the base is `N`, but *which* base is
///   read does, so a forward-only check bans the wrong chimeras.
/// - **The read `N` ban is unconditional.** `bR>3` sits outside the
///   `filter.genomicN` guard, so a read carrying `N` across the junction is
///   rejected even under `--chimFilter None`.
fn junction_scan_is_clean(c: &ChimericAlignment, index: &GenomeIndex, ban_genomic_n: bool) -> bool {
    // STAR works in ORIGINAL-read coordinates, converting a reverse segment's
    // aligned offset back: `roStart = Str==0 ? EX_R : Lread - EX_R - EX_L`.
    // Our `read_start`/`read_end` are the aligned-orientation offsets, so the
    // same conversion is needed or the scan walks the wrong span for a reverse
    // segment -- and the span is what decides how far into a gap it reaches.
    let read_len = c.read_seq.len();
    // STAR takes roStart from the scanned EXON, not the whole segment:
    // `roStart = Str==0 ? exons[e][EX_R] : Lread - exons[e][EX_R] - exons[e][EX_L]`.
    let ro_of = |ex: &ExonSpan, is_reverse: bool| -> usize {
        if is_reverse {
            read_len.saturating_sub(ex.read_end)
        } else {
            ex.read_start
        }
    };
    // Ordering uses the SEGMENT's roStart, as STAR does: it compares
    // `trChim[].roStart`, a transcript-level field, and only picks `e0`/`e1`
    // afterwards. Using an exon-derived offset here would be circular, since
    // which exon applies depends on the ordering.
    let ro_seg = |seg: &ChimericSegment| -> usize {
        if seg.is_reverse {
            read_len.saturating_sub(seg.read_end)
        } else {
            seg.read_start
        }
    };
    // STAR's scan assumes segment 0 is the one that leads in the read: `jRmax`
    // is derived as `roStart1 + L1 - roStart0 - 1`, which clamps to zero if the
    // two are the other way round, silently skipping the scan. Our donor and
    // acceptor are ordered by the junction, not by read position, so order them
    // here rather than trusting the field names.
    let (seg0, seg1) = if ro_seg(&c.donor) <= ro_seg(&c.acceptor) {
        (&c.donor, &c.acceptor)
    } else {
        (&c.acceptor, &c.donor)
    };
    // Once ordered, each segment's scanned exon follows from its role: the
    // leading segment contributes `e0`, the trailing one `e1`.
    let ex0 = seg0.junction_exon(true);
    let ex1 = seg1.junction_exon(false);
    let ro_start0 = ro_of(&ex0, seg0.is_reverse);
    let ro_start1 = ro_of(&ex1, seg1.is_reverse);
    // jRmax = roStart1 + L1 - roStart0 - 1, clamped at 0 as STAR clamps it.
    // jRmax uses the exon's length (`EX_L`), which for a gapless exon is both
    // its read and its reference length.
    let acceptor_read_end = ro_start1 + (ex1.read_end - ex1.read_start);
    let jr_max = acceptor_read_end
        .saturating_sub(ro_start0)
        .saturating_sub(1);

    // An out-of-range position is not an `N`: STAR indexes the genome directly
    // and only ever compares the value, so a position we cannot answer for must
    // not be treated as a ban.
    let base = |pos: u64| index.genome.get_base(pos);

    for jr in 0..jr_max {
        // Read base. STAR reads `Read1[0][roStart0+jR]`; ours is the same read
        // in the same orientation.
        if let Some(&b_r) = c.read_seq.get(ro_start0 + jr)
            && b_r > 3
        {
            return false; // unconditional, regardless of --chimFilter
        }

        if !ban_genomic_n {
            continue;
        }

        let jr64 = jr as u64;
        let b0 = if seg0.is_reverse {
            base(ex0.genome_end.wrapping_sub(1).wrapping_sub(jr64))
        } else {
            base(ex0.genome_start.wrapping_add(jr64))
        };
        // The acceptor is walked in the donor's read frame, so its genome
        // offset carries the (roStart0 - roStart1) shift STAR applies.
        let shift = ro_start0 as i64 + jr as i64 - ro_start1 as i64;
        let b1 = if seg1.is_reverse {
            let pos = ex1.genome_end as i64 - 1 - shift;
            if pos < 0 { None } else { base(pos as u64) }
        } else {
            let pos = ex1.genome_start as i64 + shift;
            if pos < 0 { None } else { base(pos as u64) }
        };

        if b0.is_some_and(|b| b > 3) || b1.is_some_and(|b| b > 3) {
            return false;
        }
    }
    true
}

/// `detect_chimeric_old` with an optional combined-read mate boundary (`read_length[0]`
/// in ro-space). When a candidate segment lies in a different mate than the primary
/// segment (STAR's `diffMates`), the read-gap check is waived — matching
/// `ReadAlign_chimericDetectionOld.cpp` lines 67-71.
#[allow(clippy::too_many_arguments)]
pub fn detect_chimeric_old_impl(
    all_transcripts: &[Transcript],
    tr_best: &Transcript,
    read_seq: &[u8],
    read_name: &str,
    params: &Parameters,
    index: &GenomeIndex,
    mate_boundary: Option<usize>,
) -> Result<Vec<ChimericAlignment>, Error> {
    let read_len = read_seq.len();
    let min_seg = params.chim_segment_min as usize;
    let score_min = params.chim_score_min;
    let score_drop_max = params.chim_score_drop_max;
    let score_separation = params.chim_score_separation;
    let gap_max = params.chim_segment_read_gap_max as usize;
    let main_mult_max = params.chim_main_segment_mult_nmax as usize;

    // STAR: reject if main segment is too multimapping (nTr > mainSegmentMultNmax && nTr!=2)
    let n_total = all_transcripts.len();
    if n_total > main_mult_max && n_total != 2 {
        return Ok(vec![]);
    }

    // ro coords for the best transcript
    if tr_best.exons.is_empty() {
        return Ok(vec![]);
    }
    let (ro_start1, ro_end1) = ro_coords(tr_best, read_len);
    let r_length1 = ro_end1 + 1 - ro_start1; // aligned read bases in primary

    // Main segment must be long enough
    if r_length1 < min_seg {
        return Ok(vec![]);
    }

    // There must be space for a partner segment at one end of the read
    let has_right_space = ro_end1 + min_seg < read_len;
    let has_left_space = ro_start1 >= min_seg;
    if !has_right_space && !has_left_space {
        return Ok(vec![]);
    }

    // Main segment must have no non-canonical junctions and a consistent motif strand
    use crate::align::score::SpliceMotif;
    if tr_best.junction_motifs.contains(&SpliceMotif::NonCanonical) {
        return Ok(vec![]);
    }
    let has_plus = tr_best
        .junction_motifs
        .iter()
        .any(|m| matches!(m, SpliceMotif::GtAg | SpliceMotif::GcAg | SpliceMotif::AtAc));
    let has_minus = tr_best
        .junction_motifs
        .iter()
        .any(|m| matches!(m, SpliceMotif::CtAc | SpliceMotif::CtGc | SpliceMotif::GtAt));
    if has_plus && has_minus {
        return Ok(vec![]);
    }
    // 0=undefined, 1=same as RNA (+ strand), 2=opposite to RNA (- strand)
    let chim_str1: u8 = if !has_plus && !has_minus {
        0
    } else if tr_best.is_reverse != has_plus {
        1
    } else {
        2
    };

    let score1 = tr_best.score;

    let mut chim_score_best: i32 = i32::MIN;
    let mut chim_score_next: i32 = i32::MIN;
    let mut best_tr2: Option<&Transcript> = None;
    let mut best_overlap: usize = 0;

    for tr2 in all_transcripts {
        if std::ptr::eq(tr2, tr_best) {
            continue;
        }
        if tr2.exons.is_empty() {
            continue;
        }

        // Partner must not have non-canonical junctions
        if tr2.junction_motifs.contains(&SpliceMotif::NonCanonical) {
            continue;
        }

        // Partner motif strand
        let has_plus2 = tr2
            .junction_motifs
            .iter()
            .any(|m| matches!(m, SpliceMotif::GtAg | SpliceMotif::GcAg | SpliceMotif::AtAc));
        let has_minus2 = tr2
            .junction_motifs
            .iter()
            .any(|m| matches!(m, SpliceMotif::CtAc | SpliceMotif::CtGc | SpliceMotif::GtAt));
        let chim_str2: u8 = if !has_plus2 && !has_minus2 {
            0
        } else if tr2.is_reverse != has_plus2 {
            1
        } else {
            2
        };

        // Strands must be consistent (STAR: if both defined they must match)
        if chim_str1 != 0 && chim_str2 != 0 && chim_str1 != chim_str2 {
            continue;
        }

        let (ro_start2, ro_end2) = ro_coords(tr2, read_len);

        // Overlap in read orientation coordinates
        let overlap = if ro_start2 > ro_start1 {
            if ro_start2 > ro_end1 {
                0
            } else {
                ro_end1 - ro_start2 + 1
            }
        } else if ro_end2 < ro_start1 {
            0
        } else {
            ro_end2 - ro_start1 + 1
        };

        let r_length2 = ro_end2 + 1 - ro_start2;

        // Both segments must be long enough (after subtracting overlap)
        // Same inclusive-coordinate boundary as the multimap path above
        // (`ReadAlign_chimericDetectionOld.cpp`, `chimericAlignScore`).
        if r_length1 <= min_seg + overlap + 1 || r_length2 <= min_seg + overlap + 1 {
            continue;
        }

        // Read gap check: the two segments must be close enough in read space —
        // UNLESS they come from different mates (STAR's `diffMates`), in which case
        // the inter-mate fragment gap is expected and the check is waived
        // (ReadAlign_chimericDetectionOld.cpp:67-71). `diffMates` requires a combined
        // read with a known mate boundary; it is always false for SE / per-mate pools.
        let diff_mates = mate_boundary
            .is_some_and(|b| (ro_end1 < b && ro_start2 >= b) || (ro_end2 < b && ro_start1 >= b));
        let gap_ok = diff_mates
            || ((ro_end1 + gap_max + 1 >= ro_start2) && (ro_end2 + gap_max + 1 >= ro_start1));
        if !gap_ok {
            continue;
        }

        let score2 = tr2.score;
        let chim_score = score1 + score2 - overlap as i32;

        // Track overlap of partner vs best partner (same-window case)
        let overlap_with_best: usize = if chim_score_best > i32::MIN {
            if let Some(prev) = best_tr2 {
                let (prev_s, prev_e) = ro_coords(prev, read_len);
                if ro_start2 > prev_s {
                    if ro_start2 > prev_e {
                        0
                    } else {
                        prev_e - ro_start2 + 1
                    }
                } else if ro_end2 < prev_s {
                    0
                } else {
                    ro_end2 - prev_s + 1
                }
            } else {
                0
            }
        } else {
            0
        };

        if chim_score > chim_score_best {
            best_tr2 = Some(tr2);
            if overlap_with_best == 0 {
                chim_score_next = chim_score_best;
            }
            chim_score_best = chim_score;
            best_overlap = overlap;
            let _ = chim_str2; // strand info tracked for extension later
        } else if chim_score > chim_score_next && overlap_with_best == 0 {
            chim_score_next = chim_score;
        }
    }

    // No chimeric partner found
    let Some(tr2) = best_tr2 else {
        return Ok(vec![]);
    };

    // Score filters
    if chim_score_best < score_min {
        return Ok(vec![]);
    }
    // Score-drop gate (STAR chimericDetectionOld.cpp:99, `readLength[0]+readLength[1]`).
    // `read_len` is the length of the read passed in, so this scales correctly for both
    // modes: per-mate length for the per-mate PE pools (an intra-mate chimera spans one
    // mate), and the combined length when a combined read is supplied via `mate_boundary`.
    if chim_score_best + score_drop_max < read_len as i32 {
        return Ok(vec![]);
    }
    // Uniqueness: next-best must be clearly worse
    if chim_score_next + score_separation >= chim_score_best {
        return Ok(vec![]);
    }

    let finalized = finalize_chimera(
        tr_best,
        tr2,
        best_overlap,
        chim_score_best,
        read_len,
        read_seq,
        read_name,
        params,
        index,
    )?;

    Ok(finalized
        .map(|(chim, _score)| vec![chim])
        .unwrap_or_default())
}

/// Turn an accepted segment pair into a `ChimericAlignment`, or reject it.
///
/// Everything from here on is common to both detection paths: order the two
/// segments by read position, check the junction overhang and the geometry,
/// classify the junction motif, apply the non-GTAG penalty and re-check the
/// score. `detect_chimeric_old_impl` reaches it once, with its single best
/// pair; [`detect_chimeric_mult`] reaches it for every surviving pair.
///
/// Returns the alignment together with its post-penalty score. The multimap
/// path needs that score, not `ChimericAlignment::total_score`, because the
/// latter is the plain sum of the two segment scores and knows nothing about
/// the read overlap or the motif penalty.
#[allow(clippy::too_many_arguments)]
fn finalize_chimera(
    tr1: &Transcript,
    tr2: &Transcript,
    overlap: usize,
    chim_score: i32,
    read_len: usize,
    read_seq: &[u8],
    read_name: &str,
    params: &Parameters,
    index: &GenomeIndex,
) -> Result<Option<(ChimericAlignment, i32)>, Error> {
    let score_min = params.chim_score_min;
    let score_drop_max = params.chim_score_drop_max;
    let overhang_min = params.chim_junction_overhang_min as usize;
    let non_gtag_penalty = params.chim_score_junction_non_gtag;

    // Determine donor / acceptor by read position
    let (ro_start1, ro_end1) = ro_coords(tr1, read_len);
    let (ro_start2, ro_end2) = ro_coords(tr2, read_len);
    let (tr_donor, tr_acceptor) = if ro_start1 <= ro_start2 {
        (tr1, tr2)
    } else {
        (tr2, tr1)
    };
    let (ro_donor_end, ro_acceptor_start) = if ro_start1 <= ro_start2 {
        (ro_end1, ro_start2)
    } else {
        (ro_end2, ro_start1)
    };

    // Junction overhang check (when segments don't overlap)
    if overlap == 0 {
        // Non-overlapping case: overhang = segment length at the boundary
        let donor_overhang = ro_donor_end + 1 - ro_coords(tr_donor, read_len).0;
        let acceptor_overhang = ro_coords(tr_acceptor, read_len).1 + 1 - ro_acceptor_start;
        if donor_overhang < overhang_min || acceptor_overhang < overhang_min {
            return Ok(None);
        }
    }

    // Final geometry check: must be truly chimeric (different chr/strand or far apart).
    // STAR: chimeric if chr/strand differ, OR if same-chr same-strand span > alignIntronMax.
    // (For PE inter-mate: > alignMatesGapMax; for SE we use alignIntronMax as the limit.)
    let intron_max = params.align_intron_max as u64;
    let is_chimeric = tr_donor.chr_idx != tr_acceptor.chr_idx
        || tr_donor.is_reverse != tr_acceptor.is_reverse
        || {
            let span = if tr_donor.genome_end <= tr_acceptor.genome_start {
                tr_acceptor.genome_start - tr_donor.genome_end
            } else {
                tr_donor.genome_start.saturating_sub(tr_acceptor.genome_end)
            };
            intron_max > 0 && span > intron_max
        };

    if !is_chimeric {
        return Ok(None);
    }

    // Build chimeric segments
    let donor_seg = transcript_to_segment(tr_donor)
        .map_err(|e| Error::Chimeric(format!("chimeric donor segment: {e}")))?;
    let acceptor_seg = transcript_to_segment(tr_acceptor)
        .map_err(|e| Error::Chimeric(format!("chimeric acceptor segment: {e}")))?;

    // Minimum segment length check
    if !donor_seg.meets_min_length(params.chim_segment_min)
        || !acceptor_seg.meets_min_length(params.chim_segment_min)
    {
        return Ok(None);
    }

    // Classify junction and compute repeats
    let junction_type = classify_junction_type(
        &index.genome,
        donor_seg.chr_idx,
        donor_seg.genome_end,
        donor_seg.is_reverse,
        acceptor_seg.chr_idx,
        acceptor_seg.genome_start,
        acceptor_seg.is_reverse,
    );

    // Apply non-GTAG score penalty and re-check score min
    let effective_score = if junction_type == 0 {
        chim_score + 1 + non_gtag_penalty
    } else {
        chim_score
    };
    if effective_score < score_min || effective_score + score_drop_max < read_len as i32 {
        return Ok(None);
    }

    let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
        &index.genome,
        donor_seg.chr_idx,
        donor_seg.genome_end,
        acceptor_seg.chr_idx,
        acceptor_seg.genome_start,
        20,
    );

    let chim = ChimericAlignment::new(
        donor_seg,
        acceptor_seg,
        junction_type,
        repeat_len_donor,
        repeat_len_acceptor,
        read_seq.to_vec(),
        read_name.to_string(),
    );

    Ok(Some((chim, effective_score)))
}

/// Convert a transcript to a chimeric segment
pub(crate) fn transcript_to_segment(transcript: &Transcript) -> Result<ChimericSegment, Error> {
    if transcript.exons.is_empty() {
        return Err(Error::Alignment(
            "Cannot convert empty transcript to segment".to_string(),
        ));
    }

    // Get overall bounds
    let read_start = transcript.exons[0].read_start;
    let read_end = transcript.exons.last().unwrap().read_end;

    // STAR's junction scan reads a single exon, so carry the outer two through
    // rather than only the transcript-wide span. They coincide for a gapless
    // single-exon segment and diverge as soon as there is a splice or an indel.
    let exon_span = |e: &crate::align::transcript::Exon| ExonSpan {
        genome_start: e.genome_start,
        genome_end: e.genome_end,
        read_start: e.read_start,
        read_end: e.read_end,
    };

    Ok(ChimericSegment {
        chr_idx: transcript.chr_idx,
        genome_start: transcript.genome_start,
        genome_end: transcript.genome_end,
        is_reverse: transcript.is_reverse,
        read_start,
        read_end,
        cigar: transcript.cigar.clone(),
        score: transcript.score,
        n_mismatch: transcript.n_mismatch,
        first_exon: exon_span(&transcript.exons[0]),
        last_exon: exon_span(transcript.exons.last().unwrap()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::WindowAlignment;
    use crate::align::transcript::{Exon, Transcript};
    use crate::genome::{Genome, GenomeSeq};
    use crate::index::GenomeIndex;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;
    use crate::junction::SpliceJunctionDb;
    use noodles::sam::alignment::record::cigar;

    /// Minimal two-chromosome genome for chimeric tests.
    /// Each chromosome has 200 bases of A (=0), padded to 256-byte bins.
    fn make_test_genome() -> Genome {
        let chr_len = 200u64;
        let chr_pad = 256u64;
        let n_genome = chr_pad * 2;
        let sequence = vec![0u8; 2 * n_genome as usize];
        Genome {
            transform_blocks: None,
            sequence: sequence.into(),
            n_genome,
            n_genome_real: n_genome,
            n_chr_real: 2,
            chr_name: vec!["chr0".to_string(), "chr1".to_string()],
            chr_length: vec![chr_len, chr_len],
            chr_start: vec![0, chr_pad, n_genome],
        }
    }

    /// `make_test_index`, with `N` (base 4) planted at the given absolute
    /// positions. The default fixture is all `A`, so a genomic-`N` test has to
    /// put one there deliberately -- the previous test did not, which is why it
    /// could not fail.
    fn make_test_index_with_ns(n_positions: &[u64]) -> GenomeIndex {
        let mut index = make_test_index();
        // The fixture builds an owned sequence, so rebuild it with the Ns in.
        let n = index.genome.n_genome as usize;
        let mut seq = vec![0u8; 2 * n];
        for &p in n_positions {
            seq[p as usize] = 4;
        }
        index.genome.sequence = GenomeSeq::Owned(seq);
        index
    }

    fn make_test_index() -> GenomeIndex {
        let genome = make_test_genome();
        let gstrand_bit = 32u32;
        GenomeIndex {
            genome,
            suffix_array: SuffixArray {
                data: PackedArray::new(33, 0),
                gstrand_bit,
                gstrand_mask: (1u64 << gstrand_bit) - 1,
            },
            sa_index: SaIndex {
                nbases: 0,
                genome_sa_index_start: vec![0],
                data: PackedArray::new(35, 0),
                word_length: 35,
                gstrand_bit,
            },
            junction_db: SpliceJunctionDb::empty(),
            transcriptome: None,
            prepared_junctions: Vec::new(),
            sjdb_overhang: 0,
        }
    }

    fn make_transcript(
        chr_idx: usize,
        genome_start: u64,
        genome_end: u64,
        is_reverse: bool,
    ) -> Transcript {
        use cigar::op::{Kind, Op};
        let read_len = (genome_end - genome_start) as usize;
        Transcript {
            chr_idx,
            genome_start,
            genome_end,
            is_reverse,
            exons: vec![Exon {
                genome_start,
                genome_end,
                read_start: 0,
                read_end: read_len,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, read_len)],
            score: read_len as i32,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        }
    }

    /// Helper to create a minimal SeedCluster for chimeric detection tests
    fn make_test_cluster(
        chr_idx: usize,
        genome_start: u64,
        genome_end: u64,
        is_reverse: bool,
    ) -> SeedCluster {
        SeedCluster {
            alignments: vec![WindowAlignment {
                seed_idx: 0,
                read_pos: 0,
                length: (genome_end - genome_start) as usize,
                genome_pos: genome_start,
                sa_pos: genome_start,
                n_rep: 1,
                is_anchor: true,
                mate_id: 2,
                pre_ext_score: (genome_end - genome_start) as i32,
            }],
            chr_idx,
            genome_start,
            genome_end,
            is_reverse,
            anchor_idx: 0,
            anchor_bin: 0,
        }
    }

    #[test]
    fn test_genomic_distance_same_chr() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 1200, 1300, false);

        assert_eq!(genomic_distance(&c1, &c2), 100);
        assert_eq!(genomic_distance(&c2, &c1), 100);
    }

    #[test]
    fn test_genomic_distance_overlapping() {
        let c1 = make_test_cluster(0, 1000, 1200, false);
        let c2 = make_test_cluster(0, 1100, 1300, false);

        assert_eq!(genomic_distance(&c1, &c2), 0);
    }

    #[test]
    fn test_genomic_distance_different_chr() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(1, 1000, 1100, false);

        assert_eq!(genomic_distance(&c1, &c2), u64::MAX);
    }

    #[test]
    fn test_is_chimeric_signature_different_chr() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(1, 1000, 1100, false);

        assert!(is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_strand_break() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 1200, 1300, true);

        assert!(is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_large_distance() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 2_000_000, 2_000_100, false);

        assert!(is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_close_same_strand() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 1200, 1300, false);

        assert!(!is_chimeric_signature(&c1, &c2));
    }

    // --- transcript_to_segment tests ---

    #[test]
    fn test_transcript_to_segment_basic() {
        let t = make_transcript(0, 1000, 1100, false);
        let seg = transcript_to_segment(&t).unwrap();

        assert_eq!(seg.chr_idx, 0);
        assert_eq!(seg.genome_start, 1000);
        assert_eq!(seg.genome_end, 1100);
        assert!(!seg.is_reverse);
        assert_eq!(seg.read_start, 0);
        assert_eq!(seg.read_end, 100);
        assert_eq!(seg.score, 100);
    }

    #[test]
    fn test_transcript_to_segment_empty_returns_error() {
        let t = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 0,
            is_reverse: false,
            exons: vec![],
            cigar: vec![],
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };
        assert!(transcript_to_segment(&t).is_err());
    }

    // --- detect_inter_mate_chimeric tests ---

    fn params(args: &[&str]) -> Parameters {
        let mut full_args = vec!["rustar-aligner", "--readFilesIn", "reads.fq"];
        full_args.extend_from_slice(args);
        Parameters::parse_from(full_args)
    }

    #[test]
    fn test_inter_mate_chimeric_concordant_returns_none() {
        // Normal FR pair on the same chromosome, close together → not chimeric
        let params = params(&["--chimSegmentMin", "10"]);
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 60, false); // mate1 forward
        let t2 = make_transcript(0, 80, 130, true); // mate2 reverse, same chr, close
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_none());
    }

    #[test]
    fn test_inter_mate_chimeric_different_chromosomes() {
        let params = params(&["--chimSegmentMin", "10"]);
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 60, false); // mate1 chr0
        let t2 = make_transcript(1, 10, 60, true); // mate2 chr1
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_some());
        let chim = result.unwrap();
        // Donor is the mate with earlier read_start (both 0 here; donor is t1 by read_start tie)
        assert_ne!(chim.donor.chr_idx, chim.acceptor.chr_idx);
    }

    #[test]
    fn test_inter_mate_chimeric_same_strand() {
        // Both mates forward on the same chromosome → chimeric (strand break)
        let params = params(&["--chimSegmentMin", "10"]);
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 60, false); // mate1 forward
        let t2 = make_transcript(0, 80, 130, false); // mate2 also forward (abnormal)
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_some());
    }
    /// A chimera whose forward donor scan window covers a genomic `N`.
    ///
    /// The scan runs `jR` in `0..jRmax` from the read-leading segment, so for
    /// this geometry the donor is read at `genome_start + jR` over
    /// `10..=58` and the acceptor at `genome_start + (jR - 30)` over
    /// `70..=118`. Planting an `N` anywhere in either span must ban the
    /// chimera under `banGenomicN` and leave it under `None`.
    #[test]
    fn chim_filter_bans_a_genomic_n_inside_the_scan_window() {
        let banned = params(&["--chimSegmentMin", "10"]);
        let unfiltered = params(&["--chimSegmentMin", "10", "--chimFilter", "None"]);

        // No `N` anywhere: kept either way. Establishes the control.
        let clean_index = make_test_index();
        let chim = forward_chimera();
        assert_eq!(
            apply_chim_filter(vec![chim.clone()], &banned, &clean_index).len(),
            1,
            "a clean junction must survive banGenomicN"
        );

        // `N` inside the donor's span.
        let donor_n = make_test_index_with_ns(&[15]);
        assert_eq!(
            apply_chim_filter(vec![chim.clone()], &banned, &donor_n).len(),
            0,
            "a genomic N in the donor span must ban the chimera"
        );
        assert_eq!(
            apply_chim_filter(vec![chim.clone()], &unfiltered, &donor_n).len(),
            1,
            "--chimFilter None must keep it"
        );

        // `N` inside the acceptor's span (jR = 30 maps to acceptor base 100).
        let acceptor_n = make_test_index_with_ns(&[100]);
        assert_eq!(
            apply_chim_filter(vec![chim.clone()], &banned, &acceptor_n).len(),
            0,
            "a genomic N in the acceptor span must ban the chimera"
        );

        // Outside both spans: not a ban. 200 is past the acceptor window and
        // sits in the second chromosome's padding.
        let far_n = make_test_index_with_ns(&[200]);
        assert_eq!(
            apply_chim_filter(vec![chim], &banned, &far_n).len(),
            1,
            "an N outside the scan window must not ban"
        );
    }

    /// A spliced reverse segment must be scanned from its junction-side exon,
    /// not from the whole-transcript span.
    ///
    /// STAR reads `exons[e0]` with `e0 = Str==1 ? 0 : nExons-1`, so a reverse
    /// leading segment is walked down from its *first* exon's `genome_end`. If
    /// the transcript-wide `genome_end` is used instead, every base examined is
    /// shifted by the intron length, and a gap beside the real junction is
    /// missed while an unrelated position is banned.
    #[test]
    fn chim_filter_scans_the_junction_exon_of_a_spliced_segment() {
        let banned = params(&["--chimSegmentMin", "10"]);
        let mut chim = forward_chimera();
        // Donor becomes reverse and spliced: exon1 [10,25), intron, exon2
        // [125,140). The transcript span stays [10,140), so the two disagree by
        // the 100-base intron.
        chim.donor.is_reverse = true;
        chim.donor.genome_end = 140;
        chim.donor.first_exon = ExonSpan {
            genome_start: 10,
            genome_end: 25,
            read_start: 0,
            read_end: 15,
        };
        chim.donor.last_exon = ExonSpan {
            genome_start: 125,
            genome_end: 140,
            read_start: 15,
            read_end: 30,
        };

        // Reverse + leading => first exon, so the walk starts at 25 - 1 = 24.
        let at_exon = make_test_index_with_ns(&[24]);
        assert_eq!(
            apply_chim_filter(vec![chim.clone()], &banned, &at_exon).len(),
            0,
            "the junction-side exon's end must be the first base scanned"
        );

        // 139 is the transcript-wide genome_end - 1: what the old code scanned.
        // It belongs to the far exon and must not ban.
        let at_transcript_end = make_test_index_with_ns(&[139]);
        assert_eq!(
            apply_chim_filter(vec![chim], &banned, &at_transcript_end).len(),
            1,
            "the transcript-wide span must not be used for a spliced segment"
        );
    }

    /// The strand case, which a position-blind filter gets wrong.
    ///
    /// STAR walks a reverse segment from the far end
    /// (`G[EX_G + EX_L - 1 - jR]`), so for a donor spanning `[10, 40)` the
    /// bases examined run *down* from 39, not up from 10. An `N` at 15 — which
    /// bans the forward case above — must therefore be reached at a different
    /// `jR`, and an `N` planted below the segment must not be reached at all.
    #[test]
    fn chim_filter_reads_a_reverse_donor_from_the_other_end() {
        let banned = params(&["--chimSegmentMin", "10"]);
        let mut chim = forward_chimera();
        chim.donor.is_reverse = true;

        // 39 is the reverse donor's first base (jR = 0). Banned.
        let at_far_end = make_test_index_with_ns(&[39]);
        assert_eq!(
            apply_chim_filter(vec![chim.clone()], &banned, &at_far_end).len(),
            0,
            "reverse donor must be read from genome_end - 1"
        );

        // 5 is below genome_start: the forward walk would never reach it and
        // neither does the reverse walk, which descends from 39.
        let below = make_test_index_with_ns(&[5]);
        assert_eq!(
            apply_chim_filter(vec![chim], &banned, &below).len(),
            1,
            "a position outside the reverse walk must not ban"
        );
    }

    /// The read-`N` ban is unconditional: STAR tests `bR>3` outside the
    /// `filter.genomicN` guard, so `--chimFilter None` does not disable it.
    #[test]
    fn chim_filter_bans_a_read_n_even_with_filtering_off() {
        let index = make_test_index();
        let mut chim = forward_chimera();
        chim.read_seq[5] = 4; // within the scan window, which starts at read 0

        for args in [
            vec!["--chimSegmentMin", "10"],
            vec!["--chimSegmentMin", "10", "--chimFilter", "None"],
        ] {
            let p = params(&args);
            assert_eq!(
                apply_chim_filter(vec![chim.clone()], &p, &index).len(),
                0,
                "an N in the read bans the chimera regardless of --chimFilter ({args:?})"
            );
        }
    }

    /// Shared geometry for the filter tests: donor `[10,40)` over read `[0,30)`,
    /// acceptor `[100,130)` over read `[30,50)`, both forward.
    fn forward_chimera() -> ChimericAlignment {
        ChimericAlignment::new(
            ChimericSegment {
                chr_idx: 0,
                genome_start: 10,
                genome_end: 40,
                read_start: 0,
                read_end: 30,
                is_reverse: false,
                cigar: Vec::new(),
                n_mismatch: 0,
                score: 30,
                first_exon: ExonSpan {
                    genome_start: 10,
                    genome_end: 40,
                    read_start: 0,
                    read_end: 30,
                },
                last_exon: ExonSpan {
                    genome_start: 10,
                    genome_end: 40,
                    read_start: 0,
                    read_end: 30,
                },
            },
            ChimericSegment {
                chr_idx: 0,
                genome_start: 100,
                genome_end: 130,
                read_start: 30,
                read_end: 50,
                is_reverse: false,
                cigar: Vec::new(),
                n_mismatch: 0,
                score: 30,
                first_exon: ExonSpan {
                    genome_start: 100,
                    genome_end: 130,
                    read_start: 30,
                    read_end: 50,
                },
                last_exon: ExonSpan {
                    genome_start: 100,
                    genome_end: 130,
                    read_start: 30,
                    read_end: 50,
                },
            },
            1,
            0,
            0,
            vec![0u8; 50],
            "read1".to_string(),
        )
    }

    #[test]
    fn test_inter_mate_chimeric_too_far() {
        // Opposite-strand pair but >1Mb apart → chimeric
        let params = params(&["--chimSegmentMin", "10"]);
        let index = make_test_index();

        // Use large positions — out-of-bounds for sequence but score.rs guards handle this
        let t1 = make_transcript(0, 10, 60, false);
        let t2 = make_transcript(0, 2_000_000, 2_000_050, true);
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_some());
    }

    #[test]
    fn test_inter_mate_chimeric_segment_too_short() {
        // chimSegmentMin=100 but segments are only 20bp → None
        let params = params(&["--chimSegmentMin", "100"]);
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 30, false);
        let t2 = make_transcript(1, 10, 30, true);
        let read_seq = vec![0u8; 20];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_none());
    }

    #[test]
    fn test_inter_mate_chimeric_empty_exons_returns_none() {
        let params = params(&["--chimSegmentMin", "10"]);
        let index = make_test_index();
        let read_seq = vec![0u8; 50];

        let t1 = make_transcript(0, 10, 60, false);
        let mut t2 = make_transcript(1, 10, 60, true);
        t2.exons.clear();

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_none());
    }

    // --- detect_chimeric_old tests ---

    fn make_read_seq(n: usize) -> Vec<u8> {
        vec![0u8; n]
    }

    // Build a transcript with a soft-clip at one end: the exon covers [left_clip..read_len-right_clip].
    fn make_clipped_transcript(
        chr_idx: usize,
        genome_start: u64,
        is_reverse: bool,
        read_len: usize,
        left_clip: usize,
        right_clip: usize,
    ) -> Transcript {
        use cigar::op::{Kind, Op};
        let aligned_len = read_len - left_clip - right_clip;
        let mut cigar = vec![];
        if left_clip > 0 {
            cigar.push(Op::new(Kind::SoftClip, left_clip));
        }
        cigar.push(Op::new(Kind::Match, aligned_len));
        if right_clip > 0 {
            cigar.push(Op::new(Kind::SoftClip, right_clip));
        }
        Transcript {
            chr_idx,
            genome_start,
            genome_end: genome_start + aligned_len as u64,
            is_reverse,
            exons: vec![Exon {
                genome_start,
                genome_end: genome_start + aligned_len as u64,
                read_start: left_clip,
                read_end: left_clip + aligned_len,
                i_frag: 0,
            }],
            cigar,
            score: aligned_len as i32,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        }
    }

    #[test]
    fn test_detect_chimeric_old_no_chimera_single_transcript() {
        // Only one transcript → no partner → None
        let params = params(&["--chimSegmentMin", "20"]);
        let index = make_test_index();
        let read_len = 100usize;
        let t1 = make_clipped_transcript(0, 50, false, read_len, 0, 0);
        let read_seq = make_read_seq(read_len);
        let result = detect_chimeric_old(
            std::slice::from_ref(&t1),
            &t1,
            &read_seq,
            "r",
            &params,
            &index,
        )
        .unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_detect_chimeric_old_inter_chr_pair() {
        // Primary: covers read[0..80] on chr0; secondary: covers read[80..100] on chr1.
        // With chimSegmentMin=20, scoreDropMax=20, scoreSeparation=10, this should produce a chimera.
        let params = params(&[
            "--chimSegmentMin",
            "15",
            "--chimScoreDropMax",
            "100",
            "--chimScoreSeparation",
            "10",
            "--chimJunctionOverhangMin",
            "10",
        ]);
        let index = make_test_index();
        let read_len = 100usize;
        // Primary: chr0, read[0..80], right clip = 20
        let t_main = make_clipped_transcript(0, 0, false, read_len, 0, 20);
        // Partner: chr1, read[80..100], left clip = 80
        let t_partner = make_clipped_transcript(1, 0, false, read_len, 80, 0);

        let all = vec![t_main.clone(), t_partner];
        let result =
            detect_chimeric_old(&all, &t_main, &read_seq_n(read_len), "r", &params, &index)
                .unwrap();
        // Should find a chimeric alignment
        assert_eq!(result.len(), 1);
        let chim = &result[0];
        assert_ne!(chim.donor.chr_idx, chim.acceptor.chr_idx);
    }

    fn read_seq_n(n: usize) -> Vec<u8> {
        vec![0u8; n]
    }

    #[test]
    fn test_detect_chimeric_old_segment_too_short() {
        // Segments are too short after chimSegmentMin filter
        let params = params(&["--chimSegmentMin", "50", "--chimScoreDropMax", "100"]);
        let index = make_test_index();
        let read_len = 100usize;
        let t_main = make_clipped_transcript(0, 0, false, read_len, 0, 60); // 40 bp → < 50
        let t_partner = make_clipped_transcript(1, 0, false, read_len, 60, 0); // 40 bp → < 50
        let all = vec![t_main.clone(), t_partner];
        let result =
            detect_chimeric_old(&all, &t_main, &read_seq_n(read_len), "r", &params, &index)
                .unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_detect_chimeric_old_score_drop_too_large() {
        // Score drop is too large: chimScoreDropMax=5 means combined_score + 5 >= read_len=100
        // combined_score = 50 + 50 - 0 = 100, 100 + 5 = 105 >= 100 → should pass
        // But with drop=5, score = 40+40=80, 80+5=85 < 100 → should fail
        let params = params(&[
            "--chimSegmentMin",
            "20",
            "--chimScoreDropMax",
            "5",
            "--chimScoreSeparation",
            "200", // suppress uniqueness filter
        ]);
        let index = make_test_index();
        let read_len = 100usize;
        let t_main = make_clipped_transcript(0, 0, false, read_len, 0, 40); // 60 bp aligned
        let t_partner = make_clipped_transcript(1, 0, false, read_len, 60, 0); // 40 bp aligned
        // combined_score = 60 + 40 = 100, 100 + 5 = 105 >= 100 → OK, should pass
        let all = vec![t_main.clone(), t_partner];
        let result =
            detect_chimeric_old(&all, &t_main, &read_seq_n(read_len), "r", &params, &index)
                .unwrap();
        // Score drop filter: 100 + 5 >= 100 → passes; uniqueness: score_separation=200, next=-inf → passes
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_detect_chimeric_old_diff_mates_waives_gap() {
        // Combined read len 100, mate boundary at 50: mate1=[0..50), mate2=[50..100).
        // Primary covers read[0..40] (mate1, chr0); partner covers read[60..100] (mate2, chr1).
        // The 20bp inter-segment gap exceeds chimSegmentReadGapMax (0), so the standard
        // gap check fails — but the segments are in different mates, so STAR's diffMates
        // waives it. Without a mate boundary the chimera is rejected; with one it is found.
        let params = params(&[
            "--chimSegmentMin",
            "15",
            "--chimScoreDropMax",
            "100",
            "--chimScoreSeparation",
            "200",
            "--chimJunctionOverhangMin",
            "10",
        ]);
        let index = make_test_index();
        let read_len = 100usize;
        let t_main = make_clipped_transcript(0, 0, false, read_len, 0, 60); // read[0..40], mate1
        let t_partner = make_clipped_transcript(1, 0, false, read_len, 60, 0); // read[60..100], mate2
        let all = vec![t_main.clone(), t_partner];
        let read_seq = read_seq_n(read_len);

        // No boundary (SE / per-mate pool): gap check enforced → rejected.
        let without = detect_chimeric_old(&all, &t_main, &read_seq, "r", &params, &index).unwrap();
        assert!(
            without.is_empty(),
            "gap check should reject without diffMates"
        );

        // Combined read with mate boundary at 50: diffMates waives the gap → chimera found.
        let with =
            detect_chimeric_old_impl(&all, &t_main, &read_seq, "r", &params, &index, Some(50))
                .unwrap();
        assert_eq!(with.len(), 1, "diffMates should waive the inter-mate gap");
        assert_ne!(with[0].donor.chr_idx, with[0].acceptor.chr_idx);
    }

    // --- detect_chimeric_mult tests ---

    /// A read whose 3' end maps equally well to two places: the same donor, two
    /// acceptors, identical scores.
    fn two_locus_pool(read_len: usize) -> Vec<Transcript> {
        vec![
            make_clipped_transcript(0, 0, false, read_len, 0, 30), // read[0..70] on chr0
            make_clipped_transcript(1, 0, false, read_len, 70, 0), // read[70..100] on chr1
            make_clipped_transcript(1, 150, false, read_len, 70, 0), // ...and again, elsewhere
        ]
    }

    fn mult_params(extra: &[&str]) -> Parameters {
        let mut args = vec![
            "--chimSegmentMin",
            "15",
            "--chimScoreDropMax",
            "100",
            "--chimJunctionOverhangMin",
            "10",
            "--chimNonchimScoreDropMin",
            "5",
        ];
        args.extend_from_slice(extra);
        params(&args)
    }

    /// The point of the whole path: a read with two equally good chimeric loci
    /// yields two junctions, where the old single-best path can only ever name
    /// one of them and silently discards the other.
    #[test]
    fn chim_multimap_reports_every_locus_within_the_score_range() {
        let index = make_test_index();
        let read_len = 100usize;
        let pool = two_locus_pool(read_len);
        let read_seq = read_seq_n(read_len);

        let old = detect_chimeric_old(
            &pool,
            &pool[0],
            &read_seq,
            "r",
            &mult_params(&["--chimScoreSeparation", "10"]),
            &index,
        )
        .unwrap();
        assert_eq!(
            old.len(),
            1,
            "the old path reports one locus and drops the equally good one"
        );

        let mult = detect_chimeric_mult(
            &pool,
            pool[0].score,
            &read_seq,
            "r",
            &mult_params(&["--chimMultimapNmax", "10"]),
            &index,
            None,
        )
        .unwrap();
        assert_eq!(mult.len(), 2, "both loci should be reported");
        for chim in &mult {
            assert_ne!(chim.donor.chr_idx, chim.acceptor.chr_idx);
        }
    }

    /// Past the cap STAR reports nothing at all, not the first `nmax`.
    /// A banned locus must not be counted against `--chimMultimapNmax`.
    ///
    /// STAR zeroes a candidate's score inside stitching when it meets an `N`
    /// (`ChimericAlign_chimericStitching.cpp:63-66`), so it never reaches the
    /// list the cap is measured against
    /// (`ChimericDetection_chimericDetectionMult.cpp:79-95`). Filtering after
    /// selection instead makes a banned locus push the survivor count over the
    /// cap, and STAR's "over the cap, report nothing" rule then throws away the
    /// clean locus too.
    #[test]
    fn chim_multimap_does_not_count_a_banned_locus_against_the_cap() {
        let read_len = 100usize;
        let pool = two_locus_pool(read_len);
        let read = read_seq_n(read_len);
        let cap1 = mult_params(&["--chimMultimapNmax", "1"]);

        // Both loci clean: two survivors over a cap of 1, so nothing is
        // reported. This is the control that the cap is live.
        let clean = make_test_index();
        assert!(
            detect_chimeric_mult(&pool, pool[0].score, &read, "r", &cap1, &clean, None)
                .unwrap()
                .is_empty(),
            "two clean loci over a cap of 1 must report nothing"
        );

        // Ban only the second acceptor locus. With the donor on read[0..70] and
        // the acceptors on read[70..100], the scan spans are: donor 0..98,
        // acceptor-at-0 -70..28, acceptor-at-150 80..178. So 160 sits inside the
        // second acceptor's span and outside both the donor's and the first
        // acceptor's, banning exactly one candidate. One survivor remains, which
        // is within the cap, so it must be reported rather than the whole read
        // being discarded.
        let banned_second = make_test_index_with_ns(&[160]);
        let got = detect_chimeric_mult(
            &pool,
            pool[0].score,
            &read,
            "r",
            &cap1,
            &banned_second,
            None,
        )
        .unwrap();
        assert_eq!(
            got.len(),
            1,
            "a banned locus must not count towards --chimMultimapNmax"
        );
    }

    #[test]
    fn chim_multimap_beyond_the_cap_reports_nothing() {
        let index = make_test_index();
        let read_len = 100usize;
        let pool = two_locus_pool(read_len);
        let mult = detect_chimeric_mult(
            &pool,
            pool[0].score,
            &read_seq_n(read_len),
            "r",
            &mult_params(&["--chimMultimapNmax", "1"]),
            &index,
            None,
        )
        .unwrap();
        assert!(
            mult.is_empty(),
            "2 loci over a cap of 1 must report nothing, not a truncated list"
        );
    }

    /// A locus more than `--chimMultimapScoreRange` below the best is dropped,
    /// and dropping it leaves a single unambiguous chimera.
    #[test]
    fn chim_multimap_score_range_drops_the_weaker_locus() {
        let index = make_test_index();
        let read_len = 100usize;
        let mut pool = two_locus_pool(read_len);
        // Make the second acceptor cover 5 fewer read bases, so its chimeric
        // score is 5 lower.
        pool[2] = make_clipped_transcript(1, 150, false, read_len, 70, 5);

        let mult = detect_chimeric_mult(
            &pool,
            pool[0].score,
            &read_seq_n(read_len),
            "r",
            &mult_params(&["--chimMultimapNmax", "10", "--chimMultimapScoreRange", "1"]),
            &index,
            None,
        )
        .unwrap();
        assert_eq!(mult.len(), 1, "only the best locus is within the range");
        assert_eq!(mult[0].acceptor.read_start, 70);

        // Widen the range and the weaker locus comes back.
        let wide = detect_chimeric_mult(
            &pool,
            pool[0].score,
            &read_seq_n(read_len),
            "r",
            &mult_params(&["--chimMultimapNmax", "10", "--chimMultimapScoreRange", "10"]),
            &index,
            None,
        )
        .unwrap();
        assert_eq!(wide.len(), 2);
    }

    /// `--chimNonchimScoreDropMin`: a read that already aligns linearly across
    /// its whole length is not a fusion candidate, however well the two halves
    /// score on their own.
    #[test]
    fn chim_multimap_requires_the_linear_alignment_to_leave_the_read_unexplained() {
        let index = make_test_index();
        let read_len = 100usize;
        let pool = two_locus_pool(read_len);
        let mult = detect_chimeric_mult(
            &pool,
            read_len as i32, // a full-length linear alignment
            &read_seq_n(read_len),
            "r",
            &mult_params(&["--chimMultimapNmax", "10"]),
            &index,
            None,
        )
        .unwrap();
        assert!(mult.is_empty());
    }

    /// The knob is opt-in: at its default of 0 the path is inert.
    #[test]
    fn chim_multimap_is_off_by_default() {
        let index = make_test_index();
        let read_len = 100usize;
        let pool = two_locus_pool(read_len);
        let mult = detect_chimeric_mult(
            &pool,
            pool[0].score,
            &read_seq_n(read_len),
            "r",
            &mult_params(&[]),
            &index,
            None,
        )
        .unwrap();
        assert!(mult.is_empty());
    }
}
