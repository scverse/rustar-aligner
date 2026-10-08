/// Scoring functions for alignment gaps and splice junctions
use crate::genome::{Genome, SeqView};
use crate::params::Parameters;

/// Alignment scorer with user-defined penalties
#[derive(Debug, Clone)]
pub struct AlignmentScorer {
    /// Canonical splice junction penalty (GT-AG)
    pub score_gap: i32,
    /// Non-canonical splice junction penalty
    pub score_gap_noncan: i32,
    /// GC-AG splice junction penalty
    pub score_gap_gcag: i32,
    /// AT-AC splice junction penalty
    pub score_gap_atac: i32,
    /// Deletion open penalty
    pub score_del_open: i32,
    /// Deletion extension penalty (per base)
    pub score_del_base: i32,
    /// Insertion open penalty
    pub score_ins_open: i32,
    /// Insertion extension penalty (per base)
    pub score_ins_base: i32,
    /// Minimum intron length (gaps >= this are treated as splice junctions)
    pub align_intron_min: u32,
    /// Bonus for annotated splice junctions (from GTF)
    pub sjdb_score: i32,
    /// Max mismatches for stitching SJs: [non-canonical, GT-AG, GC-AG, AT-AC]
    /// -1 means unlimited
    pub align_sj_stitch_mismatch_nmax: [i32; 4],
    /// Max absolute number of mismatches for alignment extension (outFilterMismatchNmax)
    pub n_mm_max: u32,
    /// Max ratio of mismatches to total alignment length (outFilterMismatchNoverLmax)
    pub p_mm_max: f64,
    /// Minimum overhang for splice junctions (alignSJoverhangMin, default 5)
    pub align_sj_overhang_min: u32,
    /// Minimum overhang for annotated splice junctions (alignSJDBoverhangMin, default 3)
    pub align_sjdb_overhang_min: u32,
    /// Maximum intron length. 0 = no limit (STAR-faithful: alignIntronMax=0 disables the check).
    /// Use u32::MAX when the CLI param is 0. Set to a finite value when user specifies a limit.
    pub align_intron_max: u32,
    /// Extra score log-scaled with genomic length: scale * log2(genomicLength)
    pub score_genomic_length_log2_scale: f64,
    /// Max score reduction for SJ stitching shift (scoreStitchSJshift, default 1)
    pub score_stitch_sj_shift: i32,
    /// Min mapped length of spliced mates (absolute, default 0)
    pub align_spliced_mate_map_lmin: u32,
    /// Min mapped length of spliced mates as fraction of read length (default 0.66)
    pub align_spliced_mate_map_lmin_over_lmate: f64,
    /// Minimum alignment score relative to read length (outFilterScoreMinOverLread, default 0.66)
    pub out_filter_score_min_over_lread: f64,
    /// Read-end extension policy (alignEndsType). `ext[iMate][iEnd]==true` forces
    /// full end-to-end extension (no terminal soft-clip) of that mate/end.
    pub align_ends_type: crate::params::AlignEndsType,
}

impl AlignmentScorer {
    /// Create a minimal scorer for motif detection only (used in junction recording)
    pub fn from_params_minimal() -> Self {
        Self {
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
            align_ends_type: crate::params::AlignEndsType::default(),
        }
    }

    /// Create scorer from parameters
    pub fn from_params(params: &Parameters) -> Self {
        Self {
            score_gap: params.score_gap,
            score_gap_noncan: params.score_gap_noncan,
            score_gap_gcag: params.score_gap_gcag,
            score_gap_atac: params.score_gap_atac,
            score_del_open: params.score_del_open,
            score_del_base: params.score_del_base,
            score_ins_open: params.score_ins_open,
            score_ins_base: params.score_ins_base,
            align_intron_min: params.align_intron_min,
            sjdb_score: params.sjdb_score,
            align_sj_stitch_mismatch_nmax: [
                params.align_sj_stitch_mismatch_nmax[0],
                params.align_sj_stitch_mismatch_nmax[1],
                params.align_sj_stitch_mismatch_nmax[2],
                params.align_sj_stitch_mismatch_nmax[3],
            ],
            n_mm_max: params.out_filter_mismatch_nmax,
            p_mm_max: params.out_filter_mismatch_nover_lmax,
            align_sj_overhang_min: params.align_sj_overhang_min,
            align_sjdb_overhang_min: params.align_sjdb_overhang_min,
            // STAR: when alignIntronMax==0 the check `Del>alignIntronMax && alignIntronMax>0`
            // is never true, so all intron sizes are allowed. Mirror with u32::MAX sentinel.
            align_intron_max: if params.align_intron_max == 0 {
                u32::MAX
            } else {
                params.align_intron_max
            },
            score_genomic_length_log2_scale: params.score_genomic_length_log2_scale,
            score_stitch_sj_shift: params.score_stitch_sj_shift,
            align_spliced_mate_map_lmin: params.align_spliced_mate_map_lmin,
            align_spliced_mate_map_lmin_over_lmate: params.align_spliced_mate_map_lmin_over_lmate,
            out_filter_score_min_over_lread: params.out_filter_score_min_over_lread,
            // Parsed+validated in Parameters::validate; default to Local if unset.
            align_ends_type: params.align_ends_type.parse().unwrap_or_default(),
        }
    }

    /// Compute genomic length penalty for a transcript.
    /// STAR formula: ceil(log2(genomicLength) * scale - 0.5), clamped so score >= 0.
    pub fn genomic_length_penalty(&self, genomic_span: u64) -> i32 {
        if self.score_genomic_length_log2_scale == 0.0 || genomic_span == 0 {
            return 0;
        }
        ((genomic_span as f64).log2() * self.score_genomic_length_log2_scale - 0.5).ceil() as i32
    }

    /// Annotated junctions score `sjdb_score`; unannotated junctions score `motif_score`.
    pub fn score_annotated_junction(&self, motif_score: i32, annotated: bool) -> i32 {
        if annotated {
            self.sjdb_score
        } else {
            motif_score
        }
    }

    /// Check if the number of mismatches is allowed for this junction motif type
    pub fn stitch_mismatch_allowed(&self, motif: &SpliceMotif, n_mismatch: u32) -> bool {
        let idx = match motif {
            SpliceMotif::NonCanonical => 0,
            SpliceMotif::GtAg | SpliceMotif::CtAc => 1,
            SpliceMotif::GcAg | SpliceMotif::CtGc => 2,
            SpliceMotif::AtAc | SpliceMotif::GtAt => 3,
        };
        let max_mm = self.align_sj_stitch_mismatch_nmax[idx];
        max_mm < 0 || n_mismatch <= max_mm as u32
    }

    /// Score a gap between two aligned regions
    ///
    /// # Arguments
    /// - `genome_gap`: Gap in genome coordinates (can be negative for insertions)
    /// - `read_gap`: Gap in read coordinates
    /// - `genome_pos`: Starting genomic position (for motif detection)
    /// - `genome`: Genome reference (for motif detection)
    ///
    /// # Returns
    /// `(score, gap_type)`
    pub fn score_gap(
        &self,
        genome_gap: i64,
        read_gap: i64,
        genome_pos: u64,
        genome: &Genome,
    ) -> (i32, GapType) {
        self.score_gap_with_strand(genome_gap, read_gap, genome_pos, genome, false, 0)
    }

    /// Score a gap between consecutive seeds, with strand-aware motif detection.
    ///
    /// For reverse-strand reads, `genome_pos` is a raw SA position in the RC genome.
    /// The donor position must be converted to forward genome coordinates for motif
    /// detection: `forward_donor = n_genome - rc_donor - intron_len`.
    pub fn score_gap_with_strand(
        &self,
        genome_gap: i64,
        read_gap: i64,
        genome_pos: u64,
        genome: &Genome,
        is_reverse: bool,
        n_genome: u64,
    ) -> (i32, GapType) {
        match (genome_gap, read_gap) {
            // Insertion: read advances but genome doesn't
            (0, rg) if rg > 0 => {
                let len = rg as u32;
                let score = self.score_ins_open + self.score_ins_base * len as i32;
                (score, GapType::Insertion(len))
            }
            // Deletion or splice junction: genome advances but read doesn't
            (gg, 0) if gg > 0 => {
                let len = gg as u32;
                if len >= self.align_intron_min && len <= self.align_intron_max {
                    // Splice junction — detect motif on forward genome
                    let donor = if is_reverse {
                        n_genome - genome_pos - len as u64
                    } else {
                        genome_pos
                    };
                    let motif = self.detect_splice_motif(donor, len, genome);
                    let score = self.score_splice_junction(motif);
                    (
                        score,
                        GapType::SpliceJunction {
                            intron_len: len,
                            motif,
                        },
                    )
                } else {
                    // Deletion (too short for intron, or exceeds max intron length)
                    let score = self.score_del_open + self.score_del_base * len as i32;
                    (score, GapType::Deletion(len))
                }
            }
            // Both advance: combined gap (handled by CIGAR builder in stitch.rs)
            // Score the net indel portion
            (gg, rg) if gg > 0 && rg > 0 => {
                let excess = gg - rg;
                match excess {
                    1.. => {
                        let del_len = excess as u32;
                        if del_len >= self.align_intron_min && del_len <= self.align_intron_max {
                            let rc_donor = genome_pos + rg as u64;
                            let donor = if is_reverse {
                                n_genome - rc_donor - del_len as u64
                            } else {
                                rc_donor
                            };
                            let motif = self.detect_splice_motif(donor, del_len, genome);
                            let score = self.score_splice_junction(motif);
                            (
                                score,
                                GapType::SpliceJunction {
                                    intron_len: del_len,
                                    motif,
                                },
                            )
                        } else {
                            let score = self.score_del_open + self.score_del_base * del_len as i32;
                            (score, GapType::Deletion(del_len))
                        }
                    }
                    ..=-1 => {
                        let ins_len = (-excess) as u32;
                        let score = self.score_ins_open + self.score_ins_base * ins_len as i32;
                        (score, GapType::Insertion(ins_len))
                    }
                    0 => {
                        // Equal gaps: no net indel
                        (0, GapType::Deletion(0))
                    }
                }
            }
            // Other cases (negative gaps, etc.)
            _ => (0, GapType::Deletion(0)),
        }
    }

    /// Memoized wrapper around [`AlignmentScorer::find_best_junction_position`].
    ///
    /// The scan is a pure function of its arguments. Within one window's stitch
    /// recursion, `read_seq`, the genome, `is_reverse` and `n_genome` are fixed,
    /// so the six remaining coordinates identify a scan completely.
    /// `stitchWindowAligns`' include/exclude recursion reaches the same
    /// (exon A end, seed B) pair through many different branch paths, so without
    /// a memo the identical scan is repeated thousands of times per window.
    /// Results are bit-identical to calling the uncached function.
    #[allow(clippy::too_many_arguments)]
    pub fn find_best_junction_position_cached(
        &self,
        cache: &mut JunctionScanCache,
        read_seq: &[u8],
        r_a_end: usize,
        g_a_end: u64,
        r_gap: i64,
        g_gap: i64,
        genome: &Genome,
        is_reverse: bool,
        n_genome: u64,
        prev_exon_len: usize,
        next_seed_len: usize,
    ) -> (i32, SpliceMotif, i32, u32, u32) {
        let key = JunctionScanKey {
            r_a_end,
            g_a_end,
            r_gap,
            g_gap,
            prev_exon_len,
            next_seed_len,
        };
        if let Some(hit) = cache.map.get(&key) {
            return *hit;
        }
        let val = self.find_best_junction_position(
            read_seq,
            r_a_end,
            g_a_end,
            r_gap,
            g_gap,
            genome,
            is_reverse,
            n_genome,
            prev_exon_len,
            next_seed_len,
        );
        cache.map.insert(key, val);
        val
    }

    /// Find the optimal junction boundary position by scanning all candidates.
    ///
    /// STAR's jR scanning: given a gap between seeds A and B where gGap > rGap,
    /// slide the junction boundary through all valid positions. At each position,
    /// score how well the read matches the upstream vs downstream genome, plus
    /// the splice motif quality. Return the shift that maximizes the combined score.
    ///
    /// All coordinates are in SA coordinate space (forward for fwd reads, RC genome for rev).
    ///
    /// Returns: (jr_shift, best_motif, best_motif_score)
    /// - jr_shift: how many bases to shift the junction boundary (positive = rightward)
    /// - best_motif: the splice motif at the optimal position
    /// - best_motif_score: the motif penalty at the optimal position
    #[allow(clippy::too_many_arguments)]
    pub fn find_best_junction_position(
        &self,
        read_seq: &[u8],
        r_a_end: usize, // prev.read_end (exclusive, rustar-aligner convention)
        g_a_end: u64,   // prev.genome_end (exclusive, rustar-aligner convention)
        r_gap: i64,     // read gap between seeds
        g_gap: i64,     // genome gap between seeds
        genome: &Genome,
        is_reverse: bool,
        n_genome: u64,
        prev_exon_len: usize,
        next_seed_len: usize, // length of the B seed (for scan range)
    ) -> (i32, SpliceMotif, i32, u32, u32) {
        let del = g_gap - r_gap; // net intron/deletion length (constant)
        debug_assert!(del > 0);

        // Convert to STAR-style inclusive coordinates for the scanning algorithm
        // rustar-aligner: r_a_end is exclusive (one past last base of seed A)
        // STAR:   rAend is inclusive (last base of seed A)
        let g_a_end_inc = g_a_end - 1; // last genome base of seed A (inclusive)
        let r_a_end_inc = r_a_end - 1; // last read base of seed A (inclusive)

        // gBstart1 = position in genome corresponding to the acceptor side at jR=0
        // In STAR: gBstart1 = gAend + gGap - rGap (= gAend + Del)
        // With inclusive coords: gBstart1 = g_a_end_inc + del
        let g_b_start1 = g_a_end_inc as i64 + del;

        let genome_offset: u64 = if is_reverse { n_genome } else { 0 };
        // Resolve the genome storage once: the three scans below read a base per
        // iteration, and `Genome::get_base` re-checks the `GenomeSeq` variant on
        // every one of them.
        let seq = genome.sequence.view();

        // Phase 1: Move LEFT from jR1=1, scoring mismatches
        // Find how far left we need to start scanning
        let mut jr1: i32 = 1;
        let mut score1: i32 = 0;
        loop {
            jr1 -= 1;
            let ri = r_a_end_inc as i64 + jr1 as i64;
            if ri < 0 || ri >= read_seq.len() as i64 {
                break;
            }
            let read_base = read_seq[ri as usize];

            let g_up_pos = g_a_end_inc as i64 + jr1 as i64;
            let g_dn_pos = g_b_start1 + jr1 as i64;
            if g_up_pos < 0 || g_dn_pos < 0 {
                break;
            }

            let g_up = seq.base((g_up_pos as u64 + genome_offset) as usize);
            let g_dn = seq.base((g_dn_pos as u64 + genome_offset) as usize);

            if g_up >= 4 || g_dn >= 4 {
                break;
            }
            if read_base == g_up && read_base != g_dn {
                // Moving left costs: this base matches upstream but not downstream
                score1 -= 1;
            }

            if score1 + self.score_stitch_sj_shift < 0 {
                break;
            }
            if prev_exon_len as i32 + jr1 <= 1 {
                break;
            }
        }
        // jr1 is now one past where we stopped; the scan will start from jr1

        // Phase 2: Scan RIGHT through all jR1 positions
        score1 = 0;
        let mut max_score2 = i32::MIN;
        let mut best_jr: i32 = 0;
        let mut best_motif = SpliceMotif::NonCanonical;
        let mut best_motif_score = self.score_gap_noncan;

        // `del` does not change across the scan, so whether the motif branch
        // runs is decided once rather than per iteration.
        let motif_in_range =
            del >= self.align_intron_min as i64 && del <= self.align_intron_max as i64;

        // The four motif bases sit at `donor`, `donor+1`, `donor+del-2` and
        // `donor+del-1`, and `donor` moves exactly one base per iteration, so
        // consecutive iterations share two of the four. Carrying the window
        // instead of re-fetching it turns four genome reads per position into
        // two. `window` is `None` until the first motif iteration primes it.
        let mut window: Option<MotifWindow> = None;

        loop {
            let ri = r_a_end_inc as i64 + jr1 as i64;
            if ri >= 0 && (ri as usize) < read_seq.len() {
                let read_base = read_seq[ri as usize];
                let g_up_pos = g_a_end_inc as i64 + jr1 as i64;
                let g_dn_pos = g_b_start1 + jr1 as i64;

                if g_up_pos >= 0 && g_dn_pos >= 0 {
                    let gu = seq.base((g_up_pos as u64 + genome_offset) as usize);
                    let gd = seq.base((g_dn_pos as u64 + genome_offset) as usize);

                    if gu < 4 && gd < 4 {
                        if read_base == gu && read_base != gd {
                            score1 += 1;
                        } else if read_base != gu && read_base == gd {
                            score1 -= 1;
                        }
                    }
                }
            }

            // Check splice motif at this junction position
            if motif_in_range {
                // Donor position in SA space: one past the last donor-exon base
                let donor_sa = (g_a_end_inc as i64 + jr1 as i64 + 1) as u64;
                // Convert to forward genome coordinates for motif detection
                let donor_fwd = if is_reverse {
                    n_genome - donor_sa - del as u64
                } else {
                    donor_sa
                };
                let w = match window.as_mut() {
                    Some(w) => {
                        w.slide_to(donor_fwd, del as u64, seq);
                        &*w
                    }
                    None => window.insert(MotifWindow::at(donor_fwd, del as u64, seq)),
                };
                let motif = w.motif();
                let motif_score = self.score_splice_junction(motif);
                let score2 = score1 + motif_score;

                if score2 > max_score2 {
                    max_score2 = score2;
                    best_jr = jr1;
                    best_motif = motif;
                    best_motif_score = motif_score;
                }
            } else {
                // Deletion (Del < alignIntronMin) or out-of-range gap — no motif,
                // pure positional score. STAR: jCan1=-1, jPen1=0, Score2=Score1.
                // For reverse strand, rightmost-in-RC wins ties (= leftmost in forward).
                if score1 > max_score2 || (is_reverse && score1 == max_score2) {
                    max_score2 = score1;
                    best_jr = jr1;
                    best_motif = SpliceMotif::NonCanonical;
                    best_motif_score = 0;
                }
            }

            jr1 += 1;
            // STAR: jR1 < int(rBend) - int(rAend) where rBend = rBstart + L - 1
            // This equals jR1 < rGap + L, allowing the scan to extend into the B seed.
            // This is critical for finding canonical junctions when seeds are adjacent
            // (rGap=0) but the junction needs to shift into the B seed's territory.
            if jr1 >= r_gap as i32 + next_seed_len as i32 {
                break;
            }
        }

        // Phase 3: Repeat detection around best_jr + left-flush for non-canonical/deletion
        // Count matching bases between donor and acceptor sides around best_jr
        // STAR: jjL = left repeat, jjR = right repeat
        let mut jj_l: i32 = 0;
        loop {
            let left_pos = g_a_end_inc as i64 + best_jr as i64 - jj_l as i64;
            let right_pos = g_b_start1 + best_jr as i64 - jj_l as i64;
            if left_pos < 0 || right_pos < 0 {
                break;
            }
            let gl = seq.base((left_pos as u64 + genome_offset) as usize);
            let gr = seq.base((right_pos as u64 + genome_offset) as usize);
            if gl < 4 && gl == gr {
                jj_l += 1;
            } else {
                break;
            }
            if jj_l > 255 {
                break;
            }
        }

        // Right repeat: count matching bases going right from junction
        let mut jj_r: i32 = 0;
        loop {
            let left_pos = g_a_end_inc as i64 + jj_r as i64 + best_jr as i64 + 1;
            let right_pos = g_b_start1 + jj_r as i64 + best_jr as i64 + 1;
            if left_pos < 0 || right_pos < 0 {
                break;
            }
            let gl = seq.base((left_pos as u64 + genome_offset) as usize);
            let gr = seq.base((right_pos as u64 + genome_offset) as usize);
            if gl < 4 && gl == gr {
                jj_r += 1;
            } else {
                break;
            }
            if jj_r > 255 {
                break;
            }
        }

        // Non-canonical or deletion: flush to be deterministic in repeat regions.
        // STAR operates in forward-genome space, so flush-left = leftmost in forward coords.
        // Our function operates in SA space (RC genome for reverse reads), so:
        //   - Forward reads: flush LEFT in SA space = LEFT in forward space ✓
        //   - Reverse reads: flush RIGHT in SA space = LEFT in forward space ✓
        if best_motif == SpliceMotif::NonCanonical || del < self.align_intron_min as i64 {
            if is_reverse {
                // Reverse strand: flush RIGHT in RC space = LEFT in forward space
                best_jr += jj_r;
                jj_l += jj_r;
                jj_r = 0;
            } else {
                // Forward strand: flush LEFT (STAR default)
                best_jr -= jj_l;
                jj_r += jj_l;
                jj_l = 0;
            }
            // STAR: if (int(EX_L)+jR<1) return -1000005;
            // Clamp: don't let exon A become zero-length (STAR rejects, we clamp)
            best_jr = best_jr.max(1 - prev_exon_len as i32);
            // Re-check motif at flushed position
            if del >= self.align_intron_min as i64 && del <= self.align_intron_max as i64 {
                let donor_sa = (g_a_end_inc as i64 + best_jr as i64 + 1) as u64;
                let donor_fwd = if is_reverse {
                    n_genome - donor_sa - del as u64
                } else {
                    donor_sa
                };
                best_motif = self.detect_splice_motif(donor_fwd, del as u32, genome);
                best_motif_score = self.score_splice_junction(best_motif);
            }
        }

        (
            best_jr,
            best_motif,
            best_motif_score,
            jj_l.max(0) as u32,
            jj_r.max(0) as u32,
        )
    }

    /// Detect splice junction motif (thin wrapper over the free function
    /// so `AlignmentScorer` callers keep working).
    pub fn detect_splice_motif(
        &self,
        donor_pos: u64,
        intron_len: u32,
        genome: &Genome,
    ) -> SpliceMotif {
        detect_splice_motif(donor_pos, intron_len, genome)
    }

    /// Score a splice junction based on motif
    pub(crate) fn score_splice_junction(&self, motif: SpliceMotif) -> i32 {
        match motif {
            SpliceMotif::GtAg | SpliceMotif::CtAc => self.score_gap,
            SpliceMotif::GcAg | SpliceMotif::CtGc => self.score_gap_gcag,
            SpliceMotif::AtAc | SpliceMotif::GtAt => self.score_gap_atac,
            SpliceMotif::NonCanonical => self.score_gap_noncan,
        }
    }
}

/// Detect splice junction motif from forward-strand bases at the intron
/// boundaries. Stateless — exposed as a free function so both alignment
/// scoring and `genomeGenerate` splice-junction insertion can share one
/// truth table.
///
/// `donor_pos` is the 0-based position of the intron's first base on the
/// forward strand; `intron_len` is the intron length in bases.
pub fn detect_splice_motif(donor_pos: u64, intron_len: u32, genome: &Genome) -> SpliceMotif {
    MotifWindow::at(donor_pos, intron_len as u64, genome.sequence.view()).motif()
}

/// The four bases that decide a splice motif: the intron's first two and last
/// two, on the forward strand.
///
/// Kept as a struct so the junction scan can slide it. Successive junction
/// positions differ by one base, so three of the four positions overlap the
/// previous window and only two bases have to be read from the genome.
#[derive(Clone, Copy)]
struct MotifWindow {
    donor: u64,
    d1: u8,
    d2: u8,
    a1: u8,
    a2: u8,
}

impl MotifWindow {
    #[inline]
    fn at(donor: u64, intron_len: u64, seq: SeqView<'_>) -> Self {
        Self {
            donor,
            d1: seq.base(donor as usize),
            d2: seq.base((donor + 1) as usize),
            a1: seq.base((donor + intron_len - 2) as usize),
            a2: seq.base((donor + intron_len - 1) as usize),
        }
    }

    /// Move the window to `donor`, reusing what overlaps.
    ///
    /// A step of exactly one base either way shares two of the four: moving
    /// right, the old `d2`/`a2` become the new `d1`/`a1`; moving left, the old
    /// `d1`/`a1` become the new `d2`/`a2`. Any other step is rare enough that
    /// re-reading all four is the simpler answer.
    #[inline]
    fn slide_to(&mut self, donor: u64, intron_len: u64, seq: SeqView<'_>) {
        if donor == self.donor + 1 {
            self.d1 = self.d2;
            self.a1 = self.a2;
            self.d2 = seq.base((donor + 1) as usize);
            self.a2 = seq.base((donor + intron_len - 1) as usize);
            self.donor = donor;
        } else if donor + 1 == self.donor {
            self.d2 = self.d1;
            self.a2 = self.a1;
            self.d1 = seq.base(donor as usize);
            self.a1 = seq.base((donor + intron_len - 2) as usize);
            self.donor = donor;
        } else if donor != self.donor {
            *self = Self::at(donor, intron_len, seq);
        }
    }

    /// Base encoding: A=0, C=1, G=2, T=3.
    ///
    /// A table lookup rather than a match on the four bases. The match compiled
    /// to a chain of compares whose outcome is data-dependent and close to
    /// unpredictable, which the scan pays at every junction position; the table
    /// is one load at a computed index. `|` binds tighter than `>=` in Rust, so
    /// the guard tests the union of the four bases and rejects anything that is
    /// not `A`, `C`, `G` or `T` — `N`, the chromosome-boundary byte, and the
    /// out-of-range sentinel all take that path, exactly as no match arm
    /// covered them.
    #[inline]
    fn motif(&self) -> SpliceMotif {
        let (d1, d2, a1, a2) = (self.d1, self.d2, self.a1, self.a2);
        if d1 | d2 | a1 | a2 >= 4 {
            return SpliceMotif::NonCanonical;
        }
        MOTIF_TABLE[motif_index(d1, d2, a1, a2)]
    }
}

/// Every four-base combination, packed `d1 d2 a1 a2` at two bits each.
static MOTIF_TABLE: [SpliceMotif; 256] = build_motif_table();

/// Pack the four bases into the table index, two bits each.
const fn motif_index(d1: u8, d2: u8, a1: u8, a2: u8) -> usize {
    ((d1 << 6) | (d2 << 4) | (a1 << 2) | a2) as usize
}

const fn build_motif_table() -> [SpliceMotif; 256] {
    // A=0, C=1, G=2, T=3.
    let mut t = [SpliceMotif::NonCanonical; 256];
    t[motif_index(2, 3, 0, 2)] = SpliceMotif::GtAg;
    t[motif_index(2, 1, 0, 2)] = SpliceMotif::GcAg;
    t[motif_index(0, 3, 0, 1)] = SpliceMotif::AtAc;
    t[motif_index(1, 3, 0, 1)] = SpliceMotif::CtAc;
    t[motif_index(1, 3, 2, 1)] = SpliceMotif::CtGc;
    t[motif_index(2, 3, 0, 3)] = SpliceMotif::GtAt;
    t
}

/// Key for [`JunctionScanCache`]: the arguments of
/// `find_best_junction_position` that vary within a single window's stitch
/// recursion. Everything else (`read_seq`, the genome, `is_reverse`,
/// `n_genome`) is loop-invariant there, so these six fields identify a scan
/// exactly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct JunctionScanKey {
    r_a_end: usize,
    g_a_end: u64,
    r_gap: i64,
    g_gap: i64,
    prev_exon_len: usize,
    next_seed_len: usize,
}

/// Per-window memo table for the junction-position scan.
///
/// Create one per `stitch_seeds_core` call and pass it down the recursion; it
/// must not outlive the read, genome and strand it was filled for; the
/// per-window lifetime guarantees by construction.
#[derive(Default)]
pub struct JunctionScanCache {
    map: rustc_hash::FxHashMap<JunctionScanKey, (i32, SpliceMotif, i32, u32, u32)>,
}

impl JunctionScanCache {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Splice junction motif types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpliceMotif {
    /// GT-AG (canonical, + strand)
    GtAg,
    /// CT-AC (canonical, - strand; reverse complement of GT-AG)
    CtAc,
    /// GC-AG (semi-canonical, + strand)
    GcAg,
    /// CT-GC (semi-canonical, - strand; reverse complement of GC-AG)
    CtGc,
    /// AT-AC (semi-canonical, + strand)
    AtAc,
    /// GT-AT (semi-canonical, - strand; reverse complement of AT-AC)
    GtAt,
    /// Non-canonical
    NonCanonical,
}

impl SpliceMotif {
    /// Get the implied transcript strand from this motif.
    /// Forward-strand motifs (GT-AG, GC-AG, AT-AC) → Some('+')
    /// Reverse-strand motifs (CT-AC, CT-GC, GT-AT) → Some('-')
    /// Non-canonical → None (no strand information)
    pub fn implied_strand(&self) -> Option<char> {
        match self {
            SpliceMotif::GtAg | SpliceMotif::GcAg | SpliceMotif::AtAc => Some('+'),
            SpliceMotif::CtAc | SpliceMotif::CtGc | SpliceMotif::GtAt => Some('-'),
            SpliceMotif::NonCanonical => None,
        }
    }

    /// Get the motif filter category index for outSJfilter* parameters.
    /// 0 = non-canonical, 1 = GT/AG or CT/AC, 2 = GC/AG or CT/GC, 3 = AT/AC or GT/AT
    pub fn filter_category(&self) -> usize {
        match self {
            SpliceMotif::NonCanonical => 0,
            SpliceMotif::GtAg | SpliceMotif::CtAc => 1,
            SpliceMotif::GcAg | SpliceMotif::CtGc => 2,
            SpliceMotif::AtAc | SpliceMotif::GtAt => 3,
        }
    }

    /// Get the filter category from an encoded motif value (0-6 as in SJ.out.tab).
    /// 0→0 (non-canonical), 1|2→1 (GT/AG family), 3|4→2 (GC/AG family), 5|6→3 (AT/AC family)
    pub fn filter_category_from_encoded(encoded: u8) -> usize {
        match encoded {
            1 | 2 => 1,
            3 | 4 => 2,
            5 | 6 => 3,
            _ => 0, // 0 or any unknown = non-canonical
        }
    }
}

/// Gap type between aligned regions
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GapType {
    /// Insertion in read
    Insertion(u32),
    /// Deletion in read
    Deletion(u32),
    /// Splice junction
    SpliceJunction { intron_len: u32, motif: SpliceMotif },
}

#[cfg(test)]
mod tests {
    use crate::params::Parameters;

    use super::*;

    /// The table has to agree with the match it replaced on every input the
    /// scan can present, not merely on the six motifs. That includes `N` (4),
    /// the chromosome-boundary byte (5) and the out-of-range sentinel, all of
    /// which no match arm covered and which must therefore be `NonCanonical`.
    #[test]
    fn the_motif_table_agrees_with_the_original_match_on_every_input() {
        fn original(d1: u8, d2: u8, a1: u8, a2: u8) -> SpliceMotif {
            match (d1, d2, a1, a2) {
                (2, 3, 0, 2) => SpliceMotif::GtAg,
                (2, 1, 0, 2) => SpliceMotif::GcAg,
                (0, 3, 0, 1) => SpliceMotif::AtAc,
                (1, 3, 0, 1) => SpliceMotif::CtAc,
                (1, 3, 2, 1) => SpliceMotif::CtGc,
                (2, 3, 0, 3) => SpliceMotif::GtAt,
                _ => SpliceMotif::NonCanonical,
            }
        }

        let values = [0u8, 1, 2, 3, 4, 5, crate::genome::OUT_OF_RANGE];
        for &d1 in &values {
            for &d2 in &values {
                for &a1 in &values {
                    for &a2 in &values {
                        let w = MotifWindow {
                            donor: 0,
                            d1,
                            d2,
                            a1,
                            a2,
                        };
                        assert_eq!(
                            w.motif(),
                            original(d1, d2, a1, a2),
                            "disagreement at ({d1}, {d2}, {a1}, {a2})"
                        );
                    }
                }
            }
        }
    }

    fn make_test_genome(seq: &[u8]) -> Genome {
        // Create simple genome with one chromosome
        let n_genome = ((seq.len() as u64 + 1) / 64 + 1) * 64; // Pad to 64-byte boundary
        let mut sequence = vec![5u8; (n_genome * 2) as usize];

        // Copy forward sequence
        sequence[0..seq.len()].copy_from_slice(seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        Genome {
            transform_blocks: None,
            sequence: sequence.into(),
            n_genome,
            n_genome_real: n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![seq.len() as u64],
            chr_start: vec![0, n_genome],
        }
    }

    #[test]
    fn test_detect_gtag_motif() {
        // Sequence layout (0-based):
        // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
        // A  A [G  T  C  C  C  C  C  C  C  C  A  G] A  A
        //       ^donor                    ^acceptor
        // Intron from position 2 to 14 (exclusive), length = 12
        // Donor: positions 2,3 (GT)
        // Acceptor: positions 12,13 (AG)
        // Bases: A=0, C=1, G=2, T=3
        let seq = vec![
            0, 0, // AA (positions 0-1)
            2, 3, // GT (positions 2-3, donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's (positions 4-11, intron body)
            0, 2, // AG (positions 12-13, acceptor)
            0, 0, // AA (positions 14-15)
        ];
        let genome = make_test_genome(&seq);

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        // Intron from position 2, length 12 (spans positions 2-13 inclusive)
        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::GtAg);

        let score = scorer.score_splice_junction(motif);
        assert_eq!(score, 0); // Canonical
    }

    #[test]
    fn test_detect_gcag_motif() {
        // GC-AG motif: (2,1,0,2)
        let seq = vec![
            0, 0, // AA
            2, 1, // GC (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 2, // AG (acceptor)
            0, 0, // AA
        ];
        let genome = make_test_genome(&seq);

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::GcAg);

        let score = scorer.score_splice_junction(motif);
        assert_eq!(score, -4);
    }

    #[test]
    fn test_detect_atac_motif() {
        // AT-AC motif: (0,3,0,1)
        let seq = vec![
            0, 0, // AA
            0, 3, // AT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 1, // AC (acceptor)
            0, 0, // AA
        ];
        let genome = make_test_genome(&seq);

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::AtAc);

        let score = scorer.score_splice_junction(motif);
        assert_eq!(score, -8);
    }

    #[test]
    fn test_detect_noncanonical_motif() {
        // Some random motif
        let seq = vec![
            0, 0, 1, 1, // AA CC (donor)
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 10 C's
            3, 3, 0, 0, // TT AA (acceptor)
        ];
        let genome = make_test_genome(&seq);

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        let motif = scorer.detect_splice_motif(2, 12, &genome);
        assert_eq!(motif, SpliceMotif::NonCanonical);

        let score = scorer.score_splice_junction(motif);
        assert_eq!(score, -8);
    }

    #[test]
    fn test_score_gap_insertion() {
        let genome = make_test_genome(&[0, 1, 2, 3]);
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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        let (score, gap_type) = scorer.score_gap(0, 5, 0, &genome);
        assert_eq!(score, -2 + (-2 * 5)); // open + 5*base
        assert_eq!(gap_type, GapType::Insertion(5));
    }

    #[test]
    fn test_score_gap_deletion() {
        let genome = make_test_genome(&[0, 1, 2, 3]);
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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        // Small gap (< align_intron_min) is deletion
        let (score, gap_type) = scorer.score_gap(10, 0, 0, &genome);
        assert_eq!(score, -2 + (-2 * 10));
        assert_eq!(gap_type, GapType::Deletion(10));
    }

    #[test]
    fn test_score_gap_splice_junction() {
        // Create genome with GT-AG motif
        let seq = vec![
            0, 0, 2, 3, // AA GT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 22 C's
            0, 2, 0, 0, // AG AA (acceptor)
        ];
        let genome = make_test_genome(&seq);

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        // Gap starting at position 2 (GT), length 26 (>= 21) is splice junction
        let (score, gap_type) = scorer.score_gap(26, 0, 2, &genome);
        assert_eq!(score, 0); // Canonical GT-AG
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction {
                intron_len: 26,
                motif: SpliceMotif::GtAg
            }
        ));
    }

    #[test]
    fn test_annotated_junction_bonus() {
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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        let annotated_score = scorer.score_annotated_junction(0, true);
        assert_eq!(annotated_score, 2);

        let novel_score = scorer.score_annotated_junction(0, false);
        assert_eq!(novel_score, 0);

        let annotated_noncanon = scorer.score_annotated_junction(-8, true);
        assert_eq!(annotated_noncanon, 2);
    }

    #[test]
    fn test_detect_reverse_complement_motifs() {
        // Test all 3 reverse-complement motifs on the forward genome
        // These appear at minus-strand gene splice sites

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        // CT-AC motif: (1,3,0,1) — reverse complement of GT-AG
        let seq_ctac = vec![
            0, 0, // AA
            1, 3, // CT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 1, // AC (acceptor)
            0, 0, // AA
        ];
        let genome_ctac = make_test_genome(&seq_ctac);
        let motif = scorer.detect_splice_motif(2, 12, &genome_ctac);
        assert_eq!(motif, SpliceMotif::CtAc);
        // Should score same as canonical GT-AG
        assert_eq!(scorer.score_splice_junction(motif), 0);

        // CT-GC motif: (1,3,2,1) — reverse complement of GC-AG
        let seq_ctgc = vec![
            0, 0, // AA
            1, 3, // CT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            2, 1, // GC (acceptor)
            0, 0, // AA
        ];
        let genome_ctgc = make_test_genome(&seq_ctgc);
        let motif = scorer.detect_splice_motif(2, 12, &genome_ctgc);
        assert_eq!(motif, SpliceMotif::CtGc);
        assert_eq!(scorer.score_splice_junction(motif), -4);

        // GT-AT motif: (2,3,0,3) — reverse complement of AT-AC
        let seq_gtat = vec![
            0, 0, // AA
            2, 3, // GT (donor)
            1, 1, 1, 1, 1, 1, 1, 1, // 8 C's
            0, 3, // AT (acceptor)
            0, 0, // AA
        ];
        let genome_gtat = make_test_genome(&seq_gtat);
        let motif = scorer.detect_splice_motif(2, 12, &genome_gtat);
        assert_eq!(motif, SpliceMotif::GtAt);
        assert_eq!(scorer.score_splice_junction(motif), -8);
    }

    #[test]
    fn test_align_intron_max_default() {
        // When alignIntronMax=0 (default), scorer uses u32::MAX (no limit) — STAR-faithful.
        // STAR's stitchAlignToTranscript.cpp: `if (Del>alignIntronMax && alignIntronMax>0)`
        // meaning alignIntronMax=0 disables the check entirely.
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "reads.fq"]);
        assert_eq!(params.align_intron_max, 0);
        let scorer = AlignmentScorer::from_params(&params);
        assert_eq!(scorer.align_intron_max, u32::MAX);
    }

    #[test]
    fn test_align_intron_max_custom() {
        // Custom alignIntronMax should be passed through directly
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "reads.fq",
            "--alignIntronMax",
            "100000",
        ]);
        assert_eq!(params.align_intron_max, 100_000);
        let scorer = AlignmentScorer::from_params(&params);
        assert_eq!(scorer.align_intron_max, 100_000);
    }

    #[test]
    fn test_gap_at_intron_max_is_splice_junction() {
        // A gap exactly at alignIntronMax should still be a splice junction
        // Create genome large enough with GT-AG motif at boundaries
        let mut seq = vec![0u8; 600_000]; // ~600kb genome
        // Place GT at position 100
        seq[100] = 2; // G
        seq[101] = 3; // T
        // Place AG at position 100 + 589824 - 2, 100 + 589824 - 1
        let acceptor_pos = 100 + 589_824 - 2;
        if acceptor_pos + 1 < seq.len() {
            seq[acceptor_pos] = 0; // A
            seq[acceptor_pos + 1] = 2; // G
        }
        let genome = make_test_genome(&seq);

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
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        // Gap of exactly 589824 starting at position 100 should be splice junction
        let (score, gap_type) = scorer.score_gap(589_824, 0, 100, &genome);
        assert_eq!(score, 0); // GT-AG canonical
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction {
                intron_len: 589_824,
                ..
            }
        ));
    }

    #[test]
    fn test_implied_strand() {
        // Forward-strand motifs
        assert_eq!(SpliceMotif::GtAg.implied_strand(), Some('+'));
        assert_eq!(SpliceMotif::GcAg.implied_strand(), Some('+'));
        assert_eq!(SpliceMotif::AtAc.implied_strand(), Some('+'));
        // Reverse-strand motifs
        assert_eq!(SpliceMotif::CtAc.implied_strand(), Some('-'));
        assert_eq!(SpliceMotif::CtGc.implied_strand(), Some('-'));
        assert_eq!(SpliceMotif::GtAt.implied_strand(), Some('-'));
        // Non-canonical
        assert_eq!(SpliceMotif::NonCanonical.implied_strand(), None);
    }

    #[test]
    fn test_filter_category() {
        assert_eq!(SpliceMotif::NonCanonical.filter_category(), 0);
        assert_eq!(SpliceMotif::GtAg.filter_category(), 1);
        assert_eq!(SpliceMotif::CtAc.filter_category(), 1);
        assert_eq!(SpliceMotif::GcAg.filter_category(), 2);
        assert_eq!(SpliceMotif::CtGc.filter_category(), 2);
        assert_eq!(SpliceMotif::AtAc.filter_category(), 3);
        assert_eq!(SpliceMotif::GtAt.filter_category(), 3);
    }

    #[test]
    fn test_filter_category_from_encoded() {
        assert_eq!(SpliceMotif::filter_category_from_encoded(0), 0); // non-canonical
        assert_eq!(SpliceMotif::filter_category_from_encoded(1), 1); // GT/AG
        assert_eq!(SpliceMotif::filter_category_from_encoded(2), 1); // CT/AC
        assert_eq!(SpliceMotif::filter_category_from_encoded(3), 2); // GC/AG
        assert_eq!(SpliceMotif::filter_category_from_encoded(4), 2); // CT/GC
        assert_eq!(SpliceMotif::filter_category_from_encoded(5), 3); // AT/AC
        assert_eq!(SpliceMotif::filter_category_from_encoded(6), 3); // GT/AT
        assert_eq!(SpliceMotif::filter_category_from_encoded(7), 0); // unknown → non-canonical
    }

    #[test]
    fn test_gap_exceeding_intron_max_is_deletion() {
        // A gap exceeding alignIntronMax should be treated as deletion
        let seq = vec![0u8; 100]; // Small genome, gap length doesn't need real sequence
        let genome = make_test_genome(&seq);

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
            align_intron_max: 1000, // Small max for testing
            score_genomic_length_log2_scale: -0.25,
            score_stitch_sj_shift: 1,
            align_spliced_mate_map_lmin: 0,
            align_spliced_mate_map_lmin_over_lmate: 0.66,
            out_filter_score_min_over_lread: 0.66,
            align_ends_type: crate::params::AlignEndsType::default(),
        };

        // Gap of 1001 (> 1000 max) should be deletion, not splice junction
        let (_score, gap_type) = scorer.score_gap(1001, 0, 0, &genome);
        assert!(matches!(gap_type, GapType::Deletion(1001)));

        // Gap of 1000 (== max) should still be splice junction
        let (_score, gap_type) = scorer.score_gap(1000, 0, 0, &genome);
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction {
                intron_len: 1000,
                ..
            }
        ));

        // Gap of 21 (== min) should be splice junction
        let (_score, gap_type) = scorer.score_gap(21, 0, 0, &genome);
        assert!(matches!(
            gap_type,
            GapType::SpliceJunction { intron_len: 21, .. }
        ));

        // Gap of 20 (< min) should be deletion
        let (_score, gap_type) = scorer.score_gap(20, 0, 0, &genome);
        assert!(matches!(gap_type, GapType::Deletion(20)));
    }

    fn make_scorer_for_junction_test() -> AlignmentScorer {
        AlignmentScorer {
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
            align_ends_type: crate::params::AlignEndsType::default(),
        }
    }

    #[test]
    fn test_junction_scan_finds_canonical() {
        // Genome layout (forward):
        //   pos: 0  1  2  3  4  5  6  7  8  ...  28 29 30 31
        //        A  C  G  T  A  G  T  C  C  ...   C  A  G  A
        //                       ^GT                 ^AG
        // Two seeds: A ends at genome pos 5 (exclusive), B starts at genome pos 30 (exclusive end)
        // The gap: genome_gap = 30-5 = 25, read_gap = 1
        // del = 24 (>= 21 intron min)
        // Without jR shift (jR=0): intron starts at pos 5, motif at pos 5..6 = AG (non-canonical with G,T at 5,6)
        // Wait, let me set up the motif so jR=+1 shift reveals GT-AG:
        // At jR=0: donor at pos 5 = AG... not GT
        // At jR=+1: donor at pos 6 = GT, acceptor at pos 6+24-2 = 28, 6+24-1 = 29 = AG? Need to set that up.

        // Let me design carefully:
        // Intron length = del = genome_gap - read_gap
        // We want: 1 read gap base, donor at position (g_a_end + jr_shift)
        // Let g_a_end = 5 (exclusive), read_gap = 1
        // At jR=0: donor at pos 5, check motif at fwd pos 5
        // At jR=1: donor at pos 6, check motif at fwd pos 6
        // We want GT-AG at jR=1 position.

        // Genome: build so that:
        //   pos 6,7 = G,T (GT donor at pos 6)
        //   pos 6+del-2, 6+del-1 = A,G (AG acceptor)
        // With del = genome_gap - read_gap, genome_gap = 25, read_gap = 1, del = 24
        //   acceptor at pos 6+24-2=28, 6+24-1=29 → set pos 28=A(0), 29=G(2)

        let mut seq = vec![0u8; 100]; // A=0 fill
        // Seed A region (pos 0-4): matching bases
        seq[0] = 0; // A
        seq[1] = 1; // C
        seq[2] = 2; // G
        seq[3] = 3; // T
        seq[4] = 0; // A
        // At pos 5: not GT (make it A,C so jR=0 finds non-canonical)
        seq[5] = 0; // A
        // At pos 6,7: GT donor
        seq[6] = 2; // G
        seq[7] = 3; // T
        // Fill intron body
        seq[8..28].fill(1); // C
        // AG acceptor for jR=1 position (donor at pos 6, del=24, acceptor at 28,29)
        seq[28] = 0; // A
        seq[29] = 2; // G
        // Seed B region (pos 30+)
        seq[30] = 0; // A
        seq[31] = 1; // C

        let genome = make_test_genome(&seq);
        let scorer = make_scorer_for_junction_test();

        // Read: seed A covers read[0..5], gap has 1 base, seed B covers read[6..8]
        // The gap base at read[5] should match genome[5] (=A) on upstream side
        // and genome[5+24] = genome[29] (=G) on downstream side
        // Read[5] = A matches upstream but not downstream → score1 += 1 for jR=1
        let read_seq = vec![0, 1, 2, 3, 0, 0, 0, 1]; // ACGTAACG...

        let (jr_shift, motif, _motif_score, _, _) = scorer.find_best_junction_position(
            &read_seq,
            5,  // r_a_end (exclusive)
            5,  // g_a_end (exclusive)
            1,  // read_gap
            25, // genome_gap
            &genome,
            false, // forward strand
            genome.n_genome,
            5, // prev_exon_len
            3, // next_seed_len
        );

        // jR=1 should find GT-AG motif
        assert_eq!(motif, SpliceMotif::GtAg);
        assert_eq!(jr_shift, 1);
    }

    #[test]
    fn test_junction_scan_no_shift_needed() {
        // GT-AG motif is already at the natural junction boundary (jR=0)
        // Genome: donor GT at pos 5,6 and acceptor AG at pos 5+24-2=27, 5+24-1=28
        let mut seq = vec![1u8; 100]; // C fill
        seq[0] = 0;
        seq[1] = 1;
        seq[2] = 2;
        seq[3] = 3;
        seq[4] = 0;
        // GT donor at pos 5,6
        seq[5] = 2; // G
        seq[6] = 3; // T
        // AG acceptor at pos 27,28 (donor at 5, del=24: 5+24-2=27, 5+24-1=28)
        seq[27] = 0; // A
        seq[28] = 2; // G
        // Seed B
        seq[29] = 0;
        seq[30] = 1;

        let genome = make_test_genome(&seq);
        let scorer = make_scorer_for_junction_test();

        // Read gap = 0, genome gap = 24 (pure splice junction)
        let read_seq = vec![0, 1, 2, 3, 0, 0, 1]; // ACGTAAC

        let (jr_shift, motif, _motif_score, _, _) = scorer.find_best_junction_position(
            &read_seq,
            5,  // r_a_end (exclusive)
            5,  // g_a_end (exclusive)
            0,  // read_gap = 0 (pure splice)
            24, // genome_gap = intron length
            &genome,
            false,
            genome.n_genome,
            5,
            2, // next_seed_len
        );

        assert_eq!(motif, SpliceMotif::GtAg);
        assert_eq!(jr_shift, 0);
    }

    #[test]
    fn test_junction_scan_left_flush_noncanonical() {
        // Non-canonical junction in a repeat region should be flushed left
        // Build genome where no canonical motif exists, and there's a repeat
        // at the junction boundary
        let mut seq = vec![1u8; 100]; // C fill
        // Seed A region
        seq[0] = 0;
        seq[1] = 1;
        seq[2] = 2;
        seq[3] = 3;
        seq[4] = 0;
        // At pos 5: non-canonical (CC at donor, CC at acceptor)
        seq[5] = 1; // C
        seq[6] = 1; // C
        // Make a repeat: genome[5] == genome[5+24] (both C)
        // del = 24, acceptor at 5+24-2=27, 5+24-1=28
        seq[27] = 1; // C (same as donor d1 at jR=0)
        seq[28] = 1; // C (same as donor d2 at jR=0)
        // Also make genome[4] == genome[4+24] = genome[28] (both 0=A)
        // Already seq[4]=0, seq[28]=1. Not equal, so repeat only extends 0 left.
        // Make them equal: set seq[28]=0
        // Actually, let me set up a 2-base repeat:
        seq[4] = 1; // C at pos 4 (end of seed A)
        seq[28] = 1; // C at pos 28 = pos 4+del
        seq[5] = 1; // C at pos 5
        seq[29] = 1; // C at pos 29 = pos 5+del

        let genome = make_test_genome(&seq);
        let scorer = make_scorer_for_junction_test();

        // read_gap = 0
        let read_seq = vec![0, 1, 2, 3, 1, 1, 1]; // the gap base

        let (jr_shift, motif, _, _, _) = scorer.find_best_junction_position(
            &read_seq,
            5,  // r_a_end
            5,  // g_a_end
            0,  // read_gap
            24, // genome_gap
            &genome,
            false,
            genome.n_genome,
            5,
            2, // next_seed_len
        );

        // Non-canonical junction should be flushed left
        assert_eq!(motif, SpliceMotif::NonCanonical);
        // jr_shift should be <= 0 (flushed left)
        assert!(
            jr_shift <= 0,
            "Non-canonical should flush left, got jr_shift={jr_shift}"
        );
    }
}
