/// Alignment statistics tracking and reporting
use log::info;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::align::transcript::{CigarOp, Transcript};
use crate::junction::encode_motif;

/// Reason a read could not be mapped
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnmappedReason {
    /// No seeds or no clusters found
    Other,
    /// Transcripts generated but all filtered out (score, match count, length)
    TooShort,
    /// All transcripts filtered by mismatch count/rate
    TooManyMismatches,
    /// Too many multi-mapping loci (nTr > outFilterMultimapNmax)
    TooManyLoci,
}

/// Tracks alignment statistics for a read mapping run
/// Thread-safe using atomic counters
#[derive(Debug)]
pub struct AlignmentStats {
    /// Total number of reads processed
    pub total_reads: AtomicU64,
    /// Reads that mapped uniquely (exactly 1 locus)
    pub uniquely_mapped: AtomicU64,
    /// Reads that mapped to multiple loci (2-N loci)
    pub multi_mapped: AtomicU64,
    /// Reads that did not map
    pub unmapped: AtomicU64,
    /// Reads that mapped to too many loci (exceeds outFilterMultimapNmax)
    pub too_many_loci: AtomicU64,

    // --- Log.final.out fields ---
    /// Sum of all input read lengths (for average input read length)
    pub read_bases: AtomicU64,

    /// Sum of Match/Equal/Diff lengths from unique mappers (exon-aligned bases)
    pub mapped_bases: AtomicU64,
    /// Sum of n_mismatch from unique mappers
    pub mapped_mismatches: AtomicU64,
    /// Count of Ins CIGAR ops from unique mappers
    pub mapped_ins_count: AtomicU64,
    /// Sum of Ins op lengths from unique mappers
    pub mapped_ins_bases: AtomicU64,
    /// Count of Del CIGAR ops from unique mappers
    pub mapped_del_count: AtomicU64,
    /// Sum of Del op lengths from unique mappers
    pub mapped_del_bases: AtomicU64,

    /// Per-motif splice event counts (unique mappers only).
    /// Index matches encode_motif(): [0]=NonCanonical, [1]=GT/AG, [2]=CT/AC,
    /// [3]=GC/AG, [4]=CT/GC, [5]=AT/AC, [6]=GT/AT
    pub splices_by_motif: [AtomicU64; 7],
    /// Splice events at annotated junctions (unique mappers only)
    pub splices_annotated: AtomicU64,

    /// Unmapped reads: too many mismatches
    pub unmapped_mismatches: AtomicU64,
    /// Unmapped reads: too short / low score / filtered
    pub unmapped_short: AtomicU64,
    /// Unmapped reads: no seeds, no clusters, other
    pub unmapped_other: AtomicU64,

    /// Number of chimeric reads
    pub chimeric_reads: AtomicU64,

    /// Number of half-mapped pairs (one mate mapped, other unmapped/rescue failed)
    pub half_mapped_pairs: AtomicU64,
}

impl Default for AlignmentStats {
    fn default() -> Self {
        Self {
            total_reads: AtomicU64::new(0),
            uniquely_mapped: AtomicU64::new(0),
            multi_mapped: AtomicU64::new(0),
            unmapped: AtomicU64::new(0),
            too_many_loci: AtomicU64::new(0),
            read_bases: AtomicU64::new(0),
            mapped_bases: AtomicU64::new(0),
            mapped_mismatches: AtomicU64::new(0),
            mapped_ins_count: AtomicU64::new(0),
            mapped_ins_bases: AtomicU64::new(0),
            mapped_del_count: AtomicU64::new(0),
            mapped_del_bases: AtomicU64::new(0),
            splices_by_motif: [
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
            ],
            splices_annotated: AtomicU64::new(0),
            unmapped_mismatches: AtomicU64::new(0),
            unmapped_short: AtomicU64::new(0),
            unmapped_other: AtomicU64::new(0),
            chimeric_reads: AtomicU64::new(0),
            half_mapped_pairs: AtomicU64::new(0),
        }
    }
}

impl AlignmentStats {
    /// Create new statistics tracker
    pub fn new() -> Self {
        Self::default()
    }

    /// Record an alignment result (thread-safe)
    ///
    /// # Arguments
    /// * `n_alignments` - Number of valid alignments found for the read
    /// * `max_multimaps` - Maximum allowed multi-map count (outFilterMultimapNmax)
    pub fn record_alignment(&self, n_alignments: usize, max_multimaps: usize) {
        self.total_reads.fetch_add(1, Ordering::Relaxed);
        match n_alignments {
            0 => {
                self.unmapped.fetch_add(1, Ordering::Relaxed);
            }
            1 => {
                self.uniquely_mapped.fetch_add(1, Ordering::Relaxed);
            }
            n if n <= max_multimaps => {
                self.multi_mapped.fetch_add(1, Ordering::Relaxed);
            }
            _ => {
                self.too_many_loci.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Record the length of an input read (for average input read length)
    pub fn record_read_bases(&self, len: u64) {
        self.read_bases.fetch_add(len, Ordering::Relaxed);
    }

    /// Record detailed stats from a uniquely mapped transcript.
    /// Walks the CIGAR to extract mapped bases, ins/del counts, and splice motif counts.
    /// Only call this for unique mappers (n_alignments == 1).
    pub fn record_transcript_stats(&self, transcript: &Transcript) {
        // Walk CIGAR for mapped bases, ins/del
        for op in &transcript.cigar {
            match op {
                CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                    self.mapped_bases.fetch_add(*len as u64, Ordering::Relaxed);
                }
                CigarOp::Ins(len) => {
                    self.mapped_ins_count.fetch_add(1, Ordering::Relaxed);
                    self.mapped_ins_bases
                        .fetch_add(*len as u64, Ordering::Relaxed);
                }
                CigarOp::Del(len) => {
                    self.mapped_del_count.fetch_add(1, Ordering::Relaxed);
                    self.mapped_del_bases
                        .fetch_add(*len as u64, Ordering::Relaxed);
                }
                CigarOp::RefSkip(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
            }
        }

        // Mismatches
        self.mapped_mismatches
            .fetch_add(transcript.n_mismatch as u64, Ordering::Relaxed);

        // Splice motif counts
        for motif in &transcript.junction_motifs {
            let idx = encode_motif(*motif) as usize;
            self.splices_by_motif[idx].fetch_add(1, Ordering::Relaxed);
        }

        // Annotated junction counts
        for annotated in &transcript.junction_annotated {
            if *annotated {
                self.splices_annotated.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Record the reason a read was unmapped
    pub fn record_unmapped_reason(&self, reason: UnmappedReason) {
        match reason {
            UnmappedReason::TooManyMismatches => {
                self.unmapped_mismatches.fetch_add(1, Ordering::Relaxed);
            }
            UnmappedReason::TooShort => {
                self.unmapped_short.fetch_add(1, Ordering::Relaxed);
            }
            UnmappedReason::TooManyLoci => {
                // Tracked via record_alignment's too_many_loci path, not here
            }
            UnmappedReason::Other => {
                self.unmapped_other.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Record a chimeric read
    pub fn record_chimeric(&self) {
        self.chimeric_reads.fetch_add(1, Ordering::Relaxed);
    }

    /// Record a half-mapped pair (one mate mapped, other unmapped)
    pub fn record_half_mapped(&self) {
        self.half_mapped_pairs.fetch_add(1, Ordering::Relaxed);
    }

    /// Print summary statistics to log
    pub fn print_summary(&self) {
        // Load all atomics once at start
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        let unmapped = self.unmapped.load(Ordering::Relaxed);
        let too_many_loci = self.too_many_loci.load(Ordering::Relaxed);

        if total_reads == 0 {
            info!("No reads processed");
            return;
        }

        info!("=== Alignment Summary ===");
        info!("Number of input reads: {}", total_reads);
        info!(
            "Uniquely mapped reads: {} ({:.2}%)",
            uniquely_mapped,
            100.0 * uniquely_mapped as f64 / total_reads as f64
        );
        info!(
            "Multi-mapped reads: {} ({:.2}%)",
            multi_mapped,
            100.0 * multi_mapped as f64 / total_reads as f64
        );
        info!(
            "Unmapped reads: {} ({:.2}%)",
            unmapped,
            100.0 * unmapped as f64 / total_reads as f64
        );
        let half_mapped = self.half_mapped_pairs.load(Ordering::Relaxed);
        if half_mapped > 0 {
            info!(
                "Half-mapped pairs: {} ({:.2}%)",
                half_mapped,
                100.0 * half_mapped as f64 / total_reads as f64
            );
        }
        if too_many_loci > 0 {
            info!(
                "Reads with too many loci: {} ({:.2}%)",
                too_many_loci,
                100.0 * too_many_loci as f64 / total_reads as f64
            );
        }

        // Calculate mapped percentage
        let mapped = uniquely_mapped + multi_mapped;
        info!(
            "Total mapped: {} ({:.2}%)",
            mapped,
            100.0 * mapped as f64 / total_reads as f64
        );
    }

    /// Get percentage of uniquely mapped reads
    pub fn unique_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * uniquely_mapped as f64 / total_reads as f64
        }
    }

    /// Get percentage of multi-mapped reads
    pub fn multi_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * multi_mapped as f64 / total_reads as f64
        }
    }

    /// Get percentage of unmapped reads
    pub fn unmapped_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let unmapped = self.unmapped.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * unmapped as f64 / total_reads as f64
        }
    }

    /// Get total mapped reads (unique + multi)
    pub fn total_mapped(&self) -> u64 {
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        uniquely_mapped + multi_mapped
    }

    /// Get percentage of mapped reads
    pub fn mapped_percent(&self) -> f64 {
        let total_reads = self.total_reads.load(Ordering::Relaxed);
        if total_reads == 0 {
            0.0
        } else {
            100.0 * self.total_mapped() as f64 / total_reads as f64
        }
    }

    /// Get total number of reads processed
    pub fn total_reads(&self) -> u64 {
        self.total_reads.load(Ordering::Relaxed)
    }

    /// Undo a mapped read record for BySJout filtering.
    /// Moves one read from uniquely_mapped (or multi_mapped) to unmapped.
    /// Since we don't know which category the read was in, we try unique first
    /// (most reads are uniquely mapped), then multi.
    pub fn undo_mapped_record_bysj(&self) {
        // Try to decrement uniquely_mapped first
        let mut current = self.uniquely_mapped.load(Ordering::Relaxed);
        loop {
            if current == 0 {
                break;
            }
            match self.uniquely_mapped.compare_exchange_weak(
                current,
                current - 1,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    self.unmapped.fetch_add(1, Ordering::Relaxed);
                    return;
                }
                Err(x) => current = x,
            }
        }

        // If no unique reads, try multi_mapped
        let mut current = self.multi_mapped.load(Ordering::Relaxed);
        loop {
            if current == 0 {
                break;
            }
            match self.multi_mapped.compare_exchange_weak(
                current,
                current - 1,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    self.unmapped.fetch_add(1, Ordering::Relaxed);
                    return;
                }
                Err(x) => current = x,
            }
        }
    }

    /// Write STAR-compatible Log.final.out file
    pub fn write_log_final(
        &self,
        path: &Path,
        time_start: chrono::DateTime<chrono::Local>,
        time_map_start: chrono::DateTime<chrono::Local>,
        time_finish: chrono::DateTime<chrono::Local>,
    ) -> std::io::Result<()> {
        use std::io::Write;

        let total_reads = self.total_reads.load(Ordering::Relaxed);
        let uniquely_mapped = self.uniquely_mapped.load(Ordering::Relaxed);
        let multi_mapped = self.multi_mapped.load(Ordering::Relaxed);
        let too_many_loci = self.too_many_loci.load(Ordering::Relaxed);
        let read_bases = self.read_bases.load(Ordering::Relaxed);
        let mapped_bases = self.mapped_bases.load(Ordering::Relaxed);
        let mapped_mismatches = self.mapped_mismatches.load(Ordering::Relaxed);
        let mapped_ins_count = self.mapped_ins_count.load(Ordering::Relaxed);
        let mapped_ins_bases = self.mapped_ins_bases.load(Ordering::Relaxed);
        let mapped_del_count = self.mapped_del_count.load(Ordering::Relaxed);
        let mapped_del_bases = self.mapped_del_bases.load(Ordering::Relaxed);
        let unmapped_mismatches = self.unmapped_mismatches.load(Ordering::Relaxed);
        let unmapped_short = self.unmapped_short.load(Ordering::Relaxed);
        let unmapped_other = self.unmapped_other.load(Ordering::Relaxed);
        let chimeric_reads = self.chimeric_reads.load(Ordering::Relaxed);

        // Splice counts
        let splices: Vec<u64> = (0..7)
            .map(|i| self.splices_by_motif[i].load(Ordering::Relaxed))
            .collect();
        let splices_annotated = self.splices_annotated.load(Ordering::Relaxed);
        let splices_total: u64 = splices.iter().sum();
        // GT/AG = motif[1] + motif[2], GC/AG = motif[3] + motif[4], AT/AC = motif[5] + motif[6]
        let splices_gtag = splices[1] + splices[2];
        let splices_gcag = splices[3] + splices[4];
        let splices_atac = splices[5] + splices[6];
        let splices_noncanonical = splices[0];

        // Computed fields
        let avg_input_read_length = read_bases.checked_div(total_reads).unwrap_or(0);

        let avg_mapped_length = if uniquely_mapped > 0 {
            mapped_bases as f64 / uniquely_mapped as f64
        } else {
            0.0
        };

        let elapsed_hours = {
            let elapsed = time_finish - time_map_start;
            elapsed.num_milliseconds() as f64 / 3_600_000.0
        };
        let mapping_speed = if elapsed_hours > 0.0 {
            total_reads as f64 / elapsed_hours / 1_000_000.0
        } else {
            0.0
        };

        let mismatch_rate = if mapped_bases > 0 {
            mapped_mismatches as f64 / mapped_bases as f64 * 100.0
        } else {
            0.0
        };
        let del_rate = if mapped_bases > 0 {
            mapped_del_bases as f64 / mapped_bases as f64 * 100.0
        } else {
            0.0
        };
        let del_avg_len = if mapped_del_count > 0 {
            mapped_del_bases as f64 / mapped_del_count as f64
        } else {
            0.0
        };
        let ins_rate = if mapped_bases > 0 {
            mapped_ins_bases as f64 / mapped_bases as f64 * 100.0
        } else {
            0.0
        };
        let ins_avg_len = if mapped_ins_count > 0 {
            mapped_ins_bases as f64 / mapped_ins_count as f64
        } else {
            0.0
        };

        // Percentages
        let pct = |n: u64| -> f64 {
            if total_reads > 0 {
                n as f64 / total_reads as f64 * 100.0
            } else {
                0.0
            }
        };

        let time_fmt = "%b %d %H:%M:%S";

        let mut f = std::fs::File::create(path)?;

        // Timestamps and speed
        writeln!(
            f,
            "{:>47} |\t{}",
            "Started job on",
            time_start.format(time_fmt)
        )?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Started mapping on",
            time_map_start.format(time_fmt)
        )?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Finished on",
            time_finish.format(time_fmt)
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}",
            "Mapping speed, Million of reads per hour", mapping_speed
        )?;

        // Input stats
        writeln!(f)?;
        writeln!(f, "{:>47} |\t{}", "Number of input reads", total_reads)?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Average input read length", avg_input_read_length
        )?;

        // UNIQUE READS
        writeln!(f, "{:>49}", "UNIQUE READS:")?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Uniquely mapped reads number", uniquely_mapped
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "Uniquely mapped reads %",
            pct(uniquely_mapped)
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}",
            "Average mapped length", avg_mapped_length
        )?;
        writeln!(f, "{:>47} |\t{}", "Number of splices: Total", splices_total)?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of splices: Annotated (sjdb)", splices_annotated
        )?;
        writeln!(f, "{:>47} |\t{}", "Number of splices: GT/AG", splices_gtag)?;
        writeln!(f, "{:>47} |\t{}", "Number of splices: GC/AG", splices_gcag)?;
        writeln!(f, "{:>47} |\t{}", "Number of splices: AT/AC", splices_atac)?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of splices: Non-canonical", splices_noncanonical
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "Mismatch rate per base, %", mismatch_rate
        )?;
        writeln!(f, "{:>47} |\t{:.2}%", "Deletion rate per base", del_rate)?;
        writeln!(f, "{:>47} |\t{:.2}", "Deletion average length", del_avg_len)?;
        writeln!(f, "{:>47} |\t{:.2}%", "Insertion rate per base", ins_rate)?;
        writeln!(
            f,
            "{:>47} |\t{:.2}",
            "Insertion average length", ins_avg_len
        )?;

        // MULTI-MAPPING READS
        writeln!(f, "{:>49}", "MULTI-MAPPING READS:")?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of reads mapped to multiple loci", multi_mapped
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "% of reads mapped to multiple loci",
            pct(multi_mapped)
        )?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of reads mapped to too many loci", too_many_loci
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "% of reads mapped to too many loci",
            pct(too_many_loci)
        )?;

        // UNMAPPED READS
        writeln!(f, "{:>49}", "UNMAPPED READS:")?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of reads unmapped: too many mismatches", unmapped_mismatches
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "% of reads unmapped: too many mismatches",
            pct(unmapped_mismatches)
        )?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of reads unmapped: too short", unmapped_short
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "% of reads unmapped: too short",
            pct(unmapped_short)
        )?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of reads unmapped: other", unmapped_other
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "% of reads unmapped: other",
            pct(unmapped_other)
        )?;

        // CHIMERIC READS
        writeln!(f, "{:>49}", "CHIMERIC READS:")?;
        writeln!(
            f,
            "{:>47} |\t{}",
            "Number of chimeric reads", chimeric_reads
        )?;
        writeln!(
            f,
            "{:>47} |\t{:.2}%",
            "% of chimeric reads",
            pct(chimeric_reads)
        )?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_default() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 0);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.too_many_loci.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_unique() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_unmapped() {
        let stats = AlignmentStats::new();
        stats.record_alignment(0, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_record_multi() {
        let stats = AlignmentStats::new();
        stats.record_alignment(5, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_too_many() {
        let stats = AlignmentStats::new();
        stats.record_alignment(15, 10);
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.too_many_loci.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_multiple_reads() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(0, 10); // unmapped
        stats.record_alignment(5, 10); // multi
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(15, 10); // too many

        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 5);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 2);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.too_many_loci.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_percentages() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(0, 10); // unmapped
        stats.record_alignment(5, 10); // multi
        stats.record_alignment(1, 10); // unique

        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 4);
        assert!((stats.unique_percent() - 50.0).abs() < 0.01);
        assert!((stats.multi_percent() - 25.0).abs() < 0.01);
        assert!((stats.unmapped_percent() - 25.0).abs() < 0.01);
        assert!((stats.mapped_percent() - 75.0).abs() < 0.01);
    }

    #[test]
    fn test_empty_stats() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.unique_percent(), 0.0);
        assert_eq!(stats.multi_percent(), 0.0);
        assert_eq!(stats.unmapped_percent(), 0.0);
        assert_eq!(stats.mapped_percent(), 0.0);
        assert_eq!(stats.total_mapped(), 0);
    }

    #[test]
    fn test_undo_mapped_record_bysj_unique() {
        let stats = AlignmentStats::new();
        stats.record_alignment(1, 10); // unique
        stats.record_alignment(1, 10); // unique

        stats.undo_mapped_record_bysj();

        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 2);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_undo_mapped_record_bysj_multi() {
        let stats = AlignmentStats::new();
        stats.record_alignment(5, 10); // multi

        stats.undo_mapped_record_bysj();

        // No unique reads, so multi should be decremented
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_undo_mapped_record_bysj_noop_when_empty() {
        let stats = AlignmentStats::new();
        stats.record_alignment(0, 10); // unmapped

        stats.undo_mapped_record_bysj();

        // Should be a no-op (no mapped reads to undo)
        assert_eq!(stats.unmapped.load(Ordering::Relaxed), 1);
        assert_eq!(stats.uniquely_mapped.load(Ordering::Relaxed), 0);
        assert_eq!(stats.multi_mapped.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_record_transcript_stats() {
        use crate::align::score::SpliceMotif;
        use crate::align::transcript::Exon;

        let stats = AlignmentStats::new();

        // Build a transcript with known CIGAR: 5S 45M 3I 2M 100N 50M 2D 5M
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1204,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1204,
                read_start: 0,
                read_end: 110,
                i_frag: 0,
            }],
            cigar: vec![
                CigarOp::SoftClip(5),
                CigarOp::Match(45),
                CigarOp::Ins(3),
                CigarOp::Match(2),
                CigarOp::RefSkip(100),
                CigarOp::Match(50),
                CigarOp::Del(2),
                CigarOp::Match(5),
            ],
            score: 100,
            n_mismatch: 3,
            n_gap: 1,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![true],
            read_seq: vec![0; 110],
        };

        stats.record_transcript_stats(&transcript);

        // mapped_bases = 45 + 2 + 50 + 5 = 102
        assert_eq!(stats.mapped_bases.load(Ordering::Relaxed), 102);
        // mismatches = 3
        assert_eq!(stats.mapped_mismatches.load(Ordering::Relaxed), 3);
        // 1 insertion of 3 bases
        assert_eq!(stats.mapped_ins_count.load(Ordering::Relaxed), 1);
        assert_eq!(stats.mapped_ins_bases.load(Ordering::Relaxed), 3);
        // 1 deletion of 2 bases
        assert_eq!(stats.mapped_del_count.load(Ordering::Relaxed), 1);
        assert_eq!(stats.mapped_del_bases.load(Ordering::Relaxed), 2);
        // 1 GT/AG splice (encode_motif(GtAg) = 1)
        assert_eq!(stats.splices_by_motif[1].load(Ordering::Relaxed), 1);
        // 1 annotated junction
        assert_eq!(stats.splices_annotated.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_record_unmapped_reason() {
        let stats = AlignmentStats::new();

        stats.record_unmapped_reason(UnmappedReason::Other);
        stats.record_unmapped_reason(UnmappedReason::TooShort);
        stats.record_unmapped_reason(UnmappedReason::TooShort);
        stats.record_unmapped_reason(UnmappedReason::TooManyMismatches);

        assert_eq!(stats.unmapped_other.load(Ordering::Relaxed), 1);
        assert_eq!(stats.unmapped_short.load(Ordering::Relaxed), 2);
        assert_eq!(stats.unmapped_mismatches.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_splice_motif_aggregation() {
        use crate::align::score::SpliceMotif;
        use crate::align::transcript::Exon;

        let stats = AlignmentStats::new();

        // Transcript with multiple junction types
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 2000,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 2000,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 4,
            junction_motifs: vec![
                SpliceMotif::GtAg,         // motif[1]
                SpliceMotif::CtAc,         // motif[2]
                SpliceMotif::GcAg,         // motif[3]
                SpliceMotif::NonCanonical, // motif[0]
            ],
            junction_annotated: vec![true, false, true, false],
            read_seq: vec![0; 100],
        };

        stats.record_transcript_stats(&transcript);

        // GT/AG = motif[1] + motif[2] = 1 + 1 = 2
        let gtag = stats.splices_by_motif[1].load(Ordering::Relaxed)
            + stats.splices_by_motif[2].load(Ordering::Relaxed);
        assert_eq!(gtag, 2);

        // GC/AG = motif[3] + motif[4] = 1 + 0 = 1
        let gcag = stats.splices_by_motif[3].load(Ordering::Relaxed)
            + stats.splices_by_motif[4].load(Ordering::Relaxed);
        assert_eq!(gcag, 1);

        // Non-canonical = motif[0] = 1
        assert_eq!(stats.splices_by_motif[0].load(Ordering::Relaxed), 1);

        // Annotated = 2 (GtAg + GcAg)
        assert_eq!(stats.splices_annotated.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_write_log_final_format() {
        use chrono::TimeZone;

        let stats = AlignmentStats::new();

        // Set up some stats
        stats.total_reads.store(10000, Ordering::Relaxed);
        stats.uniquely_mapped.store(8300, Ordering::Relaxed);
        stats.multi_mapped.store(500, Ordering::Relaxed);
        stats.too_many_loci.store(50, Ordering::Relaxed);
        stats.unmapped.store(1150, Ordering::Relaxed);
        stats.read_bases.store(1500000, Ordering::Relaxed);
        stats.mapped_bases.store(1200000, Ordering::Relaxed);
        stats.mapped_mismatches.store(4800, Ordering::Relaxed);
        stats.mapped_ins_count.store(100, Ordering::Relaxed);
        stats.mapped_ins_bases.store(150, Ordering::Relaxed);
        stats.mapped_del_count.store(80, Ordering::Relaxed);
        stats.mapped_del_bases.store(160, Ordering::Relaxed);
        stats.splices_by_motif[1].store(200, Ordering::Relaxed); // GtAg
        stats.splices_by_motif[2].store(50, Ordering::Relaxed); // CtAc
        stats.splices_annotated.store(180, Ordering::Relaxed);
        stats.unmapped_short.store(1000, Ordering::Relaxed);
        stats.unmapped_other.store(150, Ordering::Relaxed);
        stats.chimeric_reads.store(5, Ordering::Relaxed);

        let time_start = chrono::Local
            .with_ymd_and_hms(2025, 2, 10, 17, 10, 0)
            .unwrap();
        let time_map_start = chrono::Local
            .with_ymd_and_hms(2025, 2, 10, 17, 11, 0)
            .unwrap();
        let time_finish = chrono::Local
            .with_ymd_and_hms(2025, 2, 10, 17, 12, 0)
            .unwrap();

        let dir = tempfile::tempdir().unwrap();
        let log_path = dir.path().join("Log.final.out");

        stats
            .write_log_final(&log_path, time_start, time_map_start, time_finish)
            .unwrap();

        let content = std::fs::read_to_string(&log_path).unwrap();

        // Check field names are right-justified in 50 chars with " |\t" separator
        assert!(content.contains("Started job on |\t"));
        assert!(content.contains("Started mapping on |\t"));
        assert!(content.contains("Finished on |\t"));
        assert!(content.contains("Number of input reads |\t10000"));
        assert!(content.contains("Average input read length |\t150"));
        assert!(content.contains("UNIQUE READS:"));
        assert!(content.contains("Uniquely mapped reads number |\t8300"));
        assert!(content.contains("Uniquely mapped reads % |\t83.00%"));
        assert!(content.contains("Number of splices: Total |\t250"));
        assert!(content.contains("Number of splices: Annotated (sjdb) |\t180"));
        assert!(content.contains("Number of splices: GT/AG |\t250"));
        assert!(content.contains("MULTI-MAPPING READS:"));
        assert!(content.contains("Number of reads mapped to multiple loci |\t500"));
        assert!(content.contains("UNMAPPED READS:"));
        assert!(content.contains("Number of reads unmapped: too short |\t1000"));
        assert!(content.contains("CHIMERIC READS:"));
        assert!(content.contains("Number of chimeric reads |\t5"));
    }

    #[test]
    fn test_log_final_multiqc_fields() {
        use chrono::TimeZone;

        let stats = AlignmentStats::new();
        stats.total_reads.store(100, Ordering::Relaxed);
        stats.uniquely_mapped.store(80, Ordering::Relaxed);
        stats.read_bases.store(15000, Ordering::Relaxed);
        stats.mapped_bases.store(12000, Ordering::Relaxed);

        let t = chrono::Local.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let log_path = dir.path().join("Log.final.out");
        stats.write_log_final(&log_path, t, t, t).unwrap();

        let content = std::fs::read_to_string(&log_path).unwrap();

        // All 23 MultiQC-required field names (from MultiQC STAR module)
        let required_fields = [
            "Number of input reads",
            "Average input read length",
            "Uniquely mapped reads number",
            "Uniquely mapped reads %",
            "Average mapped length",
            "Number of splices: Total",
            "Number of splices: Annotated (sjdb)",
            "Number of splices: GT/AG",
            "Number of splices: GC/AG",
            "Number of splices: AT/AC",
            "Number of splices: Non-canonical",
            "Mismatch rate per base, %",
            "Deletion rate per base",
            "Deletion average length",
            "Insertion rate per base",
            "Insertion average length",
            "Number of reads mapped to multiple loci",
            "% of reads mapped to multiple loci",
            "Number of reads mapped to too many loci",
            "% of reads mapped to too many loci",
            "Number of reads unmapped: too many mismatches",
            "Number of reads unmapped: too short",
            "Number of reads unmapped: other",
        ];

        for field in &required_fields {
            assert!(
                content.contains(field),
                "Missing MultiQC field: '{}'",
                field
            );
        }
    }

    #[test]
    fn test_record_chimeric() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.chimeric_reads.load(Ordering::Relaxed), 0);
        stats.record_chimeric();
        stats.record_chimeric();
        assert_eq!(stats.chimeric_reads.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_half_mapped_counter() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.half_mapped_pairs.load(Ordering::Relaxed), 0);
        stats.record_half_mapped();
        assert_eq!(stats.half_mapped_pairs.load(Ordering::Relaxed), 1);
        stats.record_half_mapped();
        stats.record_half_mapped();
        assert_eq!(stats.half_mapped_pairs.load(Ordering::Relaxed), 3);
    }
}
