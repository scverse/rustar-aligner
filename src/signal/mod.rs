//! Coverage signal / `Signal.*.out.bg` bedGraph output (STAR `signalFromBAM.cpp`,
//! `--outWigType bedGraph`, stranded).
//!
//! Ported from STAR-rs `crates/star-align/src/star_wig.rs`. Per strand (str1 =
//! forward, str2 = reverse) two tracks accumulate: `Unique` (only `NH == 1`
//! reads, +1 per covered base) and `UniqueMultiple` (every reported alignment,
//! `+1/NH`). bedGraph intervals are emitted wherever the coverage changes.
//!
//! Only `--outWigType bedGraph` (stranded, full-length) is implemented, matching
//! STAR-rs's own align-time scope: `wiggle`, `--outWigStrand Unstranded`, and the
//! `--outWigType` 2nd word (`read1_5p`/`read2`) selector are only available via
//! STAR's `--runMode inputAlignmentsFromBAM`, which neither STAR-rs nor
//! rustar-aligner implement at align time. `Parameters::validate` rejects them
//! loudly rather than silently emitting the wrong tracks.

use crate::align::transcript::Transcript;
use crate::genome::Genome;

/// Per-chromosome coverage for the two strands, `Unique` and `UniqueMultiple` tracks.
pub struct Signal {
    /// `[chr][strand] = (unique, unique_multiple)`, each `chr_len + 1` bases (STAR's trailing zero).
    per_chr: Vec<[(Vec<f64>, Vec<f64>); 2]>,
    chr_len: Vec<usize>,
    chr_name: Vec<String>,
}

impl Signal {
    pub fn new(chr_name: &[String], chr_length: &[u64]) -> Self {
        let chr_len: Vec<usize> = chr_length.iter().map(|&l| l as usize).collect();
        let per_chr = chr_len
            .iter()
            .map(|&l| {
                [
                    (vec![0.0; l + 1], vec![0.0; l + 1]),
                    (vec![0.0; l + 1], vec![0.0; l + 1]),
                ]
            })
            .collect();
        Self {
            per_chr,
            chr_len,
            chr_name: chr_name.to_vec(),
        }
    }

    /// Accumulate covered genomic `M`-blocks (chr-local `(start, len)`) on one strand
    /// with `NH = n_tr`: `+1` to `Unique` iff `n_tr == 1`, always `+1/n_tr` to
    /// `UniqueMultiple`. Introns and deletions are simply the gaps between blocks.
    pub fn add_blocks(
        &mut self,
        chr: usize,
        strand: usize,
        blocks: &[(usize, usize)],
        n_tr: usize,
    ) {
        let (uniq, umult) = {
            let s = &mut self.per_chr[chr][strand];
            (&mut s.0, &mut s.1)
        };
        let w = 1.0 / n_tr as f64;
        for &(g0, len) in blocks {
            for p in g0..g0 + len {
                if n_tr == 1 {
                    uniq[p] += 1.0;
                }
                umult[p] += w;
            }
        }
    }

    /// Accumulate one reported alignment (`tr`) with `NH = n_tr`. Only exon blocks
    /// are covered; introns/deletions between exons are skipped by the block's
    /// genomic jump. `second_mate` selects STAR's `iStrand` rule (`signalFromBAM.cpp`):
    /// `iStrand = is_reverse == is_not_second_mate` — mate1 (and SE) use the
    /// alignment's own strand directly, mate2's is flipped.
    pub fn add_transcript(
        &mut self,
        genome: &Genome,
        tr: &Transcript,
        n_tr: usize,
        second_mate: bool,
    ) {
        let cs = genome.chr_start[tr.chr_idx];
        let blocks: Vec<(usize, usize)> = tr
            .exons
            .iter()
            .map(|ex| {
                (
                    (ex.genome_start - cs) as usize,
                    (ex.genome_end - ex.genome_start) as usize,
                )
            })
            .collect();
        let strand = usize::from(tr.is_reverse != second_mate);
        self.add_blocks(tr.chr_idx, strand, &blocks, n_tr);
    }

    /// bedGraph text for one track (`unique = true` for the `Unique` track, else
    /// `UniqueMultiple`) on one strand (0 = str1, 1 = str2), un-normalized
    /// (`--outWigNorm None`). STAR emits `chr start end value` per constant-coverage run.
    pub fn bedgraph(&self, unique: bool, strand: usize) -> String {
        self.bedgraph_fmt(unique, strand, fmt_g6)
    }

    /// bedGraph text normalized to reads-per-million (`--outWigNorm RPM`): each
    /// coverage value is scaled by `norm` (`1e6 / nReadsUnique`) and printed with 5
    /// decimals, matching STAR.
    pub fn bedgraph_rpm(&self, unique: bool, strand: usize, norm: f64) -> String {
        self.bedgraph_fmt(unique, strand, |v| format!("{:.5}", v * norm))
    }

    /// The shared bedGraph run-length emitter; `fmt` renders one coverage value.
    fn bedgraph_fmt(&self, unique: bool, strand: usize, fmt: impl Fn(f64) -> String) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        for chr in 0..self.chr_len.len() {
            let track = if unique {
                &self.per_chr[chr][strand].0
            } else {
                &self.per_chr[chr][strand].1
            };
            let mut prev = 0.0f64;
            for ig in 0..=self.chr_len[chr] {
                let new = if ig < track.len() { track[ig] } else { 0.0 };
                // Exact comparison is intentional: `new`/`prev` are read back
                // unmodified from the accumulator, so a run boundary is a real
                // change, not float noise (clippy::float_cmp).
                #[allow(clippy::float_cmp)]
                let changed = new != prev;
                if changed {
                    if prev != 0.0 {
                        let _ = writeln!(out, "{}\t{}", ig, fmt(prev));
                    }
                    if new != 0.0 {
                        let _ = write!(out, "{}\t{}\t", self.chr_name[chr], ig);
                    }
                    prev = new;
                }
            }
        }
        out
    }
}

/// Format a coverage value like C++'s default `ostream` for `--outWigNorm None`:
/// `%g` with 6 significant digits (NOT 6 decimal places). Integers print without a
/// decimal point; for `|v| >= 1` the number of fractional digits is
/// `5 - floor(log10|v|)` (so a total of 6 significant digits), then trailing zeros
/// and a trailing `.` are stripped. Values `|v| < 1` keep 6 fractional digits (the
/// leading `0.` does not count toward significant digits at these magnitudes).
/// Coverage sums are small and positive, so the scientific branch of `%g` never
/// triggers here.
pub(crate) fn fmt_g6(v: f64) -> String {
    // Exact comparison is intentional: this checks whether `v` is a whole number
    // (to print it without a decimal point), not an approximate-equality test
    // (clippy::float_cmp).
    #[allow(clippy::float_cmp)]
    let is_integer = v == v.trunc();
    if is_integer && v.abs() < 1e15 {
        return format!("{}", v as i64);
    }
    let decimals = if v.abs() >= 1.0 {
        (5 - v.abs().log10().floor() as i64).max(0) as usize
    } else {
        6
    };
    let mut s = format!("{v:.decimals$}");
    while s.ends_with('0') {
        s.pop();
    }
    if s.ends_with('.') {
        s.pop();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::Exon;
    use noodles::sam::alignment::record::cigar::{Op, op::Kind};

    fn tr_fwd(chr_idx: usize, genome_start: u64, len: u64) -> Transcript {
        Transcript {
            chr_idx,
            genome_start,
            genome_end: genome_start + len,
            is_reverse: false,
            exons: vec![Exon {
                genome_start,
                genome_end: genome_start + len,
                read_start: 0,
                read_end: len as usize,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, len as usize)],
            score: len as i32 - 1,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        }
    }

    #[test]
    fn fmt_g6_matches_cpp_percent_g() {
        assert_eq!(fmt_g6(1.0), "1");
        assert_eq!(fmt_g6(5.0), "5");
        assert_eq!(fmt_g6(0.5), "0.5");
        assert_eq!(fmt_g6(1.0 + 1.0 / 3.0), "1.33333");
        assert_eq!(fmt_g6(1.0 / 3.0), "0.333333");
        assert_eq!(fmt_g6(12.0 + 5.0 / 6.0), "12.8333");
        assert_eq!(fmt_g6(100.0 + 1.0 / 3.0), "100.333");
    }

    #[test]
    fn add_blocks_unique_and_multi() {
        let mut sig = Signal::new(&["chrI".to_string()], &[100]);
        sig.add_blocks(0, 0, &[(10, 5)], 1);
        assert_eq!(sig.bedgraph(true, 0), "chrI\t10\t15\t1\n");
        assert_eq!(sig.bedgraph(false, 0), "chrI\t10\t15\t1\n");
    }

    #[test]
    fn add_blocks_multimapper_splits_weight() {
        let mut sig = Signal::new(&["chrI".to_string()], &[100]);
        sig.add_blocks(0, 0, &[(10, 2)], 2);
        // Not unique -> no Unique track contribution, UniqueMultiple gets 1/2.
        assert_eq!(sig.bedgraph(true, 0), "");
        assert_eq!(sig.bedgraph(false, 0), "chrI\t10\t12\t0.5\n");
    }

    #[test]
    fn add_transcript_mate2_strand_is_flipped() {
        let mut sig = Signal::new(&["chrI".to_string()], &[100]);
        let genome = crate::genome::Genome {
            transform_blocks: None,
            sequence: vec![].into(),
            n_genome: 100,
            n_genome_real: 100,
            n_chr_real: 1,
            chr_name: vec!["chrI".to_string()],
            chr_length: vec![100],
            chr_start: vec![0],
        };
        let tr = tr_fwd(0, 10, 5); // forward alignment
        // SE / mate1: forward alignment -> str1 (strand 0).
        sig.add_transcript(&genome, &tr, 1, false);
        assert_eq!(sig.bedgraph(true, 0), "chrI\t10\t15\t1\n");
        assert_eq!(sig.bedgraph(true, 1), "");
        // mate2: forward alignment -> flipped to str2 (strand 1).
        let mut sig2 = Signal::new(&["chrI".to_string()], &[100]);
        sig2.add_transcript(&genome, &tr, 1, true);
        assert_eq!(sig2.bedgraph(true, 0), "");
        assert_eq!(sig2.bedgraph(true, 1), "chrI\t10\t15\t1\n");
    }
}
