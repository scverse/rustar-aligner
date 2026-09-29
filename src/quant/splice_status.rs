//! Bulk spliced / unspliced / ambiguous gene quantification
//! (`--quantMode GeneSplicing`, a rustar-aligner extension).
//!
//! Provenance: the classification is a port of STARsolo's Velocyto feature;
//! only the STAR source file names below still carry that name.
//!
//! STARsolo's `--soloFeatures Velocyto` classifies every read against every
//! annotated transcript that contains it, then collapses those per-transcript
//! calls into one of three categories per gene. STAR only exposes this for
//! single-cell runs. This module ports the same two steps so they can run on
//! bulk data, where each read (or read pair) plays the role of one UMI:
//!
//! 1. Per transcript: [`align_to_transcript_min_overlap`], a port of STAR's
//!    `alignToTranscriptMinOverlap` (`Transcriptome_classifyAlign.cpp`), with
//!    the hard-coded `minOverlapMinusOne = 6` (velocyto.py's `MIN_FLANK = 5`)
//!    and the 1 Mb intron cap. The transcript loop mirrors
//!    `Transcriptome::classifyAlign`: only transcripts that fully contain the
//!    alignment are tested, and the strand filter follows `--soloStrand`
//!    semantics (`Forward` = read 1 on the transcript strand).
//! 2. Per gene: [`collapse_gene_category`], a port of the per-UMI collapse in
//!    `SoloFeature_countVelocyto.cpp`: all transcripts must belong to one gene
//!    (else the read is multi-gene and dropped, as STAR does); only-exonic
//!    models give `spliced`, only-intronic / spanning models give
//!    `unspliced`, anything mixed gives `ambiguous`.
//!
//! As in STAR, "spliced" means "compatible only with mature mRNA" (the read
//! need not cross a junction), "unspliced" means "requires pre-mRNA", and
//! "ambiguous" means "compatible with both", which is where reads falling in
//! a retained intron that is exonic in another isoform end up.
//!
//! The per-gene table has the three categories for each of the three strand
//! conventions of `ReadsPerGene.out.tab` (unstranded, forward, reverse), so a
//! single run serves any library type.
use std::io::Write as _;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use noodles::sam::alignment::record::cigar::op::Kind;

use crate::align::read_align::PairedAlignment;
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::quant::transcriptome::{TrExon, TranscriptomeIndex};

/// STAR's hard-coded `minOverlapMinusOne` for Velocyto (`classifyAlign`:
/// "6 is the hard code minOverlapMinusOne, to agree with velocyto's
/// MIN_FLANK=5").
pub const MIN_OVERLAP_MINUS_ONE: u64 = 6;

/// STAR's intron size cap in `alignToTranscriptMinOverlap`: a read that is
/// intronic in an intron longer than this is not called intronic for that
/// transcript ("prevents large introns from swallowing small genes").
pub const MAX_INTRONIC_INTRON: u64 = 1_000_000;

/// STAR's `AlignVsTranscript` enum (`AlignVsTranscript.h`). The discriminants
/// are the bit positions used in the per-transcript type bitset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum AlignVsTranscript {
    /// Purely intronic.
    Intron = 0,
    /// Exonic and intronic blocks, none crossing an exon/intron boundary.
    ExonIntron = 1,
    /// At least one block crosses an exon/intron boundary.
    ExonIntronSpan = 2,
    /// Purely exonic.
    Concordant = 3,
}

impl AlignVsTranscript {
    /// STAR's `reAnn1` bitset for one transcript: the status bit, plus the
    /// `Intron` and `Concordant` bits for a span ("span is also considered
    /// intronic ... also considered exonic").
    pub fn type_bits(self) -> u8 {
        let mut bits = 1u8 << (self as u8);
        if self == Self::ExonIntronSpan {
            bits |= 1 << (Self::Intron as u8);
            bits |= 1 << (Self::Concordant as u8);
        }
        bits
    }
}

/// Splicing status of one read for one gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpliceStatus {
    /// Compatible only with exonic (mature) models.
    Spliced = 0,
    /// Requires an intronic (pre-mRNA) model.
    Unspliced = 1,
    /// Compatible with both.
    Ambiguous = 2,
}

/// Outcome of classifying one read under one strand convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadSplicing {
    /// No transcript on the selected strand contains the read (or every
    /// containing transcript rejected it).
    NoFeature,
    /// Containing transcripts belong to more than one gene.
    MultiGene,
    /// One gene, one category.
    Gene(u32, SpliceStatus),
}

/// Strand conventions, in output column order (same as `ReadsPerGene.out.tab`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandMode {
    /// Strand ignored.
    Unstranded = 0,
    /// Read 1 on the transcript strand (STARsolo `--soloStrand Forward`).
    Forward = 1,
    /// Read 1 opposite to the transcript strand (`--soloStrand Reverse`,
    /// dUTP / TruSeq Stranded libraries).
    Reverse = 2,
}

impl StrandMode {
    /// All three modes in output order.
    pub const ALL: [StrandMode; 3] = [Self::Unstranded, Self::Forward, Self::Reverse];

    fn name(self) -> &'static str {
        match self {
            Self::Unstranded => "unstranded",
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }

    /// STAR's filter in `classifyAlign`:
    /// `(trStr==1 ? aG.Str : 1-aG.Str) != pSolo.strand` rejects the transcript.
    fn keeps(self, tr_is_reverse: bool, read_is_reverse: bool) -> bool {
        match self {
            Self::Unstranded => true,
            Self::Forward => tr_is_reverse == read_is_reverse,
            Self::Reverse => tr_is_reverse != read_is_reverse,
        }
    }
}

/// A genomic alignment reduced to what the classifier reads: its
/// aligned blocks with indels merged (STAR expands a block across
/// `canonSJ` -1/-2), whether it has a splice junction (`sjYes`), and the
/// strand of read 1 (`aG.Str`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignBlocks {
    /// Chromosome index.
    pub chr_idx: usize,
    /// Blocks as absolute 0-based inclusive `(start, end)`, sorted by start.
    pub blocks: Vec<(u64, u64)>,
    /// Any mate crosses a splice junction.
    pub has_junction: bool,
    /// Strand of read 1.
    pub is_reverse: bool,
}

impl AlignBlocks {
    /// Blocks of a single alignment (SE read, or one mate).
    pub fn from_transcript(t: &Transcript) -> Self {
        let mut blocks = Vec::new();
        push_blocks(t, &mut blocks);
        blocks.sort_unstable();
        AlignBlocks {
            chr_idx: t.chr_idx,
            blocks,
            has_junction: t.n_junction > 0,
            is_reverse: t.is_reverse,
        }
    }

    /// Blocks of a read pair, both mates pooled like STAR's combined
    /// paired-end `Transcript`; the strand is mate 1's.
    pub fn from_pair(pair: &PairedAlignment) -> Self {
        let mut blocks = Vec::new();
        push_blocks(&pair.mate1_transcript, &mut blocks);
        push_blocks(&pair.mate2_transcript, &mut blocks);
        blocks.sort_unstable();
        AlignBlocks {
            chr_idx: pair.mate1_transcript.chr_idx,
            blocks,
            has_junction: pair.mate1_transcript.n_junction + pair.mate2_transcript.n_junction > 0,
            is_reverse: pair.mate1_transcript.is_reverse,
        }
    }

    /// Leftmost aligned base (0-based).
    fn start(&self) -> u64 {
        self.blocks.first().map_or(0, |b| b.0)
    }

    /// Rightmost aligned base (0-based, inclusive).
    fn end_incl(&self) -> u64 {
        self.blocks.iter().map(|b| b.1).max().unwrap_or(0)
    }
}

/// Append the indel-merged blocks of `t`: a new block starts only at a splice
/// junction (`N`). Falls back to the exon list when no CIGAR is attached.
fn push_blocks(t: &Transcript, out: &mut Vec<(u64, u64)>) {
    let Some(first) = t.exons.first() else {
        return;
    };
    if t.cigar.is_empty() {
        // No CIGAR (synthetic transcripts): exons are the blocks; merge
        // across insertions (read gap, no genome gap).
        let mut cur = (first.genome_start, first.genome_end);
        for e in &t.exons[1..] {
            if e.genome_start == cur.1 {
                cur.1 = e.genome_end;
            } else {
                out.push((cur.0, cur.1 - 1));
                cur = (e.genome_start, e.genome_end);
            }
        }
        out.push((cur.0, cur.1 - 1));
        return;
    }
    let mut pos = first.genome_start;
    let mut block_start = pos;
    for op in &t.cigar {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                pos += op.len() as u64;
            }
            Kind::Skip => {
                if pos > block_start {
                    out.push((block_start, pos - 1));
                }
                pos += op.len() as u64;
                block_start = pos;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if pos > block_start {
        out.push((block_start, pos - 1));
    }
}

/// Port of STAR's `alignToTranscriptMinOverlap`
/// (`Transcriptome_classifyAlign.cpp`). `exons` are the transcript's exons in
/// ascending order; the alignment must be contained in the transcript span
/// (the caller checks, as in STAR). Returns `None` for STAR's `-1`: a spliced
/// alignment touching an intron, or an intronic call inside an intron longer
/// than [`MAX_INTRONIC_INTRON`].
pub fn align_to_transcript_min_overlap(
    blocks: &[(u64, u64)],
    has_junction: bool,
    exons: &[TrExon],
    min_overlap_minus_one: u64,
) -> Option<AlignVsTranscript> {
    let m = min_overlap_minus_one;
    let n_ex = exons.len();
    let mut intronic = false;
    let mut exonic = false;
    let mut span = false;

    for &(bs, be) in blocks {
        // STAR `binarySearch1` over the interleaved exon start/end array,
        // halved: the last exon whose start is <= the block start.
        let ex1 = exons
            .partition_point(|e| e.genome_start <= bs)
            .checked_sub(1)?;
        if ex1 == n_ex - 1 {
            // Reached the last exon: with the alignment inside the transcript
            // (and blocks sorted), everything from here on is exonic.
            exonic = true;
            break;
        }
        if be - bs < m {
            continue; // block too short to call
        }
        let e_e = exons[ex1].genome_end - 1; // exon1 end (inclusive)
        let en_s = exons[ex1 + 1].genome_start; // exon2 start
        let en_e = exons[ex1 + 1].genome_end - 1; // exon2 end (inclusive)

        if bs + m <= e_e {
            // start is certainly in exon1
            if be <= e_e + m {
                exonic = true;
            } else {
                span = true;
            }
        } else if bs + m < en_s {
            // start is in intron1
            if be >= en_s + m {
                span = true;
            } else if be > e_e + m {
                if en_s - e_e > MAX_INTRONIC_INTRON {
                    return None;
                }
                intronic = true;
            }
        } else if be > en_e + m {
            // start too close to exon2 start; end certainly in intron2
            span = true;
        } else if be >= en_s + m {
            exonic = true;
        }

        if has_junction && (intronic || span) {
            return None; // a spliced alignment cannot overlap an intron
        }
    }

    Some(if span {
        AlignVsTranscript::ExonIntronSpan
    } else if !intronic {
        AlignVsTranscript::Concordant
    } else if exonic {
        AlignVsTranscript::ExonIntron
    } else {
        AlignVsTranscript::Intron
    })
}

/// Per-transcript type bits for every transcript that contains the
/// alignment, regardless of strand: `(transcript index, type bits)`. Mirrors
/// the transcript walk of STAR's `classifyAlign` (the strand filter is applied
/// later, per output column, by [`classify_read`]).
pub fn transcript_types(align: &AlignBlocks, idx: &TranscriptomeIndex, out: &mut Vec<(usize, u8)>) {
    out.clear();
    if align.blocks.is_empty() || idx.n_transcripts() == 0 {
        return;
    }
    let a_start = align.start();
    let a_end_excl = align.end_incl() + 1;
    let upper = idx.tr_starts_sorted.partition_point(|&s| s <= a_start);
    if upper == 0 {
        return;
    }
    let mut i = upper - 1;
    loop {
        if idx.tr_end_max_sorted[i] < a_end_excl {
            break;
        }
        let tr = idx.tr_order[i];
        if idx.tr_chr_idx[tr] == align.chr_idx
            && idx.tr_start[tr] <= a_start
            && idx.tr_end[tr] >= a_end_excl
            && let Some(status) = align_to_transcript_min_overlap(
                &align.blocks,
                align.has_junction,
                &idx.tr_exons[tr],
                MIN_OVERLAP_MINUS_ONE,
            )
        {
            out.push((tr, status.type_bits()));
        }
        if i == 0 {
            break;
        }
        i -= 1;
    }
}

/// Port of the per-UMI collapse in STAR's `SoloFeature_countVelocyto.cpp`
/// for one read: `types` are the `(transcript, type bits)` kept under one
/// strand convention.
pub fn collapse_gene_category(
    types: impl IntoIterator<Item = (usize, u8)>,
    idx: &TranscriptomeIndex,
) -> ReadSplicing {
    const INTRON: u8 = 1 << AlignVsTranscript::Intron as u8;
    const EXON_INTRON: u8 = 1 << AlignVsTranscript::ExonIntron as u8;
    const SPAN: u8 = 1 << AlignVsTranscript::ExonIntronSpan as u8;
    const CONCORDANT: u8 = 1 << AlignVsTranscript::Concordant as u8;

    let mut gene: Option<u32> = None;
    let mut exon_model = false;
    let mut intron_model = false;
    let mut span_model = true;
    let mut mixed_model = false;
    for (tr, ty) in types {
        let g = idx.tr_gene_idx[tr];
        match gene {
            None => gene = Some(g),
            Some(g0) if g0 != g => return ReadSplicing::MultiGene,
            Some(_) => {}
        }
        let has = |bit: u8| ty & bit != 0;
        mixed_model |= ((has(INTRON) && has(CONCORDANT)) || has(EXON_INTRON)) && !has(SPAN);
        span_model &= has(SPAN);
        exon_model |= has(CONCORDANT) && !has(INTRON) && !has(EXON_INTRON);
        intron_model |= has(INTRON) && !has(EXON_INTRON) && !has(CONCORDANT);
    }
    let Some(g) = gene else {
        return ReadSplicing::NoFeature;
    };
    let cat = if exon_model && !intron_model && !mixed_model {
        SpliceStatus::Spliced
    } else if span_model || ((intron_model || mixed_model) && !exon_model) {
        SpliceStatus::Unspliced
    } else {
        SpliceStatus::Ambiguous
    };
    ReadSplicing::Gene(g, cat)
}

/// Classify one uniquely mapped read (or pair) under the three strand
/// conventions, in [`StrandMode::ALL`] order.
pub fn classify_read(align: &AlignBlocks, idx: &TranscriptomeIndex) -> [ReadSplicing; 3] {
    let mut types = Vec::new();
    transcript_types(align, idx, &mut types);
    // STAR's per-UMI intersection works on transcript-sorted lists; for a
    // single read the order does not change the result, but sort anyway so
    // the multi-gene check sees the same first gene as STAR.
    types.sort_unstable_by_key(|&(tr, _)| tr);
    StrandMode::ALL.map(|mode| {
        collapse_gene_category(
            types
                .iter()
                .copied()
                .filter(|&(tr, _)| mode.keeps(idx.tr_strand[tr] == 2, align.is_reverse)),
            idx,
        )
    })
}

// ---------------------------------------------------------------------------
// Counters + output
// ---------------------------------------------------------------------------

/// Thread-safe counters for `ReadsPerGeneSplicing.out.tab`.
pub struct SplicingCounts {
    /// Per gene: `[strand mode][category]`, flattened as `mode * 3 + cat`.
    per_gene: Vec<[AtomicU64; 9]>,
    /// Reads that did not map (including too-many-loci), as in GeneCounts.
    pub n_unmapped: AtomicU64,
    /// Multi-mapping reads (not classified, as STAR's Velocyto requires a
    /// single alignment).
    pub n_multimapping: AtomicU64,
    /// Per strand mode: unique reads contained in no transcript.
    pub n_no_feature: [AtomicU64; 3],
    /// Per strand mode: unique reads whose transcripts span several genes.
    pub n_multi_gene: [AtomicU64; 3],
}

impl SplicingCounts {
    /// Zeroed counters for `n_genes` genes.
    pub fn new(n_genes: usize) -> Self {
        SplicingCounts {
            per_gene: (0..n_genes)
                .map(|_| std::array::from_fn(|_| AtomicU64::new(0)))
                .collect(),
            n_unmapped: AtomicU64::new(0),
            n_multimapping: AtomicU64::new(0),
            n_no_feature: std::array::from_fn(|_| AtomicU64::new(0)),
            n_multi_gene: std::array::from_fn(|_| AtomicU64::new(0)),
        }
    }

    fn record(&self, classes: &[ReadSplicing; 3]) {
        for (mode, class) in classes.iter().enumerate() {
            match *class {
                ReadSplicing::NoFeature => {
                    self.n_no_feature[mode].fetch_add(1, Ordering::Relaxed);
                }
                ReadSplicing::MultiGene => {
                    self.n_multi_gene[mode].fetch_add(1, Ordering::Relaxed);
                }
                ReadSplicing::Gene(g, cat) => {
                    self.per_gene[g as usize][mode * 3 + cat as usize]
                        .fetch_add(1, Ordering::Relaxed);
                }
            }
        }
    }

    /// Count a single-end read, with GeneCounts' unmapped / multimapping
    /// rules (`GeneCounts::count_se_read`).
    pub fn count_se_read(&self, transcripts: &[Transcript], idx: &TranscriptomeIndex) {
        match transcripts.len() {
            0 => {
                self.n_unmapped.fetch_add(1, Ordering::Relaxed);
            }
            1 => self.record(&classify_read(
                &AlignBlocks::from_transcript(&transcripts[0]),
                idx,
            )),
            _ => {
                self.n_multimapping.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Count a read pair, with GeneCounts' rules (`GeneCounts::count_pe_read`:
    /// half-mapped pairs count as unmapped).
    pub fn count_pe_read(
        &self,
        both_mapped: &[&PairedAlignment],
        unmapped: bool,
        idx: &TranscriptomeIndex,
    ) {
        if unmapped || both_mapped.is_empty() {
            self.n_unmapped.fetch_add(1, Ordering::Relaxed);
        } else if both_mapped.len() > 1 {
            self.n_multimapping.fetch_add(1, Ordering::Relaxed);
        } else {
            self.record(&classify_read(&AlignBlocks::from_pair(both_mapped[0]), idx));
        }
    }

    fn get(&self, g: usize, mode: usize, cat: usize) -> u64 {
        self.per_gene[g][mode * 3 + cat].load(Ordering::Relaxed)
    }

    /// Totals per `[mode][category]`.
    fn totals(&self) -> [[u64; 3]; 3] {
        let mut t = [[0u64; 3]; 3];
        for g in 0..self.per_gene.len() {
            for (mode, row) in t.iter_mut().enumerate() {
                for (cat, v) in row.iter_mut().enumerate() {
                    *v += self.get(g, mode, cat);
                }
            }
        }
        t
    }

    /// Write `ReadsPerGeneSplicing.out.tab`: a header line, then one line per
    /// gene (`geneInfo.tab` order) with spliced / unspliced / ambiguous
    /// counts for each strand convention.
    pub fn write_table(&self, path: &Path, idx: &TranscriptomeIndex) -> Result<(), Error> {
        let mut out = String::new();
        out.push_str("gene_id");
        for mode in StrandMode::ALL {
            for cat in ["spliced", "unspliced", "ambiguous"] {
                out.push('\t');
                out.push_str(mode.name());
                out.push('_');
                out.push_str(cat);
            }
        }
        out.push('\n');
        for (g, gene_id) in idx.gene_ids.iter().enumerate() {
            out.push_str(gene_id);
            for mode in 0..3 {
                for cat in 0..3 {
                    out.push('\t');
                    out.push_str(&self.get(g, mode, cat).to_string());
                }
            }
            out.push('\n');
        }
        std::fs::write(path, out).map_err(|e| Error::io(e, path))
    }

    /// Write `ReadsPerGeneSplicing.summary.tsv`: read accounting and the
    /// spliced / unspliced / ambiguous shares for each strand convention.
    pub fn write_summary(&self, path: &Path) -> Result<(), Error> {
        let mut f = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
        let io = |e| Error::io(e, path);
        let nu = self.n_unmapped.load(Ordering::Relaxed);
        let nm = self.n_multimapping.load(Ordering::Relaxed);
        let t = self.totals();
        let nf: [u64; 3] = std::array::from_fn(|m| self.n_no_feature[m].load(Ordering::Relaxed));
        let ng: [u64; 3] = std::array::from_fn(|m| self.n_multi_gene[m].load(Ordering::Relaxed));

        writeln!(f, "metric\tunstranded\tforward\treverse").map_err(io)?;
        let row = |f: &mut std::fs::File, name: &str, v: [u64; 3]| {
            writeln!(f, "{name}\t{}\t{}\t{}", v[0], v[1], v[2])
        };
        row(&mut f, "N_unmapped", [nu; 3]).map_err(io)?;
        row(&mut f, "N_multimapping", [nm; 3]).map_err(io)?;
        row(&mut f, "N_noFeature", nf).map_err(io)?;
        row(&mut f, "N_multiGene", ng).map_err(io)?;
        for (cat, name) in ["N_spliced", "N_unspliced", "N_ambiguous"]
            .iter()
            .enumerate()
        {
            row(&mut f, name, std::array::from_fn(|m| t[m][cat])).map_err(io)?;
        }
        for (cat, name) in [
            "fraction_spliced",
            "fraction_unspliced",
            "fraction_ambiguous",
        ]
        .iter()
        .enumerate()
        {
            let frac: [String; 3] = std::array::from_fn(|m| {
                let assigned: u64 = t[m].iter().sum();
                if assigned == 0 {
                    "NA".to_string()
                } else {
                    format!("{:.4}", t[m][cat] as f64 / assigned as f64)
                }
            });
            writeln!(f, "{name}\t{}\t{}\t{}", frac[0], frac[1], frac[2]).map_err(io)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::Exon;
    use crate::genome::Genome;
    use crate::junction::gtf::GtfRecord;
    use noodles::sam::alignment::record::cigar::Op;
    use std::collections::HashMap;

    const V_I: u8 = 1 << AlignVsTranscript::Intron as u8;
    const V_EI: u8 = 1 << AlignVsTranscript::ExonIntron as u8;
    const V_C: u8 = 1 << AlignVsTranscript::Concordant as u8;

    fn exons(v: &[(u64, u64)]) -> Vec<TrExon> {
        // (start, end_exclusive)
        let mut cum = 0u32;
        v.iter()
            .map(|&(s, e)| {
                let ex = TrExon {
                    genome_start: s,
                    genome_end: e,
                    ex_len_cum: cum,
                };
                cum += (e - s) as u32;
                ex
            })
            .collect()
    }

    /// Three-exon model: [100,200) [300,400) [500,600).
    fn model() -> Vec<TrExon> {
        exons(&[(100, 200), (300, 400), (500, 600)])
    }

    fn call(blocks: &[(u64, u64)], sj: bool) -> Option<AlignVsTranscript> {
        align_to_transcript_min_overlap(blocks, sj, &model(), MIN_OVERLAP_MINUS_ONE)
    }

    #[test]
    fn min_overlap_exonic_intronic_span() {
        use AlignVsTranscript::*;
        assert_eq!(call(&[(120, 169)], false), Some(Concordant));
        assert_eq!(call(&[(220, 269)], false), Some(Intron));
        // Crosses the exon1 / intron1 boundary by more than 6 bases.
        assert_eq!(call(&[(170, 219)], false), Some(ExonIntronSpan));
        // Crosses intron1 / exon2.
        assert_eq!(call(&[(270, 319)], false), Some(ExonIntronSpan));
        // Exonic block plus an intronic block (pair), no span.
        assert_eq!(call(&[(120, 169), (220, 269)], false), Some(ExonIntron));
        // Overhang of <= 6 bases into the intron is tolerated (MIN_FLANK).
        assert_eq!(call(&[(156, 205)], false), Some(Concordant));
        // Reaching the last exon short-circuits to exonic.
        assert_eq!(call(&[(520, 569)], false), Some(Concordant));
    }

    #[test]
    fn min_overlap_spliced_alignment_touching_intron_is_rejected() {
        // Spliced read whose second block runs into intron2.
        assert_eq!(call(&[(170, 199), (300, 449)], true), None);
        // Spliced read with both blocks exonic stays concordant.
        assert_eq!(
            call(&[(150, 199), (300, 349)], true),
            Some(AlignVsTranscript::Concordant)
        );
    }

    #[test]
    fn min_overlap_huge_intron_is_not_intronic() {
        let ex = exons(&[(0, 100), (2_000_000, 2_000_100)]);
        assert_eq!(
            align_to_transcript_min_overlap(&[(500_000, 500_049)], false, &ex, 6),
            None
        );
    }

    #[test]
    fn span_type_bits_include_intron_and_concordant() {
        let b = AlignVsTranscript::ExonIntronSpan.type_bits();
        assert_eq!(b, 0b1101);
        assert_eq!(AlignVsTranscript::Intron.type_bits(), V_I);
    }

    // --- Gene-level collapse and full classification on a GTF ---------------

    fn genome() -> Genome {
        // chr1 (10 kb), chrM (1 kb), chrX (2 kb), chrY (2 kb)
        let lens = [10_000u64, 1_000, 2_000, 2_000];
        let mut starts = vec![0u64];
        for l in lens {
            starts.push(starts.last().unwrap() + l);
        }
        let total = *starts.last().unwrap();
        Genome {
            transform_blocks: None,
            sequence: vec![0u8; total as usize].into(),
            n_genome: total,
            n_genome_real: total,
            n_chr_real: 4,
            chr_start: starts,
            chr_length: lens.to_vec(),
            chr_name: ["chr1", "chrM", "chrX", "chrY"].map(String::from).to_vec(),
        }
    }

    fn rec(chr: &str, s: u64, e: u64, strand: char, gene: &str, tr: &str) -> GtfRecord {
        let mut attributes = HashMap::new();
        attributes.insert("gene_id".to_string(), gene.to_string());
        attributes.insert("transcript_id".to_string(), tr.to_string());
        GtfRecord {
            seqname: chr.to_string(),
            feature: "exon".to_string(),
            start: s,
            end: e,
            strand,
            attributes,
        }
    }

    /// Synthetic annotation (GTF 1-based inclusive):
    /// - G1 (+): T1a exons 1001-1200, 1501-1700, 2001-2200 (fully spliced);
    ///   T1b exons 1001-1700 (retains intron1), 2001-2200.
    /// - G2 (-): antisense gene inside G1's intron 2: 1751-1950, single exon.
    /// - G3 (+): 5001-5100, 5201-5300 (boundary tests).
    /// - MT1 (+) on chrM: single-exon 101-600.
    /// - PAR gene on chrX and its GENCODE `_PAR_Y` copy on chrY.
    fn index() -> TranscriptomeIndex {
        let recs = vec![
            rec("chr1", 1001, 1200, '+', "G1", "T1a"),
            rec("chr1", 1501, 1700, '+', "G1", "T1a"),
            rec("chr1", 2001, 2200, '+', "G1", "T1a"),
            rec("chr1", 1001, 1700, '+', "G1", "T1b"),
            rec("chr1", 2001, 2200, '+', "G1", "T1b"),
            rec("chr1", 1751, 1950, '-', "G2", "T2"),
            rec("chr1", 5001, 5100, '+', "G3", "T3"),
            rec("chr1", 5201, 5300, '+', "G3", "T3"),
            rec("chrM", 101, 600, '+', "MT1", "TM"),
            rec("chrX", 101, 300, '+', "PAR1", "TP"),
            rec("chrX", 501, 700, '+', "PAR1", "TP"),
            rec("chrY", 101, 300, '+', "PAR1_PAR_Y", "TP_PAR_Y"),
            rec("chrY", 501, 700, '+', "PAR1_PAR_Y", "TP_PAR_Y"),
        ];
        TranscriptomeIndex::from_gtf_exons(&recs, &genome()).unwrap()
    }

    fn gene(idx: &TranscriptomeIndex, id: &str) -> u32 {
        idx.gene_ids.iter().position(|g| g == id).unwrap() as u32
    }

    /// SE transcript from 0-based chr-relative blocks `[s, e)`.
    fn aln(g: &Genome, chr: usize, blocks: &[(u64, u64)], rev: bool) -> Transcript {
        use noodles::sam::alignment::record::cigar::op::Kind as K;
        let off = g.chr_start[chr];
        let mut cigar = Vec::new();
        let mut ex = Vec::new();
        let mut rpos = 0usize;
        for (i, &(s, e)) in blocks.iter().enumerate() {
            if i > 0 {
                cigar.push(Op::new(K::Skip, (s - blocks[i - 1].1) as usize));
            }
            cigar.push(Op::new(K::Match, (e - s) as usize));
            ex.push(Exon {
                genome_start: off + s,
                genome_end: off + e,
                read_start: rpos,
                read_end: rpos + (e - s) as usize,
                i_frag: 0,
            });
            rpos += (e - s) as usize;
        }
        Transcript {
            chr_idx: chr,
            genome_start: off + blocks[0].0,
            genome_end: off + blocks.last().unwrap().1,
            is_reverse: rev,
            exons: ex,
            cigar,
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: (blocks.len() - 1) as u32,
            junction_motifs: vec![],
            junction_annotated: vec![],
        }
    }

    fn classify(t: &Transcript, idx: &TranscriptomeIndex) -> [ReadSplicing; 3] {
        classify_read(&AlignBlocks::from_transcript(t), idx)
    }

    #[test]
    fn constitutive_exon_read_is_spliced() {
        let g = genome();
        let idx = index();
        let g1 = gene(&idx, "G1");
        let r = classify(&aln(&g, 0, &[(2050, 2100)], false), &idx);
        assert_eq!(r[0], ReadSplicing::Gene(g1, SpliceStatus::Spliced));
        assert_eq!(r[1], ReadSplicing::Gene(g1, SpliceStatus::Spliced));
        // Reverse library: a + read is antisense to G1, no feature.
        assert_eq!(r[2], ReadSplicing::NoFeature);
    }

    #[test]
    fn junction_read_is_spliced() {
        let g = genome();
        let idx = index();
        let g1 = gene(&idx, "G1");
        let r = classify(&aln(&g, 0, &[(1650, 1700), (2000, 2050)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g1, SpliceStatus::Spliced));
    }

    #[test]
    fn retained_intron_read_is_ambiguous() {
        // Inside intron1 of T1a, exonic in T1b: compatible with both.
        let g = genome();
        let idx = index();
        let g1 = gene(&idx, "G1");
        let r = classify(&aln(&g, 0, &[(1300, 1350)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g1, SpliceStatus::Ambiguous));
    }

    #[test]
    fn constitutive_intron_read_is_unspliced() {
        // Intron2 of G1 (1700..2000), sense strand, clear of G2 (1750..1950).
        let g = genome();
        let idx = index();
        let g1 = gene(&idx, "G1");
        let r = classify(&aln(&g, 0, &[(1955, 1990)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g1, SpliceStatus::Unspliced));
        // Exon/intron boundary read of a constitutive exon: unspliced.
        let r = classify(&aln(&g, 0, &[(1680, 1730)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g1, SpliceStatus::Unspliced));
    }

    #[test]
    fn antisense_overlapping_gene_is_split_by_strand() {
        // A read in G2's exon, which lies in G1's intron 2 on the other strand.
        let g = genome();
        let idx = index();
        let g1 = gene(&idx, "G1");
        let g2 = gene(&idx, "G2");
        // Read on - strand = sense for G2.
        let r = classify(&aln(&g, 0, &[(1800, 1850)], true), &idx);
        assert_eq!(r[0], ReadSplicing::MultiGene); // unstranded: G1 intron + G2 exon
        assert_eq!(r[1], ReadSplicing::Gene(g2, SpliceStatus::Spliced));
        assert_eq!(r[2], ReadSplicing::Gene(g1, SpliceStatus::Unspliced));
    }

    #[test]
    fn read_outside_transcripts_is_no_feature() {
        let g = genome();
        let idx = index();
        let r = classify(&aln(&g, 0, &[(8000, 8050)], false), &idx);
        assert_eq!(r, [ReadSplicing::NoFeature; 3]);
        // Protruding past a transcript end: not contained, no feature.
        let r = classify(&aln(&g, 0, &[(5280, 5330)], false), &idx);
        assert_eq!(r, [ReadSplicing::NoFeature; 3]);
    }

    #[test]
    fn boundary_exon_reads() {
        let g = genome();
        let idx = index();
        let g3 = gene(&idx, "G3");
        // First base of the first exon to within the exon.
        let r = classify(&aln(&g, 0, &[(5000, 5050)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g3, SpliceStatus::Spliced));
        // Last bases of the last exon.
        let r = classify(&aln(&g, 0, &[(5250, 5300)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g3, SpliceStatus::Spliced));
    }

    #[test]
    fn single_exon_chrm_gene_is_spliced() {
        let g = genome();
        let idx = index();
        let mt = gene(&idx, "MT1");
        let r = classify(&aln(&g, 1, &[(200, 250)], false), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(mt, SpliceStatus::Spliced));
    }

    #[test]
    fn gencode_par_y_copies_are_separate_genes() {
        let g = genome();
        let idx = index();
        let x = gene(&idx, "PAR1");
        let y = gene(&idx, "PAR1_PAR_Y");
        let rx = classify(&aln(&g, 2, &[(350, 400)], false), &idx);
        let ry = classify(&aln(&g, 3, &[(350, 400)], false), &idx);
        assert_eq!(rx[1], ReadSplicing::Gene(x, SpliceStatus::Unspliced));
        assert_eq!(ry[1], ReadSplicing::Gene(y, SpliceStatus::Unspliced));
    }

    #[test]
    fn pair_blocks_are_pooled() {
        // Mate 1 exonic in exon2, mate 2 intronic in intron2: ExonIntron in
        // every G1 model -> unspliced.
        let g = genome();
        let idx = index();
        let g1 = gene(&idx, "G1");
        let pair = PairedAlignment {
            mate1_transcript: aln(&g, 0, &[(1550, 1600)], false),
            mate2_transcript: aln(&g, 0, &[(1955, 1990)], true),
            mate1_region: (0, 50),
            mate2_region: (0, 35),
            is_proper_pair: true,
            insert_size: 440,
            combined_wt_score: 0,
            combined_n_match: 85,
        };
        let r = classify_read(&AlignBlocks::from_pair(&pair), &idx);
        assert_eq!(r[1], ReadSplicing::Gene(g1, SpliceStatus::Unspliced));
    }

    #[test]
    fn indels_do_not_split_blocks() {
        use noodles::sam::alignment::record::cigar::op::Kind as K;
        let g = genome();
        let mut t = aln(&g, 0, &[(2050, 2100)], false);
        t.cigar = vec![
            Op::new(K::SoftClip, 3),
            Op::new(K::Match, 20),
            Op::new(K::Deletion, 2),
            Op::new(K::Match, 10),
            Op::new(K::Insertion, 1),
            Op::new(K::Match, 18),
        ];
        let b = AlignBlocks::from_transcript(&t);
        assert_eq!(b.blocks, vec![(2050, 2099)]);
        assert!(!b.has_junction);
    }

    #[test]
    fn collapse_matches_star_rules() {
        let idx = index();
        let t1a = idx.tr_ids.iter().position(|t| t == "T1a").unwrap();
        let t1b = idx.tr_ids.iter().position(|t| t == "T1b").unwrap();
        let t2 = idx.tr_ids.iter().position(|t| t == "T2").unwrap();
        let g1 = gene(&idx, "G1");
        let span = AlignVsTranscript::ExonIntronSpan.type_bits();
        let cases: &[(&[(usize, u8)], ReadSplicing)] = &[
            (&[], ReadSplicing::NoFeature),
            (&[(t1a, V_C), (t2, V_C)], ReadSplicing::MultiGene),
            (
                &[(t1a, V_C), (t1b, V_C)],
                ReadSplicing::Gene(g1, SpliceStatus::Spliced),
            ),
            (
                &[(t1a, V_I), (t1b, V_I)],
                ReadSplicing::Gene(g1, SpliceStatus::Unspliced),
            ),
            (
                &[(t1a, V_EI)],
                ReadSplicing::Gene(g1, SpliceStatus::Unspliced),
            ),
            (
                &[(t1a, span), (t1b, span)],
                ReadSplicing::Gene(g1, SpliceStatus::Unspliced),
            ),
            (
                &[(t1a, V_I), (t1b, V_C)],
                ReadSplicing::Gene(g1, SpliceStatus::Ambiguous),
            ),
            // STAR: a span in one model plus an only-exonic model is "spliced"
            // (the span sets neither exonModel, intronModel nor mixedModel).
            (
                &[(t1a, span), (t1b, V_C)],
                ReadSplicing::Gene(g1, SpliceStatus::Spliced),
            ),
        ];
        for (types, want) in cases {
            assert_eq!(
                collapse_gene_category(types.iter().copied(), &idx),
                *want,
                "{types:?}"
            );
        }
    }

    #[test]
    fn counts_table_and_summary() {
        let g = genome();
        let idx = index();
        let counts = SplicingCounts::new(idx.gene_ids.len());
        // spliced (constitutive exon), unspliced (intron 2), ambiguous
        // (retained intron), unmapped, multimapper, no feature.
        counts.count_se_read(&[aln(&g, 0, &[(2050, 2100)], false)], &idx);
        counts.count_se_read(&[aln(&g, 0, &[(1955, 1990)], false)], &idx);
        counts.count_se_read(&[aln(&g, 0, &[(1300, 1350)], false)], &idx);
        counts.count_se_read(&[], &idx);
        let t = aln(&g, 0, &[(2050, 2100)], false);
        counts.count_se_read(&[t.clone(), t], &idx);
        counts.count_se_read(&[aln(&g, 0, &[(8000, 8050)], false)], &idx);

        let dir = tempfile::tempdir().unwrap();
        let tab = dir.path().join("t.tab");
        counts.write_table(&tab, &idx).unwrap();
        let tab = std::fs::read_to_string(tab).unwrap();
        let mut lines = tab.lines();
        assert_eq!(
            lines.next().unwrap(),
            "gene_id\tunstranded_spliced\tunstranded_unspliced\tunstranded_ambiguous\t\
             forward_spliced\tforward_unspliced\tforward_ambiguous\t\
             reverse_spliced\treverse_unspliced\treverse_ambiguous"
        );
        assert_eq!(lines.next().unwrap(), "G1\t1\t1\t1\t1\t1\t1\t0\t0\t0");
        assert_eq!(lines.count(), idx.gene_ids.len() - 1);

        let sum = dir.path().join("s.tsv");
        counts.write_summary(&sum).unwrap();
        let sum = std::fs::read_to_string(sum).unwrap();
        let get = |k: &str| {
            sum.lines()
                .find(|l| l.split('\t').next() == Some(k))
                .unwrap()
                .split('\t')
                .skip(1)
                .map(str::to_string)
                .collect::<Vec<_>>()
        };
        assert_eq!(get("N_unmapped"), ["1", "1", "1"]);
        assert_eq!(get("N_multimapping"), ["1", "1", "1"]);
        assert_eq!(get("N_noFeature"), ["1", "1", "4"]);
        assert_eq!(get("N_unspliced"), ["1", "1", "0"]);
        assert_eq!(get("fraction_unspliced"), ["0.3333", "0.3333", "NA"]);
    }
}
