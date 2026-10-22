//! Per-read gene assignment for the STARsolo `Gene` feature (Phase 14.3).
//!
//! A read is assigned to a gene by intersecting the gene model with the read's
//! alignment(s). Following STARsolo's `Gene` feature under the default
//! `--soloMultiMappers Unique`, the read's gene set is the UNION of genes
//! concordant with any of its alignments (strand-filtered by `--soloStrand`):
//! exactly one gene → assigned; zero → no feature; more than one → ambiguous.
//! A multi-locus read whose loci all fall in the same gene is therefore still
//! gene-unique, unlike `--quantMode GeneCounts` which drops all multimappers.

use crate::align::transcript::Transcript;
use crate::quant::GeneAnnotation;
use crate::quant::transcriptome::{BlockGap, TranscriptomeIndex, block_gaps};
use std::cell::RefCell;
use std::str::FromStr;

/// `--soloStrand`: orientation of the cDNA read relative to its gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SoloStrand {
    /// Read maps to the sense (same) strand as the gene (10x 3'/5', default).
    #[default]
    Forward,
    /// Read maps to the antisense (opposite) strand.
    Reverse,
    /// Strand is ignored.
    Unstranded,
}

impl FromStr for SoloStrand {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Forward" => Ok(Self::Forward),
            "Reverse" => Ok(Self::Reverse),
            "Unstranded" => Ok(Self::Unstranded),
            _ => Err(format!(
                "unknown soloStrand '{s}'; expected Forward, Reverse, or Unstranded"
            )),
        }
    }
}

/// A STARsolo `--soloFeatures` value that quantifies genes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SoloFeature {
    /// Exonic counting: a read counts toward a gene only if it overlaps an exon.
    Gene,
    /// Full gene-body counting (CellRanger `include-introns`): a read counts if
    /// it overlaps the gene locus, including purely intronic reads.
    GeneFull,
}

impl SoloFeature {
    /// Output sub-directory name (`Solo.out/<dir>/raw/`).
    pub fn dir_name(self) -> &'static str {
        match self {
            SoloFeature::Gene => "Gene",
            SoloFeature::GeneFull => "GeneFull",
        }
    }
}

impl FromStr for SoloFeature {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Gene" => Ok(Self::Gene),
            "GeneFull" => Ok(Self::GeneFull),
            _ => Err(format!(
                "unsupported soloFeature '{s}'; supported: Gene, GeneFull"
            )),
        }
    }
}

/// Outcome of assigning a read to a gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneAssignment {
    /// Concordant with exactly one gene (the assigned gene index).
    Gene(u32),
    /// Mapped but overlaps no gene on the selected strand.
    NoFeature,
    /// Overlaps more than one gene → not uniquely assignable.
    Ambiguous,
    /// Read did not map (no transcripts / too many loci).
    Unmapped,
}

/// Whether gene `g` is kept for read alignment `tr` under `strand`.
#[inline]
fn strand_keeps(strand: SoloStrand, gene_is_reverse: bool, read_is_reverse: bool) -> bool {
    match strand {
        SoloStrand::Unstranded => true,
        SoloStrand::Forward => gene_is_reverse == read_is_reverse,
        SoloStrand::Reverse => gene_is_reverse != read_is_reverse,
    }
}

/// The transcript models of STARsolo's `Gene` feature, each transcript's gene
/// given as an index into the solo gene list (`GeneAnnotation`, the order of
/// `features.tsv`).
#[derive(Default)]
pub struct GeneModel {
    pub transcripts: TranscriptomeIndex,
    /// Solo gene index per transcript; `u32::MAX` if the gene is not listed.
    tr_gene: Vec<u32>,
}

impl GeneModel {
    /// Solo gene index of transcript `tr`; `u32::MAX` if its gene is not listed.
    pub fn tr_gene(&self, tr: usize) -> u32 {
        self.tr_gene[tr]
    }

    /// A model that knows only each transcript's gene, for tests of the
    /// per-molecule steps.
    #[cfg(test)]
    pub(crate) fn with_tr_genes(tr_gene: Vec<u32>) -> Self {
        Self {
            transcripts: TranscriptomeIndex::default(),
            tr_gene,
        }
    }

    pub fn new(transcripts: TranscriptomeIndex, gene_ann: &GeneAnnotation) -> Self {
        let pos: std::collections::HashMap<&str, u32> = gene_ann
            .gene_ids
            .iter()
            .enumerate()
            .map(|(i, g)| (g.as_str(), i as u32))
            .collect();
        let tr_gene = transcripts
            .tr_gene_idx
            .iter()
            .map(|&g| {
                transcripts
                    .gene_ids
                    .get(g as usize)
                    .and_then(|id| pos.get(id.as_str()).copied())
                    .unwrap_or(u32::MAX)
            })
            .collect();
        GeneModel {
            transcripts,
            tr_gene,
        }
    }
}

/// STAR's `Gene` annotation of a read (`ReadAnnotFeature`): the genes of all
/// its alignments (`fSet`) and of each alignment (`fAlign`), both sorted.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GeneSets {
    pub genes: Vec<u32>,
    pub per_align: Vec<Vec<u32>>,
}

/// STAR's `binarySearch_leLeft` (`serviceFuns.cpp:212-235`) over the
/// transcript-relative exon starts and inclusive ends `[s0, e0, s1, e1, ...]`:
/// `None` if `x` is outside `[X[0], X[N-1]]`, else the index of the element
/// `<= x`, the leftmost one if several are equal. Like STAR's loop it never
/// returns `N-1`.
fn exon_bound_le_left(x: u64, bound: impl Fn(usize) -> u64, n: usize) -> Option<usize> {
    if x > bound(n - 1) || x < bound(0) {
        return None;
    }
    let (mut i1, mut i2) = (0usize, n - 1);
    while i2 > i1 + 1 {
        let i3 = usize::midpoint(i1, i2);
        if bound(i3) > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while i1 > 0 && x == bound(i1 - 1) {
        i1 -= 1;
    }
    Some(i1)
}

/// STAR's `alignToTranscript` (`Transcriptome_classifyAlign.cpp:8-90`), only
/// as far as the `Gene` feature needs it: is the alignment `Concordant` with
/// transcript `tr`, i.e. every block inside one of its exons and every splice
/// junction exactly one of its introns.
fn concordant_with_transcript(
    aln: &Transcript,
    gaps: &[BlockGap],
    idx: &TranscriptomeIndex,
    tr: usize,
) -> bool {
    let tr_s = idx.tr_start[tr];
    let exons = &idx.tr_exons[tr];
    let n_bounds = 2 * exons.len();
    // exSE1: exon starts and inclusive ends, relative to the transcript start.
    let bound = |k: usize| {
        let e = &exons[k / 2];
        if k.is_multiple_of(2) {
            e.genome_start - tr_s
        } else {
            e.genome_end - 1 - tr_s
        }
    };
    let (mut align_intronic, mut align_spans_exon_intr, mut sj_concordant) = (false, false, true);
    let (mut ex1, mut b_e, mut e_e, mut en_s) = (0usize, 0u64, 0u64, 0u64);
    for (iab, block) in aln.exons.iter().enumerate() {
        let b_e_prev = b_e;
        if block.genome_start < tr_s {
            return false; // a 2nd mate protruding to the left of the 1st
        }
        let b_s = block.genome_start - tr_s;
        b_e = block.genome_end - 1 - tr_s;
        let mate_start = iab == 0 || block.i_frag != aln.exons[iab - 1].i_frag; // canonSJ -3
        if mate_start {
            match exon_bound_le_left(b_s, bound, n_bounds) {
                Some(i) => ex1 = i / 2,
                None => return false,
            }
        } else if gaps.get(iab - 1) == Some(&BlockGap::Junction) {
            if b_e_prev == e_e && b_s == en_s {
                ex1 += 1; // the junction is this transcript's intron
            } else {
                sj_concordant = false;
                match exon_bound_le_left(b_s, bound, n_bounds) {
                    Some(i) => ex1 = i / 2,
                    None => return false,
                }
            }
        }
        e_e = bound(2 * ex1 + 1);
        en_s = if ex1 + 1 < exons.len() {
            bound(2 * (ex1 + 1))
        } else {
            0
        };
        if b_s <= e_e {
            if b_e > e_e {
                align_spans_exon_intr = true;
            }
        } else {
            if b_e >= en_s {
                align_spans_exon_intr = true;
            }
            align_intronic = true;
        }
    }
    sj_concordant && !align_spans_exon_intr && !align_intronic
}

/// STAR's `Transcriptome::classifyAlign` for the `Gene` feature
/// (`Transcriptome_classifyAlign.cpp:173-275`): for each alignment, walk back
/// from the last transcript starting at or before it while the running
/// maximum of transcript ends still reaches the alignment's end, and collect
/// the genes of the transcripts on the `--soloStrand` strand that contain it
/// and are concordant with it.
pub fn classify_align(
    transcripts: &[Transcript],
    model: &GeneModel,
    strand: SoloStrand,
) -> GeneSets {
    let idx = &model.transcripts;
    let mut sets = GeneSets {
        genes: Vec::new(),
        per_align: vec![Vec::new(); transcripts.len()],
    };
    if idx.tr_starts_sorted.is_empty() {
        return sets;
    }
    for (iag, aln) in transcripts.iter().enumerate() {
        let (Some(first), Some(last)) = (aln.exons.first(), aln.exons.last()) else {
            continue;
        };
        // binarySearch1a: the last transcript whose start is <= the align start.
        let pp = idx
            .tr_starts_sorted
            .partition_point(|&s| s <= first.genome_start);
        if pp == 0 {
            continue; // before every transcript
        }
        let a_g_end = last.genome_end - 1;
        let gaps = block_gaps(aln);
        let mut k = pp - 1;
        loop {
            let tr = idx.tr_order[k];
            // `trStr==1 ? aG.Str : 1-aG.Str` against `pSolo.strand`.
            let read_str = u32::from(aln.is_reverse);
            let rel_str = if idx.tr_strand[tr] == 1 {
                read_str
            } else {
                1 - read_str
            };
            let strand_ok = match strand {
                SoloStrand::Unstranded => true,
                SoloStrand::Forward => rel_str == 0,
                SoloStrand::Reverse => rel_str == 1,
            };
            if a_g_end < idx.tr_end[tr]
                && strand_ok
                && concordant_with_transcript(aln, &gaps, idx, tr)
                && model.tr_gene[tr] != u32::MAX
            {
                let g = model.tr_gene[tr];
                sets.genes.push(g);
                sets.per_align[iag].push(g);
            }
            // `while (trEmax[tr1]>=aGend && tr1>0)`
            if !(idx.tr_end_max_sorted[k] > a_g_end && k > 0) {
                break;
            }
            k -= 1;
        }
        sets.per_align[iag].sort_unstable();
        sets.per_align[iag].dedup();
    }
    sets.genes.sort_unstable();
    sets.genes.dedup();
    sets
}

/// STAR's `GeneFull` annotation of a read (`Transcriptome_geneFullAlignOverlap.cpp`):
/// per alignment, the genes whose body overlaps one of its blocks on the strand
/// `--soloStrand` keeps.
pub fn gene_full_sets(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    strand: SoloStrand,
) -> GeneSets {
    let mut sets = GeneSets {
        genes: Vec::new(),
        per_align: vec![Vec::new(); transcripts.len()],
    };
    let mut raw = Vec::new();
    for (iag, tr) in transcripts.iter().enumerate() {
        gene_ann.overlapping_genes_full_into(tr, &mut raw);
        for &g in &raw {
            if strand_keeps(strand, gene_ann.gene_is_reverse[g], tr.is_reverse) {
                sets.per_align[iag].push(g as u32);
                sets.genes.push(g as u32);
            }
        }
    }
    sets.genes.sort_unstable();
    sets.genes.dedup();
    sets
}

/// The Velocyto annotation of one uniquely mapped alignment (STAR's
/// `trVelocytoType`, `Transcriptome_classifyAlign.cpp:227-258`; strand
/// filter `:207-209`): for every transcript that contains it on a strand
/// `--soloStrand` keeps, the transcript (an index into `model.transcripts`)
/// and STAR's `reAnn1` type bits from `alignToTranscriptMinOverlap` with
/// velocyto's 5-base flank.
/// A spliced alignment leaves out a transcript whose intron one of its blocks
/// runs into (its junctions are not compared with the transcript's). Sorted
/// by transcript, as `SoloReadFeature_record.cpp:164-165` sorts it.
///
/// The alignment's blocks are its exons merged across indels, mate by mate,
/// in transcript order (STAR expands a block over `canonSJ` -1/-2, not the
/// mate gap -3); `sjYes` is whether any mate crosses a junction.
pub fn velocyto_transcript_types(
    aln: &Transcript,
    model: &GeneModel,
    strand: SoloStrand,
) -> Vec<(u32, u8)> {
    use crate::quant::splice_status::{MIN_OVERLAP_MINUS_ONE, align_to_transcript_min_overlap};
    let idx = &model.transcripts;
    let mut out = Vec::new();
    let (Some(first), Some(last)) = (aln.exons.first(), aln.exons.last()) else {
        return out;
    };
    if idx.tr_starts_sorted.is_empty() {
        return out;
    }
    // Blocks: exons merged across indels within a mate.
    let gaps = block_gaps(aln);
    let mut blocks: Vec<(u64, u64)> = Vec::with_capacity(aln.exons.len());
    let mut has_junction = false;
    for (i, e) in aln.exons.iter().enumerate() {
        let joins = i > 0
            && e.i_frag == aln.exons[i - 1].i_frag
            && gaps.get(i - 1) != Some(&BlockGap::Junction);
        if i > 0 && e.i_frag == aln.exons[i - 1].i_frag && !joins {
            has_junction = true;
        }
        match blocks.last_mut() {
            Some(b) if joins => b.1 = e.genome_end - 1,
            _ => blocks.push((e.genome_start, e.genome_end - 1)),
        }
    }
    // STAR's transcript walk (`binarySearch1a`, then back while `trEmax` covers
    // the alignment end), with the `--soloStrand` filter applied before the
    // Velocyto call, as in `classifyAlign`.
    let pp = idx
        .tr_starts_sorted
        .partition_point(|&s| s <= first.genome_start);
    if pp == 0 {
        return out;
    }
    let a_g_end = last.genome_end - 1;
    let read_str = u32::from(aln.is_reverse);
    let mut k = pp - 1;
    loop {
        let tr = idx.tr_order[k];
        let rel_str = if idx.tr_strand[tr] == 1 {
            read_str
        } else {
            1 - read_str
        };
        let strand_ok = match strand {
            SoloStrand::Unstranded => true,
            SoloStrand::Forward => rel_str == 0,
            SoloStrand::Reverse => rel_str == 1,
        };
        if a_g_end < idx.tr_end[tr]
            && strand_ok
            && let Some(status) = align_to_transcript_min_overlap(
                &blocks,
                has_junction,
                &idx.tr_exons[tr],
                MIN_OVERLAP_MINUS_ONE,
            )
        {
            out.push((tr as u32, status.type_bits()));
        }
        if !(idx.tr_end_max_sorted[k] > a_g_end && k > 0) {
            break;
        }
        k -= 1;
    }
    out.sort_unstable_by_key(|&(tr, _)| tr);
    out
}

/// CellRanger-style positional region of a uniquely-mapped read (independent of
/// strand): which genomic region the read falls in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    /// Overlaps ≥1 annotated exon.
    Exonic,
    /// Overlaps a gene body but no exon (purely intronic).
    Intronic,
    /// Overlaps no gene body.
    Intergenic,
}

/// Everything one read's alignment set tells us, computed in a single pass over
/// the gene model (the two overlap queries are shared between the per-feature
/// gene assignment and the region classification, so this costs no more than the
/// old two `assign_gene_se` calls).
#[derive(Debug, Clone)]
pub struct ReadClass {
    /// Sense-strand exonic gene assignment (the `Gene` feature). `Unmapped` if
    /// exon overlap was not requested.
    pub gene: GeneAssignment,
    /// Sense-strand gene-body assignment (the `GeneFull` feature). `Unmapped` if
    /// body overlap was not requested.
    pub gene_full: GeneAssignment,
    /// Positional region (only when both exon + body overlap were computed).
    pub region: Option<Region>,
    /// Read maps to a gene body on the antisense strand and to none on the sense
    /// strand (CellRanger's "Reads Mapped Antisense to Gene").
    pub antisense: bool,
    /// Multi-gene set for the `Gene` feature (the sense exon genes), populated
    /// only when `want_multi` and the read is gene-ambiguous (>1 gene). Used by
    /// `--soloMultiMappers` to distribute the read; empty otherwise.
    pub gene_multi: Vec<u32>,
    /// Multi-gene set for the `GeneFull` feature (sense body genes).
    pub gene_full_multi: Vec<u32>,
}

fn assignment_of(sense_genes: &[usize]) -> GeneAssignment {
    match sense_genes.len() {
        0 => GeneAssignment::NoFeature,
        1 => GeneAssignment::Gene(sense_genes[0] as u32),
        _ => GeneAssignment::Ambiguous,
    }
}

/// Classify a read in one pass: sense-strand `Gene`/`GeneFull` assignments plus
/// the CellRanger-style positional region + antisense flag. `want_exon` /
/// `want_body` skip the corresponding overlap query when a feature is not needed.
pub fn classify_read(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    model: &GeneModel,
    strand: SoloStrand,
    want_exon: bool,
    want_body: bool,
    want_multi: bool,
) -> ReadClass {
    if transcripts.is_empty() {
        return ReadClass {
            gene: GeneAssignment::Unmapped,
            gene_full: GeneAssignment::Unmapped,
            region: None,
            antisense: false,
            gene_multi: Vec::new(),
            gene_full_multi: Vec::new(),
        };
    }

    thread_local! {
        static RAW: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) };
        static EXON_S: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) };
        static BODY_S: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) };
    }

    RAW.with(|rb| {
        EXON_S.with(|eb| {
            BODY_S.with(|bb| {
                let mut raw = rb.borrow_mut();
                let mut exon_s = eb.borrow_mut();
                let mut body_s = bb.borrow_mut();
                exon_s.clear();
                body_s.clear();
                // `*_any` track positional (either-strand) overlap for the region;
                // `body_anti_any` tracks an antisense-only body hit.
                let (mut exon_any, mut body_any, mut body_anti_any) = (false, false, false);

                // The `Gene` feature's genes: those whose transcripts the
                // alignments are concordant with (STAR's `classifyAlign`).
                if want_exon {
                    exon_s.extend(
                        classify_align(transcripts, model, strand)
                            .genes
                            .iter()
                            .map(|&g| g as usize),
                    );
                }
                for tr in transcripts {
                    if want_exon {
                        gene_ann.overlapping_genes_into(tr, &mut raw);
                        exon_any |= !raw.is_empty();
                    }
                    if want_body {
                        gene_ann.overlapping_genes_full_into(tr, &mut raw);
                        for &g in raw.iter() {
                            body_any = true;
                            if strand_keeps(strand, gene_ann.gene_is_reverse[g], tr.is_reverse) {
                                body_s.push(g);
                            } else {
                                body_anti_any = true;
                            }
                        }
                    }
                }
                exon_s.sort_unstable();
                exon_s.dedup();
                body_s.sort_unstable();
                body_s.dedup();

                let region = if want_exon && want_body {
                    Some(if exon_any {
                        Region::Exonic
                    } else if body_any {
                        Region::Intronic
                    } else {
                        Region::Intergenic
                    })
                } else {
                    None
                };

                // Capture the multi-gene sets only when requested and ambiguous,
                // for --soloMultiMappers distribution.
                let gene_multi = if want_multi && want_exon && exon_s.len() > 1 {
                    exon_s.iter().map(|&g| g as u32).collect()
                } else {
                    Vec::new()
                };
                let gene_full_multi = if want_multi && want_body && body_s.len() > 1 {
                    body_s.iter().map(|&g| g as u32).collect()
                } else {
                    Vec::new()
                };

                ReadClass {
                    gene: if want_exon {
                        assignment_of(&exon_s)
                    } else {
                        GeneAssignment::Unmapped
                    },
                    gene_full: if want_body {
                        assignment_of(&body_s)
                    } else {
                        GeneAssignment::Unmapped
                    },
                    region,
                    antisense: body_anti_any && body_s.is_empty(),
                    gene_multi,
                    gene_full_multi,
                }
            })
        })
    })
}

/// Assign a single-end (cDNA) read to a gene from its alignment set, using the
/// `Gene` (exonic) or `GeneFull` (gene-body, intron-inclusive) overlap basis.
/// Thin wrapper over [`classify_read`] for the single-feature case (and tests).
pub fn assign_gene_se(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    model: &GeneModel,
    strand: SoloStrand,
    feature: SoloFeature,
) -> GeneAssignment {
    let want_exon = feature == SoloFeature::Gene;
    let class = classify_read(
        transcripts,
        gene_ann,
        model,
        strand,
        want_exon,
        !want_exon,
        false,
    );
    match feature {
        SoloFeature::Gene => class.gene,
        SoloFeature::GeneFull => class.gene_full,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{Exon, Transcript};
    use crate::genome::Genome;
    use crate::junction::gtf::GtfRecord;
    use std::collections::HashMap;

    fn genome() -> Genome {
        Genome {
            transform_blocks: None,
            sequence: vec![0u8; 2000].into(),
            n_genome: 2000,
            n_genome_real: 2000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        }
    }

    fn gtf_exon(start: u64, end: u64, strand: char, gene: &str) -> GtfRecord {
        let mut attrs = HashMap::new();
        attrs.insert("gene_id".to_string(), gene.to_string());
        attrs.insert("transcript_id".to_string(), format!("{gene}_t1"));
        GtfRecord {
            seqname: "chr1".to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes: attrs,
        }
    }

    /// G1 (+) at 100-200, G2 (-) at 300-400.
    fn annotation_exons() -> Vec<GtfRecord> {
        vec![gtf_exon(100, 200, '+', "G1"), gtf_exon(300, 400, '-', "G2")]
    }

    fn annotation() -> GeneAnnotation {
        GeneAnnotation::from_gtf_exons(&annotation_exons(), &genome())
    }

    /// Transcript models built from the same GTF records as the annotation.
    fn model_of(exons: &[GtfRecord], ann: &GeneAnnotation) -> GeneModel {
        GeneModel::new(
            TranscriptomeIndex::from_gtf_exons(exons, &genome()).unwrap(),
            ann,
        )
    }

    fn read_at(start: u64, end: u64, is_reverse: bool) -> Transcript {
        Transcript {
            output_rank: 0,
            chr_idx: 0,
            genome_start: start,
            genome_end: end,
            is_reverse,
            exons: vec![Exon {
                genome_start: start,
                genome_end: end,
                read_start: 0,
                read_end: (end - start) as usize,
                i_frag: 0,
            }],
            cigar: Vec::new(),
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            junction_strand: vec![],
        }
    }

    #[test]
    fn unmapped_when_no_transcripts() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        assert_eq!(
            assign_gene_se(&[], &ann, &model, SoloStrand::Forward, SoloFeature::Gene),
            GeneAssignment::Unmapped
        );
    }

    #[test]
    fn forward_sense_assigns_g1() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        // Read on + strand overlapping G1 (a + gene).
        let tr = read_at(120, 180, false);
        match assign_gene_se(&[tr], &ann, &model, SoloStrand::Forward, SoloFeature::Gene) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1, got {other:?}"),
        }
    }

    #[test]
    fn forward_antisense_is_no_feature() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        // Read on - strand overlapping G1 (+): wrong strand under Forward.
        let tr = read_at(120, 180, true);
        assert_eq!(
            assign_gene_se(&[tr], &ann, &model, SoloStrand::Forward, SoloFeature::Gene),
            GeneAssignment::NoFeature
        );
    }

    #[test]
    fn reverse_strand_picks_antisense() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        // Read on - strand overlapping G1 (+): kept under Reverse.
        let tr = read_at(120, 180, true);
        match assign_gene_se(&[tr], &ann, &model, SoloStrand::Reverse, SoloFeature::Gene) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1 under Reverse, got {other:?}"),
        }
    }

    #[test]
    fn no_overlap_is_no_feature() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        let tr = read_at(500, 600, false);
        assert_eq!(
            assign_gene_se(
                &[tr],
                &ann,
                &model,
                SoloStrand::Unstranded,
                SoloFeature::Gene
            ),
            GeneAssignment::NoFeature
        );
    }

    #[test]
    fn multilocus_same_gene_is_unique() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        // Two loci both inside G1 → still gene-unique.
        let a = read_at(110, 150, false);
        let b = read_at(150, 190, false);
        match assign_gene_se(
            &[a, b],
            &ann,
            &model,
            SoloStrand::Forward,
            SoloFeature::Gene,
        ) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1, got {other:?}"),
        }
    }

    #[test]
    fn two_genes_unstranded_is_ambiguous() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        // One locus in G1 (+), one in G2 (-); unstranded sees both.
        let a = read_at(120, 180, false);
        let b = read_at(320, 380, true);
        assert_eq!(
            assign_gene_se(
                &[a, b],
                &ann,
                &model,
                SoloStrand::Unstranded,
                SoloFeature::Gene
            ),
            GeneAssignment::Ambiguous
        );
    }

    /// A spliced read: two blocks joined by an `N` (no CIGAR attached, so the
    /// genome gap is read as a junction).
    fn spliced_read(a: (u64, u64), b: (u64, u64)) -> Transcript {
        let mut t = read_at(a.0, b.1, false);
        let la = (a.1 - a.0) as usize;
        t.exons = vec![
            Exon {
                genome_start: a.0,
                genome_end: a.1,
                read_start: 0,
                read_end: la,
                i_frag: 0,
            },
            Exon {
                genome_start: b.0,
                genome_end: b.1,
                read_start: la,
                read_end: la + (b.1 - b.0) as usize,
                i_frag: 0,
            },
        ];
        t.n_junction = 1;
        t
    }

    /// STAR's `Gene` feature counts an alignment only if it is concordant with
    /// one of the gene's transcripts (`alignToTranscript`): a junction must be
    /// one of its introns, and no block may run into an intron. Being inside
    /// the gene's exons is not enough.
    #[test]
    fn gene_needs_a_concordant_transcript() {
        // G3 (+): exons [500,600) and [800,900).
        let g = genome();
        let exons = vec![gtf_exon(501, 600, '+', "G3"), gtf_exon(801, 900, '+', "G3")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        let model = model_of(&exons, &ann);
        let gene = |t: Transcript| {
            assign_gene_se(&[t], &ann, &model, SoloStrand::Forward, SoloFeature::Gene)
        };
        // The annotated junction 600 -> 800.
        assert_eq!(
            gene(spliced_read((550, 600), (800, 850))),
            GeneAssignment::Gene(0)
        );
        // A junction 580 -> 820: both blocks are exonic, but it is not the
        // transcript's intron.
        assert_eq!(
            gene(spliced_read((550, 580), (820, 850))),
            GeneAssignment::NoFeature
        );
        // An unspliced block running from the exon into the intron.
        assert_eq!(gene(read_at(590, 610, false)), GeneAssignment::NoFeature);
        // Wholly in one exon.
        assert_eq!(gene(read_at(820, 880, false)), GeneAssignment::Gene(0));
    }

    /// `fSet` collects the genes of all alignments, `fAlign` each alignment's,
    /// which is what GX/GN are written from.
    #[test]
    fn classify_align_keeps_per_alignment_genes() {
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        let sets = classify_align(
            &[read_at(120, 180, false), read_at(600, 650, false)],
            &model,
            SoloStrand::Forward,
        );
        assert_eq!(sets.genes, vec![0]);
        assert_eq!(sets.per_align, vec![vec![0], vec![]]);
    }

    /// A mate pair is one alignment (`pair_alignment`): concordance asks one
    /// transcript to hold both mates, a deletion inside a mate is not taken
    /// for a junction, and the pair's strand is mate 1's.
    #[test]
    fn classify_align_takes_a_mate_pair_as_one_alignment() {
        use noodles::sam::alignment::record::cigar::{Op, op::Kind};
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        let deleted = |a: (u64, u64), b: (u64, u64), is_reverse: bool| {
            let mut t = read_at(a.0, b.1, is_reverse);
            let la = (a.1 - a.0) as usize;
            t.exons = vec![
                Exon {
                    genome_start: a.0,
                    genome_end: a.1,
                    read_start: 0,
                    read_end: la,
                    i_frag: 0,
                },
                Exon {
                    genome_start: b.0,
                    genome_end: b.1,
                    read_start: la,
                    read_end: la + (b.1 - b.0) as usize,
                    i_frag: 0,
                },
            ];
            t.cigar = vec![
                Op::new(Kind::Match, la),
                Op::new(Kind::Deletion, (b.0 - a.1) as usize),
                Op::new(Kind::Match, (b.1 - b.0) as usize),
            ];
            t
        };
        let cigar_of = |t: &mut Transcript| {
            t.cigar = vec![Op::new(
                Kind::Match,
                (t.genome_end - t.genome_start) as usize,
            )];
        };
        // Forward pair inside G1's exon, mate 2 with a 1-base deletion.
        let mut m1 = read_at(120, 150, false);
        cigar_of(&mut m1);
        let m2 = deleted((155, 170), (171, 190), true);
        let pair = crate::solo::pair_alignment(&m1, &m2);
        assert_eq!(
            pair.exons.iter().map(|e| e.i_frag).collect::<Vec<_>>(),
            [0, 1, 1]
        );
        let sets = classify_align(std::slice::from_ref(&pair), &model, SoloStrand::Forward);
        assert_eq!(sets.per_align, vec![vec![0]]);

        // Reverse pair: mate 2 is the left mate and the pair is antisense to G1.
        let mut r1 = read_at(160, 190, true);
        cigar_of(&mut r1);
        let mut r2 = read_at(110, 140, false);
        cigar_of(&mut r2);
        let pair = crate::solo::pair_alignment(&r1, &r2);
        assert_eq!(pair.exons[0].genome_start, 110, "left mate first");
        assert!(pair.is_reverse);
        let fwd = classify_align(std::slice::from_ref(&pair), &model, SoloStrand::Forward);
        assert_eq!(fwd.genes, Vec::<u32>::new());
        let rev = classify_align(std::slice::from_ref(&pair), &model, SoloStrand::Reverse);
        assert_eq!(rev.genes, vec![0]);

        // One mate past the transcript's end: not concordant.
        let mut o2 = read_at(190, 230, true);
        cigar_of(&mut o2);
        let pair = crate::solo::pair_alignment(&m1, &o2);
        let sets = classify_align(std::slice::from_ref(&pair), &model, SoloStrand::Forward);
        assert_eq!(sets.per_align, vec![Vec::<u32>::new()]);
    }

    /// Velocyto's per-read transcript types: only transcripts on the strand
    /// `--soloStrand` keeps, each with STAR's type bits.
    #[test]
    fn velocyto_types_follow_the_strand_and_the_exons() {
        use crate::quant::splice_status::AlignVsTranscript as A;
        let ann = annotation();
        let model = model_of(&annotation_exons(), &ann);
        // G1 (+) is a single exon [100, 200): a read inside it is concordant.
        let sense = read_at(120, 180, false);
        let types = velocyto_transcript_types(&sense, &model, SoloStrand::Forward);
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].1, A::Concordant.type_bits());
        // The same read on the other strand: kept only when unstranded.
        let anti = read_at(120, 180, true);
        assert_eq!(
            velocyto_transcript_types(&anti, &model, SoloStrand::Forward),
            []
        );
        assert_eq!(
            velocyto_transcript_types(&anti, &model, SoloStrand::Unstranded).len(),
            1
        );
        // Outside every transcript: nothing.
        let far = read_at(600, 650, false);
        assert_eq!(
            velocyto_transcript_types(&far, &model, SoloStrand::Unstranded),
            []
        );
    }

    /// A spliced alignment leaves out a transcript only if one of its blocks
    /// runs into that transcript's intron: STAR does not ask its junctions to
    /// be the transcript's (`Transcriptome_classifyAlign.cpp:103-158`). A mate
    /// pair's gap is not a junction.
    #[test]
    fn velocyto_types_drop_transcripts_whose_intron_a_spliced_block_enters() {
        use crate::quant::splice_status::AlignVsTranscript as A;
        // G3 (+): T1 spliced, exons [500,600) and [800,900); T2 retains the
        // intron, [500,900).
        let g = genome();
        let mut retained = gtf_exon(501, 900, '+', "G3");
        retained
            .attributes
            .insert("transcript_id".to_string(), "G3_t2".to_string());
        let exons = vec![
            gtf_exon(501, 600, '+', "G3"),
            gtf_exon(801, 900, '+', "G3"),
            retained,
        ];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        let model = model_of(&exons, &ann);
        let bits = |t: &Transcript| {
            let mut b: Vec<u8> = velocyto_transcript_types(t, &model, SoloStrand::Forward)
                .into_iter()
                .map(|(_, ty)| ty)
                .collect();
            b.sort_unstable();
            b
        };
        let exon = A::Concordant.type_bits();
        let intron = A::Intron.type_bits();
        // T1's junction, and one neither transcript has: both blocks are
        // exonic in both transcripts, so both are kept.
        assert_eq!(bits(&spliced_read((550, 600), (800, 850))), [exon, exon]);
        assert_eq!(bits(&spliced_read((550, 580), (820, 850))), [exon, exon]);
        // A spliced read whose first block runs into T1's intron: T1 is out.
        assert_eq!(bits(&spliced_read((550, 620), (820, 850))), [exon]);
        // Unspliced, the same block spans T1's exon and intron.
        let span = A::ExonIntronSpan.type_bits();
        let mut both = vec![exon, span];
        both.sort_unstable();
        assert_eq!(bits(&read_at(550, 620, false)), both);
        // In the intron: intronic in T1, exonic in T2.
        let mut both = vec![exon, intron];
        both.sort_unstable();
        assert_eq!(bits(&read_at(650, 700, false)), both);
        // A pair with a mate in each exon: the mates' gap is no junction, so
        // T2 keeps it, and both transcripts hold it.
        use noodles::sam::alignment::record::cigar::{Op, op::Kind};
        let mate = |start: u64, end: u64, is_reverse: bool| {
            let mut t = read_at(start, end, is_reverse);
            t.cigar = vec![Op::new(Kind::Match, (end - start) as usize)];
            t
        };
        let pair = crate::solo::pair_alignment(&mate(550, 590, false), &mate(810, 850, true));
        assert_eq!(bits(&pair), [exon, exon]);
    }

    /// `GeneFull` gene sets: per alignment, the sense genes whose body it
    /// overlaps.
    #[test]
    fn gene_full_sets_are_per_alignment() {
        let ann = annotation();
        let sets = gene_full_sets(
            &[
                read_at(120, 180, false),
                read_at(120, 180, true),
                read_at(600, 650, false),
            ],
            &ann,
            SoloStrand::Forward,
        );
        assert_eq!(sets.genes, vec![0]);
        assert_eq!(sets.per_align, vec![vec![0], vec![], vec![]]);
    }

    #[test]
    fn genefull_counts_intronic_read() {
        // Two-exon gene G3 (+): exons [500,600) and [800,900) → gene body
        // [500,900) with an intron at [600,800).
        let g = genome();
        let exons = vec![gtf_exon(501, 600, '+', "G3"), gtf_exon(801, 900, '+', "G3")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        let model = model_of(&exons, &ann);
        // A read entirely inside the intron overlaps no exon...
        assert_eq!(
            assign_gene_se(
                &[read_at(650, 700, false)],
                &ann,
                &model,
                SoloStrand::Forward,
                SoloFeature::Gene
            ),
            GeneAssignment::NoFeature
        );
        // ...but does overlap the gene body, so GeneFull counts it.
        match assign_gene_se(
            &[read_at(650, 700, false)],
            &ann,
            &model,
            SoloStrand::Forward,
            SoloFeature::GeneFull,
        ) {
            GeneAssignment::Gene(gi) => assert_eq!(ann.gene_ids[gi as usize], "G3"),
            other => panic!("expected G3 under GeneFull, got {other:?}"),
        }
    }

    /// `gx`/`gn` are per alignment, so two alignments of one read report their
    /// own genes; `sF` reports the read-level overlap type and gene count, and
    /// falls back to `(-1, -1)` on an alignment with no gene of its own.
    /// A read inside an intron is intronic-sense, and carries the gene only on
    /// the `GeneFull` basis.
    #[test]
    fn classify_read_regions_and_antisense() {
        // Ga (+): exons [100,200) and [400,500) → body [100,500), intron [200,400).
        let g = genome();
        let exons = vec![gtf_exon(101, 200, '+', "Ga"), gtf_exon(401, 500, '+', "Ga")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        let model = model_of(&exons, &ann);
        let cls = |start, end, rev| {
            classify_read(
                &[read_at(start, end, rev)],
                &ann,
                &model,
                SoloStrand::Forward,
                true,
                true,
                false,
            )
        };

        // In an exon, sense strand → Exonic, not antisense.
        let c = cls(120, 180, false);
        assert_eq!(c.region, Some(Region::Exonic));
        assert!(!c.antisense);
        assert!(matches!(c.gene, GeneAssignment::Gene(_)));

        // Entirely within the intron → Intronic (body but no exon).
        assert_eq!(cls(250, 350, false).region, Some(Region::Intronic));

        // Outside the gene → Intergenic.
        assert_eq!(cls(700, 800, false).region, Some(Region::Intergenic));

        // Exonic position but read on the opposite strand of a (+) gene:
        // positionally Exonic, flagged antisense, no sense gene assignment.
        let c = cls(120, 180, true);
        assert_eq!(c.region, Some(Region::Exonic));
        assert!(c.antisense);
        assert_eq!(c.gene, GeneAssignment::NoFeature);

        // No region computed when only one side requested.
        assert_eq!(
            classify_read(
                &[read_at(120, 180, false)],
                &ann,
                &model,
                SoloStrand::Forward,
                true,
                false,
                false
            )
            .region,
            None
        );
    }
}
