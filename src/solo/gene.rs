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
pub struct GeneModel {
    pub transcripts: TranscriptomeIndex,
    /// Solo gene index per transcript; `u32::MAX` if the gene is not listed.
    tr_gene: Vec<u32>,
}

impl GeneModel {
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

/// RNA-velocity read category (Sullivan et al. 2025 mature/nascent/ambiguous,
/// reported as scVelo's spliced/unspliced/ambiguous).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VelocytoCategory {
    /// Spans an exon–exon junction → processed (mature) mRNA.
    Spliced,
    /// No junction, but a block extends into an intron → nascent mRNA.
    Unspliced,
    /// No junction, all blocks wholly within exons → origin indistinguishable.
    Ambiguous,
}

/// Classify a uniquely-mapped read (assigned to gene `g` by gene-body overlap)
/// into its velocity category from the alignment: a splice in the CIGAR means
/// the read is mature; otherwise an aligned block that leaves the exons (into an
/// intron) means nascent; a wholly-exonic block is ambiguous.
pub fn velocyto_category(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    g: u32,
) -> VelocytoCategory {
    if transcripts.iter().any(|t| t.n_junction > 0) {
        return VelocytoCategory::Spliced;
    }
    let g = g as usize;
    let all_exonic = transcripts.iter().all(|t| {
        t.exons
            .iter()
            .all(|e| gene_ann.block_is_exonic(g, e.genome_start, e.genome_end))
    });
    if all_exonic {
        VelocytoCategory::Ambiguous
    } else {
        VelocytoCategory::Unspliced
    }
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
