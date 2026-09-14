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

                for tr in transcripts {
                    if want_exon {
                        gene_ann.overlapping_genes_into(tr, &mut raw);
                        for &g in raw.iter() {
                            exon_any = true;
                            // STAR's Gene feature requires the alignment to be
                            // *concordant* with the gene, not merely to overlap an
                            // exon: every aligned block must lie wholly within the
                            // gene's exons. A read whose block extends past an exon
                            // boundary into an intron (no matching splice) is NOT a
                            // Gene read (it may still be GeneFull/intronic). Pure
                            // exon-overlap over-assigns such boundary reads — the
                            // measured rustar-vs-STAR divergence.
                            let concordant = tr
                                .exons
                                .iter()
                                .all(|b| gene_ann.block_is_exonic(g, b.genome_start, b.genome_end));
                            if concordant
                                && strand_keeps(strand, gene_ann.gene_is_reverse[g], tr.is_reverse)
                            {
                                exon_s.push(g);
                            }
                        }
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

/// STAR's `ReadAnnotFeature::overlapTypes` (`ReadAnnotations.h:13`), the first
/// value of the `sF` SAM tag. The read takes the lowest-numbered type any of its
/// alignments reaches.
///
/// `ExonicSense50p` / `Exonic50pAntisense` are STAR's `GeneFull_Ex50pAS` types
/// and never arise here, since that feature is not implemented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum OverlapType {
    None = 0,
    ExonicSense = 1,
    ExonicAntisense = 2,
    ExonicSense50p = 3,
    Exonic50pAntisense = 4,
    IntronicSense = 5,
    IntronicAntisense = 6,
    Intergenic = 7,
}

impl OverlapType {
    /// The sense-strand types, for which STAR reports `sF:B:i,-1,-1` on an
    /// alignment that carries no gene of its own
    /// (`ReadAlign_alignBAM.cpp:441-447`).
    pub fn is_sense(self) -> bool {
        matches!(
            self,
            Self::ExonicSense | Self::ExonicSense50p | Self::IntronicSense
        )
    }
}

/// The genes each individual alignment of a read belongs to, for the `gx`/`gn`
/// SAM tags (STAR's `readAnnot.annotFeatures[f].fAlign`), plus the read-level
/// gene set size and overlap type behind `sF`.
#[derive(Debug, Clone, Default)]
pub struct AlignGenes {
    /// One sorted gene-index list per alignment, parallel to the transcripts.
    pub per_align: Vec<Vec<u32>>,
    /// Union over alignments (STAR's `fSet`).
    pub gene_set: Vec<u32>,
    /// STAR's `ovType`.
    pub ov_type: i32,
}

/// Per-alignment gene assignment for `gx`/`gn`/`sF`, on the `Gene` (exon,
/// transcript-concordant) or `GeneFull` (gene body) basis.
///
/// Mirrors the per-transcript half of [`classify_read`]: the same concordance
/// and strand rules, kept per alignment instead of unioned over the read.
pub fn align_genes(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    strand: SoloStrand,
    feature: SoloFeature,
) -> AlignGenes {
    let want_exon = feature == SoloFeature::Gene;
    let mut out = AlignGenes {
        per_align: Vec::with_capacity(transcripts.len()),
        ..Default::default()
    };
    if transcripts.is_empty() {
        out.ov_type = OverlapType::None as i32;
        return out;
    }
    let (mut exon_sense, mut exon_anti) = (false, false);
    let (mut body_sense, mut body_anti) = (false, false);
    let mut raw: Vec<usize> = Vec::new();

    for tr in transcripts {
        let mut genes: Vec<u32> = Vec::new();
        // Exon overlap decides the exonic/intronic half of `ovType` even for the
        // GeneFull basis, which is what makes an intronic GeneFull read report
        // `intronic` rather than `exonic`.
        gene_ann.overlapping_genes_into(tr, &mut raw);
        for &g in &raw {
            let concordant = tr
                .exons
                .iter()
                .all(|b| gene_ann.block_is_exonic(g, b.genome_start, b.genome_end));
            if !concordant {
                continue;
            }
            if strand_keeps(strand, gene_ann.gene_is_reverse[g], tr.is_reverse) {
                exon_sense = true;
                if want_exon {
                    genes.push(g as u32);
                }
            } else {
                exon_anti = true;
            }
        }
        gene_ann.overlapping_genes_full_into(tr, &mut raw);
        for &g in &raw {
            if strand_keeps(strand, gene_ann.gene_is_reverse[g], tr.is_reverse) {
                body_sense = true;
                if !want_exon {
                    genes.push(g as u32);
                }
            } else {
                body_anti = true;
            }
        }
        genes.sort_unstable();
        genes.dedup();
        out.gene_set.extend_from_slice(&genes);
        out.per_align.push(genes);
    }
    out.gene_set.sort_unstable();
    out.gene_set.dedup();

    // Lower types win, as in STAR's `otFinal` scan.
    out.ov_type = if exon_sense {
        OverlapType::ExonicSense
    } else if exon_anti {
        OverlapType::ExonicAntisense
    } else if body_sense {
        OverlapType::IntronicSense
    } else if body_anti {
        OverlapType::IntronicAntisense
    } else {
        OverlapType::Intergenic
    } as i32;
    out
}

/// Assign a single-end (cDNA) read to a gene from its alignment set, using the
/// `Gene` (exonic) or `GeneFull` (gene-body, intron-inclusive) overlap basis.
/// Thin wrapper over [`classify_read`] for the single-feature case (and tests).
pub fn assign_gene_se(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    strand: SoloStrand,
    feature: SoloFeature,
) -> GeneAssignment {
    let want_exon = feature == SoloFeature::Gene;
    let class = classify_read(transcripts, gene_ann, strand, want_exon, !want_exon, false);
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
    fn annotation() -> GeneAnnotation {
        let exons = vec![gtf_exon(100, 200, '+', "G1"), gtf_exon(300, 400, '-', "G2")];
        GeneAnnotation::from_gtf_exons(&exons, &genome())
    }

    fn read_at(start: u64, end: u64, is_reverse: bool) -> Transcript {
        Transcript {
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
        }
    }

    #[test]
    fn unmapped_when_no_transcripts() {
        let ann = annotation();
        assert_eq!(
            assign_gene_se(&[], &ann, SoloStrand::Forward, SoloFeature::Gene),
            GeneAssignment::Unmapped
        );
    }

    #[test]
    fn forward_sense_assigns_g1() {
        let ann = annotation();
        // Read on + strand overlapping G1 (a + gene).
        let tr = read_at(120, 180, false);
        match assign_gene_se(&[tr], &ann, SoloStrand::Forward, SoloFeature::Gene) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1, got {other:?}"),
        }
    }

    #[test]
    fn forward_antisense_is_no_feature() {
        let ann = annotation();
        // Read on - strand overlapping G1 (+): wrong strand under Forward.
        let tr = read_at(120, 180, true);
        assert_eq!(
            assign_gene_se(&[tr], &ann, SoloStrand::Forward, SoloFeature::Gene),
            GeneAssignment::NoFeature
        );
    }

    #[test]
    fn reverse_strand_picks_antisense() {
        let ann = annotation();
        // Read on - strand overlapping G1 (+): kept under Reverse.
        let tr = read_at(120, 180, true);
        match assign_gene_se(&[tr], &ann, SoloStrand::Reverse, SoloFeature::Gene) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1 under Reverse, got {other:?}"),
        }
    }

    #[test]
    fn no_overlap_is_no_feature() {
        let ann = annotation();
        let tr = read_at(500, 600, false);
        assert_eq!(
            assign_gene_se(&[tr], &ann, SoloStrand::Unstranded, SoloFeature::Gene),
            GeneAssignment::NoFeature
        );
    }

    #[test]
    fn multilocus_same_gene_is_unique() {
        let ann = annotation();
        // Two loci both inside G1 → still gene-unique.
        let a = read_at(110, 150, false);
        let b = read_at(150, 190, false);
        match assign_gene_se(&[a, b], &ann, SoloStrand::Forward, SoloFeature::Gene) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1, got {other:?}"),
        }
    }

    #[test]
    fn two_genes_unstranded_is_ambiguous() {
        let ann = annotation();
        // One locus in G1 (+), one in G2 (-); unstranded sees both.
        let a = read_at(120, 180, false);
        let b = read_at(320, 380, true);
        assert_eq!(
            assign_gene_se(&[a, b], &ann, SoloStrand::Unstranded, SoloFeature::Gene),
            GeneAssignment::Ambiguous
        );
    }

    #[test]
    fn genefull_counts_intronic_read() {
        // Two-exon gene G3 (+): exons [500,600) and [800,900) → gene body
        // [500,900) with an intron at [600,800).
        let g = genome();
        let exons = vec![gtf_exon(501, 600, '+', "G3"), gtf_exon(801, 900, '+', "G3")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        // A read entirely inside the intron overlaps no exon...
        assert_eq!(
            assign_gene_se(
                &[read_at(650, 700, false)],
                &ann,
                SoloStrand::Forward,
                SoloFeature::Gene
            ),
            GeneAssignment::NoFeature
        );
        // ...but does overlap the gene body, so GeneFull counts it.
        match assign_gene_se(
            &[read_at(650, 700, false)],
            &ann,
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
    #[test]
    fn align_genes_are_per_alignment_with_a_read_level_overlap_type() {
        // Ga (+) exons [100,200); Gb (+) exons [400,500). Intergenic elsewhere.
        let g = genome();
        let exons = vec![gtf_exon(101, 200, '+', "Ga"), gtf_exon(401, 500, '+', "Gb")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);

        let hits = vec![read_at(110, 150, false), read_at(410, 450, false)];
        let out = align_genes(&hits, &ann, SoloStrand::Forward, SoloFeature::Gene);
        assert_eq!(out.per_align, vec![vec![0], vec![1]]);
        assert_eq!(out.gene_set, vec![0, 1]);
        assert_eq!(out.ov_type, OverlapType::ExonicSense as i32);

        // Second alignment intergenic: it keeps no gene, but the read is still
        // exonic-sense from the first.
        let mixed = vec![read_at(110, 150, false), read_at(900, 950, false)];
        let out = align_genes(&mixed, &ann, SoloStrand::Forward, SoloFeature::Gene);
        assert_eq!(out.per_align, vec![vec![0], Vec::<u32>::new()]);
        assert_eq!(out.ov_type, OverlapType::ExonicSense as i32);
        assert!(OverlapType::ExonicSense.is_sense());

        // Antisense read over the same exon: exonicAS, no gene kept.
        let anti = vec![read_at(110, 150, true)];
        let out = align_genes(&anti, &ann, SoloStrand::Forward, SoloFeature::Gene);
        assert_eq!(out.per_align, vec![Vec::<u32>::new()]);
        assert_eq!(out.ov_type, OverlapType::ExonicAntisense as i32);
        assert!(!OverlapType::ExonicAntisense.is_sense());

        // No gene anywhere near: intergenic.
        let far = vec![read_at(900, 950, false)];
        let out = align_genes(&far, &ann, SoloStrand::Forward, SoloFeature::Gene);
        assert_eq!(out.ov_type, OverlapType::Intergenic as i32);
        assert!(out.gene_set.is_empty());
    }

    /// A read inside an intron is intronic-sense, and carries the gene only on
    /// the `GeneFull` basis.
    #[test]
    fn intronic_reads_report_the_intronic_overlap_type() {
        // One gene, two exons: body [100,500), intron [200,400).
        let g = genome();
        let exons = vec![gtf_exon(101, 200, '+', "Ga"), gtf_exon(401, 500, '+', "Ga")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        let intronic = vec![read_at(250, 300, false)];

        let full = align_genes(&intronic, &ann, SoloStrand::Forward, SoloFeature::GeneFull);
        assert_eq!(full.per_align, vec![vec![0]]);
        assert_eq!(full.ov_type, OverlapType::IntronicSense as i32);

        let gene = align_genes(&intronic, &ann, SoloStrand::Forward, SoloFeature::Gene);
        assert_eq!(gene.per_align, vec![Vec::<u32>::new()]);
        assert_eq!(gene.ov_type, OverlapType::IntronicSense as i32);
    }

    #[test]
    fn classify_read_regions_and_antisense() {
        // Ga (+): exons [100,200) and [400,500) → body [100,500), intron [200,400).
        let g = genome();
        let exons = vec![gtf_exon(101, 200, '+', "Ga"), gtf_exon(401, 500, '+', "Ga")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &g);
        let cls = |start, end, rev| {
            classify_read(
                &[read_at(start, end, rev)],
                &ann,
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
