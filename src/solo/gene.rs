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

/// Assign a single-end (cDNA) read to a gene from its alignment set.
pub fn assign_gene_se(
    transcripts: &[Transcript],
    gene_ann: &GeneAnnotation,
    strand: SoloStrand,
) -> GeneAssignment {
    if transcripts.is_empty() {
        return GeneAssignment::Unmapped;
    }

    let mut genes: Vec<usize> = Vec::new();
    for tr in transcripts {
        for g in gene_ann.overlapping_genes(tr) {
            if strand_keeps(strand, gene_ann.gene_is_reverse[g], tr.is_reverse) {
                genes.push(g);
            }
        }
    }
    genes.sort_unstable();
    genes.dedup();

    match genes.len() {
        0 => GeneAssignment::NoFeature,
        1 => GeneAssignment::Gene(genes[0] as u32),
        _ => GeneAssignment::Ambiguous,
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
            sequence: vec![0u8; 2000],
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
            read_seq: Vec::new(),
        }
    }

    #[test]
    fn unmapped_when_no_transcripts() {
        let ann = annotation();
        assert_eq!(
            assign_gene_se(&[], &ann, SoloStrand::Forward),
            GeneAssignment::Unmapped
        );
    }

    #[test]
    fn forward_sense_assigns_g1() {
        let ann = annotation();
        // Read on + strand overlapping G1 (a + gene).
        let tr = read_at(120, 180, false);
        match assign_gene_se(&[tr], &ann, SoloStrand::Forward) {
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
            assign_gene_se(&[tr], &ann, SoloStrand::Forward),
            GeneAssignment::NoFeature
        );
    }

    #[test]
    fn reverse_strand_picks_antisense() {
        let ann = annotation();
        // Read on - strand overlapping G1 (+): kept under Reverse.
        let tr = read_at(120, 180, true);
        match assign_gene_se(&[tr], &ann, SoloStrand::Reverse) {
            GeneAssignment::Gene(g) => assert_eq!(ann.gene_ids[g as usize], "G1"),
            other => panic!("expected G1 under Reverse, got {other:?}"),
        }
    }

    #[test]
    fn no_overlap_is_no_feature() {
        let ann = annotation();
        let tr = read_at(500, 600, false);
        assert_eq!(
            assign_gene_se(&[tr], &ann, SoloStrand::Unstranded),
            GeneAssignment::NoFeature
        );
    }

    #[test]
    fn multilocus_same_gene_is_unique() {
        let ann = annotation();
        // Two loci both inside G1 → still gene-unique.
        let a = read_at(110, 150, false);
        let b = read_at(150, 190, false);
        match assign_gene_se(&[a, b], &ann, SoloStrand::Forward) {
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
            assign_gene_se(&[a, b], &ann, SoloStrand::Unstranded),
            GeneAssignment::Ambiguous
        );
    }
}
