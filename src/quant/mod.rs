/// Gene-level read quantification (`--quantMode GeneCounts`).
///
/// Implements HTSeq-style "union" counting: a read counts toward a gene if any
/// of its exons overlap any exon of that gene.  Output: `ReadsPerGene.out.tab`
/// with four columns — gene_id, unstranded, strand1 (same strand as gene),
/// strand2 (opposite strand).
///
/// Submodules:
/// - `transcriptome` — transcript-level alignment projection for
///   `--quantMode TranscriptomeSAM` (Salmon / RSEM input).
pub mod transcriptome;

use std::io::Write as _;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::align::read_align::PairedAlignment;
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::genome::Genome;
use crate::junction::gtf::GtfRecord;

// ---------------------------------------------------------------------------
// GeneAnnotation — interval index for overlap queries
// ---------------------------------------------------------------------------

/// Per-gene annotation built from GTF exon records.
pub struct GeneAnnotation {
    /// gene_id strings in GTF-file order (index = gene_idx).
    pub gene_ids: Vec<String>,
    /// Strand per gene: true = reverse/minus strand.
    pub gene_is_reverse: Vec<bool>,
    /// Per-chromosome exon interval list, sorted by (start, end).
    /// Each entry: (start_0based_incl, end_0based_excl, gene_idx).
    pub chr_exons: Vec<Vec<(u64, u64, usize)>>,
}

impl GeneAnnotation {
    /// Build from GTF exon records using `gene_tag` as the gene grouping attribute.
    ///
    /// `gene_tag` is STAR's `sjdbGTFtagExonParentGene` (default `"gene_id"`).
    pub fn from_gtf_exons_configured(exons: &[GtfRecord], genome: &Genome, gene_tag: &str) -> Self {
        let n_chrs = genome.n_chr_real;
        let mut gene_ids: Vec<String> = Vec::new();
        let mut gene_is_reverse: Vec<bool> = Vec::new();
        let mut gene_id_to_idx: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        let mut chr_exons: Vec<Vec<(u64, u64, usize)>> = vec![Vec::new(); n_chrs];

        for exon in exons {
            let gene_id = match exon.attributes.get(gene_tag) {
                Some(id) => id.clone(),
                None => continue,
            };

            let gene_idx = if let Some(&idx) = gene_id_to_idx.get(&gene_id) {
                idx
            } else {
                let idx = gene_ids.len();
                gene_id_to_idx.insert(gene_id.clone(), idx);
                let is_rev = exon.strand == '-';
                gene_is_reverse.push(is_rev);
                gene_ids.push(gene_id);
                idx
            };

            let chr_idx = match genome.chr_name.iter().position(|n| n == &exon.seqname) {
                Some(i) => i,
                None => continue,
            };
            if chr_idx >= n_chrs {
                continue;
            }

            // GTF is 1-based inclusive; convert to absolute 0-based [start, end)
            // matching Transcript.exon.genome_start / genome_end (concatenated genome coords).
            let chr_offset = genome.chr_start[chr_idx];
            let start = chr_offset + exon.start.saturating_sub(1);
            let end = chr_offset + exon.end;

            chr_exons[chr_idx].push((start, end, gene_idx));
        }

        for exons in &mut chr_exons {
            exons.sort_unstable_by_key(|&(s, e, _)| (s, e));
            exons.dedup();
        }

        GeneAnnotation {
            gene_ids,
            gene_is_reverse,
            chr_exons,
        }
    }

    /// Build from GTF exon records using default `"gene_id"` attribute (backward-compatible).
    pub fn from_gtf_exons(exons: &[GtfRecord], genome: &Genome) -> Self {
        Self::from_gtf_exons_configured(exons, genome, "gene_id")
    }

    pub fn n_genes(&self) -> usize {
        self.gene_ids.len()
    }

    /// Return indices of all genes whose exons overlap any exon of `transcript`.
    /// Result is sorted and deduplicated.
    pub fn overlapping_genes(&self, transcript: &Transcript) -> Vec<usize> {
        if transcript.chr_idx >= self.chr_exons.len() {
            return Vec::new();
        }
        let chr = &self.chr_exons[transcript.chr_idx];
        if chr.is_empty() {
            return Vec::new();
        }

        let mut genes: Vec<usize> = Vec::new();

        for exon in &transcript.exons {
            let rs = exon.genome_start;
            let re = exon.genome_end;
            if re <= rs {
                continue;
            }
            // All gene exons with start < re are candidates.
            let upper = chr.partition_point(|&(gs, _, _)| gs < re);
            for &(_, ge, gene_idx) in &chr[..upper] {
                // Overlap condition: ge > rs (start already guaranteed < re by upper bound).
                if ge > rs {
                    genes.push(gene_idx);
                }
            }
        }

        genes.sort_unstable();
        genes.dedup();
        genes
    }
}

// ---------------------------------------------------------------------------
// GeneCounts — thread-safe per-gene counters
// ---------------------------------------------------------------------------

/// Atomic counters for `ReadsPerGene.out.tab`.
///
/// STAR runs THREE INDEPENDENT counting passes per read — one per output column:
///   col1 (unstranded): all overlapping genes regardless of strand
///   col2 (strand1):    only genes on the SAME strand as the read
///   col3 (strand2):    only genes on the OPPOSITE strand from the read
///
/// N_noFeature and N_ambiguous are independent per column.
/// N_unmapped and N_multimapping are shared (same value in all columns).
pub struct GeneCounts {
    /// Column 1: unstranded count (any overlap regardless of strand).
    pub unstranded: Vec<AtomicU64>,
    /// Column 2: reads whose strand matches the gene strand.
    pub strand1: Vec<AtomicU64>,
    /// Column 3: reads whose strand is opposite to the gene strand.
    pub strand2: Vec<AtomicU64>,
    /// N_unmapped — same for all 3 columns (includes too-many-loci reads).
    pub n_unmapped: AtomicU64,
    /// N_multimapping — same for all 3 columns.
    pub n_multimapping: AtomicU64,
    /// N_noFeature per column (independent).
    pub n_no_feature: [AtomicU64; 3],
    /// N_ambiguous per column (independent).
    pub n_ambiguous: [AtomicU64; 3],
}

impl GeneCounts {
    pub fn new(n_genes: usize) -> Self {
        GeneCounts {
            unstranded: (0..n_genes).map(|_| AtomicU64::new(0)).collect(),
            strand1: (0..n_genes).map(|_| AtomicU64::new(0)).collect(),
            strand2: (0..n_genes).map(|_| AtomicU64::new(0)).collect(),
            n_unmapped: AtomicU64::new(0),
            n_multimapping: AtomicU64::new(0),
            n_no_feature: [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)],
            n_ambiguous: [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)],
        }
    }

    /// Count a single-end read.
    ///
    /// - `transcripts.is_empty() && n_for_mapq > 0` → too-many-loci → N_unmapped (STAR behavior)
    /// - `transcripts.is_empty()` → genuinely unmapped → N_unmapped
    /// - `transcripts.len() > 1` → multi-mapper → N_multimapping
    /// - `transcripts.len() == 1` → unique; run 3 independent overlap passes
    pub fn count_se_read(
        &self,
        transcripts: &[Transcript],
        _n_for_mapq: usize,
        gene_ann: &GeneAnnotation,
    ) {
        if transcripts.is_empty() {
            // Both genuine unmapped AND too-many-loci (n_for_mapq > 0) go to N_unmapped.
            // STAR counts too-many-loci in N_unmapped for GeneCounts (not N_multimapping).
            self.n_unmapped.fetch_add(1, Ordering::Relaxed);
            return;
        }
        if transcripts.len() > 1 {
            self.n_multimapping.fetch_add(1, Ordering::Relaxed);
            return;
        }

        self.count_unique_hit(&transcripts[0], transcripts[0].is_reverse, gene_ann);
    }

    /// Count a paired-end fragment (both mates uniquely mapped).
    ///
    /// - `unmapped` or empty with no half-mapped → N_unmapped
    /// - `both_mapped.len() > 1` → N_multimapping
    /// - `both_mapped.len() == 1` → unique; union of exon overlaps from both mates,
    ///   strand determined by mate1
    pub fn count_pe_read(
        &self,
        both_mapped: &[&PairedAlignment],
        unmapped: bool,
        half_mapped: bool,
        gene_ann: &GeneAnnotation,
    ) {
        if unmapped || (both_mapped.is_empty() && half_mapped) {
            self.n_unmapped.fetch_add(1, Ordering::Relaxed);
            return;
        }
        if both_mapped.len() > 1 {
            self.n_multimapping.fetch_add(1, Ordering::Relaxed);
            return;
        }
        if both_mapped.is_empty() {
            self.n_unmapped.fetch_add(1, Ordering::Relaxed);
            return;
        }

        let pair: &PairedAlignment = both_mapped[0];
        // Union of all genes overlapping either mate.
        let mut all_genes = gene_ann.overlapping_genes(&pair.mate1_transcript);
        let genes2 = gene_ann.overlapping_genes(&pair.mate2_transcript);
        all_genes.extend_from_slice(&genes2);
        all_genes.sort_unstable();
        all_genes.dedup();

        // Strand is determined by mate1.
        let read_is_rev = pair.mate1_transcript.is_reverse;
        self.apply_three_columns(&all_genes, read_is_rev, gene_ann);
    }

    /// Write `ReadsPerGene.out.tab` in STAR's format.
    pub fn write_output(&self, path: &Path, gene_ann: &GeneAnnotation) -> Result<(), Error> {
        let mut file = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;

        macro_rules! wl {
            ($($arg:tt)*) => {
                writeln!(file, $($arg)*).map_err(|e| Error::io(e, path))?;
            };
        }

        let nu = self.n_unmapped.load(Ordering::Relaxed);
        let nm = self.n_multimapping.load(Ordering::Relaxed);
        wl!("N_unmapped\t{nu}\t{nu}\t{nu}");
        wl!("N_multimapping\t{nm}\t{nm}\t{nm}");
        wl!(
            "N_noFeature\t{}\t{}\t{}",
            self.n_no_feature[0].load(Ordering::Relaxed),
            self.n_no_feature[1].load(Ordering::Relaxed),
            self.n_no_feature[2].load(Ordering::Relaxed),
        );
        wl!(
            "N_ambiguous\t{}\t{}\t{}",
            self.n_ambiguous[0].load(Ordering::Relaxed),
            self.n_ambiguous[1].load(Ordering::Relaxed),
            self.n_ambiguous[2].load(Ordering::Relaxed),
        );

        for (i, gene_id) in gene_ann.gene_ids.iter().enumerate() {
            wl!(
                "{}\t{}\t{}\t{}",
                gene_id,
                self.unstranded[i].load(Ordering::Relaxed),
                self.strand1[i].load(Ordering::Relaxed),
                self.strand2[i].load(Ordering::Relaxed),
            );
        }

        Ok(())
    }

    /// Run 3 independent counting passes on a pre-computed set of overlapping genes.
    ///
    /// `all_genes`: sorted+deduplicated gene indices overlapping the read (any strand).
    /// `read_is_rev`: strand of the read (mate1 for PE).
    fn apply_three_columns(
        &self,
        all_genes: &[usize],
        read_is_rev: bool,
        gene_ann: &GeneAnnotation,
    ) {
        // Col 1 — unstranded: use all_genes directly.
        self.apply_column(all_genes, 0, &self.unstranded);

        // Col 2 — strand1: genes whose strand matches the read.
        let same: Vec<usize> = all_genes
            .iter()
            .copied()
            .filter(|&g| gene_ann.gene_is_reverse[g] == read_is_rev)
            .collect();
        self.apply_column(&same, 1, &self.strand1);

        // Col 3 — strand2: genes whose strand is opposite to the read.
        let opp: Vec<usize> = all_genes
            .iter()
            .copied()
            .filter(|&g| gene_ann.gene_is_reverse[g] != read_is_rev)
            .collect();
        self.apply_column(&opp, 2, &self.strand2);
    }

    fn apply_column(&self, genes: &[usize], col: usize, per_gene: &[AtomicU64]) {
        match genes.len() {
            0 => {
                self.n_no_feature[col].fetch_add(1, Ordering::Relaxed);
            }
            1 => {
                per_gene[genes[0]].fetch_add(1, Ordering::Relaxed);
            }
            _ => {
                self.n_ambiguous[col].fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Count a unique single-transcript hit through all 3 columns.
    fn count_unique_hit(
        &self,
        transcript: &Transcript,
        read_is_reverse: bool,
        gene_ann: &GeneAnnotation,
    ) {
        let all_genes = gene_ann.overlapping_genes(transcript);
        self.apply_three_columns(&all_genes, read_is_reverse, gene_ann);
    }
}

// ---------------------------------------------------------------------------
// QuantContext — top-level bundle (passed as Arc to alignment loops)
// ---------------------------------------------------------------------------

/// Bundles GeneAnnotation + GeneCounts for cheap Arc sharing across threads.
pub struct QuantContext {
    pub gene_ann: GeneAnnotation,
    pub counts: GeneCounts,
}

impl QuantContext {
    /// Build from a GTF file.  Call once before alignment.
    pub fn build(
        gtf_path: &Path,
        genome: &Genome,
        feature_exon: &str,
        chr_prefix: &str,
        gene_tag: &str,
    ) -> Result<Self, Error> {
        let exons = crate::junction::gtf::parse_gtf_configured(gtf_path, feature_exon, chr_prefix)?;
        let gene_ann = GeneAnnotation::from_gtf_exons_configured(&exons, genome, gene_tag);
        let n = gene_ann.n_genes();
        log::info!("quantMode GeneCounts: {} genes loaded from GTF", n);
        let counts = GeneCounts::new(n);
        Ok(QuantContext { gene_ann, counts })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{Exon, Transcript};
    use crate::genome::Genome;
    use crate::junction::gtf::GtfRecord;

    fn make_genome() -> Genome {
        Genome {
            sequence: vec![0u8; 2000],
            n_genome: 2000,
            n_chr_real: 2,
            chr_start: vec![0, 1000, 2000],
            chr_length: vec![1000, 1000],
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        }
    }

    fn make_gtf_exon(
        seqname: &str,
        start: u64,
        end: u64,
        strand: char,
        gene_id: &str,
    ) -> GtfRecord {
        let mut attrs = std::collections::HashMap::new();
        attrs.insert("gene_id".to_string(), gene_id.to_string());
        attrs.insert("transcript_id".to_string(), "T1".to_string());
        GtfRecord {
            seqname: seqname.to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes: attrs,
        }
    }

    fn make_transcript(chr_idx: usize, gs: u64, ge: u64, is_reverse: bool) -> Transcript {
        Transcript {
            chr_idx,
            genome_start: gs,
            genome_end: ge,
            is_reverse,
            exons: vec![Exon {
                genome_start: gs,
                genome_end: ge,
                read_start: 0,
                read_end: (ge - gs) as usize,
                i_frag: 0,
            }],
            cigar: vec![],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        }
    }

    #[test]
    fn test_gene_annotation_basic() {
        let genome = make_genome();
        let exons = vec![
            make_gtf_exon("chr1", 101, 200, '+', "G1"), // 0-based: [100, 200)
            make_gtf_exon("chr1", 301, 400, '+', "G1"),
            make_gtf_exon("chr1", 501, 600, '-', "G2"),
        ];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);
        assert_eq!(ann.n_genes(), 2);
        assert_eq!(ann.gene_ids[0], "G1");
        assert_eq!(ann.gene_ids[1], "G2");
        assert!(!ann.gene_is_reverse[0]); // G1 is +
        assert!(ann.gene_is_reverse[1]); // G2 is -
        assert_eq!(ann.chr_exons[0].len(), 3);
    }

    #[test]
    fn test_overlapping_genes_exact() {
        let genome = make_genome();
        let exons = vec![
            make_gtf_exon("chr1", 101, 200, '+', "G1"),
            make_gtf_exon("chr1", 501, 600, '+', "G2"),
        ];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);

        // Read fully inside G1's exon
        let t = make_transcript(0, 110, 150, false);
        let genes = ann.overlapping_genes(&t);
        assert_eq!(genes, vec![0]); // G1

        // Read fully inside G2's exon
        let t2 = make_transcript(0, 510, 550, false);
        let genes2 = ann.overlapping_genes(&t2);
        assert_eq!(genes2, vec![1]); // G2
    }

    #[test]
    fn test_overlapping_genes_none() {
        let genome = make_genome();
        let exons = vec![make_gtf_exon("chr1", 101, 200, '+', "G1")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);

        // Read in gap between exons
        let t = make_transcript(0, 250, 290, false);
        assert!(ann.overlapping_genes(&t).is_empty());
    }

    #[test]
    fn test_overlapping_genes_ambiguous() {
        let genome = make_genome();
        let exons = vec![
            make_gtf_exon("chr1", 101, 300, '+', "G1"),
            make_gtf_exon("chr1", 201, 400, '+', "G2"),
        ];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);

        // Read overlaps both G1 and G2
        let t = make_transcript(0, 220, 280, false);
        let genes = ann.overlapping_genes(&t);
        assert_eq!(genes.len(), 2);
    }

    #[test]
    fn test_gene_counts_unique() {
        let genome = make_genome();
        let exons = vec![make_gtf_exon("chr1", 101, 200, '+', "G1")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);
        let counts = GeneCounts::new(ann.n_genes());

        // chr_start[0]=0, so gene exon is [100,200) absolute. Transcript at [110,150).
        let t = make_transcript(0, 110, 150, false); // forward, same as G1
        counts.count_se_read(&[t], 1, &ann);

        assert_eq!(counts.unstranded[0].load(Ordering::Relaxed), 1);
        assert_eq!(counts.strand1[0].load(Ordering::Relaxed), 1); // same strand
        assert_eq!(counts.strand2[0].load(Ordering::Relaxed), 0); // opposite: no gene
        assert_eq!(counts.n_unmapped.load(Ordering::Relaxed), 0);
        // col3: gene G1 is on +, read is on + → no opposite-strand genes → N_noFeature[2]++
        assert_eq!(counts.n_no_feature[2].load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_gene_counts_multimapper() {
        let genome = make_genome();
        let exons = vec![make_gtf_exon("chr1", 101, 200, '+', "G1")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);
        let counts = GeneCounts::new(ann.n_genes());

        let t1 = make_transcript(0, 110, 150, false);
        let t2 = make_transcript(0, 110, 150, false);
        counts.count_se_read(&[t1, t2], 2, &ann);

        assert_eq!(counts.n_multimapping.load(Ordering::Relaxed), 1);
        assert_eq!(counts.unstranded[0].load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_gene_counts_unmapped() {
        let genome = make_genome();
        let ann = GeneAnnotation::from_gtf_exons(&[], &genome);
        let counts = GeneCounts::new(0);

        counts.count_se_read(&[], 0, &ann);
        assert_eq!(counts.n_unmapped.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_gene_counts_too_many_loci() {
        let genome = make_genome();
        let ann = GeneAnnotation::from_gtf_exons(&[], &genome);
        let counts = GeneCounts::new(0);

        // too-many-loci: transcripts empty but n_for_mapq > 0 → N_unmapped (STAR behavior)
        counts.count_se_read(&[], 5, &ann);
        assert_eq!(counts.n_unmapped.load(Ordering::Relaxed), 1);
        assert_eq!(counts.n_multimapping.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_gene_counts_no_feature() {
        let genome = make_genome();
        let exons = vec![make_gtf_exon("chr1", 101, 200, '+', "G1")];
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);
        let counts = GeneCounts::new(ann.n_genes());

        let t = make_transcript(0, 800, 850, false); // No overlap with G1
        counts.count_se_read(&[t], 1, &ann);

        // All 3 columns: no gene → N_noFeature for each
        assert_eq!(counts.n_no_feature[0].load(Ordering::Relaxed), 1);
        assert_eq!(counts.n_no_feature[1].load(Ordering::Relaxed), 1);
        assert_eq!(counts.n_no_feature[2].load(Ordering::Relaxed), 1);
        assert_eq!(counts.unstranded[0].load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_gene_counts_strand2() {
        let genome = make_genome();
        let exons = vec![make_gtf_exon("chr1", 101, 200, '+', "G1")]; // Gene on + strand
        let ann = GeneAnnotation::from_gtf_exons(&exons, &genome);
        let counts = GeneCounts::new(ann.n_genes());

        let t = make_transcript(0, 110, 150, true); // Read on - strand → opposite to G1
        counts.count_se_read(&[t], 1, &ann);

        // col1 (unstranded): G1 overlaps → count
        assert_eq!(counts.unstranded[0].load(Ordering::Relaxed), 1);
        // col2 (strand1 = same as read = -): G1 is on + → no same-strand gene → N_noFeature[1]
        assert_eq!(counts.strand1[0].load(Ordering::Relaxed), 0);
        assert_eq!(counts.n_no_feature[1].load(Ordering::Relaxed), 1);
        // col3 (strand2 = opposite = +): G1 is on + → G1 counts
        assert_eq!(counts.strand2[0].load(Ordering::Relaxed), 1);
    }
}
