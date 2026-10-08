//! Per-gene, per-exonic-base coverage accumulation (`--quantMode GeneCoverage`).
//!
//! **Not a STAR feature.** This is the input half of the DegNorm degradation
//! normalization pipeline (see [`crate::degnorm`]), captured during alignment so
//! that no sorted BAM has to be written, indexed, and re-read to obtain gene
//! coverage curves.
//!
//! Coverage lives in *transcript space*: a gene's merged, sorted exons are
//! concatenated, and position `j` of gene `g` is the `j`-th exonic base of that
//! gene. The DegNorm model fits a rank-one envelope over exactly these
//! coordinates.

use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use crate::align::read_align::PairedAlignment;
use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::quant::GeneAnnotation;

const COV_MAGIC: &[u8; 8] = b"RSDGNCOV";
const COV_VERSION: u32 = 1;

/// Lock-free per-gene coverage in transcript space.
pub struct GeneCoverage {
    /// Prefix sums of per-gene transcript lengths; `len == n_genes + 1`.
    offsets: Vec<u64>,
    /// Flat coverage array, `len == offsets[n_genes]`.
    cov: Vec<AtomicU32>,
    /// Per-gene raw unique read/fragment counts (the `GeneCounts` column-1 rule).
    counts: Vec<AtomicU32>,
    /// Total reads/fragments counted into any gene (library size).
    n_counted: AtomicU64,
}

impl GeneCoverage {
    pub fn new(ann: &GeneAnnotation) -> Self {
        let n = ann.n_genes();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut acc: u64 = 0;
        offsets.push(0);
        for exons in ann.gene_exons.iter().take(n) {
            acc += exons.iter().map(|&(s, e)| e - s).sum::<u64>();
            offsets.push(acc);
        }
        GeneCoverage {
            offsets,
            cov: (0..acc).map(|_| AtomicU32::new(0)).collect(),
            counts: (0..n).map(|_| AtomicU32::new(0)).collect(),
            n_counted: AtomicU64::new(0),
        }
    }

    pub fn n_genes(&self) -> usize {
        self.counts.len()
    }

    pub fn gene_len(&self, g: usize) -> u32 {
        (self.offsets[g + 1] - self.offsets[g]) as u32
    }

    pub fn total_len(&self) -> u64 {
        self.offsets[self.offsets.len() - 1]
    }

    pub fn n_counted(&self) -> u64 {
        self.n_counted.load(Ordering::Relaxed)
    }

    pub fn gene_count(&self, g: usize) -> u32 {
        self.counts[g].load(Ordering::Relaxed)
    }

    /// Increment gene `g`'s coverage over the genomic interval `[start, end)`,
    /// restricted to the gene's merged exons.
    pub fn add_block(&self, g: usize, start: u64, end: u64, ann: &GeneAnnotation) {
        let mut tx_off: u64 = 0;
        for &(es, ee) in &ann.gene_exons[g] {
            if end > es && start < ee {
                let s = start.max(es);
                let e = end.min(ee);
                let base = self.offsets[g] + tx_off + (s - es);
                for i in 0..(e - s) {
                    self.cov[(base + i) as usize].fetch_add(1, Ordering::Relaxed);
                }
            }
            tx_off += ee - es;
        }
    }

    /// Single-end: accumulate coverage for a uniquely mapped read assigned to
    /// exactly one gene.
    ///
    /// Mirrors the `GeneCounts` column-1 rule so counts and coverage always
    /// agree; multimappers and gene-ambiguous reads are skipped, matching
    /// DegNorm's default (unique alignments only).
    pub fn count_se_read(&self, transcripts: &[Transcript], ann: &GeneAnnotation) {
        if transcripts.len() != 1 {
            return;
        }
        let t = &transcripts[0];
        let mut genes = Vec::new();
        ann.overlapping_genes_into(t, &mut genes);
        if genes.len() != 1 {
            return;
        }
        let g = genes[0];
        for ex in &t.exons {
            self.add_block(g, ex.genome_start, ex.genome_end, ann);
        }
        self.counts[g].fetch_add(1, Ordering::Relaxed);
        self.n_counted.fetch_add(1, Ordering::Relaxed);
    }

    /// Paired-end: the fragment is one observation. Mate blocks are merged
    /// before accumulation, so an overlapping pair contributes 1 per base, which
    /// is what DegNorm's paired-read coverage does.
    pub fn count_pe_read(&self, both_mapped: &[&PairedAlignment], ann: &GeneAnnotation) {
        if both_mapped.len() != 1 {
            return;
        }
        let pair = both_mapped[0];
        let mut genes = Vec::new();
        let mut genes2 = Vec::new();
        ann.overlapping_genes_into(&pair.mate1_transcript, &mut genes);
        ann.overlapping_genes_into(&pair.mate2_transcript, &mut genes2);
        genes.extend_from_slice(&genes2);
        genes.sort_unstable();
        genes.dedup();
        if genes.len() != 1 {
            return;
        }
        let g = genes[0];

        let mut blocks: Vec<(u64, u64)> = pair
            .mate1_transcript
            .exons
            .iter()
            .chain(pair.mate2_transcript.exons.iter())
            .map(|e| (e.genome_start, e.genome_end))
            .collect();
        blocks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(blocks.len());
        for (s, e) in blocks {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        for (s, e) in merged {
            self.add_block(g, s, e, ann);
        }
        self.counts[g].fetch_add(1, Ordering::Relaxed);
        self.n_counted.fetch_add(1, Ordering::Relaxed);
    }

    /// Read back one gene's coverage vector.
    pub fn gene_slice(&self, g: usize) -> Vec<u32> {
        (self.offsets[g]..self.offsets[g + 1])
            .map(|i| self.cov[i as usize].load(Ordering::Relaxed))
            .collect()
    }

    /// Write `GeneCoverage.out.bin`: a gzip stream holding a fixed header, a
    /// per-gene table, the gene id block, and the flat coverage array.
    pub fn write_file(
        &self,
        path: &Path,
        ann: &GeneAnnotation,
        sample_id: &str,
        paired: bool,
    ) -> Result<(), Error> {
        let file = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
        let mut w = GzEncoder::new(BufWriter::new(file), Compression::new(6));
        let io = |e: std::io::Error| Error::io(e, path);

        w.write_all(COV_MAGIC).map_err(io)?;
        w.write_all(&COV_VERSION.to_le_bytes()).map_err(io)?;
        w.write_all(&u32::from(paired).to_le_bytes()).map_err(io)?;
        w.write_all(&(self.n_genes() as u64).to_le_bytes())
            .map_err(io)?;
        w.write_all(&self.total_len().to_le_bytes()).map_err(io)?;
        w.write_all(&self.n_counted().to_le_bytes()).map_err(io)?;

        let sid = sample_id.as_bytes();
        w.write_all(&(sid.len() as u16).to_le_bytes()).map_err(io)?;
        w.write_all(sid).map_err(io)?;

        for g in 0..self.n_genes() {
            w.write_all(&self.gene_len(g).to_le_bytes()).map_err(io)?;
            w.write_all(&self.gene_count(g).to_le_bytes()).map_err(io)?;
        }
        for id in ann.gene_ids.iter().take(self.n_genes()) {
            let b = id.as_bytes();
            w.write_all(&(b.len() as u16).to_le_bytes()).map_err(io)?;
            w.write_all(b).map_err(io)?;
        }
        // Buffer the coverage block in chunks to avoid a syscall per base.
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 16);
        for c in &self.cov {
            buf.extend_from_slice(&c.load(Ordering::Relaxed).to_le_bytes());
            if buf.len() >= (1 << 16) {
                w.write_all(&buf).map_err(io)?;
                buf.clear();
            }
        }
        if !buf.is_empty() {
            w.write_all(&buf).map_err(io)?;
        }
        w.finish().map_err(io)?;
        Ok(())
    }
}

/// A loaded `GeneCoverage.out.bin`.
pub struct CoverageFile {
    pub sample_id: String,
    pub paired: bool,
    pub n_counted: u64,
    pub gene_ids: Vec<String>,
    pub gene_lens: Vec<u32>,
    pub counts: Vec<u32>,
    /// Prefix sums of `gene_lens`; `len == n_genes + 1`.
    pub offsets: Vec<u64>,
    pub cov: Vec<u32>,
}

impl CoverageFile {
    pub fn read(path: &Path) -> Result<Self, Error> {
        let file = std::fs::File::open(path).map_err(|e| Error::io(e, path))?;
        let mut r = GzDecoder::new(BufReader::new(file));
        let io = |e: std::io::Error| Error::io(e, path);

        let mut magic = [0u8; 8];
        r.read_exact(&mut magic).map_err(io)?;
        if &magic != COV_MAGIC {
            return Err(Error::Parameter(format!(
                "{} is not a rustar-aligner GeneCoverage file",
                path.display()
            )));
        }
        let mut b4 = [0u8; 4];
        let mut b8 = [0u8; 8];
        let mut b2 = [0u8; 2];

        r.read_exact(&mut b4).map_err(io)?;
        let version = u32::from_le_bytes(b4);
        if version != COV_VERSION {
            return Err(Error::Parameter(format!(
                "{} has GeneCoverage version {version}, expected {COV_VERSION}",
                path.display()
            )));
        }
        r.read_exact(&mut b4).map_err(io)?;
        let paired = u32::from_le_bytes(b4) != 0;
        r.read_exact(&mut b8).map_err(io)?;
        let n_genes = u64::from_le_bytes(b8) as usize;
        r.read_exact(&mut b8).map_err(io)?;
        let total_len = u64::from_le_bytes(b8) as usize;
        r.read_exact(&mut b8).map_err(io)?;
        let n_counted = u64::from_le_bytes(b8);

        r.read_exact(&mut b2).map_err(io)?;
        let mut sid = vec![0u8; u16::from_le_bytes(b2) as usize];
        r.read_exact(&mut sid).map_err(io)?;
        let sample_id = String::from_utf8_lossy(&sid).into_owned();

        let mut gene_lens = Vec::with_capacity(n_genes);
        let mut counts = Vec::with_capacity(n_genes);
        for _ in 0..n_genes {
            r.read_exact(&mut b4).map_err(io)?;
            gene_lens.push(u32::from_le_bytes(b4));
            r.read_exact(&mut b4).map_err(io)?;
            counts.push(u32::from_le_bytes(b4));
        }
        let mut gene_ids = Vec::with_capacity(n_genes);
        for _ in 0..n_genes {
            r.read_exact(&mut b2).map_err(io)?;
            let mut b = vec![0u8; u16::from_le_bytes(b2) as usize];
            r.read_exact(&mut b).map_err(io)?;
            gene_ids.push(String::from_utf8_lossy(&b).into_owned());
        }
        let mut raw = vec![0u8; total_len * 4];
        r.read_exact(&mut raw).map_err(io)?;
        let cov: Vec<u32> = raw
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect();

        let mut offsets = Vec::with_capacity(n_genes + 1);
        let mut acc = 0u64;
        offsets.push(0);
        for &l in &gene_lens {
            acc += u64::from(l);
            offsets.push(acc);
        }

        Ok(CoverageFile {
            sample_id,
            paired,
            n_counted,
            gene_ids,
            gene_lens,
            counts,
            offsets,
            cov,
        })
    }

    pub fn gene(&self, g: usize) -> &[u32] {
        &self.cov[self.offsets[g] as usize..self.offsets[g + 1] as usize]
    }

    pub fn n_genes(&self) -> usize {
        self.gene_ids.len()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{Exon, Transcript};

    /// One two-exon gene on chr1: `[100,110)` and `[200,205)`, transcript length 15.
    ///
    /// Built through the real GTF path so the overlap segment trees are present.
    fn two_exon_ann() -> GeneAnnotation {
        let genome = crate::genome::Genome {
            transform_blocks: None,
            sequence: vec![0u8; 2000].into(),
            n_genome: 2000,
            n_genome_real: 2000,
            n_chr_real: 2,
            chr_start: vec![0, 1000, 2000],
            chr_length: vec![1000, 1000],
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        };
        let exon = |start: u64, end: u64| {
            let mut attrs = std::collections::HashMap::new();
            attrs.insert("gene_id".to_string(), "G1".to_string());
            attrs.insert("transcript_id".to_string(), "T1".to_string());
            crate::junction::gtf::GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start,
                end,
                strand: '+',
                attributes: attrs,
            }
        };
        // GTF is 1-based inclusive: 101..110 -> [100, 110), 201..205 -> [200, 205).
        GeneAnnotation::from_gtf_exons(&[exon(101, 110), exon(201, 205)], &genome)
    }

    fn tr(blocks: &[(u64, u64)]) -> Transcript {
        Transcript {
            chr_idx: 0,
            genome_start: blocks[0].0,
            genome_end: blocks[blocks.len() - 1].1,
            is_reverse: false,
            exons: blocks
                .iter()
                .map(|&(s, e)| Exon {
                    genome_start: s,
                    genome_end: e,
                    read_start: 0,
                    read_end: (e - s) as usize,
                    i_frag: 0,
                })
                .collect(),
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
    fn block_inside_first_exon_maps_to_transcript_prefix() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        assert_eq!(cov.gene_len(0), 15);
        cov.add_block(0, 102, 105, &ann);
        assert_eq!(
            cov.gene_slice(0),
            vec![0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
    }

    #[test]
    fn block_spanning_intron_covers_both_exons_only() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        cov.add_block(0, 105, 202, &ann);
        assert_eq!(
            cov.gene_slice(0),
            vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
        );
    }

    #[test]
    fn block_outside_all_exons_is_ignored() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        cov.add_block(0, 120, 150, &ann);
        assert_eq!(cov.gene_slice(0).iter().sum::<u32>(), 0);
    }

    #[test]
    fn multimapper_is_not_counted() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        cov.count_se_read(&[tr(&[(100, 105)]), tr(&[(100, 105)])], &ann);
        assert_eq!(cov.gene_slice(0).iter().sum::<u32>(), 0);
        assert_eq!(cov.gene_count(0), 0);
        assert_eq!(cov.n_counted(), 0);
    }

    #[test]
    fn unique_read_increments_coverage_and_count() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        cov.count_se_read(&[tr(&[(100, 105)])], &ann);
        assert_eq!(cov.gene_slice(0)[..5], [1, 1, 1, 1, 1]);
        assert_eq!(cov.gene_count(0), 1);
        assert_eq!(cov.n_counted(), 1);
    }

    #[test]
    fn overlapping_mates_are_counted_once_per_base() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        let pair = PairedAlignment {
            mate1_transcript: tr(&[(100, 106)]),
            mate2_transcript: tr(&[(104, 110)]),
            mate1_region: (0, 6),
            mate2_region: (0, 6),
            is_proper_pair: true,
            insert_size: 10,
            combined_wt_score: 0,
        };
        cov.count_pe_read(&[&pair], &ann);
        assert_eq!(cov.gene_slice(0)[..10], [1u32; 10]);
        assert_eq!(cov.gene_count(0), 1);
        assert_eq!(cov.n_counted(), 1);
    }

    #[test]
    fn coverage_file_round_trips() {
        let ann = two_exon_ann();
        let cov = GeneCoverage::new(&ann);
        cov.count_se_read(&[tr(&[(100, 105)])], &ann);

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("GeneCoverage.out.bin");
        cov.write_file(&path, &ann, "sampleA", false).unwrap();

        let f = CoverageFile::read(&path).unwrap();
        assert_eq!(f.sample_id, "sampleA");
        assert!(!f.paired);
        assert_eq!(f.n_counted, 1);
        assert_eq!(f.gene_ids, vec!["G1".to_string()]);
        assert_eq!(f.gene_lens, vec![15]);
        assert_eq!(f.counts, vec![1]);
        assert_eq!(f.gene(0), &cov.gene_slice(0)[..]);
    }

    #[test]
    fn bad_magic_is_an_error_not_a_panic() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("junk.bin");
        std::fs::write(&path, b"not a coverage file").unwrap();
        assert!(CoverageFile::read(&path).is_err());
    }
}
