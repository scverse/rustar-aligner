//! STARsolo single-cell support (Phase 14).
//!
//! Phase 14.1 covers barcode-read input plumbing: parsing the cell barcode (CB)
//! and unique molecular identifier (UMI) out of the barcode read for
//! `--soloType CB_UMI_Simple` (droplet 10x-style geometry). Whitelist
//! correction (14.2), gene assignment (14.3), UMI deduplication and matrix
//! output (14.4+) build on the structures defined here.
//!
//! The barcode read is the SECOND `--readFilesIn` file (STAR convention:
//! `--readFilesIn cDNA_read barcode_read`). It is never aligned — only parsed.

pub mod count;
pub mod gene;
pub mod whitelist;

pub use count::{UmiDedup, UmiFiltering, write_gene_matrix};
pub use gene::{GeneAssignment, SoloStrand, assign_gene_se};
pub use whitelist::{
    CbCandidate, CbMatch, CbMatchStats, CbMatchType, CbWhitelist, UmiCheck, check_umi, pack_barcode,
};

use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::io::fastq::{EncodedRead, FastqReader, decode_base};
use crate::params::{Parameters, SoloType};
use crate::quant::GeneAnnotation;
use std::path::Path;
use std::sync::Mutex;

/// Fixed-position cell-barcode + UMI geometry for `CB_UMI_Simple`.
///
/// All offsets are stored 0-based (converted from STAR's 1-based
/// `--soloCBstart` / `--soloUMIstart`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SoloBarcodeLayout {
    /// 0-based start of the cell barcode in the barcode read.
    pub cb_start: usize,
    /// Cell-barcode length in bases.
    pub cb_len: usize,
    /// 0-based start of the UMI in the barcode read.
    pub umi_start: usize,
    /// UMI length in bases.
    pub umi_len: usize,
}

impl SoloBarcodeLayout {
    /// Build the layout from CLI parameters, converting 1-based starts to
    /// 0-based offsets.
    pub fn from_params(params: &Parameters) -> Self {
        Self {
            cb_start: (params.solo_cb_start.max(1) - 1) as usize,
            cb_len: params.solo_cb_len as usize,
            umi_start: (params.solo_umi_start.max(1) - 1) as usize,
            umi_len: params.solo_umi_len as usize,
        }
    }

    /// Minimum barcode-read length required to extract both CB and UMI.
    pub fn min_read_len(&self) -> usize {
        (self.cb_start + self.cb_len).max(self.umi_start + self.umi_len)
    }

    /// Extract the CB and UMI from one barcode read. Returns `None` if the
    /// read is shorter than [`Self::min_read_len`] (the read is then treated
    /// as having no valid barcode).
    pub fn extract(&self, barcode_read: &EncodedRead) -> Option<CellBarcode> {
        let seq = &barcode_read.sequence;
        let qual = &barcode_read.quality;
        if seq.len() < self.min_read_len() {
            return None;
        }
        let cb_seq = seq[self.cb_start..self.cb_start + self.cb_len].to_vec();
        let umi_seq = seq[self.umi_start..self.umi_start + self.umi_len].to_vec();
        // Quality vectors track the FASTQ length; guard in case quality is
        // shorter than sequence (malformed record) by clamping.
        let cb_qual = slice_or_empty(qual, self.cb_start, self.cb_len);
        let umi_qual = slice_or_empty(qual, self.umi_start, self.umi_len);
        Some(CellBarcode {
            cb_seq,
            cb_qual,
            umi_seq,
            umi_qual,
        })
    }
}

fn slice_or_empty(data: &[u8], start: usize, len: usize) -> Vec<u8> {
    if start + len <= data.len() {
        data[start..start + len].to_vec()
    } else {
        Vec::new()
    }
}

/// A cell barcode + UMI extracted from one barcode read.
///
/// Sequences are stored in genome encoding (0=A, 1=C, 2=G, 3=T, 4=N) to match
/// the rest of the pipeline; qualities are raw Phred+33 ASCII bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CellBarcode {
    pub cb_seq: Vec<u8>,
    pub cb_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
}

impl CellBarcode {
    /// True if the cell barcode contains an `N` (encoded 4) — such barcodes
    /// cannot match a whitelist exactly.
    pub fn cb_has_n(&self) -> bool {
        self.cb_seq.contains(&4)
    }

    /// True if the UMI contains an `N`. STARsolo discards reads whose UMI has
    /// any ambiguous base.
    pub fn umi_has_n(&self) -> bool {
        self.umi_seq.contains(&4)
    }

    /// Decode the cell barcode to an ASCII `ACGTN` string (for CB SAM tags and
    /// `barcodes.tsv`).
    pub fn cb_string(&self) -> String {
        decode_seq(&self.cb_seq)
    }

    /// Decode the UMI to an ASCII `ACGTN` string (for UB SAM tags).
    pub fn umi_string(&self) -> String {
        decode_seq(&self.umi_seq)
    }
}

fn decode_seq(encoded: &[u8]) -> String {
    encoded.iter().map(|&b| decode_base(b) as char).collect()
}

/// Reads cDNA reads and their paired barcode reads in lockstep from two FASTQ
/// files. The cDNA read flows into the normal alignment path; the barcode read
/// is parsed into a [`CellBarcode`] (or `None` when too short).
pub struct SoloReadReader {
    cdna: FastqReader,
    barcode: FastqReader,
    layout: SoloBarcodeLayout,
}

/// One cDNA read paired with its (optional) extracted barcode.
pub struct SoloRead {
    pub cdna: EncodedRead,
    /// `None` when the barcode read was too short to extract CB+UMI.
    pub barcode: Option<CellBarcode>,
}

impl SoloReadReader {
    /// Open the cDNA and barcode FASTQ files for a solo run.
    pub fn open(
        cdna_path: &Path,
        barcode_path: &Path,
        layout: SoloBarcodeLayout,
        decompress_cmd: Option<&str>,
    ) -> Result<Self, Error> {
        Ok(Self {
            cdna: FastqReader::open(cdna_path, decompress_cmd)?,
            barcode: FastqReader::open(barcode_path, decompress_cmd)?,
            layout,
        })
    }

    /// Fetch the next paired (cDNA, barcode) read. Errors if the two files
    /// have different lengths.
    pub fn next_read(&mut self) -> Result<Option<SoloRead>, Error> {
        let cdna_opt = self.cdna.next_encoded()?;
        let barcode_opt = self.barcode.next_encoded()?;
        match (cdna_opt, barcode_opt) {
            (Some(cdna), Some(bc)) => {
                let barcode = self.layout.extract(&bc);
                Ok(Some(SoloRead { cdna, barcode }))
            }
            (None, None) => Ok(None),
            (Some(_), None) => Err(Error::from(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "solo: cDNA read file has more reads than the barcode read file",
            ))),
            (None, Some(_)) => Err(Error::from(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "solo: barcode read file has more reads than the cDNA read file",
            ))),
        }
    }

    /// Read up to `batch_size` paired reads for parallel processing.
    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<SoloRead>, Error> {
        let mut batch = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            match self.next_read()? {
                Some(read) => batch.push(read),
                None => break,
            }
        }
        Ok(batch)
    }
}

/// Build a [`SoloReadReader`] from parameters, resolving the cDNA/barcode files
/// from `--readFilesIn`. Returns an error if solo is enabled but the read files
/// are missing (validation should have caught this earlier).
pub fn open_reader(params: &Parameters) -> Result<SoloReadReader, Error> {
    debug_assert!(params.solo_type == SoloType::CbUmiSimple);
    let cdna = params.cdna_read_file().ok_or_else(|| {
        Error::from(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "solo: missing cDNA read file",
        ))
    })?;
    let barcode = params.barcode_read_file().ok_or_else(|| {
        Error::from(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "solo: missing barcode read file",
        ))
    })?;
    let layout = SoloBarcodeLayout::from_params(params);
    SoloReadReader::open(cdna, barcode, layout, params.read_files_command.as_deref())
}

// ---------------------------------------------------------------------------
// CellRanger4 adapter clipping (--clipAdapterType CellRanger4)
// ---------------------------------------------------------------------------

/// The 10x template-switch oligo (TSO), clipped from the 5' of the cDNA read
/// under `--clipAdapterType CellRanger4`. Encoded 0=A,1=C,2=G,3=T.
const TSO_SEQ: &[u8] = b"AAGCAGTGGTATCAACGCAGAGTACATGGG";

/// Clip the 10x TSO from the 5' end and trim a 3' polyA tail of the cDNA read,
/// matching `--clipAdapterType CellRanger4`. Operates on encoded bases
/// (0=A..3=T,4=N) with parallel quality bytes. Returns the clipped read.
///
/// Conservative thresholds (full-length TSO match ≤ 3 mismatches at the 5'
/// anchor; trailing polyA run ≥ 8) keep this a no-op on adapter-free reads.
pub fn clip_adapter_cr4(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut start = 0usize;
    let mut end = seq.len();

    // 5' TSO: compare the read prefix against the full TSO; clip on a match.
    if seq.len() >= TSO_SEQ.len() {
        let tso: Vec<u8> = TSO_SEQ
            .iter()
            .map(|&b| crate::io::fastq::encode_base(b))
            .collect();
        let mismatches = seq[..tso.len()]
            .iter()
            .zip(&tso)
            .filter(|(a, b)| a != b)
            .count();
        if mismatches <= 3 {
            start = tso.len();
        }
    }

    // 3' polyA: trim a trailing run of A (encoded 0) of length >= 8.
    let mut run = 0usize;
    while end > start && seq[end - 1] == 0 {
        run += 1;
        end -= 1;
    }
    if run < 8 {
        end += run; // not a real polyA tail; keep those bases
    }

    if start == 0 && end == seq.len() {
        return (seq.to_vec(), qual.to_vec());
    }
    (
        seq[start..end].to_vec(),
        qual.get(start..end.min(qual.len()))
            .map(<[u8]>::to_vec)
            .unwrap_or_default(),
    )
}

// ---------------------------------------------------------------------------
// Solo counting context + per-read processing (Phase 14.3)
// ---------------------------------------------------------------------------

/// A fully-resolved per-read count record: one (cell, UMI, gene) observation.
/// These are collapsed by UMI per (cell, gene) into the count matrix (14.4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SoloCountRecord {
    /// Sorted whitelist index of the cell barcode.
    pub cb: u32,
    /// 2-bit packed UMI.
    pub umi: u64,
    /// Assigned gene index.
    pub gene: u32,
}

/// A read whose cell barcode matched multiple whitelist entries by 1MM
/// (`1MM_multi`). Resolution to a single CB needs the global exact-count table
/// and is deferred to the collation stage (Phase 14.4).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SoloMultiRecord {
    /// Candidate whitelist barcodes + mismatch quality.
    pub candidates: Vec<CbCandidate>,
    pub umi: u64,
    pub gene: u32,
}

/// Thread-safe sink for the records produced during alignment.
#[derive(Default)]
pub struct SoloRecorder {
    pub records: Mutex<Vec<SoloCountRecord>>,
    pub multi_records: Mutex<Vec<SoloMultiRecord>>,
}

impl SoloRecorder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a batch's records (called from the sequential write phase).
    pub fn extend(&self, recs: Vec<SoloCountRecord>, multi: Vec<SoloMultiRecord>) {
        if !recs.is_empty() {
            self.records.lock().unwrap().extend(recs);
        }
        if !multi.is_empty() {
            self.multi_records.lock().unwrap().extend(multi);
        }
    }

    /// Number of fully-resolved count records collected so far.
    pub fn n_records(&self) -> usize {
        self.records.lock().unwrap().len()
    }

    /// Number of deferred multi-CB records collected so far.
    pub fn n_multi_records(&self) -> usize {
        self.multi_records.lock().unwrap().len()
    }
}

/// Everything the alignment loop needs to quantify a solo run, shared as an
/// `Arc` across rayon threads. The gene model is built from `--sjdbGTFfile`;
/// the whitelist and stats are read concurrently (interior atomics).
pub struct SoloContext {
    pub layout: SoloBarcodeLayout,
    pub whitelist: CbWhitelist,
    pub match_type: CbMatchType,
    pub strand: SoloStrand,
    pub gene_ann: GeneAnnotation,
    pub stats: CbMatchStats,
    pub recorder: SoloRecorder,
}

/// What happened to one solo read — drives the produced record(s) and stats.
#[derive(Debug, Default)]
pub struct SoloReadOutcome {
    /// A resolved count record, if the read was fully assignable.
    pub record: Option<SoloCountRecord>,
    /// A deferred multi-CB record, if the CB was an unresolved 1MM_multi.
    pub multi: Option<SoloMultiRecord>,
}

impl SoloContext {
    /// Build the solo context from parameters: load the whitelist and build the
    /// gene model from `--sjdbGTFfile`. Call once before alignment.
    pub fn build(params: &Parameters, genome: &crate::genome::Genome) -> Result<Self, Error> {
        let whitelist = match params.solo_cb_whitelist_path() {
            Some(path) => {
                log::info!(
                    "STARsolo: loading cell-barcode whitelist from {}",
                    path.display()
                );
                let wl = CbWhitelist::load(&path)?;
                log::info!("STARsolo: {} whitelist barcodes loaded", wl.len());
                wl
            }
            None => CbWhitelist::NoWhitelist {
                len: params.solo_cb_len as usize,
            },
        };

        // Gene model from the GTF (validated to be present for Gene/GeneFull).
        let gtf_path = params.sjdb_gtf_file.as_ref().ok_or_else(|| {
            Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "STARsolo Gene feature requires --sjdbGTFfile",
            ))
        })?;
        let exons = crate::junction::gtf::parse_gtf_configured(
            gtf_path,
            &params.sjdb_gtf_feature_exon,
            &params.sjdb_gtf_chr_prefix,
        )?;
        let gene_ann = GeneAnnotation::from_gtf_exons_configured(
            &exons,
            genome,
            &params.sjdb_gtf_tag_exon_parent_gene,
        );
        log::info!(
            "STARsolo: {} genes loaded from {}",
            gene_ann.n_genes(),
            gtf_path.display()
        );

        let strand: SoloStrand = params.solo_strand.parse().map_err(|e: String| {
            Error::from(std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
        })?;

        Ok(Self {
            layout: SoloBarcodeLayout::from_params(params),
            whitelist,
            match_type: params.solo_cb_match_type(),
            strand,
            gene_ann,
            stats: CbMatchStats::new(),
            recorder: SoloRecorder::new(),
        })
    }

    /// Process one solo read: match the cell barcode, validate the UMI, assign
    /// a gene, and (on success) produce a count record. Stats are recorded
    /// here; the returned records are appended to the recorder by the caller.
    pub fn process_read(
        &self,
        cdna_transcripts: &[Transcript],
        barcode: Option<&CellBarcode>,
    ) -> SoloReadOutcome {
        let mut out = SoloReadOutcome::default();

        // No barcode read (too short) → nothing to count.
        let Some(bc) = barcode else {
            return out;
        };

        // Cell-barcode match.
        let cb_match = self
            .whitelist
            .match_cb(&bc.cb_seq, &bc.cb_qual, self.match_type);
        self.stats.record_cb(&cb_match);

        let cb_resolved: Option<u32> = match &cb_match {
            CbMatch::Exact(idx) | CbMatch::Corrected(idx) => Some(*idx),
            CbMatch::Multi(_) => None, // deferred to collation
            CbMatch::NoMatch | CbMatch::NinCb | CbMatch::MultMatchRejected => return out,
        };

        // UMI validity.
        let umi = match check_umi(&bc.umi_seq) {
            UmiCheck::Ok(packed) => {
                self.stats.record_umi(&UmiCheck::Ok(packed));
                packed
            }
            rejected => {
                self.stats.record_umi(&rejected);
                return out;
            }
        };

        // Gene assignment (only counted reads produce records).
        let gene = match assign_gene_se(cdna_transcripts, &self.gene_ann, self.strand) {
            GeneAssignment::Gene(g) => g,
            GeneAssignment::NoFeature | GeneAssignment::Ambiguous | GeneAssignment::Unmapped => {
                return out;
            }
        };

        match (cb_resolved, &cb_match) {
            (Some(cb), _) => {
                out.record = Some(SoloCountRecord { cb, umi, gene });
            }
            (None, CbMatch::Multi(cands)) => {
                out.multi = Some(SoloMultiRecord {
                    candidates: cands.clone(),
                    umi,
                    gene,
                });
            }
            (None, _) => unreachable!("non-multi unresolved CB returned early"),
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::encode_base;

    fn encoded_read(name: &str, seq: &str, qual: &str) -> EncodedRead {
        EncodedRead {
            name: name.to_string(),
            sequence: seq.bytes().map(encode_base).collect(),
            quality: qual.bytes().collect(),
        }
    }

    fn v2_layout() -> SoloBarcodeLayout {
        // 10x v2: CB at 1..16 (16 bp), UMI at 17..26 (10 bp).
        SoloBarcodeLayout {
            cb_start: 0,
            cb_len: 16,
            umi_start: 16,
            umi_len: 10,
        }
    }

    #[test]
    fn layout_from_params_converts_to_zero_based() {
        let params = Parameters::try_parse_from([
            "rustar-aligner",
            "--soloType",
            "CB_UMI_Simple",
            "--readFilesIn",
            "cdna.fq",
            "bc.fq",
            "--sjdbGTFfile",
            "genes.gtf",
            "--soloCBwhitelist",
            "wl.txt",
        ])
        .unwrap();
        let layout = SoloBarcodeLayout::from_params(&params);
        assert_eq!(layout.cb_start, 0);
        assert_eq!(layout.cb_len, 16);
        assert_eq!(layout.umi_start, 16);
        assert_eq!(layout.umi_len, 10);
        assert_eq!(layout.min_read_len(), 26);
    }

    #[test]
    fn extract_v2_barcode() {
        let layout = v2_layout();
        // 16bp CB = AAAAAAAACCCCCCCC, 10bp UMI = GGGGGTTTTT.
        let read = encoded_read(
            "bc1",
            "AAAAAAAACCCCCCCCGGGGGTTTTT",
            "IIIIIIIIIIIIIIIIJJJJJJJJJJ",
        );
        let bc = layout.extract(&read).expect("should extract");
        assert_eq!(bc.cb_string(), "AAAAAAAACCCCCCCC");
        assert_eq!(bc.umi_string(), "GGGGGTTTTT");
        assert_eq!(bc.cb_qual.len(), 16);
        assert_eq!(bc.umi_qual.len(), 10);
        assert!(!bc.cb_has_n());
        assert!(!bc.umi_has_n());
    }

    #[test]
    fn extract_too_short_returns_none() {
        let layout = v2_layout();
        let read = encoded_read("short", "AAAAAAAACCCC", "IIIIIIIIIIII");
        assert!(layout.extract(&read).is_none());
    }

    #[test]
    fn detects_n_in_cb_and_umi() {
        let layout = v2_layout();
        let read = encoded_read(
            "bcN",
            "AAAAAAAANCCCCCCCGGGGGTTTTN",
            "IIIIIIIIIIIIIIIIJJJJJJJJJJ",
        );
        let bc = layout.extract(&read).unwrap();
        assert!(bc.cb_has_n());
        assert!(bc.umi_has_n());
    }

    #[test]
    fn reader_pairs_cdna_and_barcode() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut cdna = NamedTempFile::new().unwrap();
        writeln!(cdna, "@r1\nACGTACGTAC\n+\nIIIIIIIIII").unwrap();
        writeln!(cdna, "@r2\nTTTTGGGGCC\n+\nIIIIIIIIII").unwrap();
        cdna.flush().unwrap();

        let mut bc = NamedTempFile::new().unwrap();
        writeln!(
            bc,
            "@r1\nAAAAAAAACCCCCCCCGGGGGTTTTT\n+\nIIIIIIIIIIIIIIIIJJJJJJJJJJ"
        )
        .unwrap();
        writeln!(
            bc,
            "@r2\nGGGGGGGGTTTTTTTTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIJJJJJJJJJJ"
        )
        .unwrap();
        bc.flush().unwrap();

        let mut reader = SoloReadReader::open(cdna.path(), bc.path(), v2_layout(), None).unwrap();
        let batch = reader.read_batch(10).unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0].cdna.name, "r1");
        assert_eq!(
            batch[0].barcode.as_ref().unwrap().cb_string(),
            "AAAAAAAACCCCCCCC"
        );
        assert_eq!(
            batch[1].barcode.as_ref().unwrap().umi_string(),
            "ACGTACGTAC"
        );
    }

    #[test]
    fn reader_length_mismatch_errors() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut cdna = NamedTempFile::new().unwrap();
        writeln!(cdna, "@r1\nACGT\n+\nIIII").unwrap();
        writeln!(cdna, "@r2\nTTTT\n+\nIIII").unwrap();
        cdna.flush().unwrap();

        let mut bc = NamedTempFile::new().unwrap();
        writeln!(
            bc,
            "@r1\nAAAAAAAACCCCCCCCGGGGGTTTTT\n+\nIIIIIIIIIIIIIIIIJJJJJJJJJJ"
        )
        .unwrap();
        bc.flush().unwrap();

        let mut reader = SoloReadReader::open(cdna.path(), bc.path(), v2_layout(), None).unwrap();
        assert!(reader.read_batch(10).is_err());
    }
}
