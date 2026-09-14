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

pub mod cell_reads;
pub mod count;
pub mod gene;
pub mod libcxx_rng;
pub mod sgt;
pub mod smartseq;
pub mod transcript3p;
pub mod whitelist;

pub use count::{UmiDedup, UmiFiltering, write_gene_matrix};
pub use gene::{
    GeneAssignment, Region, SoloFeature, SoloStrand, VelocytoCategory, assign_gene_se,
    classify_read, velocyto_category,
};
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
use std::sync::atomic::{AtomicU64, Ordering};

/// Cell-barcode + UMI read geometry. `Simple` is a single fixed-position CB +
/// UMI (`CB_UMI_Simple`); `Complex` assembles the CB from several fixed-position
/// segments (`CB_UMI_Complex`). All offsets are 0-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SoloBarcodeLayout {
    Simple {
        cb_start: usize,
        cb_len: usize,
        umi_start: usize,
        umi_len: usize,
    },
    /// Multi-segment CB: each `(start, len)` is one segment, concatenated in
    /// order to form the cell barcode; `umi = (start, len)`.
    Complex {
        cb_segments: Vec<(usize, usize)>,
        umi: (usize, usize),
    },
}

/// Parse a `--soloCBposition`/`--soloUMIposition` spec
/// (`startAnchor_startDist_endAnchor_endDist`) into a 0-based `(start, len)`.
/// Only read-start anchoring (`anchor = 0`) is supported.
fn parse_position(spec: &str) -> Result<(usize, usize), Error> {
    let f: Vec<&str> = spec.split('_').collect();
    if f.len() != 4 {
        return Err(invalid_pos(
            spec,
            "expected startAnchor_startDist_endAnchor_endDist",
        ));
    }
    let (sa, sd, ea, ed) = (
        f[0].parse::<i64>().ok(),
        f[1].parse::<i64>().ok(),
        f[2].parse::<i64>().ok(),
        f[3].parse::<i64>().ok(),
    );
    match (sa, sd, ea, ed) {
        (Some(0), Some(sd), Some(0), Some(ed)) if sd >= 0 && ed >= sd => {
            Ok((sd as usize, (ed - sd + 1) as usize))
        }
        (Some(0), _, Some(0), _) => Err(invalid_pos(spec, "end < start")),
        _ => Err(invalid_pos(
            spec,
            "only read-start anchoring (anchor=0) is supported",
        )),
    }
}

fn invalid_pos(spec: &str, why: &str) -> Error {
    Error::from(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!("invalid position spec '{spec}': {why}"),
    ))
}

impl SoloBarcodeLayout {
    /// Build the layout from CLI parameters. `CB_UMI_Complex` parses
    /// `--soloCBposition`/`--soloUMIposition`; otherwise fixed Simple geometry.
    pub fn from_params(params: &Parameters) -> Self {
        if params.solo_type == SoloType::CbUmiComplex && !params.solo_cb_position.is_empty() {
            let cb_segments = params
                .solo_cb_position
                .iter()
                .filter_map(|s| parse_position(s).ok())
                .collect();
            let umi = parse_position(&params.solo_umi_position).unwrap_or((0, 0));
            return Self::Complex { cb_segments, umi };
        }
        Self::Simple {
            cb_start: (params.solo_cb_start.max(1) - 1) as usize,
            cb_len: params.solo_cb_len as usize,
            umi_start: (params.solo_umi_start.max(1) - 1) as usize,
            umi_len: params.solo_umi_len as usize,
        }
    }

    /// Minimum barcode-read length required to extract the CB and UMI.
    pub fn min_read_len(&self) -> usize {
        match self {
            Self::Simple {
                cb_start,
                cb_len,
                umi_start,
                umi_len,
            } => (cb_start + cb_len).max(umi_start + umi_len),
            Self::Complex { cb_segments, umi } => cb_segments
                .iter()
                .map(|&(s, l)| s + l)
                .chain(std::iter::once(umi.0 + umi.1))
                .max()
                .unwrap_or(0),
        }
    }

    /// Extract the CB (concatenating segments for `Complex`) and UMI from one
    /// barcode read. `None` if the read is shorter than [`Self::min_read_len`].
    pub fn extract(&self, barcode_read: &EncodedRead) -> Option<CellBarcode> {
        let seq = &barcode_read.sequence;
        let qual = &barcode_read.quality;
        if seq.len() < self.min_read_len() {
            return None;
        }
        match self {
            Self::Simple {
                cb_start,
                cb_len,
                umi_start,
                umi_len,
            } => Some(CellBarcode {
                cb_seq: seq[*cb_start..cb_start + cb_len].to_vec(),
                cb_qual: slice_or_empty(qual, *cb_start, *cb_len),
                umi_seq: seq[*umi_start..umi_start + umi_len].to_vec(),
                umi_qual: slice_or_empty(qual, *umi_start, *umi_len),
            }),
            Self::Complex { cb_segments, umi } => {
                let mut cb_seq = Vec::new();
                let mut cb_qual = Vec::new();
                for &(s, l) in cb_segments {
                    cb_seq.extend_from_slice(&seq[s..s + l]);
                    cb_qual.extend_from_slice(&slice_or_empty(qual, s, l));
                }
                Some(CellBarcode {
                    cb_seq,
                    cb_qual,
                    umi_seq: seq[umi.0..umi.0 + umi.1].to_vec(),
                    umi_qual: slice_or_empty(qual, umi.0, umi.1),
                })
            }
        }
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

/// Decoded per-read values behind the STARsolo barcode SAM tags, built once per
/// read and borrowed by [`crate::io::sam::add_solo_barcode_tags`].
#[derive(Debug, Default)]
pub struct SoloTagStrings {
    cb_seq: String,
    cb_qual: String,
    umi_seq: String,
    umi_qual: String,
    barcode_seq: String,
    barcode_qual: String,
    cb_corrected: Option<String>,
    cb_match: i32,
}

impl SoloTagStrings {
    /// Decode one read's barcode into tag values.
    ///
    /// `barcode` is `None` when the barcode read was too short to carry a
    /// CB+UMI: STAR pads such a read with `N`s, which always scores `cbMatch=-2`
    /// (Ns in the barcode), so that is what `sM` reports, but there is no
    /// barcode sequence to put in `CR`/`CY`/`UR`/`UY`, and those are left off.
    pub fn build(
        barcode: Option<&CellBarcode>,
        barcode_read: Option<&EncodedRead>,
        tags: Option<SoloBarcodeTags>,
        cb_corrected: Option<String>,
    ) -> Self {
        let mut out = Self {
            cb_match: tags.map_or(-2, |t| t.cb_match),
            cb_corrected,
            ..Default::default()
        };
        if let Some(bc) = barcode {
            out.cb_seq = bc.cb_string();
            out.cb_qual = String::from_utf8_lossy(&bc.cb_qual).into_owned();
            out.umi_seq = bc.umi_string();
            out.umi_qual = String::from_utf8_lossy(&bc.umi_qual).into_owned();
        }
        if let Some(read) = barcode_read {
            out.barcode_seq = decode_seq(&read.sequence);
            out.barcode_qual = String::from_utf8_lossy(&read.quality).into_owned();
        }
        out
    }

    /// Borrowed view for the SAM writer.
    pub fn as_values(&self) -> crate::io::sam::SoloBarcodeTagValues<'_> {
        crate::io::sam::SoloBarcodeTagValues {
            cb_seq: &self.cb_seq,
            cb_qual: &self.cb_qual,
            umi_seq: &self.umi_seq,
            umi_qual: &self.umi_qual,
            cb_corrected: self.cb_corrected.as_deref(),
            cb_match: self.cb_match,
            barcode_seq: &self.barcode_seq,
            barcode_qual: &self.barcode_qual,
        }
    }
}

/// STAR's `--soloBarcodeReadLength` handling for a separate barcode read
/// (`SoloReadBarcode_getCBandUMI.cpp:228-241`).
///
/// The default (`1`) means the barcode read must be exactly CB+UMI long, and any
/// other length is a fatal input error. `0` turns the check off, and a read
/// shorter than CB+UMI is then padded with `N` (quality `H`) so it still parses
/// (and, having Ns, is never counted).
#[derive(Debug, Clone, Copy, Default)]
pub struct BarcodeReadLength {
    expected: Option<usize>,
    cbumi_len: usize,
}

impl BarcodeReadLength {
    pub fn from_params(params: &Parameters) -> Self {
        // The barcode inside a cDNA mate has the mate's length, so STAR turns the
        // check off there (`ParametersSolo.cpp:143`), as it does for the
        // variable-geometry Complex chemistry.
        if params.solo_barcode_on_mate1() || params.solo_type == SoloType::CbUmiComplex {
            return Self::default();
        }
        let cbumi_len = params.solo_cb_len as usize + params.solo_umi_len as usize;
        let expected = match params.solo_barcode_read_length {
            0 => None,
            1 => Some(cbumi_len),
            n if n > 0 => Some(n as usize),
            _ => None,
        };
        Self {
            expected,
            cbumi_len,
        }
    }

    /// Check the barcode read's length, padding it up to CB+UMI when checking is
    /// off and it falls short.
    pub fn apply(&self, read: &mut EncodedRead) -> Result<(), Error> {
        if let Some(expected) = self.expected {
            if read.sequence.len() != expected {
                return Err(Error::from(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "solo: barcode read '{}' is {} bases, not the expected {}. \
                         If the CB+UMI length is not the barcode read length, set \
                         --soloBarcodeReadLength to that length, or 0 to skip the check.",
                        read.name,
                        read.sequence.len(),
                        expected
                    ),
                )));
            }
        } else if read.sequence.len() < self.cbumi_len {
            // 4 = N in the genome encoding; 'H' is STAR's filler quality.
            read.sequence.resize(self.cbumi_len, 4);
            read.quality.resize(self.cbumi_len, b'H');
        }
        Ok(())
    }
}

/// Reads cDNA reads and their paired barcode reads in lockstep from two FASTQ
/// files. The cDNA read flows into the normal alignment path; the barcode read
/// is parsed into a [`CellBarcode`] (or `None` when too short).
pub struct SoloReadReader {
    cdna: FastqReader,
    barcode: FastqReader,
    layout: SoloBarcodeLayout,
    keep_barcode_read: bool,
    barcode_len: BarcodeReadLength,
}

/// One cDNA read paired with its (optional) extracted barcode.
pub struct SoloRead {
    pub cdna: EncodedRead,
    /// `None` when the barcode read was too short to extract CB+UMI.
    pub barcode: Option<CellBarcode>,
    /// The whole barcode read, kept only for the `sS`/`sQ` SAM tags.
    pub barcode_read: Option<EncodedRead>,
}

impl SoloReadReader {
    /// Open the cDNA and barcode FASTQ files for a solo run.
    pub fn open(
        cdna_path: &Path,
        barcode_path: &Path,
        layout: SoloBarcodeLayout,
        decompress_cmd: Option<&str>,
        keep_barcode_read: bool,
        barcode_len: BarcodeReadLength,
    ) -> Result<Self, Error> {
        Ok(Self {
            cdna: FastqReader::open(cdna_path, decompress_cmd)?,
            barcode: FastqReader::open(barcode_path, decompress_cmd)?,
            layout,
            keep_barcode_read,
            barcode_len,
        })
    }

    /// Fetch the next paired (cDNA, barcode) read. Errors if the two files
    /// have different lengths.
    pub fn next_read(&mut self) -> Result<Option<SoloRead>, Error> {
        let cdna_opt = self.cdna.next_encoded()?;
        let barcode_opt = self.barcode.next_encoded()?;
        match (cdna_opt, barcode_opt) {
            (Some(cdna), Some(mut bc)) => {
                self.barcode_len.apply(&mut bc)?;
                let barcode = self.layout.extract(&bc);
                let barcode_read = self.keep_barcode_read.then_some(bc);
                Ok(Some(SoloRead {
                    cdna,
                    barcode,
                    barcode_read,
                }))
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
    debug_assert!(matches!(
        params.solo_type,
        SoloType::CbUmiSimple | SoloType::CbUmiComplex | SoloType::CbSamTagOut
    ));
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
    SoloReadReader::open(
        cdna,
        barcode,
        layout,
        params.read_files_command.as_deref(),
        params.solo_keeps_barcode_read(),
        BarcodeReadLength::from_params(params),
    )
}

/// One paired-end solo read for `--soloBarcodeMate 1` (5' 10x): both mates carry
/// cDNA, and the cell barcode (CB+UMI) is a prefix of mate 1.
pub struct SoloPairedRead {
    pub mate1: EncodedRead,
    pub mate2: EncodedRead,
    /// `None` when the barcode read (or mate 1) was too short to extract CB+UMI.
    pub barcode: Option<CellBarcode>,
    /// The whole barcode read, kept only for the `sS`/`sQ` SAM tags. `None` for
    /// a `--soloBarcodeMate 1` run, where mate 1 itself is the barcode read.
    pub barcode_read: Option<EncodedRead>,
}

/// Reads a solo cDNA mate pair in lockstep, taking the barcode either from the
/// start of mate 1 (`--soloBarcodeMate 1`) or from a third barcode-read file
/// (`--readFilesIn cDNA_read1 cDNA_read2 barcode_read`).
pub struct SoloPairedReader {
    mate1: FastqReader,
    mate2: FastqReader,
    /// The separate barcode read, when the run has one.
    barcode: Option<FastqReader>,
    layout: SoloBarcodeLayout,
    keep_barcode_read: bool,
    barcode_len: BarcodeReadLength,
}

impl SoloPairedReader {
    pub fn open(
        mate1_path: &Path,
        mate2_path: &Path,
        barcode_path: Option<&Path>,
        layout: SoloBarcodeLayout,
        decompress_cmd: Option<&str>,
        keep_barcode_read: bool,
        barcode_len: BarcodeReadLength,
    ) -> Result<Self, Error> {
        Ok(Self {
            mate1: FastqReader::open(mate1_path, decompress_cmd)?,
            mate2: FastqReader::open(mate2_path, decompress_cmd)?,
            barcode: barcode_path
                .map(|p| FastqReader::open(p, decompress_cmd))
                .transpose()?,
            layout,
            keep_barcode_read,
            barcode_len,
        })
    }

    /// Fetch the next (mate1, mate2) pair with its barcode.
    pub fn next_read(&mut self) -> Result<Option<SoloPairedRead>, Error> {
        match (self.mate1.next_encoded()?, self.mate2.next_encoded()?) {
            (Some(mate1), Some(mate2)) => {
                let (barcode, barcode_read) = match &mut self.barcode {
                    // Separate barcode read: must stay in lockstep with the mates.
                    Some(reader) => {
                        let mut bc = reader.next_encoded()?.ok_or_else(|| {
                            Error::from(std::io::Error::new(
                                std::io::ErrorKind::UnexpectedEof,
                                "solo: barcode read file has fewer reads than the cDNA mate files",
                            ))
                        })?;
                        self.barcode_len.apply(&mut bc)?;
                        let extracted = self.layout.extract(&bc);
                        (extracted, self.keep_barcode_read.then_some(bc))
                    }
                    // `--soloBarcodeMate 1`: the barcode is a prefix of mate 1.
                    None => (self.layout.extract(&mate1), None),
                };
                Ok(Some(SoloPairedRead {
                    mate1,
                    mate2,
                    barcode,
                    barcode_read,
                }))
            }
            (None, None) => Ok(None),
            _ => Err(Error::from(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "solo: mate 1 and mate 2 read files have different lengths",
            ))),
        }
    }

    /// Read up to `batch_size` pairs for parallel processing.
    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<SoloPairedRead>, Error> {
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

/// Build a [`SoloPairedReader`] for `--soloBarcodeMate 1`, resolving the two cDNA
/// mate files from `--readFilesIn`.
pub fn open_paired_reader(params: &Parameters) -> Result<SoloPairedReader, Error> {
    let (mate1, mate2) = params.solo_cdna_mate_files().ok_or_else(|| {
        Error::from(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "solo: a paired-end solo run requires two --readFilesIn cDNA mate files",
        ))
    })?;
    let layout = SoloBarcodeLayout::from_params(params);
    SoloPairedReader::open(
        mate1,
        mate2,
        params.barcode_read_file().map(std::path::PathBuf::as_path),
        layout,
        params.read_files_command.as_deref(),
        params.solo_keeps_barcode_read(),
        BarcodeReadLength::from_params(params),
    )
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
/// Returns `(clipped_seq, clipped_qual, clip5p, clip3p)` — the CR4-clipped read plus
/// the bases trimmed from the 5' (TSO) and 3' (polyA) ends, so the caller can soft-clip
/// them (STARsolo keeps them in SEQ as soft-clips, e.g. `60M30S`, not dropped).
pub fn clip_adapter_cr4(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>, usize, usize) {
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
        return (seq.to_vec(), qual.to_vec(), 0, 0);
    }
    (
        seq[start..end].to_vec(),
        qual.get(start..end.min(qual.len()))
            .map(<[u8]>::to_vec)
            .unwrap_or_default(),
        start,
        seq.len() - end,
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
    /// Input read index (STAR's `iReadAll`), or [`NO_READ_INDEX`] when the
    /// `CB`/`UB` SAM tags were not requested. Fits in the struct's existing
    /// padding, so tracking it costs no memory.
    pub read_index: u32,
}

/// `read_index` of a record from a run that does not need STAR's readInfo.
pub const NO_READ_INDEX: u32 = u32::MAX;

/// One (cell, UMI, splice-junction) observation for the `SJ` feature. The
/// junction is identified by its absolute intron coordinates; it is mapped to a
/// matrix row (the `SJ.out.tab` order) at output time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SjCountRecord {
    pub cb: u32,
    pub umi: u64,
    pub intron_start: u64,
    pub intron_end: u64,
}

/// One (cell, UMI, gene) observation for the `Velocyto` feature, tagged with the
/// read's spliced/unspliced/ambiguous category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VelocytoRecord {
    pub cb: u32,
    pub umi: u64,
    pub gene: u32,
    pub category: VelocytoCategory,
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
    /// Input read index, as in [`SoloCountRecord`].
    pub read_index: u32,
}

/// A read that mapped to multiple genes (gene-ambiguous). Distributed across its
/// gene set by `--soloMultiMappers` into the `UniqueAndMult-*.mtx` matrices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiGeneRecord {
    pub cb: u32,
    pub umi: u64,
    pub genes: Vec<u32>,
}

/// Thread-safe sink for the records produced during alignment.
#[derive(Default)]
pub struct SoloRecorder {
    pub records: Mutex<Vec<SoloCountRecord>>,
    pub multi_records: Mutex<Vec<SoloMultiRecord>>,
    /// Gene-ambiguous reads for `--soloMultiMappers` (resolved CB only).
    pub multi_gene: Mutex<Vec<MultiGeneRecord>>,
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
    /// Quantified features (`Gene`, `GeneFull`, …), each with its own recorder
    /// and `Solo.out/<feature>/raw/` output. Parallel to `recorders`.
    pub features: Vec<SoloFeature>,
    pub recorders: Vec<SoloRecorder>,
    /// Reads uniquely assigned to a gene per feature (parallel to `features`),
    /// among valid-barcode reads — the STARsolo "Reads Mapped to <feature>:
    /// Unique" metric.
    pub feature_reads: Vec<AtomicU64>,
    /// CellRanger-style positional mapping funnel over uniquely-mapped reads
    /// (independent of barcode), populated only when both `Gene` and `GeneFull`
    /// features run.
    pub region_stats: RegionStats,
    /// `--soloFeatures SJ`: collect per-cell splice-junction counts.
    pub sj_enabled: bool,
    /// (cell, UMI, junction) observations for the SJ feature.
    pub sj_records: Mutex<Vec<SjCountRecord>>,
    /// `--soloFeatures Velocyto`: collect spliced/unspliced/ambiguous counts.
    pub velocyto_enabled: bool,
    /// (cell, UMI, gene, category) observations for the Velocyto feature.
    pub velocyto_records: Mutex<Vec<VelocytoRecord>>,
    /// `--soloMultiMappers` includes a non-`Unique` method → capture gene-
    /// ambiguous reads for distribution into `UniqueAndMult-*.mtx`.
    pub want_multi: bool,
    /// `--soloCellReadStats CB`: the per-cell read summary, or `None` when the
    /// flag is off. Behind a mutex like the other per-read collections; the
    /// lock is only taken when the flag asks for the file.
    pub cell_read_stats: Option<Mutex<crate::solo::cell_reads::CellReadStats>>,
    /// Chromosome indices named by `--genomeChrSetMitochondrial`, for the
    /// `mito` column.
    pub mito_chr: std::collections::HashSet<usize>,
    /// `--soloFeatures Transcript3p`: the transcriptome to assign reads to, and
    /// the accumulated per-read records. Both `None`/absent unless asked for,
    /// since building the transcriptome costs a GTF pass.
    pub transcriptome: Option<crate::quant::transcriptome::TranscriptomeIndex>,
    pub transcript3p: Option<Mutex<crate::solo::transcript3p::Transcript3pAcc>>,
    /// STAR's `readInfo` (`ParametersSolo.cpp:418-435`): the cell and corrected
    /// UMI each read ended up counted under, filled by UMI collapsing and read
    /// back when the sorted BAM is written to fill `CB`/`UB`. `None` unless one
    /// of those tags was requested.
    pub read_info: Option<Mutex<Vec<ReadInfo>>>,
}

/// One `readInfo` entry: what a read was counted as, once collapsing has run.
/// Both fields keep STAR's "undefined" sentinel (all ones), which surfaces as
/// `CB:Z:-` / `UB:Z:-`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadInfo {
    pub cb: u32,
    pub umi: u64,
}

impl Default for ReadInfo {
    fn default() -> Self {
        Self {
            cb: u32::MAX,
            umi: u64::MAX,
        }
    }
}

/// Per-region read tallies for the `Summary.csv` mapping funnel (uniquely-mapped
/// reads, mirroring CellRanger's "confidently mapped to ... regions").
#[derive(Default)]
pub struct RegionStats {
    pub exonic: AtomicU64,
    pub intronic: AtomicU64,
    pub intergenic: AtomicU64,
    pub antisense: AtomicU64,
}

/// What happened to one solo read — one `(record, multi)` per quantified
/// feature, parallel to [`SoloContext::features`].
#[derive(Debug, Default)]
pub struct SoloReadOutcome {
    pub per_feature: Vec<FeatureOutcome>,
    /// SJ-feature records for this read (one per crossed junction); empty unless
    /// `--soloFeatures SJ` and the read is uniquely mapped with a resolved CB.
    pub sj: Vec<SjCountRecord>,
    /// Velocyto record for this read (resolved CB, gene-assigned), if enabled.
    pub velocyto: Option<VelocytoRecord>,
    /// Barcode facts for the per-read SAM tags. `None` when the barcode read was
    /// too short to extract a CB+UMI at all.
    pub barcode: Option<SoloBarcodeTags>,
}

impl SoloReadOutcome {
    /// Stamp the input read index onto every count record this read produced, so
    /// UMI collapsing can fill STAR's readInfo (the `CB`/`UB` SAM tags).
    pub fn set_read_index(&mut self, read_index: u32) {
        for fo in &mut self.per_feature {
            if let Some(r) = &mut fo.record {
                r.read_index = read_index;
            }
            if let Some(m) = &mut fo.multi {
                m.read_index = read_index;
            }
        }
    }
}

/// The `gx`/`gn`/`sF` tag values of one alignment of a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignGeneTag {
    pub gx: String,
    pub gn: String,
    pub sf: [i32; 2],
}

/// What the barcode of one read resolved to, for the STARsolo SAM tags.
///
/// `cb_match` is STAR's `cbMatch` code (the `sM` tag); `cb_index` is the
/// whitelist entry the barcode corrected to, which also seeds the `CB` tag of a
/// `CB_samTagOut` run and the readInfo entry that fills `CB`/`UB` in the sorted
/// BAM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SoloBarcodeTags {
    pub cb_match: i32,
    pub cb_index: Option<u32>,
}

/// The record(s) one read produces for a single feature.
#[derive(Debug, Default)]
pub struct FeatureOutcome {
    /// A resolved count record, if the read was fully assignable.
    pub record: Option<SoloCountRecord>,
    /// A deferred multi-CB record, if the CB was an unresolved 1MM_multi.
    pub multi: Option<SoloMultiRecord>,
    /// A gene-ambiguous record (resolved CB), for `--soloMultiMappers`.
    pub multi_gene: Option<MultiGeneRecord>,
}

impl SoloContext {
    /// Build the solo context from parameters: load the whitelist and build the
    /// gene model from `--sjdbGTFfile`. Call once before alignment.
    pub fn build(params: &Parameters, genome: &crate::genome::Genome) -> Result<Self, Error> {
        let whitelist = if params.solo_type == SoloType::CbUmiComplex {
            // One whitelist per CB segment → combined cartesian-product whitelist.
            let paths: Vec<std::path::PathBuf> = params
                .solo_cb_whitelist
                .iter()
                .map(std::path::PathBuf::from)
                .collect();
            log::info!(
                "STARsolo CB_UMI_Complex: combining {} segment whitelists",
                paths.len()
            );
            let wl = CbWhitelist::load_complex(&paths)?;
            log::info!("STARsolo: {} combined whitelist barcodes", wl.len());
            wl
        } else {
            match params.solo_cb_whitelist_path() {
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
            }
        };

        // Gene model from the GTF (validated to be present for Gene/GeneFull).
        // `CB_samTagOut` quantifies nothing (`Solo.cpp:13` builds no
        // SoloFeature), so it runs without a gene model.
        let tag_out_only = params.solo_type == SoloType::CbSamTagOut;
        let (exons, gene_ann) = if tag_out_only {
            (Vec::new(), GeneAnnotation::default())
        } else {
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
            (exons, gene_ann)
        };

        let strand: SoloStrand = params.solo_strand.parse().map_err(|e: String| {
            Error::from(std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
        })?;

        // Quantified gene features (Gene, GeneFull). Validation guarantees these
        // parse; default to Gene if somehow empty.
        let features: Vec<SoloFeature> = if tag_out_only {
            Vec::new()
        } else {
            let parsed: Vec<SoloFeature> = params
                .solo_features
                .iter()
                .filter_map(|f| f.parse().ok())
                .collect();
            if parsed.is_empty() {
                vec![SoloFeature::Gene]
            } else {
                parsed
            }
        };
        let recorders = features.iter().map(|_| SoloRecorder::new()).collect();
        let feature_reads = features.iter().map(|_| AtomicU64::new(0)).collect();
        let sj_enabled = params.solo_features.iter().any(|f| f == "SJ");
        let velocyto_enabled = params.solo_features.iter().any(|f| f == "Velocyto");
        let transcript3p = params.solo_features.iter().any(|f| f == "Transcript3p");
        let want_multi = params.solo_multi_mappers.iter().any(|m| m != "Unique");

        Ok(Self {
            layout: SoloBarcodeLayout::from_params(params),
            whitelist,
            match_type: params.solo_cb_match_type(),
            strand,
            gene_ann,
            stats: CbMatchStats::new(),
            features,
            recorders,
            feature_reads,
            region_stats: RegionStats::default(),
            sj_enabled,
            sj_records: Mutex::new(Vec::new()),
            velocyto_enabled,
            velocyto_records: Mutex::new(Vec::new()),
            want_multi,
            cell_read_stats: (params.solo_cell_read_stats == "CB")
                .then(|| Mutex::new(crate::solo::cell_reads::CellReadStats::new())),
            mito_chr: params
                .genome_chr_set_mitochondrial
                .iter()
                .filter(|n| n.as_str() != "-")
                .filter_map(|n| genome.chr_name.iter().position(|c| c == n))
                .collect(),
            transcriptome: transcript3p
                .then(|| {
                    crate::quant::transcriptome::TranscriptomeIndex::from_gtf_exons_configured(
                        &exons,
                        genome,
                        &params.sjdb_gtf_tag_exon_parent_transcript,
                        &params.sjdb_gtf_tag_exon_parent_gene,
                        &params.sjdb_gtf_tag_exon_parent_gene_name,
                        &params.sjdb_gtf_tag_exon_parent_gene_type,
                    )
                })
                .transpose()?,
            transcript3p: transcript3p
                .then(|| Mutex::new(crate::solo::transcript3p::Transcript3pAcc::new())),
            // Sized once the read count is known (`reserve_read_info`).
            read_info: params
                .solo_read_info_needed()
                .then(|| Mutex::new(Vec::new())),
        })
    }

    /// Whether this run has to track STAR's readInfo (the `CB`/`UB` SAM tags).
    pub fn read_info_enabled(&self) -> bool {
        self.read_info.is_some()
    }

    /// `--soloType CB_samTagOut`: match the barcode to the whitelist and stop
    /// there. Returns the `sM` facts plus the corrected barcode for the `CB`
    /// tag, the whitelist entry for an Exact/1MM hit, the raw barcode when
    /// there is no whitelist to correct against, and `"-"` otherwise
    /// (`SoloReadBarcode_getCBandUMI.cpp:311-328`).
    pub fn tag_barcode(&self, bc: &CellBarcode) -> (SoloBarcodeTags, String) {
        let cb_match = self
            .whitelist
            .match_cb(&bc.cb_seq, &bc.cb_qual, self.match_type);
        self.stats.record_cb(&cb_match);
        let corrected = match cb_match.resolved_index() {
            Some(idx) => self
                .whitelist
                .barcode_string(idx)
                // No whitelist: STAR passes the barcode through uncorrected.
                .unwrap_or_else(|| bc.cb_string()),
            None => "-".to_string(),
        };
        (
            SoloBarcodeTags {
                cb_match: cb_match.star_code(),
                cb_index: cb_match.resolved_index(),
            },
            corrected,
        )
    }

    /// Size the readInfo array for `n_reads` input reads, once the alignment
    /// pass has counted them. No-op unless `CB`/`UB` were requested.
    pub fn reserve_read_info(&self, n_reads: usize) {
        if let Some(info) = &self.read_info {
            let mut info = info.lock().unwrap();
            info.clear();
            info.resize(n_reads, ReadInfo::default());
        }
    }

    /// Process one solo read: match the cell barcode, validate the UMI, assign
    /// a gene, and (on success) produce a count record. Stats are recorded
    /// here; the returned records are appended to the recorder by the caller.
    /// GX/GN SAM tag strings for a read's `Gene`-feature assignment:
    /// `(gene_id, gene_name)` when uniquely assigned, else `("-", "-")`
    /// (STARsolo convention). Drives `--outSAMattributes GX GN`.
    pub fn gene_tags<'a>(&'a self, transcripts: &[Transcript]) -> (&'a str, &'a str) {
        match assign_gene_se(transcripts, &self.gene_ann, self.strand, self.tag_feature()) {
            GeneAssignment::Gene(g) => (
                self.gene_ann.gene_ids[g as usize].as_str(),
                self.gene_ann.gene_names[g as usize].as_str(),
            ),
            _ => ("-", "-"),
        }
    }

    /// `GX`/`GN` for a paired-end solo read: the pair's single gene, or `("-",
    /// "-")` when it has none or several. Owned strings, since the pair's
    /// effective transcripts are built here.
    pub fn gene_tags_pe(&self, pairs: &[(&Transcript, &Transcript)]) -> (String, String) {
        let mut eff: Vec<Transcript> = Vec::with_capacity(pairs.len() * 2);
        for (m1, m2) in pairs {
            let mut m2c = (*m2).clone();
            m2c.is_reverse = m1.is_reverse;
            eff.push((*m1).clone());
            eff.push(m2c);
        }
        let (gx, gn) = self.gene_tags(&eff);
        (gx.to_string(), gn.to_string())
    }

    /// The feature the gene SAM tags are computed from: STAR's `samAttrFeature`,
    /// which is the first entry of `--soloFeatures`
    /// (`ParametersSolo.cpp:423`).
    fn tag_feature(&self) -> SoloFeature {
        self.features.first().copied().unwrap_or(SoloFeature::Gene)
    }

    /// Per-alignment `gx`/`gn`/`sF` values for one read.
    ///
    /// `gx`/`gn` list every gene of *that* alignment (`;`-joined, `"-"` when it
    /// has none), unlike `GX`/`GN`, which name the read's single gene. `sF`
    /// carries `(overlap type, genes for the read)`, or `(-1, -1)` when the read
    /// overlaps a sense-strand feature but this particular alignment has no gene
    /// (`ReadAlign_alignBAM.cpp:441-473`).
    pub fn align_gene_tags(&self, transcripts: &[Transcript]) -> Vec<AlignGeneTag> {
        let genes = crate::solo::gene::align_genes(
            transcripts,
            &self.gene_ann,
            self.strand,
            self.tag_feature(),
        );
        self.gene_tag_values(&genes, 1)
    }

    /// The same for a paired-end solo read: an alignment is a mate pair, so the
    /// two mates' genes are merged and both records carry the pair's tags.
    pub fn align_gene_tags_pe(&self, pairs: &[(&Transcript, &Transcript)]) -> Vec<AlignGeneTag> {
        // Both mates evaluated against the pair's (mate 1's) strand, as in
        // `process_read_pe`.
        let mut eff: Vec<Transcript> = Vec::with_capacity(pairs.len() * 2);
        for (m1, m2) in pairs {
            let mut m2c = (*m2).clone();
            m2c.is_reverse = m1.is_reverse;
            eff.push((*m1).clone());
            eff.push(m2c);
        }
        let genes =
            crate::solo::gene::align_genes(&eff, &self.gene_ann, self.strand, self.tag_feature());
        self.gene_tag_values(&genes, 2)
    }

    /// Render `AlignGenes` into per-alignment tag strings, merging every
    /// `mates_per_align` consecutive entries into one alignment.
    fn gene_tag_values(
        &self,
        genes: &crate::solo::gene::AlignGenes,
        mates_per_align: usize,
    ) -> Vec<AlignGeneTag> {
        let sense = matches!(
            genes.ov_type,
            x if x == crate::solo::gene::OverlapType::ExonicSense as i32
                || x == crate::solo::gene::OverlapType::ExonicSense50p as i32
                || x == crate::solo::gene::OverlapType::IntronicSense as i32
        );
        let n_genes = genes.gene_set.len() as i32;
        genes
            .per_align
            .chunks(mates_per_align.max(1))
            .map(|mates| {
                let mut gs: Vec<u32> = mates.concat();
                gs.sort_unstable();
                gs.dedup();
                let join = |names: &[String]| {
                    if gs.is_empty() {
                        "-".to_string()
                    } else {
                        gs.iter()
                            .map(|&g| names[g as usize].as_str())
                            .collect::<Vec<_>>()
                            .join(";")
                    }
                };
                AlignGeneTag {
                    gx: join(&self.gene_ann.gene_ids),
                    gn: join(&self.gene_ann.gene_names),
                    sf: if sense && gs.is_empty() {
                        [-1, -1]
                    } else {
                        [genes.ov_type, n_genes]
                    },
                }
            })
            .collect()
    }

    pub fn process_read(
        &self,
        cdna_transcripts: &[Transcript],
        n_loci: usize,
        barcode: Option<&CellBarcode>,
        junctions: &[(u64, u64)],
    ) -> SoloReadOutcome {
        let mut out = SoloReadOutcome::default();

        // One-pass classification: the two overlap queries are shared between the
        // per-feature gene assignment and the CellRanger-style mapping funnel, so
        // this is no more work than the old per-feature `assign_gene_se` calls.
        let want_exon = self.features.contains(&SoloFeature::Gene);
        // Velocyto assigns its gene by gene-body overlap, so it needs `want_body`.
        let want_body = self.features.contains(&SoloFeature::GeneFull) || self.velocyto_enabled;
        let class = classify_read(
            cdna_transcripts,
            &self.gene_ann,
            self.strand,
            want_exon,
            want_body,
            self.want_multi,
        );

        // Mapping funnel: count uniquely-mapped reads by region (CellRanger's
        // "confidently mapped" = MAPQ 255 ≈ a single alignment), independent of
        // barcode validity. `n_loci` is the number of genomic loci the read (or
        // pair) maps to — for SE this equals `cdna_transcripts.len()`.
        if n_loci == 1 {
            match class.region {
                Some(Region::Exonic) => {
                    self.region_stats.exonic.fetch_add(1, Ordering::Relaxed);
                }
                Some(Region::Intronic) => {
                    self.region_stats.intronic.fetch_add(1, Ordering::Relaxed);
                }
                Some(Region::Intergenic) => {
                    self.region_stats.intergenic.fetch_add(1, Ordering::Relaxed);
                }
                None => {}
            }
            if class.antisense {
                self.region_stats.antisense.fetch_add(1, Ordering::Relaxed);
            }
        }

        // No barcode read (too short) → nothing to count (region already tallied).
        let Some(bc) = barcode else {
            return out;
        };

        // Cell-barcode match.
        let cb_match = self
            .whitelist
            .match_cb(&bc.cb_seq, &bc.cb_qual, self.match_type);
        self.stats.record_cb(&cb_match);

        // Per-read SAM-tag facts, recorded before any early return so that a
        // rejected barcode still tags its alignments. STAR reports a UMI
        // rejection in place of the CB code (`getCBandUMI.cpp:304`).
        let umi_status = check_umi(&bc.umi_seq);
        out.barcode = Some(SoloBarcodeTags {
            cb_match: umi_status
                .star_code()
                .unwrap_or_else(|| cb_match.star_code()),
            cb_index: cb_match.resolved_index(),
        });

        let cb_resolved: Option<u32> = match &cb_match {
            CbMatch::Exact(idx) | CbMatch::Corrected(idx) => Some(*idx),
            CbMatch::Multi(_) => None, // deferred to collation
            CbMatch::NoMatch | CbMatch::NinCb | CbMatch::MultMatchRejected => return out,
        };

        // UMI validity.
        let umi = match umi_status {
            UmiCheck::Ok(packed) => {
                self.stats.record_umi(&UmiCheck::Ok(packed));
                packed
            }
            rejected => {
                self.stats.record_umi(&rejected);
                return out;
            }
        };

        // SJ feature: record (cell, UMI, junction) for each crossed junction.
        // Only for resolved CBs (1MM_multi deferral is not applied to SJ).
        if self.sj_enabled
            && !junctions.is_empty()
            && let Some(cb) = cb_resolved
        {
            out.sj = junctions
                .iter()
                .map(|&(intron_start, intron_end)| SjCountRecord {
                    cb,
                    umi,
                    intron_start,
                    intron_end,
                })
                .collect();
        }

        // Velocyto feature: gene from gene-body overlap, then classify the read
        // spliced/unspliced/ambiguous. Resolved CB only.
        if self.velocyto_enabled
            && let Some(cb) = cb_resolved
            && let GeneAssignment::Gene(gene) = class.gene_full
        {
            out.velocyto = Some(VelocytoRecord {
                cb,
                umi,
                gene,
                category: velocyto_category(cdna_transcripts, &self.gene_ann, gene),
            });
        }

        // The CB match + UMI are shared across features; reuse the cached
        // per-feature gene assignment from `classify_read`. One outcome/feature.
        out.per_feature = self
            .features
            .iter()
            .enumerate()
            .map(|(fi, &feature)| {
                let mut fo = FeatureOutcome::default();
                let assignment = match feature {
                    SoloFeature::Gene => class.gene,
                    SoloFeature::GeneFull => class.gene_full,
                };
                let gene = match assignment {
                    GeneAssignment::Gene(g) => g,
                    GeneAssignment::Ambiguous => {
                        // Gene-ambiguous read: record its gene set for
                        // --soloMultiMappers distribution (resolved CB only).
                        if let Some(cb) = cb_resolved {
                            let genes = match feature {
                                SoloFeature::Gene => &class.gene_multi,
                                SoloFeature::GeneFull => &class.gene_full_multi,
                            };
                            if !genes.is_empty() {
                                fo.multi_gene = Some(MultiGeneRecord {
                                    cb,
                                    umi,
                                    genes: genes.clone(),
                                });
                            }
                        }
                        return fo;
                    }
                    GeneAssignment::NoFeature | GeneAssignment::Unmapped => return fo,
                };
                // Reads uniquely mapped to a gene under this feature, among
                // valid-barcode reads (STARsolo "Reads Mapped to <feature>").
                self.feature_reads[fi].fetch_add(1, Ordering::Relaxed);
                match (cb_resolved, &cb_match) {
                    (Some(cb), _) => {
                        fo.record = Some(SoloCountRecord {
                            cb,
                            umi,
                            gene,
                            read_index: NO_READ_INDEX,
                        });
                    }
                    (None, CbMatch::Multi(cands)) => {
                        fo.multi = Some(SoloMultiRecord {
                            candidates: cands.clone(),
                            umi,
                            gene,
                            read_index: NO_READ_INDEX,
                        });
                    }
                    (None, _) => unreachable!("non-multi unresolved CB returned early"),
                }
                fo
            })
            .collect();

        self.record_cell_read(
            cb_resolved,
            &cb_match,
            n_loci,
            &class,
            cdna_transcripts,
            &out,
        );
        // Transcript3p: every transcript this read is concordant with, and how
        // far its 3' end sits from each transcript's. Uniquely-mapped reads
        // only — a read at several genomic loci says nothing about isoforms.
        if let (Some(acc), Some(tx), Some(cb)) =
            (&self.transcript3p, &self.transcriptome, cb_resolved)
            && n_loci == 1
            && let Some(align) = cdna_transcripts.first()
        {
            let hits = crate::solo::transcript3p::concordant_transcripts(
                align,
                tx,
                align.read_length() as u32,
            );
            if !hits.is_empty() {
                acc.lock().unwrap().add(cb, umi, hits);
            }
        }

        out
    }

    /// Fold one read into `CellReads.stats`, when `--soloCellReadStats CB`
    /// asked for it.
    ///
    /// Called at the end of read processing so the counted flags reflect what
    /// the read actually produced, rather than what it looked eligible for.
    #[allow(clippy::too_many_arguments)]
    fn record_cell_read(
        &self,
        cb_resolved: Option<u32>,
        cb_match: &CbMatch,
        n_loci: usize,
        class: &crate::solo::gene::ReadClass,
        transcripts: &[Transcript],
        out: &SoloReadOutcome,
    ) {
        let Some(stats) = &self.cell_read_stats else {
            return;
        };
        let counted_u = out.per_feature.iter().any(|f| f.record.is_some());
        let counted_m = out.per_feature.iter().any(|f| f.multi_gene.is_some());
        let feature_u = counted_u || out.per_feature.iter().any(|f| f.multi.is_some());
        let flag = crate::solo::cell_reads::CellReadFlag {
            cb_perfect: matches!(cb_match, CbMatch::Exact(_)),
            cb_mm_unique: matches!(cb_match, CbMatch::Corrected(_)),
            cb_mm_multiple: matches!(cb_match, CbMatch::Multi(_)),
            genome_u: n_loci == 1,
            genome_m: n_loci > 1,
            feature_u,
            feature_m: counted_m,
            // The region columns split by strand: STAR reports an antisense
            // read under `exonicAS`/`intronicAS`, not under `exonic`/`intronic`.
            exonic: !class.antisense && class.region == Some(Region::Exonic),
            intronic: !class.antisense && class.region == Some(Region::Intronic),
            exonic_as: class.antisense && class.region == Some(Region::Exonic),
            intronic_as: class.antisense && class.region == Some(Region::Intronic),
            mito: !self.mito_chr.is_empty()
                && transcripts
                    .iter()
                    .any(|t| self.mito_chr.contains(&t.chr_idx)),
            counted_u,
            counted_m,
        };
        let mut stats = stats.lock().unwrap();
        match cb_resolved {
            Some(cb) => stats.add_cell(cb, &flag),
            None => stats.add_no_cb(&flag),
        }
    }

    /// Process one 5' paired-end solo read (`--soloBarcodeMate 1`): the barcode is
    /// from mate 1, and both mates align as a pair. Genes are assigned from the
    /// union of both mates evaluated against the pair's (mate 1's) transcription
    /// strand — matching STAR's PE quantification (cf. [`crate::quant`]'s
    /// `count_pe_read`). Each aligned pair counts as one locus.
    ///
    /// `pairs` are the `(mate1, mate2)` transcripts of every BothMapped pair (>1
    /// pair ⇒ multimapper). `junctions` are the pair's crossed junctions.
    pub fn process_read_pe(
        &self,
        pairs: &[(&Transcript, &Transcript)],
        barcode: Option<&CellBarcode>,
        junctions: &[(u64, u64)],
    ) -> SoloReadOutcome {
        // Effective transcripts: give both mates of a pair the pair's strand
        // (mate 1's), so `classify_read`'s per-transcript strand filter treats them
        // as one stranded observation. Genomic overlap is unaffected by is_reverse.
        let mut eff: Vec<Transcript> = Vec::with_capacity(pairs.len() * 2);
        for (m1, m2) in pairs {
            let mut m2c = (*m2).clone();
            m2c.is_reverse = m1.is_reverse;
            eff.push((*m1).clone());
            eff.push(m2c);
        }
        self.process_read(&eff, pairs.len(), barcode, junctions)
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
        SoloBarcodeLayout::Simple {
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
        assert_eq!(
            layout,
            SoloBarcodeLayout::Simple {
                cb_start: 0,
                cb_len: 16,
                umi_start: 16,
                umi_len: 10,
            }
        );
        assert_eq!(layout.min_read_len(), 26);
    }

    #[test]
    fn complex_layout_assembles_segments() {
        // Two CB segments [0..2] + [4..6] (skipping a 2bp linker), UMI [6..8].
        let layout = SoloBarcodeLayout::Complex {
            cb_segments: vec![(0, 2), (4, 2)],
            umi: (6, 2),
        };
        let read = encoded_read("r", "AACCGGTT", "IIIIIIII");
        let bc = layout.extract(&read).unwrap();
        // CB = bases [0,1] ++ [4,5] = "AA" ++ "GG"; UMI = [6,7] = "TT".
        assert_eq!(
            bc.cb_seq,
            "AAGG".bytes().map(encode_base).collect::<Vec<_>>()
        );
        assert_eq!(
            bc.umi_seq,
            "TT".bytes().map(encode_base).collect::<Vec<_>>()
        );
    }

    #[test]
    fn parse_position_read_start() {
        assert_eq!(parse_position("0_0_0_7").unwrap(), (0, 8));
        assert_eq!(parse_position("0_8_0_15").unwrap(), (8, 8));
        assert!(parse_position("2_0_2_7").is_err()); // adapter anchor unsupported
        assert!(parse_position("0_5_0_2").is_err()); // end < start
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
    fn extract_barcode_from_mate1_with_cdna_tail() {
        // `--soloBarcodeMate 1`: mate 1 = [CB+UMI][cDNA]. The barcode is read from
        // the fixed prefix; the trailing cDNA is ignored (same CB/UMI as a
        // barcode-only read).
        let layout = v2_layout();
        let prefix = "AAAAAAAACCCCCCCCGGGGGTTTTT"; // 16 CB + 10 UMI
        let tail = "TTTTGGGGCCCCAAAATTTTGGGG"; // 24 bp of cDNA
        let read = encoded_read(
            "m1",
            &format!("{prefix}{tail}"),
            &"I".repeat(prefix.len() + tail.len()),
        );
        let bc = layout
            .extract(&read)
            .expect("should extract despite cDNA tail");
        assert_eq!(bc.cb_string(), "AAAAAAAACCCCCCCC");
        assert_eq!(bc.umi_string(), "GGGGGTTTTT");
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

    /// STAR checks the barcode read's length by default and refuses anything
    /// else; with the check off, a short read is padded with `N` (which then
    /// scores as an N-containing barcode rather than being silently dropped).
    #[test]
    fn barcode_read_length_is_checked_then_padded() {
        use crate::io::fastq::encode_base;

        let read = |seq: &str| EncodedRead {
            name: "r1".to_string(),
            sequence: seq.bytes().map(encode_base).collect(),
            quality: vec![b'I'; seq.len()],
        };
        // Default: exactly CB+UMI = 26 bases.
        let checked = BarcodeReadLength {
            expected: Some(26),
            cbumi_len: 26,
        };
        assert!(checked.apply(&mut read("A".repeat(26).as_str())).is_ok());
        let err = checked
            .apply(&mut read("A".repeat(20).as_str()))
            .unwrap_err()
            .to_string();
        assert!(err.contains("20 bases, not the expected 26"), "{err}");
        // Longer is refused too.
        assert!(checked.apply(&mut read("A".repeat(30).as_str())).is_err());

        // --soloBarcodeReadLength 0: pad up to CB+UMI with N (encoded 4) / 'H'.
        let unchecked = BarcodeReadLength {
            expected: None,
            cbumi_len: 26,
        };
        let mut short = read("ACGT");
        unchecked.apply(&mut short).unwrap();
        assert_eq!(short.sequence.len(), 26);
        assert!(short.sequence[4..].iter().all(|&b| b == 4));
        assert!(short.quality[4..].iter().all(|&q| q == b'H'));
        // A read at or past the length is left alone.
        let mut long = read("A".repeat(30).as_str());
        unchecked.apply(&mut long).unwrap();
        assert_eq!(long.sequence.len(), 30);
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

        let mut reader = SoloReadReader::open(
            cdna.path(),
            bc.path(),
            v2_layout(),
            None,
            false,
            BarcodeReadLength::default(),
        )
        .unwrap();
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

        let mut reader = SoloReadReader::open(
            cdna.path(),
            bc.path(),
            v2_layout(),
            None,
            false,
            BarcodeReadLength::default(),
        )
        .unwrap();
        assert!(reader.read_batch(10).is_err());
    }
}
