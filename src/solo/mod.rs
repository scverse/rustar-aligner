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

pub use count::{UmiDedup, UmiFiltering, write_barcode_stats, write_gene_matrix};
pub use gene::{
    GeneAssignment, GeneModel, GeneSets, ReadClass, Region, SoloFeature, SoloStrand,
    assign_gene_se, classify_align, classify_read, gene_full_sets, velocyto_transcript_types,
};
pub use whitelist::{
    BarcodeStats, CbCandidate, CbMatch, CbMatchType, CbWhitelist, UmiCheck, check_umi, pack_barcode,
};

use crate::align::transcript::Transcript;
use crate::error::Error;
use crate::io::fastq::{EncodedRead, FastqReader, decode_base};
use crate::params::{Parameters, SoloType};
use crate::quant::GeneAnnotation;
use std::path::PathBuf;
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
    /// Multi-segment CB whose positions are measured from an adapter rather
    /// than from the read ends.
    ///
    /// The adapter moves from read to read, so unlike every other layout the
    /// offsets cannot be resolved once at startup; each read is searched. A
    /// read where the adapter is not found within the mismatch budget yields no
    /// barcode at all, which is what STAR does — a barcode read at an unknown
    /// offset is not salvageable.
    Anchored {
        cb_segments: Vec<(Anchored, Anchored)>,
        umi: (Anchored, Anchored),
        adapter: Vec<u8>,
        max_mismatches: usize,
    },
}

/// Parse a `--soloCBposition`/`--soloUMIposition` spec
/// (`startAnchor_startDist_endAnchor_endDist`) into a 0-based `(start, len)`.
/// Only read-start anchoring (`anchor = 0`) is supported.
/// One end of a `--soloCBposition` / `--soloUMIposition` field: which anchor
/// it is measured from, and how far.
///
/// STAR's anchor codes: 0 = read start, 1 = read end, 2 = adapter start,
/// 3 = adapter end. Codes 2 and 3 make the position depend on where the
/// adapter turns out to be in each read, which is why resolution happens per
/// read rather than once at startup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchored {
    pub anchor: u8,
    pub dist: i64,
}

impl Anchored {
    /// Position of this anchor plus its distance, as STAR's `extractBarcode`
    /// computes it (`SoloBarcode_extractBarcode.cpp:7-23`), for a read of
    /// length `read_len` whose adapter starts at `adapter_start` and is
    /// `adapter_len` long. Bounds are checked by the caller, as in STAR.
    fn position(self, read_len: usize, adapter_start: usize, adapter_len: usize) -> i64 {
        let base = match self.anchor {
            0 => 0i64,
            1 => read_len as i64 - 1,
            2 => adapter_start as i64,
            _ => adapter_start as i64 + adapter_len as i64 - 1,
        };
        base + self.dist
    }
}

/// STAR's `extractBarcode` bounds check and `substr`
/// (`SoloBarcode_extractBarcode.cpp:25-32`): refused when the start is
/// negative, the end is past the read length or the start is after the end.
/// An end equal to the read length is allowed and `substr` clamps it.
fn barcode_span(start: i64, end: i64, read_len: usize) -> Option<(usize, usize)> {
    if start < 0 || end > read_len as i64 || start > end {
        return None;
    }
    let s0 = start as usize;
    Some((s0, (end as usize + 1).min(read_len) - s0))
}

/// STAR's `localAlignHammingDist` (`SequenceFuns.cpp:341-360`): the
/// fewest-mismatch offset of `adapter` in `seq` and its mismatch count. An `N`
/// in the adapter never counts as a mismatch; ties go to the leftmost offset
/// (the incumbent is replaced only on a strictly smaller distance). The
/// distance starts at the adapter length with the offset at 0, so when no
/// offset does better the result is `(adapter length, 0)`; a read shorter than
/// the adapter gives `(read length + 1, 0)`.
fn local_align_hamming_dist(seq: &[u8], adapter: &[u8]) -> (usize, usize) {
    const N_CODE: u8 = 4;
    if seq.len() < adapter.len() {
        return (seq.len() + 1, 0);
    }
    let (mut dist_best, mut pos) = (adapter.len(), 0usize);
    for start in 0..=(seq.len() - adapter.len()) {
        let dist = seq[start..start + adapter.len()]
            .iter()
            .zip(adapter)
            .filter(|&(a, b)| *b != N_CODE && a != b)
            .count();
        if dist < dist_best {
            dist_best = dist;
            pos = start;
        }
    }
    (dist_best, pos)
}

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

/// Parse a full `startAnchor_startDist_endAnchor_endDist` field, keeping the
/// anchors rather than collapsing them to read-start offsets.
fn parse_position_anchored(spec: &str) -> Result<(Anchored, Anchored), Error> {
    let f: Vec<&str> = spec.split('_').collect();
    if f.len() != 4 {
        return Err(invalid_pos(
            spec,
            "expected startAnchor_startDist_endAnchor_endDist",
        ));
    }
    let num = |x: &str| x.parse::<i64>().ok();
    match (num(f[0]), num(f[1]), num(f[2]), num(f[3])) {
        (Some(sa), Some(sd), Some(ea), Some(ed))
            if (0..=3).contains(&sa) && (0..=3).contains(&ea) =>
        {
            Ok((
                Anchored {
                    anchor: sa as u8,
                    dist: sd,
                },
                Anchored {
                    anchor: ea as u8,
                    dist: ed,
                },
            ))
        }
        _ => Err(invalid_pos(spec, "anchor must be 0, 1, 2 or 3")),
    }
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

            // An adapter turns the position specs into per-read offsets, so it
            // selects a different layout rather than being a tweak to this one.
            if params.solo_adapter_sequence != "-" {
                let adapter: Vec<u8> = params
                    .solo_adapter_sequence
                    .bytes()
                    .map(crate::io::fastq::encode_base)
                    .collect();
                let anchored_cb: Vec<_> = params
                    .solo_cb_position
                    .iter()
                    .filter_map(|s| parse_position_anchored(s).ok())
                    .collect();
                if let Ok(anchored_umi) = parse_position_anchored(&params.solo_umi_position) {
                    return Self::Anchored {
                        cb_segments: anchored_cb,
                        umi: anchored_umi,
                        adapter,
                        max_mismatches: params.solo_adapter_mismatches_nmax,
                    };
                }
            }
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
            // Not knowable without the read: the adapter's position decides.
            // `extract` does its own bounds checking instead.
            Self::Anchored { .. } => 0,
            Self::Complex { cb_segments, umi } => cb_segments
                .iter()
                .map(|&(s, l)| s + l)
                .chain(std::iter::once(umi.0 + umi.1))
                .max()
                .unwrap_or(0),
        }
    }

    /// Extract the CB (concatenating segments for `Complex`) and UMI from one
    /// barcode read, or why it cannot be (STAR's `getCBandUMI`: no adapter,
    /// then no UMI, then no CB).
    pub fn extract(&self, barcode_read: &EncodedRead) -> Result<CellBarcode, BarcodeMiss> {
        let seq = &barcode_read.sequence;
        let qual = &barcode_read.quality;
        match self {
            Self::Simple {
                cb_start,
                cb_len,
                umi_start,
                umi_len,
            } => {
                if seq.len() < self.min_read_len() {
                    return Err(BarcodeMiss::NoCb);
                }
                Ok(CellBarcode {
                    cb_seq: seq[*cb_start..cb_start + cb_len].to_vec(),
                    cb_qual: slice_or_empty(qual, *cb_start, *cb_len),
                    umi_seq: seq[*umi_start..umi_start + umi_len].to_vec(),
                    umi_qual: slice_or_empty(qual, *umi_start, *umi_len),
                })
            }
            Self::Anchored {
                cb_segments,
                umi,
                adapter,
                max_mismatches,
            } => {
                // `SoloReadBarcode_getCBandUMI.cpp:336-374`: no barcode when
                // the adapter is not found within the budget (`noNoAdapter`),
                // whatever the anchors; then the UMI, then each CB segment.
                let (dist, adapter_start) = local_align_hamming_dist(seq, adapter);
                if dist > *max_mismatches {
                    return Err(BarcodeMiss::NoAdapter);
                }
                let span = |(a, b): &(Anchored, Anchored)| {
                    barcode_span(
                        a.position(seq.len(), adapter_start, adapter.len()),
                        b.position(seq.len(), adapter_start, adapter.len()),
                        seq.len(),
                    )
                };
                let (us, ul) = span(umi).ok_or(BarcodeMiss::NoUmi)?;
                let mut cb_seq = Vec::new();
                let mut cb_qual = Vec::new();
                for segment in cb_segments {
                    let (s0, l) = span(segment).ok_or(BarcodeMiss::NoCb)?;
                    cb_seq.extend_from_slice(&seq[s0..s0 + l]);
                    cb_qual.extend_from_slice(&slice_or_empty(qual, s0, l));
                }
                Ok(CellBarcode {
                    cb_seq,
                    cb_qual,
                    umi_seq: seq[us..us + ul].to_vec(),
                    umi_qual: slice_or_empty(qual, us, ul),
                })
            }
            Self::Complex { cb_segments, umi } => {
                if seq.len() < umi.0 + umi.1 {
                    return Err(BarcodeMiss::NoUmi);
                }
                if seq.len() < self.min_read_len() {
                    return Err(BarcodeMiss::NoCb);
                }
                let mut cb_seq = Vec::new();
                let mut cb_qual = Vec::new();
                for &(s, l) in cb_segments {
                    cb_seq.extend_from_slice(&seq[s..s + l]);
                    cb_qual.extend_from_slice(&slice_or_empty(qual, s, l));
                }
                Ok(CellBarcode {
                    cb_seq,
                    cb_qual,
                    umi_seq: seq[umi.0..umi.0 + umi.1].to_vec(),
                    umi_qual: slice_or_empty(qual, umi.0, umi.1),
                })
            }
        }
    }
}

/// Why no barcode could be read, as STAR counts it in `Barcodes.stats`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarcodeMiss {
    /// No adapter within `--soloAdapterMismatchesNmax` (`noNoAdapter`).
    NoAdapter,
    /// The UMI cannot be extracted (`noNoUMI`).
    NoUmi,
    /// A cell barcode cannot be extracted (`noNoCB`, `cbMatch` -11).
    NoCb,
}

/// The `(Q30 bases, bases)` of a barcode read that STAR adds to its CB+UMI
/// quality histogram (`SoloReadBarcode_getCBandUMI.cpp:243-266`): the CB and
/// UMI qualities for `CB_UMI_Simple`, else the whole barcode read (its first
/// `--soloBarcodeReadLength` bases when that is set), extracted or not.
pub fn barcode_q30(
    layout: &SoloBarcodeLayout,
    read: &EncodedRead,
    bc: Result<&CellBarcode, BarcodeMiss>,
    read_len: Option<usize>,
) -> (u64, u64) {
    let q30 = |q: &[u8]| {
        (
            q.iter().filter(|&&b| b >= 33 + 30).count() as u64,
            q.len() as u64,
        )
    };
    match (layout, bc) {
        (SoloBarcodeLayout::Simple { .. }, Ok(bc)) => {
            let (a, b) = q30(&bc.cb_qual);
            let (c, d) = q30(&bc.umi_qual);
            (a + c, b + d)
        }
        (SoloBarcodeLayout::Simple { .. }, Err(_)) => (0, 0),
        _ => {
            let n = read_len
                .unwrap_or(read.quality.len())
                .min(read.quality.len());
            q30(&read.quality[..n])
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
                    // `SoloReadBarcode_getCBandUMI.cpp:228-236`.
                    format!(
                        "EXITING because of FATAL ERROR in input read file: the total length of \
                         barcode sequence is {} not equal to expected {}\n\
                         Read ID={} ;  Sequence={}\n\
                         SOLUTION: check the formatting of input read files.\n\
                         If UMI+CB length is not equal to the barcode read length, specify \
                         barcode read length with --soloBarcodeReadLength\n\
                         To avoid checking of barcode read length, specify \
                         --soloBarcodeReadLength 0",
                        read.sequence.len(),
                        expected,
                        read.name,
                        read.sequence
                            .iter()
                            .map(|&b| crate::io::fastq::decode_base(b) as char)
                            .collect::<String>()
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
    /// The extracted barcode, or why there is none.
    pub barcode: Result<CellBarcode, BarcodeMiss>,
    /// The barcode read's `(Q30 bases, bases)` for `Summary.csv` ([`barcode_q30`]).
    pub barcode_q30: (u64, u64),
    /// The whole barcode read, kept only for the `sS`/`sQ` SAM tags.
    pub barcode_read: Option<EncodedRead>,
}

impl SoloReadReader {
    /// Open the cDNA and barcode FASTQ files for a solo run.
    pub fn open(
        cdna_path: &[PathBuf],
        barcode_path: &[PathBuf],
        layout: SoloBarcodeLayout,
        decompress_cmd: Option<&str>,
        keep_barcode_read: bool,
        barcode_len: BarcodeReadLength,
    ) -> Result<Self, Error> {
        Ok(Self {
            cdna: FastqReader::open_files(cdna_path, decompress_cmd)?,
            barcode: FastqReader::open_files(barcode_path, decompress_cmd)?,
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
                let barcode_q30 = barcode_q30(
                    &self.layout,
                    &bc,
                    barcode.as_ref().map_err(|m| *m),
                    self.barcode_len.expected,
                );
                let barcode_read = self.keep_barcode_read.then_some(bc);
                Ok(Some(SoloRead {
                    cdna,
                    barcode,
                    barcode_q30,
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
    /// The extracted barcode, or why there is none.
    pub barcode: Result<CellBarcode, BarcodeMiss>,
    /// The barcode read's `(Q30 bases, bases)` for `Summary.csv` ([`barcode_q30`]).
    pub barcode_q30: (u64, u64),
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
        mate1_path: &[PathBuf],
        mate2_path: &[PathBuf],
        barcode_path: Option<&[PathBuf]>,
        layout: SoloBarcodeLayout,
        decompress_cmd: Option<&str>,
        keep_barcode_read: bool,
        barcode_len: BarcodeReadLength,
    ) -> Result<Self, Error> {
        Ok(Self {
            mate1: FastqReader::open_files(mate1_path, decompress_cmd)?,
            mate2: FastqReader::open_files(mate2_path, decompress_cmd)?,
            barcode: barcode_path
                .map(|p| FastqReader::open_files(p, decompress_cmd))
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
                // A separate barcode read must stay in lockstep with the mates;
                // with `--soloBarcodeMate 1` the barcode is a prefix of mate 1.
                let (barcode, barcode_q30, barcode_read) = if let Some(reader) = &mut self.barcode {
                    let mut bc = reader.next_encoded()?.ok_or_else(|| {
                        Error::from(std::io::Error::new(
                            std::io::ErrorKind::UnexpectedEof,
                            "solo: barcode read file has fewer reads than the cDNA mate files",
                        ))
                    })?;
                    self.barcode_len.apply(&mut bc)?;
                    let extracted = self.layout.extract(&bc);
                    let q = barcode_q30(
                        &self.layout,
                        &bc,
                        extracted.as_ref().map_err(|m| *m),
                        self.barcode_len.expected,
                    );
                    (extracted, q, self.keep_barcode_read.then_some(bc))
                } else {
                    let extracted = self.layout.extract(&mate1);
                    let q = barcode_q30(
                        &self.layout,
                        &mate1,
                        extracted.as_ref().map_err(|m| *m),
                        None,
                    );
                    (extracted, q, None)
                };
                Ok(Some(SoloPairedRead {
                    mate1,
                    mate2,
                    barcode,
                    barcode_q30,
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
        params.barcode_read_file(),
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

/// STAR aligns the TSO against the first 91 bases of the read (`ClipCR4::readLen`).
const CR4_READ_LEN: usize = 91;

/// CR4 overlap-alignment scoring, matching `ClipCR4()`: match +1, mismatch -2,
/// any-vs-N -2, N-vs-N 0; gap open = gap extend = 2. Alphabet A,C,G,T,N = 0..4.
fn cr4_scoring() -> &'static hyalite::Scoring {
    static SCORING: std::sync::OnceLock<hyalite::Scoring> = std::sync::OnceLock::new();
    SCORING.get_or_init(|| {
        #[rustfmt::skip]
        let matrix = vec![
             1, -2, -2, -2, -2,
            -2,  1, -2, -2, -2,
            -2, -2,  1, -2, -2,
            -2, -2, -2,  1, -2,
            -2, -2, -2, -2,  0,
        ];
        hyalite::Scoring::new(5, matrix, 2, 2).expect("valid CR4 scoring matrix")
    })
}

/// The TSO query, encoded once (0=A..3=T).
fn tso_query() -> &'static [u8] {
    static Q: std::sync::OnceLock<Vec<u8>> = std::sync::OnceLock::new();
    Q.get_or_init(|| {
        TSO_SEQ
            .iter()
            .map(|&b| crate::io::fastq::encode_base(b))
            .collect()
    })
}

/// STAR-faithful 5' TSO clip length via a hyalite overlap alignment, replacing
/// STAR's Opal call (`ClipCR4::opalAlign` + `ClipMate_clipChunk.cpp:41-49`).
///
/// `read` is encoded (0=A..4=N). The target is the first `min(len, 91)` bases,
/// N-padded to 91 (exactly `ClipCR4::opalFillOneSeq`); query is the TSO; overlap
/// mode with end tracking. The clip is `L = target_end + 1`, rejected to 0 by
/// STAR's score gate: `S < 20 || (S == 20 && L > 26) || (S == 21 && L > 30)`.
fn tso_clip_len_cr4(read: &[u8]) -> usize {
    if read.is_empty() {
        return 0;
    }
    let target = cr4_target(read);

    let Ok(hit) = hyalite::align_pair(
        tso_query(),
        &target,
        cr4_scoring(),
        hyalite::Mode::Ov,
        hyalite::SearchType::ScoreEnd,
    ) else {
        return 0;
    };

    cr4_gate(&hit, read.len())
}

/// Apply STAR's clip gate to one hyalite hit, yielding the 5' clip length.
///
/// `S < 20 || (S == 20 && L > 26) || (S == 21 && L > 30)` rejects to 0
/// (`ClipMate_clipChunk.cpp:44-47`), then STAR takes `min(clip, Lread)`
/// (`ClipMate_clip.cpp:53`).
fn cr4_gate(hit: &hyalite::BestHit, read_len: usize) -> usize {
    let s = hit.score;
    let l = hit.target_end.map_or(0, |e| e + 1) as i32; // ScoreEnd ⇒ Some
    let reject = s < 20 || (s == 20 && l > 26) || (s == 21 && l > 30);
    if reject {
        0
    } else {
        (l as usize).min(read_len)
    }
}

/// Build the N-padded 91-base target for one read, exactly `ClipCR4::opalFillOneSeq`.
fn cr4_target(read: &[u8]) -> Vec<u8> {
    let take = read.len().min(CR4_READ_LEN);
    let mut target = Vec::with_capacity(CR4_READ_LEN);
    target.extend_from_slice(&read[..take]);
    target.resize(CR4_READ_LEN, 4); // pad with N
    target
}

/// How many reads share one hyalite `Database` in the batched TSO scan. STAR uses
/// 64 (`ClipCR4::dbN`); throughput here is flat from 128 to 10k, so this is sized
/// to give rayon enough chunks to spread across cores rather than to tune SIMD.
const CR4_SCAN_CHUNK: usize = 512;

/// Batched 5' TSO clip lengths for a whole read batch — the SIMD counterpart to
/// [`tso_clip_len_cr4`], via hyalite's per-target `Database::scan_all`.
///
/// STAR does the same thing: `ClipMate::clipChunk` fills an Opal "database" with
/// a chunk of reads and aligns the adapter against all of them in one call. The
/// results are bit-identical to the per-read path on every backend (see the
/// `cr4_tso_scan_all_batch_matches_scalar_and_is_backend_stable` test), and ~7x
/// faster than per-read `align_pair` on AVX2 including the per-chunk build.
///
/// Falls back to the per-read scalar path for any chunk whose `Database` fails to
/// build, so this can never change results — only how fast they are computed.
pub fn tso_clip_lens_cr4_batch(reads: &[&[u8]]) -> Vec<usize> {
    use rayon::prelude::*;

    reads
        .par_chunks(CR4_SCAN_CHUNK)
        .flat_map_iter(|chunk| {
            let targets: Vec<Vec<u8>> = chunk.iter().map(|r| cr4_target(r)).collect();
            let db = hyalite::Database::builder()
                .sequences(&targets)
                .scoring(cr4_scoring().clone())
                .mode(hyalite::Mode::Ov)
                .search_type(hyalite::SearchType::ScoreEnd)
                .max_query_len(TSO_SEQ.len())
                .build();

            let Ok(db) = db else {
                // Should not happen for this fixed shape; degrade to scalar rather
                // than change results.
                return chunk
                    .iter()
                    .map(|r| tso_clip_len_cr4(r))
                    .collect::<Vec<_>>();
            };

            let mut scratch = hyalite::Scratch::new(&db);
            let mut hits = Vec::with_capacity(chunk.len());
            db.scan_all(&mut scratch, tso_query(), &mut hits);
            // `scan_all` returns one hit per sequence, in order; the zip relies on it.
            debug_assert_eq!(hits.len(), chunk.len());

            chunk
                .iter()
                .zip(&hits)
                .map(|(read, hit)| {
                    if read.is_empty() {
                        0
                    } else {
                        cr4_gate(hit, read.len())
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Number of 3' bases to trim as a CellRanger4 polyA tail — STAR
/// `ClipCR4::polyTail3p`.
///
/// Walks in from the 3' end scoring +1 per `A` (encoded 0) and -2 per non-`A`,
/// remembering the longest prefix of that walk whose running score still clears
/// a 70% density threshold (`score * 10 >= ib * 7`). It gives up once the score
/// has fallen more than 27 behind the position, and returns nothing unless the
/// remembered score reached 20. Reads shorter than 20 bases are never trimmed.
///
/// Note the scan does **not** stop at the tail boundary, so A-rich sequence
/// upstream of a clean run legitimately extends the trim. That reads like an
/// off-by-one but is STAR's behaviour; the oracle pins it.
///
/// Ported from STAR; the Rust transcription is [Benjamin Demaille's in PR
/// #148](https://github.com/Psy-Fer/rustar-aligner/pull/148), validated against STAR's own compiled `polyTail3p` — see
/// `cr4_tso_matches_star_opal_oracle`.
pub(crate) fn poly_tail_3p(seq: &[u8]) -> usize {
    let seq_len = seq.len();
    if seq_len < 20 {
        return 0;
    }
    let mut best_len: i64 = seq_len as i64 - 1;
    let mut score: i64 = 0;
    let mut best_score: i64 = 0;
    for ib in 1..=seq_len as i64 {
        if seq[seq_len - ib as usize] == 0 {
            score += 1;
            if score * 10 >= ib * 7 {
                best_len = ib;
                best_score = score;
            }
        } else {
            score -= 2;
            if ib - score > 27 {
                break;
            }
        }
    }
    if best_score < 20 {
        0
    } else {
        best_len as usize
    }
}

/// Clip the 10x TSO from the 5' end and trim a 3' polyA tail of the cDNA read,
/// matching `--clipAdapterType CellRanger4`. Operates on encoded bases
/// (0=A..3=T,4=N) with parallel quality bytes. Returns the clipped read.
///
/// The 5' TSO clip is a STAR-faithful overlap alignment (see [`tso_clip_len_cr4`]);
/// the 3' polyA trim is STAR's scored scan (see [`poly_tail_3p`]).
/// Returns `(clipped_seq, clipped_qual, clip5p, clip3p)` — the CR4-clipped read plus
/// the bases trimmed from the 5' (TSO) and 3' (polyA) ends, so the caller can soft-clip
/// them (STARsolo keeps them in SEQ as soft-clips, e.g. `60M30S`, not dropped).
/// No fixed clips; runs go through [`crate::clip::clip_mate_cr4`].
pub fn clip_adapter_cr4(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let (c5, c3) = crate::clip::clip_mate_cr4(
        seq,
        &crate::clip::ClipParams::default(),
        tso_clip_len_cr4(seq),
    );
    let end = seq.len() - c3;
    (
        seq[c5..end].to_vec(),
        qual.get(c5..end.min(qual.len()))
            .map(<[u8]>::to_vec)
            .unwrap_or_default(),
        c5,
        c3,
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

/// `read_index` of a record from a run that does not need STAR's readInfo
/// (before [`SoloReadOutcome::set_read_index`] stamps it).
pub const NO_READ_INDEX: u32 = u32::MAX;

/// `gene` of a no-feature record ([`FeatureOutcome::no_feature`]), STAR's
/// feature `-1`.
pub const NO_FEATURE: u32 = u32::MAX;

/// The SJ feature's records: one per junction a uniquely mapped read crosses,
/// as STAR writes one line per junction (`outputReadCB`,
/// `SoloReadFeature_record.cpp:253-262`).
#[derive(Debug, Default)]
pub struct SjRecords {
    pub unique: Vec<SjCountRecord>,
    /// Records of 1MM_multi barcodes, resolved by the posterior when the matrix
    /// is built (`SoloReadFeature_inputRecords.cpp:66-104`).
    pub multi: Vec<SjMultiRecord>,
}

impl SjRecords {
    pub fn is_empty(&self) -> bool {
        self.unique.is_empty() && self.multi.is_empty()
    }

    pub fn append(&mut self, other: SjRecords) {
        self.unique.extend(other.unique);
        self.multi.extend(other.multi);
    }
}

/// An SJ record whose barcode matched several whitelist entries by 1MM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SjMultiRecord {
    pub candidates: Vec<CbCandidate>,
    pub umi: u64,
    pub intron_start: u64,
    pub intron_end: u64,
}

/// One (cell, UMI, splice-junction) observation for the `SJ` feature. The
/// junction is identified by its absolute intron coordinates; it is mapped to a
/// matrix row (the `SJ.out.tab` order) at output time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SjCountRecord {
    pub cb: u32,
    pub umi: u64,
    pub intron_start: u64,
    pub intron_end: u64,
    /// The read's barcode matched the whitelist exactly (`cbMatch` 0), for
    /// `yessubWLmatchExact`.
    pub exact: bool,
}

/// One read's record for the `Velocyto` feature (STAR writes `iRead nTr tr
/// type ...`, `SoloReadFeature_record.cpp:160-171`): the read, and the
/// transcripts its unique alignment is compatible with, each with STAR's type
/// bits ([`velocyto_transcript_types`]). Its cell and UMI come later, from the
/// Gene feature's readInfo.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VelocytoRecord {
    pub read_index: u32,
    pub types: Vec<(u32, u8)>,
    /// The whitelist entries the read's barcode matched (several for an
    /// unresolved 1MM_multi barcode): STAR counts each as a cell this feature
    /// detected (`cbReadCount`, `SoloReadFeature_record.cpp:187-189`), which
    /// decides the cells of the filtered matrices.
    pub cbs: Vec<u32>,
}

/// The Velocyto feature's records for the whole run, kept compact: STAR
/// writes them to a temporary file, and a run can have hundreds of millions.
/// Each read is its index and transcript count, and each transcript its index
/// and type bits packed in one `u32` (`tr << 4 | bits`; STAR's transcript
/// index is a `uint32`, and the type bits are the four of
/// [`AlignVsTranscript`](crate::quant::splice_status::AlignVsTranscript)).
/// The cells the records' barcodes matched are a bitset over the barcode
/// index.
#[derive(Debug, Default)]
pub struct VelocytoRecords {
    reads: Vec<(u32, u32)>,
    types: Vec<u32>,
    cells: Vec<u64>,
}

impl VelocytoRecords {
    pub fn push(&mut self, r: &VelocytoRecord) {
        self.reads.push((r.read_index, r.types.len() as u32));
        self.types.extend(r.types.iter().map(|&(tr, ty)| {
            debug_assert!(tr < 1 << 28 && ty < 1 << 4);
            tr << 4 | u32::from(ty)
        }));
        for &cb in &r.cbs {
            let word = cb as usize / 64;
            if word >= self.cells.len() {
                self.cells.resize(word + 1, 0);
            }
            self.cells[word] |= 1 << (cb % 64);
        }
    }

    /// Each read's index and its `(transcript, type bits)`, sorted by
    /// transcript.
    pub fn iter(&self) -> impl Iterator<Item = (u32, impl Iterator<Item = (u32, u8)> + '_)> {
        let mut start = 0;
        self.reads.iter().map(move |&(read_index, n)| {
            let packed = &self.types[start..start + n as usize];
            start += n as usize;
            (
                read_index,
                packed.iter().map(|&p| (p >> 4, (p & 0xf) as u8)),
            )
        })
    }

    /// Whether a record's barcode matched cell `cb`.
    pub fn detected(&self, cb: u32) -> bool {
        self.cells
            .get(cb as usize / 64)
            .is_some_and(|w| w & (1 << (cb % 64)) != 0)
    }
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
    /// The whitelist candidates of a 1MM_multi barcode, resolved at counting
    /// time as for a unique-gene read (`cb` is meaningless until then); empty
    /// when the barcode was resolved per read.
    pub candidates: Vec<CbCandidate>,
}

/// Thread-safe sink for the records produced during alignment.
#[derive(Default)]
pub struct SoloRecorder {
    pub records: Mutex<Vec<SoloCountRecord>>,
    pub multi_records: Mutex<Vec<SoloMultiRecord>>,
    /// Gene-ambiguous reads for `--soloMultiMappers` (resolved CB only).
    pub multi_gene: Mutex<Vec<MultiGeneRecord>>,
    /// No-feature reads, for readInfo only ([`FeatureOutcome::no_feature`]).
    pub no_feature: Mutex<Vec<SoloCountRecord>>,
    pub no_feature_multi: Mutex<Vec<SoloMultiRecord>>,
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
    /// Transcript models for the `Gene` feature and GX/GN (`classifyAlign`).
    pub gene_model: GeneModel,
    /// `Barcodes.stats` and the CB+UMI Q30 tally.
    pub stats: BarcodeStats,
    /// `Features.stats` per quantified feature, parallel to `features`.
    pub feature_stats: Vec<FeatureStats>,
    /// `Features.stats` of the SJ and Velocyto features.
    pub sj_stats: FeatureStats,
    pub velocyto_stats: FeatureStats,
    /// `(Q30 bases, bases)` of the cDNA reads' clipped sequence, for
    /// "Q30 Bases in RNA read" (`ReadAlign_oneRead.cpp:27-33`).
    pub rna_q30: [AtomicU64; 2],
    /// Quantified features (`Gene`, `GeneFull`, …), each with its own recorder
    /// and `Solo.out/<feature>/raw/` output. Parallel to `recorders`.
    pub features: Vec<SoloFeature>,
    pub recorders: Vec<SoloRecorder>,
    /// CellRanger-style positional mapping funnel over uniquely-mapped reads
    /// (independent of barcode), kept only for `--soloCellRangerSummary Yes`
    /// and written only when both `Gene` and `GeneFull` features run.
    pub region_stats: Option<RegionStats>,
    /// `--soloFeatures SJ`: collect per-cell splice-junction counts.
    pub sj_enabled: bool,
    /// (cell, UMI, junction) observations for the SJ feature.
    pub sj_records: Mutex<SjRecords>,
    /// `--soloFeatures Velocyto`: collect spliced/unspliced/ambiguous counts.
    pub velocyto_enabled: bool,
    /// The Velocyto feature's per-read records.
    pub velocyto_records: Mutex<VelocytoRecords>,
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
    ///
    /// One slot per feature, as STAR keeps a readInfo per feature: the first
    /// feature's for `CB`/`UB` (`samAttrFeature`), the Gene feature's for
    /// Velocyto (`ParametersSolo.cpp:420-422`).
    pub read_info: Vec<Option<Mutex<Vec<ReadInfo>>>>,
    /// The slot of `read_info` the `CB`/`UB` SAM tags read (the first
    /// feature's), when those tags were requested.
    pub tag_read_info_slot: Option<usize>,
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

/// STAR's `SoloReadFeatureStats` (`SoloReadFeatureStats.h`) for one feature,
/// written to `Solo.out/<feature>/Features.stats` and behind `Summary.csv`.
#[derive(Debug, Default)]
pub struct FeatureStats {
    counts: [AtomicU64; FeatureStats::NAMES.len()],
}

impl FeatureStats {
    /// STAR's names, in its order.
    pub const NAMES: [&'static str; 11] = [
        "noUnmapped",
        "noNoFeature",
        "MultiFeature",
        "subMultiFeatureMultiGenomic",
        "noTooManyWLmatches",
        "noMMtoWLwithoutExact",
        "yesWLmatch",
        "yessubWLmatchExact",
        "yessubWLmatch_UniqueFeature",
        "yesCellBarcodes",
        "yesUMIs",
    ];
    pub const NO_UNMAPPED: usize = 0;
    pub const NO_NO_FEATURE: usize = 1;
    pub const MULTI_FEATURE: usize = 2;
    pub const SUB_MULTI_FEATURE_MULTI_GENOMIC: usize = 3;
    pub const NO_TOO_MANY_WL_MATCHES: usize = 4;
    pub const NO_MM_TO_WL_WITHOUT_EXACT: usize = 5;
    pub const YES_WL_MATCH: usize = 6;
    pub const YES_SUB_WL_MATCH_EXACT: usize = 7;
    pub const YES_SUB_WL_MATCH_UNIQUE_FEATURE: usize = 8;
    pub const YES_CELL_BARCODES: usize = 9;
    pub const YES_UMIS: usize = 10;

    pub fn add(&self, i: usize, n: u64) {
        self.counts[i].fetch_add(n, Ordering::Relaxed);
    }

    pub fn get(&self, i: usize) -> u64 {
        self.counts[i].load(Ordering::Relaxed)
    }

    /// STAR's `numInvalidBarcodes`: reads dropped for their barcode at counting.
    pub fn invalid(&self) -> u64 {
        self.get(Self::NO_TOO_MANY_WL_MATCHES) + self.get(Self::NO_MM_TO_WL_WITHOUT_EXACT)
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
    /// `--soloFeatures SJ` and the read is uniquely mapped with a whitelisted
    /// barcode and a valid UMI.
    pub sj: SjRecords,
    /// Velocyto record for this read (resolved CB, gene-assigned), if enabled.
    pub velocyto: Option<VelocytoRecord>,
    /// Barcode facts for the per-read SAM tags. `None` when the barcode read was
    /// too short to extract a CB+UMI at all.
    pub barcode: Option<SoloBarcodeTags>,
}

impl SoloReadOutcome {
    /// Stamp the input read index onto every count record this read produced, so
    /// UMI collapsing can fill STAR's readInfo (for the `CB`/`UB` SAM tags and
    /// Velocyto).
    pub fn set_read_index(&mut self, read_index: u32) {
        for fo in &mut self.per_feature {
            if let Some(r) = &mut fo.record {
                r.read_index = read_index;
            }
            if let Some(m) = &mut fo.multi {
                m.read_index = read_index;
            }
            if let Some(r) = &mut fo.no_feature {
                r.read_index = read_index;
            }
            if let Some(m) = &mut fo.no_feature_multi {
                m.read_index = read_index;
            }
        }
        if let Some(v) = &mut self.velocyto {
            v.read_index = read_index;
        }
    }
}

/// One paired-end alignment as STAR holds it: a single `Transcript` with both
/// mates' blocks, the left mate first (`i_frag` 0) and the right mate after it
/// (`i_frag` 1), on mate 1's strand. STAR's read is mate 1 followed by mate 2
/// reverse-complemented, so mate 1 is the left mate on the forward strand and
/// mate 2 on the reverse strand. Only the blocks, CIGAR, strand and chromosome
/// are meaningful; this feeds gene annotation, not output.
pub fn pair_alignment(m1: &Transcript, m2: &Transcript) -> Transcript {
    let (left, right) = if m1.is_reverse { (m2, m1) } else { (m1, m2) };
    let mut exons = Vec::with_capacity(left.exons.len() + right.exons.len());
    for (frag, mate) in [(0, left), (1, right)] {
        exons.extend(mate.exons.iter().cloned().map(|mut e| {
            e.i_frag = frag;
            e
        }));
    }
    // The mates' own CIGARs, so `block_gaps` tells a deletion from a junction,
    // with a skip between them for the mate boundary (never consulted:
    // `classify_align` finds mate starts from `i_frag`).
    let mut cigar = left.cigar.clone();
    cigar.push(noodles::sam::alignment::record::cigar::Op::new(
        noodles::sam::alignment::record::cigar::op::Kind::Skip,
        1,
    ));
    cigar.extend(right.cigar.iter().copied());
    Transcript {
        exons,
        cigar,
        genome_start: m1.genome_start.min(m2.genome_start),
        genome_end: m1.genome_end.max(m2.genome_end),
        is_reverse: m1.is_reverse,
        ..m1.clone()
    }
}

/// The `sF` tag of one alignment (`ReadAlign_alignBAM.cpp:444-452`): the
/// read's `(ovType, number of genes)`, or `(-1, -1)` for an alignment with no
/// gene of its own when `ovType` is a sense type (exonic, exonic50p,
/// intronic: 1, 3, 5).
fn sf_tag(ov_type: i32, n_genes: i32, align_has_no_gene: bool) -> [i32; 2] {
    if matches!(ov_type, 1 | 3 | 5) && align_has_no_gene {
        [-1, -1]
    } else {
        [ov_type, n_genes]
    }
}

/// The gene tag values of one alignment of a read
/// ([`SoloContext::gene_tag_values`]).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignGeneTag {
    /// `GX`/`GN`: the read's single gene, `-` otherwise.
    pub gx_read: String,
    pub gn_read: String,
    /// `gx`/`gn`: this alignment's genes.
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
    /// STAR's no-feature record (`SoloReadFeature_record.cpp:180-182`): a read
    /// with a valid barcode that this feature does not count to a single gene,
    /// kept only to fill readInfo (for the `CB`/`UB` tags or Velocyto). `gene`
    /// is [`NO_FEATURE`]; a 1MM_multi barcode goes in `no_feature_multi`
    /// instead.
    pub no_feature: Option<SoloCountRecord>,
    pub no_feature_multi: Option<SoloMultiRecord>,
}

impl SoloContext {
    /// Build the solo context from parameters: load the whitelist and build the
    /// gene model from `--sjdbGTFfile`. Call once before alignment.
    pub fn build(
        params: &Parameters,
        genome: &crate::genome::Genome,
        transcriptome: Option<&crate::quant::transcriptome::TranscriptomeIndex>,
    ) -> Result<Self, Error> {
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

        // Gene model: --sjdbGTFfile if given, else the index's annotation tables.
        // `CB_samTagOut` quantifies nothing (`Solo.cpp:13` builds no
        // SoloFeature), so it runs without a gene model.
        let tag_out_only = params.solo_type == SoloType::CbSamTagOut;
        let (gene_ann, gene_model) = if tag_out_only {
            (GeneAnnotation::default(), GeneModel::default())
        } else {
            let gene_ann = crate::quant::resolve_gene_annotation(params, genome, transcriptome)?;
            log::info!("STARsolo: {} genes loaded", gene_ann.n_genes());
            // The Gene feature tests alignments against transcript models
            // (`Transcriptome::classifyAlign`). The index's transcriptome is the
            // mapping-time GTF's when one is given, else the index tables.
            let gene_model = GeneModel::new(
                transcriptome.cloned().ok_or_else(|| {
                    Error::Index(
                        "STARsolo needs transcript models: the index's transcriptInfo.tab \
                     or --sjdbGTFfile"
                            .into(),
                    )
                })?,
                &gene_ann,
            );
            (gene_ann, gene_model)
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
        let sj_enabled = params.solo_features.iter().any(|f| f == "SJ");
        let velocyto_enabled = params.solo_features.iter().any(|f| f == "Velocyto");
        let transcript3p = params.solo_features.iter().any(|f| f == "Transcript3p");
        let want_multi = params.solo_multi_mappers.iter().any(|m| m != "Unique");

        // STAR's readInfo, per feature; sized once the read count is known
        // (`reserve_read_info`).
        let read_info = features
            .iter()
            .enumerate()
            .map(|(fi, &f)| {
                let tags = fi == 0 && params.solo_read_info_needed();
                let velocyto = velocyto_enabled && f == SoloFeature::Gene;
                (tags || velocyto).then(|| Mutex::new(Vec::new()))
            })
            .collect();
        Ok(Self {
            layout: SoloBarcodeLayout::from_params(params),
            whitelist,
            match_type: params.solo_cb_match_type(),
            strand,
            gene_ann,
            gene_model,
            stats: BarcodeStats::new(),
            feature_stats: features.iter().map(|_| FeatureStats::default()).collect(),
            sj_stats: FeatureStats::default(),
            velocyto_stats: FeatureStats::default(),
            rna_q30: [AtomicU64::new(0), AtomicU64::new(0)],
            features,
            recorders,
            region_stats: params.solo_cellranger_summary().then(RegionStats::default),
            sj_enabled,
            sj_records: Mutex::new(SjRecords::default()),
            velocyto_enabled,
            velocyto_records: Mutex::new(VelocytoRecords::default()),
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
                .then(|| -> Result<_, Error> {
                    // GTF at mapping time wins (STAR); else the index tables.
                    if let Some(gtf_path) = params.sjdb_gtf_file.as_ref() {
                        crate::quant::transcriptome::TranscriptomeIndex::from_mapping_gtf(
                            params, gtf_path, genome,
                        )
                    } else {
                        transcriptome.cloned().ok_or_else(|| {
                            Error::Index(
                                "Solo transcript 3' features need the index's transcriptInfo.tab \
                                 or --sjdbGTFfile"
                                    .into(),
                            )
                        })
                    }
                })
                .transpose()?,
            transcript3p: transcript3p
                .then(|| Mutex::new(crate::solo::transcript3p::Transcript3pAcc::new())),
            read_info,
            tag_read_info_slot: params.solo_read_info_needed().then_some(0),
        })
    }

    /// Whether this run has to track STAR's readInfo (for the `CB`/`UB` SAM
    /// tags or Velocyto).
    pub fn read_info_enabled(&self) -> bool {
        self.read_info.iter().any(Option::is_some)
    }

    /// The readInfo the `CB`/`UB` tags are filled from (the first feature's),
    /// when those tags were requested.
    pub fn tag_read_info(&self) -> Option<&Mutex<Vec<ReadInfo>>> {
        self.read_info.get(self.tag_read_info_slot?)?.as_ref()
    }

    /// The Gene feature's readInfo, which Velocyto takes its cells and UMIs
    /// from.
    pub fn gene_read_info(&self) -> Option<&Mutex<Vec<ReadInfo>>> {
        let fi = self.features.iter().position(|&f| f == SoloFeature::Gene)?;
        self.read_info[fi].as_ref()
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
        self.record_barcode_code(cb_match.star_code(), &cb_match);
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
        for info in self.read_info.iter().flatten() {
            let mut info = info.lock().unwrap();
            info.clear();
            info.resize(n_reads, ReadInfo::default());
        }
    }

    /// The feature the gene SAM tags are computed from: STAR's `samAttrFeature`,
    /// the first entry of `--soloFeatures` (`ParametersSolo.cpp:423`).
    fn tag_feature(&self) -> SoloFeature {
        self.features.first().copied().unwrap_or(SoloFeature::Gene)
    }

    /// The gene tags of one read, one entry per alignment in `alignments`
    /// order, each a mate pair for paired-end (see [`pair_alignment`]).
    ///
    /// They come from the tag feature's read annotation (`ReadAnnotFeature`):
    /// `Gene` from `classifyAlign` (whose `ovType` is exonic when it found a
    /// gene, else none, `Transcriptome_classifyAlign.cpp:262-263`); `GeneFull`
    /// from gene-body overlap (`Transcriptome_geneFullAlignOverlap.cpp`, which
    /// sets no `ovType`). Then, as `ReadAlign_alignBAM.cpp:99-106, 420-452`:
    /// - `GX`/`GN`: the alignment's gene when the read's set and the
    ///   alignment's each hold exactly one gene, else `-`;
    /// - `gx`/`gn`: the alignment's genes, `;`-joined, `-` when none;
    /// - `sF`: `(ovType, genes of the read)`, or `(-1, -1)` for an alignment
    ///   with no gene when `ovType` is a sense type (1, 3 or 5).
    pub fn gene_tag_values(&self, alignments: &[Transcript]) -> Vec<AlignGeneTag> {
        let (sets, ov_type) = match self.tag_feature() {
            SoloFeature::GeneFull => (gene_full_sets(alignments, &self.gene_ann, self.strand), 0),
            SoloFeature::Gene => {
                let sets = classify_align(alignments, &self.gene_model, self.strand);
                let ov = i32::from(!sets.genes.is_empty());
                (sets, ov)
            }
        };
        let n_genes = sets.genes.len() as i32;
        let join = |genes: &[u32], names: &[String]| {
            if genes.is_empty() {
                "-".to_string()
            } else {
                genes
                    .iter()
                    .map(|&g| names[g as usize].as_str())
                    .collect::<Vec<_>>()
                    .join(";")
            }
        };
        sets.per_align
            .iter()
            .map(|genes| {
                let unique = match (sets.genes.len(), genes.as_slice()) {
                    (1, &[g]) => Some(g as usize),
                    _ => None,
                };
                AlignGeneTag {
                    gx_read: unique
                        .map_or("-", |g| &self.gene_ann.gene_ids[g])
                        .to_string(),
                    gn_read: unique
                        .map_or("-", |g| &self.gene_ann.gene_names[g])
                        .to_string(),
                    gx: join(genes, &self.gene_ann.gene_ids),
                    gn: join(genes, &self.gene_ann.gene_names),
                    sf: sf_tag(ov_type, n_genes, genes.is_empty()),
                }
            })
            .collect()
    }

    /// Count a read with no barcode in `Barcodes.stats`: a missing adapter or
    /// UMI directly (`SoloReadBarcode_getCBandUMI.cpp:340-350`), a missing CB
    /// as `cbMatch` -11.
    pub fn record_barcode_miss(&self, miss: BarcodeMiss) {
        match miss {
            BarcodeMiss::NoAdapter => self.stats.no_adapter(),
            BarcodeMiss::NoUmi => self.stats.no_umi(),
            BarcodeMiss::NoCb => self.add_stats(-11),
        }
    }

    /// STAR's `addStats`: a read's bin by its final `cbMatch`, none without a
    /// whitelist (`SoloReadBarcode_getCBandUMI.cpp:95-96`).
    fn add_stats(&self, code: i32) {
        if !matches!(self.whitelist, CbWhitelist::NoWhitelist { .. }) {
            self.stats.record(code);
        }
    }

    /// Add one read's quality tallies: its barcode read's (`barcode_q30`) and
    /// its cDNA mates' clipped qualities (STAR's `qualHist`, bases at or above
    /// `readQualityScoreBase + 30`).
    pub fn record_qualities(&self, barcode_q30: (u64, u64), mates: &[&[u8]]) {
        self.stats.add_q30(barcode_q30.0, barcode_q30.1);
        for q in mates {
            let q30 = q.iter().filter(|&&b| b >= 33 + 30).count() as u64;
            self.rna_q30[0].fetch_add(q30, Ordering::Relaxed);
            self.rna_q30[1].fetch_add(q.len() as u64, Ordering::Relaxed);
        }
    }

    /// Count a read by its final `cbMatch` code (`addStats`), and for an exact
    /// match with a valid UMI (code 0), its barcode's exact count.
    fn record_barcode_code(&self, code: i32, cb_match: &CbMatch) {
        self.add_stats(code);
        if code == 0
            && let CbMatch::Exact(idx) = cb_match
        {
            self.whitelist.bump_exact(*idx);
        }
    }

    /// STAR's per-read feature counts for a read with a valid barcode and UMI
    /// (`SoloReadFeature_record.cpp:75-176`), plus `yessubWLmatchExact` for
    /// its counted records (`SoloReadFeature_inputRecords.cpp:111-113`; an
    /// exact barcode always survives the later filters).
    fn record_feature_stats(
        &self,
        n_loci: usize,
        class: &ReadClass,
        out: &SoloReadOutcome,
        exact: bool,
        junctions: &[(u64, u64)],
    ) {
        use FeatureStats as F;
        for (fi, (&feature, st)) in self.features.iter().zip(&self.feature_stats).enumerate() {
            if n_loci == 0 {
                st.add(F::NO_UNMAPPED, 1);
                continue;
            }
            let assignment = match feature {
                SoloFeature::Gene => class.gene,
                SoloFeature::GeneFull => class.gene_full,
            };
            match assignment {
                GeneAssignment::Gene(_) => {
                    if exact {
                        st.add(F::YES_SUB_WL_MATCH_EXACT, 1);
                    }
                }
                GeneAssignment::Ambiguous => {
                    st.add(F::MULTI_FEATURE, 1);
                    if n_loci > 1 {
                        st.add(F::SUB_MULTI_FEATURE_MULTI_GENOMIC, 1);
                    }
                    let counted = out
                        .per_feature
                        .get(fi)
                        .is_some_and(|fo| fo.multi_gene.is_some());
                    if exact && counted {
                        st.add(F::YES_SUB_WL_MATCH_EXACT, 1);
                    }
                }
                GeneAssignment::NoFeature | GeneAssignment::Unmapped => {
                    st.add(F::NO_NO_FEATURE, 1);
                }
            }
        }
        if self.sj_enabled {
            let st = &self.sj_stats;
            if n_loci == 0 {
                st.add(F::NO_UNMAPPED, 1);
            } else if n_loci > 1 {
                st.add(F::SUB_MULTI_FEATURE_MULTI_GENOMIC, 1);
                st.add(F::MULTI_FEATURE, 1);
            } else if junctions.is_empty() {
                st.add(F::NO_NO_FEATURE, 1);
            }
            // The rest are counted per junction record that has a matrix row,
            // when the matrix is built (`SoloReadFeature_inputRecords.cpp:23-27`).
        }
        // Velocyto: STAR's `record` counts its reads too, but `countVelocyto`
        // never sums the per-thread stats into the feature's (only
        // `countCBgeneUMI` calls `addStats`, `SoloFeature_countCBgeneUMI.cpp:47`),
        // so Velocyto's Features.stats carries only the UMIs and cells.
    }

    /// Process one solo read: match the cell barcode, validate the UMI, assign
    /// a gene, and (on success) produce a count record. Stats are recorded
    /// here; the returned records are appended to the recorder by the caller.
    pub fn process_read(
        &self,
        cdna_transcripts: &[Transcript],
        n_loci: usize,
        barcode: Result<&CellBarcode, BarcodeMiss>,
        junctions: &[(u64, u64)],
    ) -> SoloReadOutcome {
        let mut out = SoloReadOutcome::default();

        // One-pass classification: the two overlap queries are shared between the
        // per-feature gene assignment and the CellRanger-style mapping funnel, so
        // this is no more work than the old per-feature `assign_gene_se` calls.
        let want_exon = self.features.contains(&SoloFeature::Gene);
        let want_body = self.features.contains(&SoloFeature::GeneFull);
        let class = classify_read(
            cdna_transcripts,
            &self.gene_ann,
            &self.gene_model,
            self.strand,
            want_exon,
            want_body,
            self.want_multi,
        );

        // Mapping funnel: count uniquely-mapped reads by region (CellRanger's
        // "confidently mapped" = MAPQ 255 ≈ a single alignment), independent of
        // barcode validity. `n_loci` is the number of genomic loci the read (or
        // pair) maps to — for SE this equals `cdna_transcripts.len()`.
        if let Some(funnel) = &self.region_stats
            && n_loci == 1
        {
            match class.region {
                Some(Region::Exonic) => {
                    funnel.exonic.fetch_add(1, Ordering::Relaxed);
                }
                Some(Region::Intronic) => {
                    funnel.intronic.fetch_add(1, Ordering::Relaxed);
                }
                Some(Region::Intergenic) => {
                    funnel.intergenic.fetch_add(1, Ordering::Relaxed);
                }
                None => {}
            }
            if class.antisense {
                funnel.antisense.fetch_add(1, Ordering::Relaxed);
            }
        }

        // No barcode → nothing to count.
        let bc = match barcode {
            Ok(bc) => bc,
            Err(miss) => {
                self.record_barcode_miss(miss);
                return out;
            }
        };

        // Cell-barcode match.
        let cb_match = self
            .whitelist
            .match_cb(&bc.cb_seq, &bc.cb_qual, self.match_type);

        // Per-read SAM-tag facts, recorded before any early return so that a
        // rejected barcode still tags its alignments. STAR reports a UMI
        // rejection in place of the CB code (`getCBandUMI.cpp:304`).
        let umi_status = check_umi(&bc.umi_seq);
        let code = umi_status
            .star_code()
            .unwrap_or_else(|| cb_match.star_code());
        self.record_barcode_code(code, &cb_match);
        out.barcode = Some(SoloBarcodeTags {
            cb_match: code,
            cb_index: cb_match.resolved_index(),
        });

        let cb_resolved: Option<u32> = match &cb_match {
            CbMatch::Exact(idx) | CbMatch::Corrected(idx) => Some(*idx),
            CbMatch::Multi(_) => None, // deferred to collation
            CbMatch::NoMatch | CbMatch::NinCb | CbMatch::MultMatchRejected => return out,
        };

        // UMI validity.
        let UmiCheck::Ok(umi) = umi_status else {
            return out;
        };

        // SJ feature: record (cell, UMI, junction) for each crossed junction;
        // a 1MM_multi barcode's cell is resolved when the matrix is built.
        if self.sj_enabled && !junctions.is_empty() {
            match (cb_resolved, &cb_match) {
                (Some(cb), _) => {
                    out.sj.unique = junctions
                        .iter()
                        .map(|&(intron_start, intron_end)| SjCountRecord {
                            cb,
                            umi,
                            intron_start,
                            intron_end,
                            exact: code == 0,
                        })
                        .collect();
                }
                (None, CbMatch::Multi(cands)) => {
                    out.sj.multi = junctions
                        .iter()
                        .map(|&(intron_start, intron_end)| SjMultiRecord {
                            candidates: cands.clone(),
                            umi,
                            intron_start,
                            intron_end,
                        })
                        .collect();
                }
                (None, _) => {}
            }
        }

        // Velocyto: the transcripts a uniquely mapped read is compatible with
        // (`classifyAlign` computes them only for `nAlignG==1`), recorded for
        // any read whose barcode and UMI pass (`record` returns before
        // `cbMatch<0`), counted or not; its cell and UMI come from the Gene
        // feature's readInfo when the molecules are collapsed.
        if self.velocyto_enabled
            && n_loci == 1
            && let [aln] = cdna_transcripts
        {
            let types = velocyto_transcript_types(aln, &self.gene_model, self.strand);
            if !types.is_empty() {
                let cbs = match &cb_match {
                    CbMatch::Exact(i) | CbMatch::Corrected(i) => vec![*i],
                    CbMatch::Multi(cands) => cands.iter().map(|c| c.wl_index).collect(),
                    CbMatch::NoMatch | CbMatch::NinCb | CbMatch::MultMatchRejected => Vec::new(),
                };
                out.velocyto = Some(VelocytoRecord {
                    read_index: NO_READ_INDEX,
                    types,
                    cbs,
                });
            }
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
                // STAR's no-feature record, for the features whose readInfo
                // is kept (the first for `CB`/`UB`, Gene for Velocyto).
                let no_feature = |fo: &mut FeatureOutcome| {
                    if self.read_info.get(fi).is_none_or(Option::is_none) {
                        return;
                    }
                    match (cb_resolved, &cb_match) {
                        (Some(cb), _) => {
                            fo.no_feature = Some(SoloCountRecord {
                                cb,
                                umi,
                                gene: NO_FEATURE,
                                read_index: NO_READ_INDEX,
                            });
                        }
                        (None, CbMatch::Multi(cands)) => {
                            fo.no_feature_multi = Some(SoloMultiRecord {
                                candidates: cands.clone(),
                                umi,
                                gene: NO_FEATURE,
                                read_index: NO_READ_INDEX,
                            });
                        }
                        (None, _) => {}
                    }
                };
                let gene = match assignment {
                    GeneAssignment::Gene(g) => g,
                    GeneAssignment::Ambiguous => {
                        // Gene-ambiguous read: record its gene set for
                        // --soloMultiMappers distribution. A 1MM_multi barcode
                        // goes along with its candidates: STAR records these
                        // reads too and resolves the barcode when it counts
                        // (`SoloReadFeature_inputRecords.cpp`).
                        let genes = match feature {
                            SoloFeature::Gene => &class.gene_multi,
                            SoloFeature::GeneFull => &class.gene_full_multi,
                        };
                        let cb_and_candidates = match (cb_resolved, &cb_match) {
                            (Some(cb), _) => Some((cb, Vec::new())),
                            (None, CbMatch::Multi(cands)) => Some((0, cands.clone())),
                            (None, _) => None,
                        };
                        if let Some((cb, candidates)) = cb_and_candidates
                            && !genes.is_empty()
                        {
                            fo.multi_gene = Some(MultiGeneRecord {
                                cb,
                                umi,
                                genes: genes.clone(),
                                candidates,
                            });
                        }
                        // readInfo: the resolved cell and the raw UMI, whether
                        // the read is left uncounted or goes to the multi-gene
                        // matrices (`SoloFeature_collapseUMIall.cpp:279-284`).
                        no_feature(&mut fo);
                        return fo;
                    }
                    GeneAssignment::NoFeature | GeneAssignment::Unmapped => {
                        no_feature(&mut fo);
                        return fo;
                    }
                };
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

        // The read's `cbMatch` is >= 0 here: STAR's `SoloReadFeature::record`
        // counts it for every feature.
        self.record_feature_stats(n_loci, &class, &out, code == 0, junctions);
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
        barcode: Result<&CellBarcode, BarcodeMiss>,
        junctions: &[(u64, u64)],
    ) -> SoloReadOutcome {
        // A mate pair is one alignment in STAR (one `Transcript`, both mates'
        // blocks), so the Gene feature asks one transcript to contain the pair.
        let alns: Vec<Transcript> = pairs
            .iter()
            .map(|(m1, m2)| pair_alignment(m1, m2))
            .collect();
        self.process_read(&alns, pairs.len(), barcode, junctions)
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
            file_index: 0,
        }
    }

    fn enc(s: &str) -> Vec<u8> {
        s.bytes().map(encode_base).collect()
    }

    /// `sF` as STAR writes it: `Gene` gives `ovType` 1 when the read has a
    /// gene, else 0; `GeneFull` always 0.
    #[test]
    fn sf_tag_follows_star() {
        assert_eq!(sf_tag(1, 2, false), [1, 2]);
        // A multimapper with a gene elsewhere but none on this alignment.
        assert_eq!(sf_tag(1, 1, true), [-1, -1]);
        // No gene anywhere (Gene), and GeneFull, whose ovType is never set.
        assert_eq!(sf_tag(0, 0, true), [0, 0]);
        assert_eq!(sf_tag(0, 1, true), [0, 1]);
    }

    #[test]
    fn cr4_tso_clip_removes_exact_tso_prefix() {
        // Read = exact TSO (30 nt) + 61 nt of cDNA (all C, so the 3' polyA trim
        // is a no-op). Overlap of the 30-nt TSO against target[0..30] scores 30,
        // ends at target index 29 ⇒ L = 30, well above STAR's gate.
        let cdna: String = std::iter::repeat_n('C', 61).collect();
        let read = format!("AAGCAGTGGTATCAACGCAGAGTACATGGG{cdna}");
        let seq = enc(&read);
        assert_eq!(seq.len(), 91);
        assert_eq!(tso_clip_len_cr4(&seq), 30);

        let qual = vec![b'I'; seq.len()];
        let (clipped, _q, c5, c3) = clip_adapter_cr4(&seq, &qual);
        assert_eq!(c5, 30, "5' TSO clipped");
        assert_eq!(c3, 0, "no polyA tail");
        assert_eq!(clipped, enc(&cdna));
    }

    /// `clip_adapter_cr4` must run [`poly_tail_3p`] on the read *after* the 5'
    /// TSO clip (STAR clips type 10 before type 11), and must use STAR's scored
    /// scan rather than the old "trailing run of A >= 8" rule.
    #[test]
    fn cr4_polya_trim_follows_star_not_a_run_length_rule() {
        let qual = |n: usize| vec![b'I'; n];

        // 10 clean A's: the old >=8 run rule trimmed these; STAR does not,
        // because polyTail3p's best score (10) never reaches the 20 floor.
        let seq = enc(&format!("{}{}", "CGT".repeat(20), "A".repeat(10)));
        let (_, _, c5, c3) = clip_adapter_cr4(&seq, &qual(seq.len()));
        assert_eq!(
            (c5, c3),
            (0, 0),
            "a 10 base tail is below STAR's score floor"
        );

        // 25 clean A's clears the floor and is trimmed whole.
        let seq = enc(&format!("{}{}", "CGT".repeat(20), "A".repeat(25)));
        let (_, _, c5, c3) = clip_adapter_cr4(&seq, &qual(seq.len()));
        assert_eq!((c5, c3), (0, 25));

        // An interrupted tail survives at 70% density — the run-length rule
        // would have stopped at the first non-A and trimmed only 15.
        let seq = enc(&format!(
            "{}{}C{}",
            "CGT".repeat(20),
            "A".repeat(15),
            "A".repeat(15)
        ));
        let (_, _, _, c3) = clip_adapter_cr4(&seq, &qual(seq.len()));
        assert_eq!(c3, 31, "scored scan spans the interruption");

        // Both halves together: the polyA scan sees the post-TSO-clip read, so
        // the 5' and 3' clips are independent and both apply.
        let cdna = "CGT".repeat(15); // 45 nt, A-free
        let seq = enc(&format!(
            "AAGCAGTGGTATCAACGCAGAGTACATGGG{cdna}{}",
            "A".repeat(30)
        ));
        let (clipped, _, c5, c3) = clip_adapter_cr4(&seq, &qual(seq.len()));
        assert_eq!((c5, c3), (30, 30), "TSO off the 5', polyA off the 3'");
        assert_eq!(clipped, enc(&cdna));
    }

    #[test]
    fn cr4_tso_clip_is_noop_without_adapter() {
        // Homopolymer read with no TSO similarity ⇒ score below STAR's gate ⇒ 0.
        assert_eq!(tso_clip_len_cr4(&[3u8; 91]), 0);
    }

    #[test]
    fn cr4_tso_clip_handles_short_and_empty_reads() {
        assert_eq!(tso_clip_len_cr4(&[]), 0);
        // Shorter than the TSO: target is N-padded to 91; no qualifying match.
        assert_eq!(tso_clip_len_cr4(&enc("ACGTACGT")), 0);
    }

    /// Backend-forcing twin of [`tso_clip_lens_cr4_batch`], so the determinism
    /// test can pin a specific backend (the production function uses runtime
    /// dispatch and offers no such knob).
    fn tso_clip_lens_batch(reads: &[Vec<u8>], backend: hyalite::BackendChoice) -> Vec<usize> {
        use hyalite::{Database, Mode, Scratch, SearchType};
        let targets: Vec<Vec<u8>> = reads
            .iter()
            .map(|r| {
                let take = r.len().min(CR4_READ_LEN);
                let mut t = r[..take].to_vec();
                t.resize(CR4_READ_LEN, 4);
                t
            })
            .collect();
        let db = Database::builder()
            .sequences(&targets)
            .scoring(cr4_scoring().clone())
            .mode(Mode::Ov)
            .search_type(SearchType::ScoreEnd)
            .max_query_len(TSO_SEQ.len())
            .backend(backend)
            .build()
            .expect("build cr4 database");
        let mut scratch = Scratch::new(&db);
        let mut out = Vec::new();
        db.scan_all(&mut scratch, tso_query(), &mut out);
        reads
            .iter()
            .zip(&out)
            .map(|(r, hit)| {
                let s = hit.score;
                let l = hit.target_end.map_or(0, |e| e + 1) as i32;
                let reject = s < 20 || (s == 20 && l > 26) || (s == 21 && l > 30);
                if reject { 0 } else { (l as usize).min(r.len()) }
            })
            .collect()
    }

    /// A deterministic, diverse CR4 read set (TSO with mutations/shifts + random),
    /// all length 91, for exercising the batched path in-tree without fixtures.
    fn synthetic_cr4_reads() -> Vec<Vec<u8>> {
        let mut st = 0x51ED_2701_u64;
        let mut lcg = || {
            st = st
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            st >> 33
        };
        let tso = tso_query().to_vec();
        let mut reads = Vec::new();
        // 700 reads in all, more than one CR4_SCAN_CHUNK, so the batch tests
        // cross a chunk boundary.
        for _ in 0..550 {
            let mut r = tso.clone();
            let k = (lcg() % 13) as usize; // 0..12 mismatches
            for _ in 0..k {
                let p = (lcg() as usize) % r.len();
                r[p] = ((r[p] as u64 + 1 + lcg() % 3) % 4) as u8;
            }
            let shift = (lcg() % 6) as usize;
            let mut read: Vec<u8> = (0..shift).map(|_| (lcg() % 4) as u8).collect();
            read.extend(r);
            while read.len() < CR4_READ_LEN {
                read.push((lcg() % 4) as u8);
            }
            read.truncate(CR4_READ_LEN);
            reads.push(read);
        }
        for _ in 0..150 {
            reads.push((0..CR4_READ_LEN).map(|_| (lcg() % 4) as u8).collect());
        }
        assert!(reads.len() > CR4_SCAN_CHUNK);
        reads
    }

    /// Confirms the patched per-target `scan_all` (SIMD) path serves CR4: it
    /// agrees with the scalar `align_pair` clip we ship, is identical across
    /// backends, and (when the oracle fixture is provided) matches STAR's Opal.
    #[test]
    fn cr4_tso_scan_all_batch_matches_scalar_and_is_backend_stable() {
        let reads = synthetic_cr4_reads();
        let scalar_ref: Vec<usize> = reads.iter().map(|r| tso_clip_len_cr4(r)).collect();

        // Batched scan_all on Scalar must equal the per-read scalar path.
        let batch_scalar = tso_clip_lens_batch(
            &reads,
            hyalite::BackendChoice::Force(hyalite::Backend::Scalar),
        );
        assert_eq!(batch_scalar, scalar_ref, "scan_all(Scalar) vs align_pair");

        // And every available SIMD backend must match it bit-for-bit.
        for b in [hyalite::Backend::Sse41, hyalite::Backend::Avx2] {
            if b.is_available() {
                let batch = tso_clip_lens_batch(&reads, hyalite::BackendChoice::Force(b));
                assert_eq!(batch, scalar_ref, "scan_all({b:?}) diverged");
            }
        }

        // The production batch entry point (runtime-dispatched backend, rayon
        // chunking) must land on the same clips as the per-read scalar path.
        let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        assert_eq!(
            tso_clip_lens_cr4_batch(&refs),
            scalar_ref,
            "tso_clip_lens_cr4_batch vs align_pair"
        );
    }

    /// Faithfulness gate: our CR4 clip — both halves — vs STAR's own C++.
    ///
    /// `tests/data/cr4_opal_oracle.tsv` (columns: read, clip, score, L, polyA)
    /// is generated by `tests/data/cr4_opal_oracle.cpp`, which links STAR's own
    /// `opal.cpp` and reproduces `ClipCR4` + `ClipMate::clipChunk` verbatim.
    /// 900 reads straddle the 5' gate `S<20 / S==20&&L>26 / S==21&&L>30` (TSO
    /// with 0-14 mismatches, 5' shifts, truncations, indels, embedded Ns, polyA,
    /// random, short reads); a further 38 straddle `polyTail3p`'s own boundaries
    /// (the score>=20 floor, the 70% density rule, the `ib-score>27` give-up,
    /// the `seqLen<20` early return, and A-rich sequence upstream of the tail).
    ///
    /// The polyA column is computed on the read *after* the 5' clip, matching
    /// `ClipMate::clip`'s ordering. Committing the expected output lets this run
    /// in CI with no Opal build; regenerate with the recipe at the top of the
    /// `.cpp` and override the path via `RUSTAR_CR4_ORACLE`.
    #[test]
    fn cr4_tso_matches_star_opal_oracle() {
        let path = std::env::var("RUSTAR_CR4_ORACLE").unwrap_or_else(|_| {
            format!(
                "{}/tests/data/cr4_opal_oracle.tsv",
                env!("CARGO_MANIFEST_DIR")
            )
        });
        let data = std::fs::read_to_string(&path).expect("read oracle tsv");
        let mut reads: Vec<Vec<u8>> = Vec::new();
        let mut expected: Vec<usize> = Vec::new();
        let mut expected_polya: Vec<usize> = Vec::new();
        for line in data.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            reads.push(f[0].bytes().map(encode_base).collect());
            expected.push(f[1].parse().expect("clip int"));
            expected_polya.push(f[4].parse().expect("polyA int"));
        }
        assert!(!reads.is_empty(), "oracle file was empty");

        // 3' polyA half: STAR's polyTail3p on the post-5'-clip read.
        let mut polya_mism = 0usize;
        for (i, read) in reads.iter().enumerate() {
            let got = poly_tail_3p(&read[expected[i].min(read.len())..]);
            if got != expected_polya[i] {
                eprintln!("POLYA MISMATCH star={} ours={got}", expected_polya[i]);
                polya_mism += 1;
            }
        }
        assert_eq!(
            polya_mism,
            0,
            "{polya_mism}/{} reads disagreed with STAR's polyTail3p",
            reads.len()
        );

        // Scalar per-read path (what we ship).
        let scalar: Vec<usize> = reads.iter().map(|r| tso_clip_len_cr4(r)).collect();
        // Batched SIMD per-target path (scan_all), on the best available backend.
        let backend = if hyalite::Backend::Avx2.is_available() {
            hyalite::Backend::Avx2
        } else if hyalite::Backend::Sse41.is_available() {
            hyalite::Backend::Sse41
        } else {
            hyalite::Backend::Scalar
        };
        let batch = tso_clip_lens_batch(&reads, hyalite::BackendChoice::Force(backend));
        // The production entry point, exactly as the solo read loop calls it.
        let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        let prod = tso_clip_lens_cr4_batch(&refs);

        let mut mism = 0usize;
        for i in 0..reads.len() {
            if scalar[i] != expected[i] || batch[i] != expected[i] || prod[i] != expected[i] {
                eprintln!(
                    "MISMATCH star={} scalar={} simd={} prod={}",
                    expected[i], scalar[i], batch[i], prod[i]
                );
                mism += 1;
            }
        }
        assert_eq!(
            mism,
            0,
            "{mism}/{} reads disagreed with STAR's Opal (scalar, {backend:?}, or production batch)",
            reads.len()
        );
    }

    /// Kicks hyalite's headline "bit-identical across backends" guarantee on the
    /// CR4 workload shape, and cross-checks the SIMD `Database` path against the
    /// scalar `align_pair` path we actually ship in `tso_clip_len_cr4`.
    ///
    /// NB: `Database::scan` returns the single best hit over the DB (with
    /// `db_index`), not a per-sequence result, so this validates the best-hit
    /// only. STAR's CR4 needs a per-read result, which is what
    /// [`tso_clip_lens_cr4_batch`] gets from `scan_all` (0.2.0); this test keeps
    /// the single-best `scan` path covered alongside it.
    #[test]
    fn cr4_tso_database_is_deterministic_across_backends_and_matches_scalar() {
        use hyalite::{Backend, BackendChoice, Database, Mode, Scratch, SearchType};

        // Reads-as-targets, TSO-as-query (STAR's orientation → target_end is the
        // read coordinate). Read 0 is TSO-led, so it is the unambiguous winner.
        let mut reads: Vec<Vec<u8>> = Vec::new();
        let mut tso_led = tso_query().to_vec();
        tso_led.extend(std::iter::repeat_n(1u8, CR4_READ_LEN - tso_led.len())); // + C's → 91
        reads.push(tso_led);
        reads.push(enc(
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
        )); // 91
        reads.push(vec![3u8; CR4_READ_LEN]); // homopolymer T

        let mut backends = vec![Backend::Scalar];
        for b in [Backend::Sse41, Backend::Avx2, Backend::Neon] {
            if b.is_available() {
                backends.push(b);
            }
        }

        let mut per_backend: Vec<(i32, usize, usize)> = Vec::new(); // (score, db_index, target_end)
        for backend in &backends {
            let Ok(db) = Database::builder()
                .sequences(&reads)
                .scoring(cr4_scoring().clone())
                .mode(Mode::Ov)
                .search_type(SearchType::ScoreEnd)
                .max_query_len(TSO_SEQ.len())
                .backend(BackendChoice::Force(*backend))
                .build()
            else {
                continue; // backend not applicable to this DB; skip
            };
            let mut scratch = Scratch::new(&db);
            let hit = db.scan(&mut scratch, tso_query());
            per_backend.push((hit.score, hit.db_index, hit.target_end.unwrap()));
        }

        // Every available backend agrees, bit for bit.
        assert_ne!(per_backend.len(), 0);
        let first = per_backend[0];
        for (i, r) in per_backend.iter().enumerate() {
            assert_eq!(*r, first, "backend {:?} diverged", backends[i]);
        }

        // The SIMD/Database best hit matches the scalar align_pair path on the
        // winning read — same score, same end coordinate.
        let (score, db_index, target_end) = first;
        let scalar = hyalite::align_pair(
            tso_query(),
            &reads[db_index],
            cr4_scoring(),
            Mode::Ov,
            SearchType::ScoreEnd,
        )
        .unwrap();
        assert_eq!(scalar.score, score);
        assert_eq!(scalar.target_end.unwrap(), target_end);

        // And that winner is the TSO-led read, clipped at the full TSO.
        assert_eq!(db_index, 0);
        assert_eq!(target_end + 1, TSO_SEQ.len());
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
    fn extract_too_short_is_no_cb() {
        let layout = v2_layout();
        let read = encoded_read("short", "AAAAAAAACCCC", "IIIIIIIIIIII");
        assert_eq!(layout.extract(&read), Err(BarcodeMiss::NoCb));
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
            file_index: 0,
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
        assert!(
            err.contains("the total length of barcode sequence is 20 not equal to expected 26"),
            "{err}"
        );
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
            &[cdna.path().to_path_buf()],
            &[bc.path().to_path_buf()],
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
            &[cdna.path().to_path_buf()],
            &[bc.path().to_path_buf()],
            v2_layout(),
            None,
            false,
            BarcodeReadLength::default(),
        )
        .unwrap();
        assert!(reader.read_batch(10).is_err());
    }

    #[test]
    fn adapter_anchoring_finds_the_adapter_and_positions_around_it() {
        // Adapter ACGT sits at offset 4. CB is the 4 bases before it
        // (anchor 2, the adapter start), UMI the 4 after (anchor 3, adapter
        // end). Both move with the adapter, which is the whole point.
        let adapter = vec![0u8, 1, 2, 3];
        let layout = SoloBarcodeLayout::Anchored {
            cb_segments: vec![(
                Anchored {
                    anchor: 2,
                    dist: -4,
                },
                Anchored {
                    anchor: 2,
                    dist: -1,
                },
            )],
            umi: (
                Anchored { anchor: 3, dist: 1 },
                Anchored { anchor: 3, dist: 4 },
            ),
            adapter: adapter.clone(),
            max_mismatches: 1,
        };

        // GGGG ACGT TTTT  →  CB = GGGG, UMI = TTTT
        let read = EncodedRead {
            name: "r".into(),
            sequence: vec![2, 2, 2, 2, 0, 1, 2, 3, 3, 3, 3, 3],
            quality: vec![b'I'; 12],
            file_index: 0,
        };
        let bc = layout.extract(&read).expect("adapter present");
        assert_eq!(bc.cb_seq, vec![2, 2, 2, 2]);
        assert_eq!(bc.umi_seq, vec![3, 3, 3, 3]);

        // Shift the adapter two bases right; both fields follow it.
        let read = EncodedRead {
            name: "r".into(),
            sequence: vec![1, 1, 2, 2, 2, 2, 0, 1, 2, 3, 3, 3, 3, 3],
            quality: vec![b'I'; 14],
            file_index: 0,
        };
        let bc = layout.extract(&read).expect("adapter present, shifted");
        assert_eq!(bc.cb_seq, vec![2, 2, 2, 2]);
        assert_eq!(bc.umi_seq, vec![3, 3, 3, 3]);

        // One mismatch inside the adapter is within budget.
        let read = EncodedRead {
            name: "r".into(),
            sequence: vec![2, 2, 2, 2, 0, 1, 2, 0, 3, 3, 3, 3],
            quality: vec![b'I'; 12],
            file_index: 0,
        };
        assert!(layout.extract(&read).is_ok());

        // No adapter at all: no barcode. A barcode read at an unknown offset
        // is not salvageable, so guessing would be worse than declining.
        let read = EncodedRead {
            name: "r".into(),
            sequence: vec![2; 12],
            quality: vec![b'I'; 12],
            file_index: 0,
        };
        assert_eq!(layout.extract(&read), Err(BarcodeMiss::NoAdapter));
    }

    /// The discriminating case: an early offset that fits the budget against a
    /// later one that fits it better. STAR scores every offset and keeps the
    /// best, so the later exact match wins. Taking the first acceptable one
    /// instead would anchor the barcode four bases off.
    #[test]
    fn adapter_search_takes_the_best_scoring_offset_not_the_leftmost() {
        // ACGT at 0 with one mismatch, ACGT exactly at 5.
        let seq = [0u8, 1, 2, 0, 9, 0, 1, 2, 3];
        assert_eq!(local_align_hamming_dist(&seq, &[0, 1, 2, 3]), (0, 5));
    }

    /// Among equally good offsets STAR keeps the first, because it replaces its
    /// best only on a strictly smaller distance.
    #[test]
    fn adapter_search_breaks_ties_leftmost() {
        let seq = [0u8, 1, 2, 3, 9, 0, 1, 2, 3];
        assert_eq!(local_align_hamming_dist(&seq, &[0, 1, 2, 3]), (0, 0));
    }

    /// An `N` in the adapter is a wildcard in STAR's query, not a base that
    /// fails to match, so it costs nothing against any read base.
    #[test]
    fn an_n_in_the_adapter_is_not_a_mismatch() {
        assert_eq!(
            local_align_hamming_dist(&[0u8, 1, 2, 3], &[0, 4, 2, 3]),
            (0, 0)
        );
    }

    /// STAR's starting values: the distance is the adapter length at offset 0
    /// until an offset does better, and a read shorter than the adapter gives
    /// its length + 1. So with a budget of at least the adapter length an
    /// adapter matching nowhere is still "found" at offset 0, as in STAR.
    #[test]
    fn adapter_search_keeps_stars_starting_values() {
        assert_eq!(local_align_hamming_dist(&[2u8; 8], &[0, 1, 2, 3]), (3, 0));
        assert_eq!(local_align_hamming_dist(&[2u8; 8], &[0, 0, 0, 0]), (4, 0));
        assert_eq!(local_align_hamming_dist(&[0u8, 1], &[0, 1, 2, 3]), (3, 0));
    }

    /// STAR's `extractBarcode` bounds: the end may equal the read length (its
    /// `substr` clamps), but not exceed it; a negative start or a start after
    /// the end is refused.
    #[test]
    fn barcode_span_follows_extract_barcode() {
        assert_eq!(barcode_span(2, 5, 10), Some((2, 4)));
        assert_eq!(barcode_span(8, 10, 10), Some((8, 2)));
        assert_eq!(barcode_span(8, 11, 10), None);
        assert_eq!(barcode_span(-1, 3, 10), None);
        assert_eq!(barcode_span(4, 3, 10), None);
    }

    /// STAR drops a read whose adapter is not found within the budget
    /// (`noNoAdapter`), even when no position is measured from the adapter.
    #[test]
    fn no_adapter_means_no_barcode_whatever_the_anchors() {
        let layout = SoloBarcodeLayout::Anchored {
            cb_segments: vec![(
                Anchored { anchor: 0, dist: 0 },
                Anchored { anchor: 0, dist: 3 },
            )],
            umi: (
                Anchored { anchor: 0, dist: 4 },
                Anchored { anchor: 0, dist: 7 },
            ),
            adapter: vec![0, 1, 2, 3],
            max_mismatches: 0,
        };
        let read = EncodedRead {
            name: "r".into(),
            sequence: vec![2; 12],
            quality: vec![b'I'; 12],
            file_index: 0,
        };
        assert_eq!(layout.extract(&read), Err(BarcodeMiss::NoAdapter));
    }

    #[test]
    fn anchors_outside_the_read_yield_no_barcode() {
        // Anchor 2 with a distance that runs off the front.
        let layout = SoloBarcodeLayout::Anchored {
            cb_segments: vec![(
                Anchored {
                    anchor: 2,
                    dist: -99,
                },
                Anchored {
                    anchor: 2,
                    dist: -1,
                },
            )],
            umi: (
                Anchored { anchor: 3, dist: 1 },
                Anchored { anchor: 3, dist: 2 },
            ),
            adapter: vec![0, 1, 2, 3],
            max_mismatches: 0,
        };
        let read = EncodedRead {
            name: "r".into(),
            sequence: vec![0, 1, 2, 3, 3, 3, 3],
            quality: vec![b'I'; 7],
            file_index: 0,
        };
        assert_eq!(layout.extract(&read), Err(BarcodeMiss::NoCb));
    }
}
