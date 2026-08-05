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
pub mod smartseq;
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
    debug_assert!(matches!(
        params.solo_type,
        SoloType::CbUmiSimple | SoloType::CbUmiComplex
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
    SoloReadReader::open(cdna, barcode, layout, params.read_files_command.as_deref())
}

/// One paired-end solo read for `--soloBarcodeMate 1` (5' 10x): both mates carry
/// cDNA, and the cell barcode (CB+UMI) is a prefix of mate 1.
pub struct SoloPairedRead {
    pub mate1: EncodedRead,
    pub mate2: EncodedRead,
    /// `None` when mate 1 was too short to extract CB+UMI.
    pub barcode: Option<CellBarcode>,
}

/// Reads the two cDNA mate files in lockstep for a `--soloBarcodeMate 1` run,
/// extracting the barcode from the start of mate 1.
pub struct SoloPairedReader {
    mate1: FastqReader,
    mate2: FastqReader,
    layout: SoloBarcodeLayout,
}

impl SoloPairedReader {
    pub fn open(
        mate1_path: &Path,
        mate2_path: &Path,
        layout: SoloBarcodeLayout,
        decompress_cmd: Option<&str>,
    ) -> Result<Self, Error> {
        Ok(Self {
            mate1: FastqReader::open(mate1_path, decompress_cmd)?,
            mate2: FastqReader::open(mate2_path, decompress_cmd)?,
            layout,
        })
    }

    /// Fetch the next (mate1, mate2) pair with the barcode extracted from mate 1.
    pub fn next_read(&mut self) -> Result<Option<SoloPairedRead>, Error> {
        match (self.mate1.next_encoded()?, self.mate2.next_encoded()?) {
            (Some(mate1), Some(mate2)) => {
                let barcode = self.layout.extract(&mate1);
                Ok(Some(SoloPairedRead {
                    mate1,
                    mate2,
                    barcode,
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
            "solo: --soloBarcodeMate 1 requires two --readFilesIn cDNA mate files",
        ))
    })?;
    let layout = SoloBarcodeLayout::from_params(params);
    SoloPairedReader::open(mate1, mate2, layout, params.read_files_command.as_deref())
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
/// #148](https://github.com/Psy-Fer/rustar-aligner/pull/148), moved here onto
/// the solo path (the only place `--clipAdapterType CellRanger4` is reached)
/// and validated against STAR's own compiled `polyTail3p` — see
/// `cr4_tso_matches_star_opal_oracle`.
fn poly_tail_3p(seq: &[u8]) -> usize {
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
pub fn clip_adapter_cr4(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>, usize, usize) {
    // 5' TSO: STAR-faithful overlap alignment (ClipCR4's Opal call → hyalite).
    clip_adapter_cr4_with_tso(seq, qual, tso_clip_len_cr4(seq))
}

/// [`clip_adapter_cr4`] with the 5' TSO clip length already computed — the entry
/// point for the batched scan, which resolves every read's clip up front via
/// [`tso_clip_lens_cr4_batch`].
pub fn clip_adapter_cr4_with_tso(
    seq: &[u8],
    qual: &[u8],
    start: usize,
) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let start = start.min(seq.len());

    // 3' polyA, on the read as it stands after the 5' clip — STAR clips type 10
    // (which shifts the sequence and shrinks Lread) before type 11.
    //
    // NB: with a non-zero --clip5pNbases/--clip3pNbases alongside CellRanger4,
    // STAR's window differs (it applies those fixed trims first). Narrowing the
    // window to match made agreement *worse*, because that pairing has a larger
    // pre-existing +/-clip5pNbases position bug that swamps it — see TODO.md.
    // Left matching the default (fixed clips = 0) case, which is validated.
    let end = seq.len() - poly_tail_3p(&seq[start..]);

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

        // Quantified gene features (Gene, GeneFull). Validation guarantees these
        // parse; default to Gene if somehow empty.
        let features: Vec<SoloFeature> = params
            .solo_features
            .iter()
            .filter_map(|f| f.parse().ok())
            .collect();
        let features = if features.is_empty() {
            vec![SoloFeature::Gene]
        } else {
            features
        };
        let recorders = features.iter().map(|_| SoloRecorder::new()).collect();
        let feature_reads = features.iter().map(|_| AtomicU64::new(0)).collect();
        let sj_enabled = params.solo_features.iter().any(|f| f == "SJ");
        let velocyto_enabled = params.solo_features.iter().any(|f| f == "Velocyto");
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
        })
    }

    /// Process one solo read: match the cell barcode, validate the UMI, assign
    /// a gene, and (on success) produce a count record. Stats are recorded
    /// here; the returned records are appended to the recorder by the caller.
    /// GX/GN SAM tag strings for a read's `Gene`-feature assignment:
    /// `(gene_id, gene_name)` when uniquely assigned, else `("-", "-")`
    /// (STARsolo convention). Drives `--outSAMattributes GX GN`.
    pub fn gene_tags<'a>(&'a self, transcripts: &[Transcript]) -> (&'a str, &'a str) {
        match assign_gene_se(transcripts, &self.gene_ann, self.strand, SoloFeature::Gene) {
            GeneAssignment::Gene(g) => (
                self.gene_ann.gene_ids[g as usize].as_str(),
                self.gene_ann.gene_names[g as usize].as_str(),
            ),
            _ => ("-", "-"),
        }
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
                    (Some(cb), _) => fo.record = Some(SoloCountRecord { cb, umi, gene }),
                    (None, CbMatch::Multi(cands)) => {
                        fo.multi = Some(SoloMultiRecord {
                            candidates: cands.clone(),
                            umi,
                            gene,
                        });
                    }
                    (None, _) => unreachable!("non-multi unresolved CB returned early"),
                }
                fo
            })
            .collect();
        out
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

    fn enc(s: &str) -> Vec<u8> {
        s.bytes().map(encode_base).collect()
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
        for _ in 0..400 {
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
        for _ in 0..100 {
            reads.push((0..CR4_READ_LEN).map(|_| (lcg() % 4) as u8).collect());
        }
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
        assert!(!per_backend.is_empty());
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
