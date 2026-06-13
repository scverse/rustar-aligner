//! Cell-barcode whitelist loading and read-stage CB/UMI matching (Phase 14.2).
//!
//! Faithful port of STAR's `SoloReadBarcode_getCBandUMI.cpp` read stage:
//! barcodes are 2-bit packed (seq[0] in the high bits) into a `u64` and the
//! whitelist is a sorted array searched by binary search. Exact match,
//! single-N correction, and 1-mismatch (1MM / 1MM_multi) correction follow
//! STAR's enumeration exactly.
//!
//! The 1MM_multi *posterior* resolution (count + quality weighted) is a
//! collation-stage concern and is deferred to Phase 14.4 — here a multi-match
//! read records all candidate whitelist indices plus the mismatch-position
//! quality, exactly as STAR's `cbMatchString`.

use crate::error::Error;
use crate::io::fastq::{decode_base, encode_base};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;
use std::sync::atomic::{AtomicU64, Ordering};

/// Maximum barcode length representable in a `u64` (32 × 2-bit bases).
pub const CB_LEN_MAX: usize = 32;

// ---------------------------------------------------------------------------
// Barcode packing
// ---------------------------------------------------------------------------

/// Result of packing an encoded barcode into a `u64`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackResult {
    /// No ambiguous bases; the packed value.
    NoN(u64),
    /// Exactly one `N`; `packed` has `A` (0) at the N position.
    OneN { packed: u64, pos: usize },
    /// More than one `N` — uncorrectable.
    ManyN,
}

/// 2-bit pack an encoded barcode (`0=A,1=C,2=G,3=T,4=N`) with `seq[0]` in the
/// high bits, matching STAR's `convertNuclStrToInt64`.
pub fn pack_barcode(seq: &[u8]) -> PackResult {
    let len = seq.len();
    let mut packed: u64 = 0;
    let mut n_pos: Option<usize> = None;
    let mut n_count = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        let shift = 2 * (len - 1 - i);
        if b >= 4 {
            n_count += 1;
            if n_count > 1 {
                return PackResult::ManyN;
            }
            n_pos = Some(i);
            // leave 0 (A) at this position; correction substitutes all 4 bases
        } else {
            packed |= (b as u64) << shift;
        }
    }
    match n_pos {
        None => PackResult::NoN(packed),
        Some(pos) => PackResult::OneN { packed, pos },
    }
}

/// Unpack a `u64` of `len` 2-bit bases back to an ASCII `ACGT` string
/// (`seq[0]` from the high bits).
pub fn unpack_barcode(packed: u64, len: usize) -> String {
    (0..len)
        .map(|i| {
            let shift = 2 * (len - 1 - i);
            decode_base(((packed >> shift) & 0b11) as u8) as char
        })
        .collect()
}

/// Bit shift for the base at sequence index `pos` in a `len`-base packing.
#[inline]
fn shift_for(pos: usize, len: usize) -> u32 {
    (2 * (len - 1 - pos)) as u32
}

// ---------------------------------------------------------------------------
// Match-type configuration (--soloCBmatchWLtype)
// ---------------------------------------------------------------------------

/// Flags decoded from `--soloCBmatchWLtype`. Mirrors STAR's `CBmatchWL`
/// boolean fields one-for-one, so the multiple bools are intentional.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct CbMatchType {
    /// Allow a single mismatch to the whitelist.
    pub mm1: bool,
    /// Keep multiple 1MM candidates for posterior resolution.
    pub mm1_multi: bool,
    /// Allow multiple matches for the N-substitution path.
    pub mm1_multi_nbase: bool,
    /// Add pseudocounts in posterior resolution (collation stage).
    pub pseudocounts: bool,
}

impl FromStr for CbMatchType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Exact" => Ok(Self {
                mm1: false,
                mm1_multi: false,
                mm1_multi_nbase: false,
                pseudocounts: false,
            }),
            "1MM" => Ok(Self {
                mm1: true,
                mm1_multi: false,
                mm1_multi_nbase: false,
                pseudocounts: false,
            }),
            "1MM_multi" => Ok(Self {
                mm1: true,
                mm1_multi: true,
                mm1_multi_nbase: false,
                pseudocounts: false,
            }),
            "1MM_multi_pseudocounts" => Ok(Self {
                mm1: true,
                mm1_multi: true,
                mm1_multi_nbase: false,
                pseudocounts: true,
            }),
            "1MM_multi_Nbase_pseudocounts" => Ok(Self {
                mm1: true,
                mm1_multi: true,
                mm1_multi_nbase: true,
                pseudocounts: true,
            }),
            _ => Err(format!(
                "unknown soloCBmatchWLtype '{s}'; expected Exact, 1MM, 1MM_multi, 1MM_multi_pseudocounts, or 1MM_multi_Nbase_pseudocounts"
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// Match result
// ---------------------------------------------------------------------------

/// One candidate whitelist barcode reachable by a single edit, plus the quality
/// of the mismatched base (for posterior resolution at collation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbCandidate {
    /// Index into the sorted whitelist.
    pub wl_index: u32,
    /// 0-based mismatch position in the read barcode.
    pub mismatch_pos: usize,
    /// Raw Phred+33 quality byte at the mismatch position.
    pub mismatch_qual: u8,
}

/// Outcome of matching one cell barcode to the whitelist. The negative STAR
/// `cbMatch` codes map to the rejection variants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CbMatch {
    /// Exact whitelist hit (cbMatch=0); carries the sorted whitelist index.
    Exact(u32),
    /// Unambiguous single-edit correction (cbMatch=1).
    Corrected(u32),
    /// Multiple 1MM candidates kept for later posterior resolution (cbMatch>1).
    Multi(Vec<CbCandidate>),
    /// No whitelist match within one edit (cbMatch=-1).
    NoMatch,
    /// More than one `N` in the barcode (cbMatch=-2).
    NinCb,
    /// >1 whitelist match but `mm1_multi` not enabled (cbMatch=-3).
    MultMatchRejected,
}

// ---------------------------------------------------------------------------
// UMI validity (matches STAR umiCheck=-23 / -24)
// ---------------------------------------------------------------------------

/// Outcome of validating a UMI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiCheck {
    /// Valid UMI; carries the packed value.
    Ok(u64),
    /// Contains an `N` (cbMatch=-23).
    NinUmi,
    /// Exact homopolymer, e.g. all-A (cbMatch=-24).
    Homopolymer,
}

/// Validate a UMI: reject any `N`, then reject exact homopolymers.
pub fn check_umi(umi_seq: &[u8]) -> UmiCheck {
    match pack_barcode(umi_seq) {
        PackResult::ManyN | PackResult::OneN { .. } => UmiCheck::NinUmi,
        PackResult::NoN(packed) => {
            if is_homopolymer(umi_seq) {
                UmiCheck::Homopolymer
            } else {
                UmiCheck::Ok(packed)
            }
        }
    }
}

fn is_homopolymer(seq: &[u8]) -> bool {
    match seq.first() {
        None => false,
        Some(&first) => seq.iter().all(|&b| b == first),
    }
}

// ---------------------------------------------------------------------------
// Whitelist
// ---------------------------------------------------------------------------

/// Cell-barcode whitelist. `List` is an explicit, sorted, de-duplicated set of
/// packed barcodes; `NoWhitelist` corresponds to `--soloCBwhitelist None`.
pub enum CbWhitelist {
    List {
        /// Sorted unique packed barcodes (binary-search target).
        sorted: Vec<u64>,
        /// `orig_index[k]` = line number of `sorted[k]` in the whitelist file,
        /// for `barcodes.tsv` column ordering (Phase 14.4).
        orig_index: Vec<u32>,
        /// Per-whitelist-index exact-match read counts (posterior prior).
        exact_counts: Vec<AtomicU64>,
        /// Barcode length in bases.
        len: usize,
    },
    /// `--soloCBwhitelist None`: keep every valid (N-free) barcode as observed.
    NoWhitelist { len: usize },
}

impl CbWhitelist {
    /// Number of whitelist barcodes (0 for `NoWhitelist`).
    pub fn len(&self) -> usize {
        match self {
            Self::List { sorted, .. } => sorted.len(),
            Self::NoWhitelist { .. } => 0,
        }
    }

    /// True if the whitelist has no barcodes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Barcode length in bases.
    pub fn barcode_len(&self) -> usize {
        match self {
            Self::List { len, .. } | Self::NoWhitelist { len } => *len,
        }
    }

    /// Decode the whitelist barcode at sorted index `idx` to an ASCII string.
    pub fn barcode_string(&self, idx: u32) -> Option<String> {
        match self {
            Self::List { sorted, len, .. } => {
                sorted.get(idx as usize).map(|&p| unpack_barcode(p, *len))
            }
            Self::NoWhitelist { .. } => None,
        }
    }

    /// Load a whitelist from a file (plain or gzip). One barcode per line;
    /// blank lines ignored. Barcodes are encoded, packed, sorted, de-duplicated.
    pub fn load(path: &Path) -> Result<Self, Error> {
        let reader = open_maybe_gzip(path)?;
        let mut packed: Vec<u64> = Vec::new();
        let mut len: usize = 0;
        for (lineno, line) in reader.lines().enumerate() {
            let line = line.map_err(Error::from)?;
            let bc = line.trim();
            if bc.is_empty() {
                continue;
            }
            // STARsolo whitelists may carry a second column (e.g. translated
            // barcodes for multi-ome); take the first whitespace token.
            let bc = bc.split_whitespace().next().unwrap_or("");
            if bc.is_empty() {
                continue;
            }
            if len == 0 {
                len = bc.len();
                if len == 0 || len > CB_LEN_MAX {
                    return Err(Error::from(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("whitelist barcode length {len} out of range (1..={CB_LEN_MAX})"),
                    )));
                }
            } else if bc.len() != len {
                return Err(Error::from(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "whitelist barcode on line {} has length {} (expected {len})",
                        lineno + 1,
                        bc.len()
                    ),
                )));
            }
            let encoded: Vec<u8> = bc.bytes().map(encode_base).collect();
            match pack_barcode(&encoded) {
                PackResult::NoN(p) => packed.push(p),
                _ => {
                    return Err(Error::from(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("whitelist barcode '{bc}' on line {} contains N", lineno + 1),
                    )));
                }
            }
        }
        if packed.is_empty() {
            return Err(Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "whitelist is empty",
            )));
        }
        // Sort by packed value, carrying the original line index; de-duplicate.
        let mut indexed: Vec<(u64, u32)> = packed
            .into_iter()
            .enumerate()
            .map(|(i, p)| (p, i as u32))
            .collect();
        indexed.sort_unstable_by_key(|&(p, _)| p);
        indexed.dedup_by_key(|&mut (p, _)| p);
        let sorted: Vec<u64> = indexed.iter().map(|&(p, _)| p).collect();
        let orig_index: Vec<u32> = indexed.iter().map(|&(_, i)| i).collect();
        let exact_counts = (0..sorted.len()).map(|_| AtomicU64::new(0)).collect();
        Ok(Self::List {
            sorted,
            orig_index,
            exact_counts,
            len,
        })
    }

    /// Binary-search the sorted whitelist for `packed`; returns the sorted index.
    fn search(&self, packed: u64) -> Option<u32> {
        match self {
            Self::List { sorted, .. } => sorted.binary_search(&packed).ok().map(|i| i as u32),
            Self::NoWhitelist { .. } => None,
        }
    }

    /// Increment the exact-match count for sorted whitelist index `idx`.
    fn bump_exact(&self, idx: u32) {
        if let Self::List { exact_counts, .. } = self {
            exact_counts[idx as usize].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Snapshot of exact-match counts per sorted whitelist index (for the
    /// Phase 14.4 posterior). Empty for `NoWhitelist`.
    pub fn exact_count_snapshot(&self) -> Vec<u64> {
        match self {
            Self::List { exact_counts, .. } => exact_counts
                .iter()
                .map(|c| c.load(Ordering::Relaxed))
                .collect(),
            Self::NoWhitelist { .. } => Vec::new(),
        }
    }

    /// Match one cell barcode against the whitelist following STAR's read stage.
    ///
    /// `cb_seq` is encoded (`0..=4`); `cb_qual` is raw Phred+33 (parallel to
    /// `cb_seq`). On an exact hit the whitelist's exact-count is incremented.
    pub fn match_cb(&self, cb_seq: &[u8], cb_qual: &[u8], match_type: CbMatchType) -> CbMatch {
        let len = cb_seq.len();
        match self {
            Self::NoWhitelist { .. } => match pack_barcode(cb_seq) {
                // No whitelist: every N-free barcode is its own "cell". We
                // cannot return a stable index without a whitelist, so callers
                // treat NoWhitelist specially; report NoMatch for N-containing.
                PackResult::NoN(_) => CbMatch::Exact(0),
                _ => CbMatch::NinCb,
            },
            Self::List { .. } => match pack_barcode(cb_seq) {
                PackResult::ManyN => CbMatch::NinCb,
                PackResult::NoN(packed) => {
                    if let Some(idx) = self.search(packed) {
                        self.bump_exact(idx);
                        return CbMatch::Exact(idx);
                    }
                    if !match_type.mm1 {
                        return CbMatch::NoMatch;
                    }
                    // 1MM: every position × the 3 alternate bases.
                    let mut candidates: Vec<CbCandidate> = Vec::new();
                    for pos in 0..len {
                        let shift = shift_for(pos, len);
                        let orig = (packed >> shift) & 0b11;
                        for alt in 0u64..4 {
                            if alt == orig {
                                continue;
                            }
                            let cand = (packed & !(0b11 << shift)) | (alt << shift);
                            if let Some(idx) = self.search(cand) {
                                candidates.push(CbCandidate {
                                    wl_index: idx,
                                    mismatch_pos: pos,
                                    mismatch_qual: qual_at(cb_qual, pos),
                                });
                            }
                        }
                    }
                    Self::resolve(candidates, match_type.mm1_multi)
                }
                PackResult::OneN { packed, pos } => {
                    if !match_type.mm1 {
                        return CbMatch::NoMatch;
                    }
                    // Substitute all 4 bases at the single N position.
                    let shift = shift_for(pos, len);
                    let mut candidates: Vec<CbCandidate> = Vec::new();
                    for base in 0u64..4 {
                        let cand = (packed & !(0b11 << shift)) | (base << shift);
                        if let Some(idx) = self.search(cand) {
                            candidates.push(CbCandidate {
                                wl_index: idx,
                                mismatch_pos: pos,
                                mismatch_qual: qual_at(cb_qual, pos),
                            });
                        }
                    }
                    Self::resolve(candidates, match_type.mm1_multi_nbase)
                }
            },
        }
    }

    /// Turn a candidate list into a [`CbMatch`], honoring the multi flag.
    fn resolve(candidates: Vec<CbCandidate>, allow_multi: bool) -> CbMatch {
        match candidates.len() {
            0 => CbMatch::NoMatch,
            1 => CbMatch::Corrected(candidates[0].wl_index),
            _ => {
                if allow_multi {
                    CbMatch::Multi(candidates)
                } else {
                    CbMatch::MultMatchRejected
                }
            }
        }
    }
}

#[inline]
fn qual_at(qual: &[u8], pos: usize) -> u8 {
    qual.get(pos).copied().unwrap_or(b'!') // '!' = Phred 0
}

/// Open a file, transparently decompressing `.gz`.
fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>, Error> {
    let file = File::open(path).map_err(|e| Error::io(e, path))?;
    let is_gz = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    if is_gz {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

// ---------------------------------------------------------------------------
// Stats (STAR cbMatch categories)
// ---------------------------------------------------------------------------

/// Per-run barcode-matching statistics, mirroring STAR's `SoloReadBarcodeStats`.
#[derive(Debug, Default)]
pub struct CbMatchStats {
    pub yes_exact: AtomicU64,
    pub yes_one_mm: AtomicU64,
    pub yes_mult_mm: AtomicU64,
    pub no_match: AtomicU64,
    pub n_in_cb: AtomicU64,
    pub mult_rejected: AtomicU64,
    pub n_in_umi: AtomicU64,
    pub umi_homopolymer: AtomicU64,
}

impl CbMatchStats {
    pub fn new() -> Self {
        Self::default()
    }

    /// Record one CB match outcome.
    pub fn record_cb(&self, m: &CbMatch) {
        let c = match m {
            CbMatch::Exact(_) => &self.yes_exact,
            CbMatch::Corrected(_) => &self.yes_one_mm,
            CbMatch::Multi(_) => &self.yes_mult_mm,
            CbMatch::NoMatch => &self.no_match,
            CbMatch::NinCb => &self.n_in_cb,
            CbMatch::MultMatchRejected => &self.mult_rejected,
        };
        c.fetch_add(1, Ordering::Relaxed);
    }

    /// Record one UMI check outcome (only the rejection cases are counted).
    pub fn record_umi(&self, u: &UmiCheck) {
        match u {
            UmiCheck::NinUmi => {
                self.n_in_umi.fetch_add(1, Ordering::Relaxed);
            }
            UmiCheck::Homopolymer => {
                self.umi_homopolymer.fetch_add(1, Ordering::Relaxed);
            }
            UmiCheck::Ok(_) => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn enc(s: &str) -> Vec<u8> {
        s.bytes().map(encode_base).collect()
    }

    fn write_wl(barcodes: &[&str]) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        for b in barcodes {
            writeln!(f, "{b}").unwrap();
        }
        f.flush().unwrap();
        f
    }

    #[test]
    fn pack_roundtrip() {
        let s = "ACGTACGT";
        match pack_barcode(&enc(s)) {
            PackResult::NoN(p) => assert_eq!(unpack_barcode(p, 8), s),
            _ => panic!("should pack cleanly"),
        }
    }

    #[test]
    fn pack_detects_one_and_many_n() {
        assert!(matches!(
            pack_barcode(&enc("ACNT")),
            PackResult::OneN { pos: 2, .. }
        ));
        assert_eq!(pack_barcode(&enc("ANNT")), PackResult::ManyN);
    }

    #[test]
    fn exact_match_and_count() {
        let f = write_wl(&["AAAA", "ACGT", "TTTT"]);
        let wl = CbWhitelist::load(f.path()).unwrap();
        let t = CbMatchType::from_str("1MM_multi").unwrap();
        let m = wl.match_cb(&enc("ACGT"), b"IIII", t);
        match m {
            CbMatch::Exact(idx) => assert_eq!(wl.barcode_string(idx).unwrap(), "ACGT"),
            other => panic!("expected exact, got {other:?}"),
        }
        let counts = wl.exact_count_snapshot();
        assert_eq!(counts.iter().sum::<u64>(), 1);
    }

    #[test]
    fn single_mismatch_correction() {
        let f = write_wl(&["AAAA", "ACGT", "TTTT"]);
        let wl = CbWhitelist::load(f.path()).unwrap();
        let t = CbMatchType::from_str("1MM").unwrap();
        // ACGA differs from ACGT at last position only.
        let m = wl.match_cb(&enc("ACGA"), b"IIII", t);
        match m {
            CbMatch::Corrected(idx) => assert_eq!(wl.barcode_string(idx).unwrap(), "ACGT"),
            other => panic!("expected corrected, got {other:?}"),
        }
    }

    #[test]
    fn ambiguous_multi_match_behavior() {
        // AAAA and CAAA both within 1MM of NAAA-ish read "GAAA"? Use TAAA read:
        // candidates AAAA (pos0 T->A) and CAAA (pos0 T->C). Both in WL.
        let f = write_wl(&["AAAA", "CAAA"]);
        let wl = CbWhitelist::load(f.path()).unwrap();

        // 1MM (no multi): rejected as ambiguous.
        let rej = wl.match_cb(&enc("TAAA"), b"IIII", CbMatchType::from_str("1MM").unwrap());
        assert_eq!(rej, CbMatch::MultMatchRejected);

        // 1MM_multi: both candidates kept for later resolution.
        let multi = wl.match_cb(
            &enc("TAAA"),
            b"IIII",
            CbMatchType::from_str("1MM_multi").unwrap(),
        );
        match multi {
            CbMatch::Multi(c) => assert_eq!(c.len(), 2),
            other => panic!("expected multi, got {other:?}"),
        }
    }

    #[test]
    fn no_match_when_too_far() {
        let f = write_wl(&["AAAA", "TTTT"]);
        let wl = CbWhitelist::load(f.path()).unwrap();
        let t = CbMatchType::from_str("1MM_multi").unwrap();
        // GGGG is >1 edit from both.
        assert_eq!(wl.match_cb(&enc("GGGG"), b"IIII", t), CbMatch::NoMatch);
    }

    #[test]
    fn n_correction_single() {
        let f = write_wl(&["AAAA", "ACGT"]);
        let wl = CbWhitelist::load(f.path()).unwrap();
        let t = CbMatchType::from_str("1MM_multi").unwrap();
        // ACGN → only ACGT matches among the 4 substitutions.
        let m = wl.match_cb(&enc("ACGN"), b"IIII", t);
        match m {
            CbMatch::Corrected(idx) => assert_eq!(wl.barcode_string(idx).unwrap(), "ACGT"),
            other => panic!("expected corrected, got {other:?}"),
        }
    }

    #[test]
    fn many_n_rejected() {
        let f = write_wl(&["AAAA"]);
        let wl = CbWhitelist::load(f.path()).unwrap();
        let t = CbMatchType::from_str("1MM_multi").unwrap();
        assert_eq!(wl.match_cb(&enc("NNAA"), b"IIII", t), CbMatch::NinCb);
    }

    #[test]
    fn exact_only_mode_no_correction() {
        let f = write_wl(&["ACGT"]);
        let wl = CbWhitelist::load(f.path()).unwrap();
        let t = CbMatchType::from_str("Exact").unwrap();
        assert_eq!(wl.match_cb(&enc("ACGA"), b"IIII", t), CbMatch::NoMatch);
    }

    #[test]
    fn umi_checks() {
        assert!(matches!(check_umi(&enc("ACGTAC")), UmiCheck::Ok(_)));
        assert_eq!(check_umi(&enc("ACGTNC")), UmiCheck::NinUmi);
        assert_eq!(check_umi(&enc("AAAAAA")), UmiCheck::Homopolymer);
        assert_eq!(check_umi(&enc("TTTTTT")), UmiCheck::Homopolymer);
    }

    #[test]
    fn whitelist_length_mismatch_errors() {
        let f = write_wl(&["AAAA", "TTT"]);
        assert!(CbWhitelist::load(f.path()).is_err());
    }

    #[test]
    fn whitelist_gzip_load() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let f = tempfile::Builder::new().suffix(".gz").tempfile().unwrap();
        let mut enc = GzEncoder::new(f.as_file(), Compression::default());
        writeln!(enc, "AAAA\nACGT\nTTTT").unwrap();
        enc.finish().unwrap();
        let wl = CbWhitelist::load(f.path()).unwrap();
        assert_eq!(wl.len(), 3);
    }

    #[test]
    fn match_type_parsing() {
        assert!(!CbMatchType::from_str("Exact").unwrap().mm1);
        assert!(CbMatchType::from_str("1MM").unwrap().mm1);
        assert!(!CbMatchType::from_str("1MM").unwrap().mm1_multi);
        assert!(CbMatchType::from_str("1MM_multi").unwrap().mm1_multi);
        let n = CbMatchType::from_str("1MM_multi_Nbase_pseudocounts").unwrap();
        assert!(n.mm1_multi_nbase && n.pseudocounts);
        assert!(CbMatchType::from_str("bogus").is_err());
    }
}
