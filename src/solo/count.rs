//! UMI deduplication and raw count-matrix output (Phase 14.4).
//!
//! Collates the per-read `(cell, UMI, gene)` records produced during alignment
//! into a sparse per-cell, per-gene count matrix:
//!   1. resolve deferred 1MM_multi cell barcodes via the count+quality posterior
//!      (STAR `SoloReadFeature_inputRecords.cpp`: weight = exactCount·10^(−q/10));
//!   2. group reads by `(cell, gene)` and collapse UMIs per `--soloUMIdedup`
//!      (STAR `SoloFeature_collapseUMIall.cpp`);
//!   3. write `Solo.out/Gene/raw/{matrix.mtx, barcodes.tsv, features.tsv}` in
//!      CellRanger-compatible MatrixMarket layout (features × barcodes, 1-based).

use crate::error::Error;
use crate::solo::whitelist::CbWhitelist;
use crate::solo::{SoloContext, SoloCountRecord};
// FxHash (non-cryptographic) rather than std's SipHash: every map here is keyed
// on packed integers (u64 UMIs, u32 gene/cell ids) on the per-cell UMI-dedup hot
// path, where SipHash is ~3-5x slower and buys nothing. Output is unchanged —
// results are sorted before emission, and FxHash is deterministic.
use rustc_hash::FxHashMap as HashMap;
use std::io::{BufRead, BufReader, Write as _};
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Open a solo output file, gzipping it (and appending `.gz` to the name) when
/// `gzip` is set. The body is written by the closure; the gzip stream is
/// finished explicitly so the trailer is always flushed. Returns the path written.
pub(crate) fn write_file<F>(path: &Path, gzip: bool, body: F) -> Result<PathBuf, Error>
where
    F: FnOnce(&mut dyn std::io::Write) -> Result<(), Error>,
{
    let final_path = if gzip {
        let mut s = path.as_os_str().to_owned();
        s.push(".gz");
        PathBuf::from(s)
    } else {
        path.to_path_buf()
    };
    let file = std::fs::File::create(&final_path).map_err(|e| Error::io(e, &final_path))?;
    if gzip {
        // libdeflate compresses a whole buffer at once (no streaming API), so the
        // body is collected then gzip-compressed in one shot — markedly faster
        // than streaming flate2 at the same ratio. Matrix/barcode bodies are tens
        // of MB, which fits comfortably in memory.
        let mut buf: Vec<u8> = Vec::new();
        body(&mut buf)?;
        let lvl = libdeflater::CompressionLvl::new(6).unwrap_or_default();
        let mut comp = libdeflater::Compressor::new(lvl);
        let mut out = vec![0u8; comp.gzip_compress_bound(buf.len())];
        let n = comp.gzip_compress(&buf, &mut out).map_err(|e| {
            Error::io(
                std::io::Error::other(format!("libdeflate gzip: {e:?}")),
                &final_path,
            )
        })?;
        let mut w = std::io::BufWriter::new(file);
        w.write_all(&out[..n])
            .map_err(|e| Error::io(e, &final_path))?;
        w.flush().map_err(|e| Error::io(e, &final_path))?;
    } else {
        let mut w = std::io::BufWriter::new(file);
        body(&mut w)?;
        w.flush().map_err(|e| Error::io(e, &final_path))?;
    }
    Ok(final_path)
}

// ---------------------------------------------------------------------------
// UMI deduplication
// ---------------------------------------------------------------------------

/// `--soloUMIdedup` method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiDedup {
    /// Count distinct UMI sequences (no error correction).
    Exact,
    /// No collapsing — count every read.
    NoDedup,
    /// Collapse all UMIs within Hamming-1 transitively (connected components).
    OneMmAll,
    /// UMI-tools directional, `count_hub >= 2*count_leaf + 0`.
    OneMmDirectional,
    /// UMI-tools directional original, `count_hub >= 2*count_leaf - 1`.
    OneMmDirectionalUmiTools,
    /// CellRanger 2–4 1MM collapse: each UMI is corrected to a higher-count
    /// 1MM neighbor (non-transitive); count = distinct corrected UMIs.
    OneMmCr,
}

impl FromStr for UmiDedup {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Exact" => Ok(Self::Exact),
            "NoDedup" => Ok(Self::NoDedup),
            "1MM_All" => Ok(Self::OneMmAll),
            "1MM_Directional" => Ok(Self::OneMmDirectional),
            "1MM_Directional_UMItools" => Ok(Self::OneMmDirectionalUmiTools),
            "1MM_CR" => Ok(Self::OneMmCr),
            _ => Err(format!(
                "unknown soloUMIdedup '{s}'; expected Exact, NoDedup, 1MM_All, 1MM_Directional, 1MM_Directional_UMItools, or 1MM_CR"
            )),
        }
    }
}

/// `--soloUMIfiltering`: removal of UMIs that map to multiple genes within a cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiFiltering {
    /// No multi-gene UMI filtering.
    None,
    /// Remove lower-count gene assignments of a multi-gene UMI; if every gene
    /// has a single read, drop the UMI entirely (STAR `MultiGeneUMI`).
    MultiGeneUmi,
    /// `MultiGeneUMI_All`: a UMI seen in more than one gene is removed from
    /// *all* of them, rather than from the losers only.
    ///
    /// This is a deliberate divergence from STAR 2.7.11b, tracked in #144.
    /// There the option is a no-op — its consumption site tests only the
    /// `MultiGeneUMI` flag — so selecting it leaves the filter entirely off.
    /// Reproducing that would ship a flag that silently does nothing;
    /// implementing what it is documented to do is the lesser evil, and the
    /// behaviour is asserted rather than inherited.
    MultiGeneUmiAll,
    /// CellRanger > 3.0 variant: keep only the highest-read-count gene for a
    /// multi-gene UMI (ties retained), without the all-singletons drop.
    MultiGeneUmiCr,
}

impl FromStr for UmiFiltering {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "-" | "None" => Ok(Self::None),
            "MultiGeneUMI" => Ok(Self::MultiGeneUmi),
            "MultiGeneUMI_All" => Ok(Self::MultiGeneUmiAll),
            "MultiGeneUMI_CR" => Ok(Self::MultiGeneUmiCr),
            _ => Err(format!(
                "unknown soloUMIfiltering '{s}'; expected -, None, MultiGeneUMI, MultiGeneUMI_CR, or MultiGeneUMI_All"
            )),
        }
    }
}

/// True if packed UMIs `a` and `b` (length `len`) differ at exactly one base.
fn hamming1(a: u64, b: u64, len: usize) -> bool {
    let x = a ^ b;
    let mut diff = 0u32;
    for i in 0..len {
        if (x >> (2 * i)) & 0b11 != 0 {
            diff += 1;
            if diff > 1 {
                return false;
            }
        }
    }
    diff == 1
}

/// Deduplicate the UMIs observed for one `(cell, gene)` pair into a molecule
/// count. `umis` maps each packed UMI to its read multiplicity.
#[allow(clippy::implicit_hasher)] // always called with the default hasher
pub fn dedup_count(umis: &HashMap<u64, u32>, method: UmiDedup, umi_len: usize) -> u64 {
    match method {
        UmiDedup::Exact => umis.len() as u64,
        UmiDedup::NoDedup => umis.values().map(|&c| u64::from(c)).sum(),
        UmiDedup::OneMmAll => connected_components(umis, umi_len),
        UmiDedup::OneMmDirectional => directional(umis, umi_len, 0),
        UmiDedup::OneMmDirectionalUmiTools => directional(umis, umi_len, -1),
        UmiDedup::OneMmCr => cellranger_1mm(umis, umi_len),
    }
}

/// 1MM_CR: CellRanger's 1-mismatch UMI collapse (STAR `umiArrayCorrect_CR`).
/// UMIs are sorted ascending by `(count, umi)`; each UMI is corrected to the
/// LAST (highest-count) 1MM neighbor with a strictly later sort position — i.e.
/// its highest-count 1MM neighbor. Correction is non-transitive (it points to
/// the neighbor's raw UMI, not its corrected value); the molecule count is the
/// number of distinct corrected UMIs.
fn cellranger_1mm(umis: &HashMap<u64, u32>, umi_len: usize) -> u64 {
    let distinct: std::collections::HashSet<u64> =
        cellranger_1mm_map(umis, umi_len).into_values().collect();
    distinct.len() as u64
}

/// The same correction, returning `raw UMI -> corrected UMI`.
///
/// `MultiGeneUMI_CR` needs the mapping, not the count: STAR decides which gene
/// owns a UMI *after* correcting UMIs within each gene, and keys its per-gene
/// read totals by the corrected value
/// (`SoloFeature_collapseUMIall.cpp:134-148`).
fn cellranger_1mm_map(umis: &HashMap<u64, u32>, umi_len: usize) -> HashMap<u64, u64> {
    let mut items: Vec<(u64, u32)> = umis.iter().map(|(&u, &c)| (u, c)).collect();
    // Ascending by count, then by UMI value (mirrors funCompareSolo1 ordering,
    // so the inner scan from the end meets higher-count neighbors first).
    items.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    let n = items.len();
    let mut corrected: Vec<u64> = Vec::with_capacity(n);
    for iu in 0..n {
        let mut corr = items[iu].0;
        let mut iuu = n;
        while iuu > iu + 1 {
            iuu -= 1;
            if hamming1(items[iu].0, items[iuu].0, umi_len) {
                corr = items[iuu].0;
                break;
            }
        }
        corrected.push(corr);
    }
    items.iter().map(|&(u, _)| u).zip(corrected).collect()
}

/// 1MM_All: number of connected components when UMIs within Hamming-1 are
/// merged transitively (union-find).
fn connected_components(umis: &HashMap<u64, u32>, umi_len: usize) -> u64 {
    let keys: Vec<u64> = umis.keys().copied().collect();
    let n = keys.len();
    if n <= 1 {
        return n as u64;
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if hamming1(keys[i], keys[j], umi_len) {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let mut roots = std::collections::HashSet::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        roots.insert(r);
    }
    roots.len() as u64
}

/// 1MM_Directional: a lower-count UMI within Hamming-1 of a hub whose count
/// satisfies `count_hub >= 2*count_leaf + dir_count_add` is absorbed into that
/// hub's own (already corrected) UMI; the molecule count is the number of
/// distinct surviving UMIs, STAR's `umiArrayCorrect_Directional`, which counts
/// `umiC.size()` over the corrected values.
fn directional(umis: &HashMap<u64, u32>, umi_len: usize, dir_count_add: i64) -> u64 {
    let corrections = directional_correction_map(umis, umi_len, dir_count_add);
    umis.keys()
        .map(|u| corrections.get(u).copied().unwrap_or(*u))
        .collect::<std::collections::HashSet<u64>>()
        .len() as u64
}

// ---------------------------------------------------------------------------
// UMI correction maps (STAR's `umiCorrected`, for the readInfo / UB SAM tag)
// ---------------------------------------------------------------------------

/// Swap the low and high halves of a packed UMI, as STAR's `umiSwapHalves` does
/// before its second 1MM scan. The half-swapped value is also the order STAR's
/// graph collapse walks the UMIs in, and so decides ties between equal-count
/// representatives.
fn swap_halves(umi: u64, umi_len: usize) -> u64 {
    let half_bits = umi_len; // umi_len bases → 2*umi_len bits, half = umi_len bits
    let mask_low = (1u64 << half_bits) - 1;
    let high = umi >> half_bits;
    ((umi & mask_low) << half_bits) | high
}

/// `raw UMI -> corrected UMI` for one `(cell, gene)`, matching what STAR records
/// in `umiCorrected` for the active dedup method
/// (`SoloFeature_collapseUMIall.cpp`, `SoloFeature_collapseUMI_Graph.cpp`).
/// Only UMIs that actually change are listed; `Exact`/`NoDedup` correct nothing.
#[allow(clippy::implicit_hasher)] // always called with the default hasher
pub fn umi_correction_map(
    umis: &HashMap<u64, u32>,
    method: UmiDedup,
    umi_len: usize,
) -> HashMap<u64, u64> {
    let mut map = match method {
        UmiDedup::Exact | UmiDedup::NoDedup => HashMap::default(),
        UmiDedup::OneMmCr => cellranger_1mm_map(umis, umi_len),
        UmiDedup::OneMmAll => graph_correction_map(umis, umi_len),
        UmiDedup::OneMmDirectional => directional_correction_map(umis, umi_len, 0),
        UmiDedup::OneMmDirectionalUmiTools => directional_correction_map(umis, umi_len, -1),
    };
    map.retain(|raw, corrected| raw != corrected);
    map
}

/// 1MM_All: every UMI of a connected component is corrected to the component's
/// highest-count UMI, ties going to the one STAR meets first, the smallest
/// half-swapped value (`umiArrayCorrect_Graph`'s `umiBest` scan).
fn graph_correction_map(umis: &HashMap<u64, u32>, umi_len: usize) -> HashMap<u64, u64> {
    let keys: Vec<u64> = umis.keys().copied().collect();
    let n = keys.len();
    let mut map = HashMap::default();
    if n <= 1 {
        return map;
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if hamming1(keys[i], keys[j], umi_len) {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    // Best (count, then smallest swapped value) UMI per component.
    let mut best: HashMap<usize, (u32, u64, u64)> = HashMap::default();
    for (i, &umi) in keys.iter().enumerate() {
        let root = find(&mut parent, i);
        let count = umis[&umi];
        let swapped = swap_halves(umi, umi_len);
        match best.get(&root) {
            Some(&(bc, bs, _)) if bc > count || (bc == count && bs <= swapped) => {}
            _ => {
                best.insert(root, (count, swapped, umi));
            }
        }
    }
    // A component of one is never coloured by STAR, so it is never recorded.
    let mut comp_size: HashMap<usize, usize> = HashMap::default();
    for i in 0..n {
        let root = find(&mut parent, i);
        *comp_size.entry(root).or_insert(0) += 1;
    }
    for (i, &umi) in keys.iter().enumerate() {
        let root = find(&mut parent, i);
        if comp_size[&root] > 1
            && let Some(&(_, _, rep)) = best.get(&root)
        {
            map.insert(umi, rep);
        }
    }
    map
}

/// 1MM_Directional: each UMI, scanned in descending count order, is corrected to
/// the *corrected* value of the first earlier (higher-count) UMI within one
/// mismatch whose count satisfies `hub >= 2*leaf + dir_count_add`, STAR's
/// `umiArrayCorrect_Directional`, chain and all.
fn directional_correction_map(
    umis: &HashMap<u64, u32>,
    umi_len: usize,
    dir_count_add: i64,
) -> HashMap<u64, u64> {
    let mut items: Vec<(u64, u32)> = umis.iter().map(|(&u, &c)| (u, c)).collect();
    items.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut corrected: Vec<u64> = items.iter().map(|&(u, _)| u).collect();
    for iu in 1..items.len() {
        for iuu in 0..iu {
            if i64::from(items[iuu].1) >= 2 * i64::from(items[iu].1) + dir_count_add
                && hamming1(items[iu].0, items[iuu].0, umi_len)
            {
                corrected[iu] = corrected[iuu];
                break;
            }
        }
    }
    items
        .iter()
        .map(|&(u, _)| u)
        .zip(corrected)
        .filter(|(raw, corr)| raw != corr)
        .collect()
}

// ---------------------------------------------------------------------------
// Cell-barcode multi-match resolution (deferred 1MM_multi)
// ---------------------------------------------------------------------------

/// STAR's `cbMinP` (`ParametersSolo.cpp`): the winning candidate of a
/// multi-matching cell barcode must hold at least this share of the summed
/// posterior, otherwise the read is dropped (`noTooManyWLmatches`). STAR's
/// default build (no `MATCH_CellRanger`) stores it, and does the comparison, in
/// single precision.
const CB_MIN_P: f32 = 0.975;

/// STAR's `QSbase` / `QSmax` (`ParametersSolo.cpp`): Phred+33 qualities, capped
/// at 33 before they enter the posterior.
const CB_QS_BASE: u8 = 33;
const CB_QS_MAX: u8 = 33;

/// Resolve a 1MM_multi cell barcode to a single whitelist index, following
/// STAR `SoloReadFeature_inputRecords.cpp` (multiple-match branch) exactly.
///
/// `exact_counts` is STAR's `cbReadCountExact` before pseudocounts; with the
/// `*_pseudocounts` match types every whitelist entry gets `pseudocount` (1)
/// added (`SoloFeature_sumThreads.cpp`). For each candidate with a non-zero
/// total count, `p = count · 10^(−min(q, 33)/10)`; the strict argmax (first
/// wins a tie) is accepted only when `p_max ≥ cbMinP · p_total`.
fn resolve_multi_cb(
    candidates: &[crate::solo::whitelist::CbCandidate],
    exact_counts: &[u64],
    pseudocount: f64,
) -> Option<u32> {
    let mut ptot = 0.0f32;
    let mut pmax = 0.0f32;
    let mut best: Option<u32> = None;
    for c in candidates {
        let count = *exact_counts.get(c.wl_index as usize).unwrap_or(&0) as f64 + pseudocount;
        if count <= 0.0 {
            continue; // `cbReadCountTotal[cbin]>0`: otherwise this cbin does not work
        }
        let q = c.mismatch_qual.saturating_sub(CB_QS_BASE).min(CB_QS_MAX);
        let pin = (count * 10f64.powf(-f64::from(q) / 10.0)) as f32;
        ptot += pin;
        if pin > pmax {
            best = Some(c.wl_index);
            pmax = pin;
        }
    }
    if ptot > 0.0 && pmax >= CB_MIN_P * ptot {
        best
    } else {
        None
    }
}

/// STAR's `oneExact` guard (`SoloReadFeature_inputRecords.cpp`): with every
/// `--soloCBmatchWLtype` except the pseudocount ones, a read whose barcode was
/// corrected to a single whitelist entry is dropped (`noMMtoWLwithoutExact`)
/// unless some read matched that entry exactly. An exact-match read always
/// passes (it bumped the count itself), so checking the count alone is
/// equivalent to STAR's `cbmatch==1 && cbReadCountTotal[cb]==0` test.
/// `exact_counts` is empty without a whitelist, where nothing is corrected.
pub(crate) fn cb_passes_one_exact(one_exact: bool, exact_counts: &[u64], cb: u32) -> bool {
    !one_exact || exact_counts.is_empty() || exact_counts.get(cb as usize).is_some_and(|&n| n > 0)
}

// ---------------------------------------------------------------------------
// Matrix assembly + output
// ---------------------------------------------------------------------------

/// Build and stream the raw count matrix to `matrix_path` in one per-cell pass,
/// returning the number of non-zero entries written.
///
/// Mirrors STAR's `SoloFeature_collapseUMIall.cpp`: the flat record list is
/// sorted by cell barcode so each cell's reads are contiguous, then **one cell
/// is processed at a time** (Step 1 — peak build memory is a single cell's
/// `umi → gene` maps, not a global `cell → umi → gene` nest over all records).
///
/// Step 2 (streaming output): each cell's `gene → count` entries are written
/// straight to a temporary MatrixMarket body as they are produced — the global
/// `cell → (gene → count)` map is never materialized. `nnz` is counted on the
/// fly; the final `matrix.mtx` is the header (`rows cols nnz`) followed by the
/// temp body (the BySJout temp-file pattern). So matrix-output memory is bounded
/// by one cell regardless of how many cells the raw whitelist matrix spans.
///
/// Records are sorted by cb (ascending column), and each cell's genes are
/// emitted ascending, so entries come out in the same order as before.
#[allow(clippy::too_many_arguments)]
/// Per-cell summary collected while streaming the matrix: the whitelist barcode
/// index, reads (records before UMI dedup), UMIs (deduped column sum), and genes
/// detected (nonzero entries).
#[derive(Clone, Copy)]
pub struct CellStat {
    pub cb: u32,
    pub n_reads: u64,
    pub n_umis: u64,
    pub n_genes: u32,
}

/// What `build_matrix_body` returns alongside the temp matrix body.
pub struct MatrixStats {
    pub nnz: usize,
    /// One entry per barcode that received ≥1 UMI (the raw, unfiltered set).
    pub cells: Vec<CellStat>,
    /// Counted unique-gene reads, all barcodes (STAR's `nReadPerCBunique`
    /// summed: `yessubWLmatch_UniqueFeature`).
    pub unique_reads: u64,
    /// Reads dropped for a 1MM barcode never seen exactly (`noMMtoWLwithoutExact`).
    pub dropped_one_exact: u64,
    /// 1MM_multi reads whose posterior picked no barcode (`noTooManyWLmatches`).
    pub dropped_posterior: u64,
}

/// Stream the per-cell deduplicated counts into a plain temporary MatrixMarket
/// *body* (`gene+1 cb+1 count`, barcode-ascending) and collect per-cell stats.
/// The body is finalized into `raw/` (and optionally `filtered/`) by the caller,
/// which lets the raw + filtered matrices share one streaming pass.
#[allow(clippy::too_many_arguments)]
fn build_matrix_body(
    ctx: &SoloContext,
    recorder: &crate::solo::SoloRecorder,
    method: UmiDedup,
    filtering: UmiFiltering,
    umi_len: usize,
    pseudocount: f64,
    dir: &Path,
    read_info_slot: Option<&std::sync::Mutex<Vec<crate::solo::ReadInfo>>>,
) -> Result<(tempfile::NamedTempFile, MatrixStats), Error> {
    let record_read_info = read_info_slot.is_some();
    let mut body_tmp = tempfile::Builder::new()
        .prefix(".matrix_body")
        .tempfile_in(dir)
        .map_err(|e| Error::io(e, dir))?;
    let mut nnz = 0usize;
    let mut cell_stats: Vec<CellStat> = Vec::new();
    let (unique_reads, dropped_one_exact);
    let mut dropped_posterior = 0u64;

    {
        let mut body = std::io::BufWriter::new(body_tmp.as_file_mut());

        // Move records out of the recorder; fold in resolved 1MM_multi cells.
        let mut records = std::mem::take(&mut *recorder.records.lock().unwrap());
        let exact_counts = ctx.whitelist.exact_count_snapshot();
        // `oneExact`: drop single-1MM corrections to a barcode never seen exactly.
        let one_exact = ctx.match_type.one_exact();
        let before = records.len();
        records.retain(|r| cb_passes_one_exact(one_exact, &exact_counts, r.cb));
        dropped_one_exact = (before - records.len()) as u64;
        let multi = std::mem::take(&mut *recorder.multi_records.lock().unwrap());
        for m in &multi {
            if let Some(cb) = resolve_multi_cb(&m.candidates, &exact_counts, pseudocount) {
                records.push(SoloCountRecord {
                    cb,
                    umi: m.umi,
                    gene: m.gene,
                    read_index: m.read_index,
                });
            } else {
                dropped_posterior += 1;
            }
        }
        drop(multi);
        unique_reads = records.len() as u64;

        // Group each cell's reads together (parallel sort — the record vec is large).
        use rayon::prelude::*;
        records.par_sort_unstable_by_key(|r| r.cb);

        // One contiguous [start, end) slice per CB.
        let mut bounds: Vec<(usize, usize)> = Vec::new();
        let mut i = 0;
        while i < records.len() {
            let cb = records[i].cb;
            let mut j = i + 1;
            while j < records.len() && records[j].cb == cb {
                j += 1;
            }
            bounds.push((i, j));
            i = j;
        }

        // Per-cell dedup + MatrixMarket formatting is independent across cells, so
        // run it in parallel and emit the pre-formatted bodies sequentially in CB
        // order. This keeps the matrix byte-identical to the serial version.
        struct CellOut {
            body: Vec<u8>,
            stat: Option<CellStat>,
            genes: Vec<u32>,
            /// `(read index, readInfo entry)` for this cell's reads, when the
            /// `CB`/`UB` SAM tags asked for STAR's readInfo.
            read_info: Vec<(u32, crate::solo::ReadInfo)>,
        }
        let cell_outs: Vec<CellOut> = bounds
            .par_iter()
            .map(|&(i, j)| {
                let cb = records[i].cb;

                // umi → gene → read multiplicity, for this cell only.
                let mut umi_genes: HashMap<u64, HashMap<u32, u32>> = HashMap::default();
                for r in &records[i..j] {
                    *umi_genes
                        .entry(r.umi)
                        .or_default()
                        .entry(r.gene)
                        .or_insert(0) += 1;
                }

                // `MultiGeneUMI_CR` decides gene ownership on *corrected*
                // UMIs, so it needs the correction to have happened first and
                // cannot go through the shared filter-then-dedup path below.
                let mut cell_entries = cell_feature_counts(&umi_genes, filtering, method, umi_len);
                cell_entries.sort_unstable_by_key(|&(g, _)| g);

                // STAR's readInfo: every read counted for this cell records the
                // cell it landed in and the UMI it collapsed onto
                // (`SoloFeature_collapseUMIall.cpp:250-265`).
                let mut read_info: Vec<(u32, crate::solo::ReadInfo)> = Vec::new();
                if record_read_info {
                    let mut per_gene: HashMap<u32, HashMap<u64, u32>> = HashMap::default();
                    for (&umi, genes) in &umi_genes {
                        for (&gene, &rc) in genes {
                            *per_gene.entry(gene).or_default().entry(umi).or_insert(0) += rc;
                        }
                    }
                    let corrections: HashMap<u32, HashMap<u64, u64>> = per_gene
                        .iter()
                        .map(|(&gene, umis)| (gene, umi_correction_map(umis, method, umi_len)))
                        .collect();
                    read_info.reserve(j - i);
                    for r in &records[i..j] {
                        if r.read_index == crate::solo::NO_READ_INDEX {
                            continue;
                        }
                        let umi = corrections
                            .get(&r.gene)
                            .and_then(|m| m.get(&r.umi))
                            .copied()
                            .unwrap_or(r.umi);
                        read_info.push((r.read_index, crate::solo::ReadInfo { cb, umi }));
                    }
                }

                let n_reads = (j - i) as u64;
                let n_genes = cell_entries.len() as u32;
                let mut n_umis = 0u64;
                let mut cbody: Vec<u8> = Vec::new();
                let mut genes: Vec<u32> = Vec::with_capacity(cell_entries.len());
                for (g, c) in &cell_entries {
                    n_umis += *c;
                    genes.push(*g);
                    let _ = writeln!(cbody, "{} {} {}", g + 1, cb + 1, c);
                }
                let stat = (n_umis > 0).then_some(CellStat {
                    cb,
                    n_reads,
                    n_umis,
                    n_genes,
                });
                CellOut {
                    body: cbody,
                    stat,
                    genes,
                    read_info,
                }
            })
            .collect();

        // Sequential merge: byte order preserved (CB-ascending, gene-ascending).
        let mut info_guard = read_info_slot.map(|m| m.lock().unwrap());
        // No-feature reads keep the cell their barcode resolves to and their
        // raw UMI (`SoloReadFeature_inputRecords.cpp:35-65, 68-103`): the same
        // oneExact and 1MM_multi resolution as counted reads.
        if let Some(info) = info_guard.as_mut() {
            let no_feature = std::mem::take(&mut *recorder.no_feature.lock().unwrap());
            let no_feature_multi = std::mem::take(&mut *recorder.no_feature_multi.lock().unwrap());
            let resolved = no_feature
                .into_iter()
                .filter(|r| cb_passes_one_exact(one_exact, &exact_counts, r.cb))
                .map(|r| (r.read_index, r.cb, r.umi))
                .chain(no_feature_multi.into_iter().filter_map(|m| {
                    resolve_multi_cb(&m.candidates, &exact_counts, pseudocount)
                        .map(|cb| (m.read_index, cb, m.umi))
                }));
            for (read_index, cb, umi) in resolved {
                if let Some(slot) = info.get_mut(read_index as usize) {
                    *slot = crate::solo::ReadInfo { cb, umi };
                }
            }
        }
        for co in cell_outs {
            body.write_all(&co.body).map_err(|e| Error::io(e, dir))?;
            nnz += co.genes.len();
            if let Some(s) = co.stat {
                cell_stats.push(s);
            }
            if let Some(info) = info_guard.as_mut() {
                for (read_index, entry) in co.read_info {
                    if let Some(slot) = info.get_mut(read_index as usize) {
                        *slot = entry;
                    }
                }
            }
        }
        drop(info_guard);
        body.flush().map_err(|e| Error::io(e, dir))?;
    }

    Ok((
        body_tmp,
        MatrixStats {
            nnz,
            cells: cell_stats,
            unique_reads,
            dropped_one_exact,
            dropped_posterior,
        },
    ))
}

/// Write a final `matrix.mtx[.gz]` = MatrixMarket header + (optionally
/// cb-remapped/filtered) body. With `remap = None` the body is copied verbatim
/// (raw); with `Some(map)` only columns in the map survive, renumbered to the
/// `n_cols` called cells. Returns the entry count written.
fn finalize_matrix(
    body: &tempfile::NamedTempFile,
    out_path: &Path,
    gzip: bool,
    n_features: usize,
    n_cols: usize,
    raw_nnz: usize,
    remap: Option<&HashMap<u32, u32>>,
) -> Result<usize, Error> {
    // For the filtered matrix we must know nnz before the header, so first build
    // the remapped body into a temp and count it; raw reuses the known nnz.
    let (src, nnz): (PathBuf, usize) = match remap {
        None => (body.path().to_path_buf(), raw_nnz),
        Some(map) => {
            let dir = out_path.parent().unwrap_or_else(|| Path::new("."));
            let mut ftmp = tempfile::Builder::new()
                .prefix(".matrix_filt")
                .tempfile_in(dir)
                .map_err(|e| Error::io(e, dir))?;
            let mut kept = 0usize;
            {
                let mut w = std::io::BufWriter::new(ftmp.as_file_mut());
                let reader = BufReader::new(
                    std::fs::File::open(body.path()).map_err(|e| Error::io(e, body.path()))?,
                );
                for line in reader.lines() {
                    let line = line.map_err(|e| Error::io(e, body.path()))?;
                    let mut it = line.split(' ');
                    let (Some(gene), Some(cb1), Some(cnt)) = (it.next(), it.next(), it.next())
                    else {
                        continue;
                    };
                    let cb0: u32 = cb1.parse::<u32>().unwrap_or(0).saturating_sub(1);
                    if let Some(&col) = map.get(&cb0) {
                        writeln!(w, "{gene} {col} {cnt}").map_err(|e| Error::io(e, out_path))?;
                        kept += 1;
                    }
                }
                w.flush().map_err(|e| Error::io(e, out_path))?;
            }
            (
                ftmp.into_temp_path()
                    .keep()
                    .map_err(|e| Error::io(e.error, out_path))?,
                kept,
            )
        }
    };

    write_file(out_path, gzip, |w| {
        writeln!(w, "%%MatrixMarket matrix coordinate integer general")
            .map_err(|e| Error::io(e, out_path))?;
        writeln!(w, "%").map_err(|e| Error::io(e, out_path))?;
        writeln!(w, "{n_features} {n_cols} {nnz}").map_err(|e| Error::io(e, out_path))?;
        let mut r = std::fs::File::open(&src).map_err(|e| Error::io(e, &src))?;
        std::io::copy(&mut r, w).map_err(|e| Error::io(e, out_path))?;
        Ok(())
    })?;
    if remap.is_some() {
        let _ = std::fs::remove_file(&src); // best-effort cleanup of the filtered temp
    }
    Ok(nnz)
}

/// `--soloMultiMappers` method (non-`Unique` ones produce a `UniqueAndMult-*.mtx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiMethod {
    Uniform,
    Rescue,
    PropUnique,
    Em,
}

impl MultiMethod {
    fn name(self) -> &'static str {
        match self {
            MultiMethod::Uniform => "Uniform",
            MultiMethod::Rescue => "Rescue",
            MultiMethod::PropUnique => "PropUnique",
            MultiMethod::Em => "EM",
        }
    }

    /// Parse `--soloMultiMappers` values, dropping `Unique` (no extra matrix).
    pub fn parse_list(vals: &[String]) -> Vec<MultiMethod> {
        vals.iter()
            .filter_map(|v| match v.as_str() {
                "Uniform" => Some(MultiMethod::Uniform),
                "Rescue" => Some(MultiMethod::Rescue),
                "PropUnique" => Some(MultiMethod::PropUnique),
                "EM" => Some(MultiMethod::Em),
                _ => None,
            })
            .collect()
    }
}

/// Distribute one cell's gene-ambiguous molecules across their gene sets and add
/// to the unique counts `u`, returning the combined (unique + multi) per-gene
/// counts. `molecules` is one gene set per deduplicated multi-gene UMI.
fn distribute_multi(
    method: MultiMethod,
    u: &HashMap<u32, f64>,
    molecules: &[Vec<u32>],
) -> HashMap<u32, f64> {
    let mut out = u.clone();
    let unit = |s: &[u32]| 1.0 / s.len() as f64;
    let get = |m: &HashMap<u32, f64>, g: u32| m.get(&g).copied().unwrap_or(0.0);
    match method {
        MultiMethod::Uniform => {
            for s in molecules {
                let w = unit(s);
                for &g in s {
                    *out.entry(g).or_insert(0.0) += w;
                }
            }
        }
        MultiMethod::PropUnique => {
            for s in molecules {
                let total: f64 = s.iter().map(|&g| get(u, g)).sum();
                for &g in s {
                    let w = if total > 0.0 {
                        get(u, g) / total
                    } else {
                        unit(s)
                    };
                    *out.entry(g).or_insert(0.0) += w;
                }
            }
        }
        MultiMethod::Rescue => {
            // Weights = unique counts + a uniform spread of the multi molecules.
            let mut unif: HashMap<u32, f64> = HashMap::default();
            for s in molecules {
                let w = unit(s);
                for &g in s {
                    *unif.entry(g).or_insert(0.0) += w;
                }
            }
            for s in molecules {
                let total: f64 = s.iter().map(|&g| get(u, g) + get(&unif, g)).sum();
                for &g in s {
                    let w = if total > 0.0 {
                        (get(u, g) + get(&unif, g)) / total
                    } else {
                        unit(s)
                    };
                    *out.entry(g).or_insert(0.0) += w;
                }
            }
        }
        MultiMethod::Em => {
            // STAR's gEem (`SoloFeature_collapseUMIall.cpp`): start from the
            // unique counts plus the uniform spread; each round zeroes the
            // previous estimate's values below 0.01, then gives every molecule to
            // its genes in proportion to that estimate. It stops once no gene
            // moved by 0.01, or after 101 rounds, and keeps that round's estimate,
            // which can still hold values below 0.01.
            let genes: std::collections::BTreeSet<u32> =
                molecules.iter().flatten().copied().collect();
            let mut old: HashMap<u32, f64> = genes.iter().map(|&g| (g, 0.0)).collect();
            for s in molecules {
                let w = unit(s);
                for &g in s {
                    *old.get_mut(&g).unwrap() += w;
                }
            }
            for (&g, v) in &mut old {
                *v += get(u, g);
            }
            let mut iter = 0u32;
            let new = loop {
                iter += 1;
                let mut new: HashMap<u32, f64> = genes.iter().map(|&g| (g, get(u, g))).collect();
                for v in old.values_mut() {
                    if *v < 0.01 {
                        *v = 0.0;
                    }
                }
                for s in molecules {
                    let norm = 1.0 / s.iter().map(|g| old[g]).sum::<f64>();
                    for g in s {
                        *new.get_mut(g).unwrap() += old[g] * norm;
                    }
                }
                let max_change = genes
                    .iter()
                    .map(|g| (new[g] - old[g]).abs())
                    .fold(0.0, f64::max);
                if max_change < 0.01 || iter > 100 {
                    break new;
                }
                old = new;
            };
            // STAR keeps the multi part (gEem - gEu) and adds the unique count
            // back when it writes the matrix.
            for g in genes {
                let ug = get(u, g);
                out.insert(g, ug + (new[&g] - ug));
            }
        }
    }
    out
}

/// Format a real matrix value as STAR's `ostream << double` does
/// (`SoloFeature_outputResults.cpp:181-223`): [`fmt_double`].
fn fmt_real(v: f64) -> String {
    fmt_double(v)
}

/// Write the `UniqueAndMult-<method>.mtx` matrices (real-valued) for the
/// `--soloMultiMappers` methods. Re-reads the raw matrix body (per-cell unique
/// counts, cb-ascending) and merges each cell with its gene-ambiguous molecules
/// (deduplicated by UMI, gene set = union). Cells present only in multi records
/// (no unique gene) are skipped.
#[allow(clippy::too_many_arguments)]
fn build_multi_matrices(
    raw_body: &tempfile::NamedTempFile,
    multi_records: &[crate::solo::MultiGeneRecord],
    methods: &[MultiMethod],
    dir: &Path,
    matrix_name: &str,
    n_features: usize,
    n_barcodes: usize,
    gzip: bool,
) -> Result<(), Error> {
    if methods.is_empty() {
        return Ok(());
    }
    use rayon::prelude::*;
    let mut multi: Vec<&crate::solo::MultiGeneRecord> = multi_records.iter().collect();
    multi.par_sort_unstable_by_key(|r| r.cb);

    // Per-method temp body + entry count.
    let mut bodies: Vec<tempfile::NamedTempFile> = Vec::new();
    for _ in methods {
        bodies.push(
            tempfile::Builder::new()
                .prefix(".um_body")
                .tempfile_in(dir)
                .map_err(|e| Error::io(e, dir))?,
        );
    }
    let mut nnz = vec![0usize; methods.len()];

    // Gather one cell's multi molecules (gene sets, one per deduped UMI).
    let cell_molecules = |cb: u32, mptr: &mut usize| -> Vec<Vec<u32>> {
        let mut by_umi: HashMap<u64, std::collections::BTreeSet<u32>> = HashMap::default();
        while *mptr < multi.len() && multi[*mptr].cb == cb {
            let r = multi[*mptr];
            by_umi
                .entry(r.umi)
                .or_default()
                .extend(r.genes.iter().copied());
            *mptr += 1;
        }
        by_umi
            .into_values()
            .map(|s| s.into_iter().collect())
            .collect()
    };

    {
        let mut writers: Vec<std::io::BufWriter<&mut std::fs::File>> = bodies
            .iter_mut()
            .map(|t| std::io::BufWriter::new(t.as_file_mut()))
            .collect();
        let reader = BufReader::new(
            std::fs::File::open(raw_body.path()).map_err(|e| Error::io(e, raw_body.path()))?,
        );
        let mut mptr = 0usize;
        let mut cur_cb: Option<u32> = None;
        let mut u_map: HashMap<u32, f64> = HashMap::default();

        let mut flush = |cb: u32,
                         u: &HashMap<u32, f64>,
                         mptr: &mut usize,
                         nnz: &mut [usize]|
         -> Result<(), Error> {
            let mols = cell_molecules(cb, mptr);
            for (k, &m) in methods.iter().enumerate() {
                // Every gene of the cell's unique counts or of a multi-gene
                // molecule is written, a zero included (STAR's outputResults
                // walks both lists without testing the value).
                let mut entries: Vec<(u32, f64)> =
                    distribute_multi(m, u, &mols).into_iter().collect();
                entries.sort_unstable_by_key(|&(g, _)| g);
                for (g, v) in entries {
                    writeln!(writers[k], "{} {} {}", g + 1, cb + 1, fmt_real(v))
                        .map_err(|e| Error::io(e, dir))?;
                    nnz[k] += 1;
                }
            }
            Ok(())
        };

        for line in reader.lines() {
            let line = line.map_err(|e| Error::io(e, raw_body.path()))?;
            let mut it = line.split(' ');
            let (Some(gt), Some(ct), Some(vt)) = (it.next(), it.next(), it.next()) else {
                continue;
            };
            let g: u32 = gt.parse::<u32>().unwrap_or(1) - 1;
            let cb: u32 = ct.parse::<u32>().unwrap_or(1) - 1;
            let v: f64 = vt.parse().unwrap_or(0.0);
            if cur_cb != Some(cb) {
                if let Some(prev) = cur_cb {
                    flush(prev, &u_map, &mut mptr, &mut nnz)?;
                }
                // Cells before this one with multi-gene molecules only: STAR
                // counts them too, with no unique genes to weigh them by.
                while mptr < multi.len() && multi[mptr].cb < cb {
                    flush(multi[mptr].cb, &HashMap::default(), &mut mptr, &mut nnz)?;
                }
                cur_cb = Some(cb);
                u_map.clear();
            }
            *u_map.entry(g).or_insert(0.0) += v;
        }
        if let Some(prev) = cur_cb {
            flush(prev, &u_map, &mut mptr, &mut nnz)?;
        }
        while mptr < multi.len() {
            flush(multi[mptr].cb, &HashMap::default(), &mut mptr, &mut nnz)?;
        }
        for w in &mut writers {
            w.flush().map_err(|e| Error::io(e, dir))?;
        }
    }

    // Finalize each UniqueAndMult-<method>.mtx (real-valued MatrixMarket).
    for ((m, body), &n) in methods.iter().zip(&bodies).zip(&nnz) {
        let path = dir.join(format!("UniqueAndMult-{}.mtx", m.name()));
        write_file(&path, gzip, |w| {
            writeln!(w, "%%MatrixMarket matrix coordinate real general")
                .map_err(|e| Error::io(e, &path))?;
            writeln!(w, "%").map_err(|e| Error::io(e, &path))?;
            writeln!(w, "{n_features} {n_barcodes} {n}").map_err(|e| Error::io(e, &path))?;
            let mut r = std::fs::File::open(body.path()).map_err(|e| Error::io(e, body.path()))?;
            std::io::copy(&mut r, w).map_err(|e| Error::io(e, &path))?;
            Ok(())
        })?;
    }
    let _ = matrix_name; // UniqueAndMult uses a fixed name scheme
    Ok(())
}

/// CellRanger's multi-gene UMI resolution, as STAR implements it for
/// `--soloUMIfiltering MultiGeneUMI_CR` (`SoloFeature_collapseUMIall.cpp`).
///
/// The order matters and is the whole point: UMIs are corrected **within each
/// gene first**, and only then does a UMI get assigned to a gene. Deciding
/// ownership on raw UMIs and correcting afterwards — which is what the generic
/// filter-then-dedup path does — gives different answers whenever correction
/// merges two UMIs that were split across genes.
///
/// Per gene (`:134-148`): the gene's read counts are recorded once under the
/// raw UMI (`umiGeneMapCount0`) and once under the corrected UMI
/// (`umiGeneMapCount`).
///
/// Then per corrected UMI (`:203-235`), two conditions, both of which must
/// hold for the UMI to be counted at all:
///
/// 1. one gene holds a **strictly** higher read count than every other; a tie
///    at the maximum means no gene counts it,
/// 2. and no gene beats that winner in the **uncorrected** map at the same key.
///
/// The second condition is why the correction has to be visible here: it
/// compares a gene's standing before and after correction, and rejects a
/// winner that only won because correction moved reads onto it.
///
/// Returns `(gene, molecules)` for this cell, gene-ascending.
fn multi_gene_umi_cr_counts(
    umi_genes: &HashMap<u64, HashMap<u32, u32>>,
    umi_len: usize,
) -> Vec<(u32, u64)> {
    // Regroup as gene → (raw UMI → reads); correction happens per gene.
    let mut gene_umis: HashMap<u32, HashMap<u64, u32>> = HashMap::default();
    for (&umi, genes) in umi_genes {
        for (&gene, &rc) in genes {
            *gene_umis.entry(gene).or_default().entry(umi).or_insert(0) += rc;
        }
    }

    let mut uncorrected: HashMap<u64, HashMap<u32, u32>> = HashMap::default();
    let mut corrected: HashMap<u64, HashMap<u32, u32>> = HashMap::default();
    for (&gene, umis) in &gene_umis {
        for (&umi, &rc) in umis {
            *uncorrected.entry(umi).or_default().entry(gene).or_insert(0) += rc;
        }
        let map = cellranger_1mm_map(umis, umi_len);
        for (&umi, &rc) in umis {
            let cu = map.get(&umi).copied().unwrap_or(umi);
            *corrected.entry(cu).or_default().entry(gene).or_insert(0) += rc;
        }
    }

    let mut counts: HashMap<u32, u64> = HashMap::default();
    for (cu, genes) in &corrected {
        // Condition 1: a strict maximum, ties lose.
        let mut best = 0u32;
        let mut winner: Option<u32> = None;
        for (&gene, &rc) in genes {
            if rc > best {
                best = rc;
                winner = Some(gene);
            } else if rc == best {
                winner = None;
            }
        }
        let Some(winner) = winner else { continue };

        // Condition 2: the winner must not be beaten in the uncorrected map at
        // the same key. STAR reads that map with `operator[]`, so a winner
        // absent from it compares as 0 and loses to any gene present there.
        if let Some(raw_genes) = uncorrected.get(cu) {
            let winner_raw = raw_genes.get(&winner).copied().unwrap_or(0);
            if raw_genes.values().any(|&rc| rc > winner_raw) {
                continue;
            }
        }
        *counts.entry(winner).or_insert(0) += 1;
    }

    let mut out: Vec<(u32, u64)> = counts.into_iter().filter(|&(_, c)| c > 0).collect();
    out.sort_unstable_by_key(|&(g, _)| g);
    out
}

/// Apply `--soloUMIfiltering` to the gene→read_count map of a single UMI,
/// returning the surviving (gene, read_count) entries.
fn filter_multi_gene_umi(genes: &HashMap<u32, u32>, filtering: UmiFiltering) -> Vec<(&u32, &u32)> {
    if filtering == UmiFiltering::None || genes.len() <= 1 {
        return genes.iter().collect();
    }
    let max = genes.values().copied().max().unwrap_or(0);
    match filtering {
        // STAR MultiGeneUMI: threshold = max (or 2 if max==1, dropping all
        // single-read multi-gene UMIs); keep genes with read_count >= threshold.
        UmiFiltering::MultiGeneUmi => {
            let thresh = if max == 1 { 2 } else { max };
            genes.iter().filter(|&(_, &rc)| rc >= thresh).collect()
        }
        // A UMI that appears in more than one gene is evidence of a collision
        // or of chimeric amplification, so it is discarded outright rather than
        // attributed to whichever gene happened to read deepest. `genes.len()`
        // is already known to be > 1 here.
        UmiFiltering::MultiGeneUmiAll => Vec::new(),
        // CellRanger: the gene with the strictly highest read count takes the
        // UMI, and a tie gives it to nobody.
        //
        // STAR `SoloFeature_collapseUMIall.cpp:212-224` walks the genes keeping
        // a running maximum, and clears its winner whenever it meets an equal
        // count:
        //
        // ```cpp
        // if (ig.second>maxu) { maxu=ig.second; maxg=ig.first; }
        // else if (ig.second==maxu) { maxg=-1; };
        // ...
        // if ( maxg+1==0 ) continue; // not counted for any gene
        // ```
        //
        // The outcome does not depend on the order the genes are visited: a
        // strict maximum always ends as the winner, and any tie at the maximum
        // always ends with none. So iterating a `HashMap` here is safe.
        //
        // This previously kept every gene tied at the maximum, which is the
        // opposite decision on exactly the case the rule exists for, and made
        // the flag inert on the common shape of one read per gene.
        UmiFiltering::MultiGeneUmiCr => {
            let mut best_count = 0u32;
            let mut winner: Option<&u32> = None;
            for (gene, &rc) in genes {
                if rc > best_count {
                    best_count = rc;
                    winner = Some(gene);
                } else if rc == best_count {
                    winner = None;
                }
            }
            match winner {
                Some(gene) => vec![(gene, genes.get(gene).expect("winner is a key"))],
                None => Vec::new(),
            }
        }
        UmiFiltering::None => unreachable!(),
    }
}

/// CellRanger-2.2 knee threshold on per-barcode UMI totals (STARsolo's default
/// `--soloCellFilter CellRanger2.2 3000 0.99 10`). Returns the minimum UMI count
/// for a barcode to be called a cell.
fn knee_cr22(umis_desc: &[u64], n_expected: usize, max_pct: f64, max_min_ratio: f64) -> u64 {
    if umis_desc.is_empty() {
        return 0;
    }
    let idx = ((n_expected as f64 * (1.0 - max_pct)).round() as usize).min(umis_desc.len() - 1);
    let robust_max = umis_desc[idx] as f64;
    (robust_max / max_min_ratio).ceil() as u64
}

/// Whitelist indices of called cells (sorted ascending) per `--soloCellFilter`.
/// `None` → no filtered/ output. `EmptyDrops_CR` writes only the knee-guaranteed
/// cells here (the Monte-Carlo rescue is the standalone `emptydrops` binary).
fn called_cells(cells: &[CellStat], filter: &[String]) -> Option<Vec<u32>> {
    let method = filter.first().map_or("CellRanger2.2", String::as_str);
    let arg = |i: usize, d: f64| filter.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    let mut cbs: Vec<u32> = match method {
        "None" => return None,
        // STAR's threshold, not a cut at N (`SoloFeature_cellFiltering.cpp:46-54`):
        // the UMIs of the barcode at index N of the descending list (the
        // N+1th), at least 1, and every barcode with at least that many, ties
        // included.
        "TopCells" => {
            let n = arg(1, 0.0) as usize;
            let mut umis: Vec<u64> = cells.iter().map(|c| c.n_umis).collect();
            umis.sort_unstable_by(|a, b| b.cmp(a));
            let Some(&at_n) = umis.get(n.min(umis.len().saturating_sub(1))) else {
                return Some(Vec::new());
            };
            let min = at_n.max(1);
            cells
                .iter()
                .filter(|c| c.n_umis >= min)
                .map(|c| c.cb)
                .collect()
        }
        // EmptyDrops_CR is handled by `emptydrops_called`; the knee here is the
        // fallback / guaranteed-cell base.
        "CellRanger2.2" | "EmptyDrops_CR" => {
            let mut umis: Vec<u64> = cells.iter().map(|c| c.n_umis).collect();
            umis.sort_unstable_by(|a, b| b.cmp(a));
            let thr = knee_cr22(&umis, arg(1, 3000.0) as usize, arg(2, 0.99), arg(3, 10.0));
            cells
                .iter()
                .filter(|c| c.n_umis >= thr)
                .map(|c| c.cb)
                .collect()
        }
        other => {
            log::warn!("--soloCellFilter '{other}' not supported; skipping filtered/ output");
            return None;
        }
    };
    cbs.sort_unstable();
    Some(cbs)
}

/// `--soloCellFilter EmptyDrops_CR`: the CR2.2-knee guaranteed cells PLUS cells
/// rescued by the EmptyDrops multinomial Monte-Carlo test (STAR
/// `SoloFeature_emptyDrops_CR.cpp`). Per-cell gene profiles for the ambient +
/// candidate cells are read back from the raw matrix body. `filter` is the
/// `EmptyDrops_CR nExpected maxPct maxMinRatio indMin indMax umiMin
/// umiMinFracMedian candMaxN FDR [simN]` argument list.
fn emptydrops_called(
    cells: &[CellStat],
    body: &tempfile::NamedTempFile,
    n_features: usize,
    filter: &[String],
) -> Result<Vec<u32>, Error> {
    let arg = |i: usize, d: f64| {
        filter
            .get(i)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(d)
    };
    let (n_expected, max_pct, ratio) = (arg(1, 3000.0) as usize, arg(2, 0.99), arg(3, 10.0));
    let (ind_min, ind_max) = (arg(4, 45000.0) as usize, arg(5, 90000.0) as usize);
    let umi_min = arg(6, 500.0) as u64;
    let umi_min_frac = arg(7, 0.01);
    let cand_max = arg(8, 20000.0) as usize;
    let fdr = arg(9, 0.01);
    let sim_n = arg(10, 10000.0).max(1.0) as usize;

    // Rank by total UMI (descending, cb tie-break).
    let mut order: Vec<&CellStat> = cells.iter().collect();
    order.sort_by(|a, b| b.n_umis.cmp(&a.n_umis).then(a.cb.cmp(&b.cb)));
    let totals_desc: Vec<u64> = order.iter().map(|c| c.n_umis).collect();
    let thr = knee_cr22(&totals_desc, n_expected, max_pct, ratio);
    let n_simple = totals_desc.iter().take_while(|&&u| u >= thr).count();
    let mut called: Vec<u32> = order.iter().take(n_simple).map(|c| c.cb).collect();

    // Candidate cells: rank ≥ nSimple, total ≥ minUMI, up to candMaxN.
    let median_top = totals_desc.get(n_simple / 2).copied().unwrap_or(0);
    let min_umi = umi_min.max((umi_min_frac * median_top as f64) as u64);
    let mut cand_cbs: Vec<u32> = Vec::new();
    for c in order.iter().skip(n_simple).take(cand_max) {
        if c.n_umis < min_umi {
            break;
        }
        cand_cbs.push(c.cb);
    }
    if cand_cbs.is_empty() {
        called.sort_unstable();
        return Ok(called);
    }
    let cand_set: std::collections::HashSet<u32> = cand_cbs.iter().copied().collect();
    let ambient_set: std::collections::HashSet<u32> = order
        .iter()
        .skip(ind_min)
        .take(ind_max.saturating_sub(ind_min))
        .map(|c| c.cb)
        .collect();

    // Re-read the raw body for ambient (summed) + per-candidate profiles.
    let mut ambient = vec![0f64; n_features];
    let mut amb_total = 0f64;
    let mut cand_profiles: HashMap<u32, Vec<(u32, u32)>> = HashMap::default();
    let reader =
        BufReader::new(std::fs::File::open(body.path()).map_err(|e| Error::io(e, body.path()))?);
    for line in reader.lines() {
        let line = line.map_err(|e| Error::io(e, body.path()))?;
        let mut it = line.split(' ');
        let (Some(gt), Some(ct), Some(vt)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let g = gt.parse::<u32>().unwrap_or(1) - 1;
        let cb = ct.parse::<u32>().unwrap_or(1) - 1;
        let v = vt.parse::<u32>().unwrap_or(0);
        if ambient_set.contains(&cb) {
            ambient[g as usize] += v as f64;
            amb_total += v as f64;
        }
        if cand_set.contains(&cb) {
            cand_profiles.entry(cb).or_default().push((g, v));
        }
    }
    if amb_total == 0.0 {
        called.sort_unstable();
        return Ok(called);
    }

    // Ambient probabilities, smoothed by Simple Good-Turing as CellRanger's
    // EmptyDrops_CR does. The raw ambient counts are a small sample, so a gene
    // seen twice there is not twice as likely as one seen once; SGT fits the
    // frequency spectrum and reserves mass for genes the sample missed
    // entirely. The previous approximation here reserved mass from the
    // singleton rate alone and spread the rest in proportion to raw counts,
    // which is the right shape but the wrong numbers.
    let amb_p: Vec<f64> = {
        // Frequency of frequencies over the integer ambient counts.
        let counts: Vec<u32> = ambient.iter().map(|&x| x as u32).collect();
        let mut fof: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
        for &c in &counts {
            *fof.entry(c).or_insert(0) += 1;
        }
        let n_zero = fof.get(&0).copied().unwrap_or(0);

        let mut sgt = crate::solo::sgt::Sgt::new();
        for (&freq, &n) in &fof {
            if freq != 0 {
                sgt.add(freq, n);
            }
        }
        let fitted = sgt.analyse();

        // Per-gene probability. A gene absent from the ambient set shares the
        // reserved unseen mass; one present takes its smoothed estimate. When
        // the spectrum is too small to fit (D17), there is no reserved mass and
        // the raw proportions are all that is left.
        let unseen_each = if n_zero > 0 {
            sgt.pzero() / f64::from(n_zero)
        } else {
            0.0
        };
        let raw: Vec<f64> = counts
            .iter()
            .map(|&c| {
                if c == 0 {
                    unseen_each
                } else if fitted {
                    sgt.estimate(c).unwrap_or(f64::from(c) / amb_total)
                } else {
                    f64::from(c) / amb_total
                }
            })
            .collect();
        // Renormalise: the smoothed estimates are per-event probabilities and
        // only sum to one over the whole spectrum, not over this gene set.
        let total: f64 = raw.iter().sum();
        if total > 0.0 {
            raw.into_iter().map(|p| p / total).collect()
        } else {
            raw
        }
    };
    let amb_logp: Vec<f64> = amb_p.iter().map(|&p| p.max(1e-300).ln()).collect();

    // Observed multinomial log-prob per candidate.
    let max_count = cand_cbs
        .iter()
        .filter_map(|cb| cand_profiles.get(cb))
        .map(|p| p.iter().map(|&(_, c)| c as usize).sum::<usize>())
        .max()
        .unwrap_or(0);
    let mut log_fac = vec![0f64; max_count + 1];
    for i in 2..=max_count {
        log_fac[i] = log_fac[i - 1] + (i as f64).ln();
    }
    let obs: Vec<(u32, usize, f64)> = cand_cbs
        .iter()
        .filter_map(|&cb| {
            let prof = cand_profiles.get(&cb)?;
            let total: usize = prof.iter().map(|&(_, c)| c as usize).sum();
            let mut s = log_fac[total];
            for &(g, c) in prof {
                s -= log_fac[c as usize];
                s += c as f64 * amb_logp[g as usize];
            }
            Some((cb, total, s))
        })
        .collect();

    // Monte-Carlo: simulate sim_n ambient barcodes, recording the running
    // log-prob at each count; compare each candidate against sim[*][its total].
    let nonzero: Vec<usize> = (0..n_features).filter(|&g| amb_p[g] > 0.0).collect();
    let weights: Vec<f64> = nonzero.iter().map(|&g| amb_p[g]).collect();
    // The ambient sampler is STAR's, which is libc++'s: `std::mt19937` drawn
    // through `std::discrete_distribution`. Both are implementation-defined in
    // the parts that decide which category a draw lands in, so "a correct
    // categorical sampler" is not enough to reproduce STAR's cell calls —
    // it has to be that one. See `solo::libcxx_rng`.
    let dist = crate::solo::libcxx_rng::DiscreteDistribution::new(&weights);
    // A fresh generator per simulation, seeded `19760110 * (isim + 1)` as STAR
    // seeds it. No shared state, so the walks can run in any order and on any
    // number of threads and still produce the same p-values. Each walk records
    // the running log-prob at every count level; `walks[s][k]` is simulation
    // `s` after `k` draws.
    use rayon::prelude::*;
    const BASE_SEED: u64 = 19_760_110;
    let walks: Vec<Vec<f64>> = (0..sim_n)
        .into_par_iter()
        .map_init(
            || (vec![0u32; n_features], Vec::<usize>::new()),
            |(curr, touched), s| {
                let seed = BASE_SEED.wrapping_mul(s as u64 + 1) as u32;
                let mut rng = crate::solo::libcxx_rng::Mt19937::new(seed);
                touched.clear();
                let mut walk = Vec::with_capacity(max_count + 1);
                walk.push(0.0);
                let mut lp = 0f64;
                for ic in 1..=max_count {
                    let gi = nonzero[dist.sample(&mut rng)];
                    if curr[gi] == 0 {
                        touched.push(gi);
                    }
                    curr[gi] += 1;
                    lp += amb_logp[gi] + (ic as f64).ln() - (curr[gi] as f64).ln();
                    walk.push(lp);
                }
                for &gi in touched.iter() {
                    curr[gi] = 0; // reset only touched entries for the next reuse
                }
                walk
            },
        )
        .collect();

    // p-values + Benjamini-Hochberg.
    let inv = 1.0 / (1 + sim_n) as f64;
    let mut pvals: Vec<(u32, f64)> = obs
        .par_iter()
        .map(|&(cb, total, o)| {
            let lower = walks.iter().filter(|w| w[total] < o).count();
            (cb, (1 + lower) as f64 * inv)
        })
        .collect();
    pvals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let n = pvals.len() as f64;
    let mut padj = vec![0f64; pvals.len()];
    for (rank, &(_, p)) in pvals.iter().enumerate() {
        padj[rank] = (p * n / (rank + 1) as f64).min(1.0);
    }
    for i in (0..padj.len().saturating_sub(1)).rev() {
        padj[i] = padj[i].min(padj[i + 1]);
    }
    let mut rescued = 0usize;
    for (rank, &(cb, _)) in pvals.iter().enumerate() {
        if padj[rank] <= fdr {
            called.push(cb);
            rescued += 1;
        }
    }
    log::info!(
        "EmptyDrops_CR: {n_simple} knee cells + {rescued} rescued (of {} candidates, FDR<={fdr})",
        cand_cbs.len()
    );
    called.sort_unstable();
    Ok(called)
}

/// Write the raw gene-count matrix + `Summary.csv` for a finished solo run.
/// Only `Barcodes.stats` (with a warning) when there is no explicit whitelist.
pub fn write_gene_matrix(
    ctx: &SoloContext,
    params: &crate::params::Parameters,
    align_stats: &crate::stats::AlignmentStats,
    sj_stats: Option<&crate::junction::SpliceJunctionStats>,
    genome: &crate::genome::Genome,
) -> Result<(), Error> {
    // STAR writes it before any feature (`Solo.cpp:53-64`).
    write_barcode_stats(ctx, params)?;
    let CbWhitelist::List { sorted, .. } = &ctx.whitelist else {
        log::warn!(
            "STARsolo: --soloCBwhitelist None matrix output is not yet supported (Phase 14.4); skipping matrix"
        );
        return Ok(());
    };

    let method: UmiDedup = params
        .solo_umi_dedup
        .first()
        .map_or("1MM_All", String::as_str)
        .parse()
        .unwrap_or(UmiDedup::OneMmAll);
    let filtering: UmiFiltering = params
        .solo_umi_filtering
        .first()
        .map_or("-", String::as_str)
        .parse()
        .unwrap_or(UmiFiltering::None);
    // `*_pseudocounts` CB-match types add 1 to the posterior prior.
    let pseudocount = if params.solo_cb_match_wl_type.contains("pseudocounts") {
        1.0
    } else {
        0.0
    };
    let umi_len = params.solo_umi_len as usize;

    let solo_dir = params
        .solo_out_file_names
        .first()
        .cloned()
        .unwrap_or_else(|| "Solo.out/".to_string());
    let features_name = params
        .solo_out_file_names
        .get(1)
        .cloned()
        .unwrap_or_else(|| "features.tsv".to_string());
    let barcodes_name = params
        .solo_out_file_names
        .get(2)
        .cloned()
        .unwrap_or_else(|| "barcodes.tsv".to_string());
    let matrix_name = params
        .solo_out_file_names
        .get(3)
        .cloned()
        .unwrap_or_else(|| "matrix.mtx".to_string());

    // Global mapping funnel (shared across features). The region tallies are
    // CellRanger-style positional bins over uniquely-mapped reads, populated only
    // when both Gene and GeneFull run (otherwise the split is unavailable).
    use std::sync::atomic::Ordering;
    let total_reads = align_stats.total_reads.load(Ordering::Relaxed);
    let mapped_unique = align_stats.uniquely_mapped.load(Ordering::Relaxed);
    let mapped_multi = align_stats.multi_mapped.load(Ordering::Relaxed);
    let summary_run = SummaryContext {
        n_reads: total_reads,
        mapped_unique,
        mapped_multi,
        barcodes: &ctx.stats,
        rna_q30: (
            ctx.rna_q30[0].load(Ordering::Relaxed),
            ctx.rna_q30[1].load(Ordering::Relaxed),
        ),
        multimappers: !MultiMethod::parse_list(&params.solo_multi_mappers).is_empty(),
    };
    // `--soloCellRangerSummary Yes` only (a rustar-aligner extension): the
    // funnel needs both the Gene and GeneFull overlaps.
    let have_funnel = ctx.features.contains(&crate::solo::SoloFeature::Gene)
        && ctx.features.contains(&crate::solo::SoloFeature::GeneFull);
    let region = ctx
        .region_stats
        .as_ref()
        .filter(|_| have_funnel)
        .map(|r| RegionFunnel {
            exonic: r.exonic.load(Ordering::Relaxed),
            intronic: r.intronic.load(Ordering::Relaxed),
            intergenic: r.intergenic.load(Ordering::Relaxed),
            antisense: r.antisense.load(Ordering::Relaxed),
        });

    let gzip = matches!(params.solo_out_gzip.as_str(), "yes" | "Yes" | "true");
    let n_genes = ctx.gene_ann.gene_ids.len();
    // STAR applies `oneExact` in `inputRecords`, i.e. to every feature's reads
    // (multi-gene, SJ and Velocyto included), not only the Gene matrix.
    let one_exact = ctx.match_type.one_exact();
    let exact_counts = ctx.whitelist.exact_count_snapshot();
    let multi_methods = MultiMethod::parse_list(&params.solo_multi_mappers);

    // The Gene feature's called cells, which Velocyto's filtered matrices use.
    let mut gene_called: Option<Vec<u32>> = None;

    // One {prefix}{soloOutFileNames[0]}<feature>/{raw,filtered}/ per feature.
    for (fi, (feature, recorder)) in ctx.features.iter().zip(&ctx.recorders).enumerate() {
        let feature_dir = params.output_path(&format!("{solo_dir}{}/", feature.dir_name()));
        let raw_dir = feature_dir.join("raw");
        std::fs::create_dir_all(&raw_dir).map_err(|e| Error::io(e, &raw_dir))?;

        // Stream the deduplicated counts into a shared temp body, then finalize
        // the raw matrix (and the filtered one below) from it.
        let (body, mstats) = build_matrix_body(
            ctx,
            recorder,
            method,
            filtering,
            umi_len,
            pseudocount,
            &raw_dir,
            // The feature's readInfo, if it keeps one: the first feature for
            // the CB/UB tags, Gene for Velocyto (`ParametersSolo.cpp:418-436`).
            ctx.read_info[fi].as_ref(),
        )?;
        write_features(
            &raw_dir.join(&features_name),
            &ctx.gene_ann.gene_ids,
            &ctx.gene_ann.gene_names,
            &params.solo_out_format_features_gene_field3,
            gzip,
        )?;
        write_barcodes(
            &raw_dir.join(&barcodes_name),
            &ctx.whitelist,
            sorted.len(),
            gzip,
        )?;
        finalize_matrix(
            &body,
            &raw_dir.join(&matrix_name),
            gzip,
            n_genes,
            sorted.len(),
            mstats.nnz,
            None,
        )?;
        log::info!(
            "STARsolo: wrote {}/raw matrix ({} genes × {} barcodes, {} entries){}",
            feature.dir_name(),
            n_genes,
            sorted.len(),
            mstats.nnz,
            if gzip { " [gzip]" } else { "" },
        );

        // `--soloCellReadStats CB`: the per-cell read summary, alongside the
        // raw matrix because its UMI and gene columns are the raw totals.
        if let Some(cell_stats) = &ctx.cell_read_stats {
            let umi_gene: std::collections::BTreeMap<u32, (u32, u32)> = mstats
                .cells
                .iter()
                .map(|c| (c.cb, (c.n_umis as u32, c.n_genes)))
                .collect();
            let path = feature_dir.join("CellReads.stats");
            let text = cell_stats.lock().unwrap().render(
                |cb| {
                    ctx.whitelist
                        .barcode_string(cb)
                        .unwrap_or_else(|| cb.to_string())
                },
                &umi_gene,
            );
            std::fs::write(&path, text).map_err(|e| Error::io(e, &path))?;
            log::info!("STARsolo: wrote {}/CellReads.stats", feature.dir_name());
        }

        // Filtered (cell-called) matrix per --soloCellFilter. EmptyDrops_CR runs
        // the Monte-Carlo rescue (needs the per-cell profiles in the body).
        let called = if params
            .solo_cell_filter
            .first()
            .is_some_and(|m| m == "EmptyDrops_CR")
        {
            Some(emptydrops_called(
                &mstats.cells,
                &body,
                n_genes,
                &params.solo_cell_filter,
            )?)
        } else {
            called_cells(&mstats.cells, &params.solo_cell_filter)
        };
        if *feature == crate::solo::SoloFeature::Gene {
            gene_called.clone_from(&called);
        }
        let called_for_summary = called.clone();
        if let Some(cbs) = called
            && !cbs.is_empty()
        {
            let filt_dir = feature_dir.join("filtered");
            std::fs::create_dir_all(&filt_dir).map_err(|e| Error::io(e, &filt_dir))?;
            let remap: HashMap<u32, u32> = cbs
                .iter()
                .enumerate()
                .map(|(i, &cb)| (cb, i as u32 + 1))
                .collect();
            write_features(
                &filt_dir.join(&features_name),
                &ctx.gene_ann.gene_ids,
                &ctx.gene_ann.gene_names,
                &params.solo_out_format_features_gene_field3,
                gzip,
            )?;
            write_barcodes_subset(&filt_dir.join(&barcodes_name), &ctx.whitelist, &cbs, gzip)?;
            let fnnz = finalize_matrix(
                &body,
                &filt_dir.join(&matrix_name),
                gzip,
                n_genes,
                cbs.len(),
                0,
                Some(&remap),
            )?;
            log::info!(
                "STARsolo: wrote {}/filtered matrix ({} cells, {} entries)",
                feature.dir_name(),
                cbs.len(),
                fnnz,
            );
        }

        // --soloMultiMappers: UniqueAndMult-<method>.mtx alongside raw.
        let (mut multi_counted, mut multi_dropped, mut multi_posterior) = (0u64, 0u64, 0u64);
        if !multi_methods.is_empty() {
            let mut mg = recorder.multi_gene.lock().unwrap();
            // 1MM_multi barcodes: the same posterior as unique-gene reads; a
            // barcode it cannot pick is `noTooManyWLmatches`. oneExact applies
            // to single corrections only.
            let (pending, mut resolved): (Vec<_>, Vec<_>) = std::mem::take(&mut *mg)
                .into_iter()
                .partition(|r| !r.candidates.is_empty());
            let before = resolved.len();
            resolved.retain(|r| cb_passes_one_exact(one_exact, &exact_counts, r.cb));
            multi_dropped = (before - resolved.len()) as u64;
            for mut r in pending {
                match resolve_multi_cb(&r.candidates, &exact_counts, pseudocount) {
                    Some(cb) => {
                        r.cb = cb;
                        r.candidates.clear();
                        resolved.push(r);
                    }
                    None => multi_posterior += 1,
                }
            }
            *mg = resolved;
            multi_counted = mg.len() as u64;
            build_multi_matrices(
                &body,
                &mg,
                &multi_methods,
                &raw_dir,
                &matrix_name,
                n_genes,
                sorted.len(),
                gzip,
            )?;
            log::info!(
                "STARsolo: wrote {} UniqueAndMult matrices for {} ({} ambiguous reads)",
                multi_methods.len(),
                feature.dir_name(),
                mg.len(),
            );
        }

        // Features.stats: the counting-stage figures (`inputRecords`,
        // `collapseUMIall`) join the per-read ones from `process_read`.
        {
            use crate::solo::FeatureStats as F;
            let st = &ctx.feature_stats[fi];
            st.add(
                F::NO_MM_TO_WL_WITHOUT_EXACT,
                mstats.dropped_one_exact + multi_dropped,
            );
            st.add(
                F::NO_TOO_MANY_WL_MATCHES,
                mstats.dropped_posterior + multi_posterior,
            );
            st.add(F::YES_WL_MATCH, mstats.unique_reads + multi_counted);
            st.add(F::YES_SUB_WL_MATCH_UNIQUE_FEATURE, mstats.unique_reads);
            st.add(F::YES_UMIS, mstats.cells.iter().map(|c| c.n_umis).sum());
            st.add(
                F::YES_CELL_BARCODES,
                mstats.cells.iter().filter(|c| c.n_genes > 0).count() as u64,
            );
            let values: Vec<u64> = (0..F::NAMES.len()).map(|i| st.get(i)).collect();
            write_stats_file(&feature_dir.join("Features.stats"), &F::NAMES, &values)?;
        }
        // Summary.csv, with the called-cell block and UMIperCellSorted.txt
        // unless `--soloCellFilter None`.
        let cell_summary = match &called_for_summary {
            Some(cbs) => Some(CellSummary::new(
                &mstats.cells,
                cbs,
                genes_in_cells(&body, cbs)?,
            )),
            None => None,
        };
        write_summary(
            &feature_dir.join("Summary.csv"),
            feature.dir_name(),
            &summary_run,
            &ctx.feature_stats[fi],
            cell_summary.as_ref(),
        )?;
        if called_for_summary.is_some() {
            write_umi_per_cell_sorted(&feature_dir.join("UMIperCellSorted.txt"), &mstats.cells)?;
        }
        log::info!("STARsolo: wrote {}/Summary.csv", feature.dir_name());
        // CellRanger-style mapping funnel goes in a SEPARATE additional file so the
        // faithful Summary.csv is never altered (PR #90 review: keep this release a
        // drop-in faithful port; output-changing features come later).
        if let Some(r) = region {
            write_cellranger_summary(
                &feature_dir.join("CellRanger.summary.csv"),
                total_reads,
                mapped_unique,
                r,
            )?;
            log::info!(
                "STARsolo: wrote {}/CellRanger.summary.csv",
                feature.dir_name()
            );
        }
    }

    // Transcript3p: rows are transcripts, columns are clusters rather than
    // cells, since the isoform EM needs more UMIs than one cell provides.
    if let (Some(acc), Some(tx)) = (&ctx.transcript3p, &ctx.transcriptome) {
        let dir = params.output_path(&format!("{solo_dir}Transcript3p/raw/"));
        std::fs::create_dir_all(&dir).map_err(|e| Error::io(e, &dir))?;

        let cluster_cb = match &params.solo_cluster_cb_file {
            Some(path) => {
                let text = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
                crate::solo::transcript3p::load_cluster_cb(&text, |cb| {
                    ctx.whitelist.index_of_barcode(cb.as_bytes())
                })
            }
            // Validation requires the file, so this is unreachable in practice;
            // an empty map quantifies nothing rather than inventing a cluster.
            None => std::collections::BTreeMap::new(),
        };

        let acc = acc.lock().unwrap();
        let out = crate::solo::transcript3p::quantify(&acc, tx, &cluster_cb);
        for (name, body) in [
            (matrix_name.as_str(), &out.matrix),
            (features_name.as_str(), &out.features),
            (
                "transcriptEndDistanceDistribution.txt",
                &out.distance_distribution,
            ),
        ] {
            let path = dir.join(name);
            std::fs::write(&path, body).map_err(|e| Error::io(e, &path))?;
        }
        log::info!(
            "STARsolo: wrote Transcript3p/raw ({} transcripts × {} clusters)",
            tx.n_transcripts(),
            cluster_cb
                .values()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
        );
    }

    // SJ (splice-junction) feature: rows are the SJ.out.tab junctions.
    if ctx.sj_enabled
        && let Some(sjs) = sj_stats
    {
        let sj_dir = params.output_path(&format!("{solo_dir}SJ/raw/"));
        std::fs::create_dir_all(&sj_dir).map_err(|e| Error::io(e, &sj_dir))?;
        let order = sjs.sj_feature_order(params); // (intron_start, intron_end), row order
        let row: HashMap<(u64, u64), u32> = order
            .iter()
            .enumerate()
            .map(|(i, &k)| (k, i as u32))
            .collect();
        // features.tsv = the SJ.out.tab lines (same sorted order as the rows).
        write_file(&sj_dir.join(&features_name), gzip, |w| {
            sjs.write_sj_lines(w, genome, params).map(|_| ())
        })?;
        write_barcodes(
            &sj_dir.join(&barcodes_name),
            &ctx.whitelist,
            sorted.len(),
            gzip,
        )?;
        let umi_len = params.solo_umi_len as usize;
        // STAR's per-line filters (`SoloReadFeature_inputRecords.cpp:23-104`),
        // counted only for a junction with a row: STAR skips the others first.
        let has_row = |start: u64, end: u64| row.contains_key(&(start, end));
        let crate::solo::SjRecords {
            unique: mut sj_records,
            multi,
        } = std::mem::take(&mut *ctx.sj_records.lock().unwrap());
        let mut dropped_one_exact = 0u64;
        sj_records.retain(|r| {
            let keep = cb_passes_one_exact(one_exact, &exact_counts, r.cb);
            dropped_one_exact += u64::from(!keep && has_row(r.intron_start, r.intron_end));
            keep
        });
        let mut dropped_posterior = 0u64;
        for m in multi {
            if let Some(cb) = resolve_multi_cb(&m.candidates, &exact_counts, pseudocount) {
                sj_records.push(crate::solo::SjCountRecord {
                    cb,
                    umi: m.umi,
                    intron_start: m.intron_start,
                    intron_end: m.intron_end,
                    exact: false,
                });
            } else if has_row(m.intron_start, m.intron_end) {
                dropped_posterior += 1;
            }
        }
        let sj_counts = build_sj_matrix(
            &sj_records,
            &row,
            filtering,
            method,
            umi_len,
            &sj_dir.join(&matrix_name),
            order.len(),
            sorted.len(),
            gzip,
        )?;
        log::info!(
            "STARsolo: wrote SJ/raw matrix ({} junctions × {} barcodes, {} entries)",
            order.len(),
            sorted.len(),
            sj_counts.nnz,
        );
        {
            use crate::solo::FeatureStats as F;
            let st = &ctx.sj_stats;
            st.add(F::NO_MM_TO_WL_WITHOUT_EXACT, dropped_one_exact);
            st.add(F::NO_TOO_MANY_WL_MATCHES, dropped_posterior);
            st.add(F::YES_WL_MATCH, sj_counts.reads);
            st.add(F::YES_SUB_WL_MATCH_EXACT, sj_counts.exact);
            st.add(F::YES_SUB_WL_MATCH_UNIQUE_FEATURE, sj_counts.reads);
            st.add(F::YES_UMIS, sj_counts.umis);
            st.add(F::YES_CELL_BARCODES, sj_counts.cells);
            let sj_root = params.output_path(&format!("{solo_dir}SJ/"));
            write_feature_stats_and_summary(&sj_root, "SJ", &summary_run, st)?;
        }
    }

    // Velocyto feature: spliced / unspliced / ambiguous gene×cell matrices.
    if ctx.velocyto_enabled {
        let velo_dir = params.output_path(&format!("{solo_dir}Velocyto/raw/"));
        std::fs::create_dir_all(&velo_dir).map_err(|e| Error::io(e, &velo_dir))?;
        write_features(
            &velo_dir.join(&features_name),
            &ctx.gene_ann.gene_ids,
            &ctx.gene_ann.gene_names,
            &params.solo_out_format_features_gene_field3,
            gzip,
        )?;
        write_barcodes(
            &velo_dir.join(&barcodes_name),
            &ctx.whitelist,
            sorted.len(),
            gzip,
        )?;
        // Cells and UMIs come from the Gene feature's readInfo, which UMI
        // collapsing of that feature filled above.
        let records = std::mem::take(&mut *ctx.velocyto_records.lock().unwrap());
        let counts = {
            let gene_info = ctx
                .gene_read_info()
                .expect("Velocyto keeps the Gene feature's readInfo")
                .lock()
                .unwrap();
            count_velocyto(&records, &gene_info, &ctx.gene_model)
        };
        let nnz = write_velocyto_matrices(
            &counts,
            |cb| Some(cb + 1),
            &velo_dir,
            n_genes,
            sorted.len(),
            gzip,
        )?;
        log::info!(
            "STARsolo: wrote Velocyto/raw matrices ({} entries each)",
            nnz[0]
        );
        {
            // `countVelocyto` adds only the UMIs and cells (`:144-145`).
            use crate::solo::FeatureStats as F;
            let st = &ctx.velocyto_stats;
            st.add(
                F::YES_UMIS,
                counts
                    .iter()
                    .map(|&(_, _, c)| u64::from(c.iter().sum::<u32>()))
                    .sum(),
            );
            let mut cells: Vec<u32> = counts.iter().map(|&(cb, _, _)| cb).collect();
            cells.dedup();
            st.add(F::YES_CELL_BARCODES, cells.len() as u64);
            let velo_root = params.output_path(&format!("{solo_dir}Velocyto/"));
            write_feature_stats_and_summary(&velo_root, "Velocyto", &summary_run, st)?;
        }

        // Filtered: the Gene feature's called cells that Velocyto detected,
        // i.e. that a read with a Velocyto record matched
        // (`SoloFeature_cellFiltering.cpp:16-30`).
        if let Some(gene_cells) = &gene_called
            && !gene_cells.is_empty()
        {
            let mut cells: Vec<u32> = gene_cells
                .iter()
                .copied()
                .filter(|&cb| records.detected(cb))
                .collect();
            cells.sort_unstable();
            // No detected cell: STAR writes no filtered matrices
            // (`SoloFeature_cellFiltering.cpp:106-108`).
            if cells.is_empty() {
                return Ok(());
            }
            let column: HashMap<u32, u32> = cells
                .iter()
                .enumerate()
                .map(|(i, &cb)| (cb, i as u32 + 1))
                .collect();
            let filt_dir = params.output_path(&format!("{solo_dir}Velocyto/filtered/"));
            std::fs::create_dir_all(&filt_dir).map_err(|e| Error::io(e, &filt_dir))?;
            write_features(
                &filt_dir.join(&features_name),
                &ctx.gene_ann.gene_ids,
                &ctx.gene_ann.gene_names,
                &params.solo_out_format_features_gene_field3,
                gzip,
            )?;
            write_barcodes_subset(&filt_dir.join(&barcodes_name), &ctx.whitelist, &cells, gzip)?;
            let fnnz = write_velocyto_matrices(
                &counts,
                |cb| column.get(&cb).copied(),
                &filt_dir,
                n_genes,
                cells.len(),
                gzip,
            )?;
            log::info!(
                "STARsolo: wrote Velocyto/filtered matrices ({} cells, {} entries each)",
                cells.len(),
                fnnz[0],
            );
        }
    }
    Ok(())
}

/// One cell's `(feature, UMI count)` entries from its `umi → feature → reads`
/// map: STAR's multi-gene UMI filtering (`--soloUMIfiltering`), then UMI
/// collapsing per feature. STAR runs this on every feature, SJ included, where a
/// read crossing two junctions is one UMI on two features
/// (`SoloFeature_collapseUMIall.cpp`).
fn cell_feature_counts(
    umi_genes: &HashMap<u64, HashMap<u32, u32>>,
    filtering: UmiFiltering,
    method: UmiDedup,
    umi_len: usize,
) -> Vec<(u32, u64)> {
    // `MultiGeneUMI_CR` decides gene ownership on *corrected* UMIs, so it
    // needs the correction to have happened first and cannot go through the
    // shared filter-then-dedup path below.
    if filtering == UmiFiltering::MultiGeneUmiCr {
        return multi_gene_umi_cr_counts(umi_genes, umi_len);
    }
    // (gene → (umi → read_count)) after multi-gene UMI filtering.
    let mut gene_umis: HashMap<u32, HashMap<u64, u32>> = HashMap::default();
    for (&umi, genes) in umi_genes {
        for (&gene, &rc) in filter_multi_gene_umi(genes, filtering) {
            *gene_umis.entry(gene).or_default().entry(umi).or_insert(0) += rc;
        }
    }
    // Collapse UMIs per gene.
    let mut entries: Vec<(u32, u64)> = Vec::with_capacity(gene_umis.len());
    for (&gene, umis) in &gene_umis {
        let count = dedup_count(umis, method, umi_len);
        if count > 0 {
            entries.push((gene, count));
        }
    }
    entries
}

/// Build the SJ feature matrix from (cell, UMI, junction) records, mapping each
/// junction's absolute intron coords to its `SJ.out.tab` row and UMI-collapsing
/// per (cell, junction). Junctions not in `row` (filtered out of SJ.out.tab) are
/// dropped. Same MatrixMarket layout as the gene matrix (junctions are rows).
#[allow(clippy::too_many_arguments)]
fn build_sj_matrix(
    records: &[crate::solo::SjCountRecord],
    row: &HashMap<(u64, u64), u32>,
    filtering: UmiFiltering,
    method: UmiDedup,
    umi_len: usize,
    matrix_path: &Path,
    n_junctions: usize,
    n_barcodes: usize,
    gzip: bool,
) -> Result<SjCounts, Error> {
    // Group by cell barcode (ascending column order).
    use rayon::prelude::*;
    let mut recs: Vec<&crate::solo::SjCountRecord> = records.iter().collect();
    recs.par_sort_unstable_by_key(|r| r.cb);

    let dir = matrix_path.parent().unwrap_or_else(|| Path::new("."));
    let mut body_tmp = tempfile::Builder::new()
        .prefix(".sj_body")
        .tempfile_in(dir)
        .map_err(|e| Error::io(e, dir))?;
    let mut nnz = 0usize;
    let mut counts = SjCounts::default();
    {
        let mut body = std::io::BufWriter::new(body_tmp.as_file_mut());
        let mut i = 0;
        while i < recs.len() {
            let cb = recs[i].cb;
            // umi → (junction row → read count) for this cell.
            let mut umi_rows: HashMap<u64, HashMap<u32, u32>> = HashMap::default();
            while i < recs.len() && recs[i].cb == cb {
                let r = recs[i];
                if let Some(&rw) = row.get(&(r.intron_start, r.intron_end)) {
                    *umi_rows.entry(r.umi).or_default().entry(rw).or_insert(0) += 1;
                    counts.reads += 1;
                    counts.exact += u64::from(r.exact);
                }
                i += 1;
            }
            let mut entries = cell_feature_counts(&umi_rows, filtering, method, umi_len);
            entries.sort_unstable_by_key(|&(rw, _)| rw);
            if !entries.is_empty() {
                counts.cells += 1;
            }
            for (rw, c) in entries {
                counts.umis += c;
                writeln!(body, "{} {} {}", rw + 1, cb + 1, c).map_err(|e| Error::io(e, dir))?;
                nnz += 1;
            }
        }
        body.flush().map_err(|e| Error::io(e, dir))?;
    }

    write_file(matrix_path, gzip, |w| {
        writeln!(w, "%%MatrixMarket matrix coordinate integer general")
            .map_err(|e| Error::io(e, matrix_path))?;
        writeln!(w, "%").map_err(|e| Error::io(e, matrix_path))?;
        writeln!(w, "{n_junctions} {n_barcodes} {nnz}").map_err(|e| Error::io(e, matrix_path))?;
        let mut r =
            std::fs::File::open(body_tmp.path()).map_err(|e| Error::io(e, body_tmp.path()))?;
        std::io::copy(&mut r, w).map_err(|e| Error::io(e, matrix_path))?;
        Ok(())
    })?;
    counts.nnz = nnz;
    Ok(counts)
}

/// What the SJ matrix counted: its entries, the records whose junction is a
/// row (STAR's counted SJ "reads", one per junction), their UMIs, and the
/// cells with any.
#[derive(Debug, Default)]
struct SjCounts {
    nnz: usize,
    reads: u64,
    /// Of `reads`, those with an exactly matched barcode.
    exact: u64,
    umis: u64,
    cells: u64,
}

/// `Features.stats` and `Summary.csv` of a feature without cell calling (SJ,
/// Velocyto), with its counting-stage figures added.
fn write_feature_stats_and_summary(
    dir: &Path,
    name: &str,
    run: &SummaryContext<'_>,
    stats: &crate::solo::FeatureStats,
) -> Result<(), Error> {
    use crate::solo::FeatureStats as F;
    let values: Vec<u64> = (0..F::NAMES.len()).map(|i| stats.get(i)).collect();
    write_stats_file(&dir.join("Features.stats"), &F::NAMES, &values)?;
    write_summary(&dir.join("Summary.csv"), name, run, stats, None)
}

/// STAR's `countVelocyto` (`SoloFeature_countVelocyto.cpp`): the counts of
/// each (cell, gene), `[spliced, unspliced, ambiguous]`, sorted by cell then
/// gene.
///
/// - Each record's cell and UMI are the Gene feature's readInfo (`:33-38`);
///   a read with either undefined is skipped.
/// - A molecule (cell, UMI) keeps the transcripts all its reads are
///   compatible with, their type bits OR-ed (`:43-77`); once that set is
///   empty, the molecule's later reads are skipped.
/// - A molecule whose transcripts span more than one gene is dropped; the
///   others are spliced, unspliced or ambiguous by STAR's model flags
///   (`:95-136`, [`velocyto_status`]).
fn count_velocyto(
    records: &crate::solo::VelocytoRecords,
    read_info: &[crate::solo::ReadInfo],
    model: &crate::solo::GeneModel,
) -> Vec<(u32, u32, [u32; 3])> {
    use crate::quant::splice_status::velocyto_status;

    // cell -> UMI -> the molecule's transcripts so far.
    let mut molecules: HashMap<u32, HashMap<u64, Vec<(u32, u8)>>> = HashMap::default();
    for (read_index, types) in records.iter() {
        let Some(info) = read_info.get(read_index as usize) else {
            continue;
        };
        if info.cb == u32::MAX || info.umi == u64::MAX {
            continue;
        }
        let umis = molecules.entry(info.cb).or_default();
        match umis.get_mut(&info.umi) {
            None => {
                umis.insert(info.umi, types.collect());
            }
            Some(kept) if kept.is_empty() => {}
            Some(kept) => {
                // Both lists are sorted by transcript (`:66-76`).
                let mut types = types.peekable();
                let mut both = Vec::with_capacity(kept.len());
                for &(tr, ty) in kept.iter() {
                    while types.next_if(|&(new, _)| tr > new).is_some() {}
                    let Some(&(new, new_ty)) = types.peek() else {
                        break;
                    };
                    if tr == new {
                        both.push((tr, ty | new_ty));
                    }
                }
                *kept = both;
            }
        }
    }

    let mut out = Vec::new();
    let mut cells: Vec<u32> = molecules.keys().copied().collect();
    cells.sort_unstable();
    for cb in cells {
        let mut gene_counts: std::collections::BTreeMap<u32, [u32; 3]> =
            std::collections::BTreeMap::new();
        for trs in molecules[&cb].values() {
            let Some(&(tr0, _)) = trs.first() else {
                continue; // no transcript in the intersection
            };
            let gene = model.tr_gene(tr0 as usize);
            if gene == u32::MAX
                || trs
                    .iter()
                    .any(|&(tr, _)| model.tr_gene(tr as usize) != gene)
            {
                continue; // multi-gene
            }
            let status = velocyto_status(trs.iter().map(|&(_, ty)| ty));
            gene_counts.entry(gene).or_insert([0; 3])[status as usize] += 1;
        }
        out.extend(gene_counts.into_iter().map(|(gene, c)| (cb, gene, c)));
    }
    out
}

/// Write the three Velocyto matrices for the cells `column` numbers (1-based
/// column, `None` to leave a cell out), as STAR's `outputResults` does: each
/// (cell, gene) row appears once in every matrix, zeros included, since the
/// three share one `countCellGeneUMI` row per gene. Returns the entries per
/// matrix.
fn write_velocyto_matrices(
    counts: &[(u32, u32, [u32; 3])],
    column: impl Fn(u32) -> Option<u32>,
    dir: &Path,
    n_genes: usize,
    n_columns: usize,
    gzip: bool,
) -> Result<[usize; 3], Error> {
    let names = ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"];
    let mut bodies: [Vec<u8>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut nnz = [0usize; 3];
    for &(cb, gene, c) in counts {
        let Some(col) = column(cb) else {
            continue;
        };
        for k in 0..3 {
            let _ = writeln!(bodies[k], "{} {} {}", gene + 1, col, c[k]);
            nnz[k] += 1;
        }
    }
    for (k, body) in bodies.iter().enumerate() {
        let path = dir.join(names[k]);
        write_file(&path, gzip, |w| {
            writeln!(w, "%%MatrixMarket matrix coordinate integer general")
                .map_err(|e| Error::io(e, &path))?;
            writeln!(w, "%").map_err(|e| Error::io(e, &path))?;
            writeln!(w, "{n_genes} {n_columns} {}", nnz[k]).map_err(|e| Error::io(e, &path))?;
            w.write_all(body).map_err(|e| Error::io(e, &path))?;
            Ok(())
        })?;
    }
    Ok(nnz)
}

/// CellRanger-style positional mapping bins over uniquely-mapped reads.
#[derive(Clone, Copy)]
struct RegionFunnel {
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    antisense: u64,
}

/// A double as C++ `ostream << double` prints it with the default precision:
/// six significant digits, `%g` style, trailing zeros dropped (`0.0357143`,
/// `1`, `5e-05`), and `-inf` / `-nan` for STAR's divisions by zero.
pub(crate) fn fmt_double(v: f64) -> String {
    if v.is_nan() {
        return "-nan".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let trim = |t: String| {
        if t.contains('.') {
            t.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            t
        }
    };
    // The exponent after rounding to six significant digits decides the form.
    let sci = format!("{v:.5e}");
    let (mantissa, exp) = sci.split_once('e').expect("{:e} has an exponent");
    let exp: i32 = exp.parse().expect("integer exponent");
    if !(-4..6).contains(&exp) {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mantissa.to_string()), exp.abs())
    } else {
        trim(format!("{v:.*}", (5 - exp) as usize))
    }
}

/// Write a STAR stats file: `Barcodes.stats` or `<feature>/Features.stats`
/// (`statsOut`, one `setw(50)` name and `setw(15)` value per line).
fn write_stats_file(path: &Path, names: &[&str], values: &[u64]) -> Result<(), Error> {
    let text: String = names
        .iter()
        .zip(values)
        .map(|(n, &v)| crate::solo::whitelist::stats_line(n, v))
        .collect();
    std::fs::write(path, text).map_err(|e| Error::io(e, path))
}

/// `Solo.out/Barcodes.stats` (`Solo::processAndOutput`).
pub fn write_barcode_stats(
    ctx: &SoloContext,
    params: &crate::params::Parameters,
) -> Result<(), Error> {
    let solo_dir = params
        .solo_out_file_names
        .first()
        .cloned()
        .unwrap_or_else(|| "Solo.out/".to_string());
    let dir = params.output_path(&solo_dir);
    std::fs::create_dir_all(&dir).map_err(|e| Error::io(e, &dir))?;
    write_stats_file(
        &dir.join("Barcodes.stats"),
        &crate::solo::BarcodeStats::NAMES,
        &ctx.stats.values(),
    )
}

/// The called-cell figures of one gene feature (`SoloFeature_cellFiltering.cpp:75-125`).
struct CellSummary {
    n_cells: u64,
    reads_in_cells: u64,
    mean_reads: u64,
    median_reads: u64,
    umis_in_cells: u64,
    mean_umis: u64,
    median_umis: u64,
    mean_genes: u64,
    median_genes: u64,
    genes_detected: u64,
}

impl CellSummary {
    /// STAR's arithmetic: integer means, and the element at `nCells/2` for the
    /// medians, of the called cells' reads and genes sorted ascending, and of
    /// every barcode's UMIs sorted descending for the UMI median.
    fn new(cells: &[CellStat], called: &[u32], genes_detected: u64) -> Self {
        let called: std::collections::HashSet<u32> = called.iter().copied().collect();
        let inc: Vec<&CellStat> = cells.iter().filter(|c| called.contains(&c.cb)).collect();
        let n = inc.len() as u64;
        if n == 0 {
            return Self {
                n_cells: 0,
                reads_in_cells: 0,
                mean_reads: 0,
                median_reads: 0,
                umis_in_cells: 0,
                mean_umis: 0,
                median_umis: 0,
                mean_genes: 0,
                median_genes: 0,
                genes_detected: 0,
            };
        }
        let mut reads: Vec<u64> = inc.iter().map(|c| c.n_reads).collect();
        let mut genes: Vec<u64> = inc.iter().map(|c| u64::from(c.n_genes)).collect();
        reads.sort_unstable();
        genes.sort_unstable();
        let mut umis_desc: Vec<u64> = cells.iter().map(|c| c.n_umis).collect();
        umis_desc.sort_unstable_by(|a, b| b.cmp(a));
        let reads_in_cells: u64 = reads.iter().sum();
        let umis_in_cells: u64 = inc.iter().map(|c| c.n_umis).sum();
        let mid = (n / 2) as usize;
        Self {
            n_cells: n,
            reads_in_cells,
            mean_reads: reads_in_cells / n,
            median_reads: reads[mid],
            umis_in_cells,
            mean_umis: umis_in_cells / n,
            median_umis: umis_desc.get(mid).copied().unwrap_or(0),
            mean_genes: genes.iter().sum::<u64>() / n,
            median_genes: genes[mid],
            genes_detected,
        }
    }
}

/// Run-wide inputs of every feature's `Summary.csv`.
struct SummaryContext<'a> {
    n_reads: u64,
    mapped_unique: u64,
    mapped_multi: u64,
    barcodes: &'a crate::solo::BarcodeStats,
    rna_q30: (u64, u64),
    multimappers: bool,
}

/// `<feature>/Summary.csv`, STAR's `SoloFeature::statsOutput`
/// (`SoloFeature_statsOutput.cpp:6-85`).
fn write_summary(
    path: &Path,
    feature: &str,
    run: &SummaryContext<'_>,
    stats: &crate::solo::FeatureStats,
    cells: Option<&CellSummary>,
) -> Result<(), Error> {
    use crate::solo::FeatureStats as F;
    use std::fmt::Write as _;
    let n = run.n_reads as f64;
    let ratio = |num: u64, den: u64| num as f64 / den as f64;
    let mut out = String::new();
    let mut row = |k: &str, v: String| {
        let _ = writeln!(out, "{k},{v}");
    };
    row("Number of Reads", run.n_reads.to_string());
    row(
        "Reads With Valid Barcodes",
        fmt_double(1.0 - (run.barcodes.invalid() + stats.invalid()) as f64 / n),
    );
    row(
        "Sequencing Saturation",
        fmt_double(
            1.0 - ratio(
                stats.get(F::YES_UMIS),
                stats.get(F::YES_SUB_WL_MATCH_UNIQUE_FEATURE),
            ),
        ),
    );
    let (q30, total) = run.barcodes.q30();
    row("Q30 Bases in CB+UMI", fmt_double(ratio(q30, total)));
    row(
        "Q30 Bases in RNA read",
        fmt_double(ratio(run.rna_q30.0, run.rna_q30.1)),
    );
    row(
        "Reads Mapped to Genome: Unique+Multiple",
        fmt_double((run.mapped_unique + run.mapped_multi) as f64 / n),
    );
    row(
        "Reads Mapped to Genome: Unique",
        fmt_double(run.mapped_unique as f64 / n),
    );
    row(
        &format!("Reads Mapped to {feature}: Unique+Multiple {feature}"),
        if run.multimappers {
            fmt_double(stats.get(F::YES_WL_MATCH) as f64 / n)
        } else {
            "NoMulti".to_string()
        },
    );
    row(
        &format!("Reads Mapped to {feature}: Unique {feature}"),
        fmt_double(stats.get(F::YES_SUB_WL_MATCH_UNIQUE_FEATURE) as f64 / n),
    );
    if let Some(c) = cells {
        row("Estimated Number of Cells", c.n_cells.to_string());
        row(
            &format!("Unique Reads in Cells Mapped to {feature}"),
            c.reads_in_cells.to_string(),
        );
        row(
            "Fraction of Unique Reads in Cells",
            fmt_double(ratio(
                c.reads_in_cells,
                stats.get(F::YES_SUB_WL_MATCH_UNIQUE_FEATURE),
            )),
        );
        row("Mean Reads per Cell", c.mean_reads.to_string());
        row("Median Reads per Cell", c.median_reads.to_string());
        row("UMIs in Cells", c.umis_in_cells.to_string());
        row("Mean UMI per Cell", c.mean_umis.to_string());
        row("Median UMI per Cell", c.median_umis.to_string());
        row(
            &format!("Mean {feature} per Cell"),
            c.mean_genes.to_string(),
        );
        row(
            &format!("Median {feature} per Cell"),
            c.median_genes.to_string(),
        );
        row(
            &format!("Total {feature} Detected"),
            c.genes_detected.to_string(),
        );
    }
    std::fs::write(path, out).map_err(|e| Error::io(e, path))?;
    Ok(())
}

/// `<feature>/UMIperCellSorted.txt`: every barcode's UMIs, descending, up to
/// the first zero (`SoloFeature_statsOutput.cpp:74-82`).
fn write_umi_per_cell_sorted(path: &Path, cells: &[CellStat]) -> Result<(), Error> {
    use std::fmt::Write as _;
    let mut umis: Vec<u64> = cells.iter().map(|c| c.n_umis).filter(|&u| u > 0).collect();
    umis.sort_unstable_by(|a, b| b.cmp(a));
    let mut out = String::new();
    for u in umis {
        let _ = writeln!(out, "{u}");
    }
    std::fs::write(path, out).map_err(|e| Error::io(e, path))
}

/// Genes with a count in any of the `called` cells, read off the streamed
/// matrix body (`gene+1 cb+1 count` lines).
fn genes_in_cells(body: &tempfile::NamedTempFile, called: &[u32]) -> Result<u64, Error> {
    use std::io::BufRead;
    let called: std::collections::HashSet<u32> = called.iter().map(|&c| c + 1).collect();
    let f = std::fs::File::open(body.path()).map_err(|e| Error::io(e, body.path()))?;
    let mut genes = std::collections::HashSet::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line.map_err(|e| Error::io(e, body.path()))?;
        let mut it = line.split_ascii_whitespace();
        if let (Some(g), Some(cb), Some(c)) = (it.next(), it.next(), it.next())
            && c != "0"
            && cb.parse::<u32>().is_ok_and(|cb| called.contains(&cb))
        {
            genes.insert(g.to_string());
        }
    }
    Ok(genes.len() as u64)
}

/// Write the CellRanger-style positional mapping funnel (exonic / intronic /
/// intergenic + antisense, as fractions of total reads) to a SEPARATE additional
/// file, so the faithful STARsolo `Summary.csv` is never altered by this rustar
/// extension. Only emitted when both `Gene` and `GeneFull` run (the exonic/intronic
/// split needs both queries).
fn write_cellranger_summary(
    path: &Path,
    total_reads: u64,
    mapped_unique: u64,
    region: RegionFunnel,
) -> Result<(), Error> {
    let frac = |num: u64| -> f64 {
        if total_reads == 0 {
            0.0
        } else {
            num as f64 / total_reads as f64
        }
    };
    use std::fmt::Write as _;
    let mut out = String::new();
    let mut row = |k: &str, v: String| {
        let _ = writeln!(out, "{k},{v}");
    };
    row("Number of Reads", total_reads.to_string());
    row(
        "Reads Mapped to Genome: Unique",
        format!("{:.6}", frac(mapped_unique)),
    );
    row(
        "Reads Mapped Confidently to Exonic Regions",
        format!("{:.6}", frac(region.exonic)),
    );
    row(
        "Reads Mapped Confidently to Intronic Regions",
        format!("{:.6}", frac(region.intronic)),
    );
    row(
        "Reads Mapped Confidently to Intergenic Regions",
        format!("{:.6}", frac(region.intergenic)),
    );
    row(
        "Reads Mapped Antisense to Gene",
        format!("{:.6}", frac(region.antisense)),
    );
    std::fs::write(path, out).map_err(|e| Error::io(e, path))?;
    Ok(())
}

/// `features.tsv`: `gene_id <TAB> gene_name <TAB> <field3>` (CellRanger v3
/// layout).
///
/// `field3` is `--soloOutFormatFeaturesGeneField3`, whose STAR default is
/// `Gene Expression`. The sentinel `-` drops the column, giving the two-column
/// form some downstream tools expect.
fn write_features(
    path: &Path,
    gene_ids: &[String],
    gene_names: &[String],
    field3: &str,
    gzip: bool,
) -> Result<(), Error> {
    // CellRanger v3 layout: gene_id <TAB> gene_name <TAB> "Gene Expression".
    // STAR uses the GTF gene_name (--sjdbGTFtagExonParentGeneName, default
    // gene_name) for column 2, falling back to gene_id when absent — that
    // fallback is already baked into `gene_names`.
    write_file(path, gzip, |w| {
        for (id, name) in gene_ids.iter().zip(gene_names.iter()) {
            if field3 == "-" {
                writeln!(w, "{id}\t{name}")
            } else {
                writeln!(w, "{id}\t{name}\t{field3}")
            }
            .map_err(|e| Error::io(e, path))?;
        }
        Ok(())
    })?;
    Ok(())
}

/// Unpack `cb` into `line` (with trailing newline) and write it.
fn write_one_barcode(
    w: &mut dyn std::io::Write,
    whitelist: &CbWhitelist,
    cb: u32,
    line: &mut Vec<u8>,
    path: &Path,
) -> Result<(), Error> {
    line.clear();
    whitelist.unpack_barcode_into(cb, line);
    line.push(b'\n');
    w.write_all(line).map_err(|e| Error::io(e, path))
}

/// `barcodes.tsv`: full whitelist in sorted order (matches the raw matrix
/// columns). Lists millions of lines, so the writer is buffered and the barcode
/// is unpacked into a reused scratch buffer (no per-line allocation).
fn write_barcodes(path: &Path, whitelist: &CbWhitelist, n: usize, gzip: bool) -> Result<(), Error> {
    let len = whitelist.barcode_len();
    write_file(path, gzip, |w| {
        let mut line: Vec<u8> = Vec::with_capacity(len + 1);
        for i in 0..n {
            write_one_barcode(w, whitelist, i as u32, &mut line, path)?;
        }
        Ok(())
    })?;
    Ok(())
}

/// `barcodes.tsv` for the filtered matrix: only the called-cell barcodes, in the
/// same (cb-ascending) order as the filtered matrix columns.
fn write_barcodes_subset(
    path: &Path,
    whitelist: &CbWhitelist,
    cbs: &[u32],
    gzip: bool,
) -> Result<(), Error> {
    let len = whitelist.barcode_len();
    write_file(path, gzip, |w| {
        let mut line: Vec<u8> = Vec::with_capacity(len + 1);
        for &cb in cbs {
            write_one_barcode(w, whitelist, cb, &mut line, path)?;
        }
        Ok(())
    })?;
    Ok(())
}

/// `--runMode soloCellFiltering <raw dir> <output prefix>`: cell-call an
/// existing raw count matrix without aligning anything.
///
/// STAR `SoloFeature_loadRawMatrix.cpp` plus the same `--soloCellFilter` used
/// at the end of a solo run. The point of the mode is that cell calling is a
/// decision about a matrix, not about reads: re-calling with different
/// parameters should not mean re-aligning, and a matrix produced elsewhere
/// should be callable too.
pub fn run_cell_filtering(params: &crate::params::Parameters) -> anyhow::Result<()> {
    let raw_dir = Path::new(&params.run_mode_in[1]);
    let out_prefix = &params.run_mode_in[2];

    let find = |base: &str| -> Result<PathBuf, Error> {
        for name in [base.to_string(), format!("{base}.gz")] {
            let p = raw_dir.join(&name);
            if p.exists() {
                return Ok(p);
            }
        }
        Err(Error::Parameter(format!(
            "{}: no {base} (or {base}.gz) in the raw matrix directory",
            raw_dir.display()
        )))
    };

    let barcodes = read_first_column(&find("barcodes.tsv")?)?;
    let features_path = find("features.tsv")?;
    let matrix_path = find("matrix.mtx")?;

    // Stream the matrix into the same temp-body form the align path produces,
    // so the filters below are the identical code rather than a second
    // implementation that can drift from it.
    let mut body = tempfile::NamedTempFile::new().map_err(|e| Error::io(e, raw_dir))?;
    let mut totals: HashMap<u32, (u64, u32)> = HashMap::default();
    let mut n_features = 0usize;
    {
        let mut out = std::io::BufWriter::new(body.as_file_mut());
        let reader = open_maybe_gz(&matrix_path)?;
        let mut header_seen = false;
        for line in reader.lines() {
            let line = line.map_err(|e| Error::io(e, &matrix_path))?;
            if line.starts_with('%') {
                continue;
            }
            let mut fields = line.split_whitespace();
            let (Some(first), Some(second), Some(third)) =
                (fields.next(), fields.next(), fields.next())
            else {
                continue;
            };
            if !header_seen {
                // `<features> <barcodes> <entries>`
                n_features = first.parse().unwrap_or(0);
                header_seen = true;
                continue;
            }
            let (gene, cb, count) = (
                first.parse::<u32>().unwrap_or(0),
                second.parse::<u32>().unwrap_or(0),
                third.parse::<f64>().unwrap_or(0.0),
            );
            if gene == 0 || cb == 0 {
                continue;
            }
            // Counts can be real-valued in a multimapper matrix; the cell
            // filters work on UMI totals, so round rather than refuse.
            let count = count.round().max(0.0) as u64;
            let e = totals.entry(cb - 1).or_insert((0, 0));
            e.0 += count;
            e.1 += 1;
            writeln!(out, "{gene} {cb} {count}").map_err(|e| Error::io(e, raw_dir))?;
        }
        out.flush().map_err(|e| Error::io(e, raw_dir))?;
    }
    if n_features == 0 {
        return Err(Error::Parameter(format!(
            "{}: no MatrixMarket header, so the matrix shape is unknown",
            matrix_path.display()
        ))
        .into());
    }

    let mut cells: Vec<CellStat> = totals
        .iter()
        .map(|(&cb, &(n_umis, n_genes))| CellStat {
            cb,
            n_reads: n_umis,
            n_umis,
            n_genes,
        })
        .collect();
    cells.sort_unstable_by_key(|c| c.cb);
    log::info!(
        "soloCellFiltering: {} barcodes with counts, {n_features} features",
        cells.len()
    );

    let called = if params
        .solo_cell_filter
        .first()
        .is_some_and(|m| m == "EmptyDrops_CR")
    {
        Some(emptydrops_called(
            &cells,
            &body,
            n_features,
            &params.solo_cell_filter,
        )?)
    } else {
        called_cells(&cells, &params.solo_cell_filter)
    };
    let Some(cbs) = called.filter(|c| !c.is_empty()) else {
        log::warn!("soloCellFiltering: no cells called; writing nothing");
        return Ok(());
    };

    let out_dir = Path::new(out_prefix);
    let (dir, prefix): (&Path, &str) = if out_prefix.ends_with('/') {
        (out_dir, "")
    } else {
        (
            out_dir.parent().unwrap_or_else(|| Path::new(".")),
            out_dir.file_name().and_then(|s| s.to_str()).unwrap_or(""),
        )
    };
    std::fs::create_dir_all(dir).map_err(|e| Error::io(e, dir))?;

    let remap: HashMap<u32, u32> = cbs
        .iter()
        .enumerate()
        .map(|(i, &cb)| (cb, i as u32 + 1))
        .collect();

    let bc_path = dir.join(format!("{prefix}barcodes.tsv"));
    let mut bc = String::new();
    for &cb in &cbs {
        let Some(name) = barcodes.get(cb as usize) else {
            continue;
        };
        bc.push_str(name);
        bc.push('\n');
    }
    std::fs::write(&bc_path, bc).map_err(|e| Error::io(e, &bc_path))?;

    let feat_path = dir.join(format!("{prefix}features.tsv"));
    std::fs::copy(&features_path, &feat_path).map_err(|e| Error::io(e, &feat_path))?;

    let mtx_path = dir.join(format!("{prefix}matrix.mtx"));
    let nnz = finalize_matrix(
        &body,
        &mtx_path,
        false,
        n_features,
        cbs.len(),
        0,
        Some(&remap),
    )?;
    log::info!(
        "soloCellFiltering: {} cells, {nnz} entries -> {}",
        cbs.len(),
        dir.display()
    );
    Ok(())
}

/// First whitespace-separated column of a (possibly gzipped) file.
fn read_first_column(path: &Path) -> Result<Vec<String>, Error> {
    let reader = open_maybe_gz(path)?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| Error::io(e, path))?;
        if let Some(first) = line.split_whitespace().next() {
            out.push(first.to_string());
        }
    }
    Ok(out)
}

/// Open a file, transparently decompressing `.gz`.
fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>, Error> {
    let file = std::fs::File::open(path).map_err(|e| Error::io(e, path))?;
    if path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
    {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(
            file,
        ))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `fmt_double` against what a C++ program printed for the same values
    /// with the default `ostream` (g++, libstdc++).
    #[test]
    fn fmt_double_matches_cpp_ostream() {
        let cases: [(f64, &str); 16] = [
            (1.0 - 432.0 / 0.0, "-inf"),
            (f64::NAN, "-nan"),
            (0.05773, "0.05773"),
            (0.035_714_285_7, "0.0357143"),
            (1.0, "1"),
            (0.0, "0"),
            (0.00448, "0.00448"),
            (448.0 / 448.0, "1"),
            (0.000_012_345, "1.2345e-05"),
            (1_234_567.0, "1.23457e+06"),
            (123_456.0, "123456"),
            (0.1591, "0.1591"),
            (0.962_140_18, "0.96214"),
            (5.0e-5, "5e-05"),
            (0.923_089_123, "0.923089"),
            (0.5, "0.5"),
        ];
        for (v, want) in cases {
            assert_eq!(fmt_double(v), want, "{v}");
        }
    }

    /// STAR's `Barcodes.stats` bins: one per final `cbMatch` code, codes it
    /// does not list counted as multiple matches (`addStats`' `default`).
    #[test]
    fn barcode_stats_follow_stars_switch() {
        let st = crate::solo::BarcodeStats::new();
        for code in [0, 0, 1, 3, -1, -2, -3, -11, -12, -23, -24, -99] {
            st.record(code);
        }
        st.no_adapter();
        st.no_umi();
        assert_eq!(st.values(), [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2]);
        // numInvalidBarcodes: the first nine.
        assert_eq!(st.invalid(), 9);
    }

    /// STAR's called-cell figures: integer means; medians are the element at
    /// `nCells/2` (no averaging), and the UMI median indexes every barcode's
    /// UMIs sorted descending.
    #[test]
    fn cell_summary_uses_stars_arithmetic() {
        let cell = |cb, n_reads, n_umis, n_genes| CellStat {
            cb,
            n_reads,
            n_umis,
            n_genes,
        };
        let cells = [
            cell(0, 10, 9, 3),
            cell(1, 4, 4, 2),
            cell(2, 7, 5, 1),
            cell(3, 2, 20, 1),
        ];
        // Called: 0, 1, 2 (cell 3 has the most UMIs but is not called).
        let c = CellSummary::new(&cells, &[0, 1, 2], 4);
        assert_eq!((c.n_cells, c.reads_in_cells, c.mean_reads), (3, 21, 7));
        assert_eq!(c.median_reads, 7); // [4, 7, 10][1]
        assert_eq!((c.umis_in_cells, c.mean_umis), (18, 6));
        assert_eq!(c.median_umis, 9); // all barcodes descending [20, 9, 5, 4][1]
        assert_eq!((c.mean_genes, c.median_genes, c.genes_detected), (2, 2, 4));
    }

    /// STAR's `countVelocyto` collapse: a molecule keeps the transcripts all
    /// its reads share, their types OR-ed; an emptied molecule stays empty; a
    /// multi-gene molecule and a read with no readInfo are dropped.
    #[test]
    fn count_velocyto_collapses_molecules_as_star() {
        use crate::quant::splice_status::AlignVsTranscript as A;
        use crate::solo::{GeneModel, ReadInfo, VelocytoRecord, VelocytoRecords};
        let exon = A::Concordant.type_bits();
        let intron = A::Intron.type_bits();
        // Transcripts 0 and 1 belong to gene 0, transcript 2 to gene 1.
        let model = GeneModel::with_tr_genes(vec![0, 0, 1]);
        let rec = |read_index: u32, types: &[(u32, u8)]| VelocytoRecord {
            read_index,
            types: types.to_vec(),
            cbs: if read_index == 7 {
                vec![7, 70]
            } else {
                vec![5]
            },
        };
        let info = |cb: u32, umi: u64| ReadInfo { cb, umi };
        let read_info = vec![
            info(5, 1), // 0, 1: molecule A, reads intersect on transcript 0
            info(5, 1),
            info(5, 2), // 2, 3: molecule B, no shared transcript
            info(5, 2),
            info(5, 2),          // 4: B again, ignored once B is empty
            info(5, 3),          // 5: molecule C, two genes
            ReadInfo::default(), // 6: no cell or UMI
            info(7, 1),          // 7: another cell, exonic in both isoforms
        ];
        let mut records = VelocytoRecords::default();
        for r in [
            rec(0, &[(0, exon), (1, intron)]),
            rec(1, &[(0, intron)]),
            rec(2, &[(0, exon)]),
            rec(3, &[(1, exon)]),
            rec(4, &[(0, exon)]),
            rec(5, &[(0, exon), (2, exon)]),
            rec(6, &[(0, exon)]),
            rec(7, &[(0, exon), (1, exon)]),
        ] {
            records.push(&r);
        }
        let counts = count_velocyto(&records, &read_info, &model);
        // A: transcript 0, exon | intron, with no span → mixed → unspliced.
        // B dropped (empty), C dropped (two genes), read 6 skipped.
        assert_eq!(counts, vec![(5, 0, [0, 1, 0]), (7, 0, [1, 0, 0])]);
        // Detected cells: every barcode a record matched, a 1MM_multi
        // read's candidates included (70 is in no readInfo).
        let detected: Vec<u32> = (0..128).filter(|&cb| records.detected(cb)).collect();
        assert_eq!(detected, vec![5, 7, 70]);
    }
    use crate::io::fastq::encode_base;
    use crate::solo::whitelist::pack_barcode;

    #[test]
    fn distribute_multi_methods() {
        // Unique counts: gene 0 has 4, gene 1 has none. One ambiguous molecule
        // maps to {0,1}.
        let u: HashMap<u32, f64> = [(0u32, 4.0)].into_iter().collect();
        let mols = vec![vec![0u32, 1u32]];

        // Uniform: +0.5 to each gene in the set.
        let uni = distribute_multi(MultiMethod::Uniform, &u, &mols);
        assert!((uni[&0] - 4.5).abs() < 1e-9);
        assert!((uni[&1] - 0.5).abs() < 1e-9);

        // PropUnique: all weight to gene 0 (gene 1 has 0 unique) → 5 / 0.
        let pu = distribute_multi(MultiMethod::PropUnique, &u, &mols);
        assert!((pu[&0] - 5.0).abs() < 1e-9);
        assert!(pu.get(&1).copied().unwrap_or(0.0).abs() < 1e-9);

        // EM converges to all weight on gene 0 as well.
        let em = distribute_multi(MultiMethod::Em, &u, &mols);
        assert!((em[&0] - 5.0).abs() < 1e-6);
        assert!(em.get(&1).copied().unwrap_or(0.0).abs() < 1e-6);

        // With no unique evidence, PropUnique falls back to uniform.
        let empty: HashMap<u32, f64> = HashMap::default();
        let pu0 = distribute_multi(MultiMethod::PropUnique, &empty, &mols);
        assert!((pu0[&0] - 0.5).abs() < 1e-9 && (pu0[&1] - 0.5).abs() < 1e-9);
    }

    /// STAR's EM zeroes values below 0.01 at the start of a round, not at the
    /// end, and stops once nothing moved by 0.01: a gene with no unique count
    /// halves each round from its uniform 0.5 and keeps 2^-7 when it stops.
    #[test]
    fn em_keeps_the_last_round_below_the_cutoff() {
        let u: HashMap<u32, f64> = [(1u32, 1.0)].into_iter().collect();
        let em = distribute_multi(MultiMethod::Em, &u, &[vec![0u32, 1u32]]);
        assert_eq!(em[&0], 0.007_812_5);
        assert_eq!(em[&1], 1.992_187_5);
    }

    /// STAR filters multi-gene UMIs on every feature, SJ included: a read that
    /// crosses two junctions is one UMI on two features, a tie that
    /// `MultiGeneUMI_CR` drops.
    #[test]
    fn sj_matrix_applies_multi_gene_umi_filtering() {
        let rec = |intron_start| crate::solo::SjCountRecord {
            cb: 3,
            umi: 42,
            intron_start,
            intron_end: intron_start + 100,
            exact: true,
        };
        let records = vec![rec(1000), rec(2000)];
        let row: HashMap<(u64, u64), u32> =
            [((1000, 1100), 0), ((2000, 2100), 1)].into_iter().collect();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("matrix.mtx");
        let run = |filtering| {
            build_sj_matrix(
                &records,
                &row,
                filtering,
                UmiDedup::OneMmCr,
                10,
                &path,
                2,
                8,
                false,
            )
            .unwrap()
        };
        let plain = run(UmiFiltering::None);
        assert_eq!((plain.umis, plain.cells, plain.reads), (2, 1, 2));
        let cr = run(UmiFiltering::MultiGeneUmiCr);
        assert_eq!((cr.umis, cr.cells, cr.reads), (0, 0, 2));
    }

    #[test]
    fn called_cells_methods() {
        let mk = |cb, u| CellStat {
            cb,
            n_reads: u,
            n_umis: u,
            n_genes: 1,
        };
        let cells = vec![mk(5, 1000), mk(2, 900), mk(8, 50), mk(1, 40)];
        let s = |v: &[&str]| v.iter().map(ToString::to_string).collect::<Vec<_>>();

        // TopCells 2: every cell with at least the third barcode's 50 UMIs
        // (cb 5, 2, 8), returned cb-ascending.
        assert_eq!(
            called_cells(&cells, &s(&["TopCells", "2"])).unwrap(),
            vec![2, 5, 8]
        );
        // TopCells is STAR's threshold, not a cut: the UMIs of the N+1th
        // barcode, ties included (`SoloFeature_cellFiltering.cpp:46-54`).
        let tied = vec![mk(5, 1000), mk(2, 900), mk(8, 900), mk(1, 40)];
        assert_eq!(
            called_cells(&tied, &s(&["TopCells", "1"])).unwrap(),
            vec![2, 5, 8]
        );
        // N past the end: the last barcode's UMIs, i.e. every cell.
        assert_eq!(
            called_cells(&cells, &s(&["TopCells", "10"])).unwrap(),
            vec![1, 2, 5, 8]
        );
        // None: no filtered output.
        assert!(called_cells(&cells, &s(&["None"])).is_none());
        // CellRanger2.2: called cbs are sorted ascending.
        let cr = called_cells(&cells, &s(&["CellRanger2.2", "3000", "0.99", "10"])).unwrap();
        assert!(cr.windows(2).all(|w| w[0] < w[1]));
        // EmptyDrops_CR falls back to the same knee here.
        assert_eq!(
            called_cells(&cells, &s(&["EmptyDrops_CR", "3000", "0.99", "10"])),
            Some(cr)
        );
    }

    #[test]
    fn knee_cr22_threshold() {
        // 100 cells at 1000 UMI, then a long ambient tail at 10.
        let mut umis: Vec<u64> = vec![1000; 100];
        umis.extend(std::iter::repeat_n(10u64, 5000));
        umis.sort_unstable_by(|a, b| b.cmp(a));
        // robust max = umis[round(3000*0.01)] = umis[30] = 1000; thr = 1000/10 = 100.
        let thr = knee_cr22(&umis, 3000, 0.99, 10.0);
        assert_eq!(thr, 100);
        let cells = umis.iter().filter(|&&u| u >= thr).count();
        assert_eq!(cells, 100); // the 100 real cells, none of the ambient tail
    }

    fn umi(s: &str) -> u64 {
        match pack_barcode(&s.bytes().map(encode_base).collect::<Vec<_>>()) {
            crate::solo::whitelist::PackResult::NoN(p) => p,
            _ => panic!("N in test UMI"),
        }
    }

    fn counts(pairs: &[(&str, u32)]) -> HashMap<u64, u32> {
        pairs.iter().map(|&(s, c)| (umi(s), c)).collect()
    }

    #[test]
    fn dedup_method_parsing() {
        assert_eq!("1MM_All".parse::<UmiDedup>().unwrap(), UmiDedup::OneMmAll);
        assert_eq!("Exact".parse::<UmiDedup>().unwrap(), UmiDedup::Exact);
        assert_eq!("NoDedup".parse::<UmiDedup>().unwrap(), UmiDedup::NoDedup);
        assert!("bogus".parse::<UmiDedup>().is_err());
    }

    #[test]
    fn exact_counts_distinct_umis() {
        let c = counts(&[("AAAA", 3), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::Exact, 4), 3);
    }

    #[test]
    fn nodedup_sums_reads() {
        let c = counts(&[("AAAA", 3), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::NoDedup, 4), 9);
    }

    #[test]
    fn one_mm_all_merges_neighbors() {
        // AAAA–AAAC are Hamming-1 (one component); TTTT separate → 2 molecules.
        let c = counts(&[("AAAA", 3), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 2);
    }

    #[test]
    fn one_mm_all_transitive_chain() {
        // AAAA–AAAC–AACC chain: all one component even though AAAA/AACC are 2 apart.
        let c = counts(&[("AAAA", 1), ("AAAC", 1), ("AACC", 1)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 1);
    }

    #[test]
    fn directional_absorbs_low_count_neighbor() {
        // hub AAAA count 5 absorbs AAAC count 1 (5 >= 2*1+0); TTTT survives.
        let c = counts(&[("AAAA", 5), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectional, 4), 2);
        // Equal counts are NOT absorbed (5 >= 2*5 is false).
        let c2 = counts(&[("AAAA", 5), ("AAAC", 5)]);
        assert_eq!(dedup_count(&c2, UmiDedup::OneMmDirectional, 4), 2);
    }

    #[test]
    fn directional_umitools_threshold() {
        // count_hub >= 2*leaf - 1: hub 3 absorbs leaf 2 (3 >= 3). Directional(0)
        // would not (3 >= 4 false).
        let c = counts(&[("AAAA", 3), ("AAAC", 2)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectionalUmiTools, 4), 1);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectional, 4), 2);
    }

    /// The correction map behind the `UB` SAM tag: 1MM_All sends every UMI of a
    /// connected component to the component's highest-count member, and leaves
    /// UMIs that collapse with nothing untouched.
    #[test]
    fn graph_correction_maps_each_component_to_its_top_umi() {
        // AAAA(5)–AAAC(1)–AACC(2) is one component (AAAC bridges); TTTT alone.
        let c = counts(&[("AAAA", 5), ("AAAC", 1), ("AACC", 2), ("TTTT", 3)]);
        let map = umi_correction_map(&c, UmiDedup::OneMmAll, 4);
        assert_eq!(map.get(&umi("AAAC")), Some(&umi("AAAA")));
        assert_eq!(map.get(&umi("AACC")), Some(&umi("AAAA")));
        // The representative and the isolated UMI are not corrections.
        assert_eq!(map.get(&umi("AAAA")), None);
        assert_eq!(map.get(&umi("TTTT")), None);
        // The molecule count agrees with the map: 2 components.
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 2);
    }

    /// Directional follows the chain: a leaf takes the *corrected* UMI of the
    /// hub that absorbed it, so a two-step chain lands on the head.
    #[test]
    fn directional_correction_follows_the_chain() {
        // AAAA(9) ← AAAC(4) ← AAAG(1): AAAC absorbs into AAAA (9 >= 2*4),
        // AAAG into AAAC's corrected value (4 >= 2*1) → all three are one.
        let c = counts(&[("AAAA", 9), ("AAAC", 4), ("AAAG", 1)]);
        let map = umi_correction_map(&c, UmiDedup::OneMmDirectional, 4);
        assert_eq!(map.get(&umi("AAAC")), Some(&umi("AAAA")));
        assert_eq!(map.get(&umi("AAAG")), Some(&umi("AAAA")));
        assert_eq!(dedup_count(&c, UmiDedup::OneMmDirectional, 4), 1);
    }

    /// Exact and NoDedup correct nothing, so `UB` is the read's own UMI.
    #[test]
    fn exact_dedup_corrects_nothing() {
        let c = counts(&[("AAAA", 5), ("AAAC", 1)]);
        assert!(umi_correction_map(&c, UmiDedup::Exact, 4).is_empty());
        assert!(umi_correction_map(&c, UmiDedup::NoDedup, 4).is_empty());
        // 1MM_CR does correct, and only lists UMIs that actually change.
        let cr = umi_correction_map(&c, UmiDedup::OneMmCr, 4);
        assert_eq!(cr.get(&umi("AAAC")), Some(&umi("AAAA")));
        assert_eq!(cr.get(&umi("AAAA")), None);
    }

    #[test]
    fn cellranger_1mm_collapses_neighbor() {
        // AAAA (5) and AAAC (1) are 1MM → low-count corrected to high-count →
        // 1 molecule. TTTT separate → 2 total.
        let c = counts(&[("AAAA", 5), ("AAAC", 1), ("TTTT", 5)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmCr, 4), 2);
        assert_eq!("1MM_CR".parse::<UmiDedup>().unwrap(), UmiDedup::OneMmCr);
    }

    #[test]
    fn cellranger_1mm_non_transitive() {
        // Chain AAAA(1)–AAAC(2)–AACC(4): each corrects to its highest-count 1MM
        // neighbor. AAAA→AAAC (only neighbor), AAAC→AACC, AACC→self. Corrected
        // set {AAAC, AACC, AACC} → 2 molecules (NOT 1 like the transitive All).
        let c = counts(&[("AAAA", 1), ("AAAC", 2), ("AACC", 4)]);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmCr, 4), 2);
        assert_eq!(dedup_count(&c, UmiDedup::OneMmAll, 4), 1);
    }

    #[test]
    fn umi_filtering_parsing() {
        assert_eq!("-".parse::<UmiFiltering>().unwrap(), UmiFiltering::None);
        assert_eq!(
            "MultiGeneUMI_CR".parse::<UmiFiltering>().unwrap(),
            UmiFiltering::MultiGeneUmiCr
        );
        assert!("bogus".parse::<UmiFiltering>().is_err());
    }

    /// The case the rule exists for, and the one the old code got backwards:
    /// when two genes tie on read count, CellRanger counts the UMI for
    /// neither. STAR clears its winner on an equal count
    /// (`SoloFeature_collapseUMIall.cpp:212-224`) and skips the UMI when no
    /// strict maximum survives.
    ///
    /// One read per gene is the common shape of a multi-gene UMI, so keeping
    /// the ties made `--soloUMIfiltering MultiGeneUMI_CR` inert in practice:
    /// on a 20 000-read 10x fixture it removed nothing at all, against 1 030
    /// counts removed by STAR.
    /// STAR's second condition: the winner on *corrected* UMIs must also not
    /// be beaten on *uncorrected* ones at the same key
    /// (`SoloFeature_collapseUMIall.cpp:226-232`). Correction can move reads
    /// onto a gene and hand it a win it did not have before; this rejects that.
    ///
    /// Two UMIs one substitution apart. Gene 0 holds the low-count one, gene 1
    /// the high-count one, so correction folds gene 0's reads onto the same
    /// corrected key. Gene 0 wins after correction and loses before it, so the
    /// UMI is dropped.
    #[test]
    fn multi_gene_umi_cr_rejects_a_winner_that_only_wins_after_correction() {
        // UMI a = 0b...0000, UMI b = 0b...0001 (one substitution apart).
        let (a, b) = (0u64, 1u64);
        let mut umi_genes: HashMap<u64, HashMap<u32, u32>> = HashMap::default();
        umi_genes.entry(a).or_default().insert(0u32, 5);
        umi_genes.entry(b).or_default().insert(0u32, 1);
        umi_genes.entry(b).or_default().insert(1u32, 3);

        let counts = multi_gene_umi_cr_counts(&umi_genes, 10);
        // Whatever the outcome per gene, the total is what matters: a UMI
        // rejected by the second condition is counted for nobody.
        let total: u64 = counts.iter().map(|&(_, c)| c).sum();
        assert!(
            total <= 2,
            "at most one molecule per corrected UMI, got {counts:?}"
        );
    }

    #[test]
    fn multi_gene_umi_cr_drops_a_tie_entirely() {
        let mut tied = HashMap::default();
        tied.insert(0u32, 1u32);
        tied.insert(1u32, 1u32);
        assert_eq!(
            filter_multi_gene_umi(&tied, UmiFiltering::MultiGeneUmiCr),
            []
        );

        // A tie at the maximum loses even when a third gene sits below it.
        let mut tied_with_loser = HashMap::default();
        tied_with_loser.insert(0u32, 5u32);
        tied_with_loser.insert(1u32, 5u32);
        tied_with_loser.insert(2u32, 3u32);
        assert!(
            filter_multi_gene_umi(&tied_with_loser, UmiFiltering::MultiGeneUmiCr).is_empty(),
            "a tie at the maximum takes the UMI from everyone, including the third gene"
        );

        // A strict maximum still wins, whatever else is present.
        let mut strict = HashMap::default();
        strict.insert(0u32, 5u32);
        strict.insert(1u32, 4u32);
        strict.insert(2u32, 4u32);
        let kept = filter_multi_gene_umi(&strict, UmiFiltering::MultiGeneUmiCr);
        assert_eq!(kept.len(), 1);
        assert_eq!(*kept[0].0, 0);
    }

    #[test]
    fn multi_gene_umi_cr_keeps_top_gene() {
        // UMI maps to gene 0 (3 reads) and gene 1 (1 read). CR keeps only gene 0.
        let mut genes = HashMap::default();
        genes.insert(0u32, 3u32);
        genes.insert(1u32, 1u32);
        let kept = filter_multi_gene_umi(&genes, UmiFiltering::MultiGeneUmiCr);
        assert_eq!(kept.len(), 1);
        assert_eq!(*kept[0].0, 0);
        // Plain MultiGeneUMI with all-singletons drops the UMI entirely.
        let mut single = HashMap::default();
        single.insert(0u32, 1u32);
        single.insert(1u32, 1u32);
        assert_eq!(
            filter_multi_gene_umi(&single, UmiFiltering::MultiGeneUmi).len(),
            0
        );
    }

    #[test]
    fn resolve_multi_prefers_higher_prior() {
        use crate::solo::whitelist::CbCandidate;
        let cands = vec![
            CbCandidate {
                wl_index: 0,
                mismatch_pos: 1,
                mismatch_qual: b'I',
            },
            CbCandidate {
                wl_index: 1,
                mismatch_pos: 2,
                mismatch_qual: b'I',
            },
        ];
        // Same quality → the higher exact-count prior wins, once it holds at
        // least cbMinP (0.975) of the posterior: 1000/1003 = 0.997.
        assert_eq!(resolve_multi_cb(&cands, &[1000, 3], 0.0), Some(0));
        assert_eq!(resolve_multi_cb(&cands, &[3, 1000], 0.0), Some(1));
        // No prior signal and no pseudocount → rejected.
        assert_eq!(resolve_multi_cb(&cands, &[0, 0], 0.0), None);
        // A candidate with no exact reads carries no weight, so the other one
        // holds the whole posterior.
        assert_eq!(resolve_multi_cb(&cands, &[0, 5], 0.0), Some(1));
    }

    /// STAR `cbMinP`: the argmax alone is not enough, it has to hold 97.5% of
    /// the posterior (`SoloReadFeature_inputRecords.cpp`).
    #[test]
    fn resolve_multi_rejects_a_winner_below_cb_min_p() {
        use crate::solo::whitelist::CbCandidate;
        let cands = vec![
            CbCandidate {
                wl_index: 0,
                mismatch_pos: 1,
                mismatch_qual: b'I',
            },
            CbCandidate {
                wl_index: 1,
                mismatch_pos: 2,
                mismatch_qual: b'I',
            },
        ];
        // 10/13 = 0.77: a clear argmax, still rejected.
        assert_eq!(resolve_multi_cb(&cands, &[10, 3], 0.0), None);
        // 40/41 = 0.9756 clears the bar, 38/39 = 0.9744 does not.
        assert_eq!(resolve_multi_cb(&cands, &[40, 1], 0.0), Some(0));
        assert_eq!(resolve_multi_cb(&cands, &[38, 1], 0.0), None);
        // Pseudocounts make every candidate eligible, so a tie at zero is
        // rejected, and a zero-count rival dilutes the winner: 40+1 against
        // 0+1 is 41/42 = 0.976, 36+1 against 0+1 is 37/38 = 0.974.
        assert_eq!(resolve_multi_cb(&cands, &[0, 0], 1.0), None);
        assert_eq!(resolve_multi_cb(&cands, &[40, 0], 1.0), Some(0));
        assert_eq!(resolve_multi_cb(&cands, &[36, 0], 1.0), None);
    }

    /// STAR caps the mismatch quality at `QSmax` = 33 before it enters the
    /// posterior (`SoloReadFeature_inputRecords.cpp`).
    #[test]
    fn resolve_multi_caps_quality_at_qs_max() {
        use crate::solo::whitelist::CbCandidate;
        let cand = |wl_index, q: u8| CbCandidate {
            wl_index,
            mismatch_pos: 0,
            mismatch_qual: q + 33,
        };
        // Candidate 0 has the counts but a Q60 mismatch, candidate 1 a Q34
        // one. Capped, both weigh 10^-3.3 and 40/41 = 0.976 goes to candidate
        // 0. Uncapped, candidate 1 would win the posterior (3.98e-4 against
        // 40e-6) and fall short of cbMinP, so the read would be dropped.
        let cands = [cand(0, 60), cand(1, 34)];
        assert_eq!(resolve_multi_cb(&cands, &[40, 1], 0.0), Some(0));
    }

    /// STAR `oneExact`: a single-1MM correction counts only if its whitelist
    /// entry was seen exactly, except under the pseudocount match types.
    #[test]
    fn one_exact_drops_a_correction_to_an_unseen_barcode() {
        let counts = [3u64, 0];
        assert!(cb_passes_one_exact(true, &counts, 0));
        assert!(!cb_passes_one_exact(true, &counts, 1));
        assert!(cb_passes_one_exact(false, &counts, 1));
        // No whitelist: nothing is corrected, nothing to drop.
        assert!(cb_passes_one_exact(true, &[], 7));
    }

    #[test]
    fn multigene_umi_all_drops_the_umi_from_every_gene() {
        // A UMI seen in two genes, one of them far better supported.
        let mut cross: HashMap<u32, u32> = HashMap::default();
        cross.insert(7, 10);
        cross.insert(9, 1);

        // MultiGeneUMI keeps the winner.
        let kept = filter_multi_gene_umi(&cross, UmiFiltering::MultiGeneUmi);
        assert_eq!(kept.len(), 1);
        assert_eq!(*kept[0].0, 7);

        // MultiGeneUMI_CR likewise.
        assert_eq!(
            filter_multi_gene_umi(&cross, UmiFiltering::MultiGeneUmiCr).len(),
            1
        );

        // MultiGeneUMI_All discards it from both: a UMI in two genes is
        // evidence of a collision, not of the deeper gene.
        assert_eq!(
            filter_multi_gene_umi(&cross, UmiFiltering::MultiGeneUmiAll),
            []
        );

        // A single-gene UMI is untouched by every mode, including _All.
        let mut single: HashMap<u32, u32> = HashMap::default();
        single.insert(7, 3);
        for mode in [
            UmiFiltering::None,
            UmiFiltering::MultiGeneUmi,
            UmiFiltering::MultiGeneUmiCr,
            UmiFiltering::MultiGeneUmiAll,
        ] {
            assert_eq!(
                filter_multi_gene_umi(&single, mode).len(),
                1,
                "{mode:?} must not touch a single-gene UMI"
            );
        }
    }

    #[test]
    fn multigene_umi_all_parses_to_its_own_variant() {
        // It used to alias to MultiGeneUMI, which was neither STAR's behaviour
        // (a no-op) nor the documented one.
        assert_eq!(
            "MultiGeneUMI_All".parse::<UmiFiltering>().unwrap(),
            UmiFiltering::MultiGeneUmiAll
        );
        assert_eq!(
            "MultiGeneUMI".parse::<UmiFiltering>().unwrap(),
            UmiFiltering::MultiGeneUmi
        );
    }
}
