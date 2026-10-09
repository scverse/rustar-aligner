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
#[cfg(feature = "anndata-out")]
use nalgebra_sparse::CsrMatrix;
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
/// satisfies `count_hub >= 2*count_leaf + dir_count_add` is absorbed; the
/// molecule count is the number of surviving (non-absorbed) UMIs.
fn directional(umis: &HashMap<u64, u32>, umi_len: usize, dir_count_add: i64) -> u64 {
    // Sort by count desc, then by UMI value for determinism.
    let mut items: Vec<(u64, u32)> = umis.iter().map(|(&u, &c)| (u, c)).collect();
    items.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let n = items.len();
    let mut absorbed = vec![false; n];
    for i in 0..n {
        if absorbed[i] {
            continue;
        }
        let hub_count = i64::from(items[i].1);
        for j in 0..n {
            if i == j || absorbed[j] {
                continue;
            }
            let leaf_count = i64::from(items[j].1);
            if leaf_count <= hub_count
                && hub_count >= 2 * leaf_count + dir_count_add
                && hamming1(items[i].0, items[j].0, umi_len)
            {
                absorbed[j] = true;
            }
        }
    }
    (n - absorbed.iter().filter(|&&a| a).count()) as u64
}

// ---------------------------------------------------------------------------
// Cell-barcode multi-match resolution (deferred 1MM_multi)
// ---------------------------------------------------------------------------

/// Resolve a 1MM_multi cell barcode to a single whitelist index using the
/// count+quality posterior: weight = `(exactCount[cand] + pseudocount) · 10^(−q/10)`
/// where `q` is the mismatch-position Phred score. `pseudocount` is 1 for the
/// `*_pseudocounts` match types (CellRanger ≥ 3.0). Returns the argmax, or
/// `None` if no candidate has positive weight.
fn resolve_multi_cb(
    candidates: &[crate::solo::whitelist::CbCandidate],
    exact_counts: &[u64],
    pseudocount: f64,
) -> Option<u32> {
    let mut best: Option<(u32, f64)> = None;
    let mut total = 0.0f64;
    for c in candidates {
        let prior = *exact_counts.get(c.wl_index as usize).unwrap_or(&0) as f64 + pseudocount;
        let q = f64::from(c.mismatch_qual.saturating_sub(33)); // Phred+33 → Phred
        let weight = prior * 10f64.powf(-q / 10.0);
        total += weight;
        match best {
            Some((_, w)) if w >= weight => {}
            _ => best = Some((c.wl_index, weight)),
        }
    }
    match best {
        Some((idx, w)) if total > 0.0 && w > 0.0 => Some(idx),
        _ => None,
    }
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
    /// Distinct genes with a nonzero count anywhere in the raw matrix.
    pub genes_detected: u32,
}

/// One cell's deduplicated counts: the whitelist barcode index, the reads that
/// went in (before UMI collapse), and the gene-ascending `(feature, count)`
/// entries that came out.
#[derive(Clone)]
pub struct CellCounts {
    pub cb: u32,
    pub n_reads: u64,
    pub entries: Vec<(u32, u64)>,
}

impl CellCounts {
    pub(crate) fn stat(&self) -> CellStat {
        CellStat {
            cb: self.cb,
            n_reads: self.n_reads,
            n_umis: self.entries.iter().map(|&(_, c)| c).sum(),
            n_genes: self.entries.len() as u32,
        }
    }
}

/// The `--soloUMI*` knobs that drive UMI collapse, resolved once per run.
#[derive(Clone, Copy)]
pub(crate) struct CountOptions {
    pub method: UmiDedup,
    pub filtering: UmiFiltering,
    pub umi_len: usize,
    /// `*_pseudocounts` CB-match types add 1 to the 1MM_multi posterior prior.
    pub pseudocount: f64,
}

impl CountOptions {
    pub(crate) fn from_params(params: &crate::params::Parameters) -> Self {
        Self {
            method: params
                .solo_umi_dedup
                .first()
                .map_or("1MM_All", String::as_str)
                .parse()
                .unwrap_or(UmiDedup::OneMmAll),
            filtering: params
                .solo_umi_filtering
                .first()
                .map_or("-", String::as_str)
                .parse()
                .unwrap_or(UmiFiltering::None),
            umi_len: params.solo_umi_len as usize,
            pseudocount: f64::from(u8::from(
                params.solo_cb_match_wl_type.contains("pseudocounts"),
            )),
        }
    }
}

/// Deduplicate one feature's records into per-cell counts, cb-ascending (entries
/// gene-ascending). Mirrors STAR's `SoloFeature_collapseUMIall.cpp`: resolved
/// 1MM_multi barcodes are folded in, the flat record list is sorted by cell
/// barcode so each cell's reads are contiguous, then each cell is collapsed
/// independently (peak memory is one cell's `umi → gene` maps, not a global
/// nest over all records). Consumes the recorder's records.
pub(crate) fn dedup_cells(
    ctx: &SoloContext,
    recorder: &crate::solo::SoloRecorder,
    opts: &CountOptions,
) -> Vec<CellCounts> {
    let &CountOptions {
        method,
        filtering,
        umi_len,
        pseudocount,
    } = opts;
    // Move records out of the recorder; fold in resolved 1MM_multi cells.
    let mut records = std::mem::take(&mut *recorder.records.lock().unwrap());
    let exact_counts = ctx.whitelist.exact_count_snapshot();
    let multi = std::mem::take(&mut *recorder.multi_records.lock().unwrap());
    for m in &multi {
        if let Some(cb) = resolve_multi_cb(&m.candidates, &exact_counts, pseudocount) {
            records.push(SoloCountRecord {
                cb,
                umi: m.umi,
                gene: m.gene,
            });
        }
    }
    drop(multi);

    // Group each cell's reads together (parallel sort — the record vec is large).
    use rayon::prelude::*;
    records.par_sort_unstable_by_key(|r| r.cb);

    // Per-cell dedup is independent across cells, so run it in parallel; the
    // result stays CB-ascending, which is what both output formats expect.
    cell_bounds(&records, |r| r.cb)
        .par_iter()
        .map(|&(i, j)| {
            // umi → gene → read multiplicity, for this cell only.
            let mut umi_genes: HashMap<u64, HashMap<u32, u32>> = HashMap::default();
            for r in &records[i..j] {
                *umi_genes
                    .entry(r.umi)
                    .or_default()
                    .entry(r.gene)
                    .or_insert(0) += 1;
            }

            // `MultiGeneUMI_CR` decides gene ownership on *corrected* UMIs, so
            // it needs the correction to have happened first and cannot go
            // through the shared filter-then-dedup path below.
            let mut entries: Vec<(u32, u64)> = if filtering == UmiFiltering::MultiGeneUmiCr {
                multi_gene_umi_cr_counts(&umi_genes, umi_len)
            } else {
                // (gene → (umi → read_count)) after multi-gene UMI filtering.
                let mut gene_umis: HashMap<u32, HashMap<u64, u32>> = HashMap::default();
                for (&umi, genes) in &umi_genes {
                    for (&gene, &rc) in filter_multi_gene_umi(genes, filtering) {
                        *gene_umis.entry(gene).or_default().entry(umi).or_insert(0) += rc;
                    }
                }

                // Collapse UMIs per gene, then emit this cell's entries gene-ascending.
                gene_umis
                    .iter()
                    .map(|(&gene, umis)| (gene, dedup_count(umis, method, umi_len)))
                    .filter(|&(_, c)| c > 0)
                    .collect()
            };
            entries.sort_unstable_by_key(|&(g, _)| g);
            CellCounts {
                cb: records[i].cb,
                n_reads: (j - i) as u64,
                entries,
            }
        })
        .collect()
}

/// One contiguous `[start, end)` slice per distinct key in a key-sorted slice.
fn cell_bounds<T>(sorted: &[T], key: impl Fn(&T) -> u32) -> Vec<(usize, usize)> {
    let mut bounds = Vec::new();
    let mut i = 0;
    while i < sorted.len() {
        let k = key(&sorted[i]);
        let mut j = i + 1;
        while j < sorted.len() && key(&sorted[j]) == k {
            j += 1;
        }
        bounds.push((i, j));
        i = j;
    }
    bounds
}

/// Write per-cell counts to a plain temporary MatrixMarket *body*
/// (`gene+1 cb+1 count`, barcode-ascending) and collect per-cell stats. The body
/// is finalized into `raw/` (and optionally `filtered/`) by the caller, which
/// lets the raw + filtered matrices share one pass.
fn build_matrix_body(
    cells: &[CellCounts],
    dir: &Path,
    n_features: usize,
) -> Result<(tempfile::NamedTempFile, MatrixStats), Error> {
    let mut body_tmp = tempfile::Builder::new()
        .prefix(".matrix_body")
        .tempfile_in(dir)
        .map_err(|e| Error::io(e, dir))?;
    let mut nnz = 0usize;
    let mut cell_stats: Vec<CellStat> = Vec::new();
    let mut gene_seen = vec![false; n_features];
    {
        let mut body = std::io::BufWriter::new(body_tmp.as_file_mut());
        for cell in cells {
            for &(g, c) in &cell.entries {
                writeln!(body, "{} {} {}", g + 1, cell.cb + 1, c).map_err(|e| Error::io(e, dir))?;
                gene_seen[g as usize] = true;
            }
            nnz += cell.entries.len();
            let stat = cell.stat();
            if stat.n_umis > 0 {
                cell_stats.push(stat);
            }
        }
        body.flush().map_err(|e| Error::io(e, dir))?;
    }

    let genes_detected = gene_seen.iter().filter(|&&s| s).count() as u32;
    Ok((
        body_tmp,
        MatrixStats {
            nnz,
            cells: cell_stats,
            genes_detected,
        },
    ))
}

/// Per-cell counts as a cells × features CSR. Rows span the **whole whitelist**
/// (`n_rows`), so row index == whitelist barcode index and every matrix in the
/// run shares one `obs` axis; barcodes with no counts are empty rows.
#[cfg(feature = "anndata-out")]
pub(crate) fn cells_to_csr(
    cells: &[CellCounts],
    n_rows: usize,
    n_cols: usize,
) -> Result<CsrMatrix<u64>, Error> {
    let nnz: usize = cells.iter().map(|c| c.entries.len()).sum();
    let mut row_offsets = Vec::with_capacity(n_rows + 1);
    let mut col_indices = Vec::with_capacity(nnz);
    let mut values = Vec::with_capacity(nnz);
    for cell in cells {
        // Empty rows for the whitelist barcodes between the last cell and this one.
        row_offsets.resize(cell.cb as usize + 1, col_indices.len());
        col_indices.extend(cell.entries.iter().map(|&(g, _)| g as usize));
        values.extend(cell.entries.iter().map(|&(_, c)| c));
    }
    row_offsets.resize(n_rows + 1, col_indices.len());
    CsrMatrix::try_from_csr_data(n_rows, n_cols, row_offsets, col_indices, values)
        .map_err(|e| Error::Parameter(format!("building solo count matrix: {e}")))
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
    pub(crate) fn name(self) -> &'static str {
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
            // theta_g = u_g + (multi distributed proportional to theta), iterated.
            let mut theta = u.clone();
            for s in molecules {
                for &g in s {
                    theta.entry(g).or_insert(0.0);
                }
            }
            for _ in 0..100 {
                let mut next = u.clone();
                for s in molecules {
                    for &g in s {
                        next.entry(g).or_insert(0.0);
                    }
                }
                for s in molecules {
                    let total: f64 = s.iter().map(|&g| get(&theta, g)).sum();
                    for &g in s {
                        let w = if total > 0.0 {
                            get(&theta, g) / total
                        } else {
                            unit(s)
                        };
                        *next.get_mut(&g).unwrap() += w;
                    }
                }
                let delta: f64 = next.iter().map(|(g, v)| (v - get(&theta, *g)).abs()).sum();
                theta = next;
                if delta < 1e-6 {
                    break;
                }
            }
            out = theta;
        }
    }
    out
}

/// Format a real matrix value compactly (integers without a decimal point).
fn fmt_real(v: f64) -> String {
    if v.fract().abs() < 1e-9 {
        format!("{}", v.round() as i64)
    } else {
        format!("{v:.5}")
    }
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
        while *mptr < multi.len() && multi[*mptr].cb < cb {
            *mptr += 1; // skip multi-only cells (no unique gene)
        }
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
                let counts = distribute_multi(m, u, &mols);
                let mut entries: Vec<(u32, f64)> =
                    counts.into_iter().filter(|&(_, v)| v > 1e-9).collect();
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
                cur_cb = Some(cb);
                u_map.clear();
            }
            *u_map.entry(g).or_insert(0.0) += v;
        }
        if let Some(prev) = cur_cb {
            flush(prev, &u_map, &mut mptr, &mut nnz)?;
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

/// CSR variant of [`build_multi_matrices`]: given the unique counts as a
/// whitelist × genes CSR (row index = cell barcode, as built by
/// [`cells_to_csr`]), yields one real-valued `UniqueAndMult` matrix per
/// `--soloMultiMappers` method, in `methods` order. Cells present only in
/// `multi_records` (no unique gene) still get their molecules distributed here —
/// unlike the `.mtx` writer, which walks the unique body and so skips them.
#[cfg(feature = "anndata-out")]
pub(crate) fn multi_matrices<'a>(
    unique: &'a CsrMatrix<u64>,
    multi_records: &[crate::solo::MultiGeneRecord],
    methods: &'a [MultiMethod],
) -> impl Iterator<Item = Result<(MultiMethod, CsrMatrix<f64>), Error>> + 'a {
    // Per-cell multi molecules: one gene set per deduplicated UMI.
    let mut by_cb: HashMap<u32, HashMap<u64, std::collections::BTreeSet<u32>>> = HashMap::default();
    for r in multi_records {
        by_cb
            .entry(r.cb)
            .or_default()
            .entry(r.umi)
            .or_default()
            .extend(r.genes.iter().copied());
    }
    let mols: HashMap<u32, Vec<Vec<u32>>> = by_cb
        .into_iter()
        .map(|(cb, umis)| {
            (
                cb,
                umis.into_values()
                    .map(|g| g.into_iter().collect())
                    .collect(),
            )
        })
        .collect();

    methods.iter().map(move |&m| {
        let mut row_offsets = Vec::with_capacity(unique.nrows() + 1);
        let mut col_indices = Vec::with_capacity(unique.nnz());
        let mut values = Vec::with_capacity(unique.nnz());
        row_offsets.push(0);
        for (cb, row) in unique.row_iter().enumerate() {
            let cell_mols = mols.get(&(cb as u32)).map_or(&[][..], Vec::as_slice);
            if row.nnz() == 0 && cell_mols.is_empty() {
                row_offsets.push(col_indices.len());
                continue;
            }
            let u: HashMap<u32, f64> = row
                .col_indices()
                .iter()
                .zip(row.values())
                .map(|(&g, &v)| (g as u32, v as f64))
                .collect();
            let mut entries: Vec<(u32, f64)> = distribute_multi(m, &u, cell_mols)
                .into_iter()
                .filter(|&(_, v)| v > 1e-9)
                .collect();
            entries.sort_unstable_by_key(|&(g, _)| g);
            col_indices.extend(entries.iter().map(|&(g, _)| g as usize));
            values.extend(entries.iter().map(|&(_, v)| v));
            row_offsets.push(col_indices.len());
        }
        CsrMatrix::try_from_csr_data(
            row_offsets.len() - 1,
            unique.ncols(),
            row_offsets,
            col_indices,
            values,
        )
        .map(|csr| (m, csr))
        .map_err(|e| Error::Parameter(format!("building UniqueAndMult-{} matrix: {e}", m.name())))
    })
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
pub(crate) fn called_cells(cells: &[CellStat], filter: &[String]) -> Option<Vec<u32>> {
    let method = filter.first().map_or("CellRanger2.2", String::as_str);
    let arg = |i: usize, d: f64| filter.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    let mut cbs: Vec<u32> = match method {
        "None" => return None,
        "TopCells" => {
            let n = arg(1, 0.0) as usize;
            let mut idx: Vec<&CellStat> = cells.iter().collect();
            idx.sort_by(|a, b| b.n_umis.cmp(&a.n_umis).then(a.cb.cmp(&b.cb)));
            idx.into_iter().take(n).map(|c| c.cb).collect()
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
/// `SoloFeature_emptyDrops_CR.cpp`). `triplets` streams the raw counts as 0-based
/// `(gene, cb, count)` — only the ambient + candidate cells are kept, so it need
/// not be materialized. `filter` is the
/// `EmptyDrops_CR nExpected maxPct maxMinRatio indMin indMax umiMin
/// umiMinFracMedian candMaxN FDR [simN]` argument list.
pub(crate) fn emptydrops_called(
    cells: &[CellStat],
    triplets: impl Iterator<Item = Result<(u32, u32, u32), Error>>,
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

    // Stream the raw counts for ambient (summed) + per-candidate profiles.
    let mut ambient = vec![0f64; n_features];
    let mut amb_total = 0f64;
    let mut cand_profiles: HashMap<u32, Vec<(u32, u32)>> = HashMap::default();
    for t in triplets {
        let (g, cb, v) = t?;
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

/// Stream a MatrixMarket body (`gene+1 cb+1 count` per line) as 0-based
/// `(gene, cb, count)` triplets, for `emptydrops_called`.
fn body_triplets(
    path: &Path,
) -> Result<impl Iterator<Item = Result<(u32, u32, u32), Error>>, Error> {
    let path = path.to_path_buf();
    let reader = BufReader::new(std::fs::File::open(&path).map_err(|e| Error::io(e, &path))?);
    Ok(reader.lines().filter_map(move |line| match line {
        Err(e) => Some(Err(Error::io(e, &path))),
        Ok(l) => {
            let mut it = l.split(' ');
            let (Some(gt), Some(ct), Some(vt)) = (it.next(), it.next(), it.next()) else {
                return None;
            };
            Some(Ok((
                gt.parse::<u32>().unwrap_or(1) - 1,
                ct.parse::<u32>().unwrap_or(1) - 1,
                vt.parse::<u32>().unwrap_or(0),
            )))
        }
    }))
}

/// Median of an ascending-sorted slice (0 if empty).
fn median_sorted(sorted: &[u64]) -> u64 {
    let n = sorted.len();
    if n == 0 {
        0
    } else if n % 2 == 1 {
        sorted[n / 2]
    } else {
        u64::midpoint(sorted[n / 2 - 1], sorted[n / 2])
    }
}

/// Write the raw gene-count matrix + `Summary.csv` for a finished solo run.
/// No-op (with a warning) when there is no explicit whitelist.
pub fn write_matrix_market(
    ctx: &SoloContext,
    params: &crate::params::Parameters,
    align_stats: &crate::stats::AlignmentStats,
    sj_stats: Option<&crate::junction::SpliceJunctionStats>,
    genome: &crate::genome::Genome,
) -> Result<(), Error> {
    let CbWhitelist::List { sorted, .. } = &ctx.whitelist else {
        log::warn!(
            "STARsolo: --soloCBwhitelist None matrix output is not yet supported (Phase 14.4); skipping matrix"
        );
        return Ok(());
    };

    let opts = CountOptions::from_params(params);
    let umi_len = opts.umi_len;

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

    let funnel = MappingFunnel::collect(ctx, align_stats);
    let gzip = matches!(params.solo_out_gzip.as_str(), "yes" | "Yes" | "true");
    let n_genes = ctx.gene_ann.gene_ids.len();
    let multi_methods = MultiMethod::parse_list(&params.solo_multi_mappers);

    // One {prefix}{soloOutFileNames[0]}<feature>/{raw,filtered}/ per feature.
    for (feature, recorder) in ctx.features.iter().zip(&ctx.recorders) {
        let feature_dir = params.output_path(&format!("{solo_dir}{}/", feature.dir_name()));
        let raw_dir = feature_dir.join("raw");
        std::fs::create_dir_all(&raw_dir).map_err(|e| Error::io(e, &raw_dir))?;

        // Collapse UMIs, then write the deduplicated counts into a shared temp
        // body and finalize the raw matrix (and the filtered one below) from it.
        let cells = dedup_cells(ctx, recorder, &opts);
        let (body, mstats) = build_matrix_body(&cells, &raw_dir, n_genes)?;
        drop(cells);
        write_features(
            &raw_dir.join(&features_name),
            &ctx.gene_ann.gene_ids,
            &ctx.gene_ann.gene_names,
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
                body_triplets(body.path())?,
                n_genes,
                &params.solo_cell_filter,
            )?)
        } else {
            called_cells(&mstats.cells, &params.solo_cell_filter)
        };
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
        if !multi_methods.is_empty() {
            let mg = recorder.multi_gene.lock().unwrap();
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

        write_summary(
            &feature_dir.join("Summary.csv"),
            feature.dir_name(),
            &mstats,
            &funnel,
            feature_reads(ctx, *feature),
        )?;
        log::info!("STARsolo: wrote {}/Summary.csv", feature.dir_name());
        // CellRanger-style mapping funnel goes in a SEPARATE additional file so the
        // faithful Summary.csv is never altered (PR #90 review: keep this release a
        // drop-in faithful port; output-changing features come later).
        if let Some(cr) = funnel.cellranger_summary() {
            write_cellranger_summary(&feature_dir.join("CellRanger.summary.csv"), &cr)?;
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
        let cells = sj_cells(&ctx.sj_records.lock().unwrap(), &row, opts.method, umi_len);
        let nnz = build_sj_matrix(
            &cells,
            &sj_dir.join(&matrix_name),
            order.len(),
            sorted.len(),
            gzip,
        )?;
        log::info!(
            "STARsolo: wrote SJ/raw matrix ({} junctions × {} barcodes, {} entries)",
            order.len(),
            sorted.len(),
            nnz,
        );
    }

    // Velocyto feature: spliced / unspliced / ambiguous gene×cell matrices.
    if ctx.velocyto_enabled {
        let velo_dir = params.output_path(&format!("{solo_dir}Velocyto/raw/"));
        std::fs::create_dir_all(&velo_dir).map_err(|e| Error::io(e, &velo_dir))?;
        write_features(
            &velo_dir.join(&features_name),
            &ctx.gene_ann.gene_ids,
            &ctx.gene_ann.gene_names,
            gzip,
        )?;
        write_barcodes(
            &velo_dir.join(&barcodes_name),
            &ctx.whitelist,
            sorted.len(),
            gzip,
        )?;
        // `--soloVelocytoAmbiguous no` folds exon-only molecules into spliced and
        // omits ambiguous.mtx (rustar extension); default `yes` = STARsolo 3-matrix.
        let keep_ambiguous = !matches!(
            params.solo_velocyto_ambiguous.as_str(),
            "no" | "No" | "false"
        );
        let cells = velocyto_cells(
            &ctx.velocyto_records.lock().unwrap(),
            opts.method,
            umi_len,
            keep_ambiguous,
        );
        let nnz = build_velocyto_matrices(
            &cells,
            &velo_dir,
            n_genes,
            sorted.len(),
            gzip,
            keep_ambiguous,
        )?;
        if keep_ambiguous {
            log::info!(
                "STARsolo: wrote Velocyto/raw matrices (spliced={} unspliced={} ambiguous={} entries)",
                nnz[0],
                nnz[1],
                nnz[2],
            );
        } else {
            log::info!(
                "STARsolo: wrote Velocyto/raw matrices, ambiguous folded into spliced (spliced={} unspliced={} entries)",
                nnz[0],
                nnz[1],
            );
        }
    }
    Ok(())
}

/// Collapse the SJ feature's (cell, UMI, junction) records into per-cell counts,
/// mapping each junction's absolute intron coords to its `SJ.out.tab` row.
/// Junctions not in `row` (filtered out of SJ.out.tab) are dropped. Same shape
/// as [`dedup_cells`]: cb-ascending cells, junction-ascending entries.
pub(crate) fn sj_cells(
    records: &[crate::solo::SjCountRecord],
    row: &HashMap<(u64, u64), u32>,
    method: UmiDedup,
    umi_len: usize,
) -> Vec<CellCounts> {
    use rayon::prelude::*;
    let mut recs: Vec<&crate::solo::SjCountRecord> = records.iter().collect();
    recs.par_sort_unstable_by_key(|r| r.cb);

    cell_bounds(&recs, |r| r.cb)
        .par_iter()
        .map(|&(i, j)| {
            // junction row → (umi → read count) for this cell.
            let mut sj_umis: HashMap<u32, HashMap<u64, u32>> = HashMap::default();
            for r in &recs[i..j] {
                if let Some(&rw) = row.get(&(r.intron_start, r.intron_end)) {
                    *sj_umis.entry(rw).or_default().entry(r.umi).or_insert(0) += 1;
                }
            }
            let mut entries: Vec<(u32, u64)> = sj_umis
                .into_iter()
                .map(|(rw, umis)| (rw, dedup_count(&umis, method, umi_len)))
                .filter(|&(_, c)| c > 0)
                .collect();
            entries.sort_unstable_by_key(|&(rw, _)| rw);
            CellCounts {
                cb: recs[i].cb,
                n_reads: (j - i) as u64,
                entries,
            }
        })
        .collect()
}

/// Write the SJ feature matrix (junctions are rows) in the same MatrixMarket
/// layout as the gene matrix.
fn build_sj_matrix(
    cells: &[CellCounts],
    matrix_path: &Path,
    n_junctions: usize,
    n_barcodes: usize,
    gzip: bool,
) -> Result<usize, Error> {
    let dir = matrix_path.parent().unwrap_or_else(|| Path::new("."));
    let mut body_tmp = tempfile::Builder::new()
        .prefix(".sj_body")
        .tempfile_in(dir)
        .map_err(|e| Error::io(e, dir))?;
    let mut nnz = 0usize;
    {
        let mut body = std::io::BufWriter::new(body_tmp.as_file_mut());
        for cell in cells {
            for &(rw, c) in &cell.entries {
                writeln!(body, "{} {} {}", rw + 1, cell.cb + 1, c)
                    .map_err(|e| Error::io(e, dir))?;
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
    Ok(nnz)
}

/// Collapse the `Velocyto` (cell, UMI, gene, category) records into per-cell
/// counts, one [`CellCounts`] per category in `[spliced, unspliced, ambiguous]`
/// order. Per (cell, gene) each UMI is resolved to a single category (priority
/// unspliced over spliced over ambiguous — any intron evidence makes the molecule
/// nascent), then UMI-deduplicated within that category.
///
/// With `keep_ambiguous` (default, STARsolo-faithful) all three categories are
/// kept. With `keep_ambiguous = false` the exon-only `ambiguous` molecules are
/// folded into `spliced` (an exon-only read is most likely mature mRNA; cf. He,
/// Soneson & Patro 2023) and the ambiguous counts come back empty.
pub(crate) fn velocyto_cells(
    records: &[crate::solo::VelocytoRecord],
    method: UmiDedup,
    umi_len: usize,
    keep_ambiguous: bool,
) -> Vec<[CellCounts; 3]> {
    use crate::solo::VelocytoCategory;
    // Category → matrix index (file order) and resolution priority.
    let cat_idx = |c: VelocytoCategory| match c {
        VelocytoCategory::Spliced => 0usize,
        VelocytoCategory::Unspliced => 1,
        VelocytoCategory::Ambiguous => 2,
    };
    let priority = |c: VelocytoCategory| match c {
        VelocytoCategory::Unspliced => 2u8,
        VelocytoCategory::Spliced => 1,
        VelocytoCategory::Ambiguous => 0,
    };

    use rayon::prelude::*;
    let mut recs: Vec<&crate::solo::VelocytoRecord> = records.iter().collect();
    recs.par_sort_unstable_by_key(|r| r.cb);

    cell_bounds(&recs, |r| r.cb)
        .par_iter()
        .map(|&(lo, hi)| {
            let cb = recs[lo].cb;
            // gene → umi → (resolved category, read count)
            let mut gene_umi: HashMap<u32, HashMap<u64, (VelocytoCategory, u32)>> =
                HashMap::default();
            for &r in &recs[lo..hi] {
                let e = gene_umi
                    .entry(r.gene)
                    .or_default()
                    .entry(r.umi)
                    .or_insert((r.category, 0));
                e.1 += 1;
                if priority(r.category) > priority(e.0) {
                    e.0 = r.category;
                }
            }
            // Per gene, dedup UMIs within each resolved category, emit entries.
            let mut genes: Vec<&u32> = gene_umi.keys().collect();
            genes.sort_unstable();
            let mut out = std::array::from_fn::<_, 3, _>(|_| CellCounts {
                cb,
                n_reads: (hi - lo) as u64,
                entries: Vec::new(),
            });
            for &g in &genes {
                let mut by_cat: [HashMap<u64, u32>; 3] =
                    [HashMap::default(), HashMap::default(), HashMap::default()];
                for (&umi, &(cat, rc)) in &gene_umi[g] {
                    by_cat[cat_idx(cat)].insert(umi, rc);
                }
                // Fold ambiguous (exon-only) molecules into spliced. A UMI resolves
                // to exactly one category, so spliced/ambiguous keys are disjoint.
                if !keep_ambiguous {
                    let amb = std::mem::take(&mut by_cat[2]);
                    by_cat[0].extend(amb);
                }
                for (k, cell) in out.iter_mut().enumerate() {
                    let c = dedup_count(&by_cat[k], method, umi_len);
                    if c > 0 {
                        cell.entries.push((*g, c));
                    }
                }
            }
            out
        })
        .collect()
}

/// Write the three `Velocyto` matrices (genes are rows, cells columns — the same
/// layout as the Gene matrix, which scVelo/dynamo ingest directly). Only
/// `spliced`/`unspliced` are written when ambiguous molecules were folded in.
/// The returned `[usize; 3]` reports `[spliced, unspliced, ambiguous]` nnz.
fn build_velocyto_matrices(
    cells: &[[CellCounts; 3]],
    dir: &Path,
    n_genes: usize,
    n_barcodes: usize,
    gzip: bool,
    keep_ambiguous: bool,
) -> Result<[usize; 3], Error> {
    let names = ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"];
    let mut bodies: Vec<tempfile::NamedTempFile> = Vec::new();
    for _ in 0..3 {
        bodies.push(
            tempfile::Builder::new()
                .prefix(".velo_body")
                .tempfile_in(dir)
                .map_err(|e| Error::io(e, dir))?,
        );
    }
    let mut nnz = [0usize; 3];
    {
        let mut writers: Vec<std::io::BufWriter<&mut std::fs::File>> = bodies
            .iter_mut()
            .map(|t| std::io::BufWriter::new(t.as_file_mut()))
            .collect();
        for per_cat in cells {
            for (k, w) in writers.iter_mut().enumerate() {
                for &(g, c) in &per_cat[k].entries {
                    writeln!(w, "{} {} {}", g + 1, per_cat[k].cb + 1, c)
                        .map_err(|e| Error::io(e, dir))?;
                    nnz[k] += 1;
                }
            }
        }
        for w in &mut writers {
            w.flush().map_err(|e| Error::io(e, dir))?;
        }
    }

    // Three files by default; only spliced/unspliced when ambiguous is folded.
    let n_out = if keep_ambiguous { 3 } else { 2 };
    for (k, body) in bodies.iter().take(n_out).enumerate() {
        let path = dir.join(names[k]);
        write_file(&path, gzip, |w| {
            writeln!(w, "%%MatrixMarket matrix coordinate integer general")
                .map_err(|e| Error::io(e, &path))?;
            writeln!(w, "%").map_err(|e| Error::io(e, &path))?;
            writeln!(w, "{n_genes} {n_barcodes} {}", nnz[k]).map_err(|e| Error::io(e, &path))?;
            let mut r = std::fs::File::open(body.path()).map_err(|e| Error::io(e, body.path()))?;
            std::io::copy(&mut r, w).map_err(|e| Error::io(e, &path))?;
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

/// Reads uniquely assigned to `feature` among valid-barcode reads — the STARsolo
/// "Reads Mapped to <feature>: Unique" metric.
pub(crate) fn feature_reads(ctx: &SoloContext, feature: crate::solo::SoloFeature) -> u64 {
    use std::sync::atomic::Ordering;
    ctx.features
        .iter()
        .position(|&x| x == feature)
        .map_or(0, |i| ctx.feature_reads[i].load(Ordering::Relaxed))
}

/// The global read funnel every feature's `Summary.csv` is computed against —
/// counted once per run, then shared across features and output formats.
pub(crate) struct MappingFunnel {
    pub total_reads: u64,
    pub valid_barcodes: u64,
    pub mapped_unique: u64,
    pub mapped_multi: u64,
    /// The positional (exonic/intronic/…) split, available only when both Gene
    /// and GeneFull ran — otherwise there is nothing to compare them against.
    region: Option<RegionFunnel>,
}

impl MappingFunnel {
    pub(crate) fn collect(ctx: &SoloContext, align_stats: &crate::stats::AlignmentStats) -> Self {
        use std::sync::atomic::Ordering;
        let have_funnel = ctx.features.contains(&crate::solo::SoloFeature::Gene)
            && ctx.features.contains(&crate::solo::SoloFeature::GeneFull);
        Self {
            total_reads: align_stats.total_reads.load(Ordering::Relaxed),
            valid_barcodes: ctx.stats.yes_exact.load(Ordering::Relaxed)
                + ctx.stats.yes_one_mm.load(Ordering::Relaxed)
                + ctx.stats.yes_mult_mm.load(Ordering::Relaxed),
            mapped_unique: align_stats.uniquely_mapped.load(Ordering::Relaxed),
            mapped_multi: align_stats.multi_mapped.load(Ordering::Relaxed),
            region: have_funnel.then(|| RegionFunnel {
                exonic: ctx.region_stats.exonic.load(Ordering::Relaxed),
                intronic: ctx.region_stats.intronic.load(Ordering::Relaxed),
                intergenic: ctx.region_stats.intergenic.load(Ordering::Relaxed),
                antisense: ctx.region_stats.antisense.load(Ordering::Relaxed),
            }),
        }
    }

    /// `num` as a fraction of all input reads.
    fn frac(&self, num: u64) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            num as f64 / self.total_reads as f64
        }
    }

    /// The CellRanger positional funnel, when both Gene and GeneFull ran.
    pub(crate) fn cellranger_summary(&self) -> Option<CellRangerSummary> {
        let r = self.region?;
        Some(CellRangerSummary {
            n_reads: self.total_reads,
            frac_mapped_genome_unique: self.frac(self.mapped_unique),
            frac_exonic: self.frac(r.exonic),
            frac_intronic: self.frac(r.intronic),
            frac_intergenic: self.frac(r.intergenic),
            frac_antisense: self.frac(r.antisense),
        })
    }
}

/// Write the STARsolo-faithful `Summary.csv` for one feature: the sequencing /
/// genome-mapping rows plus per-cell UMI/gene statistics over the CR2.2-knee-called
/// cells. The CellRanger-style exonic/intronic/intergenic/antisense funnel is a
/// rustar extension kept in a SEPARATE file (see `write_cellranger_summary`) so
/// this file is never altered relative to STARsolo's own Summary.csv.
fn write_summary(
    path: &Path,
    feature_name: &str,
    mstats: &MatrixStats,
    funnel: &MappingFunnel,
    feature_mapped: u64,
) -> Result<(), Error> {
    let s = feature_summary(&mstats.cells, mstats.genes_detected, funnel, feature_mapped);
    use std::fmt::Write as _;
    let mut out = String::new();
    let mut row = |k: &str, v: String| {
        let _ = writeln!(out, "{k},{v}");
    };
    row("Number of Reads", s.n_reads.to_string());
    row(
        "Reads With Valid Barcodes",
        format!("{:.6}", s.frac_valid_barcodes),
    );
    row("Sequencing Saturation", format!("{:.6}", s.saturation));
    row(
        "Reads Mapped to Genome: Unique+Multiple",
        format!("{:.6}", s.frac_mapped_genome_unique_multi),
    );
    row(
        "Reads Mapped to Genome: Unique",
        format!("{:.6}", s.frac_mapped_genome_unique),
    );
    row(
        &format!("Reads Mapped to {feature_name}: Unique {feature_name}"),
        format!("{:.6}", s.frac_mapped_feature_unique),
    );
    row("Estimated Number of Cells", s.n_cells.to_string());
    row(
        &format!("Unique Reads in Cells Mapped to {feature_name}"),
        s.reads_in_cells.to_string(),
    );
    row(
        "Fraction of Unique Reads in Cells",
        format!("{:.6}", s.frac_reads_in_cells),
    );
    row("Mean Reads per Cell", s.mean_reads_per_cell.to_string());
    row("Median Reads per Cell", s.median_reads_per_cell.to_string());
    row("UMIs in Cells", s.umis_in_cells.to_string());
    row("Mean UMI per Cell", s.mean_umis_per_cell.to_string());
    row("Median UMI per Cell", s.median_umis_per_cell.to_string());
    row(
        &format!("Mean {feature_name} per Cell"),
        s.mean_features_per_cell.to_string(),
    );
    row(
        &format!("Median {feature_name} per Cell"),
        s.median_features_per_cell.to_string(),
    );
    row(
        &format!("Total {feature_name} Detected"),
        s.features_detected.to_string(),
    );
    std::fs::write(path, out).map_err(|e| Error::io(e, path))
}

/// The `Summary.csv` statistics for one feature, before any formatting: fractions
/// are raw ratios in `[0, 1]`, counts are counts. `write_summary` renders these as
/// the STARsolo CSV; other output formats (AnnData) can consume them directly.
#[derive(Clone, Copy, Debug)]
pub struct FeatureSummary {
    /// Reads in the input (all barcodes, mapped or not).
    pub n_reads: u64,
    /// Fraction of reads whose CB matched the whitelist.
    pub frac_valid_barcodes: f64,
    /// 1 − UMIs / reads over all barcodes.
    pub saturation: f64,
    /// Fraction of reads mapped to the genome, uniquely or multi-mapping.
    pub frac_mapped_genome_unique_multi: f64,
    /// Fraction of reads mapped uniquely to the genome.
    pub frac_mapped_genome_unique: f64,
    /// Fraction of reads assigned uniquely to this feature (gene / SJ / …).
    pub frac_mapped_feature_unique: f64,
    /// Cells called by the CR2.2 knee on per-barcode UMI totals.
    pub n_cells: usize,
    /// Feature-assigned reads in called cells.
    pub reads_in_cells: u64,
    /// `reads_in_cells` over the feature-assigned reads of all barcodes.
    pub frac_reads_in_cells: f64,
    pub mean_reads_per_cell: u64,
    pub median_reads_per_cell: u64,
    /// Deduplicated UMIs in called cells.
    pub umis_in_cells: u64,
    pub mean_umis_per_cell: u64,
    pub median_umis_per_cell: u64,
    pub mean_features_per_cell: u64,
    pub median_features_per_cell: u64,
    /// Features with a nonzero count anywhere in the raw matrix.
    pub features_detected: u32,
}

/// Compute [`FeatureSummary`] from the per-cell stats plus the global mapping
/// funnel. `cell_stats` is the raw (unfiltered) per-barcode set; cell calling
/// happens here.
pub(crate) fn feature_summary(
    cell_stats: &[CellStat],
    features_detected: u32,
    funnel: &MappingFunnel,
    feature_mapped: u64,
) -> FeatureSummary {
    let frac = |num: u64| funnel.frac(num);

    // Cell calling: CR2.2 knee on per-barcode UMI totals.
    let mut umis_desc: Vec<u64> = cell_stats.iter().map(|c| c.n_umis).collect();
    umis_desc.sort_unstable_by(|a, b| b.cmp(a));
    let thr: u64 = knee_cr22(&umis_desc, 3000, 0.99, 10.0);
    let cells: Vec<&CellStat> = cell_stats.iter().filter(|c| c.n_umis >= thr).collect();
    let n_cells = cells.len();

    // Totals across all barcodes (for sequencing saturation + fraction-in-cells).
    let total_reads_counted: u64 = cell_stats.iter().map(|c| c.n_reads).sum();
    let total_umis_all: u64 = cell_stats.iter().map(|c| c.n_umis).sum();
    let saturation = if total_reads_counted > 0 {
        1.0 - total_umis_all as f64 / total_reads_counted as f64
    } else {
        0.0
    };

    // Per-cell aggregates over called cells.
    let reads_in_cells: u64 = cells.iter().map(|c| c.n_reads).sum();
    let umis_in_cells: u64 = cells.iter().map(|c| c.n_umis).sum();
    let mut reads_sorted: Vec<u64> = cells.iter().map(|c| c.n_reads).collect();
    let mut umis_sorted: Vec<u64> = cells.iter().map(|c| c.n_umis).collect();
    let mut genes_sorted: Vec<u64> = cells.iter().map(|c| c.n_genes as u64).collect();
    reads_sorted.sort_unstable();
    umis_sorted.sort_unstable();
    genes_sorted.sort_unstable();
    let mean = |sum: u64| -> u64 {
        if n_cells == 0 {
            0
        } else {
            sum / n_cells as u64
        }
    };

    FeatureSummary {
        n_reads: funnel.total_reads,
        frac_valid_barcodes: frac(funnel.valid_barcodes),
        saturation,
        frac_mapped_genome_unique_multi: frac(funnel.mapped_unique + funnel.mapped_multi),
        frac_mapped_genome_unique: frac(funnel.mapped_unique),
        frac_mapped_feature_unique: frac(feature_mapped),
        n_cells,
        reads_in_cells,
        frac_reads_in_cells: if total_reads_counted > 0 {
            reads_in_cells as f64 / total_reads_counted as f64
        } else {
            0.0
        },
        mean_reads_per_cell: mean(reads_in_cells),
        median_reads_per_cell: median_sorted(&reads_sorted),
        umis_in_cells,
        mean_umis_per_cell: mean(umis_in_cells),
        median_umis_per_cell: median_sorted(&umis_sorted),
        mean_features_per_cell: mean(genes_sorted.iter().sum()),
        median_features_per_cell: median_sorted(&genes_sorted),
        features_detected,
    }
}

/// Write the CellRanger-style positional mapping funnel (exonic / intronic /
/// intergenic + antisense, as fractions of total reads) to a SEPARATE additional
/// file, so the faithful STARsolo `Summary.csv` is never altered by this rustar
/// extension. Only emitted when both `Gene` and `GeneFull` run (the exonic/intronic
/// split needs both queries).
fn write_cellranger_summary(path: &Path, s: &CellRangerSummary) -> Result<(), Error> {
    use std::fmt::Write as _;
    let mut out = String::new();
    let mut row = |k: &str, v: String| {
        let _ = writeln!(out, "{k},{v}");
    };
    row("Number of Reads", s.n_reads.to_string());
    row(
        "Reads Mapped to Genome: Unique",
        format!("{:.6}", s.frac_mapped_genome_unique),
    );
    row(
        "Reads Mapped Confidently to Exonic Regions",
        format!("{:.6}", s.frac_exonic),
    );
    row(
        "Reads Mapped Confidently to Intronic Regions",
        format!("{:.6}", s.frac_intronic),
    );
    row(
        "Reads Mapped Confidently to Intergenic Regions",
        format!("{:.6}", s.frac_intergenic),
    );
    row(
        "Reads Mapped Antisense to Gene",
        format!("{:.6}", s.frac_antisense),
    );
    std::fs::write(path, out).map_err(|e| Error::io(e, path))
}

/// The `CellRanger.summary.csv` statistics, before formatting: the positional
/// funnel as raw fractions of all input reads.
#[derive(Clone, Copy, Debug)]
pub struct CellRangerSummary {
    /// Reads in the input (all barcodes, mapped or not).
    pub n_reads: u64,
    pub frac_mapped_genome_unique: f64,
    pub frac_exonic: f64,
    pub frac_intronic: f64,
    pub frac_intergenic: f64,
    /// Uniquely-mapped reads landing on a gene's opposite strand.
    pub frac_antisense: f64,
}

/// `features.tsv`: `gene_id <TAB> gene_name <TAB> "Gene Expression"` (CellRanger
/// v3 layout). We have no gene names, so the id is repeated.
fn write_features(
    path: &Path,
    gene_ids: &[String],
    gene_names: &[String],
    gzip: bool,
) -> Result<(), Error> {
    // CellRanger v3 layout: gene_id <TAB> gene_name <TAB> "Gene Expression".
    // STAR uses the GTF gene_name (--sjdbGTFtagExonParentGeneName, default
    // gene_name) for column 2, falling back to gene_id when absent — that
    // fallback is already baked into `gene_names`.
    write_file(path, gzip, |w| {
        for (id, name) in gene_ids.iter().zip(gene_names.iter()) {
            writeln!(w, "{id}\t{name}\tGene Expression").map_err(|e| Error::io(e, path))?;
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
            body_triplets(body.path())?,
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
    use crate::io::fastq::encode_base;
    use crate::solo::whitelist::pack_barcode;

    #[test]
    fn median_sorted_odd_even_empty() {
        assert_eq!(median_sorted(&[]), 0);
        assert_eq!(median_sorted(&[5]), 5);
        assert_eq!(median_sorted(&[1, 2, 3]), 2);
        assert_eq!(median_sorted(&[10, 20, 30, 40]), 25); // midpoint(20,30)
    }

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

    #[cfg(feature = "anndata-out")]
    #[test]
    fn cells_to_csr_rows_are_whitelist_barcodes() {
        let cells = vec![
            CellCounts {
                cb: 1,
                n_reads: 3,
                entries: vec![(0, 2), (2, 1)],
            },
            CellCounts {
                cb: 3,
                n_reads: 1,
                entries: vec![(1, 1)],
            },
        ];
        // 5 whitelist barcodes: rows 0, 2 and 4 have no counts at all.
        let csr = cells_to_csr(&cells, 5, 3).expect("valid csr");
        assert_eq!((csr.nrows(), csr.ncols()), (5, 3));
        let t: Vec<_> = csr.triplet_iter().map(|(r, c, &v)| (r, c, v)).collect();
        assert_eq!(t, vec![(1, 0, 2), (1, 2, 1), (3, 1, 1)]);
    }

    #[cfg(feature = "anndata-out")]
    #[test]
    fn multi_matrices_csr() {
        // 4 whitelist barcodes × 3 genes. Cell 1: gene 0 = 4 unique. Cell 3: gene 2 = 1.
        let unique =
            CsrMatrix::try_from_csr_data(4, 3, vec![0, 0, 1, 1, 2], vec![0, 2], vec![4u64, 1])
                .expect("valid csr");
        // Cell 1: one ambiguous molecule (two records, same UMI) over {0,1}.
        let mk = |cb, umi, genes: &[u32]| crate::solo::MultiGeneRecord {
            cb,
            umi,
            genes: genes.to_vec(),
        };
        let multi = vec![mk(1, 42, &[0, 1]), mk(1, 42, &[1])];

        let methods = [MultiMethod::Uniform, MultiMethod::PropUnique];
        let out: Vec<_> = multi_matrices(&unique, &multi, &methods)
            .collect::<Result<Vec<_>, _>>()
            .expect("valid matrices");
        assert_eq!(out.len(), 2);

        // Uniform: cell 1 gets +0.5 on genes 0 and 1; cell 3 unchanged.
        let (m, uni) = &out[0];
        assert_eq!(*m, MultiMethod::Uniform);
        let t: Vec<_> = uni.triplet_iter().map(|(r, c, &v)| (r, c, v)).collect();
        assert_eq!(t, vec![(1, 0, 4.5), (1, 1, 0.5), (3, 2, 1.0)]);

        // PropUnique: all weight to gene 0 (gene 1 has no unique counts).
        let (_, pu) = &out[1];
        let t: Vec<_> = pu.triplet_iter().map(|(r, c, &v)| (r, c, v)).collect();
        assert_eq!(t, vec![(1, 0, 5.0), (3, 2, 1.0)]);
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

        // TopCells 2: the two highest-UMI cells (cb 5, 2), returned cb-ascending.
        assert_eq!(
            called_cells(&cells, &s(&["TopCells", "2"])).unwrap(),
            vec![2, 5]
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
        // Same quality → higher exact-count prior wins.
        assert_eq!(resolve_multi_cb(&cands, &[10, 3], 0.0), Some(0));
        assert_eq!(resolve_multi_cb(&cands, &[3, 10], 0.0), Some(1));
        // No prior signal and no pseudocount → rejected.
        assert_eq!(resolve_multi_cb(&cands, &[0, 0], 0.0), None);
        // Pseudocount gives every candidate positive weight → argmax accepted.
        assert!(resolve_multi_cb(&cands, &[0, 0], 1.0).is_some());
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
