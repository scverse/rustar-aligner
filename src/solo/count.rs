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
use crate::solo::SoloContext;
use crate::solo::whitelist::CbWhitelist;
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;
use std::str::FromStr;

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
    /// CellRanger > 3.0 variant: keep only the highest-read-count gene for a
    /// multi-gene UMI (ties retained), without the all-singletons drop.
    MultiGeneUmiCr,
}

impl FromStr for UmiFiltering {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "-" | "None" => Ok(Self::None),
            // MultiGeneUMI_All behaves like MultiGeneUMI for the count matrix.
            "MultiGeneUMI" | "MultiGeneUMI_All" => Ok(Self::MultiGeneUmi),
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
    let distinct: std::collections::HashSet<u64> = corrected.into_iter().collect();
    distinct.len() as u64
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

/// A sparse gene-count matrix: `cell_genes[cell] = {gene → molecule_count}`.
struct CountMatrix {
    /// Per sorted-whitelist-cell-index → (gene_idx → deduped count).
    cell_genes: HashMap<u32, HashMap<u32, u64>>,
}

impl CountMatrix {
    /// Number of non-zero (cell, gene) entries.
    fn n_entries(&self) -> usize {
        self.cell_genes.values().map(HashMap::len).sum()
    }
}

/// Build the count matrix from a solo context's collected records.
///
/// Per cell, reads are grouped as `umi → gene → read_count`. Multi-gene UMIs
/// are then resolved per `--soloUMIfiltering`, and finally the surviving UMIs
/// of each gene are collapsed per `--soloUMIdedup`.
fn build_matrix(
    ctx: &SoloContext,
    method: UmiDedup,
    filtering: UmiFiltering,
    umi_len: usize,
    pseudocount: f64,
) -> CountMatrix {
    // cell → umi → gene → read multiplicity
    let mut cells: HashMap<u32, HashMap<u64, HashMap<u32, u32>>> = HashMap::new();

    let mut push = |cb: u32, gene: u32, umi: u64| {
        *cells
            .entry(cb)
            .or_default()
            .entry(umi)
            .or_default()
            .entry(gene)
            .or_insert(0) += 1;
    };

    for r in ctx.recorder.records.lock().unwrap().iter() {
        push(r.cb, r.gene, r.umi);
    }
    // Resolve deferred 1MM_multi cell barcodes against the exact-count prior.
    let exact_counts = ctx.whitelist.exact_count_snapshot();
    for m in ctx.recorder.multi_records.lock().unwrap().iter() {
        if let Some(cb) = resolve_multi_cb(&m.candidates, &exact_counts, pseudocount) {
            push(cb, m.gene, m.umi);
        }
    }

    let mut cell_genes: HashMap<u32, HashMap<u32, u64>> = HashMap::new();
    for (cb, umi_genes) in &cells {
        // (gene → (umi → read_count)) after multi-gene UMI filtering.
        let mut gene_umis: HashMap<u32, HashMap<u64, u32>> = HashMap::new();
        for (&umi, genes) in umi_genes {
            for (&gene, &rc) in filter_multi_gene_umi(genes, filtering) {
                *gene_umis.entry(gene).or_default().entry(umi).or_insert(0) += rc;
            }
        }
        for (gene, umis) in &gene_umis {
            let count = dedup_count(umis, method, umi_len);
            if count > 0 {
                cell_genes.entry(*cb).or_default().insert(*gene, count);
            }
        }
    }

    CountMatrix { cell_genes }
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
        // CellRanger > 3.0: keep the highest-read-count gene(s); no singleton drop.
        UmiFiltering::MultiGeneUmiCr => genes.iter().filter(|&(_, &rc)| rc >= max).collect(),
        UmiFiltering::None => unreachable!(),
    }
}

/// Write the raw gene-count matrix for a finished solo run. No-op (with a
/// warning) when there is no explicit whitelist, which 14.4 does not support.
pub fn write_gene_matrix(
    ctx: &SoloContext,
    params: &crate::params::Parameters,
) -> Result<(), Error> {
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

    let matrix = build_matrix(ctx, method, filtering, umi_len, pseudocount);

    // Output directory: {prefix}{soloOutFileNames[0]}Gene/raw/
    let solo_dir = params
        .solo_out_file_names
        .first()
        .cloned()
        .unwrap_or_else(|| "Solo.out/".to_string());
    let raw_dir = params.output_path(&format!("{solo_dir}Gene/raw/"));
    std::fs::create_dir_all(&raw_dir).map_err(|e| Error::io(e, &raw_dir))?;

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

    write_features(&raw_dir.join(&features_name), &ctx.gene_ann.gene_ids)?;
    write_barcodes(&raw_dir.join(&barcodes_name), &ctx.whitelist, sorted.len())?;
    write_matrix_mtx(
        &raw_dir.join(&matrix_name),
        &matrix,
        ctx.gene_ann.gene_ids.len(),
        sorted.len(),
    )?;

    log::info!(
        "STARsolo: wrote Gene/raw matrix to {} ({} genes × {} barcodes, {} entries)",
        raw_dir.display(),
        ctx.gene_ann.gene_ids.len(),
        sorted.len(),
        matrix.n_entries(),
    );
    Ok(())
}

/// `features.tsv`: `gene_id <TAB> gene_name <TAB> "Gene Expression"` (CellRanger
/// v3 layout). We have no gene names, so the id is repeated.
fn write_features(path: &Path, gene_ids: &[String]) -> Result<(), Error> {
    let mut f = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
    for id in gene_ids {
        writeln!(f, "{id}\t{id}\tGene Expression").map_err(|e| Error::io(e, path))?;
    }
    Ok(())
}

/// `barcodes.tsv`: one barcode per line in sorted whitelist order (the same
/// order the matrix columns are indexed by).
fn write_barcodes(path: &Path, whitelist: &CbWhitelist, n: usize) -> Result<(), Error> {
    let mut f = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
    for i in 0..n {
        let bc = whitelist.barcode_string(i as u32).unwrap_or_default();
        writeln!(f, "{bc}").map_err(|e| Error::io(e, path))?;
    }
    Ok(())
}

/// `matrix.mtx`: MatrixMarket coordinate format. Header `nFeatures nBarcodes
/// nEntries`; each entry `featureIndex cellIndex count` (1-based), iterated in
/// cell (column) order for stable output.
fn write_matrix_mtx(
    path: &Path,
    matrix: &CountMatrix,
    n_features: usize,
    n_barcodes: usize,
) -> Result<(), Error> {
    let mut f = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
    writeln!(f, "%%MatrixMarket matrix coordinate integer general")
        .map_err(|e| Error::io(e, path))?;
    writeln!(f, "%").map_err(|e| Error::io(e, path))?;
    writeln!(f, "{n_features} {n_barcodes} {}", matrix.n_entries())
        .map_err(|e| Error::io(e, path))?;

    // Iterate cells in ascending sorted-whitelist order; genes ascending within.
    let mut cells: Vec<&u32> = matrix.cell_genes.keys().collect();
    cells.sort_unstable();
    for &cell in cells {
        let genes = &matrix.cell_genes[&cell];
        let mut gene_idxs: Vec<&u32> = genes.keys().collect();
        gene_idxs.sort_unstable();
        for &g in gene_idxs {
            // 1-based feature index, 1-based cell index, count.
            writeln!(f, "{} {} {}", g + 1, cell + 1, genes[&g]).map_err(|e| Error::io(e, path))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::encode_base;
    use crate::solo::whitelist::pack_barcode;

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

    #[test]
    fn multi_gene_umi_cr_keeps_top_gene() {
        // UMI maps to gene 0 (3 reads) and gene 1 (1 read). CR keeps only gene 0.
        let mut genes = HashMap::new();
        genes.insert(0u32, 3u32);
        genes.insert(1u32, 1u32);
        let kept = filter_multi_gene_umi(&genes, UmiFiltering::MultiGeneUmiCr);
        assert_eq!(kept.len(), 1);
        assert_eq!(*kept[0].0, 0);
        // Plain MultiGeneUMI with all-singletons drops the UMI entirely.
        let mut single = HashMap::new();
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
}
