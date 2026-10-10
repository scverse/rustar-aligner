//! CellRanger-style secondary analysis (`outs/analysis/`).
//!
//! Independent implementation written from 10x Genomics' public
//! documentation (output file descriptions, algorithm overview) and from
//! published methods only: dispersion-based feature selection with
//! median/MAD normalisation inside mean-expression bins, truncated SVD
//! (subspace iteration, a close cousin of IRLBA) for PCA, k-means with
//! k-means++ seeding, Louvain community detection (Blondel et al. 2008) on a
//! Jaccard-weighted k-nearest-neighbour graph, and the sSeq negative
//! binomial exact test (Yu, Huber, Vitek 2013) with Benjamini-Hochberg
//! adjustment. No 10x source code was consulted. Numerical conventions
//! (which normalisation, which bin count, ...) were chosen by comparing the
//! results with the *output files* of a reference run, so agreement is close
//! but not bit-identical; see each submodule for what is exact.
//!
//! Layout written under the `analysis/` directory:
//!
//! * `pca/gene_expression_10_components/{components,dispersion,features_selected,projection,variance}.csv`
//! * `clustering/gene_expression_graphclust/clusters.csv` and
//!   `clustering/gene_expression_kmeans_{2..10}_clusters/clusters.csv`
//! * `diffexp/<same clusterings>/differential_expression.csv`
//! * `tsne/` and `umap/` through [`embed::compute_embeddings`] (a hook, see
//!   that module; currently writes nothing).

#![allow(
    clippy::many_single_char_names,
    clippy::unreadable_literal,
    clippy::similar_names,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::needless_range_loop,
    clippy::too_many_lines,
    clippy::type_complexity,
    clippy::format_push_string,
    clippy::explicit_iter_loop,
    clippy::default_trait_access,
    clippy::needless_for_each,
    clippy::cloned_instead_of_copied,
    clippy::manual_midpoint
)]

pub mod cluster;
pub mod diffexp;
pub mod embed;
pub mod pca;
mod rng;

use crate::error::Error;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

/// Number of principal components, as CellRanger's default.
pub const N_PCS: usize = 10;
/// Largest k of the k-means series.
pub const MAX_KMEANS: usize = 10;

/// A feature-by-barcode count matrix stored cell-major (CSC: one column per
/// barcode).
#[derive(Debug, Default, Clone)]
pub struct CountMatrix {
    pub feature_ids: Vec<String>,
    pub feature_names: Vec<String>,
    pub barcodes: Vec<String>,
    /// `indptr[j]..indptr[j+1]` indexes the entries of cell `j`.
    pub indptr: Vec<usize>,
    /// Feature index of every entry.
    pub rows: Vec<u32>,
    pub vals: Vec<f64>,
}

impl CountMatrix {
    pub fn n_genes(&self) -> usize {
        self.feature_ids.len()
    }
    pub fn n_cells(&self) -> usize {
        self.barcodes.len()
    }

    /// Total counts per cell.
    pub fn cell_sizes(&self) -> Vec<f64> {
        (0..self.n_cells())
            .map(|j| self.vals[self.indptr[j]..self.indptr[j + 1]].iter().sum())
            .collect()
    }

    /// Total counts per gene.
    pub fn gene_totals(&self) -> Vec<f64> {
        let mut t = vec![0.0; self.n_genes()];
        for (r, v) in self.rows.iter().zip(&self.vals) {
            t[*r as usize] += v;
        }
        t
    }

    /// Build from `(gene, cell, value)` triplets (0-based).
    pub fn from_triplets(
        feature_ids: Vec<String>,
        feature_names: Vec<String>,
        barcodes: Vec<String>,
        trip: &[(u32, u32, f64)],
    ) -> Self {
        let nc = barcodes.len();
        let mut indptr = vec![0usize; nc + 1];
        for &(_, c, _) in trip {
            indptr[c as usize + 1] += 1;
        }
        for j in 0..nc {
            indptr[j + 1] += indptr[j];
        }
        let mut fill = indptr.clone();
        let mut rows = vec![0u32; trip.len()];
        let mut vals = vec![0f64; trip.len()];
        for &(g, c, v) in trip {
            let p = fill[c as usize];
            rows[p] = g;
            vals[p] = v;
            fill[c as usize] += 1;
        }
        // Sort every column by gene so downstream passes are deterministic.
        for j in 0..nc {
            let (s, e) = (indptr[j], indptr[j + 1]);
            let mut idx: Vec<usize> = (s..e).collect();
            idx.sort_by_key(|&p| rows[p]);
            let r: Vec<u32> = idx.iter().map(|&p| rows[p]).collect();
            let v: Vec<f64> = idx.iter().map(|&p| vals[p]).collect();
            rows[s..e].copy_from_slice(&r);
            vals[s..e].copy_from_slice(&v);
        }
        Self {
            feature_ids,
            feature_names,
            barcodes,
            indptr,
            rows,
            vals,
        }
    }
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>, Error> {
    let f = File::open(path).map_err(|e| Error::io(e, path))?;
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn bad(path: &Path, msg: impl std::fmt::Display) -> Error {
    Error::io(
        std::io::Error::new(std::io::ErrorKind::InvalidData, msg.to_string()),
        path,
    )
}

/// Resolve `dir/name`, accepting either the plain or the `.gz` spelling.
fn resolve(dir: &Path, name: &str) -> PathBuf {
    let plain = dir.join(name);
    if plain.exists() {
        return plain;
    }
    let mut gz = plain.as_os_str().to_owned();
    gz.push(".gz");
    PathBuf::from(gz)
}

/// Read a `features/barcodes/matrix` MatrixMarket directory.
pub fn read_mtx_dir(
    dir: &Path,
    features: &str,
    barcodes: &str,
    matrix: &str,
) -> Result<CountMatrix, Error> {
    let fpath = resolve(dir, features);
    let mut ids = Vec::new();
    let mut names = Vec::new();
    for line in open_text(&fpath)?.lines() {
        let line = line.map_err(|e| Error::io(e, &fpath))?;
        if line.is_empty() {
            continue;
        }
        let mut c = line.split('\t');
        let id = c.next().unwrap_or_default().to_string();
        names.push(c.next().map_or_else(|| id.clone(), str::to_string));
        ids.push(id);
    }
    let bpath = resolve(dir, barcodes);
    let mut bcs = Vec::new();
    for line in open_text(&bpath)?.lines() {
        let line = line.map_err(|e| Error::io(e, &bpath))?;
        if let Some(b) = line.split('\t').next().filter(|s| !s.is_empty()) {
            bcs.push(b.to_string());
        }
    }
    let mpath = resolve(dir, matrix);
    let mut lines = open_text(&mpath)?.lines();
    let mut dims = None;
    let mut trip: Vec<(u32, u32, f64)> = Vec::new();
    for line in &mut lines {
        let line = line.map_err(|e| Error::io(e, &mpath))?;
        if line.starts_with('%') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        let a = it.next().and_then(|s| s.parse::<usize>().ok());
        let b = it.next().and_then(|s| s.parse::<usize>().ok());
        let c = it.next();
        match dims {
            None => {
                let (Some(a), Some(b)) = (a, b) else {
                    return Err(bad(&mpath, "bad dimensions line"));
                };
                dims = Some((a, b));
                trip.reserve(c.and_then(|s| s.parse().ok()).unwrap_or(0));
            }
            Some((ng, nc)) => {
                let (Some(g), Some(cc), Some(v)) = (a, b, c.and_then(|s| s.parse::<f64>().ok()))
                else {
                    return Err(bad(&mpath, "bad entry line"));
                };
                if g == 0 || cc == 0 || g > ng || cc > nc {
                    return Err(bad(&mpath, "entry out of range"));
                }
                trip.push((g as u32 - 1, cc as u32 - 1, v));
            }
        }
    }
    let (ng, nc) = dims.ok_or_else(|| bad(&mpath, "missing dimensions line"))?;
    if ng != ids.len() || nc != bcs.len() {
        return Err(bad(
            &mpath,
            "matrix dimensions disagree with features/barcodes",
        ));
    }
    Ok(CountMatrix::from_triplets(ids, names, bcs, &trip))
}

/// Format a float the way the reference CSVs do: shortest round-trip digits,
/// `NaN` for missing values and an exponent for very small magnitudes.
pub fn fmt_f64(x: f64) -> String {
    if x.is_nan() {
        "NaN".to_string()
    } else if x == 0.0 {
        "0".to_string()
    } else if x.abs() < 1e-5 || x.abs() >= 1e16 {
        format!("{x:e}")
    } else {
        format!("{x}")
    }
}

pub(crate) fn create(path: &Path) -> Result<BufWriter<File>, Error> {
    if let Some(p) = path.parent() {
        std::fs::create_dir_all(p).map_err(|e| Error::io(e, p))?;
    }
    Ok(BufWriter::new(
        File::create(path).map_err(|e| Error::io(e, path))?,
    ))
}

/// Everything the secondary analysis produced, kept for tests and hooks.
#[derive(Debug)]
pub struct Analysis {
    pub pca: pca::PcaResult,
    /// `(directory name, 1-based labels)` for graphclust then k-means 2..10.
    pub clusterings: Vec<(String, Vec<u32>)>,
}

/// Run the whole secondary analysis for `counts` and write `analysis_dir`.
///
/// Returns `Ok(None)` (writing nothing) when there are too few cells or genes
/// for a meaningful PCA.
pub fn run(counts: &CountMatrix, analysis_dir: &Path) -> Result<Option<Analysis>, Error> {
    let n_cells = counts.n_cells();
    if n_cells < 3 || counts.n_genes() == 0 {
        log::warn!("analysis: skipped (needs at least 3 cells, got {n_cells})");
        return Ok(None);
    }
    let sizes = counts.cell_sizes();
    let present: Vec<usize> = counts
        .gene_totals()
        .iter()
        .enumerate()
        .filter(|(_, t)| **t > 0.0)
        .map(|(i, _)| i)
        .collect();

    let pca = pca::run(counts, &sizes, N_PCS)?;
    let npc = pca.n_components();
    let pdir = analysis_dir
        .join("pca")
        .join(format!("gene_expression_{npc}_components"));
    pca::write(counts, &pca, &present, &pdir)?;

    let proj = pca.projection_rows();
    let mut clusterings: Vec<(String, Vec<u32>)> = Vec::new();
    let g = cluster::graph_cluster(&proj, npc);
    clusterings.push(("gene_expression_graphclust".to_string(), g));
    for k in 2..=MAX_KMEANS.min(n_cells) {
        clusterings.push((
            format!("gene_expression_kmeans_{k}_clusters"),
            cluster::kmeans(&proj, npc, k),
        ));
    }
    for (name, labels) in &clusterings {
        let cdir = analysis_dir.join("clustering").join(name);
        cluster::write_clusters(&counts.barcodes, labels, &cdir.join("clusters.csv"))?;
    }
    let de = diffexp::Context::new(counts, &sizes, &present);
    clusterings.par_iter().try_for_each(|(name, labels)| {
        let ddir = analysis_dir.join("diffexp").join(name);
        de.write(labels, &ddir.join("differential_expression.csv"))
    })?;

    embed::compute_embeddings(&pca, &counts.barcodes, analysis_dir)?;
    log::info!("analysis: wrote {}", analysis_dir.display());
    Ok(Some(Analysis { pca, clusterings }))
}

/// Convenience used by the CellRanger layout writer: read the filtered
/// MatrixMarket directory and write `analysis_dir`. Failures are logged, not
/// fatal, because the primary outputs are already on disk.
pub fn run_from_dir(
    filtered_dir: &Path,
    features: &str,
    barcodes: &str,
    matrix: &str,
    analysis_dir: &Path,
) {
    let res = read_mtx_dir(filtered_dir, features, barcodes, matrix)
        .and_then(|m| run(&m, analysis_dir).map(|_| ()));
    if let Err(e) = res {
        log::warn!("analysis: secondary analysis failed: {e}");
    }
}

pub(crate) fn write_line(w: &mut impl Write, s: &str) -> std::io::Result<()> {
    w.write_all(s.as_bytes())?;
    w.write_all(b"\n")
}

#[allow(dead_code)]
pub(crate) fn read_all(path: &Path) -> std::io::Result<String> {
    let mut s = String::new();
    File::open(path)?.read_to_string(&mut s)?;
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_blob_matrix() -> CountMatrix {
        // 30 genes x 40 cells; the first 15 genes mark the first 20 cells and
        // the rest mark the others, with a deterministic pseudo-random jitter.
        let mut rng = rng::Rng::new(7);
        let mut trip = Vec::new();
        for c in 0..40u32 {
            for g in 0..30u32 {
                let marker = (g < 15) == (c < 20);
                let lam = if marker { 6.0 } else { 0.6 };
                let v = (lam * (0.5 + rng.f64())).round();
                if v > 0.0 {
                    trip.push((g, c, v));
                }
            }
        }
        CountMatrix::from_triplets(
            (0..30).map(|g| format!("G{g}")).collect(),
            (0..30).map(|g| format!("name{g}")).collect(),
            (0..40).map(|c| format!("BC{c}-1")).collect(),
            &trip,
        )
    }

    #[test]
    fn bh_is_monotone_and_bounded() {
        let a = diffexp::bh(&[0.01, 0.04, 0.03, 0.5]);
        assert!(a.iter().all(|&x| (0.0..=1.0).contains(&x)));
        assert!((a[0] - 0.04).abs() < 1e-12);
    }

    #[test]
    fn nb_exact_is_one_for_balanced_split() {
        let p = diffexp::nb_exact_p(50, 50, 20.0, 20.0);
        assert!((p - 1.0).abs() < 1e-9, "{p}");
        assert!(diffexp::nb_exact_p(5, 95, 20.0, 20.0) < 1e-6);
    }

    #[test]
    fn separates_two_populations() {
        let m = two_blob_matrix();
        let dir = tempfile::tempdir().unwrap();
        let a = run(&m, dir.path()).unwrap().expect("analysis ran");
        let k2 = &a
            .clusterings
            .iter()
            .find(|c| c.0.ends_with("kmeans_2_clusters"))
            .unwrap()
            .1;
        assert!(k2[..20].iter().all(|&l| l == k2[0]));
        assert!(k2[20..].iter().all(|&l| l == k2[20]));
        assert_ne!(k2[0], k2[20]);
        assert!(a.pca.variance[0] > 0.3);
        for f in [
            "pca/gene_expression_10_components/projection.csv",
            "clustering/gene_expression_graphclust/clusters.csv",
            "diffexp/gene_expression_kmeans_2_clusters/differential_expression.csv",
        ] {
            assert!(dir.path().join(f).exists(), "{f}");
        }
    }

    /// Reads the reference pbmc1k filtered matrix and writes `analysis/` into
    /// the scratch directory so the result can be compared to the reference.
    #[test]
    #[ignore = "needs /Users/benjamin/rustar-parity data"]
    fn pbmc1k_reference_matrix() {
        let src = Path::new(
            "/Users/benjamin/rustar-parity/cr/ref/pbmc1k/outs/filtered_feature_bc_matrix",
        );
        let out = Path::new("/Users/benjamin/rustar-parity/cleanroom/analysis_work/out/analysis");
        let m = read_mtx_dir(src, "features.tsv", "barcodes.tsv", "matrix.mtx").unwrap();
        let a = run(&m, out).unwrap().unwrap();
        assert_eq!(a.pca.n_components(), 10);
        assert!((a.pca.variance[0] - 0.02839).abs() < 1e-4);
    }
}
