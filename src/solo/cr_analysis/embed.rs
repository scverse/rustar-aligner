//! Hook for the 2-D embeddings (`analysis/tsne/` and `analysis/umap/`).
//!
//! t-SNE and UMAP come from [`crate::solo::cr_embed`], run on the PCA
//! projection (one `[x, y]` per cell, in the order of the projection rows), and
//! [`compute_embeddings`] writes
//! `tsne/gene_expression_2_components/projection.csv` (`Barcode,TSNE-1,TSNE-2`)
//! and `umap/gene_expression_2_components/projection.csv`
//! (`Barcode,UMAP-1,UMAP-2`).

use super::{create, fmt_f64, pca::PcaResult, write_line};
use crate::error::Error;
use std::io::Write;
use std::path::Path;

fn rows(projection: &[f64], n_components: usize) -> Vec<Vec<f64>> {
    projection
        .chunks(n_components)
        .map(<[f64]>::to_vec)
        .collect()
}

/// t-SNE of the PCA projection (`n_cells x n_components`, row-major).
pub fn tsne(projection: &[f64], n_components: usize) -> Option<Vec<[f64; 2]>> {
    (n_components > 0 && !projection.is_empty())
        .then(|| crate::solo::cr_embed::tsne_2d(&rows(projection, n_components)))
}

/// UMAP of the PCA projection (`n_cells x n_components`, row-major).
pub fn umap(projection: &[f64], n_components: usize) -> Option<Vec<[f64; 2]>> {
    (n_components > 0 && !projection.is_empty())
        .then(|| crate::solo::cr_embed::umap_2d(&rows(projection, n_components)))
}

fn write_embedding(
    barcodes: &[String],
    xy: &[[f64; 2]],
    prefix: &str,
    path: &Path,
) -> Result<(), Error> {
    let mut w = create(path)?;
    let io = |e| Error::io(e, path);
    write_line(&mut w, &format!("Barcode,{prefix}-1,{prefix}-2")).map_err(io)?;
    for (b, p) in barcodes.iter().zip(xy) {
        write_line(&mut w, &format!("{b},{},{}", fmt_f64(p[0]), fmt_f64(p[1]))).map_err(io)?;
    }
    w.flush().map_err(io)
}

/// Compute and write whichever embeddings are available.
pub fn compute_embeddings(
    pca: &PcaResult,
    barcodes: &[String],
    analysis_dir: &Path,
) -> Result<(), Error> {
    let k = pca.n_components();
    if k == 0 {
        return Ok(());
    }
    let proj = pca.projection_rows();
    if let Some(xy) = tsne(&proj, k) {
        write_embedding(
            barcodes,
            &xy,
            "TSNE",
            &analysis_dir.join("tsne/gene_expression_2_components/projection.csv"),
        )?;
    }
    if let Some(xy) = umap(&proj, k) {
        write_embedding(
            barcodes,
            &xy,
            "UMAP",
            &analysis_dir.join("umap/gene_expression_2_components/projection.csv"),
        )?;
    }
    Ok(())
}
