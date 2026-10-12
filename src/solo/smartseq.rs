//! `--soloType SmartSeq` — plate-based full-length protocols (Smart-seq2).
//!
//! There are no cell barcodes or UMIs in the reads. Each plate well is a
//! separate library given by a `--readFilesManifest` line
//! (`read1 <TAB> read2 <TAB> cellID`); the cell identity is the manifest cellID,
//! and a gene's count for a cell is the number of its uniquely-gene-assigned
//! reads (no UMI deduplication). Output mirrors the droplet path:
//! `Solo.out/Gene/raw/{matrix.mtx, barcodes.tsv (cell IDs), features.tsv}`.
//!
//! Supports both single-end manifests (`read2 = -`, read counts) and paired-end
//! (`read2` = mate-2 file, fragment counts via paired alignment).

use crate::error::Error;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// One plate-well cell from the manifest.
pub struct SmartSeqCell {
    pub read1: PathBuf,
    /// Mate-2 file for paired-end SmartSeq; `None` for single-end (`read2 = -`).
    pub read2: Option<PathBuf>,
    pub cell_id: String,
}

/// Parse a `--readFilesManifest` TSV into per-cell entries. Lines are
/// `read1 <TAB> read2 <TAB> cellID`; blank lines and `#` comments are skipped.
/// `read2 = -` is single-end; any other value is the mate-2 file (paired-end).
pub fn parse_manifest(path: &Path) -> Result<Vec<SmartSeqCell>, Error> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
    let mut cells = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            return Err(invalid(format!(
                "readFilesManifest line {}: expected 'read1<TAB>read2<TAB>cellID', got {:?}",
                lineno + 1,
                line
            )));
        }
        cells.push(SmartSeqCell {
            read1: PathBuf::from(f[0]),
            read2: (f[1] != "-").then(|| PathBuf::from(f[1])),
            cell_id: f[2].to_string(),
        });
    }
    if cells.is_empty() {
        return Err(invalid(format!(
            "readFilesManifest {} has no cell entries",
            path.display()
        )));
    }
    Ok(cells)
}

fn invalid(msg: String) -> Error {
    Error::from(std::io::Error::new(std::io::ErrorKind::InvalidInput, msg))
}

/// Per-cell, per-gene read counts for a SmartSeq run. `cells` is the manifest
/// order (the matrix column order); `counts[cell]` maps gene → read count.
pub struct SmartSeqCounts {
    pub cell_ids: Vec<String>,
    pub counts: Vec<Mutex<std::collections::HashMap<u32, u64>>>,
    pub n_genes: usize,
}

impl SmartSeqCounts {
    pub fn new(cell_ids: Vec<String>, n_genes: usize) -> Self {
        let counts = (0..cell_ids.len())
            .map(|_| Mutex::new(std::collections::HashMap::new()))
            .collect();
        Self {
            cell_ids,
            counts,
            n_genes,
        }
    }

    /// Add `+1` to (cell, gene) for one uniquely-assigned read.
    pub fn add(&self, cell: usize, gene: u32) {
        *self.counts[cell].lock().unwrap().entry(gene).or_insert(0) += 1;
    }

    /// Write `Solo.out/Gene/raw/{matrix.mtx, barcodes.tsv, features.tsv}` —
    /// genes × cells, integer read counts. `gzip` appends `.gz`.
    pub fn write_matrix(
        &self,
        raw_dir: &Path,
        gene_ids: &[String],
        gene_names: &[String],
        gzip: bool,
    ) -> Result<usize, Error> {
        std::fs::create_dir_all(raw_dir).map_err(|e| Error::io(e, raw_dir))?;

        // features.tsv (CellRanger v3 layout: id, name, "Gene Expression").
        // Column 2 is the GTF gene_name (gene_id fallback baked into gene_names).
        crate::solo::mtx_market::write_file(&raw_dir.join("features.tsv"), gzip, |w| {
            for (id, name) in gene_ids.iter().zip(gene_names.iter()) {
                writeln!(w, "{id}\t{name}\tGene Expression").map_err(|e| Error::io(e, raw_dir))?;
            }
            Ok(())
        })?;
        // barcodes.tsv = the manifest cell IDs (one per matrix column).
        crate::solo::mtx_market::write_file(&raw_dir.join("barcodes.tsv"), gzip, |w| {
            for cid in &self.cell_ids {
                writeln!(w, "{cid}").map_err(|e| Error::io(e, raw_dir))?;
            }
            Ok(())
        })?;

        // matrix.mtx — collect entries cell-ascending, gene-ascending.
        let mut nnz = 0usize;
        let path = raw_dir.join("matrix.mtx");
        // Pre-count nnz.
        for c in &self.counts {
            nnz += c.lock().unwrap().len();
        }
        crate::solo::mtx_market::write_file(&path, gzip, |w| {
            writeln!(w, "%%MatrixMarket matrix coordinate integer general")
                .map_err(|e| Error::io(e, &path))?;
            writeln!(w, "%").map_err(|e| Error::io(e, &path))?;
            writeln!(w, "{} {} {}", self.n_genes, self.cell_ids.len(), nnz)
                .map_err(|e| Error::io(e, &path))?;
            for (ci, cell) in self.counts.iter().enumerate() {
                let map = cell.lock().unwrap();
                let mut entries: Vec<(u32, u64)> = map.iter().map(|(&g, &c)| (g, c)).collect();
                entries.sort_unstable_by_key(|&(g, _)| g);
                for (g, c) in entries {
                    writeln!(w, "{} {} {}", g + 1, ci + 1, c).map_err(|e| Error::io(e, &path))?;
                }
            }
            Ok(())
        })?;
        Ok(nnz)
    }
}
