//! `--soloOutputFormat MTX`: STARsolo's MatrixMarket output — `raw/` +
//! `filtered/` `{matrix.mtx, barcodes.tsv, features.tsv}` per feature, the
//! `Summary.csv` files, and `--runMode soloCellFiltering` over an existing
//! raw matrix. Counting and cell calling live in [`crate::solo::count`].

use crate::error::Error;
use crate::solo::SoloContext;
use crate::solo::count::{
    CellCounts, CellRangerSummary, CellStat, CountOptions, MappingFunnel, MultiMethod,
    called_cells, dedup_cells, distribute_multi, emptydrops_called, feature_reads, feature_summary,
    sj_cells, velocyto_cells,
};
use crate::solo::whitelist::CbWhitelist;
use rustc_hash::FxHashMap as HashMap;
use std::io::{BufRead, BufReader, Write as _};
use std::path::{Path, PathBuf};

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

/// What `build_matrix_body` returns alongside the temp matrix body.
pub struct MatrixStats {
    pub nnz: usize,
    /// One entry per barcode that received ≥1 UMI (the raw, unfiltered set).
    pub cells: Vec<CellStat>,
    /// Distinct genes with a nonzero count anywhere in the raw matrix.
    pub genes_detected: u32,
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
        let nnz = write_as_mtx_sj_matrix(
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

/// Write the SJ feature matrix (junctions are rows) in the same MatrixMarket
/// layout as the gene matrix.
fn write_as_mtx_sj_matrix(
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
