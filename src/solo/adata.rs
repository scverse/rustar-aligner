//! `--soloOutputFormat scverse`: the whole solo run as one mudata (if junctions are present) or anndata store (i.e., just gex).
//
//! On-disk layout (scverse v3; `mod/<name>` groups are plain AnnData stores, the
//! mudata scaffolding around them is written here by hand if necessary — there is no mudata
//! crate):
//!
//! ```text
//! Solo.out/matrix.zarr/            encoding-type = mudata
//!   uns/summary/<Feature>/…        the Summary.csv statistics, unformatted
//!   uns/cell_reads_not_in_passlist CellReads.stats' CBnotInPasslist dataframe (--soloCellReadStats CB)
//!   mod/                           mod-order = [gex, sj]
//!     gex/                         cells × genes, no X
//!       layers/{Gene,GeneFull}     one per --soloFeatures gene feature
//!       layers/{spliced,unspliced,ambiguous}
//!       layers/Gene_UniqueAndMult-EM …
//!       obs (whitelist barcodes), var (gene_id + gene_name, feature_types)
//!       obsm/stats_<Feature>       per-barcode reads/UMIs/genes + is_cell
//!                                  (+ the CellReads.stats counters with --soloCellReadStats CB)
//!     sj/                          cells × junctions (X only)
//! ```
//!
//! Without `--soloFeatures SJ` there is no second modality to wrap, so the
//! store is plain AnnData: `gex`'s contents (and the root `uns`) sit at
//! `matrix.zarr/` itself, with no `mod/` and no mudata attributes.
//!
//! Three things are deliberate. **`obs` is the full whitelist**, identical in every
//! matrix and both modalities — so the modalities can share one mudata `obs` and
//! no cell-calling decision is baked into the axis. **Cell calling lives in each
//! output type's `obsm` frame** as an `is_cell` column, rather than one
//! `is_cell_<Feature>` column per feature sprayed across `obs`; the MTX writer's
//! `filtered/` directory is that same boolean, materialized. And **`gex` has no
//! `X`**: no gene feature is privileged as *the* matrix, so every one of them is
//! a named layer and the caller picks (`sj` has a single matrix, so it keeps `X`).

/// Container format for the solo count matrices (`--soloOutputFormat`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum OutputFormat {
    /// STARsolo's `raw/` + `filtered/` MatrixMarket triplet.
    Mtx,
    /// One sharded scverse v3 {mu,ann}data store.
    #[allow(non_camel_case_types)]
    scverse,
}

impl OutputFormat {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "MTX" => Some(Self::Mtx),
            "scverse" => Some(Self::scverse),
            _ => None,
        }
    }

    /// False when this binary was built without the backend's cargo feature.
    pub fn is_available(self) -> bool {
        match self {
            Self::Mtx => true,
            Self::scverse => cfg!(feature = "scverse-out"),
        }
    }

    pub fn cargo_feature(self) -> &'static str {
        match self {
            Self::Mtx => "",
            Self::scverse => "scverse-out",
        }
    }
}

#[cfg(feature = "scverse-out")]
pub use zarr::write_scverse;

#[cfg(feature = "scverse-out")]
pub use zarr::is_scverse_store;

/// `--runMode soloCellFiltering <matrix.zarr> <column>`: re-call cells in a
/// scverse store without aligning anything.
///
/// The MatrixMarket mode leaves the raw matrix alone and writes the new call
/// to a path the user names. The store keeps the same contract: counts and the
/// run's own `is_cell` are untouched, and the new call is added as the boolean
/// column `<column>` of every `obsm/stats_<Feature>` frame, so the two calls
/// can be compared side by side. Naming an existing column (e.g. `is_cell`)
/// overwrites it.
pub fn run_cell_filtering(params: &crate::params::Parameters) -> anyhow::Result<()> {
    #[cfg(feature = "scverse-out")]
    {
        zarr::recall_cells(
            std::path::Path::new(&params.run_mode_in[1]),
            &params.run_mode_in[2],
            params,
        )?;
        Ok(())
    }
    #[cfg(not(feature = "scverse-out"))]
    {
        anyhow::bail!(
            "{} is a scverse store, but this binary was built without the `scverse-out` feature",
            params.run_mode_in[1]
        )
    }
}

#[cfg(feature = "scverse-out")]
mod zarr {
    use crate::error::Error;
    use crate::solo::SoloContext;
    use crate::solo::cell_reads::CellReadStats;
    use crate::solo::count::{
        CellCounts, CellStat, FeatureSummary, MultiMethod, cells_to_csr, multi_matrices,
    };
    use crate::solo::whitelist::CbWhitelist;
    use anndata::ArrayElemOp;
    use anndata::backend::{AttributeOp, Value};
    use anndata::data::{DataFrameIndex, Mapping};
    use anndata::{
        AnnData, AnnDataOp, AxisArraysOp, Backend, ElemCollectionOp, backend::GroupOp, data::Data,
    };
    use anndata_zarr::Zarr;
    use nalgebra_sparse::CsrMatrix;
    use polars::prelude::{Column, DataFrame};
    use rustc_hash::FxHashMap as HashMap;
    use std::path::Path;

    /// mudata's on-disk format versions (`mudata.__mudataversion__` /
    /// `__anndataversion__`), which its reader checks the groups against.
    /// anndata-rs does not stamp its anndata stores with this but it is technically required..
    const MUDATA_VERSION: &str = "0.1.0";
    const ANNDATA_VERSION: &str = "0.1.0";

    /// True when `path` is the root of a zarr store (a `--soloOutputFormat scverse`
    /// output) rather than a MatrixMarket `raw/` directory.
    pub fn is_scverse_store(path: &std::path::Path) -> bool {
        let group_res = Zarr::open(path);
        group_res.is_ok_and(|store| {
            store.open_group("/").is_ok_and(|g| {
                g.get_json_attr("encoding-type").is_ok_and(|encoding_type| {
                    g.get_json_attr("encoding-version")
                        .is_ok_and(|encoding_version| match (encoding_type, encoding_version) {
                            (Value::String(t), Value::String(v)) => {
                                (v == MUDATA_VERSION && t == "mudata")
                                    || (v == ANNDATA_VERSION && t == "anndata")
                            }
                            (_, _) => false,
                        })
                })
            })
        })
    }

    /// Write the solo count matrices as one mudata store. No-op (with a warning)
    /// when there is no explicit whitelist, matching the MatrixMarket writer.
    pub fn write_scverse(
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
        let n_obs = sorted.len();
        let n_genes = ctx.gene_ann.gene_ids.len();
        let opts = crate::solo::count::CountOptions::from_params(params);
        let multi_methods = MultiMethod::parse_list(&params.solo_multi_mappers);
        let funnel = crate::solo::count::MappingFunnel::collect(ctx, align_stats);

        let path = store_path(params)?;
        // Without the sj modality there is nothing to wrap: gex is the store root.
        let should_write_mudata = ctx.sj_enabled && sj_stats.is_some();
        // The store root. `Zarr::new` wipes any previous store at this path.
        let store = Zarr::new(&path).map_err(zarr_err)?;
        if should_write_mudata {
            std::fs::create_dir_all(path.join("mod")).map_err(|e| Error::io(e, &path))?;
        }

        // The one `obs` axis every matrix and both modalities are written against.
        // This matches the use of `write_barcodes` against the full whitelist in the matrix market output.
        let mut buf = Vec::with_capacity(ctx.whitelist.barcode_len());
        let obs_index = DataFrameIndex::from(
            (0..n_obs)
                .map(|i| {
                    ctx.whitelist.unpack_barcode_into(i as u32, &mut buf);
                    String::from_utf8_lossy(&buf).into_owned()
                })
                .collect::<Vec<String>>(),
        );
        drop(buf);

        // -- gex: every gene-indexed output type shares one cells × genes axis --
        // `features` is never empty (it defaults to Gene), so gex always exists.
        let gex_path = if should_write_mudata {
            path.join("mod").join("gex")
        } else {
            path.clone()
        };
        let gex = AnnData::<Zarr>::new(&gex_path).map_err(zarr_err)?;
        gex.set_var(gene_frame(&ctx.gene_ann.gene_names)?)
            .map_err(zarr_err)?;
        gex.set_var_names(DataFrameIndex::from(ctx.gene_ann.gene_ids.clone()))
            .map_err(zarr_err)?;
        gex.set_obs_names(obs_index.clone()).map_err(zarr_err)?;

        // `--soloCellReadStats CB`: the read counters are per barcode, not per
        // feature, but the MTX writer puts `CellReads.stats` in every feature's
        // directory, so they join every gene feature's stats frame.
        let cell_reads = ctx.cell_read_stats.as_ref().map(|s| {
            let s = s.lock().unwrap();
            (s.columns(n_obs), *s.not_in_passlist())
        });

        // Each layer is written as soon as it is built and dropped right after,
        // so at most a feature's matrix plus one derived matrix is in memory.
        let mut summaries: Vec<(String, Mapping)> = Vec::new();
        for (feature, recorder) in ctx.features.iter().zip(&ctx.recorders) {
            let name = feature.dir_name().to_string();
            let cells = crate::solo::count::dedup_cells(ctx, recorder, &opts);
            let stats: Vec<CellStat> = cells
                .iter()
                .map(CellCounts::stat)
                .filter(|s| s.n_umis > 0)
                .collect();
            let matrix = cells_to_csr(&cells, n_obs, n_genes)?;
            drop(cells);
            let called = call_cells(&stats, &matrix, n_genes, params)?;

            // UniqueAndMult-<method> variants share the layer axis, so they
            // ride along as extra layers of this feature.
            let mg = recorder.multi_gene.lock().unwrap();
            for m in multi_matrices(&matrix, &mg, &multi_methods) {
                let (method, mat) = m?;
                gex.layers()
                    .add(&format!("{name}_UniqueAndMult-{}", method.name()), mat)
                    .map_err(zarr_err)?;
            }
            drop(mg);

            let summary = crate::solo::count::feature_summary(
                &stats,
                detected(&matrix, n_genes),
                &funnel,
                crate::solo::count::feature_reads(ctx, *feature),
            );
            gex.layers().add(&name, matrix).map_err(zarr_err)?;
            let mut frame = stats_frame(&stats, &called, n_obs)?;
            if let Some((cols, _)) = &cell_reads {
                for (col, values) in CellReadStats::COLUMNS.iter().zip(cols) {
                    frame
                        .with_column(Column::new((*col).into(), values))
                        .map_err(|e| Error::Parameter(format!("building solo obsm frame: {e}")))?;
                }
            }
            gex.obsm()
                .add(&format!("stats_{name}"), frame)
                .map_err(zarr_err)?;
            summaries.push((name, summary_map(&summary)));
        }

        // Velocyto: three more layers on the same gene axis, no separate stats.
        if ctx.velocyto_enabled {
            let keep_ambiguous = !matches!(
                params.solo_velocyto_ambiguous.as_str(),
                "no" | "No" | "false"
            );
            let cells = crate::solo::count::velocyto_cells(
                &ctx.velocyto_records.lock().unwrap(),
                opts.method,
                opts.umi_len,
                keep_ambiguous,
            );
            let names = ["spliced", "unspliced", "ambiguous"];
            for (k, name) in names
                .iter()
                .enumerate()
                .take(if keep_ambiguous { 3 } else { 2 })
            {
                let per_cat: Vec<CellCounts> = cells.iter().map(|c| c[k].clone()).collect();
                gex.layers()
                    .add(name, cells_to_csr(&per_cat, n_obs, n_genes)?)
                    .map_err(zarr_err)?;
            }
        }

        // Root `uns`: on the mudata root, or on gex itself when it is the root.
        let mut uns: HashMap<String, Data> = HashMap::default();
        uns.insert(
            "summary".to_string(),
            Mapping::from(
                summaries
                    .into_iter()
                    .map(|(k, v)| (k, Data::from(v)))
                    .collect::<std::collections::HashMap<_, _>>(),
            )
            .into(),
        );
        if let Some(cr) = funnel.cellranger_summary() {
            uns.insert("cellranger_summary".to_string(), cellranger_map(&cr).into());
        }
        uns.insert("run_info".to_string(), run_info().into());
        // The file's `CBnotInPasslist` row: reads that never reached a barcode,
        // so they have no `obs` row to live in.
        if let Some((_, no_cb)) = &cell_reads {
            uns.insert(
                "cell_reads_not_in_passlist".to_string(),
                Mapping::from(
                    CellReadStats::COLUMNS
                        .iter()
                        .zip(no_cb)
                        .map(|(k, &v)| ((*k).to_string(), Data::from(v)))
                        .collect::<std::collections::HashMap<_, _>>(),
                )
                .into(),
            );
        }

        if !should_write_mudata {
            for (k, v) in uns.drain() {
                gex.uns().add(&k, v).map_err(zarr_err)?;
            }
        }
        gex.close().map_err(zarr_err)?;
        let mut modalities: Vec<&str> = vec!["gex"];
        log::info!(
            "STARsolo: wrote {} ({n_obs} barcodes × {n_genes} genes)",
            gex_path.display(),
        );
        if !should_write_mudata {
            // AnnData does not tag its own root; mudata's `mod/*` get this below.
            let mut root = store.open_group("/").map_err(zarr_err)?;
            return set_encoding(&mut root, "anndata", ANNDATA_VERSION);
        }

        // -- sj: junction-indexed, so its own modality (it cannot share `var`) --
        if ctx.sj_enabled
            && let Some(sjs) = sj_stats
        {
            let junctions = sjs.sorted_junctions(params);
            let row: HashMap<(u64, u64), u32> = junctions
                .iter()
                .enumerate()
                .map(|(i, (k, _))| ((k.intron_start, k.intron_end), i as u32))
                .collect();
            let cells = crate::solo::count::sj_cells(
                &ctx.sj_records.lock().unwrap(),
                &row,
                opts.method,
                opts.umi_len,
            );
            let stats: Vec<CellStat> = cells
                .iter()
                .map(CellCounts::stat)
                .filter(|s| s.n_umis > 0)
                .collect();
            let matrix = cells_to_csr(&cells, n_obs, junctions.len())?;
            let called = call_cells(&stats, &matrix, junctions.len(), params)?;

            let sj = AnnData::<Zarr>::new(path.join("mod").join("sj")).map_err(zarr_err)?;
            let (index, var) = junction_frame(&junctions, genome)?;
            sj.set_var(var).map_err(zarr_err)?;
            sj.set_var_names(DataFrameIndex::from(index))
                .map_err(zarr_err)?;
            sj.set_obs_names(obs_index.clone()).map_err(zarr_err)?;
            sj.set_x(matrix).map_err(zarr_err)?;
            sj.obsm()
                .add("stats_SJ", stats_frame(&stats, &called, n_obs)?)
                .map_err(zarr_err)?;
            sj.close().map_err(zarr_err)?;
            modalities.push("sj");
            log::info!(
                "STARsolo: wrote {}/mod/sj ({n_obs} barcodes × {} junctions)",
                path.display(),
                junctions.len(),
            );
        }

        // -- the mudata scaffolding: root uns, group attributes, mod-order --
        anndata::data::Writable::write(
            &Mapping::from(uns.into_iter().collect::<std::collections::HashMap<_, _>>()),
            &store,
            "uns",
        )
        .map_err(zarr_err)?;

        for m in &modalities {
            let mut g = store.open_group(&format!("mod/{m}")).map_err(zarr_err)?;
            set_encoding(&mut g, "anndata", ANNDATA_VERSION)?;
        }
        let mut mod_group = store.new_group("mod").map_err(zarr_err)?;
        mod_group
            .new_json_attr(
                "mod-order",
                &Value::Array(modalities.iter().map(|m| Value::from(*m)).collect()),
            )
            .map_err(zarr_err)?;
        let mut root = store.new_group("/").map_err(zarr_err)?;
        set_encoding(&mut root, "mudata", MUDATA_VERSION)?;
        // axis = 0: the modalities share `obs` (the whitelist) and concatenate `var`.
        root.new_json_attr("axis", &Value::from(0))
            .map_err(zarr_err)?;
        Ok(())
    }

    // Lives here rather than in `smartseq.rs` to reach the store helpers.
    impl crate::solo::smartseq::SmartSeqCounts {
        /// `--soloType SmartSeq` as a plain AnnData store: manifest cells × genes,
        /// read counts in `layers/Gene` (no UMIs, so `n_umis` is the read count).
        /// Every manifest cell is a cell — the MTX writer has no `filtered/` either —
        /// so `is_cell` is all true; `soloCellFiltering` can still re-call it.
        pub fn write_scverse(
            &self,
            gene_ids: &[String],
            gene_names: &[String],
            params: &crate::params::Parameters,
        ) -> Result<std::path::PathBuf, Error> {
            let counts = self;
            let n_obs = counts.cell_ids.len();
            let cells: Vec<CellCounts> = counts
                .counts
                .iter()
                .enumerate()
                .map(|(cb, m)| {
                    let mut entries: Vec<(u32, u64)> =
                        m.lock().unwrap().iter().map(|(&g, &c)| (g, c)).collect();
                    entries.sort_unstable_by_key(|&(g, _)| g);
                    CellCounts {
                        cb: cb as u32,
                        n_reads: entries.iter().map(|&(_, c)| c).sum(),
                        entries,
                    }
                })
                .collect();
            let matrix = cells_to_csr(&cells, n_obs, counts.n_genes)?;
            drop(cells);
            let stats = row_stats(&matrix);
            let all: Vec<u32> = (0..n_obs as u32).collect();

            let path = store_path(params)?;
            let store = Zarr::new(&path).map_err(zarr_err)?;
            let adata = AnnData::<Zarr>::new(&path).map_err(zarr_err)?;
            adata.set_var(gene_frame(gene_names)?).map_err(zarr_err)?;
            adata
                .set_var_names(DataFrameIndex::from(gene_ids.to_vec()))
                .map_err(zarr_err)?;
            adata
                .set_obs_names(DataFrameIndex::from(counts.cell_ids.clone()))
                .map_err(zarr_err)?;
            adata
                .obsm()
                .add("stats_Gene", stats_frame(&stats, &all, n_obs)?)
                .map_err(zarr_err)?;
            adata.layers().add("Gene", matrix).map_err(zarr_err)?;
            adata.uns().add("run_info", run_info()).map_err(zarr_err)?;
            adata.close().map_err(zarr_err)?;
            let mut root = store.open_group("/").map_err(zarr_err)?;
            set_encoding(&mut root, "anndata", ANNDATA_VERSION)?;
            Ok(path)
        }
    }

    /// `Solo.out/<matrix stem>.zarr` under the output prefix, parent created.
    fn store_path(params: &crate::params::Parameters) -> Result<std::path::PathBuf, Error> {
        let solo_dir = params
            .solo_out_file_names
            .first()
            .cloned()
            .unwrap_or_else(|| "Solo.out/".to_string());
        let stem = params.solo_out_file_names.get(3).map_or("matrix", |m| {
            Path::new(m)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("matrix")
        });
        let path = params.output_path(&format!("{solo_dir}{stem}.zarr"));
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| Error::io(e, parent))?;
        }
        Ok(path)
    }

    /// Per-row totals of a cells × features matrix, for rows with any counts.
    /// Reads and UMIs are both the row sum (the matrix keeps no read count).
    fn row_stats(matrix: &CsrMatrix<u64>) -> Vec<CellStat> {
        matrix
            .row_iter()
            .enumerate()
            .map(|(cb, row)| {
                let n_umis = row.values().iter().sum();
                CellStat {
                    cb: cb as u32,
                    n_reads: n_umis,
                    n_umis,
                    n_genes: row.nnz() as u32,
                }
            })
            .filter(|s| s.n_umis > 0)
            .collect()
    }

    /// Columns of the stats frames that hold counts, not a cell call; the
    /// re-call must not clobber them.
    const STATS_COUNT_COLUMNS: [&str; 3] = ["n_reads", "n_umis", "n_features"];

    /// See [`super::run_cell_filtering`]. Each `stats_<Feature>` frame in obsm is
    /// re-called against its own matrix, one of the `<Feature>` layer or `X` for the
    /// single-matrix `sj` modality.
    pub(super) fn recall_cells(
        store: &Path,
        column: &str,
        params: &crate::params::Parameters,
    ) -> Result<(), Error> {
        if STATS_COUNT_COLUMNS.contains(&column) || CellReadStats::COLUMNS.contains(&column) {
            return Err(Error::Parameter(format!(
                "soloCellFiltering: `{column}` holds counts; choose another column name for the cell call"
            )));
        }
        // MuData keeps one AnnData per modality under `mod/`; plain AnnData is the root.
        let roots: Vec<std::path::PathBuf> = if store.join("mod").is_dir() {
            ["gex", "sj"]
                .iter()
                .map(|m| store.join("mod").join(m))
                .filter(|p| p.is_dir())
                .collect()
        } else {
            vec![store.to_path_buf()]
        };

        let mut n_frames = 0;
        for root in roots {
            let adata =
                AnnData::<Zarr>::open(Zarr::open_rw(&root).map_err(zarr_err)?).map_err(zarr_err)?;
            let keys: Vec<String> = adata
                .obsm()
                .keys()
                .into_iter()
                .filter(|k| k.starts_with("stats_"))
                .collect();
            for key in keys {
                let feature = &key["stats_".len()..];
                let matrix: CsrMatrix<u64> =
                    if adata.layers().keys().iter().any(|l| l == feature) {
                        adata.layers().get_item(feature)
                    } else {
                        adata.x().get()
                    }
                    .map_err(zarr_err)?
                    .ok_or_else(|| {
                        Error::Parameter(format!("{}: no matrix for obsm/{key}", root.display()))
                    })?;

                // Cell calling only reads per-barcode UMI totals.
                let stats = row_stats(&matrix);
                let called = call_cells(&stats, &matrix, matrix.ncols(), params)?;

                let mut is_cell = vec![false; matrix.nrows()];
                for &cb in &called {
                    is_cell[cb as usize] = true;
                }
                let mut frame: DataFrame = adata
                    .obsm()
                    .get_item(&key)
                    .map_err(zarr_err)?
                    .ok_or_else(|| {
                        Error::Parameter(format!("{}: obsm/{key} is empty", root.display()))
                    })?;
                frame
                    .with_column(Column::new(column.into(), is_cell))
                    .map_err(|e| Error::Parameter(format!("updating obsm/{key}: {e}")))?;
                adata.obsm().add(&key, frame).map_err(zarr_err)?;
                log::info!(
                    "soloCellFiltering: {} cells called -> {}/obsm/{key}/{column}",
                    called.len(),
                    root.display()
                );
                n_frames += 1;
            }
            adata.close().map_err(zarr_err)?;
        }
        if n_frames == 0 {
            return Err(Error::Parameter(format!(
                "{}: no obsm/stats_<Feature> frames to re-call",
                store.display()
            )));
        }
        Ok(())
    }

    fn set_encoding(
        group: &mut <Zarr as Backend>::Group,
        encoding: &str,
        version: &str,
    ) -> Result<(), Error> {
        group
            .new_json_attr("encoding-type", &Value::from(encoding))
            .map_err(zarr_err)?;
        group
            .new_json_attr("encoding-version", &Value::from(version))
            .map_err(zarr_err)?;
        Ok(())
    }

    /// The whitelist indices `--soloCellFilter` calls cells for this matrix —
    /// the `filtered/` directory of the MatrixMarket writer, as a list. Empty
    /// for `--soloCellFilter None`.
    fn call_cells(
        stats: &[CellStat],
        matrix: &nalgebra_sparse::CsrMatrix<u64>,
        n_features: usize,
        params: &crate::params::Parameters,
    ) -> Result<Vec<u32>, Error> {
        if params
            .solo_cell_filter
            .first()
            .is_some_and(|m| m == "EmptyDrops_CR")
        {
            // Rows are whitelist barcodes, so the row index *is* the barcode.
            let triplets = matrix
                .triplet_iter()
                .map(|(cb, gene, &v)| Ok((gene as u32, cb as u32, v as u32)));
            crate::solo::count::emptydrops_called(
                stats,
                triplets,
                n_features,
                &params.solo_cell_filter,
            )
        } else {
            Ok(
                crate::solo::count::called_cells(stats, &params.solo_cell_filter)
                    .unwrap_or_default(),
            )
        }
    }

    /// Genes with a nonzero count anywhere in the matrix.
    fn detected(matrix: &nalgebra_sparse::CsrMatrix<u64>, n_cols: usize) -> u32 {
        let mut seen = vec![false; n_cols];
        for &c in matrix.col_indices() {
            seen[c] = true;
        }
        seen.iter().filter(|&&s| s).count() as u32
    }

    /// One output type's per-barcode statistics, full whitelist height: the
    /// counting stats plus this output type's own cell call.
    fn stats_frame(stats: &[CellStat], called: &[u32], n_obs: usize) -> Result<DataFrame, Error> {
        let mut n_reads = vec![0u64; n_obs];
        let mut n_umis = vec![0u64; n_obs];
        let mut n_features = vec![0u32; n_obs];
        let mut is_cell = vec![false; n_obs];
        for s in stats {
            n_reads[s.cb as usize] = s.n_reads;
            n_umis[s.cb as usize] = s.n_umis;
            n_features[s.cb as usize] = s.n_genes;
        }
        for &cb in called {
            is_cell[cb as usize] = true;
        }
        DataFrame::new(
            n_obs,
            vec![
                Column::new("n_reads".into(), n_reads),
                Column::new("n_umis".into(), n_umis),
                Column::new("n_features".into(), n_features),
                Column::new("is_cell".into(), is_cell),
            ],
        )
        .map_err(|e| Error::Parameter(format!("building solo obsm frame: {e}")))
    }

    /// `var` for the gene modality (the index is set separately from `gene_ids`):
    /// `features.tsv` columns 2 and 3, named as scanpy's `read_10x_mtx` names them.
    fn gene_frame(gene_names: &[String]) -> Result<DataFrame, Error> {
        DataFrame::new(
            gene_names.len(),
            vec![
                Column::new("gene_name".into(), gene_names),
                Column::new(
                    "feature_types".into(),
                    vec!["Gene Expression"; gene_names.len()],
                ),
            ],
        )
        .map_err(|e| Error::Parameter(format!("building solo var frame: {e}")))
    }

    /// `(var_names, var)` for the SJ modality: `chr:start-end:strand` ids plus the
    /// per-junction `SJ.out.tab` columns, which STARsolo gives us for free here.
    fn junction_frame(
        junctions: &[(
            crate::junction::sj_output::SjKey,
            crate::junction::sj_output::SjRowCounts,
        )],
        genome: &crate::genome::Genome,
    ) -> Result<(Vec<String>, DataFrame), Error> {
        let strand_char = |s: u8| match s {
            1 => '+',
            2 => '-',
            _ => '.',
        };
        let mut index = Vec::with_capacity(junctions.len());
        let mut motif = Vec::with_capacity(junctions.len());
        let mut annotated = Vec::with_capacity(junctions.len());
        let mut unique = Vec::with_capacity(junctions.len());
        let mut multi = Vec::with_capacity(junctions.len());
        let mut overhang = Vec::with_capacity(junctions.len());
        for (key, counts) in junctions {
            let (chr, start, end) = crate::junction::SpliceJunctionStats::locus(key, genome)?;
            index.push(format!("{chr}:{start}-{end}:{}", strand_char(key.strand)));
            motif.push(key.motif);
            annotated.push(counts.annotated);
            unique.push(counts.unique);
            multi.push(counts.multi);
            overhang.push(counts.max_overhang);
        }
        let var = DataFrame::new(
            junctions.len(),
            vec![
                Column::new("motif".into(), motif),
                Column::new("annotated".into(), annotated),
                Column::new("n_reads_unique".into(), unique),
                Column::new("n_reads_multi".into(), multi),
                Column::new("max_overhang".into(), overhang),
            ],
        )
        .map_err(|e| Error::Parameter(format!("building solo SJ var frame: {e}")))?;
        Ok((index, var))
    }

    /// The `Summary.csv` statistics as an `uns` mapping — the same numbers, but
    /// unformatted (fractions stay fractions).
    fn summary_map(s: &FeatureSummary) -> Mapping {
        map([
            ("n_reads", Data::from(s.n_reads)),
            ("frac_valid_barcodes", s.frac_valid_barcodes.into()),
            ("saturation", s.saturation.into()),
            (
                "frac_mapped_genome_unique_multi",
                s.frac_mapped_genome_unique_multi.into(),
            ),
            (
                "frac_mapped_genome_unique",
                s.frac_mapped_genome_unique.into(),
            ),
            (
                "frac_mapped_feature_unique",
                s.frac_mapped_feature_unique.into(),
            ),
            ("n_cells", (s.n_cells as u64).into()),
            ("reads_in_cells", s.reads_in_cells.into()),
            ("frac_reads_in_cells", s.frac_reads_in_cells.into()),
            ("mean_reads_per_cell", s.mean_reads_per_cell.into()),
            ("median_reads_per_cell", s.median_reads_per_cell.into()),
            ("umis_in_cells", s.umis_in_cells.into()),
            ("mean_umis_per_cell", s.mean_umis_per_cell.into()),
            ("median_umis_per_cell", s.median_umis_per_cell.into()),
            ("mean_features_per_cell", s.mean_features_per_cell.into()),
            (
                "median_features_per_cell",
                s.median_features_per_cell.into(),
            ),
            ("features_detected", u64::from(s.features_detected).into()),
        ])
    }

    fn cellranger_map(s: &crate::solo::count::CellRangerSummary) -> Mapping {
        map([
            ("n_reads", Data::from(s.n_reads)),
            (
                "frac_mapped_genome_unique",
                s.frac_mapped_genome_unique.into(),
            ),
            ("frac_exonic", s.frac_exonic.into()),
            ("frac_intronic", s.frac_intronic.into()),
            ("frac_intergenic", s.frac_intergenic.into()),
            ("frac_antisense", s.frac_antisense.into()),
        ])
    }

    /// What produced the store: version + the exact command line.
    fn run_info() -> Mapping {
        map([
            ("aligner", Data::from("rustar-aligner".to_string())),
            ("version", env!("CARGO_PKG_VERSION").to_string().into()),
            (
                "command_line",
                shlex::try_join(
                    std::env::args()
                        .collect::<Vec<_>>()
                        .iter()
                        .map(String::as_str),
                )
                .unwrap_or_default()
                .into(),
            ),
        ])
    }

    fn map<const N: usize>(entries: [(&str, Data); N]) -> Mapping {
        Mapping::from(
            entries
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect::<std::collections::HashMap<_, _>>(),
        )
    }

    /// anndata/zarrs report `anyhow::Error`; the solo writers report [`Error`].
    #[allow(clippy::needless_pass_by_value)] // used as `.map_err(zarr_err)`
    fn zarr_err(e: anyhow::Error) -> Error {
        Error::Parameter(format!("AnnData output: {e:#}"))
    }
}
