//! `--soloOutputFormat Zarr`: the whole solo run as one MuData store (Phase 14.x,
//! a rustar extension beyond STARsolo, which only writes MatrixMarket).
//!
//! On-disk layout (Zarr v3; `mod/<name>` groups are plain AnnData stores, the
//! MuData scaffolding around them is written here by hand — there is no MuData
//! crate):
//!
//! ```text
//! Solo.out/matrix.zarr/            encoding-type = MuData
//!   uns/summary/<Feature>/…        the Summary.csv statistics, unformatted
//!   mod/                           mod-order = [gex, sj]
//!     gex/                         cells × genes, no X
//!       layers/{Gene,GeneFull}     one per --soloFeatures gene feature
//!       layers/{spliced,unspliced,ambiguous}
//!       layers/Gene_UniqueAndMult-EM …
//!       obs (whitelist barcodes), var (gene_id + gene_name)
//!       obsm/stats_<Feature>       per-barcode reads/UMIs/genes + is_cell
//!     sj/                          cells × junctions (X only)
//! ```
//!
//! Without `--soloFeatures SJ` there is no second modality to wrap, so the
//! store is plain AnnData: `gex`'s contents (and the root `uns`) sit at
//! `matrix.zarr/` itself, with no `mod/` and no MuData attributes.
//!
//! Three things are deliberate. **`obs` is the full whitelist**, identical in every
//! matrix and both modalities — so the modalities can share one MuData `obs` and
//! no cell-calling decision is baked into the axis. **Cell calling lives in each
//! output type's `obsm` frame** as an `is_cell` column, rather than one
//! `is_cell_<Feature>` column per feature sprayed across `obs`; the MTX writer's
//! `filtered/` directory is that same boolean, materialized. And **`gex` has no
//! `X`**: no gene feature is privileged as *the* matrix, so every one of them is
//! a named layer and the caller picks (`sj` has a single matrix, so it keeps `X`).
//!
//! Unlike the MatrixMarket writer, which streams each feature through a temp
//! file, `set_layers` takes materialized arrays — so every layer is held in
//! memory at once. If that becomes the ceiling on a big run, the fix is to add
//! the layers one at a time through `AxisArraysOp::add` and drop each matrix
//! after it is written.

/// Container format for the solo count matrices (`--soloOutputFormat`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum OutputFormat {
    /// STARsolo's `raw/` + `filtered/` MatrixMarket triplet.
    Mtx,
    /// One sharded Zarr v3 MuData store.
    Zarr,
}

impl OutputFormat {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "MTX" => Some(Self::Mtx),
            "Zarr" => Some(Self::Zarr),
            _ => None,
        }
    }

    /// False when this binary was built without the backend's cargo feature.
    pub fn is_available(self) -> bool {
        match self {
            Self::Mtx => true,
            Self::Zarr => cfg!(feature = "anndata-out"),
        }
    }

    pub fn cargo_feature(self) -> &'static str {
        match self {
            Self::Mtx => "",
            Self::Zarr => "anndata-out",
        }
    }
}

#[cfg(feature = "anndata-out")]
pub use zarr::write_mudata;

#[cfg(feature = "anndata-out")]
mod zarr {
    use crate::error::Error;
    use crate::solo::SoloContext;
    use crate::solo::count::{
        CellCounts, CellStat, FeatureSummary, MultiMethod, cells_to_csr, multi_matrices,
    };
    use crate::solo::whitelist::CbWhitelist;
    use anndata::backend::{AttributeOp, Value};
    use anndata::data::{DataFrameIndex, Mapping};
    use anndata::{
        AnnData, AnnDataOp, AxisArraysOp, Backend, ElemCollectionOp, backend::GroupOp,
        data::ArrayData, data::Data,
    };
    use anndata_zarr::Zarr;
    use polars::prelude::{Column, DataFrame};
    use rustc_hash::FxHashMap as HashMap;
    use std::path::Path;

    /// mudata's on-disk format versions (`mudata.__mudataversion__` /
    /// `__anndataversion__`), which its reader checks the groups against.
    const MUDATA_VERSION: &str = "0.1.0";
    const ANNDATA_VERSION: &str = "0.1.0";

    /// Write the solo count matrices as one MuData store. No-op (with a warning)
    /// when there is no explicit whitelist, matching the MatrixMarket writer.
    pub fn write_mudata(
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

        let solo_dir = params
            .solo_out_file_names
            .first()
            .cloned()
            .unwrap_or_else(|| "Solo.out/".to_string());
        let stem = params
            .solo_out_file_names
            .get(3)
            .map_or("matrix", |m| {
                Path::new(m)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("matrix")
            })
            .to_string();
        let path = params.output_path(&format!("{solo_dir}{stem}.zarr"));
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| Error::io(e, parent))?;
        }
        // Without the sj modality there is nothing to wrap: gex is the store root.
        let mudata = ctx.sj_enabled && sj_stats.is_some();
        // The store root. `Zarr::new` wipes any previous store at this path.
        let store = Zarr::new(&path).map_err(zarr_err)?;
        if mudata {
            std::fs::create_dir_all(path.join("mod")).map_err(|e| Error::io(e, &path))?;
        }

        // The one `obs` axis every matrix and both modalities are written against.
        let obs_index = DataFrameIndex::from(
            (0..n_obs)
                .map(|i| {
                    let mut buf = Vec::with_capacity(ctx.whitelist.barcode_len());
                    ctx.whitelist.unpack_barcode_into(i as u32, &mut buf);
                    String::from_utf8_lossy(&buf).into_owned()
                })
                .collect::<Vec<String>>(),
        );

        // -- gex: every gene-indexed output type shares one cells × genes axis --
        // One (layer, obsm frame, uns summary) per output type, unzipped into the
        // three collections the AnnData setters take.
        let per_feature = ctx
            .features
            .iter()
            .zip(&ctx.recorders)
            .map(|(feature, recorder)| {
                let name = feature.dir_name().to_string();
                let cells = crate::solo::count::dedup_cells(ctx, recorder, &opts);
                let stats: Vec<CellStat> = cells
                    .iter()
                    .map(CellCounts::stat)
                    .filter(|s| s.n_umis > 0)
                    .collect();
                let matrix = cells_to_csr(&cells, n_obs, n_genes)?;
                let called = call_cells(&stats, &matrix, n_genes, params)?;

                // UniqueAndMult-<method> variants share the layer axis, so they
                // ride along as extra layers of this feature.
                let mg = recorder.multi_gene.lock().unwrap();
                let mut layers: Vec<(String, ArrayData)> = Vec::new();
                for m in multi_matrices(&matrix, &mg, &multi_methods) {
                    let (method, mat) = m?;
                    layers.push((
                        format!("{name}_UniqueAndMult-{}", method.name()),
                        mat.into(),
                    ));
                }
                drop(mg);

                let summary = crate::solo::count::feature_summary(
                    &stats,
                    detected(&matrix, n_genes),
                    &funnel,
                    crate::solo::count::feature_reads(ctx, *feature),
                );
                // The feature's own counts lead, ahead of its multimapper variants.
                layers.insert(0, (name.clone(), matrix.into()));
                Ok((
                    layers,
                    (
                        (
                            format!("stats_{name}"),
                            stats_frame(&stats, &called, n_obs)?,
                        ),
                        (name, summary_map(&summary)),
                    ),
                ))
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let (layer_groups, (obsm, summaries)): (Vec<_>, (Vec<_>, Vec<_>)) =
            per_feature.into_iter().unzip();
        let mut layers: Vec<(String, ArrayData)> = layer_groups.into_iter().flatten().collect();

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
                layers.push((
                    (*name).to_string(),
                    cells_to_csr(&per_cat, n_obs, n_genes)?.into(),
                ));
            }
        }

        // Root `uns`: on the MuData root, or on gex itself when it is the root.
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

        let mut modalities: Vec<&str> = Vec::new();
        if layers.is_empty() {
            log::warn!(
                "STARsolo: no gene-indexed features to write to {}",
                path.display()
            );
        } else {
            // All layers are output as-is with keys in Layers.
            let gex_path = if mudata {
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
            gex.set_layers(layers).map_err(zarr_err)?;
            for (name, df) in obsm {
                gex.obsm().add(&name, df).map_err(zarr_err)?;
            }
            if !mudata {
                for (k, v) in uns.drain() {
                    gex.uns().add(&k, v).map_err(zarr_err)?;
                }
            }
            gex.close().map_err(zarr_err)?;
            modalities.push("gex");
            log::info!(
                "STARsolo: wrote {} ({n_obs} barcodes × {n_genes} genes)",
                gex_path.display(),
            );
        }
        if !mudata {
            // AnnData does not tag its own root; mudata's `mod/*` get this below.
            let mut root = store.open_group("/").map_err(zarr_err)?;
            return set_attrs(&mut root, "anndata", ANNDATA_VERSION);
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

        // -- the MuData scaffolding: root uns, group attributes, mod-order --
        anndata::data::Writable::write(
            &Mapping::from(uns.into_iter().collect::<std::collections::HashMap<_, _>>()),
            &store,
            "uns",
        )
        .map_err(zarr_err)?;

        for m in &modalities {
            let mut g = store.open_group(&format!("mod/{m}")).map_err(zarr_err)?;
            set_attrs(&mut g, "anndata", ANNDATA_VERSION)?;
        }
        let mut mod_group = store.new_group("mod").map_err(zarr_err)?;
        mod_group
            .new_json_attr(
                "mod-order",
                &Value::Array(modalities.iter().map(|m| Value::from(*m)).collect()),
            )
            .map_err(zarr_err)?;
        let mut root = store.new_group("/").map_err(zarr_err)?;
        set_attrs(&mut root, "MuData", MUDATA_VERSION)?;
        // axis = 0: the modalities share `obs` (the whitelist) and concatenate `var`.
        root.new_json_attr("axis", &Value::from(0))
            .map_err(zarr_err)?;
        Ok(())
    }

    fn set_attrs(
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
        group
            .new_json_attr("encoder", &Value::from("rustar-aligner"))
            .map_err(zarr_err)?;
        group
            .new_json_attr("encoder-version", &Value::from(env!("CARGO_PKG_VERSION")))
            .map_err(zarr_err)
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

    /// `var` for the gene modality (the index is set separately from `gene_ids`).
    fn gene_frame(gene_names: &[String]) -> Result<DataFrame, Error> {
        DataFrame::new(
            gene_names.len(),
            vec![Column::new("gene_name".into(), gene_names)],
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
