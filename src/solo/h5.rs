//! `--soloOutH5 yes`: CellRanger-style HDF5 count matrices (#270, a rustar
//! extension beyond STARsolo, which only writes MatrixMarket).
//!
//! Each `Gene`/`GeneFull` feature directory gains
//! `raw_feature_bc_matrix.h5` (and `filtered_feature_bc_matrix.h5` when a
//! `filtered/` matrix exists), next to the `raw/` and `filtered/` directories.
//! The files are built from the MatrixMarket triplet the run has just written,
//! so they hold exactly the same counts, barcodes and features, and the solo
//! counting code is left untouched. The layout is CellRanger v3's:
//!
//! ```text
//! /                        filetype = "matrix", version = 2, software_version
//!   matrix/
//!     barcodes             string [n_barcodes]
//!     data                 int32  [nnz]            (chunked, deflate)
//!     indices              int64  [nnz]  0-based feature row of each value
//!     indptr               int64  [n_barcodes + 1] (CSC column pointers)
//!     shape                int32  [2] = [n_features, n_barcodes]
//!     features/
//!       _all_tag_keys      ["genome"]
//!       id, name, feature_type, genome   string [n_features]
//! ```
//!
//! The writer is the pure-Rust `hdf5-pure` crate, compiled only with the
//! `hdf5-out` cargo feature, so the default build carries no HDF5 dependency.
//! `SJ`, `Velocyto` and the `UniqueAndMult-*.mtx` matrices are not converted:
//! CellRanger's `.h5` has no slot for junction rows, for several matrices over
//! one axis, or for real-valued counts.

use crate::error::Error;
use crate::params::Parameters;

/// The cargo feature that compiles the HDF5 writer in.
pub const CARGO_FEATURE: &str = "hdf5-out";

/// False when this binary was built without the `hdf5-out` cargo feature.
pub fn is_available() -> bool {
    cfg!(feature = "hdf5-out")
}

/// File names of one MatrixMarket triplet (before any `.gz` suffix).
#[derive(Debug, Clone)]
pub struct MtxNames {
    pub features: String,
    pub barcodes: String,
    pub matrix: String,
}

impl Default for MtxNames {
    fn default() -> Self {
        Self {
            features: "features.tsv".to_string(),
            barcodes: "barcodes.tsv".to_string(),
            matrix: "matrix.mtx".to_string(),
        }
    }
}

impl MtxNames {
    /// The names `--soloOutFileNames` gives the CB/UMI solo matrices.
    pub fn from_params(params: &Parameters) -> Self {
        let d = Self::default();
        let name = |i: usize, fallback: String| {
            params
                .solo_out_file_names
                .get(i)
                .cloned()
                .unwrap_or(fallback)
        };
        Self {
            features: name(1, d.features),
            barcodes: name(2, d.barcodes),
            matrix: name(3, d.matrix),
        }
    }
}

/// Write the `.h5` matrices for each feature directory `<soloOutFileNames[0]><dir>/`
/// in `feature_dirs`, from the MatrixMarket files already written there.
pub fn write_solo_h5(
    params: &Parameters,
    feature_dirs: &[&str],
    names: &MtxNames,
) -> Result<(), Error> {
    #[cfg(feature = "hdf5-out")]
    {
        imp::write_solo_h5(params, feature_dirs, names)
    }
    #[cfg(not(feature = "hdf5-out"))]
    {
        // `Parameters::validate` refuses `--soloOutH5 yes` in this build.
        let _ = (params, feature_dirs, names);
        Err(Error::Parameter(format!(
            "--soloOutH5 yes needs the `{CARGO_FEATURE}` cargo feature, which this binary was built without"
        )))
    }
}

/// `molecule_info.h5` for a CellRanger-layout run (see `--soloOutLayout`).
/// A no-op unless the run kept its molecules, which needs the `hdf5-out` build.
pub fn write_molecule_info(
    sctx: &crate::solo::SoloContext,
    params: &Parameters,
) -> Result<(), Error> {
    #[cfg(feature = "hdf5-out")]
    {
        imp::write_molecule_info_for_run(sctx, params)
    }
    #[cfg(not(feature = "hdf5-out"))]
    {
        let _ = (sctx, params);
        Ok(())
    }
}

#[cfg(feature = "hdf5-out")]
pub use imp::{CrMeta, CscMatrix, MoleculeInput, read_mtx_dir, write_10x_h5, write_10x_h5_cr};

#[cfg(feature = "hdf5-out")]
mod imp {
    use super::MtxNames;
    use crate::error::Error;
    use crate::params::Parameters;
    use hdf5_pure::{AttrValue, DatasetBuilder, FileBuilder};
    use std::io::BufRead;
    use std::path::{Path, PathBuf};

    /// Chunk length of the `data`/`indices` datasets (CellRanger's own value).
    const CHUNK_LEN: u64 = 80_000;
    /// Deflate level of the chunked datasets (h5py's `compression="gzip"` default).
    const DEFLATE_LEVEL: u32 = 4;

    /// A feature-by-barcode count matrix in compressed sparse column form, plus
    /// the axis labels, i.e. one MatrixMarket triplet held in memory.
    #[derive(Debug, Default, Clone, PartialEq, Eq)]
    pub struct CscMatrix {
        pub n_features: usize,
        pub n_barcodes: usize,
        /// Per-barcode column pointers into `indices`/`data` (`n_barcodes + 1`).
        pub indptr: Vec<i64>,
        /// 0-based feature row of each stored value.
        pub indices: Vec<i64>,
        pub data: Vec<i32>,
        pub barcodes: Vec<String>,
        pub feature_ids: Vec<String>,
        pub feature_names: Vec<String>,
        pub feature_types: Vec<String>,
    }

    pub(super) fn write_solo_h5(
        params: &Parameters,
        feature_dirs: &[&str],
        names: &MtxNames,
    ) -> Result<(), Error> {
        let solo_dir = params
            .solo_out_file_names
            .first()
            .map_or("Solo.out/", String::as_str);
        let gzip = matches!(params.solo_out_gzip.as_str(), "yes" | "Yes" | "true");
        let cr = params.solo_out_layout == "CellRanger";
        // Same directory rules as the MatrixMarket writer: CellRanger has no
        // per-feature directory when there is exactly one feature.
        let per_feature_dir = !cr || feature_dirs.len() > 1;
        let meta = cr.then(|| CrMeta::from_params(params));
        for dir in feature_dirs {
            let feature_dir = if per_feature_dir {
                params.output_path(&format!("{solo_dir}{dir}/"))
            } else {
                params.output_path(solo_dir)
            };
            let (raw_sub, filt_sub) = if cr {
                ("raw_feature_bc_matrix", "filtered_feature_bc_matrix")
            } else {
                ("raw", "filtered")
            };
            for (sub, h5_name, is_raw) in [
                (raw_sub, "raw_feature_bc_matrix.h5", true),
                (filt_sub, "filtered_feature_bc_matrix.h5", false),
            ] {
                let mtx_dir = feature_dir.join(sub);
                // `filtered/` exists only when cells were called.
                if !with_gz(&mtx_dir.join(&names.matrix), gzip).exists() {
                    continue;
                }
                let m = read_mtx_dir(&mtx_dir, names, gzip)?;
                let out = feature_dir.join(h5_name);
                match &meta {
                    Some(meta) => write_10x_h5_cr(&m, &out, meta, is_raw)?,
                    None => write_10x_h5(&m, &out)?,
                }
                log::info!(
                    "STARsolo: wrote {dir}/{h5_name} ({} features × {} barcodes, {} entries)",
                    m.n_features,
                    m.n_barcodes,
                    m.data.len(),
                );
            }
        }
        Ok(())
    }

    /// The run-level facts CellRanger records in its `.h5` files.
    #[derive(Debug, Clone)]
    pub struct CrMeta {
        pub chemistry_description: String,
        pub library_id: String,
        /// Reference genome name (`reference.json` `genomes[0]`), or empty.
        pub genome: String,
        pub software_version: String,
        pub gem_group: i64,
        /// `reference.json` text, for the hashes `molecule_info.h5` records.
        pub reference_json: String,
        pub whitelist_name: String,
        pub umi_len: usize,
        pub cb_len: usize,
    }

    impl CrMeta {
        pub fn from_params(params: &Parameters) -> Self {
            let reference_json = params
                .genome_dir
                .parent()
                .map(|p| p.join("reference.json"))
                .and_then(|p| std::fs::read_to_string(p).ok())
                .unwrap_or_default();
            let genome = json_array_first(&reference_json, "genomes").unwrap_or_default();
            let umi_len = params.solo_umi_len as usize;
            let chemistry_description = match (umi_len, params.solo_barcode_mate) {
                (12, 0) => "Single Cell 3' v3 (polyA)",
                (10, 0) => "Single Cell 3' v2",
                (_, 0) => "Single Cell 3' (custom)",
                (10, _) => "Single Cell 5' v2",
                _ => "Single Cell 5' (custom)",
            }
            .to_string();
            let whitelist_name = params
                .solo_cb_whitelist
                .first()
                .map(|w| {
                    let stem = std::path::Path::new(w)
                        .file_name()
                        .map(|f| f.to_string_lossy().into_owned())
                        .unwrap_or_default();
                    stem.trim_end_matches(".gz")
                        .trim_end_matches(".txt")
                        .to_string()
                })
                .unwrap_or_default();
            Self {
                chemistry_description,
                library_id: params.solo_out_sample_id.clone(),
                genome,
                software_version: format!("rustar-aligner {}", env!("CARGO_PKG_VERSION")),
                gem_group: 1,
                reference_json,
                whitelist_name,
                umi_len,
                cb_len: params.solo_cb_len as usize,
            }
        }
    }

    /// `"key": "value"` in flat JSON text.
    fn json_str(text: &str, key: &str) -> Option<String> {
        let at = text.find(&format!("\"{key}\""))?;
        let rest = &text[at + key.len() + 2..];
        let rest = rest.trim_start().strip_prefix(':')?.trim_start();
        let rest = rest.strip_prefix('"')?;
        Some(rest[..rest.find('"')?].to_string())
    }

    /// First string of `"key": ["a", ...]` in flat JSON text.
    fn json_array_first(text: &str, key: &str) -> Option<String> {
        let at = text.find(&format!("\"{key}\""))?;
        let rest = &text[at + key.len() + 2..];
        let rest = rest.trim_start().strip_prefix(':')?.trim_start();
        let rest = rest.strip_prefix('[')?.trim_start().strip_prefix('"')?;
        Some(rest[..rest.find('"')?].to_string())
    }

    /// `path` with `.gz` appended when the solo output is gzipped.
    fn with_gz(path: &Path, gzip: bool) -> PathBuf {
        if gzip {
            let mut s = path.as_os_str().to_owned();
            s.push(".gz");
            PathBuf::from(s)
        } else {
            path.to_path_buf()
        }
    }

    fn open(path: &Path, gzip: bool) -> Result<Box<dyn BufRead>, Error> {
        let file = std::fs::File::open(path).map_err(|e| Error::io(e, path))?;
        Ok(if gzip {
            Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(
                file,
            )))
        } else {
            Box::new(std::io::BufReader::new(file))
        })
    }

    fn bad(path: &Path, msg: impl std::fmt::Display) -> Error {
        Error::io(
            std::io::Error::new(std::io::ErrorKind::InvalidData, msg.to_string()),
            path,
        )
    }

    /// Read `dir/{features,barcodes,matrix}` (each with `.gz` when `gzip`) into
    /// a [`CscMatrix`]. The matrix must be an integer `coordinate` MatrixMarket
    /// file with features as rows and barcodes as columns, as STARsolo writes it.
    pub fn read_mtx_dir(dir: &Path, names: &MtxNames, gzip: bool) -> Result<CscMatrix, Error> {
        let features_path = with_gz(&dir.join(&names.features), gzip);
        let barcodes_path = with_gz(&dir.join(&names.barcodes), gzip);
        let matrix_path = with_gz(&dir.join(&names.matrix), gzip);

        let mut m = CscMatrix::default();
        for line in open(&features_path, gzip)?.lines() {
            let line = line.map_err(|e| Error::io(e, &features_path))?;
            if line.is_empty() {
                continue;
            }
            let mut cols = line.split('\t');
            let id = cols.next().unwrap_or_default().to_string();
            let name = cols.next().map_or_else(|| id.clone(), str::to_string);
            let kind = cols.next().unwrap_or("Gene Expression").to_string();
            m.feature_ids.push(id);
            m.feature_names.push(name);
            m.feature_types.push(kind);
        }
        for line in open(&barcodes_path, gzip)?.lines() {
            let line = line.map_err(|e| Error::io(e, &barcodes_path))?;
            if let Some(bc) = line.split('\t').next().filter(|s| !s.is_empty()) {
                m.barcodes.push(bc.to_string());
            }
        }

        let mut lines = open(&matrix_path, gzip)?.lines();
        let mut next_line = || -> Result<Option<String>, Error> {
            lines
                .next()
                .transpose()
                .map_err(|e| Error::io(e, &matrix_path))
        };
        let banner = next_line()?.unwrap_or_default();
        let banner_lc = banner.to_ascii_lowercase();
        if !banner_lc.starts_with("%%matrixmarket matrix coordinate integer") {
            return Err(bad(
                &matrix_path,
                format!("expected an integer coordinate MatrixMarket banner, got '{banner}'"),
            ));
        }
        let dims = loop {
            match next_line()? {
                Some(l) if l.starts_with('%') || l.trim().is_empty() => {}
                Some(l) => break l,
                None => return Err(bad(&matrix_path, "missing the dimensions line")),
            }
        };
        let dims: Vec<usize> = dims
            .split_whitespace()
            .map(str::parse)
            .collect::<Result<_, _>>()
            .map_err(|e| bad(&matrix_path, format!("dimensions line: {e}")))?;
        let [n_rows, n_cols, nnz] = dims[..] else {
            return Err(bad(&matrix_path, "dimensions line needs rows, cols, nnz"));
        };
        if n_rows != m.feature_ids.len() || n_cols != m.barcodes.len() {
            return Err(bad(
                &matrix_path,
                format!(
                    "matrix is {n_rows} × {n_cols} but there are {} features and {} barcodes",
                    m.feature_ids.len(),
                    m.barcodes.len()
                ),
            ));
        }

        // (column, row, value), 0-based; sorted into CSC order below.
        let mut entries: Vec<(u32, u32, i32)> = Vec::with_capacity(nnz);
        while let Some(l) = next_line()? {
            if l.trim().is_empty() {
                continue;
            }
            let mut it = l.split_whitespace();
            let (Some(r), Some(c), Some(v)) = (it.next(), it.next(), it.next()) else {
                return Err(bad(&matrix_path, format!("malformed entry '{l}'")));
            };
            let parse = |s: &str| {
                s.parse::<u64>()
                    .map_err(|e| bad(&matrix_path, format!("entry '{l}': {e}")))
            };
            let (r, c, v) = (parse(r)?, parse(c)?, parse(v)?);
            if r == 0 || r > n_rows as u64 || c == 0 || c > n_cols as u64 {
                return Err(bad(&matrix_path, format!("entry '{l}' is out of bounds")));
            }
            let v = i32::try_from(v)
                .map_err(|_| bad(&matrix_path, format!("count in '{l}' exceeds int32")))?;
            entries.push(((c - 1) as u32, (r - 1) as u32, v));
        }
        if entries.len() != nnz {
            return Err(bad(
                &matrix_path,
                format!("header declares {nnz} entries, found {}", entries.len()),
            ));
        }
        entries.sort_unstable_by_key(|&(c, r, _)| (c, r));

        m.n_features = n_rows;
        m.n_barcodes = n_cols;
        m.indptr = vec![0i64; n_cols + 1];
        for &(c, _, _) in &entries {
            m.indptr[c as usize + 1] += 1;
        }
        for i in 0..n_cols {
            m.indptr[i + 1] += m.indptr[i];
        }
        m.indices = entries.iter().map(|&(_, r, _)| i64::from(r)).collect();
        m.data = entries.iter().map(|&(_, _, v)| v).collect();
        Ok(m)
    }

    fn h5_err(path: &Path, e: impl std::fmt::Display) -> Error {
        Error::io(std::io::Error::other(format!("HDF5 write: {e}")), path)
    }

    /// A fixed-width string dataset, as CellRanger writes them: ASCII when every
    /// value is ASCII, otherwise UTF-8.
    fn strings(ds: &mut DatasetBuilder, values: &[String], path: &Path) -> Result<(), Error> {
        let refs: Vec<&str> = values.iter().map(String::as_str).collect();
        if refs.iter().all(|s| s.is_ascii()) {
            ds.with_ascii_strings(&refs).map_err(|e| h5_err(path, e))?;
        } else {
            ds.with_strings(&refs).map_err(|e| h5_err(path, e))?;
        }
        Ok(())
    }

    /// Chunk and deflate a 1-D dataset of `len` elements (an empty one stays
    /// contiguous: a zero-length chunk is not a valid HDF5 layout).
    fn compress(ds: &mut DatasetBuilder, len: usize) {
        if len > 0 {
            ds.with_chunks(&[CHUNK_LEN.min(len as u64)])
                .with_deflate(DEFLATE_LEVEL);
        }
    }

    /// Write `m` to `path` in the CellRanger v3 `feature_bc_matrix.h5` layout.
    pub fn write_10x_h5(m: &CscMatrix, path: &Path) -> Result<(), Error> {
        write_matrix(m, path, None)
    }

    /// Write `m` the way `cellranger count` does: its root attributes
    /// (chemistry, library ids, gem groups), its genome column, its string
    /// widths and its chunking. `raw` selects the raw-matrix flavour, whose
    /// string datasets CellRanger pads to 256 bytes.
    pub fn write_10x_h5_cr(
        m: &CscMatrix,
        path: &Path,
        meta: &CrMeta,
        raw: bool,
    ) -> Result<(), Error> {
        write_matrix(m, path, Some((meta, raw)))
    }

    /// Fixed-width ASCII strings of an explicit `width`.
    fn strings_w(
        ds: &mut DatasetBuilder,
        values: &[String],
        width: u32,
        path: &Path,
    ) -> Result<(), Error> {
        let refs: Vec<&str> = values.iter().map(String::as_str).collect();
        ds.with_ascii_strings_sized(&refs, width)
            .map_err(|e| h5_err(path, e))?;
        Ok(())
    }

    fn max_len(values: &[String]) -> u32 {
        values.iter().map(String::len).max().unwrap_or(1).max(1) as u32
    }

    /// Chunk, shuffle and deflate a 1-D dataset the way h5py's
    /// `compression="gzip"` does for CellRanger. `unlimited` makes the first
    /// dimension resizable, as CellRanger's growable columns are.
    fn pack(
        ds: &mut DatasetBuilder,
        len: usize,
        chunk: u64,
        shuffle: bool,
        level: u32,
        unlimited: bool,
    ) {
        if len == 0 {
            return;
        }
        ds.with_chunks(&[chunk.min(len as u64)]);
        if shuffle {
            ds.with_shuffle();
        }
        ds.with_deflate(level);
        if unlimited {
            ds.with_maxshape(&[hdf5_pure::MaxExtent::Unlimited]);
        }
    }

    fn write_matrix(m: &CscMatrix, path: &Path, cr: Option<(&CrMeta, bool)>) -> Result<(), Error> {
        let mut fb = FileBuilder::new();
        if let Some((meta, _)) = cr {
            fb.set_attr(
                "chemistry_description",
                AttrValue::VarLenString(meta.chemistry_description.clone()),
            );
            fb.set_attr("filetype", AttrValue::VarLenString("matrix".to_string()));
            fb.set_attr(
                "library_ids",
                AttrValue::ascii_string_array_sized(vec![meta.library_id.clone()], 256)
                    .map_err(|e| h5_err(path, e))?,
            );
            fb.set_attr(
                "original_gem_groups",
                AttrValue::I64Array(vec![meta.gem_group]),
            );
            fb.set_attr(
                "software_version",
                AttrValue::VarLenString(meta.software_version.clone()),
            );
            fb.set_attr("version", AttrValue::I64(2));
        } else {
            fb.set_attr("filetype", AttrValue::AsciiString("matrix".to_string()));
            fb.set_attr("version", AttrValue::I64(2));
            fb.set_attr(
                "software_version",
                AttrValue::AsciiString(format!("rustar-aligner {}", env!("CARGO_PKG_VERSION"))),
            );
        }
        let raw = cr.is_some_and(|(_, raw)| raw);

        // CellRanger's raw matrix comes from its Rust writer (gzip 1, shuffle,
        // 64Ki chunks, every dataset compressed); the filtered one from h5py
        // (gzip 4, shuffle on the numeric columns, 80 000-element chunks).
        let (chunk, level) = if raw {
            (65_536, 1)
        } else {
            (CHUNK_LEN, DEFLATE_LEVEL)
        };
        let mut matrix = fb.create_group("matrix");
        let ds = matrix.create_dataset("barcodes");
        if cr.is_some() {
            let w = max_len(&m.barcodes);
            strings_w(ds, &m.barcodes, if raw { w.max(45) } else { w }, path)?;
            if raw {
                pack(ds, m.barcodes.len(), chunk, true, level, true);
            } else {
                pack(
                    ds,
                    m.barcodes.len(),
                    (m.barcodes.len() as u64).div_ceil(2),
                    false,
                    level,
                    true,
                );
            }
        } else {
            strings(ds, &m.barcodes, path)?;
        }
        let ds = matrix.create_dataset("data");
        ds.with_i32_data(&m.data);
        if cr.is_some() {
            pack(ds, m.data.len(), chunk, true, level, true);
        } else {
            compress(ds, m.data.len());
        }
        let ds = matrix.create_dataset("indices");
        ds.with_i64_data(&m.indices);
        if cr.is_some() {
            pack(ds, m.indices.len(), chunk, true, level, true);
        } else {
            compress(ds, m.indices.len());
        }
        let ds = matrix.create_dataset("indptr");
        ds.with_i64_data(&m.indptr);
        if cr.is_some() {
            pack(
                ds,
                m.indptr.len(),
                m.indptr.len() as u64,
                true,
                level,
                false,
            );
        }
        let shape = [
            i32::try_from(m.n_features).map_err(|e| h5_err(path, e))?,
            i32::try_from(m.n_barcodes).map_err(|e| h5_err(path, e))?,
        ];
        let ds = matrix.create_dataset("shape");
        ds.with_i32_data(&shape);
        if cr.is_some() {
            pack(ds, shape.len(), shape.len() as u64, true, level, false);
        }

        let mut features = matrix.create_group("features");
        let genome = cr.map_or("", |(meta, _)| meta.genome.as_str());
        if cr.is_some() {
            let one =
                |ds: &mut DatasetBuilder, values: &[String], chunked: bool| -> Result<(), Error> {
                    if raw {
                        strings_w(ds, values, 256, path)?;
                        if chunked {
                            pack(ds, values.len(), values.len() as u64, true, 1, false);
                        }
                    } else {
                        strings_w(ds, values, max_len(values), path)?;
                        if chunked {
                            pack(
                                ds,
                                values.len(),
                                (values.len() as u64).div_ceil(32).max(1),
                                false,
                                level,
                                false,
                            );
                        }
                    }
                    Ok(())
                };
            one(
                features.create_dataset("_all_tag_keys"),
                &["genome".to_string()],
                false,
            )?;
            one(
                features.create_dataset("feature_type"),
                &m.feature_types,
                true,
            )?;
            one(
                features.create_dataset("genome"),
                &vec![genome.to_string(); m.n_features],
                true,
            )?;
            one(features.create_dataset("id"), &m.feature_ids, true)?;
            one(features.create_dataset("name"), &m.feature_names, true)?;
        } else {
            strings(
                features.create_dataset("_all_tag_keys"),
                &["genome".to_string()],
                path,
            )?;
            strings(features.create_dataset("id"), &m.feature_ids, path)?;
            strings(features.create_dataset("name"), &m.feature_names, path)?;
            strings(
                features.create_dataset("feature_type"),
                &m.feature_types,
                path,
            )?;
            // The reference name is not known here; CellRanger's readers accept an
            // empty genome and only use it to split multi-genome references.
            strings(
                features.create_dataset("genome"),
                &vec![String::new(); m.n_features],
                path,
            )?;
        }
        matrix.add_group(features.finish());
        fb.add_group(matrix.finish());

        fb.write(path).map_err(|e| h5_err(path, e))
    }

    /// Everything `molecule_info.h5` holds beyond the molecules themselves.
    pub struct MoleculeInput<'a> {
        pub meta: &'a CrMeta,
        pub table: &'a crate::solo::count::MoleculeTable,
        /// Whitelist index to barcode string (no `-1` suffix).
        pub barcode_of: &'a dyn Fn(u32) -> String,
        pub feature_ids: &'a [String],
        pub feature_names: &'a [String],
        pub feature_type: &'a str,
    }

    /// `metrics_json`: the keys `cellranger count` records, alphabetical as
    /// its serializer writes them.
    fn metrics_json(input: &MoleculeInput<'_>, feature_reads: u64, usable_reads: u64) -> String {
        let m = input.meta;
        let barcode = format!(
            r#"[{{"kind":"gel_bead","length":{},"offset":0,"read_type":"R1","whitelist":{{"name":"{}","part":null,"slide":null,"strand":"+","translation":false,"translation_whitelist_path":null}}}}]"#,
            m.cb_len, m.whitelist_name
        );
        let umi = format!(
            r#"[{{"length":{},"min_length":{},"offset":{},"read_type":"R1","whitelist":null}}]"#,
            m.umi_len,
            m.umi_len.min(10),
            m.cb_len
        );
        let rna = r#"{"length":null,"min_length":null,"offset":0,"read_type":"R2"}"#;
        let name = if m.umi_len == 12 {
            "SC3Pv3-polyA"
        } else {
            "SC3Pv2"
        };
        let desc = &m.chemistry_description;
        let rj = &m.reference_json;
        let fasta = json_str(rj, "fasta_hash").unwrap_or_default();
        let gtf_gz = json_str(rj, "gtf_hash.gz").unwrap_or_default();
        let gtf = json_str(rj, "gtf_hash").unwrap_or_default();
        let mkref = json_str(rj, "mkref_version").unwrap_or_default();
        format!(
            concat!(
                r#"{{"analysis_parameters":{{"filter_aggregates":true,"filter_high_occupancy_gems":true,"filter_probes":null,"include_introns":true}},"#,
                r#""cellranger_version":"{ver}","chemistry_barcode":{bc},"#,
                r#""chemistry_defs":{{"Gene Expression":{{"barcode":{bc},"barcode_extraction":null,"description":"{desc}","endedness":"three_prime","name":"{name}","rna":{rna},"rna2":null,"strandedness":"+","umi":{umi}}}}},"#,
                r#""chemistry_description":"{desc}","chemistry_endedness":"three_prime","chemistry_name":"{name}","chemistry_rna":{rna},"chemistry_rna2":null,"chemistry_strandedness":"+","chemistry_umi":{umi},"#,
                r#""gem_groups":{{"{gg}":{{"force_cells":null,"recovered_cells":null}}}},"#,
                r#""libraries":{{"0":{{"feature_read_pairs":{feat},"raw_read_pairs":{raw},"usable_read_pairs":{usable}}}}},"#,
                r#""molecule_info_type":"count","reference_fasta_hash":"{fasta}","reference_gtf_hash":"{gtf}","reference_gtf_hash.gz":"{gtf_gz}","reference_mkref_version":"{mkref}"}}"#,
            ),
            ver = m.software_version,
            bc = barcode,
            desc = desc,
            name = name,
            rna = rna,
            umi = umi,
            gg = m.gem_group,
            feat = feature_reads,
            raw = input.table.raw_read_pairs,
            usable = usable_reads,
            fasta = fasta,
            gtf = gtf,
            gtf_gz = gtf_gz,
            mkref = mkref,
        )
    }

    /// Write `molecule_info.h5` (file version 6) from the deduplicated
    /// molecules of the run.
    pub fn write_molecule_info(path: &Path, input: &MoleculeInput<'_>) -> Result<(), Error> {
        use hdf5_pure::MaxExtent;
        const CHUNK: u64 = 65_536;
        let table = input.table;
        let meta = input.meta;

        // Barcodes that carry a molecule, in lexicographic order (whitelist
        // indices are already lexicographic, so ascending index is enough).
        let mut cbs: Vec<u32> = table.molecules.iter().map(|m| m.cb).collect();
        cbs.sort_unstable();
        cbs.dedup();
        let barcodes: Vec<String> = cbs.iter().map(|&cb| (input.barcode_of)(cb)).collect();
        let idx_of = |cb: u32| cbs.binary_search(&cb).unwrap_or(0) as u64;

        // `(barcode, gene, umi)` order, as CellRanger writes them.
        let mut mols: Vec<&crate::solo::count::Molecule> = table.molecules.iter().collect();
        mols.sort_unstable_by_key(|m| (m.cb, m.gene, m.umi));
        let n = mols.len();
        let barcode_idx: Vec<u64> = mols.iter().map(|m| idx_of(m.cb)).collect();
        let feature_idx: Vec<u32> = mols.iter().map(|m| m.gene).collect();
        let umi: Vec<u32> = mols.iter().map(|m| m.umi as u32).collect();
        let count: Vec<u32> = mols.iter().map(|m| m.count).collect();
        let ones16 = vec![meta.gem_group as u16; n];
        let zeros16 = vec![0u16; n];
        // Bit 0 is CellRanger's "transcriptomic UMI" flag.
        let umi_type = vec![1u32; n];

        let feature_reads: u64 = count.iter().map(|&c| u64::from(c)).sum();
        let called: std::collections::HashSet<u32> = table.called.iter().copied().collect();
        let usable: u64 = mols
            .iter()
            .filter(|m| called.contains(&m.cb))
            .map(|m| u64::from(m.count))
            .sum();
        let mut pass: Vec<u64> = Vec::with_capacity(table.called.len() * 3);
        let mut pass_idx: Vec<u64> = table
            .called
            .iter()
            .filter(|cb| cbs.binary_search(cb).is_ok())
            .map(|&cb| idx_of(cb))
            .collect();
        pass_idx.sort_unstable();
        for i in &pass_idx {
            pass.extend([*i, 0, 0]);
        }

        let col = |ds: &mut DatasetBuilder| {
            if n > 0 {
                ds.with_chunks(&[CHUNK.min(n as u64)])
                    .with_shuffle()
                    .with_deflate(1)
                    .with_maxshape(&[MaxExtent::Unlimited]);
            }
        };

        let mut fb = FileBuilder::new();
        fb.set_attr("file_version", AttrValue::I64(6));
        fb.set_attr("filetype", AttrValue::VarLenString("molecule".to_string()));

        let ds = fb.create_dataset("barcode_idx");
        ds.with_u64_data(&barcode_idx);
        col(ds);

        let mut bi = fb.create_group("barcode_info");
        let ds = bi.create_dataset("genomes");
        strings_w(ds, std::slice::from_ref(&meta.genome), 256, path)?;
        ds.with_chunks(&[1])
            .with_shuffle()
            .with_deflate(1)
            .with_maxshape(&[MaxExtent::Unlimited]);
        let ds = bi.create_dataset("pass_filter");
        ds.with_u64_data(&pass)
            .with_shape(&[pass_idx.len() as u64, 3]);
        if !pass_idx.is_empty() {
            ds.with_chunks(&[CHUNK.min(pass_idx.len() as u64), 3])
                .with_shuffle()
                .with_deflate(1)
                .with_maxshape(&[MaxExtent::Unlimited, MaxExtent::Fixed(3)]);
        }
        fb.add_group(bi.finish());

        let ds = fb.create_dataset("barcodes");
        strings_w(ds, &barcodes, max_len(&barcodes), path)?;
        if !barcodes.is_empty() {
            ds.with_chunks(&[CHUNK.min(barcodes.len() as u64)])
                .with_shuffle()
                .with_deflate(1)
                .with_maxshape(&[MaxExtent::Unlimited]);
        }

        let ds = fb.create_dataset("count");
        ds.with_u32_data(&count);
        col(ds);
        let ds = fb.create_dataset("feature_idx");
        ds.with_u32_data(&feature_idx);
        col(ds);

        let mut features = fb.create_group("features");
        let nf = input.feature_ids.len();
        let tag = ["genome".to_string()];
        strings_w(features.create_dataset("_all_tag_keys"), &tag, 256, path)?;
        let packed = |ds: &mut DatasetBuilder, v: &[String]| -> Result<(), Error> {
            strings_w(ds, v, 256, path)?;
            if nf > 0 {
                ds.with_chunks(&[nf as u64]).with_shuffle().with_deflate(1);
            }
            Ok(())
        };
        packed(
            features.create_dataset("feature_type"),
            &vec![input.feature_type.to_string(); nf],
        )?;
        packed(
            features.create_dataset("genome"),
            &vec![meta.genome.clone(); nf],
        )?;
        packed(features.create_dataset("id"), input.feature_ids)?;
        packed(features.create_dataset("name"), input.feature_names)?;
        fb.add_group(features.finish());

        let ds = fb.create_dataset("gem_group");
        ds.with_u16_data(&ones16);
        col(ds);
        let ds = fb.create_dataset("library_idx");
        ds.with_u16_data(&zeros16);
        col(ds);

        let lib_info = format!(
            r#"[{{"library_id":0,"library_type":"Gene Expression","gem_group":{},"target_set_name":null}}]"#,
            meta.gem_group
        );
        strings_w(
            fb.create_dataset("library_info"),
            std::slice::from_ref(&lib_info),
            262_144,
            path,
        )?;
        // The genomes dataset above lists the genome; `library_info` repeats
        // the library table as JSON, as CellRanger does.
        fb.create_dataset("metrics_json")
            .with_vlen_strings(&[metrics_json(input, feature_reads, usable).as_str()])
            .with_shape(&[]);

        let ds = fb.create_dataset("umi");
        ds.with_u32_data(&umi);
        col(ds);
        let ds = fb.create_dataset("umi_type");
        ds.with_u32_data(&umi_type);
        col(ds);

        fb.write(path).map_err(|e| h5_err(path, e))
    }

    /// Write `<outs>/molecule_info.h5` for a CellRanger-layout run that kept
    /// its molecules. A no-op when the counting pass kept none.
    pub fn write_molecule_info_for_run(
        sctx: &crate::solo::SoloContext,
        params: &Parameters,
    ) -> Result<(), Error> {
        let Some(table) = sctx.molecules.lock().unwrap().take() else {
            return Ok(());
        };
        let meta = CrMeta::from_params(params);
        let solo_dir = params
            .solo_out_file_names
            .first()
            .map_or("Solo.out/", String::as_str);
        let out = if sctx.features.len() > 1 {
            params.output_path(&format!("{solo_dir}{}/", sctx.features[0].dir_name()))
        } else {
            params.output_path(solo_dir)
        }
        .join("molecule_info.h5");
        let wl = &sctx.whitelist;
        let barcode_of = |cb: u32| wl.barcode_string(cb).unwrap_or_default();
        write_molecule_info(
            &out,
            &MoleculeInput {
                meta: &meta,
                table: &table,
                barcode_of: &barcode_of,
                feature_ids: &sctx.gene_ann.gene_ids,
                feature_names: &sctx.gene_ann.gene_names,
                feature_type: "Gene Expression",
            },
        )?;
        log::info!(
            "STARsolo: wrote {} ({} molecules)",
            out.display(),
            table.molecules.len()
        );
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use hdf5_pure::File;

        fn write(dir: &Path, name: &str, body: &str) {
            std::fs::write(dir.join(name), body).unwrap();
        }

        /// 3 genes × 4 barcodes, entries deliberately out of column order.
        fn fixture(dir: &Path) {
            write(
                dir,
                "features.tsv",
                "G1\tAlpha\tGene Expression\nG2\tBeta\tGene Expression\nG3\tGamma\tGene Expression\n",
            );
            write(dir, "barcodes.tsv", "AAAA\nCCCC\nGGGG\nTTTT\n");
            write(
                dir,
                "matrix.mtx",
                "%%MatrixMarket matrix coordinate integer general\n%\n3 4 4\n3 4 7\n1 1 5\n2 1 1\n2 3 2\n",
            );
        }

        #[test]
        fn mtx_to_csc() {
            let tmp = tempfile::tempdir().unwrap();
            fixture(tmp.path());
            let m = read_mtx_dir(tmp.path(), &MtxNames::default(), false).unwrap();
            assert_eq!((m.n_features, m.n_barcodes), (3, 4));
            assert_eq!(m.indptr, [0, 2, 2, 3, 4]);
            assert_eq!(m.indices, [0, 1, 1, 2]);
            assert_eq!(m.data, [5, 1, 2, 7]);
            assert_eq!(m.feature_names, ["Alpha", "Beta", "Gamma"]);
            assert_eq!(m.barcodes, ["AAAA", "CCCC", "GGGG", "TTTT"]);
        }

        #[test]
        fn mtx_gzip_and_two_column_features() {
            use std::io::Write as _;
            let tmp = tempfile::tempdir().unwrap();
            let gz = |name: &str, body: &str| {
                let f = std::fs::File::create(tmp.path().join(format!("{name}.gz"))).unwrap();
                let mut w = flate2::write::GzEncoder::new(f, flate2::Compression::default());
                w.write_all(body.as_bytes()).unwrap();
                w.finish().unwrap();
            };
            gz("features.tsv", "G1\n");
            gz("barcodes.tsv", "AC\n");
            gz(
                "matrix.mtx",
                "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 9\n",
            );
            let m = read_mtx_dir(tmp.path(), &MtxNames::default(), true).unwrap();
            assert_eq!(m.feature_names, ["G1"]);
            assert_eq!(m.feature_types, ["Gene Expression"]);
            assert_eq!(m.data, [9]);
        }

        #[test]
        fn mtx_rejects_inconsistent_input() {
            let tmp = tempfile::tempdir().unwrap();
            fixture(tmp.path());
            // Wrong nnz.
            write(
                tmp.path(),
                "matrix.mtx",
                "%%MatrixMarket matrix coordinate integer general\n3 4 2\n1 1 5\n",
            );
            assert!(read_mtx_dir(tmp.path(), &MtxNames::default(), false).is_err());
            // Real-valued (UniqueAndMult-style) matrices are not converted.
            write(
                tmp.path(),
                "matrix.mtx",
                "%%MatrixMarket matrix coordinate real general\n3 4 1\n1 1 0.5\n",
            );
            assert!(read_mtx_dir(tmp.path(), &MtxNames::default(), false).is_err());
            // Row out of bounds.
            write(
                tmp.path(),
                "matrix.mtx",
                "%%MatrixMarket matrix coordinate integer general\n3 4 1\n4 1 5\n",
            );
            assert!(read_mtx_dir(tmp.path(), &MtxNames::default(), false).is_err());
        }

        #[test]
        fn h5_roundtrip_cellranger_layout() {
            let tmp = tempfile::tempdir().unwrap();
            fixture(tmp.path());
            let m = read_mtx_dir(tmp.path(), &MtxNames::default(), false).unwrap();
            let out = tmp.path().join("raw_feature_bc_matrix.h5");
            write_10x_h5(&m, &out).unwrap();

            let f = File::open(&out).unwrap();
            let root = f.root().attrs().unwrap();
            assert_eq!(root["filetype"].as_str(), Some("matrix"));
            assert_eq!(root["version"].as_i64(), Some(2));
            let ds = |p: &str| f.dataset(p).unwrap();
            assert_eq!(ds("matrix/shape").read_i32().unwrap(), [3, 4]);
            assert_eq!(ds("matrix/data").read_i32().unwrap(), m.data);
            assert_eq!(ds("matrix/indices").read_i64().unwrap(), m.indices);
            assert_eq!(ds("matrix/indptr").read_i64().unwrap(), m.indptr);
            assert_eq!(ds("matrix/barcodes").read_string().unwrap(), m.barcodes);
            assert_eq!(
                ds("matrix/features/id").read_string().unwrap(),
                m.feature_ids
            );
            assert_eq!(
                ds("matrix/features/name").read_string().unwrap(),
                m.feature_names
            );
            assert_eq!(
                ds("matrix/features/feature_type").read_string().unwrap(),
                m.feature_types
            );
            assert_eq!(
                ds("matrix/features/_all_tag_keys").read_string().unwrap(),
                ["genome"]
            );
            assert_eq!(ds("matrix/features/genome").shape().unwrap(), [3]);
        }

        fn meta() -> CrMeta {
            CrMeta {
                chemistry_description: "Single Cell 3' v3 (polyA)".into(),
                library_id: "lib".into(),
                genome: "GRCh38".into(),
                software_version: "t 1".into(),
                gem_group: 1,
                reference_json:
                    r#"{"fasta_hash": "abc", "genomes": ["GRCh38"], "gtf_hash.gz": "def"}"#.into(),
                whitelist_name: "wl".into(),
                umi_len: 12,
                cb_len: 16,
            }
        }

        #[test]
        fn cellranger_flavour_attributes_and_genome() {
            let tmp = tempfile::tempdir().unwrap();
            fixture(tmp.path());
            let m = read_mtx_dir(tmp.path(), &MtxNames::default(), false).unwrap();
            for raw in [true, false] {
                let out = tmp.path().join(format!("cr_{raw}.h5"));
                write_10x_h5_cr(&m, &out, &meta(), raw).unwrap();
                let f = File::open(&out).unwrap();
                let a = f.root().attrs().unwrap();
                assert_eq!(a["filetype"].as_str(), Some("matrix"));
                assert_eq!(
                    a["chemistry_description"].as_str(),
                    Some("Single Cell 3' v3 (polyA)")
                );
                assert!(a.contains_key("original_gem_groups"));
                let g = f
                    .dataset("matrix/features/genome")
                    .unwrap()
                    .read_string()
                    .unwrap();
                assert_eq!(g, ["GRCh38"; 3]);
                assert_eq!(
                    f.dataset("matrix/data").unwrap().read_i32().unwrap(),
                    m.data
                );
                assert_eq!(
                    f.dataset("matrix/barcodes").unwrap().read_string().unwrap(),
                    m.barcodes
                );
            }
        }

        #[test]
        fn molecule_info_roundtrip() {
            use crate::solo::count::{Molecule, MoleculeTable};
            let tmp = tempfile::tempdir().unwrap();
            let table = MoleculeTable {
                molecules: vec![
                    Molecule {
                        cb: 5,
                        gene: 1,
                        umi: 9,
                        count: 3,
                    },
                    Molecule {
                        cb: 2,
                        gene: 0,
                        umi: 7,
                        count: 1,
                    },
                    Molecule {
                        cb: 2,
                        gene: 0,
                        umi: 3,
                        count: 2,
                    },
                ],
                called: vec![2],
                raw_read_pairs: 100,
            };
            let ids = vec!["G1".to_string(), "G2".to_string()];
            let names = ids.clone();
            let bc = |cb: u32| format!("BC{cb}");
            let out = tmp.path().join("molecule_info.h5");
            let m = meta();
            write_molecule_info(
                &out,
                &MoleculeInput {
                    meta: &m,
                    table: &table,
                    barcode_of: &bc,
                    feature_ids: &ids,
                    feature_names: &names,
                    feature_type: "Gene Expression",
                },
            )
            .unwrap();
            let f = File::open(&out).unwrap();
            let u32s = |n: &str| f.dataset(n).unwrap().read_u32().unwrap();
            // Sorted (barcode, gene, umi); barcode 2 -> idx 0, barcode 5 -> idx 1.
            assert_eq!(u32s("umi"), [3, 7, 9]);
            assert_eq!(u32s("count"), [2, 1, 3]);
            assert_eq!(u32s("feature_idx"), [0, 0, 1]);
            assert_eq!(
                f.dataset("barcode_idx").unwrap().read_u64().unwrap(),
                [0, 0, 1]
            );
            assert_eq!(
                f.dataset("barcodes").unwrap().read_string().unwrap(),
                ["BC2", "BC5"]
            );
            assert_eq!(
                f.dataset("barcode_info/pass_filter")
                    .unwrap()
                    .read_u64()
                    .unwrap(),
                [0, 0, 0]
            );
            assert_eq!(f.root().attrs().unwrap()["file_version"].as_i64(), Some(6));
        }

        #[test]
        fn h5_empty_matrix() {
            let tmp = tempfile::tempdir().unwrap();
            let m = CscMatrix {
                n_features: 1,
                n_barcodes: 2,
                indptr: vec![0, 0, 0],
                barcodes: vec!["A".into(), "C".into()],
                feature_ids: vec!["G1".into()],
                feature_names: vec!["G1".into()],
                feature_types: vec!["Gene Expression".into()],
                ..CscMatrix::default()
            };
            let out = tmp.path().join("empty.h5");
            write_10x_h5(&m, &out).unwrap();
            let f = File::open(&out).unwrap();
            assert_eq!(
                f.dataset("matrix/data").unwrap().read_i32().unwrap(),
                Vec::<i32>::new()
            );
            assert_eq!(
                f.dataset("matrix/indptr").unwrap().read_i64().unwrap(),
                [0, 0, 0]
            );
        }
    }
}
