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

#[cfg(feature = "hdf5-out")]
pub use imp::{CscMatrix, read_mtx_dir, write_10x_h5};

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
        for dir in feature_dirs {
            let feature_dir = params.output_path(&format!("{solo_dir}{dir}/"));
            for (sub, h5_name) in [
                ("raw", "raw_feature_bc_matrix.h5"),
                ("filtered", "filtered_feature_bc_matrix.h5"),
            ] {
                let mtx_dir = feature_dir.join(sub);
                // `filtered/` exists only when cells were called.
                if !with_gz(&mtx_dir.join(&names.matrix), gzip).exists() {
                    continue;
                }
                let m = read_mtx_dir(&mtx_dir, names, gzip)?;
                let out = feature_dir.join(h5_name);
                write_10x_h5(&m, &out)?;
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
        let mut fb = FileBuilder::new();
        fb.set_attr("filetype", AttrValue::AsciiString("matrix".to_string()));
        fb.set_attr("version", AttrValue::I64(2));
        fb.set_attr(
            "software_version",
            AttrValue::AsciiString(format!("rustar-aligner {}", env!("CARGO_PKG_VERSION"))),
        );

        let mut matrix = fb.create_group("matrix");
        strings(matrix.create_dataset("barcodes"), &m.barcodes, path)?;
        let ds = matrix.create_dataset("data");
        ds.with_i32_data(&m.data);
        compress(ds, m.data.len());
        let ds = matrix.create_dataset("indices");
        ds.with_i64_data(&m.indices);
        compress(ds, m.indices.len());
        matrix.create_dataset("indptr").with_i64_data(&m.indptr);
        let shape = [
            i32::try_from(m.n_features).map_err(|e| h5_err(path, e))?,
            i32::try_from(m.n_barcodes).map_err(|e| h5_err(path, e))?,
        ];
        matrix.create_dataset("shape").with_i32_data(&shape);

        let mut features = matrix.create_group("features");
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
        matrix.add_group(features.finish());
        fb.add_group(matrix.finish());

        fb.write(path).map_err(|e| h5_err(path, e))
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
            assert!(
                f.dataset("matrix/data")
                    .unwrap()
                    .read_i32()
                    .unwrap()
                    .is_empty()
            );
            assert_eq!(
                f.dataset("matrix/indptr").unwrap().read_i64().unwrap(),
                [0, 0, 0]
            );
        }
    }
}
