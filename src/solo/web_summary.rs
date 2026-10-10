//! `outs/web_summary.html` for `--soloOutLayout CellRanger`: one self-contained
//! HTML report of the run.
//!
//! The page has no scripts and no external resources. Every plot is an inline
//! SVG drawn here, and the colours are CSS custom properties so the page follows
//! the reader's light or dark preference. The design is rustar's own.
//!
//! Inputs are the metrics the run computed (the same strings that go into
//! `metrics_summary.csv`), the per-barcode UMI totals for the rank plot, an
//! optional sequencing-depth curve, and, when `outs/analysis/` exists, the
//! embeddings, clusters and differential expression files found there.

#![allow(clippy::many_single_char_names)]

use crate::error::Error;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Run description shown in the page header.
#[derive(Debug, Clone, Default)]
pub struct RunInfo {
    pub sample_id: String,
    pub chemistry: String,
    pub reference: String,
    pub version: String,
}

impl RunInfo {
    /// Header facts from the run parameters and the reference folder's
    /// `reference.json` (next to the genome directory), when there is one.
    pub fn from_params(params: &crate::params::Parameters) -> Self {
        let reference_json = params
            .genome_dir
            .parent()
            .map(|p| p.join("reference.json"))
            .and_then(|p| std::fs::read_to_string(p).ok())
            .unwrap_or_default();
        let chemistry = match (params.solo_umi_len, params.solo_barcode_mate) {
            (12, 0) => "Single Cell 3' v3 (polyA)",
            (10, 0) => "Single Cell 3' v2",
            (_, 0) => "Single Cell 3' (custom)",
            (10, _) => "Single Cell 5' v2",
            _ => "Single Cell 5' (custom)",
        };
        Self {
            sample_id: params.solo_out_sample_id.clone(),
            chemistry: chemistry.to_string(),
            reference: json_array_first(&reference_json, "genomes").unwrap_or_default(),
            version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }
}

/// First string of `"key": ["a", ...]` in flat JSON text.
fn json_array_first(text: &str, key: &str) -> Option<String> {
    let at = text.find(&format!("\"{key}\""))?;
    let rest = &text[at + key.len() + 2..];
    let rest = rest.trim_start().strip_prefix(':')?.trim_start();
    let rest = rest.strip_prefix('[')?.trim_start().strip_prefix('"')?;
    Some(rest[..rest.find('"')?].to_string())
}

/// Per-barcode UMI totals for the rank plot.
#[derive(Debug, Clone, Default)]
pub struct RankData {
    /// `(umis, is_cell)` for every barcode with at least one UMI, by UMIs descending.
    pub ranked: Vec<(u64, bool)>,
    /// UMI total of each called cell, by barcode string (as in the matrix).
    pub cell_umis: HashMap<String, u64>,
}

/// One point of the sequencing-depth curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint {
    pub reads_per_cell: f64,
    /// Sequencing saturation in percent.
    pub saturation: f64,
    pub median_genes: f64,
}

/// Expected saturation and median genes per cell when only a fraction of the
/// reads is kept. Each molecule with `count` reads survives with probability
/// `1 - (1-f)^count`; the expectation is used, so no random numbers are needed.
///
/// `molecules` must be `(cb, gene, umi)`-ascending; `called` lists the cells'
/// whitelist indices ascending. Returns nothing without cells or molecules.
pub fn saturation_curve(
    molecules: &[crate::solo::count::Molecule],
    called: &[u32],
    total_reads: u64,
) -> Vec<CurvePoint> {
    if molecules.is_empty() || called.is_empty() {
        return Vec::new();
    }
    let n_cells = called.len() as f64;
    let mut out = Vec::new();
    for step in 1..=10u32 {
        let f = f64::from(step) / 10.0;
        let (mut reads, mut umis) = (0.0f64, 0.0f64);
        let mut per_cell: HashMap<u32, f64> = HashMap::new();
        let mut i = 0;
        while i < molecules.len() {
            // One (cb, gene) group: the gene is seen unless every molecule is lost.
            let (cb, gene) = (molecules[i].cb, molecules[i].gene);
            let mut miss = 1.0f64;
            while i < molecules.len() && molecules[i].cb == cb && molecules[i].gene == gene {
                let keep = 1.0 - (1.0 - f).powi(molecules[i].count as i32);
                reads += f * f64::from(molecules[i].count);
                umis += keep;
                miss *= 1.0 - keep;
                i += 1;
            }
            if called.binary_search(&cb).is_ok() {
                *per_cell.entry(cb).or_insert(0.0) += 1.0 - miss;
            }
        }
        let mut genes: Vec<f64> = per_cell.values().copied().collect();
        genes.resize(called.len(), 0.0);
        genes.sort_by(f64::total_cmp);
        let med = if genes.len() % 2 == 1 {
            genes[genes.len() / 2]
        } else {
            f64::midpoint(genes[genes.len() / 2 - 1], genes[genes.len() / 2])
        };
        out.push(CurvePoint {
            reads_per_cell: f * total_reads as f64 / n_cells,
            saturation: if reads > 0.0 {
                100.0 * (1.0 - umis / reads).max(0.0)
            } else {
                0.0
            },
            median_genes: med,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// Reading what the run already wrote
// ---------------------------------------------------------------------------

fn invalid(path: &Path, msg: &str) -> Error {
    Error::io(
        std::io::Error::new(std::io::ErrorKind::InvalidData, msg.to_string()),
        path,
    )
}

/// Open `path`, or `path.gz` / `path` without `.gz`, decompressing as needed.
fn open_any(path: &Path) -> Result<Box<dyn BufRead>, Error> {
    let mut candidates = vec![path.to_path_buf()];
    let mut gz = path.as_os_str().to_owned();
    gz.push(".gz");
    candidates.push(PathBuf::from(gz));
    if path.extension().is_some_and(|e| e == "gz") {
        candidates.push(path.with_extension(""));
    }
    for c in &candidates {
        if let Ok(f) = std::fs::File::open(c) {
            let is_gz = c.extension().is_some_and(|e| e == "gz");
            let r: Box<dyn Read> = if is_gz {
                Box::new(flate2::read::MultiGzDecoder::new(f))
            } else {
                Box::new(f)
            };
            return Ok(Box::new(BufReader::with_capacity(1 << 20, r)));
        }
    }
    Err(Error::io(
        std::io::Error::new(std::io::ErrorKind::NotFound, "file not found"),
        path,
    ))
}

/// Split one CSV line, honouring double quotes.
fn split_csv(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    for c in line.trim_end_matches(['\r', '\n']).chars() {
        match c {
            '"' => quoted = !quoted,
            ',' if !quoted => out.push(std::mem::take(&mut cur)),
            _ => cur.push(c),
        }
    }
    out.push(cur);
    out
}

/// Read `metrics_summary.csv` back into `(name, value)` pairs.
pub fn read_metrics_csv(path: &Path) -> Result<Vec<(String, String)>, Error> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
    let mut lines = text.lines();
    let names = split_csv(lines.next().ok_or_else(|| invalid(path, "empty file"))?);
    let values = split_csv(lines.next().ok_or_else(|| invalid(path, "no value row"))?);
    Ok(names.into_iter().zip(values).collect())
}

/// Per-barcode UMI totals from the raw matrix, with cells taken to be the
/// barcodes of the filtered matrix.
pub fn read_rank_data(outs_dir: &Path) -> Result<RankData, Error> {
    let raw = outs_dir.join("raw_feature_bc_matrix");
    let filt = outs_dir.join("filtered_feature_bc_matrix");
    let read_barcodes = |p: &Path| -> Result<Vec<String>, Error> {
        let mut v = Vec::new();
        for l in open_any(p)?.lines() {
            v.push(l.map_err(|e| Error::io(e, p))?.trim().to_string());
        }
        Ok(v)
    };
    let raw_bc = read_barcodes(&raw.join("barcodes.tsv.gz"))?;
    let cells: std::collections::HashSet<String> = read_barcodes(&filt.join("barcodes.tsv.gz"))
        .unwrap_or_default()
        .into_iter()
        .collect();

    let mpath = raw.join("matrix.mtx.gz");
    let mut totals = vec![0u64; raw_bc.len()];
    let mut header_seen = false;
    for l in open_any(&mpath)?.lines() {
        let l = l.map_err(|e| Error::io(e, &mpath))?;
        if l.starts_with('%') {
            continue;
        }
        if !header_seen {
            header_seen = true;
            continue;
        }
        let mut it = l.split_ascii_whitespace();
        let (_g, c, n) = (it.next(), it.next(), it.next());
        if let (Some(c), Some(n)) = (c, n)
            && let (Ok(c), Ok(n)) = (c.parse::<usize>(), n.parse::<u64>())
            && c >= 1
            && c <= totals.len()
        {
            totals[c - 1] += n;
        }
    }
    let mut ranked: Vec<(u64, bool)> = Vec::new();
    let mut cell_umis = HashMap::new();
    for (bc, &n) in raw_bc.iter().zip(&totals) {
        if n == 0 {
            continue;
        }
        let is_cell = cells.contains(bc);
        if is_cell {
            cell_umis.insert(bc.clone(), n);
        }
        ranked.push((n, is_cell));
    }
    ranked.sort_by_key(|a| std::cmp::Reverse(a.0));
    Ok(RankData { ranked, cell_umis })
}

/// Build the page from files already in `outs_dir` (no curve): reads
/// `metrics_summary.csv` and the raw and filtered matrices, then writes
/// `web_summary.html` to `out_path`.
pub fn write_web_summary_from_outs(
    outs_dir: &Path,
    out_path: &Path,
    info: &RunInfo,
    curve: Option<&[CurvePoint]>,
) -> Result<(), Error> {
    let metrics = read_metrics_csv(&outs_dir.join("metrics_summary.csv"))?;
    let ranks = read_rank_data(outs_dir)?;
    let html = render_web_summary(outs_dir, &metrics, info, &ranks, curve);
    std::fs::write(out_path, html).map_err(|e| Error::io(e, out_path))
}

/// Write `outs_dir/web_summary.html` and return its path.
pub fn write_web_summary(
    outs_dir: &Path,
    metrics: &[(String, String)],
    info: &RunInfo,
    ranks: &RankData,
    curve: Option<&[CurvePoint]>,
) -> Result<PathBuf, Error> {
    let path = outs_dir.join("web_summary.html");
    let html = render_web_summary(outs_dir, metrics, info, ranks, curve);
    std::fs::write(&path, html).map_err(|e| Error::io(e, &path))?;
    Ok(path)
}

// ---------------------------------------------------------------------------
// Analysis files
// ---------------------------------------------------------------------------

struct Embedding {
    name: &'static str,
    points: Vec<(String, f64, f64)>,
}

/// `(gene, log2 fold change, adjusted p)`.
type GeneHit = (String, f64, f64);

struct Analysis {
    clusters: HashMap<String, u32>,
    embeddings: Vec<Embedding>,
    /// Per cluster: `(gene, log2 fold change, adjusted p)`, best first.
    top_genes: Vec<(u32, Vec<GeneHit>)>,
}

fn read_csv_rows(path: &Path) -> Option<(Vec<String>, Vec<Vec<String>>)> {
    let r = open_any(path).ok()?;
    let mut lines = r.lines();
    let head = split_csv(&lines.next()?.ok()?);
    let rows = lines
        .map_while(Result::ok)
        .filter(|l| !l.is_empty())
        .map(|l| split_csv(&l))
        .collect();
    Some((head, rows))
}

fn first_subdir(dir: &Path, prefix: &str) -> Option<PathBuf> {
    let mut found: Vec<PathBuf> = std::fs::read_dir(dir)
        .ok()?
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| {
            p.is_dir()
                && p.file_name()
                    .is_some_and(|n| n.to_string_lossy().starts_with(prefix))
        })
        .collect();
    found.sort();
    found.into_iter().next()
}

fn load_analysis(outs_dir: &Path) -> Option<Analysis> {
    let a = outs_dir.join("analysis");
    if !a.is_dir() {
        return None;
    }
    let clusters: HashMap<String, u32> = read_csv_rows(
        &a.join("clustering")
            .join("gene_expression_graphclust")
            .join("clusters.csv"),
    )
    .map(|(_, rows)| {
        rows.into_iter()
            .filter(|r| r.len() >= 2)
            .filter_map(|r| Some((r[0].clone(), r[1].parse().ok()?)))
            .collect()
    })
    .unwrap_or_default();

    let mut embeddings = Vec::new();
    for (dir, name) in [("tsne", "t-SNE"), ("umap", "UMAP")] {
        let Some(sub) = first_subdir(&a.join(dir), "") else {
            continue;
        };
        if let Some((_, rows)) = read_csv_rows(&sub.join("projection.csv")) {
            let points: Vec<(String, f64, f64)> = rows
                .into_iter()
                .filter(|r| r.len() >= 3)
                .filter_map(|r| Some((r[0].clone(), r[1].parse().ok()?, r[2].parse().ok()?)))
                .collect();
            if !points.is_empty() {
                embeddings.push(Embedding { name, points });
            }
        }
    }

    let mut top_genes = Vec::new();
    if let Some((head, rows)) = read_csv_rows(
        &a.join("diffexp")
            .join("gene_expression_graphclust")
            .join("differential_expression.csv"),
    ) {
        let mut k = 1u32;
        loop {
            let lfc = head
                .iter()
                .position(|h| *h == format!("Cluster {k} Log2 fold change"));
            let pv = head
                .iter()
                .position(|h| *h == format!("Cluster {k} Adjusted p value"));
            let (Some(lfc), Some(pv)) = (lfc, pv) else {
                break;
            };
            let mut genes: Vec<GeneHit> = rows
                .iter()
                .filter_map(|r| {
                    let name = r.get(1)?.clone();
                    let l: f64 = r.get(lfc)?.parse().ok()?;
                    let p: f64 = r.get(pv)?.parse().ok()?;
                    (p < 0.05 && l > 0.0).then_some((name, l, p))
                })
                .collect();
            genes.sort_by(|a, b| b.1.total_cmp(&a.1));
            genes.truncate(5);
            top_genes.push((k, genes));
            k += 1;
        }
    }
    if embeddings.is_empty() && top_genes.is_empty() {
        return None;
    }
    Some(Analysis {
        clusters,
        embeddings,
        top_genes,
    })
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn metric<'a>(metrics: &'a [(String, String)], name: &str) -> &'a str {
    metrics
        .iter()
        .find(|(k, _)| k == name)
        .map_or("n/a", |(_, v)| v.as_str())
}

fn tick_label(v: f64) -> String {
    let a = v.abs();
    if a >= 1e6 {
        format!("{}M", trim_num(v / 1e6))
    } else if a >= 1e3 {
        format!("{}k", trim_num(v / 1e3))
    } else {
        trim_num(v)
    }
}

fn trim_num(v: f64) -> String {
    let s = format!("{v:.2}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// A plot area with margins, in SVG user units.
struct Frame {
    w: f64,
    h: f64,
    ml: f64,
    mr: f64,
    mt: f64,
    mb: f64,
}

impl Frame {
    fn x(&self, t: f64) -> f64 {
        self.ml + t * (self.w - self.ml - self.mr)
    }
    fn y(&self, t: f64) -> f64 {
        self.h - self.mb - t * (self.h - self.mt - self.mb)
    }
    fn open(&self, label: &str) -> String {
        format!(
            "<svg viewBox=\"0 0 {} {}\" role=\"img\" aria-label=\"{}\" preserveAspectRatio=\"xMidYMid meet\">",
            self.w,
            self.h,
            esc(label)
        )
    }
    fn box_and_titles(&self, xt: &str, yt: &str) -> String {
        let mut s = format!(
            "<rect class=\"frame\" x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>",
            self.ml,
            self.mt,
            self.w - self.ml - self.mr,
            self.h - self.mt - self.mb
        );
        let _ = write!(
            s,
            "<text class=\"axt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>",
            self.ml + (self.w - self.ml - self.mr) / 2.0,
            self.h - 6.0,
            esc(xt)
        );
        let _ = write!(
            s,
            "<text class=\"axt\" transform=\"translate(14 {:.1}) rotate(-90)\" text-anchor=\"middle\">{}</text>",
            self.mt + (self.h - self.mt - self.mb) / 2.0,
            esc(yt)
        );
        s
    }
}

/// Decade ticks `10^k` inside `[lo, hi]`.
fn decades(lo: f64, hi: f64) -> Vec<f64> {
    let mut v = Vec::new();
    let mut k = lo.log10().floor() as i32;
    while 10f64.powi(k) <= hi * 1.0001 {
        if 10f64.powi(k) >= lo * 0.9999 {
            v.push(10f64.powi(k));
        }
        k += 1;
    }
    v
}

/// Round-number ticks inside `[lo, hi]`, about `n` of them.
fn linear_ticks(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let span = (hi - lo).max(f64::MIN_POSITIVE);
    let raw = span / n.max(1) as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * mag)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * mag);
    let mut t = (lo / step).ceil() * step;
    let mut v = Vec::new();
    while t <= hi + step * 1e-9 {
        v.push(t);
        t += step;
    }
    v
}

fn rank_plot(ranks: &RankData) -> String {
    let n = ranks.ranked.len();
    if n < 2 {
        return String::new();
    }
    let f = Frame {
        w: 560.0,
        h: 360.0,
        ml: 56.0,
        mr: 14.0,
        mt: 12.0,
        mb: 42.0,
    };
    let max_u = ranks.ranked[0].0.max(10) as f64;
    let min_u = ranks.ranked[n - 1].0.max(1) as f64;
    let (xlo, xhi) = (0.0, (n as f64).log10().max(1.0));
    let ylo = min_u.log10().floor();
    let yhi = max_u.log10().ceil().max(ylo + 1.0);
    let px = |rank: f64| f.x((rank.log10() - xlo) / (xhi - xlo));
    let py = |u: f64| f.y((u.max(1.0).log10() - ylo) / (yhi - ylo));

    let mut s = f.open("Barcode rank plot, UMI counts against barcode rank");
    for t in decades(1.0, n as f64) {
        let x = px(t);
        let _ = write!(
            s,
            "<line class=\"grid\" x1=\"{x:.1}\" y1=\"{:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\"/><text class=\"tk\" x=\"{x:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>",
            f.mt,
            f.h - f.mb,
            f.h - f.mb + 14.0,
            tick_label(t)
        );
    }
    for t in decades(10f64.powf(ylo), 10f64.powf(yhi)) {
        let y = py(t);
        let _ = write!(
            s,
            "<line class=\"grid\" x1=\"{:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\"/><text class=\"tk\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\">{}</text>",
            f.ml,
            f.w - f.mr,
            f.ml - 5.0,
            y + 3.5,
            tick_label(t)
        );
    }
    // Log-spaced sample of ranks (all of the first 50), so 700k barcodes stay light.
    let mut picks: Vec<usize> = Vec::new();
    let mut last_bucket = i64::MIN;
    for i in 0..n {
        let bucket = if i < 50 {
            i as i64 - 1_000_000
        } else {
            (((i + 1) as f64).log10() * 120.0) as i64
        };
        if bucket != last_bucket {
            picks.push(i);
            last_bucket = bucket;
        }
    }
    if picks.last() != Some(&(n - 1)) {
        picks.push(n - 1);
    }
    // Runs of equal state, each sharing its first point with the previous run.
    let mut start = 0;
    while start + 1 < picks.len() {
        let state = ranks.ranked[picks[start + 1]].1;
        let mut end = start + 1;
        while end + 1 < picks.len() && ranks.ranked[picks[end + 1]].1 == state {
            end += 1;
        }
        let pts: Vec<String> = picks[start..=end]
            .iter()
            .map(|&i| {
                format!(
                    "{:.1},{:.1}",
                    px((i + 1) as f64),
                    py(ranks.ranked[i].0 as f64)
                )
            })
            .collect();
        let _ = write!(
            s,
            "<polyline class=\"{}\" points=\"{}\"/>",
            if state { "rk-cell" } else { "rk-bg" },
            pts.join(" ")
        );
        start = end;
    }
    s.push_str(&f.box_and_titles("Barcode rank", "UMI counts"));
    s.push_str("</svg>");
    s
}

fn curve_plot(
    points: &[CurvePoint],
    ylabel: &str,
    pick: fn(&CurvePoint) -> f64,
    ymax: Option<f64>,
) -> String {
    let f = Frame {
        w: 400.0,
        h: 260.0,
        ml: 56.0,
        mr: 14.0,
        mt: 12.0,
        mb: 42.0,
    };
    let xhi = points
        .iter()
        .map(|p| p.reads_per_cell)
        .fold(1.0f64, f64::max);
    let yhi = ymax.unwrap_or_else(|| points.iter().map(pick).fold(1.0f64, f64::max) * 1.08);
    let mut s = f.open(&format!("{ylabel} against mean reads per cell"));
    for t in linear_ticks(0.0, xhi, 4) {
        let x = f.x(t / xhi);
        let _ = write!(
            s,
            "<line class=\"grid\" x1=\"{x:.1}\" y1=\"{:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\"/><text class=\"tk\" x=\"{x:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>",
            f.mt,
            f.h - f.mb,
            f.h - f.mb + 14.0,
            tick_label(t)
        );
    }
    for t in linear_ticks(0.0, yhi, 4) {
        let y = f.y(t / yhi);
        let _ = write!(
            s,
            "<line class=\"grid\" x1=\"{:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\"/><text class=\"tk\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\">{}</text>",
            f.ml,
            f.w - f.mr,
            f.ml - 5.0,
            y + 3.5,
            tick_label(t)
        );
    }
    let pts: Vec<String> = points
        .iter()
        .map(|p| {
            format!(
                "{:.1},{:.1}",
                f.x(p.reads_per_cell / xhi),
                f.y(pick(p) / yhi)
            )
        })
        .collect();
    let _ = write!(s, "<polyline class=\"line\" points=\"{}\"/>", pts.join(" "));
    for p in &pts {
        let (x, y) = p.split_once(',').unwrap_or(("0", "0"));
        let _ = write!(s, "<circle class=\"dot\" cx=\"{x}\" cy=\"{y}\" r=\"3\"/>");
    }
    s.push_str(&f.box_and_titles("Mean reads per cell", ylabel));
    s.push_str("</svg>");
    s
}

/// Viridis-like ramp, `t` in 0..1.
fn ramp(t: f64) -> String {
    const STOPS: [(f64, f64, f64); 5] = [
        (68.0, 1.0, 84.0),
        (59.0, 82.0, 139.0),
        (33.0, 145.0, 140.0),
        (94.0, 201.0, 98.0),
        (253.0, 231.0, 37.0),
    ];
    let t = t.clamp(0.0, 1.0) * 4.0;
    let i = (t.floor() as usize).min(3);
    let u = t - i as f64;
    let c = |a: f64, b: f64| (a + (b - a) * u).round() as u8;
    let (a, b) = (STOPS[i], STOPS[i + 1]);
    format!("#{:02x}{:02x}{:02x}", c(a.0, b.0), c(a.1, b.1), c(a.2, b.2))
}

enum Colouring<'a> {
    Cluster(&'a HashMap<String, u32>),
    Umi(&'a HashMap<String, u64>),
}

fn scatter(emb: &Embedding, colouring: &Colouring<'_>, label: &str) -> String {
    let f = Frame {
        w: 400.0,
        h: 330.0,
        ml: 24.0,
        mr: 8.0,
        mt: 8.0,
        mb: 30.0,
    };
    let (mut x0, mut x1, mut y0, mut y1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for (_, x, y) in &emb.points {
        x0 = x0.min(*x);
        x1 = x1.max(*x);
        y0 = y0.min(*y);
        y1 = y1.max(*y);
    }
    let (dx, dy) = ((x1 - x0).max(1e-9), (y1 - y0).max(1e-9));
    let (px, py) = (
        |x: f64| f.x(0.02 + 0.96 * (x - x0) / dx),
        |y: f64| f.y(0.02 + 0.96 * (y - y0) / dy),
    );
    let mut order: Vec<usize> = (0..emb.points.len()).collect();
    if let Colouring::Umi(m) = colouring {
        order.sort_by_key(|&i| m.get(&emb.points[i].0).copied().unwrap_or(0));
    }
    let (lo, hi) = match colouring {
        Colouring::Umi(m) => {
            let v: Vec<f64> = m.values().map(|&u| (u.max(1) as f64).log10()).collect();
            (
                v.iter().copied().fold(f64::MAX, f64::min),
                v.iter().copied().fold(f64::MIN, f64::max),
            )
        }
        Colouring::Cluster(_) => (0.0, 1.0),
    };
    let mut s = f.open(label);
    let _ = write!(
        s,
        "<rect class=\"frame\" x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>",
        f.ml,
        f.mt,
        f.w - f.ml - f.mr,
        f.h - f.mt - f.mb
    );
    for i in order {
        let (bc, x, y) = &emb.points[i];
        match colouring {
            Colouring::Cluster(m) => match m.get(bc) {
                Some(k) => {
                    let _ = write!(
                        s,
                        "<circle class=\"k{}\" cx=\"{:.1}\" cy=\"{:.1}\" r=\"2.4\"/>",
                        (k.saturating_sub(1)) % 12,
                        px(*x),
                        py(*y)
                    );
                }
                None => {
                    let _ = write!(
                        s,
                        "<circle class=\"kx\" cx=\"{:.1}\" cy=\"{:.1}\" r=\"2.4\"/>",
                        px(*x),
                        py(*y)
                    );
                }
            },
            Colouring::Umi(m) => {
                let u = m.get(bc).copied().unwrap_or(1).max(1) as f64;
                let t = if hi > lo {
                    (u.log10() - lo) / (hi - lo)
                } else {
                    0.5
                };
                let _ = write!(
                    s,
                    "<circle fill=\"{}\" cx=\"{:.1}\" cy=\"{:.1}\" r=\"2.4\"/>",
                    ramp(t),
                    px(*x),
                    py(*y)
                );
            }
        }
    }
    let _ = write!(
        s,
        "<text class=\"axt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{} 1</text><text class=\"axt\" transform=\"translate(10 {:.1}) rotate(-90)\" text-anchor=\"middle\">{} 2</text></svg>",
        f.ml + (f.w - f.ml - f.mr) / 2.0,
        f.h - 8.0,
        esc(emb.name),
        f.mt + (f.h - f.mt - f.mb) / 2.0,
        esc(emb.name)
    );
    s
}

const CSS: &str = r#"
:root{--bg:#f6f7f9;--card:#fff;--ink:#1c2430;--mute:#5b6676;--line:#d9dee6;--faint:#b9c1cd;--accent:#1f6feb;--accent2:#d9480f;
--k0:#1f77b4;--k1:#e8710a;--k2:#2ca02c;--k3:#d62728;--k4:#8e5bc4;--k5:#8c564b;--k6:#e377c2;--k7:#6b7280;--k8:#b5a800;--k9:#17a2b8;--k10:#0b6e4f;--k11:#c2185b;--kx:#c3c9d3}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#0f141b;--card:#171e28;--ink:#e6ebf2;--mute:#99a5b6;--line:#2b3544;--faint:#4a5668;--accent:#58a6ff;--accent2:#ff9b5e;
--k0:#5aa9e6;--k1:#ffa94d;--k2:#6fcf6f;--k3:#ff6b6b;--k4:#b794f4;--k5:#c49a8c;--k6:#f28ac8;--k7:#a0a8b5;--k8:#e5d84a;--k9:#4dd0e1;--k10:#46c79e;--k11:#f06292;--kx:#3a4556}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1080px;margin:0 auto;padding:20px 16px 48px}
header{margin-bottom:20px}
h1{font-size:1.6rem;margin:0 0 4px;letter-spacing:-.01em}
h2{font-size:1.1rem;margin:28px 0 10px}
.sub{color:var(--mute);margin:0}
.meta{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:8px 20px;margin-top:14px}
.meta div{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:8px 12px}
.meta b{display:block;font-size:.72rem;font-weight:600;color:var(--mute);text-transform:uppercase;letter-spacing:.05em}
.hero{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:12px}
.hero div{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px}
.hero .v{font-size:2rem;font-weight:700;color:var(--accent);line-height:1.1;font-variant-numeric:tabular-nums}
.hero .k{color:var(--mute);font-size:.85rem}
.grp{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:12px 16px}
.card h3{font-size:.8rem;margin:0 0 6px;color:var(--mute);text-transform:uppercase;letter-spacing:.05em}
table{border-collapse:collapse;width:100%}
td,th{padding:5px 0;border-bottom:1px solid var(--line);text-align:left;font-variant-numeric:tabular-nums}
tr:last-child td{border-bottom:0}
td.n,th.n{text-align:right}
.plots{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,340px),1fr));gap:12px}
.plots svg{width:100%;height:auto;display:block}
.cap{color:var(--mute);font-size:.85rem;margin:4px 0 0}
svg .frame{fill:none;stroke:var(--line)}
svg .grid{stroke:var(--line);stroke-width:.6;stroke-dasharray:2 3}
svg .tk{fill:var(--mute);font-size:10px}
svg .axt{fill:var(--mute);font-size:11px}
svg .rk-cell{fill:none;stroke:var(--accent);stroke-width:2.6;stroke-linejoin:round}
svg .rk-bg{fill:none;stroke:var(--faint);stroke-width:2.6;stroke-linejoin:round}
svg .line{fill:none;stroke:var(--accent);stroke-width:2}
svg .dot{fill:var(--accent)}
.k0{fill:var(--k0)}.k1{fill:var(--k1)}.k2{fill:var(--k2)}.k3{fill:var(--k3)}.k4{fill:var(--k4)}.k5{fill:var(--k5)}.k6{fill:var(--k6)}.k7{fill:var(--k7)}.k8{fill:var(--k8)}.k9{fill:var(--k9)}.k10{fill:var(--k10)}.k11{fill:var(--k11)}.kx{fill:var(--kx)}
.legend{display:flex;flex-wrap:wrap;gap:4px 14px;margin:8px 0 0;font-size:.85rem;color:var(--mute)}
.legend i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px;vertical-align:baseline}
.legend .ramp{display:inline-block;width:90px;height:8px;border-radius:4px;vertical-align:middle;margin:0 6px}
.sw-cell{background:var(--accent)}.sw-bg{background:var(--faint)}
.k0i{background:var(--k0)}.k1i{background:var(--k1)}.k2i{background:var(--k2)}.k3i{background:var(--k3)}.k4i{background:var(--k4)}.k5i{background:var(--k5)}.k6i{background:var(--k6)}.k7i{background:var(--k7)}.k8i{background:var(--k8)}.k9i{background:var(--k9)}.k10i{background:var(--k10)}.k11i{background:var(--k11)}
.scroll{overflow-x:auto}
footer{margin-top:32px;color:var(--mute);font-size:.8rem}
@media (max-width:520px){.hero .v{font-size:1.6rem}h1{font-size:1.3rem}}
"#;

fn metric_group(title: &str, metrics: &[(String, String)], names: &[&str]) -> String {
    let mut s = format!("<div class=\"card\"><h3>{}</h3><table>", esc(title));
    for n in names {
        let _ = write!(
            s,
            "<tr><td>{}</td><td class=\"n\">{}</td></tr>",
            esc(n),
            esc(metric(metrics, n))
        );
    }
    s.push_str("</table></div>");
    s
}

/// Render the whole page. Analysis plots appear when `outs_dir/analysis` exists.
pub fn render_web_summary(
    outs_dir: &Path,
    metrics: &[(String, String)],
    info: &RunInfo,
    ranks: &RankData,
    curve: Option<&[CurvePoint]>,
) -> String {
    let mut h = String::with_capacity(1 << 20);
    let title = if info.sample_id.is_empty() {
        "rustar run summary".to_string()
    } else {
        format!("{} | rustar run summary", info.sample_id)
    };
    let _ = write!(
        h,
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"color-scheme\" content=\"light dark\"><title>{}</title><style>{}</style></head><body><main>",
        esc(&title),
        CSS
    );
    let _ = write!(
        h,
        "<header><h1>Run summary</h1><p class=\"sub\">Single-cell gene expression counts from rustar-aligner</p><div class=\"meta\"><div><b>Sample</b>{}</div><div><b>Chemistry</b>{}</div><div><b>Reference</b>{}</div><div><b>rustar version</b>{}</div></div></header>",
        esc(or_na(&info.sample_id)),
        esc(or_na(&info.chemistry)),
        esc(or_na(&info.reference)),
        esc(or_na(&info.version))
    );

    h.push_str("<h2>Key metrics</h2><div class=\"hero\">");
    for (k, label) in [
        ("Estimated Number of Cells", "Estimated cells"),
        ("Mean Reads per Cell", "Mean reads per cell"),
        ("Median Genes per Cell", "Median genes per cell"),
        ("Median UMI Counts per Cell", "Median UMI counts per cell"),
    ] {
        let _ = write!(
            h,
            "<div><div class=\"v\">{}</div><div class=\"k\">{}</div></div>",
            esc(metric(metrics, k)),
            label
        );
    }
    h.push_str("</div><h2>All metrics</h2><div class=\"grp\">");
    h.push_str(&metric_group(
        "Sequencing",
        metrics,
        &[
            "Number of Reads",
            "Valid Barcodes",
            "Valid UMI Sequences",
            "Sequencing Saturation",
            "Q30 Bases in Barcode",
            "Q30 Bases in RNA Read",
            "Q30 Bases in UMI",
        ],
    ));
    h.push_str(&metric_group(
        "Mapping",
        metrics,
        &[
            "Reads Mapped to Genome",
            "Reads Mapped Confidently to Genome",
            "Reads Mapped Confidently to Intergenic Regions",
            "Reads Mapped Confidently to Intronic Regions",
            "Reads Mapped Confidently to Exonic Regions",
            "Reads Mapped Confidently to Transcriptome",
            "Reads Mapped Antisense to Gene",
        ],
    ));
    h.push_str(&metric_group(
        "Cells",
        metrics,
        &[
            "Estimated Number of Cells",
            "Fraction Reads in Cells",
            "Mean Reads per Cell",
            "Median Genes per Cell",
            "Median UMI Counts per Cell",
            "Total Genes Detected",
        ],
    ));
    h.push_str("</div>");

    h.push_str("<h2>Cells</h2><div class=\"plots\">");
    let rp = rank_plot(ranks);
    if !rp.is_empty() {
        let _ = write!(
            h,
            "<div class=\"card\"><h3>Barcode rank plot</h3>{rp}<div class=\"legend\"><span><i class=\"sw-cell\"></i>Cells</span><span><i class=\"sw-bg\"></i>Background</span></div><p class=\"cap\">UMI counts per barcode, barcodes ranked from most to fewest, both axes logarithmic.</p></div>"
        );
    }
    h.push_str("</div>");

    if let Some(c) = curve.filter(|c| c.len() >= 2) {
        h.push_str("<h2>Sequencing depth</h2><div class=\"plots\">");
        let _ = write!(
            h,
            "<div class=\"card\"><h3>Sequencing saturation</h3>{}<p class=\"cap\">Expected value when only a fraction of the reads is kept.</p></div><div class=\"card\"><h3>Median genes per cell</h3>{}<p class=\"cap\">Same subsampling, genes with at least one UMI.</p></div></div>",
            curve_plot(c, "Saturation (%)", |p| p.saturation, Some(100.0)),
            curve_plot(c, "Median genes per cell", |p| p.median_genes, None)
        );
    }

    if let Some(an) = load_analysis(outs_dir) {
        h.push_str("<h2>Secondary analysis</h2>");
        let mut ks: Vec<u32> = an.clusters.values().copied().collect();
        ks.sort_unstable();
        ks.dedup();
        let legend = if ks.is_empty() {
            String::new()
        } else {
            let mut l = String::from("<div class=\"legend\">");
            for k in &ks {
                let _ = write!(
                    l,
                    "<span><i class=\"k{}i\"></i>Cluster {}</span>",
                    k.saturating_sub(1) % 12,
                    k
                );
            }
            l.push_str("</div>");
            l
        };
        let ramp_legend = format!(
            "<div class=\"legend\"><span>Low UMI<span class=\"ramp\" style=\"background:linear-gradient(90deg,{},{},{},{},{})\"></span>High UMI</span></div>",
            ramp(0.0),
            ramp(0.25),
            ramp(0.5),
            ramp(0.75),
            ramp(1.0)
        );
        h.push_str("<div class=\"plots\">");
        for emb in &an.embeddings {
            if !an.clusters.is_empty() {
                let _ = write!(
                    h,
                    "<div class=\"card\"><h3>{} by graph-based cluster</h3>{}{}</div>",
                    emb.name,
                    scatter(
                        emb,
                        &Colouring::Cluster(&an.clusters),
                        &format!("{} coloured by cluster", emb.name)
                    ),
                    legend
                );
            }
            if !ranks.cell_umis.is_empty() {
                let _ = write!(
                    h,
                    "<div class=\"card\"><h3>{} by UMI counts</h3>{}{}</div>",
                    emb.name,
                    scatter(
                        emb,
                        &Colouring::Umi(&ranks.cell_umis),
                        &format!("{} coloured by UMI counts", emb.name)
                    ),
                    ramp_legend
                );
            }
        }
        h.push_str("</div>");
        if an.top_genes.iter().any(|(_, g)| !g.is_empty()) {
            h.push_str("<h2>Top genes per cluster</h2><div class=\"card scroll\"><table><tr><th>Cluster</th><th>Genes (log2 fold change)</th></tr>");
            for (k, genes) in &an.top_genes {
                let list: Vec<String> = genes
                    .iter()
                    .map(|(n, l, _)| format!("{} ({:+.1})", esc(n), l))
                    .collect();
                let _ = write!(
                    h,
                    "<tr><td><i class=\"legend\" style=\"display:inline\"><i class=\"k{}i\"></i></i>{}</td><td>{}</td></tr>",
                    k.saturating_sub(1) % 12,
                    k,
                    list.join(", ")
                );
            }
            h.push_str("</table><p class=\"cap\">Up to five genes per cluster with the largest positive log2 fold change versus all other cells, adjusted p below 0.05.</p></div>");
        }
    }

    let _ = write!(
        h,
        "<footer>Generated by rustar-aligner {}. This page is a single file with no external resources.</footer></main></body></html>",
        esc(&info.version)
    );
    h
}

fn or_na(s: &str) -> &str {
    if s.is_empty() { "n/a" } else { s }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solo::count::Molecule;

    fn mol(cb: u32, gene: u32, umi: u64, count: u32) -> Molecule {
        Molecule {
            cb,
            gene,
            umi,
            count,
            utype: 1,
        }
    }

    #[test]
    fn csv_split_quotes() {
        assert_eq!(split_csv("a,\"1,221\",3%"), vec!["a", "1,221", "3%"]);
    }

    #[test]
    fn curve_is_monotone_and_ends_at_observed() {
        let m = vec![
            mol(0, 0, 1, 3),
            mol(0, 0, 2, 1),
            mol(0, 1, 1, 2),
            mol(1, 0, 1, 1),
        ];
        let c = saturation_curve(&m, &[0, 1], 7);
        assert_eq!(c.len(), 10);
        assert!(
            c.windows(2)
                .all(|w| w[1].saturation >= w[0].saturation - 1e-9)
        );
        // 7 reads, 4 molecules: observed saturation is 1 - 4/7.
        assert!((c[9].saturation - 100.0 * 3.0 / 7.0).abs() < 1e-6);
    }

    #[test]
    fn page_has_no_dashes_or_scripts() {
        let ranks = RankData {
            ranked: (1..200u64).rev().map(|u| (u * 10, u > 150)).collect(),
            cell_umis: HashMap::new(),
        };
        let metrics = vec![("Estimated Number of Cells".to_string(), "49".to_string())];
        let html = render_web_summary(
            Path::new("/nonexistent"),
            &metrics,
            &RunInfo::default(),
            &ranks,
            None,
        );
        assert!(!html.contains('\u{2014}'));
        assert!(!html.contains("<script"));
        assert!(html.contains("rk-cell"));
    }

    /// Build the page from existing run output; run by hand:
    /// `cargo test --release --features hdf5-out web_summary_from_outs -- --ignored`.
    #[test]
    #[cfg(unix)]
    #[ignore = "needs /Users/benjamin/rustar-parity data"]
    fn web_summary_from_outs() {
        let outs = Path::new("/Users/benjamin/rustar-parity/cr/rs7/rs_outs");
        // Stand-in analysis folder: link it in a scratch copy of the outs dir.
        let scratch = Path::new("/Users/benjamin/rustar-parity/cleanroom/web_work/outs");
        std::fs::create_dir_all(scratch).unwrap();
        for f in ["raw_feature_bc_matrix", "filtered_feature_bc_matrix"] {
            let _ = std::fs::remove_file(scratch.join(f));
            std::os::unix::fs::symlink(outs.join(f), scratch.join(f)).unwrap();
        }
        let _ = std::fs::remove_file(scratch.join("analysis"));
        std::os::unix::fs::symlink(
            "/Users/benjamin/rustar-parity/cr/ref/pbmc1k/outs/analysis",
            scratch.join("analysis"),
        )
        .unwrap();
        let info = RunInfo {
            sample_id: "pbmc1k".into(),
            chemistry: "Single Cell 3' v3 (polyA)".into(),
            reference: "GRCh38".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        };
        let metrics = read_metrics_csv(&outs.join("metrics_summary.csv")).unwrap();
        let ranks = read_rank_data(outs).unwrap();
        let html = render_web_summary(scratch, &metrics, &info, &ranks, None);
        let out = Path::new("/Users/benjamin/rustar-parity/cleanroom/web_work/web_summary.html");
        std::fs::write(out, &html).unwrap();
        assert!(html.len() < 2_000_000, "page is {} bytes", html.len());
    }
}
