//! DegNorm driver: load per-sample coverage, run the NMF-OA iterations, and
//! write DI scores plus degradation-adjusted counts.
//!
//! Port of `GeneNMFOA.run` in DegNorm's `degnorm/nmf.py`.

use std::io::Write;
use std::path::Path;

use rayon::prelude::*;

use crate::degnorm::baseline::{BaselineParams, baseline_selection};
use crate::degnorm::nmf::{Mat, ratio_svd};
use crate::error::Error;
use crate::quant::coverage::CoverageFile;

pub struct DegNormConfig {
    pub iter: usize,
    pub nmf_iter: usize,
    pub downsample_rate: usize,
    pub minimax_coverage: u32,
    pub skip_baseline: bool,
    pub bins: usize,
    pub min_high_coverage: usize,
    pub seed: u64,
}

pub struct DegNormOutput {
    pub gene_ids: Vec<String>,
    pub sample_ids: Vec<String>,
    /// DI score per included gene, per sample.
    pub rho: Vec<Vec<f64>>,
    pub raw_counts: Vec<Vec<f64>>,
    pub adjusted_counts: Vec<Vec<f64>>,
    pub scale_factors: Vec<f64>,
    /// Genes sent through baseline selection, per outer iteration.
    pub n_baseline_selected: Vec<usize>,
}

fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        f64::midpoint(s[n / 2 - 1], s[n / 2])
    }
}

/// Run the full DegNorm pipeline over `files` (one per sample).
pub fn run_degnorm(files: &[CoverageFile], cfg: &DegNormConfig) -> Result<DegNormOutput, Error> {
    let p = files.len();
    if p < 2 {
        return Err(Error::Parameter(
            "--runMode degNorm requires at least 2 coverage files: a degradation index is \
             defined across samples, not within one"
                .to_string(),
        ));
    }

    // All samples must share one gene model.
    for f in &files[1..] {
        if f.gene_ids != files[0].gene_ids || f.gene_lens != files[0].gene_lens {
            let first_diff = f
                .gene_ids
                .iter()
                .zip(&files[0].gene_ids)
                .position(|(a, b)| a != b)
                .map_or_else(
                    || "<gene length mismatch>".to_string(),
                    |i| f.gene_ids[i].clone(),
                );
            return Err(Error::Parameter(format!(
                "coverage file '{}' has a different gene model than '{}' (first difference: {}); \
                 all samples must be aligned against the same GTF",
                f.sample_id, files[0].sample_id, first_diff
            )));
        }
    }

    let sample_ids: Vec<String> = files.iter().map(|f| f.sample_id.clone()).collect();

    // Gene selection: every sample must reach the minimax coverage threshold.
    let min_cov = cfg.minimax_coverage.max(1);
    let kept: Vec<usize> = (0..files[0].n_genes())
        .filter(|&g| {
            files[0].gene_lens[g] >= 2
                && files
                    .iter()
                    .all(|f| f.gene(g).iter().copied().max().unwrap_or(0) >= min_cov)
        })
        .collect();
    if kept.is_empty() {
        return Err(Error::Parameter(
            "no gene reached --degNormMinimaxCoverage in every sample; nothing to normalize"
                .to_string(),
        ));
    }
    let n_genes = kept.len();

    // Coverage matrices, one per kept gene, shape p x L.
    let mats: Vec<Mat> = kept
        .iter()
        .map(|&g| {
            let l = files[0].gene_lens[g] as usize;
            let mut m = Mat::new(p, l);
            for (i, f) in files.iter().enumerate() {
                for (j, &v) in f.gene(g).iter().enumerate() {
                    m.set(i, j, f64::from(v));
                }
            }
            m
        })
        .collect();

    // Raw count matrix (genes x samples).
    let x: Vec<Vec<f64>> = kept
        .iter()
        .map(|&g| files.iter().map(|f| f64::from(f.counts[g])).collect())
        .collect();

    // Initialisation: one-shot over-approximation gives the starting DI scores.
    let mut rho: Vec<Vec<f64>> = mats
        .par_iter()
        .map(|m| {
            let est = ratio_svd(m);
            let fs = m.row_sums();
            let es = est.row_sums();
            (0..m.p)
                .map(|i| (1.0 - fs[i] / (es[i] + 1.0)).clamp(0.0, 0.9))
                .collect()
        })
        .collect();

    let low_di: Vec<usize> = (0..n_genes)
        .filter(|&g| rho[g].iter().copied().fold(f64::NEG_INFINITY, f64::max) < 0.1)
        .collect();
    let count_sums: Vec<f64> = (0..p)
        .map(|i| {
            if low_di.is_empty() {
                (0..n_genes).map(|g| x[g][i]).sum()
            } else {
                low_di.iter().map(|&g| x[g][i]).sum()
            }
        })
        .collect();
    let med = median(&count_sums);
    let mut norm_factors: Vec<f64> = count_sums
        .iter()
        .map(|&c| if med > 0.0 { (c / med).max(1e-12) } else { 1.0 })
        .collect();
    let mut x_weighted: Vec<Vec<f64>> = (0..n_genes)
        .map(|g| (0..p).map(|i| x[g][i] / norm_factors[i]).collect())
        .collect();
    let mut scale_factors = norm_factors.clone();
    let mut x_adj = x_weighted.clone();
    let mut n_baseline_selected = Vec::with_capacity(cfg.iter);

    let bp = BaselineParams {
        nmf_iter: cfg.nmf_iter,
        bins: cfg.bins,
        min_high_coverage: cfg.min_high_coverage,
        downsample_rate: cfg.downsample_rate,
        skip: cfg.skip_baseline,
    };

    for it in 0..cfg.iter {
        // Scale each sample's coverage curve by its depth factor, then fit.
        let results: Vec<(Vec<f64>, bool)> = mats
            .par_iter()
            .enumerate()
            .map(|(gi, m)| {
                let mut adj = Mat::new(m.p, m.l);
                for (i, factor) in scale_factors.iter().enumerate().take(m.p) {
                    let s = factor.max(1e-12);
                    for j in 0..m.l {
                        adj.set(i, j, m.get(i, j) / s);
                    }
                }
                let r = baseline_selection(&adj, &bp, cfg.seed.wrapping_add(gi as u64));
                (r.rho, r.ran_baseline)
            })
            .collect();

        let n_sel = results.iter().filter(|(_, b)| *b).count();
        n_baseline_selected.push(n_sel);
        log::info!(
            "degNorm iteration {}: {n_sel} genes through baseline selection",
            it + 1
        );
        for (g, (r, _)) in results.into_iter().enumerate() {
            rho[g] = r;
        }

        // Genes that never went through baseline selection inherit the
        // sample-average DI (DegNorm's `correct_di_scores`).
        x_adj = adjust_counts(&x_weighted, &rho, n_genes, p);
        let sample_avg_di: Vec<f64> = (0..p)
            .map(|i| {
                let w: f64 = (0..n_genes).map(|g| x_weighted[g][i]).sum();
                let a: f64 = (0..n_genes).map(|g| x_adj[g][i]).sum();
                if a > 0.0 { 1.0 - w / a } else { 0.0 }
            })
            .collect();
        for r in rho.iter_mut().take(n_genes) {
            if r.iter().copied().fold(f64::NEG_INFINITY, f64::max) <= 0.0 {
                r.clone_from(&sample_avg_di);
            }
        }

        x_adj = adjust_counts(&x_weighted, &rho, n_genes, p);

        // Fold the degradation correction back into the depth factors.
        let col_sums: Vec<f64> = (0..p)
            .map(|i| (0..n_genes).map(|g| x_adj[g][i]).sum())
            .collect();
        let med = median(&col_sums);
        norm_factors = col_sums
            .iter()
            .map(|&c| if med > 0.0 { (c / med).max(1e-12) } else { 1.0 })
            .collect();
        for row in x_weighted.iter_mut().take(n_genes) {
            for (v, nf) in row.iter_mut().zip(&norm_factors) {
                *v /= *nf;
            }
        }
        for (s, nf) in scale_factors.iter_mut().zip(&norm_factors) {
            *s *= *nf;
        }
    }

    Ok(DegNormOutput {
        gene_ids: kept.iter().map(|&g| files[0].gene_ids[g].clone()).collect(),
        sample_ids,
        rho,
        raw_counts: x,
        adjusted_counts: x_adj,
        scale_factors,
        n_baseline_selected,
    })
}

fn adjust_counts(
    x_weighted: &[Vec<f64>],
    rho: &[Vec<f64>],
    n_genes: usize,
    p: usize,
) -> Vec<Vec<f64>> {
    (0..n_genes)
        .map(|g| {
            (0..p)
                .map(|i| x_weighted[g][i] / (1.0 - rho[g][i]).max(1e-12))
                .collect()
        })
        .collect()
}

fn write_matrix(
    path: &Path,
    row_ids: &[String],
    header: &[String],
    rows: &[Vec<f64>],
) -> Result<(), Error> {
    let mut f = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
    let io = |e: std::io::Error| Error::io(e, path);
    writeln!(f, "gene\t{}", header.join("\t")).map_err(io)?;
    for (g, id) in row_ids.iter().enumerate() {
        let cells: Vec<String> = rows[g].iter().map(|v| format!("{v:.6}")).collect();
        writeln!(f, "{id}\t{}", cells.join("\t")).map_err(io)?;
    }
    Ok(())
}

/// Write `DegradationIndex.tab`, `AdjustedCounts.tab`, `RawCounts.tab`,
/// `ScaleFactors.tab`, and `Summary.txt` into `dir`.
pub fn write_outputs(out: &DegNormOutput, dir: &Path, cfg: &DegNormConfig) -> Result<(), Error> {
    std::fs::create_dir_all(dir).map_err(|e| Error::io(e, dir))?;
    write_matrix(
        &dir.join("DegradationIndex.tab"),
        &out.gene_ids,
        &out.sample_ids,
        &out.rho,
    )?;
    write_matrix(
        &dir.join("AdjustedCounts.tab"),
        &out.gene_ids,
        &out.sample_ids,
        &out.adjusted_counts,
    )?;
    write_matrix(
        &dir.join("RawCounts.tab"),
        &out.gene_ids,
        &out.sample_ids,
        &out.raw_counts,
    )?;

    let sf_path = dir.join("ScaleFactors.tab");
    let mut f = std::fs::File::create(&sf_path).map_err(|e| Error::io(e, &sf_path))?;
    {
        let io = |e: std::io::Error| Error::io(e, sf_path.as_path());
        writeln!(f, "sample\tscaleFactor").map_err(io)?;
        for (i, s) in out.sample_ids.iter().enumerate() {
            writeln!(f, "{s}\t{:.6}", out.scale_factors[i]).map_err(io)?;
        }
    }

    let sum_path = dir.join("Summary.txt");
    let mut f = std::fs::File::create(&sum_path).map_err(|e| Error::io(e, &sum_path))?;
    let io = |e: std::io::Error| Error::io(e, sum_path.as_path());
    writeln!(f, "samples\t{}", out.sample_ids.len()).map_err(io)?;
    writeln!(f, "genesIncluded\t{}", out.gene_ids.len()).map_err(io)?;
    writeln!(f, "degNormIter\t{}", cfg.iter).map_err(io)?;
    writeln!(f, "degNormNmfIter\t{}", cfg.nmf_iter).map_err(io)?;
    writeln!(f, "downsampleRate\t{}", cfg.downsample_rate).map_err(io)?;
    writeln!(f, "minimaxCoverage\t{}", cfg.minimax_coverage).map_err(io)?;
    writeln!(f, "skipBaselineSelection\t{}", cfg.skip_baseline).map_err(io)?;
    for (i, n) in out.n_baseline_selected.iter().enumerate() {
        writeln!(f, "baselineSelectedIter{}\t{}", i + 1, n).map_err(io)?;
    }
    for (i, s) in out.sample_ids.iter().enumerate() {
        let col: Vec<f64> = out.rho.iter().map(|r| r[i]).collect();
        writeln!(f, "medianDI_{}\t{:.6}", s, median(&col)).map_err(io)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an in-memory `CoverageFile` with `n_genes` genes of length `len`,
    /// coverage supplied by `f(gene, position)`.
    fn synth(
        sample: &str,
        n_genes: usize,
        len: usize,
        f: impl Fn(usize, usize) -> u32,
    ) -> CoverageFile {
        let mut cov = Vec::new();
        let mut counts = Vec::new();
        let mut offsets = vec![0u64];
        for g in 0..n_genes {
            let mut total = 0u64;
            for j in 0..len {
                let v = f(g, j);
                total += u64::from(v);
                cov.push(v);
            }
            counts.push((total / 100).max(1) as u32);
            offsets.push(offsets[g] + len as u64);
        }
        CoverageFile {
            sample_id: sample.to_string(),
            paired: false,
            n_counted: counts.iter().map(|&c| u64::from(c)).sum(),
            gene_ids: (0..n_genes).map(|g| format!("G{g}")).collect(),
            gene_lens: vec![len as u32; n_genes],
            counts,
            offsets,
            cov,
        }
    }

    fn cfg() -> DegNormConfig {
        DegNormConfig {
            iter: 2,
            nmf_iter: 50,
            downsample_rate: 1,
            minimax_coverage: 0,
            skip_baseline: false,
            bins: 20,
            min_high_coverage: 50,
            seed: 777,
        }
    }

    #[test]
    fn degraded_sample_gets_higher_di_and_upweighted_counts() {
        let len = 400;
        let a = synth("ctrl", 2, len, |_, _| 20);
        let b = synth(
            "degraded",
            2,
            len,
            |g, j| {
                if g == 0 && j >= 160 { 2 } else { 20 }
            },
        );
        let out = run_degnorm(&[a, b], &cfg()).unwrap();

        assert_eq!(
            out.sample_ids,
            vec!["ctrl".to_string(), "degraded".to_string()]
        );
        assert_eq!(out.gene_ids.len(), 2);
        assert!(
            out.rho[0][1] > out.rho[0][0] + 0.2,
            "DI of the degraded sample ({}) should exceed the control ({})",
            out.rho[0][1],
            out.rho[0][0]
        );
        assert!(out.adjusted_counts[0][1] > out.raw_counts[0][1] * 0.5);
    }

    #[test]
    fn mismatched_gene_sets_are_rejected() {
        let a = synth("a", 2, 400, |_, _| 10);
        let mut b = synth("b", 2, 400, |_, _| 10);
        b.gene_ids[1] = "OTHER".to_string();
        assert!(run_degnorm(&[a, b], &cfg()).is_err());
    }

    #[test]
    fn single_sample_is_rejected() {
        let a = synth("a", 2, 400, |_, _| 10);
        assert!(run_degnorm(&[a], &cfg()).is_err());
    }

    #[test]
    fn outputs_are_written() {
        let a = synth("a", 2, 400, |_, _| 20);
        let b = synth("b", 2, 400, |_, j| if j >= 200 { 3 } else { 20 });
        let out = run_degnorm(&[a, b], &cfg()).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("DegNorm.out");
        write_outputs(&out, &target, &cfg()).unwrap();
        for name in [
            "DegradationIndex.tab",
            "AdjustedCounts.tab",
            "RawCounts.tab",
            "ScaleFactors.tab",
            "Summary.txt",
        ] {
            assert!(target.join(name).exists(), "{name} missing");
        }
        let di = std::fs::read_to_string(target.join("DegradationIndex.tab")).unwrap();
        assert!(di.starts_with("gene\ta\tb\n"));
        assert_eq!(di.lines().count(), 3);
    }
}
