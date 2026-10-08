//! Baseline selection: DegNorm's per-gene search for a transcript region where
//! degradation is minimal, so the coverage envelope (and therefore the DI
//! scores) can be estimated from undegraded positions.
//!
//! Port of `GeneNMFOA.baseline_selection` in DegNorm's `degnorm/nmf.py`, quirks
//! included: the `+ 1` in the DI denominator, the 0.1 / 0.2 / 0.9 thresholds,
//! and the bin-dropping loop.

use std::collections::HashSet;

use crate::degnorm::nmf::{Mat, nmf_oa, outer, over_approximate};

pub struct BaselineParams {
    pub nmf_iter: usize,
    pub bins: usize,
    pub min_high_coverage: usize,
    pub downsample_rate: usize,
    pub skip: bool,
}

pub struct BaselineResult {
    /// DI score per sample for this gene.
    pub rho: Vec<f64>,
    /// Whether the bin-dropping search actually ran.
    pub ran_baseline: bool,
}

/// Columns whose sample-wise maximum exceeds 10% of the matrix maximum.
fn high_coverage_idx(f: &Mat) -> Vec<usize> {
    let thresh = 0.1 * f.max();
    let cm = f.col_max();
    (0..f.l).filter(|&j| cm[j] > thresh).collect()
}

/// Deterministic systematic sample. DegNorm draws the start offset at random;
/// here it is derived from the caller's seed and the gene index so results do
/// not depend on thread scheduling.
fn systematic_sample(n: usize, take_every: usize, seed_offset: u64) -> Vec<usize> {
    if take_every <= 1 || take_every >= n {
        return (0..n).collect();
    }
    let start = (seed_offset % take_every as u64) as usize;
    (start..n).step_by(take_every).collect()
}

/// Split `0..n` into consecutive chunks of `ceil(n / bins)` (DegNorm's
/// `split_into_chunks`, which can yield fewer than `bins` chunks).
fn split_into_chunks(n: usize, bins: usize) -> Vec<Vec<usize>> {
    let csize = n.div_ceil(bins.max(1)).max(1);
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        out.push((i..(i + csize).min(n)).collect());
        i += csize;
    }
    out
}

/// DegNorm's DI formula: `1 - rowsum(F) / (rowsum(estimate) + 1)`.
fn di_scores(f: &Mat, est: &Mat) -> Vec<f64> {
    let fs = f.row_sums();
    let es = est.row_sums();
    (0..f.p).map(|i| 1.0 - fs[i] / (es[i] + 1.0)).collect()
}

fn vmax(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn vmin(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::INFINITY, f64::min)
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

fn clamp_rho(rho: &mut [f64]) {
    for r in rho.iter_mut() {
        *r = r.clamp(0.0, 0.9);
    }
}

/// Fit gene `f` (shape `p x L`, already depth-scaled) and return its DI scores.
pub fn baseline_selection(f: &Mat, params: &BaselineParams, seed_offset: u64) -> BaselineResult {
    let zero = || BaselineResult {
        rho: vec![0.0; f.p],
        ran_baseline: false,
    };

    let mut hi = high_coverage_idx(f);
    if params.downsample_rate > 1 {
        let sampled: HashSet<usize> = systematic_sample(f.l, params.downsample_rate, seed_offset)
            .into_iter()
            .collect();
        hi.retain(|j| sampled.contains(j));
    }
    if hi.len() < params.min_high_coverage.max(2) {
        return zero();
    }

    let f_start = f.select_cols(&hi);
    if f_start.row_sums().iter().any(|&s| s <= 0.0) {
        return zero();
    }

    let (k_start, e_start) = nmf_oa(&f_start, params.nmf_iter);
    let mut f_bin = f_start.clone();
    let (mut k, _e) = (k_start.clone(), e_start.clone());
    let mut ke = outer(&k_start, &e_start);
    let mut rho = di_scores(&f_bin, &ke);

    // Exclude extreme cases where the fit did not converge (upstream check).
    let one_minus: Vec<f64> = rho.iter().map(|r| 1.0 - r).collect();
    if median(&one_minus) > 1.0 {
        return zero();
    }

    let min_gene_len = (200.0 / params.downsample_rate as f64).ceil().max(2.0) as usize;
    let min_bins = (params.bins as f64 * 0.2).ceil() as usize;
    let mut ran_baseline = false;

    let can_run = hi.len() >= min_gene_len && vmin(&rho) <= 0.2 && !params.skip;

    if can_run {
        let mut bin_segs = split_into_chunks(f_bin.l, params.bins);

        while vmax(&rho) > 0.1 {
            ran_baseline = true;

            // Per-column worst squared relative residual, averaged per bin.
            let res: Vec<f64> = (0..f_bin.l)
                .map(|j| {
                    (0..f_bin.p)
                        .map(|i| {
                            let d = (ke.get(i, j) - f_bin.get(i, j)) / (f_bin.get(i, j) + 1.0);
                            d * d
                        })
                        .fold(f64::NEG_INFINITY, f64::max)
                })
                .collect();
            let ss_r: Vec<f64> = bin_segs
                .iter()
                .map(|b| b.iter().map(|&j| res[j]).sum::<f64>() / b.len() as f64)
                .collect();
            if vmax(&ss_r) <= 0.0 {
                break;
            }

            let drop = ss_r
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map_or(0, |(i, _)| i);

            // Drop that bin's columns, then renumber the surviving bins so they
            // index the shrunken matrix (DegNorm's `shift_bins`).
            let dropped: HashSet<usize> = bin_segs[drop].iter().copied().collect();
            let keep: Vec<usize> = (0..f_bin.l).filter(|j| !dropped.contains(j)).collect();
            f_bin = f_bin.select_cols(&keep);
            bin_segs.remove(drop);
            let sizes: Vec<usize> = bin_segs.iter().map(Vec::len).collect();
            bin_segs = {
                let mut out = Vec::with_capacity(sizes.len());
                let mut i = 0;
                for s in sizes {
                    out.push((i..i + s).collect::<Vec<usize>>());
                    i += s;
                }
                out
            };

            if f_bin.l == 0 || bin_segs.is_empty() {
                break;
            }

            let (k2, e2) = nmf_oa(&f_bin, params.nmf_iter);
            k = k2;
            ke = outer(&k, &e2);
            if ke.row_sums().iter().any(|&s| s <= 0.0) {
                break;
            }
            over_approximate(&mut ke, &f_bin);
            rho = di_scores(&f_bin, &ke);

            if bin_segs.len() <= min_bins || f_bin.l < min_gene_len {
                break;
            }
        }

        if vmax(&rho) < 0.2 {
            // A baseline region was found: reuse its per-sample abundances to
            // re-derive the envelope over the whole high-coverage transcript.
            let mut kk: Vec<f64> = k.iter().map(|v| v.abs()).collect();
            let min_pos = kk
                .iter()
                .copied()
                .filter(|&v| v >= 1e-5)
                .fold(f64::INFINITY, f64::min);
            let floor = if min_pos.is_finite() { min_pos } else { 1e-5 };
            for v in &mut kk {
                if *v < 1e-5 {
                    *v = floor;
                }
            }
            let ee: Vec<f64> = (0..f_start.l)
                .map(|j| {
                    (0..f_start.p)
                        .map(|i| f_start.get(i, j) / kk[i])
                        .fold(f64::NEG_INFINITY, f64::max)
                })
                .collect();
            rho = di_scores(&f_start, &outer(&kk, &ee));

            // Long, shallow genes can produce implausibly high DI this way;
            // upstream falls back to the plain fit.
            if vmax(&rho) > 0.9 {
                let mut est = outer(&k_start, &e_start);
                over_approximate(&mut est, &f_start);
                rho = di_scores(&f_start, &est);
            }
        } else {
            let mut est = outer(&k_start, &e_start);
            over_approximate(&mut est, &f_start);
            rho = di_scores(&f_start, &est);
        }
    } else {
        // No baseline search: the DI comes from the plain over-approximation.
        let mut est = outer(&k_start, &e_start);
        over_approximate(&mut est, &f_start);
        rho = di_scores(&f_start, &est);
    }

    clamp_rho(&mut rho);
    BaselineResult { rho, ran_baseline }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> BaselineParams {
        BaselineParams {
            nmf_iter: 100,
            bins: 20,
            min_high_coverage: 50,
            downsample_rate: 1,
            skip: false,
        }
    }

    #[test]
    fn undegraded_gene_has_near_zero_di() {
        let l = 400;
        let mut f = Mat::new(2, l);
        for j in 0..l {
            f.set(0, j, 20.0);
            f.set(1, j, 20.0);
        }
        let out = baseline_selection(&f, &params(), 0);
        assert!(
            out.rho.iter().all(|&r| r.abs() < 0.05),
            "rho = {:?}",
            out.rho
        );
    }

    #[test]
    fn degraded_sample_gets_higher_di_than_control() {
        let l = 400;
        let mut f = Mat::new(2, l);
        for j in 0..l {
            f.set(0, j, 20.0);
            f.set(1, j, if j < 160 { 20.0 } else { 2.0 });
        }
        let out = baseline_selection(&f, &params(), 0);
        assert!(
            out.rho[1] > out.rho[0] + 0.2,
            "degraded DI {} should exceed control DI {}",
            out.rho[1],
            out.rho[0]
        );
        assert!(out.rho.iter().all(|&r| (0.0..=0.9).contains(&r)));
    }

    #[test]
    fn short_gene_returns_zero_di() {
        let mut f = Mat::new(2, 10);
        for j in 0..10 {
            f.set(0, j, 5.0);
            f.set(1, j, 5.0);
        }
        let out = baseline_selection(&f, &params(), 0);
        assert_eq!(out.rho, vec![0.0, 0.0]);
        assert!(!out.ran_baseline);
    }

    #[test]
    fn split_into_chunks_covers_every_index() {
        let chunks = split_into_chunks(97, 20);
        let flat: Vec<usize> = chunks.iter().flatten().copied().collect();
        assert_eq!(flat, (0..97).collect::<Vec<_>>());
    }
}
