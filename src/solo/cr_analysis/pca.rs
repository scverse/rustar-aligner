//! Feature selection by normalised dispersion and PCA.
//!
//! Method (public documentation + published techniques):
//!
//! 1. Library-size normalisation: every cell is scaled to the median total
//!    count of the cells.
//! 2. Features with fewer than 3 counts in total are dropped. The remaining
//!    features are all used for the PCA (the reference output lists them all
//!    in `features_selected.csv`), ordered by increasing normalised
//!    dispersion.
//! 3. Dispersion = variance / mean of the normalised counts. Features are
//!    binned by log mean expression (20 equal-width bins) and the dispersion
//!    is robustly standardised in each bin: `(d - median) / MAD`.
//! 4. PCA input: `ln(1 + normalised count)`, centred and scaled to unit
//!    variance per feature. The top components come from a randomised block
//!    subspace iteration on the implicit (sparse + rank one) matrix, refined
//!    with a Rayleigh-Ritz step (a truncated SVD, as IRLBA computes).
//!
//! What is exact and what is approximate versus the reference run: the set
//! of features and the PCA construction reproduce the reference projections
//! to ~1e-3 relative (equal up to sign); the dispersion values follow the
//! reference closely (rank correlation ~0.98) but the exact binning of the
//! reference is not known, so they are not identical.

use super::{CountMatrix, create, fmt_f64, rng::Rng, write_line};
use crate::error::Error;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

/// Minimum total count for a feature to take part.
const MIN_TOTAL_COUNT: f64 = 3.0;
const N_DISP_BINS: usize = 20;

#[derive(Debug)]
pub struct PcaResult {
    n_cells: usize,
    /// Feature indices that took part, in matrix order.
    pub selected: Vec<usize>,
    /// Normalised dispersion, one entry per feature (NaN when not selected).
    pub dispersion: Vec<f64>,
    /// Indices into `selected`-feature space of the dispersion ordering
    /// (ascending), stored as feature indices.
    pub order: Vec<usize>,
    /// `components[pc][i]` is the loading of `selected[i]`.
    pub components: Vec<Vec<f64>>,
    /// `scores[pc][cell]`.
    pub scores: Vec<Vec<f64>>,
    /// Proportion of total variance per component.
    pub variance: Vec<f64>,
}

impl PcaResult {
    pub fn n_components(&self) -> usize {
        self.scores.len()
    }
    /// The projection as a row-major `n_cells x n_components` matrix.
    pub fn projection_rows(&self) -> Vec<f64> {
        let k = self.n_components();
        let mut out = vec![0.0; self.n_cells * k];
        for (c, s) in self.scores.iter().enumerate() {
            for (j, v) in s.iter().enumerate() {
                out[j * k + c] = *v;
            }
        }
        out
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        f64::NAN
    } else if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Dispersion standardised inside mean-expression bins. `sel` are feature
/// indices; returns a vector over all `n_genes` (NaN elsewhere).
fn normalized_dispersion(
    counts: &CountMatrix,
    sizes: &[f64],
    med_size: f64,
    sel: &[usize],
    pos: &[u32],
) -> Vec<f64> {
    let n = counts.n_cells() as f64;
    let mut s1 = vec![0.0; sel.len()];
    let mut s2 = vec![0.0; sel.len()];
    for j in 0..counts.n_cells() {
        let f = med_size / sizes[j];
        for p in counts.indptr[j]..counts.indptr[j + 1] {
            let q = pos[counts.rows[p] as usize];
            if q != u32::MAX {
                let x = counts.vals[p] * f;
                s1[q as usize] += x;
                s2[q as usize] += x * x;
            }
        }
    }
    let mut mean = vec![0.0; sel.len()];
    let mut disp = vec![0.0; sel.len()];
    for i in 0..sel.len() {
        mean[i] = s1[i] / n;
        let var = ((s2[i] - n * mean[i] * mean[i]) / (n - 1.0)).max(0.0);
        disp[i] = var / mean[i];
    }
    let lm: Vec<f64> = mean.iter().map(|m| m.ln()).collect();
    let (lo, hi) = lm
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
            (a.min(x), b.max(x))
        });
    let width = (hi - lo).max(f64::MIN_POSITIVE);
    let bin_of = |x: f64| (((x - lo) / width * N_DISP_BINS as f64) as usize).min(N_DISP_BINS - 1);
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); N_DISP_BINS];
    for (i, &x) in lm.iter().enumerate() {
        bins[bin_of(x)].push(i);
    }
    let mut norm = vec![0.0; sel.len()];
    for members in &bins {
        if members.is_empty() {
            continue;
        }
        let mut d: Vec<f64> = members.iter().map(|&i| disp[i]).collect();
        let m = median(&mut d);
        let mut dev: Vec<f64> = members.iter().map(|&i| (disp[i] - m).abs()).collect();
        let mad = median(&mut dev);
        for &i in members {
            norm[i] = if mad > 0.0 { (disp[i] - m) / mad } else { 0.0 };
        }
    }
    let mut out = vec![f64::NAN; counts.n_genes()];
    for (i, &g) in sel.iter().enumerate() {
        out[g] = norm[i];
    }
    out
}

/// Symmetric eigendecomposition by cyclic Jacobi. `a` is `n x n` row-major;
/// returns eigenvalues (descending) and eigenvectors as columns of a
/// row-major matrix.
fn jacobi_eigen(mut a: Vec<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j] * a[i * n + j])
            .sum();
        let diag: f64 = (0..n).map(|i| a[i * n + i] * a[i * n + i]).sum();
        if off <= 1e-30 * diag.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&x, &y| a[y * n + y].partial_cmp(&a[x * n + x]).unwrap());
    let vals = idx.iter().map(|&i| a[i * n + i]).collect();
    let mut vecs = vec![0.0; n * n];
    for (newc, &oldc) in idx.iter().enumerate() {
        for r in 0..n {
            vecs[r * n + newc] = v[r * n + oldc];
        }
    }
    (vals, vecs)
}

/// Orthonormalise the columns of a row-major `rows x b` matrix (modified
/// Gram-Schmidt, twice).
fn orthonormalize(m: &mut [f64], rows: usize, b: usize) {
    for _ in 0..2 {
        for c in 0..b {
            for p in 0..c {
                let dot: f64 = (0..rows).map(|r| m[r * b + p] * m[r * b + c]).sum();
                for r in 0..rows {
                    m[r * b + c] -= dot * m[r * b + p];
                }
            }
            let nrm = (0..rows)
                .map(|r| m[r * b + c] * m[r * b + c])
                .sum::<f64>()
                .sqrt();
            if nrm > 1e-300 {
                for r in 0..rows {
                    m[r * b + c] /= nrm;
                }
            }
        }
    }
}

/// The implicit matrix `X` (cells x features): `ln(1+x)` values, centred and
/// scaled per feature.
struct Implicit {
    n_cells: usize,
    n_sel: usize,
    // CSC by cell over selected features.
    c_ptr: Vec<usize>,
    c_idx: Vec<u32>,
    c_val: Vec<f64>,
    // CSR by feature.
    g_ptr: Vec<usize>,
    g_idx: Vec<u32>,
    g_val: Vec<f64>,
    mu: Vec<f64>,
    inv_sd: Vec<f64>,
}

impl Implicit {
    /// `X V` for `V` (features x b) -> cells x b.
    fn times(&self, v: &[f64], b: usize) -> Vec<f64> {
        let mut w = vec![0.0; self.n_sel * b];
        let mut c = vec![0.0; b];
        for g in 0..self.n_sel {
            for k in 0..b {
                let x = v[g * b + k] * self.inv_sd[g];
                w[g * b + k] = x;
                c[k] += self.mu[g] * x;
            }
        }
        let mut y = vec![0.0; self.n_cells * b];
        y.par_chunks_mut(b).enumerate().for_each(|(j, row)| {
            for (k, r) in row.iter_mut().enumerate() {
                *r = -c[k];
            }
            for p in self.c_ptr[j]..self.c_ptr[j + 1] {
                let g = self.c_idx[p] as usize;
                let val = self.c_val[p];
                for k in 0..b {
                    row[k] += val * w[g * b + k];
                }
            }
        });
        y
    }

    /// `X^T Q` for `Q` (cells x b) -> features x b.
    fn t_times(&self, q: &[f64], b: usize) -> Vec<f64> {
        let mut qs = vec![0.0; b];
        for j in 0..self.n_cells {
            for k in 0..b {
                qs[k] += q[j * b + k];
            }
        }
        let mut p = vec![0.0; self.n_sel * b];
        p.par_chunks_mut(b).enumerate().for_each(|(g, row)| {
            for q_i in self.g_ptr[g]..self.g_ptr[g + 1] {
                let j = self.g_idx[q_i] as usize;
                let val = self.g_val[q_i];
                for k in 0..b {
                    row[k] += val * q[j * b + k];
                }
            }
            for (k, r) in row.iter_mut().enumerate() {
                *r = (*r - self.mu[g] * qs[k]) * self.inv_sd[g];
            }
        });
        p
    }
}

pub fn run(counts: &CountMatrix, sizes: &[f64], npc: usize) -> Result<PcaResult, Error> {
    let n_cells = counts.n_cells();
    let totals = counts.gene_totals();
    let selected: Vec<usize> = (0..counts.n_genes())
        .filter(|&g| totals[g] >= MIN_TOTAL_COUNT)
        .collect();
    let mut pos = vec![u32::MAX; counts.n_genes()];
    for (i, &g) in selected.iter().enumerate() {
        pos[g] = i as u32;
    }
    let mut ssorted = sizes.to_vec();
    let med_size = median(&mut ssorted);
    let dispersion = if selected.is_empty() {
        vec![f64::NAN; counts.n_genes()]
    } else {
        normalized_dispersion(counts, sizes, med_size, &selected, &pos)
    };
    let mut order = selected.clone();
    order.sort_by(|&a, &b| dispersion[a].partial_cmp(&dispersion[b]).unwrap());

    let k = npc.min(n_cells.saturating_sub(1)).min(selected.len());
    if k == 0 {
        return Ok(PcaResult {
            n_cells,
            selected,
            dispersion,
            order,
            components: Vec::new(),
            scores: Vec::new(),
            variance: Vec::new(),
        });
    }

    // ln(1 + x) values in both orientations.
    let n_sel = selected.len();
    let mut c_ptr = vec![0usize; n_cells + 1];
    let mut c_idx = Vec::new();
    let mut c_val = Vec::new();
    for j in 0..n_cells {
        let f = med_size / sizes[j];
        for p in counts.indptr[j]..counts.indptr[j + 1] {
            let q = pos[counts.rows[p] as usize];
            if q != u32::MAX {
                c_idx.push(q);
                c_val.push((counts.vals[p] * f).ln_1p());
            }
        }
        c_ptr[j + 1] = c_idx.len();
    }
    let mut g_ptr = vec![0usize; n_sel + 1];
    for &q in &c_idx {
        g_ptr[q as usize + 1] += 1;
    }
    for g in 0..n_sel {
        g_ptr[g + 1] += g_ptr[g];
    }
    let mut fill = g_ptr.clone();
    let mut g_idx = vec![0u32; c_idx.len()];
    let mut g_val = vec![0.0; c_idx.len()];
    for j in 0..n_cells {
        for p in c_ptr[j]..c_ptr[j + 1] {
            let g = c_idx[p] as usize;
            g_idx[fill[g]] = j as u32;
            g_val[fill[g]] = c_val[p];
            fill[g] += 1;
        }
    }
    let n = n_cells as f64;
    let mut mu = vec![0.0; n_sel];
    let mut inv_sd = vec![0.0; n_sel];
    let mut n_var = 0usize;
    for g in 0..n_sel {
        let (s, e) = (g_ptr[g], g_ptr[g + 1]);
        let s1: f64 = g_val[s..e].iter().sum();
        let s2: f64 = g_val[s..e].iter().map(|x| x * x).sum();
        mu[g] = s1 / n;
        let var = (s2 / n - mu[g] * mu[g]).max(0.0);
        if var > 1e-24 {
            inv_sd[g] = 1.0 / var.sqrt();
            n_var += 1;
        }
    }
    let x = Implicit {
        n_cells,
        n_sel,
        c_ptr,
        c_idx,
        c_val,
        g_ptr,
        g_idx,
        g_val,
        mu,
        inv_sd,
    };

    let b = (k + 10).min(n_cells).min(n_sel);
    let mut rng = Rng::new(0x5EED_0001);
    let mut q: Vec<f64> = (0..n_cells * b).map(|_| rng.f64() - 0.5).collect();
    orthonormalize(&mut q, n_cells, b);
    let mut prev = vec![0.0; k];
    let mut p;
    let mut iter = 0;
    let (lam, evec) = loop {
        p = x.t_times(&q, b);
        // C = P^T P (b x b)
        let mut c = vec![0.0; b * b];
        for r in 0..n_sel {
            for i in 0..b {
                let pi = p[r * b + i];
                for j in i..b {
                    c[i * b + j] += pi * p[r * b + j];
                }
            }
        }
        for i in 0..b {
            for j in 0..i {
                c[i * b + j] = c[j * b + i];
            }
        }
        let (lam, evec) = jacobi_eigen(c, b);
        iter += 1;
        let delta = (0..k)
            .map(|i| (lam[i] - prev[i]).abs() / lam[i].abs().max(1e-300))
            .fold(0.0, f64::max);
        prev[..k].copy_from_slice(&lam[..k]);
        if (iter > 2 && delta < 1e-13) || iter >= 300 {
            break (lam, evec);
        }
        orthonormalize(&mut p, n_sel, b);
        q = x.times(&p, b);
        orthonormalize(&mut q, n_cells, b);
    };

    let total_ss = n * n_var as f64;
    let mut components = Vec::with_capacity(k);
    let mut scores = Vec::with_capacity(k);
    let mut variance = Vec::with_capacity(k);
    for i in 0..k {
        let sigma = lam[i].max(0.0).sqrt();
        let mut load: Vec<f64> = (0..n_sel)
            .map(|g| {
                (0..b).map(|c| p[g * b + c] * evec[c * b + i]).sum::<f64>() / sigma.max(1e-300)
            })
            .collect();
        let mut sc: Vec<f64> = (0..n_cells)
            .map(|j| (0..b).map(|c| q[j * b + c] * evec[c * b + i]).sum::<f64>() * sigma)
            .collect();
        // Deterministic sign: the largest-magnitude loading is negative.
        let big = load
            .iter()
            .cloned()
            .fold(0.0f64, |a, x| if x.abs() > a.abs() { x } else { a });
        if big > 0.0 {
            load.iter_mut().for_each(|x| *x = -*x);
            sc.iter_mut().for_each(|x| *x = -*x);
        }
        components.push(load);
        scores.push(sc);
        variance.push(lam[i] / total_ss);
    }
    Ok(PcaResult {
        n_cells,
        selected,
        dispersion,
        order,
        components,
        scores,
        variance,
    })
}

/// Write the five `pca/` CSV files. `present` are the features with any count
/// (the reference lists those, and only those, in `components.csv` and
/// `dispersion.csv`).
pub fn write(
    counts: &CountMatrix,
    pca: &PcaResult,
    present: &[usize],
    dir: &Path,
) -> Result<(), Error> {
    let k = pca.n_components();
    let mut spos = vec![usize::MAX; counts.n_genes()];
    for (i, &g) in pca.selected.iter().enumerate() {
        spos[g] = i;
    }

    let path = dir.join("components.csv");
    let mut w = create(&path)?;
    let mut head = String::from("PC");
    for &g in present {
        head.push(',');
        head.push_str(&counts.feature_ids[g]);
    }
    let io = |e| Error::io(e, &path);
    write_line(&mut w, &head).map_err(io)?;
    for (c, comp) in pca.components.iter().enumerate() {
        let mut line = format!("{}", c + 1);
        for &g in present {
            line.push(',');
            line.push_str(&fmt_f64(if spos[g] == usize::MAX {
                0.0
            } else {
                comp[spos[g]]
            }));
        }
        write_line(&mut w, &line).map_err(io)?;
    }
    w.flush().map_err(io)?;

    let path = dir.join("dispersion.csv");
    let mut w = create(&path)?;
    let io = |e| Error::io(e, &path);
    write_line(&mut w, "Feature,Normalized.Dispersion").map_err(io)?;
    for &g in present {
        write_line(
            &mut w,
            &format!("{},{}", counts.feature_ids[g], fmt_f64(pca.dispersion[g])),
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)?;

    let path = dir.join("features_selected.csv");
    let mut w = create(&path)?;
    let io = |e| Error::io(e, &path);
    write_line(&mut w, "Feature").map_err(io)?;
    for (i, &g) in pca.order.iter().enumerate() {
        write_line(&mut w, &format!("{},{}", i + 1, counts.feature_ids[g])).map_err(io)?;
    }
    w.flush().map_err(io)?;

    let path = dir.join("projection.csv");
    let mut w = create(&path)?;
    let io = |e| Error::io(e, &path);
    let mut head = String::from("Barcode");
    for c in 0..k {
        head.push_str(&format!(",PC-{}", c + 1));
    }
    write_line(&mut w, &head).map_err(io)?;
    for (j, bc) in counts.barcodes.iter().enumerate() {
        let mut line = bc.clone();
        for c in 0..k {
            line.push(',');
            line.push_str(&fmt_f64(pca.scores[c][j]));
        }
        write_line(&mut w, &line).map_err(io)?;
    }
    w.flush().map_err(io)?;

    let path = dir.join("variance.csv");
    let mut w = create(&path)?;
    let io = |e| Error::io(e, &path);
    write_line(&mut w, "PC,Proportion.Variance.Explained").map_err(io)?;
    for (c, v) in pca.variance.iter().enumerate() {
        write_line(&mut w, &format!("{},{}", c + 1, fmt_f64(*v))).map_err(io)?;
    }
    w.flush().map_err(io)?;
    Ok(())
}
