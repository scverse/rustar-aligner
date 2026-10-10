//! Differential expression of every cluster against all other cells.
//!
//! For each feature and cluster the table has the mean count, the log2 fold
//! change and a Benjamini-Hochberg adjusted p value.
//!
//! * **Mean counts** = (counts in the cluster / total counts of the cluster's
//!   cells) x median cell total. This reproduces the reference values to
//!   floating-point precision.
//! * **log2 fold change** = `log2((c_in + 1) / (c_out + 1)) + log2(S_out /
//!   S_in)` with `c` the raw counts and `S` the summed cell totals of each
//!   side (a pseudocount of one count, size-normalised). The reference
//!   values follow the same formula up to a per-cluster constant that is
//!   within 0.01 of this one (0.09 for a 76-cell cluster).
//! * **p values**: negative binomial exact test on the size-factor
//!   normalised, group-summed counts (Robinson and Smyth; the sSeq test of Yu,
//!   Huber and Vitek 2013), per-feature dispersion by the method of moments
//!   on the pooled within-group variance, shrunk towards a high quantile of
//!   the dispersions. The shrinkage constants were tuned by eye against the
//!   reference output rather than derived, so p values are approximate: the
//!   ranking agrees (Spearman 0.75 to 0.92 on log10 adjusted p) and nearly
//!   every feature the reference calls significant is also significant here,
//!   but the number of features passing 0.05 differs by up to a factor 1.5.

use super::{CountMatrix, create, fmt_f64, write_line};
use crate::error::Error;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

/// Quantile of the raw dispersions used as the shrinkage target, and the
/// weight given to it.
const SHRINK_QUANTILE: f64 = 0.875;
const SHRINK_WEIGHT: f64 = 0.8;

pub struct Context<'a> {
    counts: &'a CountMatrix,
    present: &'a [usize],
    pos: Vec<u32>,
    /// Size factors (cell total / mean cell total).
    sf: Vec<f64>,
    sizes: &'a [f64],
    med_size: f64,
    /// Mean of 1 / size factor.
    inv_sf_mean: f64,
}

/// Benjamini-Hochberg adjustment.
pub fn bh(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap());
    let mut out = vec![1.0; n];
    let mut run = 1.0f64;
    for r in (0..n).rev() {
        let i = idx[r];
        run = run.min(p[i] * n as f64 / (r + 1) as f64);
        out[i] = run.min(1.0);
    }
    out
}

/// Two-sided exact test for two negative binomial sums conditioned on their
/// total `s = x1 + x2`: sizes `a` and `b` (replicates / dispersion), common
/// probability. Sums the conditional probabilities not exceeding that of the
/// observed split.
pub fn nb_exact_p(x1: u64, x2: u64, a: f64, b: f64) -> f64 {
    let s = x1 + x2;
    if s == 0 {
        return 1.0;
    }
    let sf = s as f64;
    // w(k) relative to w(x1) = 1.
    let up = |k: f64| (k + a) * (sf - k) / ((k + 1.0) * (sf - k - 1.0 + b));
    let down = |k: f64| k * (sf - k + b) / ((k - 1.0 + a) * (sf - k + 1.0));
    let mut total = 1.0;
    let mut tail = 1.0;
    let eps = 1.0 + 1e-9;
    // Upwards.
    let (mut w, mut peak) = (1.0f64, 1.0f64);
    let mut k = x1;
    while k < s {
        w *= up(k as f64);
        k += 1;
        total += w;
        if w <= eps {
            tail += w;
        }
        peak = peak.max(w);
        if w < 1e-18 * peak && up(k as f64) < 1.0 {
            break;
        }
    }
    // Downwards.
    let (mut w, mut peak) = (1.0f64, 1.0f64);
    let mut k = x1;
    while k > 0 {
        w *= down(k as f64);
        k -= 1;
        total += w;
        if w <= eps {
            tail += w;
        }
        peak = peak.max(w);
        if w < 1e-18 * peak && down(k as f64) < 1.0 {
            break;
        }
    }
    (tail / total).min(1.0)
}

impl<'a> Context<'a> {
    pub fn new(counts: &'a CountMatrix, sizes: &'a [f64], present: &'a [usize]) -> Self {
        let mut pos = vec![u32::MAX; counts.n_genes()];
        for (i, &g) in present.iter().enumerate() {
            pos[g] = i as u32;
        }
        let mean = sizes.iter().sum::<f64>() / sizes.len() as f64;
        let sf: Vec<f64> = sizes.iter().map(|s| s / mean).collect();
        let mut sorted = sizes.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = sorted.len();
        let med_size = if n % 2 == 1 {
            sorted[n / 2]
        } else {
            0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
        };
        let inv_sf_mean = sf.iter().map(|s| 1.0 / s).sum::<f64>() / n as f64;
        Self {
            counts,
            present,
            pos,
            sf,
            sizes,
            med_size,
            inv_sf_mean,
        }
    }

    /// Write `differential_expression.csv` for `labels` (1-based clusters).
    pub fn write(&self, labels: &[u32], path: &Path) -> Result<(), Error> {
        let kc = labels.iter().copied().max().unwrap_or(0) as usize;
        let ng = self.present.len();
        if kc < 2 || ng == 0 {
            return Ok(());
        }
        let n = labels.len();
        let mut raw = vec![0.0; ng * kc];
        let mut s1 = vec![0.0; ng * kc];
        let mut s2 = vec![0.0; ng * kc];
        let mut ncell = vec![0usize; kc];
        let mut ssize = vec![0.0; kc];
        for j in 0..n {
            let c = labels[j] as usize - 1;
            ncell[c] += 1;
            ssize[c] += self.sizes[j];
            for p in self.counts.indptr[j]..self.counts.indptr[j + 1] {
                let q = self.pos[self.counts.rows[p] as usize];
                if q == u32::MAX {
                    continue;
                }
                let x = self.counts.vals[p];
                let nx = x / self.sf[j];
                let i = q as usize * kc + c;
                raw[i] += x;
                s1[i] += nx;
                s2[i] += nx * nx;
            }
        }
        let tot_n: usize = ncell.iter().sum();
        let tot_size: f64 = ssize.iter().sum();

        // Per cluster: mean, log2fc, p.
        let mut mean_counts = vec![0.0; ng * kc];
        let mut l2fc = vec![0.0; ng * kc];
        let mut adj = vec![1.0; ng * kc];
        for c in 0..kc {
            let (n1, n2) = (ncell[c], tot_n - ncell[c]);
            if n1 == 0 || n2 == 0 {
                continue;
            }
            let (sz_in, sz_out) = (ssize[c], tot_size - ssize[c]);
            let tot_raw = |g: usize| (0..kc).map(|d| raw[g * kc + d]).sum::<f64>();
            let tot_s1 = |g: usize| (0..kc).map(|d| s1[g * kc + d]).sum::<f64>();
            let tot_s2 = |g: usize| (0..kc).map(|d| s2[g * kc + d]).sum::<f64>();
            // Pass 1: summaries and the method-of-moments dispersion.
            let pre: Vec<(f64, f64, f64, f64, f64)> = (0..ng)
                .into_par_iter()
                .map(|g| {
                    let i = g * kc + c;
                    let r_in = raw[i];
                    let r_out = tot_raw(g) - r_in;
                    let mean = r_in / sz_in * self.med_size;
                    let fc = ((r_in + 1.0) / (r_out + 1.0)).log2() + (sz_out / sz_in).log2();
                    let (x1, x2) = (s1[i], tot_s1(g) - s1[i]);
                    let (q1, q2) = (s2[i], tot_s2(g) - s2[i]);
                    let (f1, f2) = (n1 as f64, n2 as f64);
                    let (m1, m2) = (x1 / f1, x2 / f2);
                    let v1 = if n1 > 1 {
                        ((q1 - f1 * m1 * m1) / (f1 - 1.0)).max(0.0)
                    } else {
                        0.0
                    };
                    let v2 = if n2 > 1 {
                        ((q2 - f2 * m2 * m2) / (f2 - 1.0)).max(0.0)
                    } else {
                        0.0
                    };
                    let dfp = (f1 + f2 - 2.0).max(1.0);
                    let vp = ((f1 - 1.0).max(0.0) * v1 + (f2 - 1.0).max(0.0) * v2) / dfp;
                    let mm = ((f1 - 1.0).max(0.0) * m1 + (f2 - 1.0).max(0.0) * m2) / dfp;
                    let phi = if mm > 0.0 {
                        ((vp - mm * self.inv_sf_mean) / (mm * mm)).max(0.0)
                    } else {
                        f64::NAN
                    };
                    (mean, fc, x1, x2, phi)
                })
                .collect();
            // Shrink the dispersions towards a high quantile of their
            // distribution (sSeq-style shrinkage, parameters tuned against
            // the reference output).
            let mut phis: Vec<f64> = pre.iter().map(|r| r.4).filter(|x| x.is_finite()).collect();
            phis.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let target = if phis.is_empty() {
                0.0
            } else {
                phis[((phis.len() - 1) as f64 * SHRINK_QUANTILE) as usize]
            };
            let (f1, f2) = (n1 as f64, n2 as f64);
            let res: Vec<(f64, f64, f64)> = pre
                .par_iter()
                .map(|&(mean, fc, x1, x2, phi)| {
                    let phi = if phi.is_finite() {
                        (1.0 - SHRINK_WEIGHT) * phi + SHRINK_WEIGHT * target
                    } else {
                        target
                    }
                    .max(1e-3);
                    let p = nb_exact_p(x1.round() as u64, x2.round() as u64, f1 / phi, f2 / phi);
                    (mean, fc, p)
                })
                .collect();
            let pv: Vec<f64> = res.iter().map(|r| r.2).collect();
            let padj = bh(&pv);
            for g in 0..ng {
                mean_counts[g * kc + c] = res[g].0;
                l2fc[g * kc + c] = res[g].1;
                adj[g * kc + c] = padj[g];
            }
        }

        let mut w = create(path)?;
        let io = |e| Error::io(e, path);
        let mut head = String::from("Feature ID,Feature Name");
        for c in 1..=kc {
            head.push_str(&format!(
                ",Cluster {c} Mean Counts,Cluster {c} Log2 fold change,Cluster {c} Adjusted p value"
            ));
        }
        write_line(&mut w, &head).map_err(io)?;
        for (g, &fi) in self.present.iter().enumerate() {
            let name = &self.counts.feature_names[fi];
            let name = if name.contains(',') || name.contains('"') {
                format!("\"{}\"", name.replace('"', "\"\""))
            } else {
                name.clone()
            };
            let mut line = format!("{},{}", self.counts.feature_ids[fi], name);
            for c in 0..kc {
                line.push_str(&format!(
                    ",{},{},{}",
                    fmt_f64(mean_counts[g * kc + c]),
                    fmt_f64(l2fc[g * kc + c]),
                    fmt_f64(adj[g * kc + c])
                ));
            }
            write_line(&mut w, &line).map_err(io)?;
        }
        w.flush().map_err(io)
    }
}
