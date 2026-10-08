//! Rank-one NMF over-approximation, ported from DegNorm's `degnorm/nmf.py`
//! (NUStatBioinfo/DegNorm, GPL-free MIT-compatible research code; the algorithm
//! is described in Xiong et al., *Genome Biology* 2019).
//!
//! Pure math: no I/O, no globals, no randomness.
//!
//! DegNorm calls `scipy.sparse.linalg.svds(x, k=1)` for the leading singular
//! triplet. The number of rows here is the number of RNA-seq libraries (small),
//! so the triplet is obtained by power iteration on the `p x p` Gram matrix
//! `x x^T`, which is cheaper and deterministic.

/// A dense `p x l` matrix in row-major order.
#[derive(Clone, Debug)]
pub struct Mat {
    pub p: usize,
    pub l: usize,
    pub data: Vec<f64>,
}

impl Mat {
    pub fn new(p: usize, l: usize) -> Self {
        Mat {
            p,
            l,
            data: vec![0.0; p * l],
        }
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.l + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.l + j] = v;
    }

    pub fn row_sums(&self) -> Vec<f64> {
        (0..self.p)
            .map(|i| self.data[i * self.l..(i + 1) * self.l].iter().sum())
            .collect()
    }

    /// Sample-wise maximum per transcript position.
    pub fn col_max(&self) -> Vec<f64> {
        (0..self.l)
            .map(|j| {
                (0..self.p)
                    .map(|i| self.get(i, j))
                    .fold(f64::NEG_INFINITY, f64::max)
            })
            .collect()
    }

    pub fn max(&self) -> f64 {
        self.data.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Keep only the columns listed in `idx`, in the order given.
    #[must_use]
    pub fn select_cols(&self, idx: &[usize]) -> Mat {
        let mut out = Mat::new(self.p, idx.len());
        for i in 0..self.p {
            for (jj, &j) in idx.iter().enumerate() {
                out.set(i, jj, self.get(i, j));
            }
        }
        out
    }
}

/// Leading singular triplet as `(k, e)` with `k[i] * e[j] ~= x[i][j]`, signs
/// fixed so both factors are non-negative for a non-negative input.
pub fn rank_one(mat: &Mat) -> (Vec<f64>, Vec<f64>) {
    let n_samples = mat.p;
    if n_samples == 0 || mat.l == 0 {
        return (vec![0.0; n_samples], vec![0.0; mat.l]);
    }

    // Gram matrix `gram = mat mat^T` (n_samples x n_samples).
    let mut gram = vec![0.0f64; n_samples * n_samples];
    for row_a in 0..n_samples {
        for row_b in row_a..n_samples {
            let mut acc = 0.0;
            for j in 0..mat.l {
                acc += mat.get(row_a, j) * mat.get(row_b, j);
            }
            gram[row_a * n_samples + row_b] = acc;
            gram[row_b * n_samples + row_a] = acc;
        }
    }

    // Power iteration for the dominant eigenvector of the Gram matrix.
    let mut u_vec = vec![1.0 / (n_samples as f64).sqrt(); n_samples];
    for _ in 0..1000 {
        let mut next = vec![0.0f64; n_samples];
        for (row, slot) in next.iter_mut().enumerate() {
            *slot = (0..n_samples)
                .map(|col| gram[row * n_samples + col] * u_vec[col])
                .sum();
        }
        let norm = next.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm <= 0.0 {
            return (vec![0.0; n_samples], vec![0.0; mat.l]);
        }
        for v in &mut next {
            *v /= norm;
        }
        let delta: f64 = next.iter().zip(&u_vec).map(|(a, b)| (a - b).abs()).sum();
        u_vec = next;
        if delta < 1e-12 {
            break;
        }
    }

    // Point the dominant direction the positive way.
    if u_vec.iter().sum::<f64>() < 0.0 {
        for v in &mut u_vec {
            *v = -*v;
        }
    }

    // envelope_raw = mat^T u; sigma = ||envelope_raw||; abundance = u * sigma.
    let mut envelope: Vec<f64> = (0..mat.l)
        .map(|j| (0..n_samples).map(|i| mat.get(i, j) * u_vec[i]).sum())
        .collect();
    let sigma = envelope.iter().map(|v| v * v).sum::<f64>().sqrt();
    if sigma <= 0.0 {
        return (vec![0.0; n_samples], vec![0.0; mat.l]);
    }
    for v in &mut envelope {
        *v /= sigma;
    }
    let abundance: Vec<f64> = u_vec.iter().map(|v| v * sigma).collect();
    (abundance, envelope)
}

/// Outer product of the rank-one factors.
pub fn outer(abundance: &[f64], envelope: &[f64]) -> Mat {
    let mut m = Mat::new(abundance.len(), envelope.len());
    for (i, &a) in abundance.iter().enumerate() {
        for (j, &e) in envelope.iter().enumerate() {
            m.set(i, j, a * e);
        }
    }
    m
}

/// NMF over-approximation: faithful port of `GeneNMFOA.nmf` — dual ascent on a
/// non-negative multiplier `lambda` with step `1 / sqrt(iters)`, re-fitting the
/// rank-one factors of `x + lambda` each round.
///
/// As upstream, the returned factors are *not* clamped to dominate `x`; callers
/// that need a strict over-approximation apply [`over_approximate`] (DegNorm
/// does the same, its in-function clamp is commented out).
pub fn nmf_oa(mat: &Mat, iters: usize) -> (Vec<f64>, Vec<f64>) {
    let (mut abundance, mut envelope) = rank_one(mat);
    let mut est = outer(&abundance, &envelope);
    let mut lambda = Mat::new(mat.p, mat.l);
    let step = 1.0 / (iters.max(1) as f64).sqrt();

    let mut shifted = Mat::new(mat.p, mat.l);
    for _ in 0..iters {
        for i in 0..mat.p {
            for j in 0..mat.l {
                let residual = est.get(i, j) - mat.get(i, j);
                let dual = (lambda.get(i, j) - step * residual).max(0.0);
                lambda.set(i, j, dual);
                shifted.set(i, j, mat.get(i, j) + dual);
            }
        }
        let (a2, e2) = rank_one(&shifted);
        abundance = a2;
        envelope = e2;
        est = outer(&abundance, &envelope);
    }
    (abundance, envelope)
}

/// Raise `est` elementwise to at least `f` (DegNorm's over-approximation
/// quality control).
pub fn over_approximate(est: &mut Mat, f: &Mat) {
    for i in 0..f.p {
        for j in 0..f.l {
            if est.get(i, j) < f.get(i, j) {
                est.set(i, j, f.get(i, j));
            }
        }
    }
}

/// One-shot rank-one estimate raised to at least the input: DegNorm's
/// `ratio_svd`, used only to initialise the depth scale factors.
pub fn ratio_svd(mat: &Mat) -> Mat {
    let (abundance, envelope) = rank_one(mat);
    let mut est = outer(&abundance, &envelope);
    over_approximate(&mut est, mat);
    est
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn mat_from(rows: &[&[f64]]) -> Mat {
        let p = rows.len();
        let l = rows[0].len();
        let mut m = Mat::new(p, l);
        for (i, r) in rows.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                m.set(i, j, v);
            }
        }
        m
    }

    #[test]
    fn rank_one_recovers_a_planted_rank_one_matrix() {
        // x = k e^T with k = [1, 2, 3], e = [4, 5].
        let x = mat_from(&[&[4.0, 5.0], &[8.0, 10.0], &[12.0, 15.0]]);
        let (k, e) = rank_one(&x);
        for (i, &ki) in k.iter().enumerate() {
            for (j, &ej) in e.iter().enumerate() {
                assert!((ki * ej - x.get(i, j)).abs() < 1e-8);
            }
        }
        assert!(k.iter().all(|&v| v > 0.0));
        assert!(e.iter().all(|&v| v > 0.0));
    }

    #[test]
    fn rank_one_of_a_zero_matrix_is_zero() {
        let x = Mat::new(2, 3);
        let (k, e) = rank_one(&x);
        assert!(k.iter().all(|&v| v == 0.0));
        assert!(e.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn nmf_oa_moves_the_estimate_above_the_rank_one_fit() {
        // Sample 2 is degraded at the 3' end: a rank-one fit cannot cover it.
        let x = mat_from(&[&[10.0, 10.0, 10.0, 10.0], &[10.0, 10.0, 2.0, 1.0]]);
        let (k0, e0) = rank_one(&x);
        let (k, e) = nmf_oa(&x, 100);
        let plain = outer(&k0, &e0);
        let lifted = outer(&k, &e);
        let deficit = |m: &Mat| -> f64 {
            let mut d = 0.0;
            for i in 0..x.p {
                for j in 0..x.l {
                    d += (x.get(i, j) - m.get(i, j)).max(0.0);
                }
            }
            d
        };
        assert!(
            deficit(&lifted) < deficit(&plain),
            "NMF-OA should reduce the under-approximation deficit: {} vs {}",
            deficit(&lifted),
            deficit(&plain)
        );
    }

    #[test]
    fn ratio_svd_dominates_input() {
        let x = mat_from(&[&[5.0, 1.0], &[1.0, 5.0]]);
        let est = ratio_svd(&x);
        for i in 0..2 {
            for j in 0..2 {
                assert!(est.get(i, j) >= x.get(i, j) - 1e-9);
            }
        }
    }

    #[test]
    fn select_cols_keeps_requested_positions() {
        let x = mat_from(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]]);
        let sub = x.select_cols(&[0, 2]);
        assert_eq!(sub.l, 2);
        assert!((sub.get(1, 1) - 6.0).abs() < 1e-12);
    }
}
