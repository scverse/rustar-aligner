//! 2-D embeddings (t-SNE and UMAP) of a PCA projection, plus the CSV writers
//! for the `analysis/tsne` and `analysis/umap` projection files.
//!
//! Independent implementation from public documentation and published
//! methods: exact-gradient t-SNE (van der Maaten and Hinton 2008, with the
//! usual early exaggeration, momentum and adaptive gains) and UMAP
//! (McInnes, Healy and Melville 2018: smooth kNN distances, fuzzy simplicial
//! set union, spectral initialisation, SGD with negative sampling).
//! Parameters: t-SNE perplexity 30, 1000 iterations; UMAP n_neighbors 30,
//! min_dist 0.3, correlation metric. Both are deterministic (fixed seed).
//! The exact O(n^2) t-SNE gradient is intended for a few thousand cells.

#![allow(
    clippy::many_single_char_names,
    clippy::manual_midpoint,
    clippy::explicit_iter_loop,
    clippy::needless_for_each,
    clippy::needless_range_loop,
    clippy::similar_names
)]

use std::io::Write;
use std::path::Path;

const TSNE_PERPLEXITY: f64 = 30.0;
const TSNE_ITERS: usize = 1000;
const UMAP_NEIGHBORS: usize = 30;
const UMAP_MIN_DIST: f64 = 0.3;
const SEED: u64 = 0x5EED_1234_ABCD_0001;

/// SplitMix64 generator (deterministic, dependency free).
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

fn sq_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// First `k` principal axes projections of `x` (power iteration), used as a
/// deterministic t-SNE initialisation.
fn pca_init(x: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let n = x.len();
    let d = x[0].len();
    let mut mean = vec![0.0; d];
    for r in x {
        for (m, v) in mean.iter_mut().zip(r) {
            *m += v / n as f64;
        }
    }
    let c: Vec<Vec<f64>> = x
        .iter()
        .map(|r| r.iter().zip(&mean).map(|(v, m)| v - m).collect())
        .collect();
    let mut axes: Vec<Vec<f64>> = Vec::new();
    let mut rng = Rng(SEED);
    for _ in 0..k {
        let mut v: Vec<f64> = (0..d).map(|_| rng.uniform() - 0.5).collect();
        for _ in 0..200 {
            for a in &axes {
                let dot: f64 = v.iter().zip(a).map(|(p, q)| p * q).sum();
                v.iter_mut().zip(a).for_each(|(p, q)| *p -= dot * q);
            }
            let mut w = vec![0.0; d];
            for r in &c {
                let s: f64 = r.iter().zip(&v).map(|(p, q)| p * q).sum();
                w.iter_mut().zip(r).for_each(|(o, p)| *o += s * p);
            }
            let nrm = w.iter().map(|t| t * t).sum::<f64>().sqrt().max(1e-300);
            v = w.iter().map(|t| t / nrm).collect();
        }
        axes.push(v);
    }
    c.iter()
        .map(|r| {
            axes.iter()
                .map(|a| r.iter().zip(a).map(|(p, q)| p * q).sum())
                .collect()
        })
        .collect()
}

/// 2-D t-SNE embedding of the rows of `x` (perplexity 30, 1000 iterations).
pub fn tsne_2d(x: &[Vec<f64>]) -> Vec<[f64; 2]> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    if n < 4 || x[0].is_empty() {
        return (0..n).map(|i| [i as f64, 0.0]).collect();
    }
    let perp = TSNE_PERPLEXITY.min((n as f64 - 1.0) / 3.0).max(1.0);
    let d2: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| sq_dist(&x[i], &x[j])).collect())
        .collect();
    // Conditional probabilities by per-point binary search on the precision.
    let target = perp.ln();
    let mut p = vec![0.0f64; n * n];
    for i in 0..n {
        let (mut lo, mut hi, mut beta) = (0.0f64, f64::INFINITY, 1.0f64);
        let mut row = vec![0.0f64; n];
        for _ in 0..200 {
            let mut sum = 0.0;
            let mut sdp = 0.0;
            for j in 0..n {
                if j == i {
                    row[j] = 0.0;
                    continue;
                }
                let v = (-beta * d2[i][j]).exp();
                row[j] = v;
                sum += v;
                sdp += d2[i][j] * v;
            }
            let sum = sum.max(1e-300);
            let h = sum.ln() + beta * sdp / sum;
            let diff = h - target;
            if diff.abs() < 1e-5 {
                break;
            }
            if diff > 0.0 {
                lo = beta;
                beta = if hi.is_finite() {
                    (beta + hi) / 2.0
                } else {
                    beta * 2.0
                };
            } else {
                hi = beta;
                beta = (beta + lo) / 2.0;
            }
        }
        let s: f64 = row.iter().sum::<f64>().max(1e-300);
        for j in 0..n {
            p[i * n + j] = row[j] / s;
        }
    }
    // Symmetrise.
    let mut ps = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            ps[i * n + j] = ((p[i * n + j] + p[j * n + i]) / (2.0 * n as f64)).max(1e-12);
        }
    }
    // Initialisation: leading two PCs scaled to std 1e-4.
    let init = pca_init(x, 2);
    let sd = (init.iter().map(|r| r[0] * r[0]).sum::<f64>() / n as f64)
        .sqrt()
        .max(1e-300);
    let mut y: Vec<[f64; 2]> = init
        .iter()
        .map(|r| [r[0] / sd * 1e-4, r[1] / sd * 1e-4])
        .collect();
    let mut vel = vec![[0.0f64; 2]; n];
    let mut gains = vec![[1.0f64; 2]; n];
    let lr = (n as f64 / 12.0).max(200.0);
    let mut num = vec![0.0f64; n * n];
    for it in 0..TSNE_ITERS {
        let exag = if it < 250 { 12.0 } else { 1.0 };
        let mom = if it < 250 { 0.5 } else { 0.8 };
        let mut z = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = y[i][0] - y[j][0];
                let dy = y[i][1] - y[j][1];
                let q = 1.0 / (1.0 + dx * dx + dy * dy);
                num[i * n + j] = q;
                num[j * n + i] = q;
                z += 2.0 * q;
            }
        }
        let z = z.max(1e-300);
        let mut grad = vec![[0.0f64; 2]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let q = num[i * n + j];
                let m = (exag * ps[i * n + j] - q / z) * q;
                grad[i][0] += 4.0 * m * (y[i][0] - y[j][0]);
                grad[i][1] += 4.0 * m * (y[i][1] - y[j][1]);
            }
        }
        for i in 0..n {
            for c in 0..2 {
                let same = (grad[i][c] > 0.0) == (vel[i][c] > 0.0);
                gains[i][c] = if same {
                    (gains[i][c] * 0.8).max(0.01)
                } else {
                    gains[i][c] + 0.2
                };
                vel[i][c] = mom * vel[i][c] - lr * gains[i][c] * grad[i][c];
                y[i][c] += vel[i][c];
            }
        }
        let (mx, my) = y.iter().fold((0.0, 0.0), |a, r| {
            (a.0 + r[0] / n as f64, a.1 + r[1] / n as f64)
        });
        y.iter_mut().for_each(|r| {
            r[0] -= mx;
            r[1] -= my;
        });
    }
    y
}

/// Fit a, b of 1 / (1 + a d^(2b)) to the UMAP target curve (spread 1).
fn fit_ab(min_dist: f64) -> (f64, f64) {
    let xs: Vec<f64> = (1..300).map(|i| i as f64 * 3.0 / 300.0).collect();
    let ys: Vec<f64> = xs
        .iter()
        .map(|&x| {
            if x <= min_dist {
                1.0
            } else {
                (-(x - min_dist)).exp()
            }
        })
        .collect();
    let (mut a, mut b) = (1.5f64, 0.9f64);
    let cost = |a: f64, b: f64| -> f64 {
        xs.iter()
            .zip(&ys)
            .map(|(&x, &y)| {
                let f = 1.0 / (1.0 + a * x.powf(2.0 * b));
                (f - y) * (f - y)
            })
            .sum()
    };
    let mut lambda = 1e-3;
    let mut c = cost(a, b);
    for _ in 0..200 {
        let (mut jaa, mut jab, mut jbb, mut ga, mut gb) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (&x, &y) in xs.iter().zip(&ys) {
            let t = x.powf(2.0 * b);
            let f = 1.0 / (1.0 + a * t);
            let r = f - y;
            let da = -t * f * f;
            let db = -a * t * (2.0 * x.ln()) * f * f;
            jaa += da * da;
            jab += da * db;
            jbb += db * db;
            ga += da * r;
            gb += db * r;
        }
        let (haa, hbb) = (jaa * (1.0 + lambda), jbb * (1.0 + lambda));
        let det = haa * hbb - jab * jab;
        if det.abs() < 1e-300 {
            break;
        }
        let na = a - (hbb * ga - jab * gb) / det;
        let nb = b - (haa * gb - jab * ga) / det;
        let nc = cost(na, nb);
        if nc < c && na > 0.0 && nb > 0.0 {
            a = na;
            b = nb;
            c = nc;
            lambda *= 0.5;
        } else {
            lambda *= 4.0;
            if lambda > 1e8 {
                break;
            }
        }
    }
    (a, b)
}

/// Leading non-trivial eigenvectors of the normalised graph Laplacian by
/// subspace iteration (sparse symmetric weights), scaled to [0, 10].
fn spectral_init(n: usize, edges: &[(usize, usize, f64)], rng: &mut Rng) -> Vec<[f64; 2]> {
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(i, j, w) in edges {
        adj[i].push((j, w));
    }
    let deg: Vec<f64> = adj
        .iter()
        .map(|a| a.iter().map(|e| e.1).sum::<f64>().max(1e-12))
        .collect();
    let k = 3usize; // trivial + 2
    let mut basis: Vec<Vec<f64>> = (0..k)
        .map(|_| (0..n).map(|_| rng.uniform() - 0.5).collect())
        .collect();
    let orth = |b: &mut Vec<Vec<f64>>| {
        for a in 0..b.len() {
            for c in 0..a {
                let dot: f64 = b[a].iter().zip(&b[c]).map(|(p, q)| p * q).sum();
                let (l, r) = b.split_at_mut(a);
                r[0].iter_mut().zip(&l[c]).for_each(|(p, q)| *p -= dot * q);
            }
            let nrm = b[a].iter().map(|t| t * t).sum::<f64>().sqrt().max(1e-300);
            b[a].iter_mut().for_each(|t| *t /= nrm);
        }
    };
    orth(&mut basis);
    for _ in 0..1500 {
        for v in basis.iter_mut() {
            // (I + D^-1/2 W D^-1/2) v
            let mut out = v.clone();
            for i in 0..n {
                let mut s = 0.0;
                for &(j, w) in &adj[i] {
                    s += w * v[j] / (deg[i] * deg[j]).sqrt();
                }
                out[i] += s;
            }
            *v = out;
        }
        orth(&mut basis);
    }
    // basis[0] ~ sqrt(deg) (trivial); take basis[1], basis[2] mapped by D^-1/2.
    let mut res: Vec<[f64; 2]> = (0..n)
        .map(|i| [basis[1][i] / deg[i].sqrt(), basis[2][i] / deg[i].sqrt()])
        .collect();
    for c in 0..2 {
        let lo = res.iter().map(|r| r[c]).fold(f64::INFINITY, f64::min);
        let hi = res.iter().map(|r| r[c]).fold(f64::NEG_INFINITY, f64::max);
        let span = (hi - lo).max(1e-12);
        for r in res.iter_mut() {
            r[c] = 10.0 * (r[c] - lo) / span;
        }
    }
    res
}

/// 2-D UMAP embedding of the rows of `x` (n_neighbors 30, min_dist 0.3,
/// correlation metric).
pub fn umap_2d(x: &[Vec<f64>]) -> Vec<[f64; 2]> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    if n < 4 || x[0].is_empty() {
        return (0..n).map(|i| [i as f64, 0.0]).collect();
    }
    let kn = UMAP_NEIGHBORS.min(n - 1);
    // Correlation distance = cosine distance of row-centred vectors.
    let z: Vec<Vec<f64>> = x
        .iter()
        .map(|r| {
            let m = r.iter().sum::<f64>() / r.len() as f64;
            let c: Vec<f64> = r.iter().map(|v| v - m).collect();
            let nrm = c.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
            c.iter().map(|v| v / nrm).collect()
        })
        .collect();
    let mut knn: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut d: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let dot: f64 = z[i].iter().zip(&z[j]).map(|(p, q)| p * q).sum();
                (j, (1.0 - dot).max(0.0))
            })
            .collect();
        d.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        d.truncate(kn);
        knn.push(d);
    }
    // Smooth kNN distances -> directed membership strengths.
    let target = (kn as f64).log2();
    let mut w = vec![std::collections::HashMap::<usize, f64>::new(); n];
    for i in 0..n {
        let rho = knn[i].iter().map(|e| e.1).find(|&d| d > 0.0).unwrap_or(0.0);
        let (mut lo, mut hi, mut sigma) = (0.0f64, f64::INFINITY, 1.0f64);
        for _ in 0..64 {
            let s: f64 = knn[i]
                .iter()
                .map(|e| {
                    let d = e.1 - rho;
                    if d > 0.0 { (-d / sigma).exp() } else { 1.0 }
                })
                .sum();
            if (s - target).abs() < 1e-5 {
                break;
            }
            if s > target {
                hi = sigma;
                sigma = (lo + hi) / 2.0;
            } else {
                lo = sigma;
                sigma = if hi.is_finite() {
                    (lo + hi) / 2.0
                } else {
                    sigma * 2.0
                };
            }
        }
        for &(j, d) in &knn[i] {
            let dd = d - rho;
            let v = if dd <= 0.0 { 1.0 } else { (-dd / sigma).exp() };
            w[i].insert(j, v);
        }
    }
    // Fuzzy union: a + b - a*b.
    let mut edges: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        let mut keys: Vec<usize> = w[i].keys().copied().collect();
        keys.sort_unstable();
        for j in keys {
            let a = w[i][&j];
            let b = w[j].get(&i).copied().unwrap_or(0.0);
            edges.push((i, j, a + b - a * b));
        }
    }
    for i in 0..n {
        let mut keys: Vec<usize> = w[i].keys().copied().collect();
        keys.sort_unstable();
        for j in keys {
            if !w[j].contains_key(&i) {
                edges.push((j, i, w[i][&j]));
            }
        }
    }
    let mut rng = Rng(SEED ^ 0xA5A5);
    let mut y = spectral_init(n, &edges, &mut rng);
    for r in y.iter_mut() {
        r[0] += 1e-4 * (rng.uniform() - 0.5);
        r[1] += 1e-4 * (rng.uniform() - 0.5);
    }
    let (a, b) = fit_ab(UMAP_MIN_DIST);
    let wmax = edges.iter().map(|e| e.2).fold(0.0, f64::max);
    let epochs = if n <= 10_000 { 500 } else { 200 };
    let neg = 5;
    let clip = |v: f64| v.clamp(-4.0, 4.0);
    // Undirected edge list (each pair once, i < j) for the SGD.
    let und: Vec<(usize, usize, f64)> = edges.iter().filter(|e| e.0 < e.1).copied().collect();
    for ep in 0..epochs {
        let alpha = 1.0 - ep as f64 / epochs as f64;
        for &(i, j, wt) in &und {
            if rng.uniform() > wt / wmax {
                continue;
            }
            let dx = y[i][0] - y[j][0];
            let dy = y[i][1] - y[j][1];
            let d2 = dx * dx + dy * dy;
            if d2 > 0.0 {
                let g = -2.0 * a * b * d2.powf(b - 1.0) / (1.0 + a * d2.powf(b));
                let (gx, gy) = (clip(g * dx), clip(g * dy));
                y[i][0] += alpha * gx;
                y[i][1] += alpha * gy;
                y[j][0] -= alpha * gx;
                y[j][1] -= alpha * gy;
            }
            for _ in 0..neg {
                let k = rng.below(n);
                if k == i {
                    continue;
                }
                let dx = y[i][0] - y[k][0];
                let dy = y[i][1] - y[k][1];
                let d2 = dx * dx + dy * dy;
                let g = if d2 > 0.0 {
                    2.0 * b / ((0.001 + d2) * (1.0 + a * d2.powf(b)))
                } else {
                    0.0
                };
                let (gx, gy) = if g > 0.0 {
                    (clip(g * dx), clip(g * dy))
                } else {
                    (4.0, 4.0)
                };
                y[i][0] += alpha * gx;
                y[i][1] += alpha * gy;
            }
        }
    }
    y
}

/// Write `analysis/{tsne,umap}/.../projection.csv`. `kind` is "tsne" or
/// "umap"; the header is `Barcode,TSNE-1,TSNE-2` or `Barcode,UMAP-1,UMAP-2`.
/// Numbers use the shortest round-trip decimal representation.
pub fn write_projection_csv(
    path: &Path,
    barcodes: &[String],
    coords: &[[f64; 2]],
    kind: &str,
) -> std::io::Result<()> {
    let label = match kind {
        "tsne" => "TSNE",
        "umap" => "UMAP",
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "kind must be \"tsne\" or \"umap\"",
            ));
        }
    };
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Barcode,{label}-1,{label}-2")?;
    for (b, c) in barcodes.iter().zip(coords) {
        writeln!(f, "{},{},{}", b, c[0], c[1])?;
    }
    f.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blobs() -> (Vec<Vec<f64>>, Vec<usize>) {
        let mut rng = Rng(7);
        let mut x = Vec::new();
        let mut l = Vec::new();
        for c in 0..3 {
            for _ in 0..40 {
                x.push(
                    (0..5)
                        .map(|d| if d == c { 10.0 } else { 0.0 } + rng.uniform())
                        .collect(),
                );
                l.push(c);
            }
        }
        (x, l)
    }

    fn separated(y: &[[f64; 2]], l: &[usize]) -> bool {
        // nearest neighbour shares the label for nearly all points
        let mut ok = 0;
        for i in 0..y.len() {
            let j = (0..y.len())
                .filter(|&j| j != i)
                .min_by(|&a, &b| sq_dist(&y[i], &y[a]).total_cmp(&sq_dist(&y[i], &y[b])))
                .unwrap();
            if l[i] == l[j] {
                ok += 1;
            }
        }
        ok as f64 / y.len() as f64 > 0.95
    }

    #[test]
    fn tsne_separates_and_is_deterministic() {
        let (x, l) = blobs();
        let a = tsne_2d(&x);
        assert_eq!(a, tsne_2d(&x));
        assert!(separated(&a, &l));
    }

    #[test]
    fn umap_separates_and_is_deterministic() {
        let (x, l) = blobs();
        let a = umap_2d(&x);
        assert_eq!(a, umap_2d(&x));
        assert!(a.iter().all(|r| r[0].is_finite() && r[1].is_finite()));
        assert!(separated(&a, &l));
    }

    #[test]
    fn ab_fit_reasonable() {
        let (a, b) = fit_ab(0.3);
        assert!(
            (a - 0.99).abs() < 0.05 && (b - 1.11).abs() < 0.05,
            "{a} {b}"
        );
    }

    #[test]
    fn csv_format() {
        let p = std::env::temp_dir().join("cr_embed_test/tsne/projection.csv");
        write_projection_csv(&p, &["AAAC-1".into()], &[[1.5, -2.0]], "tsne").unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        assert_eq!(s, "Barcode,TSNE-1,TSNE-2\nAAAC-1,1.5,-2\n");
        assert!(write_projection_csv(&p, &[], &[], "x").is_err());
    }

    fn read_csv(p: &str) -> (Vec<String>, Vec<Vec<f64>>) {
        let s = std::fs::read_to_string(p).unwrap();
        let mut b = Vec::new();
        let mut v = Vec::new();
        for line in s.lines().skip(1) {
            let mut it = line.split(',');
            b.push(it.next().unwrap().to_string());
            v.push(it.map(|t| t.parse().unwrap()).collect());
        }
        (b, v)
    }

    fn knn_sets(y: &[Vec<f64>], k: usize) -> Vec<Vec<usize>> {
        (0..y.len())
            .map(|i| {
                let mut d: Vec<(usize, f64)> = (0..y.len())
                    .filter(|&j| j != i)
                    .map(|j| (j, sq_dist(&y[i], &y[j])))
                    .collect();
                d.sort_by(|a, b| a.1.total_cmp(&b.1));
                d.iter().take(k).map(|e| e.0).collect()
            })
            .collect()
    }

    fn overlap(a: &[Vec<f64>], b: &[Vec<f64>], k: usize) -> f64 {
        let (ka, kb) = (knn_sets(a, k), knn_sets(b, k));
        let s: usize = ka
            .iter()
            .zip(&kb)
            .map(|(p, q)| p.iter().filter(|e| q.contains(e)).count())
            .sum();
        s as f64 / (a.len() * k) as f64
    }

    fn trust(hi: &[Vec<f64>], lo: &[Vec<f64>], k: usize) -> f64 {
        let n = hi.len();
        let klo = knn_sets(lo, k);
        let mut pen = 0.0;
        for i in 0..n {
            let mut d: Vec<(usize, f64)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (j, sq_dist(&hi[i], &hi[j])))
                .collect();
            d.sort_by(|a, b| a.1.total_cmp(&b.1));
            let mut rank = vec![0usize; n];
            for (r, e) in d.iter().enumerate() {
                rank[e.0] = r + 1;
            }
            for &j in &klo[i] {
                if rank[j] > k {
                    pen += (rank[j] - k) as f64;
                }
            }
        }
        let nf = n as f64;
        let kf = k as f64;
        1.0 - 2.0 / (nf * kf * (2.0 * nf - 3.0 * kf - 1.0)) * pen
    }

    fn silhouette(y: &[Vec<f64>], lab: &[usize]) -> f64 {
        let n = y.len();
        let nc = lab.iter().max().unwrap() + 1;
        let mut tot = 0.0;
        for i in 0..n {
            let mut sum = vec![0.0; nc];
            let mut cnt = vec![0usize; nc];
            for j in 0..n {
                if j != i {
                    sum[lab[j]] += sq_dist(&y[i], &y[j]).sqrt();
                    cnt[lab[j]] += 1;
                }
            }
            if cnt[lab[i]] == 0 {
                continue;
            }
            let a = sum[lab[i]] / cnt[lab[i]] as f64;
            let b = (0..nc)
                .filter(|&c| c != lab[i] && cnt[c] > 0)
                .map(|c| sum[c] / cnt[c] as f64)
                .fold(f64::INFINITY, f64::min);
            tot += (b - a) / a.max(b);
        }
        tot / n as f64
    }

    #[test]
    #[ignore = "reads the cellranger reference outputs"]
    fn reference_pbmc1k() {
        let r = "/Users/benjamin/rustar-parity/cr/ref/pbmc1k/outs/analysis";
        let (bc, pca) = read_csv(&format!(
            "{r}/pca/gene_expression_10_components/projection.csv"
        ));
        let (_, rt) = read_csv(&format!(
            "{r}/tsne/gene_expression_2_components/projection.csv"
        ));
        let (_, ru) = read_csv(&format!(
            "{r}/umap/gene_expression_2_components/projection.csv"
        ));
        let (cb, cl) = read_csv(&format!(
            "{r}/clustering/gene_expression_graphclust/clusters.csv"
        ));
        assert_eq!(bc, cb);
        let lab: Vec<usize> = cl.iter().map(|v| v[0] as usize - 1).collect();
        let t: Vec<Vec<f64>> = tsne_2d(&pca).iter().map(|r| r.to_vec()).collect();
        let u: Vec<Vec<f64>> = umap_2d(&pca).iter().map(|r| r.to_vec()).collect();
        for (name, mine, refr) in [("tsne", &t, &rt), ("umap", &u, &ru)] {
            eprintln!(
                "{name}: knn15 overlap {:.3}, trust(mine) {:.4}, trust(ref) {:.4}, sil(mine) {:.3}, sil(ref) {:.3}",
                overlap(mine, refr, 15),
                trust(&pca, mine, 15),
                trust(&pca, refr, 15),
                silhouette(mine, &lab),
                silhouette(refr, &lab)
            );
        }
        eprintln!(
            "ref tsne vs ref umap knn15 overlap {:.3}",
            overlap(&rt, &ru, 15)
        );
    }
}
