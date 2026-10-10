//! Clustering of the PCA projection.
//!
//! * **k-means** (`k = 2..=10`): k-means++ seeding, Lloyd iterations, the best
//!   of several restarts by within-cluster sum of squares. The reference run
//!   is a single (not globally optimal) local minimum, so for larger `k` the
//!   partitions agree closely but not always exactly.
//! * **graph-based**: k-nearest-neighbour graph on the projection with
//!   `k = round(-230 + 120 * log10(n_cells))` (from the public `reanalyze`
//!   parameter description), edges weighted by the Jaccard index of the two
//!   neighbour sets, then Louvain modularity optimisation (Blondel et al.
//!   2008), best modularity over several randomised restarts.
//!
//! Cluster numbers are assigned by decreasing cluster size, as in the
//! reference output.

use super::{create, rng::Rng, write_line};
use crate::error::Error;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

fn sq_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// Relabel to `1..=K` by decreasing size (ties: smaller original id first).
pub fn relabel_by_size(labels: &[usize]) -> Vec<u32> {
    let kmax = labels.iter().copied().max().map_or(0, |m| m + 1);
    let mut size = vec![0usize; kmax];
    for &l in labels {
        size[l] += 1;
    }
    let mut ids: Vec<usize> = (0..kmax).filter(|&i| size[i] > 0).collect();
    ids.sort_by(|&a, &b| size[b].cmp(&size[a]).then(a.cmp(&b)));
    let mut map = vec![0u32; kmax];
    for (new, &old) in ids.iter().enumerate() {
        map[old] = new as u32 + 1;
    }
    labels.iter().map(|&l| map[l]).collect()
}

fn kmeans_once(x: &[f64], d: usize, k: usize, seed: u64) -> (Vec<usize>, f64) {
    let n = x.len() / d;
    let mut rng = Rng::new(seed);
    let row = |i: usize| &x[i * d..(i + 1) * d];
    // k-means++ seeding.
    let mut centers: Vec<f64> = Vec::with_capacity(k * d);
    centers.extend_from_slice(row(rng.below(n)));
    let mut dmin: Vec<f64> = (0..n).map(|i| sq_dist(row(i), &centers[0..d])).collect();
    while centers.len() < k * d {
        let total: f64 = dmin.iter().sum();
        let mut target = rng.f64() * total;
        let mut pick = n - 1;
        for (i, &dm) in dmin.iter().enumerate() {
            if target < dm {
                pick = i;
                break;
            }
            target -= dm;
        }
        let c0 = centers.len();
        centers.extend_from_slice(row(pick));
        for i in 0..n {
            let dd = sq_dist(row(i), &centers[c0..c0 + d]);
            if dd < dmin[i] {
                dmin[i] = dd;
            }
        }
    }
    let mut labels = vec![usize::MAX; n];
    for _ in 0..300 {
        let new: Vec<usize> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut best = (f64::INFINITY, 0);
                for c in 0..k {
                    let dd = sq_dist(row(i), &centers[c * d..(c + 1) * d]);
                    if dd < best.0 {
                        best = (dd, c);
                    }
                }
                best.1
            })
            .collect();
        let changed = new != labels;
        labels = new;
        let mut sums = vec![0.0; k * d];
        let mut cnt = vec![0usize; k];
        for i in 0..n {
            cnt[labels[i]] += 1;
            for t in 0..d {
                sums[labels[i] * d + t] += x[i * d + t];
            }
        }
        for c in 0..k {
            if cnt[c] == 0 {
                // Re-seed an empty cluster at the point farthest from its centre.
                let far = (0..n)
                    .max_by(|&a, &b| {
                        let da = sq_dist(row(a), &centers[labels[a] * d..(labels[a] + 1) * d]);
                        let db = sq_dist(row(b), &centers[labels[b] * d..(labels[b] + 1) * d]);
                        da.partial_cmp(&db).unwrap()
                    })
                    .unwrap();
                centers[c * d..(c + 1) * d].copy_from_slice(row(far));
            } else {
                for t in 0..d {
                    centers[c * d + t] = sums[c * d + t] / cnt[c] as f64;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let inertia = (0..n)
        .map(|i| sq_dist(row(i), &centers[labels[i] * d..(labels[i] + 1) * d]))
        .sum();
    (labels, inertia)
}

/// K-means on a row-major `n x d` matrix; labels `1..=k` by decreasing size.
pub fn kmeans(x: &[f64], d: usize, k: usize) -> Vec<u32> {
    let best = (0..40u64)
        .into_par_iter()
        .map(|r| kmeans_once(x, d, k, 0xC0FFEE + 7919 * r + k as u64))
        .reduce_with(|a, b| if b.1 < a.1 { b } else { a })
        .unwrap();
    relabel_by_size(&best.0)
}

/// Number of neighbours for `n` cells.
pub fn graph_k(n: usize) -> usize {
    let k = (-230.0 + 120.0 * (n as f64).log10()).round();
    (k.max(10.0) as usize).min(n.saturating_sub(1)).max(1)
}

type Adj = Vec<Vec<(u32, f64)>>;

fn knn_jaccard_graph(x: &[f64], d: usize, k: usize) -> Adj {
    let n = x.len() / d;
    let nbrs: Vec<Vec<u32>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let xi = &x[i * d..(i + 1) * d];
            let mut dist: Vec<(f64, u32)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (sq_dist(xi, &x[j * d..(j + 1) * d]), j as u32))
                .collect();
            let kk = k.min(dist.len());
            dist.select_nth_unstable_by(kk - 1, |a, b| {
                a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1))
            });
            let mut v: Vec<u32> = dist[..kk].iter().map(|p| p.1).collect();
            v.sort_unstable();
            v
        })
        .collect();
    // Undirected union of the neighbour relations.
    let mut edges: Vec<(u32, u32)> = Vec::new();
    for (i, nb) in nbrs.iter().enumerate() {
        for &j in nb {
            edges.push(((i as u32).min(j), (i as u32).max(j)));
        }
    }
    edges.sort_unstable();
    edges.dedup();
    let w: Vec<f64> = edges
        .par_iter()
        .map(|&(a, b)| {
            let (na, nb) = (&nbrs[a as usize], &nbrs[b as usize]);
            let (mut i, mut j, mut inter) = (0, 0, 0usize);
            while i < na.len() && j < nb.len() {
                match na[i].cmp(&nb[j]) {
                    std::cmp::Ordering::Less => i += 1,
                    std::cmp::Ordering::Greater => j += 1,
                    std::cmp::Ordering::Equal => {
                        inter += 1;
                        i += 1;
                        j += 1;
                    }
                }
            }
            inter as f64 / (na.len() + nb.len() - inter) as f64
        })
        .collect();
    let mut adj: Adj = vec![Vec::new(); n];
    for (&(a, b), &wt) in edges.iter().zip(&w) {
        if wt > 0.0 {
            adj[a as usize].push((b, wt));
            adj[b as usize].push((a, wt));
        }
    }
    adj
}

/// One Louvain run; returns the community of every node and the modularity.
fn louvain(adj: &Adj, rng: &mut Rng) -> (Vec<usize>, f64) {
    let n0 = adj.len();
    let m2: f64 = adj.iter().flat_map(|a| a.iter().map(|e| e.1)).sum();
    let mut member: Vec<usize> = (0..n0).collect();
    let mut g = adj.clone();
    let mut selfw = vec![0.0; n0];
    if m2 == 0.0 {
        return (member, 0.0);
    }
    loop {
        let n = g.len();
        let strength: Vec<f64> = (0..n)
            .map(|i| g[i].iter().map(|e| e.1).sum::<f64>() + selfw[i])
            .collect();
        let mut comm: Vec<usize> = (0..n).collect();
        let mut tot = strength.clone();
        let mut order: Vec<usize> = (0..n).collect();
        let mut improved = false;
        let mut w_to = vec![0.0; n];
        for _pass in 0..100 {
            rng.shuffle(&mut order);
            let mut moved = false;
            for &i in &order {
                let ci = comm[i];
                let mut touched: Vec<usize> = Vec::new();
                for &(j, w) in &g[i] {
                    let cj = comm[j as usize];
                    if w_to[cj] == 0.0 {
                        touched.push(cj);
                    }
                    w_to[cj] += w;
                }
                tot[ci] -= strength[i];
                let gain = |c: usize| w_to[c] - tot[c] * strength[i] / m2;
                let mut best_c = ci;
                let mut best_g = gain(ci);
                for &c in &touched {
                    let gg = gain(c);
                    if gg > best_g + 1e-12 {
                        best_g = gg;
                        best_c = c;
                    }
                }
                tot[best_c] += strength[i];
                if best_c != ci {
                    comm[i] = best_c;
                    moved = true;
                    improved = true;
                }
                for &c in &touched {
                    w_to[c] = 0.0;
                }
                w_to[ci] = 0.0;
            }
            if !moved {
                break;
            }
        }
        if !improved {
            break;
        }
        // Renumber communities and aggregate.
        let mut remap = vec![usize::MAX; n];
        let mut nc = 0;
        for &c in &comm {
            if remap[c] == usize::MAX {
                remap[c] = nc;
                nc += 1;
            }
        }
        for m in member.iter_mut() {
            *m = remap[comm[*m]];
        }
        let mut agg: Vec<std::collections::BTreeMap<u32, f64>> = vec![Default::default(); nc];
        let mut nself = vec![0.0; nc];
        for i in 0..n {
            let a = remap[comm[i]];
            nself[a] += selfw[i];
            for &(j, w) in &g[i] {
                let b = remap[comm[j as usize]];
                if a == b {
                    nself[a] += w;
                } else {
                    *agg[a].entry(b as u32).or_insert(0.0) += w;
                }
            }
        }
        g = agg
            .into_iter()
            .map(|m| m.into_iter().collect::<Vec<_>>())
            .collect();
        selfw = nself;
        if nc == 1 {
            break;
        }
    }
    // Modularity on the original graph.
    let nc = member.iter().copied().max().map_or(0, |m| m + 1);
    let mut inside = vec![0.0; nc];
    let mut tot = vec![0.0; nc];
    for (i, edges) in adj.iter().enumerate() {
        for &(j, w) in edges {
            tot[member[i]] += w;
            if member[i] == member[j as usize] {
                inside[member[i]] += w;
            }
        }
    }
    let q = (0..nc)
        .map(|c| inside[c] / m2 - (tot[c] / m2) * (tot[c] / m2))
        .sum();
    (member, q)
}

/// Graph-based clustering of an `n x d` projection.
pub fn graph_cluster(x: &[f64], d: usize) -> Vec<u32> {
    let n = x.len() / d;
    if n < 3 {
        return vec![1; n];
    }
    let adj = knn_jaccard_graph(x, d, graph_k(n));
    let best = (0..40u64)
        .into_par_iter()
        .map(|r| louvain(&adj, &mut Rng::new(0xFACADE + 104_729 * r)))
        .reduce_with(|a, b| if b.1 > a.1 { b } else { a })
        .unwrap();
    relabel_by_size(&best.0)
}

pub fn write_clusters(barcodes: &[String], labels: &[u32], path: &Path) -> Result<(), Error> {
    let mut w = create(path)?;
    let io = |e| Error::io(e, path);
    write_line(&mut w, "Barcode,Cluster").map_err(io)?;
    for (b, l) in barcodes.iter().zip(labels) {
        write_line(&mut w, &format!("{b},{l}")).map_err(io)?;
    }
    w.flush().map_err(io)
}
