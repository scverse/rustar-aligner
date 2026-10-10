//! Cell calling for 3' v3 libraries: an order-of-magnitude call from barcode
//! UMI totals, followed by an EmptyDrops-style rescue of barcodes whose
//! expression profile differs from the ambient profile.
//!
//! This is an independent implementation written from a behavioural
//! specification and the published methods: EmptyDrops (Lun et al. 2019),
//! Simple Good-Turing (Gale and Sampson 1995), NumPy's legacy `RandomState`
//! (MT19937) integer sampling, `default_rng` (PCG64 seeded via
//! `SeedSequence`), NumPy's default introsort `argsort`, and the Cephes
//! `lgam` used by SciPy's `gammaln`.

#![allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::manual_midpoint,
    clippy::unreadable_literal
)]

use crate::solo::libcxx_rng::Mt19937;
use rayon::prelude::*;

/// Column-compressed feature-by-barcode count matrix.
#[derive(Debug, Clone, Default)]
pub struct CscCounts {
    pub n_features: usize,
    pub col_ptr: Vec<usize>,
    pub rows: Vec<u32>,
    pub vals: Vec<u32>,
}

impl CscCounts {
    pub fn n_cols(&self) -> usize {
        self.col_ptr.len().saturating_sub(1)
    }

    /// Feature indices and counts of the nonzero entries of column `c`.
    pub fn col(&self, c: usize) -> (&[u32], &[u32]) {
        let (a, b) = (self.col_ptr[c], self.col_ptr[c + 1]);
        (&self.rows[a..b], &self.vals[a..b])
    }

    /// Total UMI count of every column.
    pub fn umis(&self) -> Vec<u64> {
        (0..self.n_cols())
            .map(|c| self.col(c).1.iter().map(|&v| u64::from(v)).sum())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const EMPTY_RANK_LO: usize = 45_000;
const EMPTY_RANK_HI: usize = 90_000;
const MIN_UMIS: u64 = 500;
const N_SIMS: usize = 100_000;
const FDR: f64 = 0.001;
const SIM_CHUNK: usize = 256;

// ---------------------------------------------------------------------------
// Random number generators
// ---------------------------------------------------------------------------

/// NumPy legacy `RandomState` restricted to `randint(0, n)`.
struct LegacyRng {
    mt: Mt19937,
}

impl LegacyRng {
    fn new(seed: u32) -> Self {
        Self {
            mt: Mt19937::new(seed),
        }
    }

    /// Uniform integer in `[0, n)` by masked rejection sampling.
    fn randint(&mut self, n: u64) -> u64 {
        let r = n - 1;
        if r == 0 {
            return 0;
        }
        if let Ok(r32) = u32::try_from(r) {
            let mut mask = r32;
            mask |= mask >> 1;
            mask |= mask >> 2;
            mask |= mask >> 4;
            mask |= mask >> 8;
            mask |= mask >> 16;
            loop {
                let v = self.mt.next_u32() & mask;
                if v <= r32 {
                    return u64::from(v);
                }
            }
        }
        let mut mask = r;
        for s in [1, 2, 4, 8, 16, 32] {
            mask |= mask >> s;
        }
        loop {
            let hi = u64::from(self.mt.next_u32());
            let lo = u64::from(self.mt.next_u32());
            let v = ((hi << 32) | lo) & mask;
            if v <= r {
                return v;
            }
        }
    }
}

/// NumPy `SeedSequence(entropy)` for a single 32-bit entropy word, producing
/// `n` 32-bit state words.
fn seed_sequence_words(entropy: u32, n: usize) -> Vec<u32> {
    const INIT_A: u32 = 0x43b0_d7e5;
    const MULT_A: u32 = 0x931e_8875;
    const INIT_B: u32 = 0x8b51_f9dd;
    const MULT_B: u32 = 0x58f3_8ded;
    const MIX_L: u32 = 0xca01_f9dd;
    const MIX_R: u32 = 0x4973_f715;
    const POOL: usize = 4;

    fn hashmix(value: u32, hc: &mut u32) -> u32 {
        let mut v = value ^ *hc;
        *hc = hc.wrapping_mul(MULT_A);
        v = v.wrapping_mul(*hc);
        v ^ (v >> 16)
    }
    fn mix(x: u32, y: u32) -> u32 {
        let r = MIX_L.wrapping_mul(x).wrapping_sub(MIX_R.wrapping_mul(y));
        r ^ (r >> 16)
    }

    let mut hc = INIT_A;
    let mut pool = [0u32; POOL];
    for (i, p) in pool.iter_mut().enumerate() {
        let e = if i == 0 { entropy } else { 0 };
        *p = hashmix(e, &mut hc);
    }
    for i_src in 0..POOL {
        for i_dst in 0..POOL {
            if i_src != i_dst {
                let h = hashmix(pool[i_src], &mut hc);
                pool[i_dst] = mix(pool[i_dst], h);
            }
        }
    }
    let mut hc = INIT_B;
    (0..n)
        .map(|i| {
            let mut d = pool[i % POOL] ^ hc;
            hc = hc.wrapping_mul(MULT_B);
            d = d.wrapping_mul(hc);
            d ^ (d >> 16)
        })
        .collect()
}

const PCG_MULT: u128 = 0x2360_ED05_1FC6_5DA4_4385_DF64_9FCC_F645;

/// NumPy's PCG64 (XSL-RR 128/64) generator.
#[derive(Debug, Clone)]
struct Pcg64 {
    state: u128,
    inc: u128,
}

impl Pcg64 {
    /// `default_rng(seed)` for a 32-bit integer seed.
    fn from_seed(seed: u32) -> Self {
        let w = seed_sequence_words(seed, 8);
        let u = |i: usize| u64::from(w[2 * i]) | (u64::from(w[2 * i + 1]) << 32);
        let initstate = (u128::from(u(0)) << 64) | u128::from(u(1));
        let initseq = (u128::from(u(2)) << 64) | u128::from(u(3));
        let mut g = Self {
            state: 0,
            inc: (initseq << 1) | 1,
        };
        g.step();
        g.state = g.state.wrapping_add(initstate);
        g.step();
        g
    }

    fn step(&mut self) {
        self.state = self.state.wrapping_mul(PCG_MULT).wrapping_add(self.inc);
    }

    fn next_u64(&mut self) -> u64 {
        self.step();
        let x = ((self.state >> 64) as u64) ^ (self.state as u64);
        x.rotate_right((self.state >> 122) as u32)
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Skip `delta` draws in O(log delta).
    fn advance(&mut self, mut delta: u128) {
        let (mut acc_mult, mut acc_plus) = (1u128, 0u128);
        let (mut cur_mult, mut cur_plus) = (PCG_MULT, self.inc);
        while delta > 0 {
            if delta & 1 == 1 {
                acc_mult = acc_mult.wrapping_mul(cur_mult);
                acc_plus = acc_plus.wrapping_mul(cur_mult).wrapping_add(cur_plus);
            }
            cur_plus = cur_mult.wrapping_add(1).wrapping_mul(cur_plus);
            cur_mult = cur_mult.wrapping_mul(cur_mult);
            delta >>= 1;
        }
        self.state = acc_mult.wrapping_mul(self.state).wrapping_add(acc_plus);
    }
}

// ---------------------------------------------------------------------------
// Rounding helpers
// ---------------------------------------------------------------------------

fn round_half_even(x: f64) -> f64 {
    let r = x.round();
    if (x - x.trunc()).abs() == 0.5 {
        2.0 * (x / 2.0).round()
    } else {
        r
    }
}

// ---------------------------------------------------------------------------
// NumPy default argsort (introsort on the index array)
// ---------------------------------------------------------------------------

fn numpy_argsort(v: &[u64]) -> Vec<usize> {
    let n = v.len();
    let mut idx: Vec<usize> = (0..n).collect();
    if n < 2 {
        return idx;
    }
    let mut pl: isize = 0;
    let mut pr: isize = n as isize - 1;
    let mut cdepth: i32 = 2 * (usize::BITS - 1 - n.leading_zeros()) as i32;
    let mut stack: Vec<(isize, isize, i32)> = Vec::new();
    let val = |idx: &Vec<usize>, p: isize| v[idx[p as usize]];
    loop {
        while pr - pl > 16 {
            let pm = pl + ((pr - pl) >> 1);
            if val(&idx, pm) < val(&idx, pl) {
                idx.swap(pm as usize, pl as usize);
            }
            if val(&idx, pr) < val(&idx, pm) {
                idx.swap(pr as usize, pm as usize);
            }
            if val(&idx, pm) < val(&idx, pl) {
                idx.swap(pm as usize, pl as usize);
            }
            let vp = val(&idx, pm);
            let mut pi = pl;
            let mut pj = pr - 1;
            idx.swap(pm as usize, pj as usize);
            loop {
                loop {
                    pi += 1;
                    if val(&idx, pi) >= vp {
                        break;
                    }
                }
                loop {
                    pj -= 1;
                    if vp >= val(&idx, pj) {
                        break;
                    }
                }
                if pi >= pj {
                    break;
                }
                idx.swap(pi as usize, pj as usize);
            }
            idx.swap(pi as usize, (pr - 1) as usize);
            cdepth -= 1;
            if pi - pl < pr - pi {
                stack.push((pi + 1, pr, cdepth));
                pr = pi - 1;
            } else {
                stack.push((pl, pi - 1, cdepth));
                pl = pi + 1;
            }
        }
        // Insertion sort of the small partition.
        for pi in (pl + 1)..=pr {
            let vi = idx[pi as usize];
            let vp = v[vi];
            let mut pj = pi;
            while pj > pl && vp < val(&idx, pj - 1) {
                idx[pj as usize] = idx[(pj - 1) as usize];
                pj -= 1;
            }
            idx[pj as usize] = vi;
        }
        match stack.pop() {
            None => break,
            Some((l, r, d)) => {
                pl = l;
                pr = r;
                cdepth = d;
            }
        }
    }
    let _ = cdepth;
    idx
}

// ---------------------------------------------------------------------------
// Stage A: order-of-magnitude call
// ---------------------------------------------------------------------------

fn bootstrap(nz: &[u64], rng: &mut LegacyRng) -> Vec<u64> {
    let n = nz.len() as u64;
    let mut s: Vec<u64> = (0..nz.len()).map(|_| nz[rng.randint(n) as usize]).collect();
    s.sort_unstable();
    s
}

/// Number of elements of ascending `s` at least a tenth of the baseline value.
fn within_ordmag(s: &[u64], b: usize) -> usize {
    let v = s[s.len() - 1 - b] as f64;
    let cutoff = round_half_even(0.1 * v).max(1.0) as u64;
    s.len() - s.partition_point(|&x| x < cutoff)
}

fn recovered_candidates() -> Vec<u64> {
    let start = 1.0f64;
    let stop = 45_000f64.log2();
    let n = 2000usize;
    let step = (stop - start) / (n - 1) as f64;
    let mut out: Vec<u64> = Vec::new();
    for i in 0..n {
        let e = if i == n - 1 {
            stop
        } else {
            start + i as f64 * step
        };
        let rc = round_half_even(e.exp2()) as u64;
        if out.last() != Some(&rc) {
            out.push(rc);
        }
    }
    out
}

fn baseline_index(x: f64) -> usize {
    round_half_even(x * (1.0 - 0.99)) as usize
}

/// Order-of-magnitude call; returns the called column indices, ascending.
pub fn ordmag_cells(umis: &[u64]) -> Vec<usize> {
    let nz: Vec<u64> = umis.iter().copied().filter(|&u| u > 0).collect();
    if nz.is_empty() {
        return Vec::new();
    }
    let mut rng = LegacyRng::new(0);
    let cands = recovered_candidates();

    let mut sum = 0.0f64;
    for _ in 0..100 {
        let s = bootstrap(&nz, &mut rng);
        let mut best = (f64::INFINITY, cands[0]);
        for &rc in &cands {
            let b = baseline_index(rc as f64).min(s.len() - 1);
            let filtered = within_ordmag(&s, b) as f64;
            let rcf = rc as f64;
            let loss = (filtered - rcf) * (filtered - rcf) / rcf;
            if loss < best.0 {
                best = (loss, rc);
            }
        }
        sum += best.1 as f64;
    }
    let recovered = round_half_even(sum / 100.0).max(50.0);
    let b = baseline_index(recovered).min(nz.len() - 1);

    let mut csum = 0.0f64;
    for _ in 0..100 {
        let s = bootstrap(&nz, &mut rng);
        csum += within_ordmag(&s, b) as f64;
    }
    let mut nbcs = round_half_even(csum / 100.0) as usize;

    if nbcs > 0 {
        let mut sorted = nz.clone();
        sorted.sort_unstable_by(|a, c| c.cmp(a));
        let nbcs0 = nbcs.min(sorted.len());
        nbcs = nbcs0;
        let cutoff = sorted[nbcs - 1];
        let mut i = nbcs - 1;
        while sorted[i] == cutoff && i + 1 < sorted.len() {
            i += 1;
            if (i + 1 - nbcs) as f64 > 0.20 * nbcs as f64 {
                break;
            }
            nbcs = i + 1;
        }
    }

    // Stable ascending sort, reversed: later columns win ties.
    let mut order: Vec<usize> = (0..umis.len()).collect();
    order.sort_by_key(|&c| umis[c]);
    order.reverse();
    order.truncate(nbcs);
    order.sort_unstable();
    order
}

// ---------------------------------------------------------------------------
// Simple Good-Turing
// ---------------------------------------------------------------------------

/// Simple Good-Turing smoothing of per-feature frequencies. Returns the
/// per-feature proportions and the mass of the unseen class, or `None` when
/// fewer than 10 distinct frequency values are present.
fn simple_good_turing(freqs: &[u64]) -> Option<(Vec<f64>, f64)> {
    let mut sorted: Vec<u64> = freqs.to_vec();
    sorted.sort_unstable();
    let mut r: Vec<f64> = Vec::new();
    let mut nr: Vec<f64> = Vec::new();
    for &f in &sorted {
        if r.last() == Some(&(f as f64)) {
            *nr.last_mut().unwrap() += 1.0;
        } else {
            r.push(f as f64);
            nr.push(1.0);
        }
    }
    let n = r.len();
    if n < 10 {
        return None;
    }
    let total: f64 = r.iter().zip(&nr).map(|(a, b)| a * b).sum();

    // Averaging transform.
    let d: Vec<f64> = (0..n)
        .map(|k| if k == 0 { 1.0 } else { r[k] - r[k - 1] })
        .collect();
    let z: Vec<f64> = (0..n)
        .map(|k| {
            let dr = if k + 1 < n {
                0.5 * (d[k + 1] + d[k])
            } else {
                d[k]
            };
            nr[k] / dr
        })
        .collect();
    let lx: Vec<f64> = r.iter().map(|x| x.ln()).collect();
    let ly: Vec<f64> = z.iter().map(|x| x.ln()).collect();
    let nf = n as f64;
    let mx = lx.iter().sum::<f64>() / nf;
    let my = ly.iter().sum::<f64>() / nf;
    let sxy: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = lx.iter().map(|a| (a - mx) * (a - mx)).sum();
    let mut slope = sxy / sxx;
    if slope >= -1.0 {
        slope = -1.0;
    }

    let lin: Vec<f64> = r
        .iter()
        .map(|&x| (1.0 + 1.0 / x).powf(1.0 + slope))
        .collect();
    let mut turing = vec![0.0f64; n];
    let mut sd = vec![1.0f64; n];
    for j in 0..n.saturating_sub(1) {
        if r[j + 1] == r[j] + 1.0 {
            turing[j] = (r[j] + 1.0) / r[j] * nr[j + 1] / nr[j];
            sd[j] = (j as f64 + 2.0) / nr[j] * (nr[j + 1] * (1.0 + nr[j + 1] / nr[j])).sqrt();
        }
    }
    let mut comb = vec![0.0f64; n];
    let mut mode = true;
    for j in 0..n {
        if !mode {
            comb[j] = lin[j];
        } else if (lin[j] - turing[j]).abs() * (1.0 + j as f64) / sd[j] > 1.65 {
            comb[j] = turing[j];
        } else {
            mode = false;
            comb[j] = lin[j];
        }
    }
    let sumpraw: f64 = (0..n).map(|j| comb[j] * r[j] * nr[j]).sum::<f64>() / total;
    let p0 = nr[0] / total;
    let rstar: Vec<f64> = (0..n)
        .map(|j| r[j] * comb[j] * (1.0 - p0) / sumpraw)
        .collect();
    let norm: f64 = (0..n).map(|j| nr[j] * rstar[j]).sum();
    let prob: Vec<f64> = freqs
        .iter()
        .map(|&f| {
            let j = r.partition_point(|&x| x < f as f64);
            (1.0 - p0) * rstar[j] / norm
        })
        .collect();
    Some((prob, p0))
}

// ---------------------------------------------------------------------------
// ln Gamma at integers (Cephes lgam, as used by SciPy gammaln)
// ---------------------------------------------------------------------------

fn lgam(x: f64) -> f64 {
    if x < 13.0 {
        let mut z = 1.0f64;
        let mut p = 0.0f64;
        let mut u = x;
        while u >= 3.0 {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while u < 2.0 {
            z /= u;
            p += 1.0;
            u = x + p;
        }
        // For integer x, u is now exactly 2.
        return z.ln();
    }
    let mut q = (x - 0.5) * x.ln() - x + 0.918_938_533_204_672_8;
    if x > 1.0e8 {
        return q;
    }
    let p = 1.0 / (x * x);
    if x >= 1000.0 {
        q += ((7.936_507_936_507_937e-4 * p - 2.777_777_777_777_778e-3) * p
            + 0.083_333_333_333_333_33)
            / x;
    } else {
        let a = [
            8.116_141_674_705_085e-4,
            -5.950_619_042_843_014e-4,
            7.936_503_404_577_169e-4,
            -2.777_777_777_300_997e-3,
            8.333_333_333_333_319e-2,
        ];
        let mut acc = 0.0;
        for c in a {
            acc = acc * p + c;
        }
        q += acc / x;
    }
    q
}

fn lgamma_table(max_k: usize) -> Vec<f64> {
    (0..=max_k.max(2)).map(|k| lgam(k as f64 + 1.0)).collect()
}

// ---------------------------------------------------------------------------
// Stage B: non-ambient barcodes
// ---------------------------------------------------------------------------

/// Ambient profile over `n_eval` evaluated features, or `None` if SGT fails.
fn ambient_profile(
    m: &CscCounts,
    pos: &[usize],
    n_eval: usize,
    ambient: &[usize],
) -> Option<Vec<f64>> {
    let mut sums = vec![0u64; n_eval];
    for &c in ambient {
        let (rows, vals) = m.col(c);
        for (&r, &v) in rows.iter().zip(vals) {
            sums[pos[r as usize]] += u64::from(v);
        }
    }
    let nz_idx: Vec<usize> = (0..n_eval).filter(|&i| sums[i] > 0).collect();
    let freqs: Vec<u64> = nz_idx.iter().map(|&i| sums[i]).collect();
    let (p, p0) = simple_good_turing(&freqs)?;
    let n0 = n_eval - nz_idx.len();
    let mut profile = vec![0.0f64; n_eval];
    if n0 == 0 {
        let s: f64 = p.iter().sum();
        for (&i, &pi) in nz_idx.iter().zip(&p) {
            profile[i] = pi / s;
        }
    } else {
        let fill = p0 / n0 as f64;
        profile.fill(fill);
        for (&i, &pi) in nz_idx.iter().zip(&p) {
            profile[i] = pi;
        }
    }
    Some(profile)
}

/// Barcodes (not in `orig_cells`) whose profile differs from the ambient one.
pub fn nonambient_cells(m: &CscCounts, feature_ids: &[String], orig_cells: &[usize]) -> Vec<usize> {
    if orig_cells.is_empty() {
        return Vec::new();
    }
    let totals = m.umis();
    let nb = totals.len();

    // 3.1 background range.
    let mut order = numpy_argsort(&totals);
    order.reverse();
    let hi = EMPTY_RANK_HI.min(nb);
    let lo = EMPTY_RANK_LO.min(hi);
    let bg = &order[lo..hi];
    let mut ambient: Vec<usize> = bg.iter().copied().filter(|&c| totals[c] > 0).collect();
    ambient.sort_unstable();
    if ambient.is_empty() {
        return Vec::new();
    }
    let max_bg = bg.iter().map(|&c| totals[c]).max().unwrap_or(0);

    // 3.2 evaluated features ordered by id.
    let mut present = vec![false; m.n_features];
    for &r in &m.rows {
        present[r as usize] = true;
    }
    let mut eval: Vec<usize> = (0..m.n_features).filter(|&f| present[f]).collect();
    eval.sort_by(|&a, &b| feature_ids[a].as_bytes().cmp(feature_ids[b].as_bytes()));
    let n_eval = eval.len();
    let mut pos = vec![usize::MAX; m.n_features];
    for (p, &f) in eval.iter().enumerate() {
        pos[f] = p;
    }

    // 3.3 ambient profile.
    let Some(profile) = ambient_profile(m, &pos, n_eval, &ambient) else {
        return Vec::new();
    };

    // 3.5 candidates.
    let min_umis = MIN_UMIS.max(1 + max_bg);
    let mut is_orig = vec![false; nb];
    for &c in orig_cells {
        is_orig[c] = true;
    }
    let cands: Vec<usize> = (0..nb)
        .filter(|&c| !is_orig[c] && totals[c] >= min_umis)
        .collect();
    if cands.is_empty() {
        return Vec::new();
    }

    // 3.6 observed log-likelihoods.
    let max_n = cands.iter().map(|&c| totals[c]).max().unwrap() as usize;
    let lg = lgamma_table(max_n);
    let ln_prof: Vec<f64> = profile.iter().map(|p| p.ln()).collect();
    let obs: Vec<f64> = cands
        .iter()
        .map(|&c| {
            let (rows, vals) = m.col(c);
            let mut l = lg[totals[c] as usize];
            for &v in vals {
                l -= lg[v as usize];
            }
            for (&r, &v) in rows.iter().zip(vals) {
                l += f64::from(v) * ln_prof[pos[r as usize]];
            }
            l
        })
        .collect();

    // 3.7 simulation.
    let mut cum = Vec::with_capacity(n_eval);
    let mut acc = 0.0f64;
    for &p in &profile {
        acc += p;
        cum.push(acc);
    }
    let ln_j: Vec<f64> = (0..=max_n).map(|j| (j as f64).ln()).collect();
    // Candidate indices grouped by total, ascending total.
    let mut by_total: Vec<usize> = (0..cands.len()).collect();
    by_total.sort_by_key(|&i| totals[cands[i]]);

    let base = Pcg64::from_seed(42);
    let n_chunks = N_SIMS.div_ceil(SIM_CHUNK);
    let counts: Vec<u32> = (0..n_chunks)
        .into_par_iter()
        .map(|ch| {
            let first = ch * SIM_CHUNK;
            let last = (first + SIM_CHUNK).min(N_SIMS);
            let mut g = base.clone();
            g.advance(first as u128 * max_n as u128);
            let mut cnt = vec![0u32; n_eval];
            let mut drawn: Vec<usize> = Vec::with_capacity(max_n);
            let mut res = vec![0u32; cands.len()];
            for _ in first..last {
                drawn.clear();
                let mut ll = 0.0f64;
                let mut ptr = 0usize;
                for j in 1..=max_n {
                    let u = g.next_f64();
                    let f = cum.partition_point(|&x| x < u).min(n_eval - 1);
                    cnt[f] += 1;
                    drawn.push(f);
                    ll += ln_j[j] - ln_j[cnt[f] as usize] + ln_prof[f];
                    while ptr < by_total.len() && totals[cands[by_total[ptr]]] as usize == j {
                        let ci = by_total[ptr];
                        if ll < obs[ci] {
                            res[ci] += 1;
                        }
                        ptr += 1;
                    }
                }
                for &f in &drawn {
                    cnt[f] = 0;
                }
            }
            res
        })
        .reduce(
            || vec![0u32; cands.len()],
            |mut a, b| {
                a.iter_mut().zip(&b).for_each(|(x, y)| *x += y);
                a
            },
        )
        .into_iter()
        .collect();

    // 3.8 Benjamini-Hochberg.
    let pvals: Vec<f64> = counts
        .iter()
        .map(|&c| (1.0 + f64::from(c)) / (1.0 + N_SIMS as f64))
        .collect();
    let m_c = pvals.len();
    let mut ord: Vec<usize> = (0..m_c).collect();
    ord.sort_by(|&a, &b| pvals[a].total_cmp(&pvals[b]));
    let mut adj = vec![0.0f64; m_c];
    let mut run = f64::INFINITY;
    for rank in (1..=m_c).rev() {
        let i = ord[rank - 1];
        run = run.min(m_c as f64 / rank as f64 * pvals[i]);
        adj[i] = run.min(1.0);
    }
    (0..m_c)
        .filter(|&i| adj[i] <= FDR)
        .map(|i| cands[i])
        .collect()
}

/// Final cell set: order-of-magnitude call united with the non-ambient
/// rescue, ascending and unique.
pub fn call_cells(m: &CscCounts, feature_ids: &[String]) -> Vec<usize> {
    let orig = ordmag_cells(&m.umis());
    let extra = nonambient_cells(m, feature_ids, &orig);
    let mut all = orig;
    all.extend(extra);
    all.sort_unstable();
    all.dedup();
    all
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufRead, BufReader};

    #[test]
    fn pcg64_first_floats() {
        let mut g = Pcg64::from_seed(42);
        assert_eq!(g.next_f64(), 0.7739560485559633);
        assert_eq!(g.next_f64(), 0.4388784397520523);
        assert_eq!(g.next_f64(), 0.8585979199113825);
    }

    #[test]
    fn pcg64_advance_equals_drawing() {
        let mut a = Pcg64::from_seed(42);
        for _ in 0..1000 {
            a.next_f64();
        }
        let mut b = Pcg64::from_seed(42);
        b.advance(1000);
        assert_eq!(a.next_f64(), 0.06206310654015523);
        assert_eq!(b.next_f64(), 0.06206310654015523);
    }

    #[test]
    fn legacy_randint() {
        let mut r = LegacyRng::new(0);
        let v: Vec<u64> = (0..3).map(|_| r.randint(10)).collect();
        assert_eq!(v, [5, 0, 3]);
        let mut r = LegacyRng::new(0);
        let v: Vec<u64> = (0..5).map(|_| r.randint(1_000_003)).collect();
        assert_eq!(v, [985772, 305711, 435829, 117952, 963395]);
    }

    #[test]
    fn lgamma_matches_scipy() {
        let t = lgamma_table(5000);
        assert_eq!(t[0], 0.0);
        assert_eq!(t[1], 0.0);
        assert!((t[5] - 120f64.ln()).abs() < 1e-15);
        for (k, want) in [
            (12, 19.98721449566189),
            (13, 22.55216385312342),
            (20, 42.335616460753485),
            (999, 5905.220423209181),
            (1000, 5912.128178488164),
            (5000, 37591.14350887677),
        ] {
            assert!((t[k] - want).abs() < 1e-9 * want.max(1.0), "k={k}");
        }
        assert!((t[20] - 2432902008176640000f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn argsort_sorts_and_keeps_small_inputs_stable() {
        // Inputs of 16 or fewer elements use a stable insertion sort.
        let small = [3u64, 1, 3, 0, 1, 2, 3, 0];
        assert_eq!(numpy_argsort(&small), [3, 7, 1, 4, 5, 0, 2, 6]);
        // Larger inputs go through the quicksort partitioning: the result is
        // a permutation ordered by value (tie order is implementation defined).
        let v: Vec<u64> = (0..5000u64).map(|i| (i * 2_654_435_761) % 97).collect();
        let idx = numpy_argsort(&v);
        let mut seen = vec![false; v.len()];
        for &i in &idx {
            assert!(!seen[i]);
            seen[i] = true;
        }
        assert!(idx.windows(2).all(|w| v[w[0]] <= v[w[1]]));
    }

    #[test]
    fn ordmag_empty() {
        assert_eq!(ordmag_cells(&[]).len(), 0);
        assert_eq!(ordmag_cells(&[0, 0]).len(), 0);
    }

    fn read_gz_lines(p: &std::path::Path) -> Vec<String> {
        let f = std::fs::File::open(p).unwrap();
        BufReader::new(flate2::read::MultiGzDecoder::new(f))
            .lines()
            .map(|l| l.unwrap())
            .collect()
    }

    #[test]
    #[ignore = "needs CR_RAW_DIR and CR_CELLS"]
    fn cell_calling_matches_cellranger_on_a_raw_matrix() {
        let dir = std::path::PathBuf::from(std::env::var("CR_RAW_DIR").unwrap());
        let expected_file = std::env::var("CR_CELLS").unwrap();
        let barcodes = read_gz_lines(&dir.join("barcodes.tsv.gz"));
        let feats: Vec<String> = read_gz_lines(&dir.join("features.tsv.gz"))
            .iter()
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect();
        let mut cols: Vec<Vec<(u32, u32)>> = vec![Vec::new(); barcodes.len()];
        let mut header_seen = false;
        for l in read_gz_lines(&dir.join("matrix.mtx.gz")) {
            if l.starts_with('%') {
                continue;
            }
            if !header_seen {
                header_seen = true;
                continue;
            }
            let mut it = l.split_whitespace();
            let g: u32 = it.next().unwrap().parse().unwrap();
            let c: usize = it.next().unwrap().parse().unwrap();
            let v: u32 = it.next().unwrap().parse().unwrap();
            cols[c - 1].push((g - 1, v));
        }
        let mut m = CscCounts {
            n_features: feats.len(),
            col_ptr: vec![0],
            ..Default::default()
        };
        for mut c in cols {
            c.sort_unstable();
            for (g, v) in c {
                m.rows.push(g);
                m.vals.push(v);
            }
            m.col_ptr.push(m.rows.len());
        }
        let got = call_cells(&m, &feats);
        let mut got: Vec<String> = got.iter().map(|&c| barcodes[c].clone()).collect();
        got.sort();
        let mut want: Vec<String> = std::fs::read_to_string(expected_file)
            .unwrap()
            .lines()
            .map(ToString::to_string)
            .collect();
        want.sort();
        assert_eq!(got.len(), want.len(), "called vs expected count");
        assert_eq!(got, want);
    }
}
