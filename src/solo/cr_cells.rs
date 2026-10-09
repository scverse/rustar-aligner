//! CellRanger's cell calling (`--soloOutLayout CellRanger`): the `ordmag` initial
//! call followed by the EmptyDrops-style rescue of barcodes whose profile is
//! not the ambient one, as `cellranger count` 10.0.0 runs them for a 3' v3
//! library (`filter_barcodes_method ordmag_nonambient`).
//!
//! Ported from `lib/python/cellranger/cell_calling_helpers.py`
//! (`filter_cellular_barcodes_ordmag`), `cell_calling.py`
//! (`find_nonambient_barcodes`), `stats.py` and `sgt.py`. CellRanger draws its
//! random numbers from NumPy, and the calls it makes depend on them, so the two
//! generators are reproduced bit for bit: `RandomState(0)` (MT19937, masked
//! rejection sampling) for the ordmag bootstrap and `default_rng(42)` (PCG64
//! seeded through `SeedSequence`) for the ambient simulations.

use crate::solo::libcxx_rng::Mt19937;
use rayon::prelude::*;

/// A feature-by-barcode count matrix in column (barcode) order.
pub struct CscCounts {
    pub n_features: usize,
    /// `col_ptr[c]..col_ptr[c + 1]` indexes `rows`/`vals` for column `c`.
    pub col_ptr: Vec<usize>,
    pub rows: Vec<u32>,
    pub vals: Vec<u32>,
}

impl CscCounts {
    pub fn n_cols(&self) -> usize {
        self.col_ptr.len() - 1
    }

    pub fn col(&self, c: usize) -> (&[u32], &[u32]) {
        let (a, b) = (self.col_ptr[c], self.col_ptr[c + 1]);
        (&self.rows[a..b], &self.vals[a..b])
    }

    pub fn umis(&self) -> Vec<u64> {
        (0..self.n_cols())
            .map(|c| self.col(c).1.iter().map(|&v| u64::from(v)).sum())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// NumPy generators
// ---------------------------------------------------------------------------

/// `numpy.random.RandomState(0).choice(a, size)`'s index draws: `randint(0, n)`
/// by masked rejection on 32-bit draws.
struct LegacyRandomState {
    mt: Mt19937,
}

impl LegacyRandomState {
    fn new(seed: u32) -> Self {
        Self {
            mt: Mt19937::new(seed),
        }
    }

    fn randint(&mut self, n: u64) -> u64 {
        let rng = n - 1;
        if rng == 0 {
            return 0;
        }
        let mask = u64::MAX >> rng.leading_zeros();
        loop {
            let v = u64::from(self.mt.next_u32()) & mask;
            if v <= rng {
                return v;
            }
        }
    }
}

const MULT_128: u128 = 0x2360_ED05_1FC6_5DA4_4385_DF64_9FCC_F645;

/// `numpy.random.default_rng(seed)` (PCG64, XSL-RR 128/64).
#[derive(Clone)]
struct Pcg64 {
    state: u128,
    inc: u128,
}

fn seed_sequence_state(entropy: u32) -> [u64; 4] {
    const INIT_A: u32 = 0x43b0_d7e5;
    const MULT_A: u32 = 0x931e_8875;
    const INIT_B: u32 = 0x8b51_f9dd;
    const MULT_B: u32 = 0x58f3_8ded;
    const MIX_L: u32 = 0xca01_f9dd;
    const MIX_R: u32 = 0x4973_f715;
    let mut hash_const = INIT_A;
    let hashmix = |value: u32, hc: &mut u32| -> u32 {
        let mut v = value ^ *hc;
        *hc = hc.wrapping_mul(MULT_A);
        v = v.wrapping_mul(*hc);
        v ^ (v >> 16)
    };
    let mix = |x: u32, y: u32| -> u32 {
        let r = MIX_L.wrapping_mul(x).wrapping_sub(MIX_R.wrapping_mul(y));
        r ^ (r >> 16)
    };
    let mut pool = [0u32; 4];
    pool[0] = hashmix(entropy, &mut hash_const);
    for p in pool.iter_mut().skip(1) {
        *p = hashmix(0, &mut hash_const);
    }
    for i_src in 0..4 {
        for i_dst in 0..4 {
            if i_src != i_dst {
                let h = hashmix(pool[i_src], &mut hash_const);
                pool[i_dst] = mix(pool[i_dst], h);
            }
        }
    }
    let mut hc = INIT_B;
    let mut words = [0u32; 8];
    for (i, w) in words.iter_mut().enumerate() {
        let mut d = pool[i % 4] ^ hc;
        hc = hc.wrapping_mul(MULT_B);
        d = d.wrapping_mul(hc);
        *w = d ^ (d >> 16);
    }
    let mut out = [0u64; 4];
    for k in 0..4 {
        out[k] = u64::from(words[2 * k]) | u64::from(words[2 * k + 1]) << 32;
    }
    out
}

impl Pcg64 {
    fn new(seed: u32) -> Self {
        let v = seed_sequence_state(seed);
        let initstate = u128::from(v[0]) << 64 | u128::from(v[1]);
        let initseq = u128::from(v[2]) << 64 | u128::from(v[3]);
        let mut r = Self {
            state: 0,
            inc: initseq << 1 | 1,
        };
        r.step();
        r.state = r.state.wrapping_add(initstate);
        r.step();
        r
    }

    fn step(&mut self) {
        self.state = self.state.wrapping_mul(MULT_128).wrapping_add(self.inc);
    }

    fn next_u64(&mut self) -> u64 {
        self.step();
        let s = self.state;
        let xsl = (s >> 64) as u64 ^ s as u64;
        xsl.rotate_right((s >> 122) as u32)
    }

    /// `Generator.random()`: 53 random bits.
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Jump `delta` draws ahead.
    fn advance(&mut self, mut delta: u128) {
        let (mut acc_mult, mut acc_plus) = (1u128, 0u128);
        let (mut cur_mult, mut cur_plus) = (MULT_128, self.inc);
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
// ordmag
// ---------------------------------------------------------------------------

const ORDMAG_BOOTSTRAPS: usize = 100;
const ORDMAG_QUANTILE: f64 = 0.99;
const MIN_RECOVERED: i64 = 50;

/// `find_within_ordmag` for each baseline index, on an ascending sample.
fn find_within_ordmag(asc: &[u64], baseline_idx: usize) -> usize {
    let baseline = asc[asc.len() - 1 - baseline_idx];
    let cutoff = ((0.1 * baseline as f64).round_ties_even() as u64).max(1);
    asc.len() - asc.partition_point(|&v| v < cutoff)
}

fn bootstrap_sorted(nonzero: &[u64], rs: &mut LegacyRandomState) -> Vec<u64> {
    let n = nonzero.len() as u64;
    let mut s: Vec<u64> = (0..nonzero.len())
        .map(|_| nonzero[rs.randint(n) as usize])
        .collect();
    s.sort_unstable();
    s
}

/// `estimate_recovered_cells_ordmag`: `(recovered_cells, loss)`.
fn estimate_recovered(asc: &[u64], max_expected: usize) -> (f64, f64) {
    // unique(round(2 ** linspace(1, log2(max), 2000)))
    let (start, stop) = (1.0f64, (max_expected as f64).log2());
    let step = (stop - start) / 1999.0;
    let mut rcs: Vec<i64> = (0..2000)
        .map(|i| {
            let y = if i == 1999 {
                stop
            } else {
                start + i as f64 * step
            };
            y.exp2().round_ties_even() as i64
        })
        .collect();
    rcs.dedup();
    let mut best = (f64::INFINITY, 0i64);
    for &rc in &rcs {
        let b = ((rc as f64) * (1.0 - ORDMAG_QUANTILE)).round_ties_even() as i64;
        let b = (b.max(0) as usize).min(asc.len() - 1);
        let filtered = find_within_ordmag(asc, b) as f64;
        let loss = (filtered - rc as f64).powi(2) / rc as f64;
        if loss < best.0 {
            best = (loss, rc);
        }
    }
    (best.1 as f64, best.0)
}

/// `filter_cellular_barcodes_ordmag` with `recovered_cells=None` and a chemistry
/// whose empty-drops range starts at 45 000: the barcode indices called.
pub fn ordmag_cells(bc_counts: &[u64]) -> Vec<usize> {
    let nonzero: Vec<u64> = bc_counts.iter().copied().filter(|&c| c > 0).collect();
    if nonzero.is_empty() {
        return Vec::new();
    }
    let mut rs = LegacyRandomState::new(0);
    // min(empty-drops range start, 2^18)
    let max_expected = 45_000usize;
    let mut sum = (0.0f64, 0.0f64);
    for _ in 0..ORDMAG_BOOTSTRAPS {
        let s = bootstrap_sorted(&nonzero, &mut rs);
        let (rc, loss) = estimate_recovered(&s, max_expected);
        sum.0 += rc;
        sum.1 += loss;
    }
    let recovered =
        ((sum.0 / ORDMAG_BOOTSTRAPS as f64).round_ties_even() as i64).max(MIN_RECOVERED);
    let baseline = ((recovered as f64 * (1.0 - ORDMAG_QUANTILE)).round_ties_even() as usize)
        .min(nonzero.len() - 1);
    let mut top_n_sum = 0u64;
    for _ in 0..ORDMAG_BOOTSTRAPS {
        let s = bootstrap_sorted(&nonzero, &mut rs);
        top_n_sum += find_within_ordmag(&s, baseline) as u64;
    }
    let mean = top_n_sum as f64 / ORDMAG_BOOTSTRAPS as f64;
    let mut nbcs = mean.round_ties_even() as usize;
    if nbcs > 0 {
        // Take every barcode tied with the last one selected, unless that grabs
        // more than 20% extra.
        let mut sorted: Vec<u64> = nonzero.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        let cutoff = sorted[nbcs - 1];
        let mut index = nbcs - 1;
        while index + 1 < sorted.len() && sorted[index] == cutoff {
            index += 1;
            if (index + 1 - nbcs) as f64 > 0.20 * nbcs as f64 {
                break;
            }
            nbcs = index + 1;
        }
    }
    // top n by count; ties resolved as a stable ascending sort read backwards.
    let mut order: Vec<usize> = (0..bc_counts.len()).collect();
    order.sort_by_key(|&i| bc_counts[i]);
    order.reverse();
    let mut top: Vec<usize> = order.into_iter().take(nbcs).collect();
    top.sort_unstable();
    top
}

// ---------------------------------------------------------------------------
// Simple Good-Turing (sgt.py)
// ---------------------------------------------------------------------------

/// `sgt_proportions` of strictly positive integer frequencies.
#[allow(
    clippy::manual_midpoint,
    clippy::many_single_char_names,
    clippy::float_cmp
)]
fn sgt_proportions(freqs: &[u64]) -> Option<(Vec<f64>, f64)> {
    let max = *freqs.iter().max()? as usize;
    let mut ff = vec![0u64; max + 1];
    for &f in freqs {
        ff[f as usize] += 1;
    }
    let use_freqs: Vec<usize> = (1..=max).filter(|&r| ff[r] > 0).collect();
    if use_freqs.len() < 10 {
        return None;
    }
    let xr: Vec<f64> = use_freqs.iter().map(|&r| r as f64).collect();
    let xnr: Vec<f64> = use_freqs.iter().map(|&r| ff[r] as f64).collect();
    let n = xr.len();
    let xn: f64 = xr.iter().zip(&xnr).map(|(a, b)| a * b).sum();

    // averaging transform
    let mut d = vec![1.0f64; n];
    for i in 1..n {
        d[i] = xr[i] - xr[i - 1];
    }
    let mut dr = vec![0.0f64; n];
    for i in 0..n - 1 {
        dr[i] = 0.5 * (d[i + 1] + d[i]);
    }
    dr[n - 1] = d[n - 1];
    let xnrz: Vec<f64> = (0..n).map(|i| xnr[i] / dr[i]).collect();

    // linear regression of log(xnrz) on log(xr)
    let lx: Vec<f64> = xr.iter().map(|v| v.ln()).collect();
    let ly: Vec<f64> = xnrz.iter().map(|v| v.ln()).collect();
    let mx = lx.iter().sum::<f64>() / n as f64;
    let my = ly.iter().sum::<f64>() / n as f64;
    let sxy: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = lx.iter().map(|a| (a - mx) * (a - mx)).sum();
    let mut slope = sxy / sxx;
    if slope >= -1.0 {
        slope = -1.0;
    }

    let xrstrel: Vec<f64> = xr
        .iter()
        .map(|&r| r * (1.0 + 1.0 / r).powf(1.0 + slope) / r)
        .collect();
    let mut xrtry = vec![false; n];
    let mut xrstarel = vec![0.0f64; n];
    for i in 0..n {
        let next = if i + 1 < n { xr[i + 1] - 1.0 } else { 0.0 };
        if xr[i] == next {
            xrtry[i] = true;
            xrstarel[i] = (xr[i] + 1.0) / xr[i] * xnr[i + 1] / xnr[i];
        }
    }
    let mut tursd = vec![1.0f64; n];
    for i in 0..n {
        if xrtry[i] {
            tursd[i] =
                (i as f64 + 2.0) / xnr[i] * (xnr[i + 1] * (1.0 + xnr[i + 1] / xnr[i])).sqrt();
        }
    }
    let mut comb = vec![0.0f64; n];
    let mut use_turing = true;
    for r in 0..n {
        if !use_turing {
            comb[r] = xrstrel[r];
        } else if (xrstrel[r] - xrstarel[r]).abs() * (1.0 + r as f64) / tursd[r] > 1.65 {
            comb[r] = xrstarel[r];
        } else {
            use_turing = false;
            comb[r] = xrstrel[r];
        }
    }
    let sumpraw: f64 = (0..n).map(|i| comb[i] * xr[i] * xnr[i] / xn).sum();
    let p0 = xnr[0] / xn;
    let rstar: Vec<f64> = (0..n)
        .map(|i| xr[i] * (comb[i] * (1.0 - xnr[0] / xn) / sumpraw))
        .collect();

    let rstar_sum: f64 = (0..n).map(|i| xnr[i] * rstar[i]).sum();
    let mut by_freq = vec![0.0f64; max + 1];
    for (k, &r) in use_freqs.iter().enumerate() {
        by_freq[r] = rstar[k];
    }
    let pstar = freqs
        .iter()
        .map(|&f| (1.0 - p0) * (by_freq[f as usize] / rstar_sum))
        .collect();
    Some((pstar, p0))
}

// ---------------------------------------------------------------------------
// Non-ambient barcodes (cell_calling.py)
// ---------------------------------------------------------------------------

const EMPTY_LOW: usize = 45_000;
const EMPTY_HIGH: usize = 90_000;
const MIN_UMIS: u64 = 500;
const NUM_SIMS: usize = 100_000;
const FDR: f64 = 0.001;

/// `scipy.special.gammaln` at the integer `n`.
fn ln_factorial_table(max: usize) -> Vec<f64> {
    // gammaln(k + 1) for k in 0..=max, by cephes' own recurrences: an exact
    // product below 13, Stirling's series above.
    (0..=max)
        .map(|k| {
            let x = (k + 1) as f64;
            if x < 13.0 {
                let mut z = 1.0f64;
                let mut u = x;
                let mut p = 0.0f64;
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
                z.ln()
            } else {
                let mut q = (x - 0.5) * x.ln() - x + 0.918_938_533_204_672_7;
                let p = 1.0 / (x * x);
                if x >= 1000.0 {
                    q += ((7.936_507_936_507_937e-4 * p - 2.777_777_777_777_778e-3) * p
                        + 0.083_333_333_333_333_33)
                        / x;
                } else {
                    const A: [f64; 5] = [
                        8.116_141_674_705_085e-4,
                        -5.950_619_042_843_014e-4,
                        7.936_503_404_577_169e-4,
                        -2.777_777_777_300_997e-3,
                        8.333_333_333_333_319e-2,
                    ];
                    let mut acc = A[0];
                    for c in &A[1..] {
                        acc = acc * p + c;
                    }
                    q += acc / x;
                }
                q
            }
        })
        .collect()
}

/// `numpy.argsort(v)` with its default kind (`npy_aquicksort`): introsort with
/// a median-of-three partition and insertion sort below 16 elements. It is not
/// stable, and which of several equal values lands where decides which
/// barcodes sit at a rank boundary, so the exact algorithm is reproduced.
#[allow(clippy::many_single_char_names)]
fn numpy_argsort(v: &[u64]) -> Vec<usize> {
    const SMALL: isize = 16;
    let n = v.len();
    let mut idx: Vec<usize> = (0..n).collect();
    if n < 2 {
        return idx;
    }
    let msb = |x: usize| (usize::BITS - 1 - x.leading_zeros()) as i32;
    let mut stack: Vec<(isize, isize, i32)> = Vec::new();
    let (mut pl, mut pr) = (0isize, n as isize - 1);
    let mut cdepth = msb(n) * 2;
    loop {
        let mut descend = true;
        if cdepth < 0 {
            // Heapsort fallback (never reached on real data): any correct sort.
            idx[pl as usize..=pr as usize].sort_by_key(|&i| v[i]);
            descend = false;
        }
        if descend {
            while pr - pl > SMALL {
                let pm = pl + ((pr - pl) >> 1);
                if v[idx[pm as usize]] < v[idx[pl as usize]] {
                    idx.swap(pm as usize, pl as usize);
                }
                if v[idx[pr as usize]] < v[idx[pm as usize]] {
                    idx.swap(pr as usize, pm as usize);
                }
                if v[idx[pm as usize]] < v[idx[pl as usize]] {
                    idx.swap(pm as usize, pl as usize);
                }
                let vp = v[idx[pm as usize]];
                let mut pi = pl;
                let mut pj = pr - 1;
                idx.swap(pm as usize, pj as usize);
                loop {
                    loop {
                        pi += 1;
                        if v[idx[pi as usize]] >= vp {
                            break;
                        }
                    }
                    loop {
                        pj -= 1;
                        if vp >= v[idx[pj as usize]] {
                            break;
                        }
                    }
                    if pi >= pj {
                        break;
                    }
                    idx.swap(pi as usize, pj as usize);
                }
                let pk = pr - 1;
                idx.swap(pi as usize, pk as usize);
                if pi - pl < pr - pi {
                    stack.push((pi + 1, pr, 0));
                    pr = pi - 1;
                } else {
                    stack.push((pl, pi - 1, 0));
                    pl = pi + 1;
                }
                cdepth -= 1;
                if let Some(top) = stack.last_mut() {
                    top.2 = cdepth;
                }
            }
            // insertion sort
            for pi in pl + 1..=pr {
                let vi = idx[pi as usize];
                let vp = v[vi];
                let mut pj = pi;
                while pj > pl && vp < v[idx[(pj - 1) as usize]] {
                    idx[pj as usize] = idx[(pj - 1) as usize];
                    pj -= 1;
                }
                idx[pj as usize] = vi;
            }
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
    idx
}

/// `find_nonambient_barcodes`: the barcode indices (not in `orig_cells`) whose
/// profile differs from the ambient one at the CellRanger FDR.
pub fn nonambient_cells(m: &CscCounts, feature_ids: &[String], orig_cells: &[usize]) -> Vec<usize> {
    let umis = m.umis();
    let n_bc = umis.len();
    // The barcodes of empty partitions: ranks [45000, 90000) from the top.
    let mut order = numpy_argsort(&umis);
    order.reverse();
    let lo = EMPTY_LOW.min(n_bc);
    let hi = EMPTY_HIGH.min(n_bc);
    let mut empty: Vec<usize> = order[lo..hi].to_vec();
    empty.sort_unstable();
    let ambient: Vec<usize> = empty.iter().copied().filter(|&i| umis[i] > 0).collect();
    if ambient.is_empty() || orig_cells.is_empty() {
        return Vec::new();
    }

    // Features with any count, ordered by id.
    let mut any = vec![false; m.n_features];
    for &r in &m.rows {
        any[r as usize] = true;
    }
    let mut eval_features: Vec<usize> = (0..m.n_features).filter(|&f| any[f]).collect();
    eval_features.sort_by(|&a, &b| feature_ids[a].cmp(&feature_ids[b]));
    let mut pos_of = vec![usize::MAX; m.n_features];
    for (k, &f) in eval_features.iter().enumerate() {
        pos_of[f] = k;
    }

    // Ambient profile by Simple Good-Turing.
    let mut profile = vec![0u64; eval_features.len()];
    for &c in &ambient {
        let (rows, vals) = m.col(c);
        for (&r, &v) in rows.iter().zip(vals) {
            profile[pos_of[r as usize]] += u64::from(v);
        }
    }
    let nz: Vec<usize> = (0..profile.len()).filter(|&i| profile[i] > 0).collect();
    let nzf: Vec<u64> = nz.iter().map(|&i| profile[i]).collect();
    let Some((p_smoothed, p0)) = sgt_proportions(&nzf) else {
        return Vec::new();
    };
    let n0 = profile.len() - nz.len();
    let mut prof_p = vec![if n0 == 0 { -1.0 } else { p0 / n0 as f64 }; profile.len()];
    if n0 == 0 {
        let total: f64 = p_smoothed.iter().sum();
        for (k, &i) in nz.iter().enumerate() {
            prof_p[i] = p_smoothed[k] / total;
        }
    } else {
        for (k, &i) in nz.iter().enumerate() {
            prof_p[i] = p_smoothed[k];
        }
    }

    // Candidates: non-cells above the minimum, ascending by index.
    let mut is_orig = vec![false; n_bc];
    for &i in orig_cells {
        is_orig[i] = true;
    }
    let max_bg = empty.iter().map(|&i| umis[i]).max().unwrap_or(0);
    let min_umis = MIN_UMIS.max(1 + max_bg);
    let eval_bcs: Vec<usize> = (0..n_bc)
        .filter(|&i| !is_orig[i] && umis[i] >= min_umis)
        .collect();
    if eval_bcs.is_empty() {
        return Vec::new();
    }
    log::info!(
        "cell calling: {} candidate barcodes, {} ambient, {} original cells",
        eval_bcs.len(),
        ambient.len(),
        orig_cells.len()
    );

    let max_n = eval_bcs.iter().map(|&i| umis[i]).max().unwrap() as usize;
    let lgam = ln_factorial_table(max_n.max(2));
    let logp: Vec<f64> = prof_p.iter().map(|p| p.ln()).collect();

    // Observed log-likelihoods.
    let obs: Vec<f64> = eval_bcs
        .iter()
        .map(|&c| {
            let (rows, vals) = m.col(c);
            let mut s = 0.0f64;
            let mut g = 0.0f64;
            for (&r, &v) in rows.iter().zip(vals) {
                g += lgam[v as usize];
                s += f64::from(v) * logp[pos_of[r as usize]];
            }
            lgam[umis[c] as usize] - g + s
        })
        .collect();
    let n_of: Vec<usize> = eval_bcs.iter().map(|&i| umis[i] as usize).collect();

    // Simulated log-likelihoods: sim `s` is the next `max_n` draws of one
    // generator, so chunks can start where the stream would have got to.
    let mut p_cum = Vec::with_capacity(prof_p.len());
    let mut acc = 0.0f64;
    for p in &prof_p {
        acc += p;
        p_cum.push(acc);
    }
    let ln_n: Vec<f64> = (1..=max_n).map(|j| (j as f64).ln()).collect();
    let chunk = 256usize;
    let n_chunks = NUM_SIMS.div_ceil(chunk);
    let lower: Vec<u32> = (0..n_chunks)
        .into_par_iter()
        .map(|ci| {
            let (s0, s1) = (ci * chunk, ((ci + 1) * chunk).min(NUM_SIMS));
            let mut rng = Pcg64::new(42);
            rng.advance(s0 as u128 * max_n as u128);
            let mut counts = vec![0u32; prof_p.len()];
            let mut lo = vec![0u32; eval_bcs.len()];
            let mut cum = vec![0.0f64; max_n];
            for _ in s0..s1 {
                let mut total = 0.0f64;
                for j in 0..max_n {
                    let u = rng.next_f64();
                    let f = p_cum.partition_point(|&c| c < u).min(prof_p.len() - 1);
                    counts[f] += 1;
                    let term = (ln_n[j] - f64::from(counts[f]).ln()) + logp[f];
                    total += term;
                    cum[j] = total;
                }
                counts.fill(0);
                for (i, &n) in n_of.iter().enumerate() {
                    if cum[n - 1] < obs[i] {
                        lo[i] += 1;
                    }
                }
            }
            lo
        })
        .reduce(
            || vec![0u32; eval_bcs.len()],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b) {
                    *x += y;
                }
                a
            },
        );
    // `reduce` above sums the per-chunk vectors, so `lower` is already per barcode.
    let pvalues: Vec<f64> = lower
        .iter()
        .map(|&l| (1.0 + f64::from(l)) / (1.0 + NUM_SIMS as f64))
        .collect();

    // Benjamini-Hochberg.
    let n = pvalues.len();
    let mut asc: Vec<usize> = (0..n).collect();
    asc.sort_by(|&a, &b| pvalues[a].partial_cmp(&pvalues[b]).unwrap());
    let mut adj = vec![0.0f64; n];
    let mut run = f64::INFINITY;
    for (k, &i) in asc.iter().enumerate().rev() {
        let scale = n as f64 / (k + 1) as f64;
        run = run.min(scale * pvalues[i]);
        adj[i] = run.min(1.0);
    }
    eval_bcs
        .iter()
        .zip(&adj)
        .filter(|&(_, &a)| a <= FDR)
        .map(|(&i, _)| i)
        .collect()
}

/// The cells `cellranger count` calls: `ordmag`, plus the non-ambient
/// barcodes, ascending.
pub fn call_cells(m: &CscCounts, feature_ids: &[String]) -> Vec<usize> {
    let umis = m.umis();
    let mut cells = ordmag_cells(&umis);
    let extra = nonambient_cells(m, feature_ids, &cells);
    log::info!(
        "cell calling: {} ordmag cells + {} non-ambient",
        cells.len(),
        extra.len()
    );
    cells.extend(extra);
    cells.sort_unstable();
    cells.dedup();
    cells
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values from NumPy 2: `default_rng(42).random(3)`.
    #[test]
    fn pcg64_matches_numpy_default_rng_42() {
        let mut r = Pcg64::new(42);
        let a = [r.next_f64(), r.next_f64(), r.next_f64()];
        assert!((a[0] - 0.773_956_048_555_963_3).abs() < 1e-15, "{a:?}");
        assert!((a[1] - 0.438_878_439_752_052_3).abs() < 1e-15, "{a:?}");
        assert!((a[2] - 0.858_597_919_911_382_5).abs() < 1e-15, "{a:?}");
    }

    #[test]
    fn advancing_equals_drawing() {
        let mut a = Pcg64::new(42);
        for _ in 0..1000 {
            a.next_u64();
        }
        let mut b = Pcg64::new(42);
        b.advance(1000);
        assert_eq!(a.next_u64(), b.next_u64());
    }

    /// `RandomState(0).randint(0, 10, 3)` is `[5, 0, 3]` in NumPy.
    #[test]
    fn legacy_randint_matches_numpy() {
        let mut r = LegacyRandomState::new(0);
        let v = [r.randint(10), r.randint(10), r.randint(10)];
        assert_eq!(v, [5, 0, 3]);
    }

    #[test]
    fn ln_factorial_matches_the_exact_values() {
        let t = ln_factorial_table(20);
        assert!((t[5] - 120f64.ln()).abs() < 1e-14);
        assert!((t[20] - 2_432_902_008_176_640_000f64.ln()).abs() < 1e-12);
        assert_eq!(t[0], 0.0);
    }

    /// Run with `CR_RAW_DIR=<dir with matrix.mtx.gz, features.tsv.gz,
    /// barcodes.tsv.gz> CR_CELLS=<expected barcodes, one per line> cargo test
    /// --release -- --ignored cell_calling_matches`.
    #[test]
    #[ignore = "needs a raw matrix and CellRanger's cell list"]
    fn cell_calling_matches_cellranger_on_a_raw_matrix() {
        use std::io::{BufRead, BufReader};
        let dir = std::env::var("CR_RAW_DIR").unwrap();
        let open = |n: &str| {
            BufReader::new(flate2::read::MultiGzDecoder::new(
                std::fs::File::open(format!("{dir}/{n}")).unwrap(),
            ))
        };
        let feats: Vec<String> = open("features.tsv.gz")
            .lines()
            .map(|l| l.unwrap().split('\t').next().unwrap().to_string())
            .collect();
        let bcs: Vec<String> = open("barcodes.tsv.gz")
            .lines()
            .map(|l| l.unwrap())
            .collect();
        let mut cols: Vec<Vec<(u32, u32)>> = vec![Vec::new(); bcs.len()];
        let mut header = true;
        for l in open("matrix.mtx.gz").lines() {
            let l = l.unwrap();
            if l.starts_with('%') {
                continue;
            }
            if header {
                header = false;
                continue;
            }
            let mut it = l.split(' ');
            let g: u32 = it.next().unwrap().parse().unwrap();
            let c: usize = it.next().unwrap().parse().unwrap();
            let v: u32 = it.next().unwrap().parse().unwrap();
            cols[c - 1].push((g - 1, v));
        }
        let (mut col_ptr, mut rows, mut vals) = (vec![0usize], Vec::new(), Vec::new());
        for mut c in cols {
            c.sort_unstable();
            for (g, v) in c {
                rows.push(g);
                vals.push(v);
            }
            col_ptr.push(rows.len());
        }
        let m = CscCounts {
            n_features: feats.len(),
            col_ptr,
            rows,
            vals,
        };
        let got: Vec<&str> = call_cells(&m, &feats)
            .into_iter()
            .map(|c| bcs[c].as_str())
            .collect();
        let want = std::fs::read_to_string(std::env::var("CR_CELLS").unwrap()).unwrap();
        let mut want: Vec<&str> = want.lines().collect();
        want.sort_unstable();
        let mut got = got;
        got.sort_unstable();
        assert_eq!(got.len(), want.len());
        assert_eq!(got, want);
    }
}
