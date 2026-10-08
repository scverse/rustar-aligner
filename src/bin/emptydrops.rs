//! Standalone EmptyDrops_CR cell caller (Rust port of STAR
//! `SoloFeature_emptyDrops_CR.cpp` / CellRanger's EmptyDrops variant).
//!
//! Reads a raw count matrix (MatrixMarket `matrix.mtx` [.gz] genes×cells +
//! `barcodes.tsv`/`features.tsv`) and writes the called cells:
//!   - guaranteed cells from the CellRanger-2.2 knee, plus
//!   - extra cells whose expression profile is significantly different from the
//!     ambient RNA profile (multinomial Monte-Carlo test, Benjamini-Hochberg).
//!
//! Output: `<out>/barcodes.tsv` (called cells) + `<out>/cells.txt` (one called
//! barcode per line) and a `<out>/emptydrops.json` summary.
//!
//! Usage:
//!   emptydrops --raw <raw_dir> --out <dir> [--seed N] [--fdr 0.01] [--sim-n 10000]
//!
//! Defaults mirror STAR `--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000`.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use rustar_aligner::rng::{SplitMix64, cumulative_weights, sample_cumulative};

struct Args {
    raw: PathBuf,
    out: PathBuf,
    seed: u64,
    fdr: f64,
    sim_n: usize,
    n_expected: usize,
    max_percentile: f64,
    max_min_ratio: f64,
    ind_min: usize,
    ind_max: usize,
    umi_min: u64,
    umi_min_frac_median: f64,
    cand_max_n: usize,
}

fn parse_args() -> Args {
    let mut a = Args {
        raw: PathBuf::new(),
        out: PathBuf::new(),
        seed: 19_760_110,
        fdr: 0.01,
        sim_n: 10_000,
        n_expected: 3000,
        max_percentile: 0.99,
        max_min_ratio: 10.0,
        ind_min: 45_000,
        ind_max: 90_000,
        umi_min: 500,
        umi_min_frac_median: 0.01,
        cand_max_n: 20_000,
    };
    let mut it = std::env::args().skip(1);
    while let Some(k) = it.next() {
        let mut v = || it.next().expect("missing value");
        match k.as_str() {
            "--raw" => a.raw = PathBuf::from(v()),
            "--out" => a.out = PathBuf::from(v()),
            "--seed" => a.seed = v().parse().unwrap(),
            "--fdr" => a.fdr = v().parse().unwrap(),
            "--sim-n" => a.sim_n = v().parse().unwrap(),
            "--n-expected" => a.n_expected = v().parse().unwrap(),
            "--cand-max-n" => a.cand_max_n = v().parse().unwrap(),
            "--ind-min" => a.ind_min = v().parse().unwrap(),
            "--ind-max" => a.ind_max = v().parse().unwrap(),
            "--umi-min" => a.umi_min = v().parse().unwrap(),
            other => panic!("unknown arg {other}"),
        }
    }
    assert!(!a.raw.as_os_str().is_empty(), "--raw required");
    assert!(!a.out.as_os_str().is_empty(), "--out required");
    a
}

fn find(d: &Path, base: &str) -> PathBuf {
    for c in [base.to_string(), format!("{base}.gz")] {
        let p = d.join(&c);
        if p.exists() {
            return p;
        }
    }
    panic!("{base}[.gz] not found in {}", d.display());
}

fn reader(p: &Path) -> Box<dyn BufRead> {
    let f = File::open(p).unwrap();
    if p.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    }
}

fn read_lines_first_col(p: &Path) -> Vec<String> {
    reader(p)
        .lines()
        .map(|l| l.unwrap().split('\t').next().unwrap().trim().to_string())
        .collect()
}

/// Per-cell sparse profile: (gene_idx, count). Plus per-cell total.
struct Matrix {
    n_genes: usize,
    barcodes: Vec<String>,
    cell_profiles: Vec<Vec<(u32, u32)>>,
    totals: Vec<u64>,
}

fn load_matrix(raw: &Path) -> Matrix {
    let barcodes = read_lines_first_col(&find(raw, "barcodes.tsv"));
    let genes = read_lines_first_col(&find(raw, "features.tsv"));
    let n_genes = genes.len();
    let n_cells = barcodes.len();

    // MatrixMarket: skip % header, then "nGenes nCells nnz", then "gene cell count".
    let mut rd = reader(&find(raw, "matrix.mtx"));
    let mut buf = String::new();
    // header
    loop {
        buf.clear();
        rd.read_line(&mut buf).unwrap();
        if !buf.starts_with('%') {
            break;
        }
    }
    let mut cell_profiles: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n_cells];
    let mut totals = vec![0u64; n_cells];
    let mut line = String::new();
    let mut content = String::new();
    rd.read_to_string(&mut content).unwrap();
    for l in content.lines() {
        line.clear();
        let mut p = l.split_whitespace();
        let g: usize = match p.next() {
            Some(x) => x.parse().unwrap(),
            None => continue,
        };
        let c: usize = p.next().unwrap().parse().unwrap();
        let v: u64 = p.next().unwrap().parse::<f64>().unwrap() as u64;
        if v == 0 {
            continue;
        }
        let gi = g - 1;
        let ci = c - 1;
        cell_profiles[ci].push((gi as u32, v as u32));
        totals[ci] += v;
    }
    Matrix {
        n_genes,
        barcodes,
        cell_profiles,
        totals,
    }
}

/// CellRanger-2.2 knee: number of guaranteed cells (top barcodes by total).
fn knee_n_cells(sorted_desc: &[u64], n_expected: usize, max_pct: f64, max_min_ratio: f64) -> usize {
    if sorted_desc.is_empty() {
        return 0;
    }
    let idx = ((n_expected as f64 * (1.0 - max_pct)).round() as usize).min(sorted_desc.len() - 1);
    let robust_max = sorted_desc[idx] as f64;
    let thr = robust_max / max_min_ratio;
    sorted_desc.iter().take_while(|&&c| c as f64 >= thr).count()
}

fn main() {
    let a = parse_args();
    eprintln!("emptydrops: loading {}", a.raw.display());
    let m = load_matrix(&a.raw);
    let n_cells = m.totals.len();

    // Rank barcodes by total UMI, descending (stable by index for ties).
    let mut order: Vec<usize> = (0..n_cells).filter(|&i| m.totals[i] > 0).collect();
    order.sort_by(|&i, &j| m.totals[j].cmp(&m.totals[i]).then(i.cmp(&j)));
    let sorted_desc: Vec<u64> = order.iter().map(|&i| m.totals[i]).collect();

    // (1) Guaranteed cells from the CR2.2 knee.
    let n_simple = knee_n_cells(
        &sorted_desc,
        a.n_expected,
        a.max_percentile,
        a.max_min_ratio,
    );
    eprintln!("emptydrops: {n_simple} guaranteed cells from CR2.2 knee");

    // (2) Ambient profile from rank [ind_min, ind_max).
    let mut amb = vec![0f64; m.n_genes];
    let mut amb_total = 0f64;
    for &cell in order
        .iter()
        .skip(a.ind_min)
        .take(a.ind_max.saturating_sub(a.ind_min))
    {
        for &(g, c) in &m.cell_profiles[cell] {
            amb[g as usize] += c as f64;
            amb_total += c as f64;
        }
    }
    if amb_total == 0.0 {
        eprintln!("emptydrops: empty ambient range — falling back to knee-only");
        write_out(&a, &m, &order[..n_simple], n_simple, 0);
        return;
    }
    // Good-Turing P0 (unseen mass) distributed over zero-count genes; seen genes
    // get proportional mass scaled by (1 - P0). Approximates STAR's SGT.
    let n1 = amb.iter().filter(|&&x| (x - 1.0).abs() < 0.5).count() as f64;
    let p0 = (n1 / amb_total).clamp(1e-12, 0.5);
    let n_zero = amb.iter().filter(|&&x| x == 0.0).count().max(1) as f64;
    let amb_prob: Vec<f64> = amb
        .iter()
        .map(|&x| {
            if x > 0.0 {
                (1.0 - p0) * x / amb_total
            } else {
                p0 / n_zero
            }
        })
        .collect();
    let amb_logp: Vec<f64> = amb_prob.iter().map(|&p| p.max(1e-300).ln()).collect();

    // (3) Candidate barcodes: rank >= n_simple, total >= minUMI, up to cand_max_n.
    let median_top = if n_simple >= 2 {
        sorted_desc[n_simple / 2]
    } else if !sorted_desc.is_empty() {
        sorted_desc[0]
    } else {
        0
    };
    let min_umi = a
        .umi_min
        .max((a.umi_min_frac_median * median_top as f64) as u64);
    let mut cands: Vec<usize> = Vec::new();
    for &cell in order.iter().skip(n_simple).take(a.cand_max_n) {
        if m.totals[cell] < min_umi {
            break;
        }
        cands.push(cell);
    }
    eprintln!(
        "emptydrops: {} candidates (minUMI={min_umi}); running {} Monte-Carlo sims",
        cands.len(),
        a.sim_n
    );
    if cands.is_empty() {
        write_out(&a, &m, &order[..n_simple], n_simple, 0);
        return;
    }

    // logFactorial up to the largest candidate total.
    let max_count = cands.iter().map(|&c| m.totals[c]).max().unwrap() as usize;
    let mut log_fac = vec![0f64; max_count + 1];
    for i in 2..=max_count {
        log_fac[i] = log_fac[i - 1] + (i as f64).ln();
    }

    // Observed multinomial log-prob per candidate.
    let obs_logp: Vec<f64> = cands
        .iter()
        .map(|&cell| {
            let total = m.totals[cell] as usize;
            let mut s = log_fac[total];
            for &(g, c) in &m.cell_profiles[cell] {
                s -= log_fac[c as usize];
                s += c as f64 * amb_logp[g as usize];
            }
            s
        })
        .collect();

    // (4/5) Monte Carlo: simulate sim_n barcodes from the ambient multinomial,
    // recording the running log-prob at every count up to max_count. Each
    // candidate of total t is compared against sim[*][t].
    let nonzero: Vec<usize> = (0..m.n_genes).filter(|&g| amb_prob[g] > 0.0).collect();
    let weights: Vec<f64> = nonzero.iter().map(|&g| amb_prob[g]).collect();
    let cumulative = cumulative_weights(&weights);
    let mut rng = SplitMix64::seed(a.seed);

    // For each count t, collect the sim log-probs (so we can compare per candidate).
    // Memory: sim_n * (max_count+1) f64 — fine for ~10k * a few-thousand.
    let mut sim_at: Vec<Vec<f64>> = vec![Vec::with_capacity(a.sim_n); max_count + 1];
    let mut curr = vec![0u32; m.n_genes];
    for _ in 0..a.sim_n {
        for v in curr.iter_mut() {
            *v = 0;
        }
        let mut lp = 0f64;
        sim_at[0].push(0.0);
        #[allow(clippy::needless_range_loop)] // ic is both index and multinomial term
        for ic in 1..=max_count {
            let gi = nonzero[sample_cumulative(&cumulative, &mut rng)];
            curr[gi] += 1;
            lp += amb_logp[gi] + (ic as f64).ln() - (curr[gi] as f64).ln();
            sim_at[ic].push(lp);
        }
    }

    // p-value: fraction of sims with LOWER log-prob than observed (more extreme).
    let mut pvals: Vec<(usize, f64)> = cands
        .iter()
        .enumerate()
        .map(|(i, &cell)| {
            let t = m.totals[cell] as usize;
            let obs = obs_logp[i];
            let n_lower = sim_at[t].iter().filter(|&&sp| sp < obs).count();
            let p = (1 + n_lower) as f64 / (1 + a.sim_n) as f64;
            (i, p)
        })
        .collect();

    // (6) Benjamini-Hochberg.
    pvals.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap());
    let n = pvals.len() as f64;
    let mut padj = vec![0f64; pvals.len()];
    for (rank, &(_, p)) in pvals.iter().enumerate() {
        padj[rank] = (p * n / (rank + 1) as f64).min(1.0);
    }
    for i in (0..padj.len() - 1).rev() {
        padj[i] = padj[i].min(padj[i + 1]);
    }

    // Called cells = guaranteed + candidates with padj <= FDR.
    let mut called: Vec<usize> = order[..n_simple].to_vec();
    let mut extra = 0usize;
    for (rank, &(ci, _)) in pvals.iter().enumerate() {
        if padj[rank] <= a.fdr {
            called.push(cands[ci]);
            extra += 1;
        }
    }
    eprintln!("emptydrops: {extra} extra cells (FDR<={})", a.fdr);
    write_out(&a, &m, &called, n_simple, extra);
}

fn write_out(a: &Args, m: &Matrix, called: &[usize], n_simple: usize, extra: usize) {
    std::fs::create_dir_all(&a.out).unwrap();
    // Stable order: by descending total then barcode.
    let mut cells: Vec<usize> = called.to_vec();
    cells.sort_by(|&i, &j| {
        m.totals[j]
            .cmp(&m.totals[i])
            .then(m.barcodes[i].cmp(&m.barcodes[j]))
    });
    cells.dedup();

    let mut bc = BufWriter::new(File::create(a.out.join("barcodes.tsv")).unwrap());
    let mut cl = BufWriter::new(File::create(a.out.join("cells.txt")).unwrap());
    for &c in &cells {
        writeln!(bc, "{}", m.barcodes[c]).unwrap();
        writeln!(cl, "{}", m.barcodes[c]).unwrap();
    }
    let summary = format!(
        "{{\"n_cells\": {}, \"n_guaranteed\": {}, \"n_emptydrops_extra\": {}, \"fdr\": {}, \"sim_n\": {}}}\n",
        cells.len(),
        n_simple,
        extra,
        a.fdr,
        a.sim_n
    );
    std::fs::write(a.out.join("emptydrops.json"), &summary).unwrap();
    println!(
        "EmptyDrops_CR: {} cells ({} guaranteed + {} EmptyDrops) -> {}",
        cells.len(),
        n_simple,
        extra,
        a.out.display()
    );
}
