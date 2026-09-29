//! `--soloFeatures Transcript3p`: quantify transcripts rather than genes, using
//! how far each read's 3' end sits from the transcript's 3' end.
//!
//! STAR `Transcriptome_classifyAlign.cpp` plus `SoloFeature_quantTranscript.cpp`.
//!
//! **STAR marks this feature as under development.** In `parametersDefault`
//! both `Transcript3p` and `--soloClusterCBfile` sit between
//! `#####UnderDevelopment_begin : not supported - do not use` and
//! `#####UnderDevelopment_end`, and `STAR --help` prints that banner around
//! them. The port follows STAR's code, so it inherits whatever that code does,
//! including its unfinished parts; it is not a stable interface either here or
//! upstream, and a differential against STAR compares two implementations of
//! something STAR itself does not support.
//!
//! In a 3'-biased assay every read lands near the transcript's 3' end, and how
//! near is informative: a read 200 bases from the end of one isoform and 4000
//! from the end of another is evidence for the first. This feature records, per
//! read, every transcript the alignment is concordant with and the spliced
//! distance from the read to that transcript's 3' end. The distribution of
//! those distances is then estimated from the data itself and used as the
//! likelihood in an EM over UMIs.
//!
//! Two things make it different from the gene features. The output is per
//! *cluster* rather than per cell (`--soloClusterCBfile` says which cell is in
//! which cluster), because a single cell has too few UMIs to run an EM over
//! isoforms. And a UMI seen on several reads contributes the *intersection* of
//! their transcript sets, not the union: reads sharing a UMI came from one
//! molecule, so a transcript missing from any of them is excluded.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use crate::align::transcript::Transcript;
use crate::quant::transcriptome::TranscriptomeIndex;

/// Size of the distance histogram (STAR `transcriptDistCount.resize(10000)`).
pub const DIST_COUNT_LEN: usize = 10_000;

/// One read's contribution: the cell, the UMI, and every concordant
/// `(transcript, distance to its 3' end)`.
pub type Record = (u32, u64, Vec<(u32, u32)>);

/// The record side of the feature, accumulated across reads.
#[derive(Debug, Default)]
pub struct Transcript3pAcc {
    /// Histogram of observed 3'-end distances, capped at [`DIST_COUNT_LEN`].
    pub dist_count: Vec<u32>,
    pub records: Vec<Record>,
}

impl Transcript3pAcc {
    pub fn new() -> Self {
        Self {
            dist_count: vec![0; DIST_COUNT_LEN],
            records: Vec::new(),
        }
    }

    /// Record one read's concordant transcripts.
    pub fn add(&mut self, cb: u32, umi: u64, hits: Vec<(u32, u32)>) {
        for &(_, dist) in &hits {
            if let Some(slot) = self.dist_count.get_mut(dist as usize) {
                *slot += 1;
            }
        }
        self.records.push((cb, umi, hits));
    }
}

/// Every transcript this alignment is concordant with, and the spliced distance
/// from the read's 3'-most base to that transcript's 3' end.
///
/// Concordance is exactly what the transcriptome projection already enforces:
/// the alignment lies inside the transcript, is purely exonic, and every splice
/// junction it crosses is one of the transcript's. A projection that survives
/// is concordant; one that does not, is not.
///
/// The projection puts the transcript's 5' end at coordinate zero for both
/// strands, so the distance is the same expression either way.
pub fn concordant_transcripts(
    align: &Transcript,
    tx: &TranscriptomeIndex,
    lread: u32,
) -> Vec<(u32, u32)> {
    crate::quant::transcriptome::align_to_transcripts(align, tx, lread)
        .into_iter()
        .filter_map(|proj| {
            let tr = proj.chr_idx; // the projection stores the transcript index here
            let tr_len = u64::from(*tx.tr_length.get(tr)?);
            let dist = tr_len.checked_sub(proj.genome_end)?;
            Some((tr as u32, u32::try_from(dist).ok()?))
        })
        .collect()
}

/// Parse `--soloClusterCBfile`: whitespace-separated `CB cluster` pairs.
///
/// A barcode not in the whitelist is skipped rather than rejected, as STAR
/// does — the file is usually produced by an external clustering run against a
/// filtered matrix, so it can legitimately name barcodes this run did not keep.
/// A trailing barcode with no cluster ends the parse, matching STAR's stream
/// extraction failing.
pub fn load_cluster_cb(
    text: &str,
    barcode_index: impl Fn(&str) -> Option<u32>,
) -> BTreeMap<u32, u32> {
    let mut out = BTreeMap::new();
    let mut tokens = text.split_whitespace();
    while let Some(cb) = tokens.next() {
        let Some(cluster) = tokens.next().and_then(|s| s.parse::<u32>().ok()) else {
            break;
        };
        if let Some(i) = barcode_index(cb) {
            out.insert(i, cluster);
        }
    }
    out
}

/// The 3'-distance distribution, estimated from the observed histogram.
///
/// Returns the normalised distribution (written out for inspection), its
/// natural log (the per-read weight table), and a per-transcript factor that
/// corrects for transcripts shorter than the distribution's support — a 300-base
/// transcript cannot produce a read 2000 bases from its end, so its abundance
/// must be scaled by the mass it can actually reach.
///
/// The histogram is smoothed with a running average and cut at the first
/// minimum past 1000, which is where the 3' peak has decayed into the body.
/// STAR's running-average divisor is `min(2N+1, ii + N)` rather than the number
/// of elements actually summed, so the first window is scaled slightly wrong;
/// that is reproduced, because the cut point and the weights both depend on it.
fn dist_function(dist_count: &[u32], tr_length: &[u32]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    const RUN_AVER_N: i64 = 50;
    let len = dist_count.len() as i64;
    let mut dist_fun = vec![0.0f64; dist_count.len()];
    let mut i = 0i64;
    while i < len - RUN_AVER_N - 1 {
        let lo = (i - RUN_AVER_N).max(0) as usize;
        let hi = (i + RUN_AVER_N + 1) as usize;
        let sum: u64 = dist_count[lo..hi].iter().map(|&x| u64::from(x)).sum();
        let divisor = (2 * RUN_AVER_N + 1).min(i + RUN_AVER_N);
        dist_fun[i as usize] = sum as f64 / divisor as f64;
        i += 1;
    }

    // Walk up to the peak past 1000, then down to the following minimum.
    let mut imax = 1000usize;
    while imax + 1 < dist_fun.len() && dist_fun[imax + 1] > dist_fun[imax] {
        imax += 1;
    }
    while imax + 1 < dist_fun.len() && dist_fun[imax + 1] < dist_fun[imax] {
        imax += 1;
    }
    dist_fun.truncate(imax);

    let norm: f64 = dist_fun.iter().sum();
    if norm > 0.0 {
        for f in &mut dist_fun {
            *f /= norm;
        }
    }
    let normalised = dist_fun.clone();

    let mut cumulative = Vec::with_capacity(dist_fun.len());
    let mut acc = 0.0f64;
    for &d in &dist_fun {
        acc += d;
        cumulative.push(acc);
    }
    let tr_factor: Vec<f64> = tr_length
        .iter()
        .map(|&l| {
            let l = l as usize;
            if l >= 1 && l < cumulative.len() && cumulative[l - 1] > 0.0 {
                -(cumulative[l - 1].ln())
            } else {
                0.0
            }
        })
        .collect();

    let log_dist: Vec<f64> = dist_fun.iter().map(|&x| x.ln()).collect();
    (normalised, log_dist, tr_factor)
}

/// The per-cluster EM over UMIs.
///
/// A UMI compatible with one transcript is evidence for it outright; a UMI
/// compatible with several is split between them in proportion to the current
/// abundance estimate, and the estimate is re-derived, until it stops moving.
/// Transcripts that fall below `1e-8` of the total, or whose estimate stops
/// changing, are frozen — otherwise the loop spends its iterations on
/// transcripts that have already decided.
fn cluster_em(umis: &BTreeMap<u64, Vec<(u32, f64)>>, n_tr: usize, tr_factor: &[f64]) -> Vec<f64> {
    let mut unique = vec![0.0f64; n_tr];
    let mut initial = vec![0.0f64; n_tr];
    let mut multi: Vec<Vec<(u32, f64)>> = Vec::new();
    let mut n_umi: u64 = 0;

    for hits in umis.values() {
        match hits.len() {
            // An empty intersection means the reads sharing this UMI agreed on
            // no transcript at all, so it is evidence for nothing.
            0 => {}
            1 => {
                unique[hits[0].0 as usize] += 1.0;
                initial[hits[0].0 as usize] += 1.0;
                n_umi += 1;
            }
            n => {
                // Shift by the maximum before exponentiating: these are log
                // weights and the raw values underflow.
                let max = hits
                    .iter()
                    .map(|&(_, w)| w)
                    .fold(f64::NEG_INFINITY, f64::max);
                let share = 1.0 / n as f64;
                let mut v = Vec::with_capacity(n);
                for &(tr, w) in hits {
                    initial[tr as usize] += share;
                    v.push((tr, (w - max).exp()));
                }
                multi.push(v);
                n_umi += 1;
            }
        }
    }

    let mut old = initial;
    let mut new = vec![0.0f64; n_tr];
    let mut converged = vec![false; n_tr];
    const DIFF_MAX: f64 = 1e-5;
    let diff_one = DIFF_MAX * 0.1;
    let expr_threshold = 1e-8 * n_umi as f64;

    for _ in 0..10_000 {
        new.copy_from_slice(&unique);
        for v in &multi {
            let denom: f64 = v.iter().map(|&(tr, w)| w * old[tr as usize]).sum();
            if denom == 0.0 {
                continue;
            }
            for &(tr, w) in v {
                if !converged[tr as usize] {
                    new[tr as usize] += w * old[tr as usize] / denom;
                }
            }
        }
        let mut worst = 0.0f64;
        for itr in 0..n_tr {
            if converged[itr] || old[itr] == 0.0 {
                continue;
            }
            let diff = (new[itr] - old[itr]).abs() / old[itr];
            worst = worst.max(diff);
            if new[itr] < expr_threshold {
                converged[itr] = true;
                unique[itr] = 0.0;
            }
            if diff < diff_one {
                converged[itr] = true;
                unique[itr] = new[itr];
            }
        }
        if worst < DIFF_MAX {
            break;
        }
        std::mem::swap(&mut new, &mut old);
    }

    // Undo the length correction and put the total back on the UMI scale, so
    // the numbers are comparable across clusters of different depth.
    let mut out = new;
    let mut norm = 0.0f64;
    for (itr, v) in out.iter_mut().enumerate() {
        *v *= tr_factor[itr].exp();
        norm += *v;
    }
    if norm > 0.0 {
        let scale = n_umi as f64 / norm;
        for v in &mut out {
            *v *= scale;
        }
    }
    out
}

/// The three output files: the cluster × transcript matrix, the transcript
/// list, and the estimated distance distribution.
pub struct Transcript3pOutput {
    pub matrix: String,
    pub features: String,
    pub distance_distribution: String,
}

/// Quantify the accumulated records.
pub fn quantify(
    acc: &Transcript3pAcc,
    tx: &TranscriptomeIndex,
    cluster_cb: &BTreeMap<u32, u32>,
) -> Transcript3pOutput {
    let n_tr = tx.n_transcripts();
    let (normalised, log_dist, tr_factor) = dist_function(&acc.dist_count, &tx.tr_length);
    let imax = log_dist.len();

    let mut distance_distribution = String::new();
    for &v in &normalised {
        distance_distribution.push_str(&fmt_cpp_g6(v));
        distance_distribution.push('\n');
    }

    // Per cluster, per UMI, the transcripts still compatible with every read
    // carrying that UMI.
    let mut per_cluster: BTreeMap<u32, BTreeMap<u64, Vec<(u32, f64)>>> = BTreeMap::new();
    for (cb, umi, hits) in &acc.records {
        let Some(&cluster) = cluster_cb.get(cb) else {
            continue; // this cell is not in any cluster
        };
        let mut weighted: Vec<(u32, f64)> = hits
            .iter()
            .filter(|&&(_, dist)| (dist as usize) < imax)
            .map(|&(tr, dist)| (tr, log_dist[dist as usize] + tr_factor[tr as usize]))
            .collect();
        if weighted.is_empty() {
            continue;
        }
        weighted.sort_by_key(|&(tr, _)| tr);

        let umi_map = per_cluster.entry(cluster).or_default();
        match umi_map.get(umi) {
            None => {
                umi_map.insert(*umi, weighted);
            }
            Some(existing) => {
                // Intersect: one molecule, so a transcript absent from either
                // read cannot be its source. Weights add, since the reads are
                // independent observations of the same molecule.
                let mut merged = Vec::new();
                let mut j = 0usize;
                for &(tr, w) in existing {
                    while j < weighted.len() && weighted[j].0 < tr {
                        j += 1;
                    }
                    if j == weighted.len() {
                        break;
                    }
                    if weighted[j].0 == tr {
                        merged.push((tr, w + weighted[j].1));
                    }
                }
                umi_map.insert(*umi, merged);
            }
        }
    }

    let expression: BTreeMap<u32, Vec<f64>> = per_cluster
        .iter()
        .map(|(&cl, umis)| (cl, cluster_em(umis, n_tr, &tr_factor)))
        .collect();

    let clusters: BTreeSet<u32> = cluster_cb.values().copied().collect();
    let n_clusters = clusters.iter().max().copied().unwrap_or(0);
    let nnz: usize = expression
        .values()
        .map(|e| e.iter().filter(|&&v| v > 0.0).count())
        .sum();

    let mut matrix = String::from("%%MatrixMarket matrix coordinate real general\n%\n");
    let _ = writeln!(matrix, "{n_tr} {n_clusters} {nnz}");
    for (&cl, expr) in &expression {
        for (itr, &v) in expr.iter().enumerate() {
            if v > 0.0 {
                let _ = writeln!(matrix, "{} {} {}", itr + 1, cl, fmt_cpp_g6(v));
            }
        }
    }

    let mut features = String::new();
    for (i, id) in tx.tr_ids.iter().enumerate() {
        let gene = tx.tr_gene_idx[i] as usize;
        let name = tx.gene_names.get(gene).map_or("-", String::as_str);
        let _ = writeln!(features, "{id}\t{name}\tTranscript3p");
    }

    Transcript3pOutput {
        matrix,
        features,
        distance_distribution,
    }
}

/// Format like C++'s default `ostream << double`: six significant digits,
/// fixed notation for exponents in `[-4, 6)` and scientific outside it, with
/// trailing zeros trimmed.
///
/// Worth the trouble because the normalised distribution runs down to ~1e-4,
/// where Rust's `{}` and C++'s default disagree on both notation and digits.
pub fn fmt_cpp_g6(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let e = v.abs().log10().floor() as i32;
    if !(-4..6).contains(&e) {
        let s = format!("{v:.5e}");
        let (mantissa, exponent) = s.split_once('e').unwrap();
        let mut mantissa = mantissa.to_string();
        if mantissa.contains('.') {
            while mantissa.ends_with('0') {
                mantissa.pop();
            }
            if mantissa.ends_with('.') {
                mantissa.pop();
            }
        }
        let exp: i32 = exponent.parse().unwrap();
        let sign = if exp < 0 { '-' } else { '+' };
        return format!("{mantissa}e{sign}{:02}", exp.abs());
    }
    // Exact comparison on purpose: C++ prints an integral value without a
    // decimal point, and "integral" there means exactly integral.
    #[allow(clippy::float_cmp)]
    let is_integral = v == v.trunc();
    if is_integral {
        return format!("{}", v as i64);
    }
    let decimals = (5 - e).max(0) as usize;
    let mut s = format!("{v:.decimals$}");
    while s.ends_with('0') {
        s.pop();
    }
    if s.ends_with('.') {
        s.pop();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A histogram with a clear 3' peak: the estimated distribution is a
    /// probability distribution, and it is cut where the peak ends rather than
    /// running to the end of the buffer.
    #[test]
    fn the_distance_distribution_is_normalised_and_cut_at_the_peak() {
        let mut counts = vec![0u32; DIST_COUNT_LEN];
        // A broad peak around 1200, decaying to nothing by 3000.
        for (i, c) in counts.iter_mut().enumerate().take(3000) {
            let d = (i as f64 - 1200.0).abs();
            *c = (1000.0 * (-d / 400.0).exp()) as u32;
        }
        let tr_len = vec![5000u32; 4];
        let (normalised, log_dist, _) = dist_function(&counts, &tr_len);

        assert_ne!(normalised, []);
        assert!(
            normalised.len() < DIST_COUNT_LEN,
            "must be cut, not full length"
        );
        let total: f64 = normalised.iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "should sum to 1, got {total}");
        assert_eq!(log_dist.len(), normalised.len());
    }

    /// The per-transcript factor exists to stop short transcripts being
    /// under-counted: one shorter than the distribution's reach gets a positive
    /// correction, one longer than it gets none.
    #[test]
    fn short_transcripts_get_a_length_correction() {
        let mut counts = vec![0u32; DIST_COUNT_LEN];
        for (i, c) in counts.iter_mut().enumerate().take(3000) {
            let d = (i as f64 - 1200.0).abs();
            *c = (1000.0 * (-d / 400.0).exp()) as u32;
        }
        let (normalised, _, tr_factor) = dist_function(&counts, &[500u32, 100_000u32]);
        assert!(
            tr_factor[0] > 0.0,
            "a transcript shorter than the distribution needs correcting"
        );
        assert_eq!(
            tr_factor[1].to_bits(),
            0.0f64.to_bits(),
            "one longer than the distribution's support needs none"
        );
        assert!(normalised.len() < 100_000);
    }

    /// A UMI compatible with one transcript is evidence for that transcript and
    /// nothing else.
    #[test]
    fn a_unique_umi_goes_entirely_to_its_transcript() {
        let mut umis = BTreeMap::new();
        umis.insert(1u64, vec![(0u32, 0.0f64)]);
        umis.insert(2u64, vec![(0u32, 0.0f64)]);
        let expr = cluster_em(&umis, 3, &[0.0, 0.0, 0.0]);
        assert!(expr[0] > 0.0);
        assert_eq!(expr[1].to_bits(), 0.0f64.to_bits());
        assert_eq!(expr[2].to_bits(), 0.0f64.to_bits());
    }

    /// An ambiguous UMI is resolved by the unambiguous ones around it: with
    /// nine UMIs pointing at transcript 0 and one split between 0 and 1, the EM
    /// gives almost all of the split one to 0.
    #[test]
    fn an_ambiguous_umi_follows_the_evidence() {
        let mut umis = BTreeMap::new();
        for u in 0..9u64 {
            umis.insert(u, vec![(0u32, 0.0f64)]);
        }
        umis.insert(9, vec![(0u32, 0.0f64), (1u32, 0.0f64)]);
        let expr = cluster_em(&umis, 2, &[0.0, 0.0]);
        assert!(
            expr[0] > expr[1] * 5.0,
            "the well-supported transcript should take the ambiguous UMI: {expr:?}"
        );
    }

    /// Reads sharing a UMI came from one molecule, so the transcript sets
    /// intersect rather than accumulate. A transcript missing from the second
    /// read is dropped, even though the first read supported it.
    #[test]
    fn a_umi_seen_twice_keeps_only_the_shared_transcripts() {
        let mut acc = Transcript3pAcc::new();
        acc.add(0, 42, vec![(0, 100), (1, 200)]);
        acc.add(0, 42, vec![(1, 150), (2, 300)]);
        let tx = tiny_index();
        let mut clusters = BTreeMap::new();
        clusters.insert(0u32, 1u32);
        let out = quantify(&acc, &tx, &clusters);
        // Transcript 1 (row 2) is the only one in both reads.
        let rows: Vec<&str> = out.matrix.lines().skip(3).collect();
        assert!(
            rows.iter().all(|r| r.starts_with("2 ")),
            "only the shared transcript should be quantified: {rows:?}"
        );
    }

    /// A cell that no cluster claims contributes nothing, rather than being
    /// silently folded into cluster 0.
    #[test]
    fn a_cell_outside_every_cluster_is_skipped() {
        let mut acc = Transcript3pAcc::new();
        acc.add(7, 1, vec![(0, 100)]);
        let tx = tiny_index();
        let out = quantify(&acc, &tx, &BTreeMap::new());
        assert_eq!(out.matrix.lines().count(), 3, "header only, no entries");
    }

    #[test]
    fn cluster_file_parsing_skips_unknown_barcodes() {
        let known = ["AAAA", "CCCC"];
        let index = |cb: &str| known.iter().position(|&k| k == cb).map(|i| i as u32);
        let map = load_cluster_cb("AAAA 1\nGGGG 2\nCCCC 3\n", index);
        assert_eq!(map.get(&0), Some(&1));
        assert_eq!(map.get(&1), Some(&3));
        assert_eq!(map.len(), 2, "the unknown barcode is skipped, not an error");
    }

    #[test]
    fn a_trailing_barcode_without_a_cluster_ends_the_parse() {
        let known = ["AAAA", "CCCC"];
        let index = |cb: &str| known.iter().position(|&k| k == cb).map(|i| i as u32);
        let map = load_cluster_cb("AAAA 1\nCCCC\n", index);
        assert_eq!(map.len(), 1);
    }

    /// The number formatting has to match C++'s default stream output, which
    /// neither Rust's `{}` nor `{:e}` does on its own.
    #[test]
    fn numbers_are_formatted_the_way_c_plus_plus_prints_them() {
        assert_eq!(fmt_cpp_g6(0.0), "0");
        assert_eq!(fmt_cpp_g6(1.0), "1");
        assert_eq!(fmt_cpp_g6(0.5), "0.5");
        assert_eq!(fmt_cpp_g6(0.000_123_456_789), "0.000123457");
        assert_eq!(fmt_cpp_g6(1.234_567_89e-7), "1.23457e-07");
        assert_eq!(fmt_cpp_g6(1.5e7), "1.5e+07");
    }

    /// Three transcripts of one gene, enough for the quantifier to run.
    fn tiny_index() -> TranscriptomeIndex {
        TranscriptomeIndex {
            tr_ids: vec!["t0".into(), "t1".into(), "t2".into()],
            tr_chr_idx: vec![0; 3],
            tr_strand: vec![1; 3],
            tr_gene_idx: vec![0; 3],
            gene_ids: vec!["g0".into()],
            gene_names: vec!["G0".into()],
            gene_biotypes: vec!["protein_coding".into()],
            tr_start: vec![0; 3],
            tr_end: vec![1000; 3],
            tr_exons: vec![Vec::new(); 3],
            tr_length: vec![1000; 3],
            tr_exi: vec![0; 3],
            tr_order: vec![0, 1, 2],
            tr_starts_sorted: vec![0; 3],
            tr_end_max_sorted: vec![1000; 3],
            retained_introns: None,
        }
    }
}
