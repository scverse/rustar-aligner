//! Is the junction scan stalled on memory, and is it the acceptor stream?
//!
//! The stack sampler gives self time, not cache misses, so this asks the
//! question by construction instead: run the same scan over the same genome,
//! changing only how far the acceptor sits from the donor. Everything else
//! (iteration count, branch pattern, instruction mix) is held fixed, because
//! `del` enters the loop only as an address offset once it is inside the
//! intron-length range.
//!
//! If time per iteration is flat across `del`, the loop is not memory-bound on
//! the acceptor and prefetching it would be wasted work. If it climbs with
//! `del`, that stream is the cost.
//!
//! Run: cargo run --release --example jrbench

use rustar_aligner::align::score::AlignmentScorer;
use rustar_aligner::genome::{Genome, GenomeSeq};
use std::time::Instant;

/// Deterministic pseudo-random bases, so the run is reproducible and the motif
/// hit rate is the same at every `del`.
fn synthetic_genome(n: usize) -> Genome {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut seq = Vec::with_capacity(n);
    for _ in 0..n {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        seq.push((state >> 33) as u8 & 3);
    }
    Genome {
        transform_blocks: None,
        sequence: GenomeSeq::Owned(seq),
        n_genome: n as u64,
        n_genome_real: n as u64,
        n_chr_real: 1,
        chr_name: vec!["chr1".to_string()],
        chr_length: vec![n as u64],
        chr_start: vec![0, n as u64],
    }
}

fn main() {
    // Large enough that the donor and acceptor streams cannot both stay in
    // cache when they are far apart.
    const N: usize = 256 << 20; // 256 MB of bases, one byte each
    const READ_LEN: usize = 150;
    const SCANS: usize = 20_000;

    let genome = synthetic_genome(N);
    let mut scorer = AlignmentScorer::from_params_minimal();
    scorer.align_intron_min = 20;
    scorer.align_intron_max = u32::MAX;

    let mut read = vec![0u8; READ_LEN];
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    for b in read.iter_mut() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        *b = (state >> 33) as u8 & 3;
    }

    println!("{:>12}  {:>10}  {:>12}", "del", "wall", "ns/iter");
    for del in [64u64, 1_024, 16_384, 262_144, 4_194_304, 67_108_864] {
        // Spread the scans over the genome so no single window stays resident.
        let stride = (N as u64 - del - 4096) / SCANS as u64;
        let iters_per_scan = 100usize; // r_gap 0 + next_seed_len 100
        let t0 = Instant::now();
        let mut sink = 0i64;
        for k in 0..SCANS {
            let g_a_end = 2048 + k as u64 * stride;
            let (jr, _motif, score, _l, _r) = scorer.find_best_junction_position(
                &read,
                75,
                g_a_end,
                0,
                del as i64,
                &genome,
                false,
                N as u64,
                75,
                iters_per_scan,
            );
            sink += jr as i64 + score as i64;
        }
        let dt = t0.elapsed();
        let iters = (SCANS * iters_per_scan) as f64;
        println!(
            "{:>12}  {:>9.3}s  {:>11.2}   (sink {})",
            del,
            dt.as_secs_f64(),
            dt.as_secs_f64() * 1e9 / iters,
            sink
        );
    }
}
