//! STAR-faithful suffix array construction via the `caps-sa` crate.
//!
//! The `caps-sa` crate sorts a generic text; this module is the bridge
//! between STAR's expectations and that generic algorithm:
//!
//! 1. **Sentinel transform.** STAR's SA is a *generalized* suffix array over
//!    `T = forward || RC` with each inter-chromosome spacer run acting as a
//!    sentinel and a position tie-break. Distinct per-run sentinels (plus a
//!    final terminal sentinel) turn STAR's order into a plain lexicographic
//!    SA — see the project plan for the correctness argument.
//! 2. **Standard SA construction** via `caps-sa`. The
//!    [`caps_sa::build_ext_mem`] path is the default for production-
//!    scale genomes (transformed text ≥ 16 MB): it streams the SA from
//!    disk-spilling buckets and keeps peak RAM bounded at ~`O(text +
//!    n/p)`. For smaller inputs — including the synthetic test fixtures
//!    where bin-padding dominates the transformed text and would make
//!    ext-mem pathologically slow — the in-memory
//!    [`caps_sa::build_in_memory`] is used instead. Override with
//!    `RUSTAR_USE_IN_MEM=1` (force in-mem) or `RUSTAR_USE_EXT_MEM={1,0}`
//!    (force ext-mem on/off); see [`use_ext_mem`] for the full
//!    decision table.
//! 3. **Filter + pack.** Each output position from caps-sa is fed into a
//!    streaming callback that keeps only ACGT (`≤ 3`) starts, applies
//!    STAR's `genomeSAsparseD` stride, and packs into the `PackedArray`
//!    with STAR's strand bit encoding. The full SA is never materialized
//!    twice — for the ext-mem path it isn't materialized at all.
//!
//! Phase 1 uses the byte-alphabet path (sentinel value `5 + n_seg ≤ 255`,
//! i.e. ≤125 chromosomes — covers yeast and the human primary assembly).
//! Highly fragmented assemblies need a `u16` text and are rejected with a
//! clear error until that fallback lands in Phase 3.

use crate::error::Error;
use crate::genome::Genome;
use crate::index::packed_array::PackedArray;
use crate::index::suffix_array::SuffixArray;

/// STAR's spacer byte. Matches `GENOME_spacingChar` in
/// `STAR/source/IncludeDefine.h` and the value `5` used throughout the
/// existing rustar-aligner code.
const SPACER: u8 = 5;

/// Base value for per-run sentinels: run 0 → 5, run 1 → 6, ….
const SENTINEL_BASE: u8 = 5;

/// Build the suffix array for `genome` using the caps-sa sample-sort
/// construction.
///
/// `genome.sequence` must already be of length `2 * genome.n_genome` (forward
/// + reverse complement laid out as `[forward | RC]`). The current call site
///   (`GenomeIndex::build` after `genome.append_sjdb`) already satisfies this.
pub fn build(genome: &Genome) -> Result<SuffixArray, Error> {
    let n_genome = genome.n_genome as usize;
    let n2 = 2 * n_genome;
    if genome.sequence.len() < n2 {
        return Err(Error::Index(format!(
            "sa_build: genome.sequence length {} < 2 * n_genome ({})",
            genome.sequence.len(),
            n2
        )));
    }

    let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
    let gstrand_mask = (1u64 << gstrand_bit) - 1;
    let word_length = gstrand_bit + 1;
    let n2_bit = 1u64 << gstrand_bit;

    // (1) Sentinel transform.
    let (t_prime, n_seg) = build_sentinel_transformed_text(&genome.sequence[..n2])?;
    log::info!(
        "sa_build: transformed text length = {}, {} per-segment sentinels \
         (alphabet max = {})",
        t_prime.len(),
        n_seg,
        SENTINEL_BASE as u32 + n_seg
    );

    // (2) Filter+pack budget: STAR's iteration in `Genome_genomeGenerate.cpp`
    //     is over the REVERSED buffer with step `gSAsparseD`, which
    //     corresponds to original-T positions `p` with
    //     `(2*n_genome - 1 - p) % gSAsparseD == 0`. With the default
    //     `gSAsparseD = 1` (rustar-aligner currently uses the default) the
    //     stride collapses to "every ACGT position." A non-1 stride is a
    //     future addition once `params.genome_sa_sparse_d` is threaded
    //     through.
    let _ = n_seg;
    let sparse_d: u64 = 1;
    let n_sa_kept = count_kept_positions(&genome.sequence[..n2], sparse_d);
    log::info!("sa_build: {n_sa_kept} entries after ACGT + sparse-d={sparse_d} filter");

    let mut data = PackedArray::new(word_length, n_sa_kept);
    let n_genome_u64 = n_genome as u64;
    let n2_minus_one = n2 as u64 - 1;
    let mut out_idx: usize = 0;

    // The packer is shared between the in-memory and ext-mem paths. It
    // consumes SA entries one at a time in lex order and writes the
    // strand-bit-encoded packed values into `data`.
    let mut pack_one = |sa_pos: u64| {
        let p = sa_pos as usize;
        if p == n2 {
            return; // terminal sentinel
        }
        if genome.sequence[p] >= 4 {
            return; // N or spacer — STAR's `G[ii] < 4` filter
        }
        if sparse_d != 1 && !(n2_minus_one - sa_pos).is_multiple_of(sparse_d) {
            return;
        }
        let packed_value = if p < n_genome {
            sa_pos
        } else {
            (sa_pos - n_genome_u64) | n2_bit
        };
        data.write(out_idx, packed_value);
        out_idx += 1;
    };

    // (3) Standard SA over the transformed text via caps-sa.
    if use_ext_mem(t_prime.len()) {
        log::info!(
            "sa_build: invoking caps-sa::build_ext_mem (text len {})",
            t_prime.len()
        );
        let opts = caps_sa::ExtMemOpts {
            work_dir: std::env::temp_dir(),
            ..Default::default()
        };
        caps_sa::build_ext_mem(&t_prime, &opts, |sa_pos| {
            pack_one(sa_pos);
            Ok(())
        })
        .map_err(|e| Error::Index(format!("caps-sa::build_ext_mem failed: {e}")))?;
        drop(t_prime);
    } else {
        log::info!(
            "sa_build: invoking caps-sa::build_in_memory (text len {})",
            t_prime.len()
        );
        let sa: Vec<u64> = caps_sa::build_in_memory(&t_prime);
        // Free the transformed text early — only `data` and `genome` are
        // read from this point on.
        drop(t_prime);
        for &sa_pos in &sa {
            pack_one(sa_pos);
        }
    }

    debug_assert_eq!(
        out_idx, n_sa_kept,
        "sa_build: packed {out_idx} entries but counted {n_sa_kept}",
    );

    Ok(SuffixArray {
        data,
        gstrand_bit,
        gstrand_mask,
    })
}

/// Select the in-memory vs. external-memory caps-sa path.
///
/// The external-memory path is the **default for genomes at or above
/// `EXT_MEM_THRESHOLD_BYTES`**: it streams the SA from disk-spilling
/// buckets and keeps peak RAM bounded at ~`O(text + n/p)` regardless of
/// genome size. Below the threshold, the in-memory path is preferred:
/// it's faster on tiny inputs, and ext-mem becomes pathological when
/// the transformed text is dominated by `genomeChrBinNbits` padding
/// (e.g. a 20 kb test fixture rounds to 256 kb of padded text, of which
/// >90% is constant spacer bytes — `caps-sa` sorts those positions and
/// then we filter them out, but the sort itself does `O(spacer_run_len²)`
/// work because every spacer-starting suffix shares a near-maximal LCP
/// with every other).
///
/// Explicit overrides (decision order):
///
/// 1. `RUSTAR_USE_IN_MEM=1` (also accepts `true`/`yes`/`on`) forces
///    in-memory.
/// 2. `RUSTAR_USE_EXT_MEM=1` / `=0` (also accepts the usual aliases)
///    forces ext-mem on / off respectively. Retained for backward
///    compatibility with the earlier opt-in flag.
///
/// Without overrides the threshold below decides. 16 MB transformed
/// text places yeast (~30 MB), human primary chromosomes (~6 MB-200 MB
/// transformed text per chromosome), and any input where the padded
/// genome would dwarf in-memory's `~4 × n × sizeof(I)` working set on
/// the ext-mem path; tiny synthetic test fixtures stay in-memory.
fn use_ext_mem(text_len: usize) -> bool {
    if let Ok(v) = std::env::var("RUSTAR_USE_IN_MEM") {
        if matches!(v.as_str(), "1" | "true" | "yes" | "on") {
            return false;
        }
    }
    if let Ok(v) = std::env::var("RUSTAR_USE_EXT_MEM") {
        if matches!(v.as_str(), "0" | "false" | "no" | "off") {
            return false;
        }
        if matches!(v.as_str(), "1" | "true" | "yes" | "on") {
            return true;
        }
    }
    const EXT_MEM_THRESHOLD_BYTES: usize = 16 * 1024 * 1024;
    text_len >= EXT_MEM_THRESHOLD_BYTES
}

/// Count positions kept by STAR's `G[ii] < 4` + `--genomeSAsparseD` filter.
fn count_kept_positions(text: &[u8], sparse_d: u64) -> usize {
    let n2 = text.len();
    if sparse_d == 1 {
        text.iter().filter(|&&b| b < 4).count()
    } else {
        let n2_minus_one = n2 as u64 - 1;
        (0..n2 as u64)
            .filter(|&p| text[p as usize] < 4 && (n2_minus_one - p).is_multiple_of(sparse_d))
            .count()
    }
}

/// Build the per-segment sentinel-transformed text `T'` of length `n + 1`.
///
/// Each maximal run of spacer bytes (value `5`) is numbered in position
/// order; the run at index `i` has every position stamped with the sentinel
/// value `SENTINEL_BASE + i`. After the last byte we append one extra
/// terminal sentinel = `SENTINEL_BASE + n_seg` (a value larger than every
/// per-run sentinel), so the final RC suffix has a sentinel-terminator and
/// the `caps-sa` "implicit smallest sentinel at the end" never affects any
/// kept suffix's order.
fn build_sentinel_transformed_text(genome: &[u8]) -> Result<(Vec<u8>, u32), Error> {
    let n = genome.len();

    // First pass: count maximal spacer runs.
    let mut n_seg: u32 = 0;
    let mut in_run = false;
    for &b in genome {
        if b == SPACER {
            if !in_run {
                n_seg += 1;
                in_run = true;
            }
        } else {
            in_run = false;
        }
    }

    let terminal = SENTINEL_BASE as u32 + n_seg;
    if terminal > 255 {
        return Err(Error::Index(format!(
            "caps-sa sentinel transform: {n_seg} spacer runs require alphabet \
             value up to {terminal}; this exceeds the byte alphabet (max 255). \
             The u16 fallback path is not yet implemented — currently supports \
             up to ~125 chromosomes (yeast + human primary assembly)."
        )));
    }

    // Second pass: emit the transform.
    let mut out = Vec::with_capacity(n + 1);
    let mut run_idx: u32 = 0;
    let mut in_run = false;
    for &b in genome {
        if b == SPACER {
            if !in_run {
                in_run = true;
            }
            out.push((SENTINEL_BASE as u32 + run_idx) as u8);
        } else {
            if in_run {
                in_run = false;
                run_idx += 1;
            }
            out.push(b);
        }
    }
    if in_run {
        run_idx += 1;
    }
    debug_assert_eq!(run_idx, n_seg);

    // Terminal sentinel — distinct from every per-run value, larger than
    // every real base, and the unique-maximum symbol of `T'`.
    out.push(terminal as u8);

    Ok((out, n_seg))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sentinel_transform_counts_runs() {
        // [a, 5, 5, c, 5, 5, 5, t]
        let input: Vec<u8> = vec![0, 5, 5, 1, 5, 5, 5, 3];
        let (out, n_seg) = build_sentinel_transformed_text(&input).unwrap();
        assert_eq!(n_seg, 2);
        // Run 0 → sentinel 5; run 1 → sentinel 6; terminal sentinel = 7.
        assert_eq!(out, vec![0, 5, 5, 1, 6, 6, 6, 3, 7]);
    }

    #[test]
    fn sentinel_transform_trailing_run() {
        // Input ending in a spacer run.
        let input: Vec<u8> = vec![0, 1, 5, 5];
        let (out, n_seg) = build_sentinel_transformed_text(&input).unwrap();
        assert_eq!(n_seg, 1);
        // Run 0 → sentinel 5; terminal sentinel = 6.
        assert_eq!(out, vec![0, 1, 5, 5, 6]);
    }

    #[test]
    fn sentinel_transform_no_spacers() {
        let input: Vec<u8> = vec![0, 1, 2, 3];
        let (out, n_seg) = build_sentinel_transformed_text(&input).unwrap();
        assert_eq!(n_seg, 0);
        assert_eq!(out, vec![0, 1, 2, 3, 5]);
    }

    #[test]
    fn sentinel_transform_rejects_byte_overflow() {
        // 251 spacer runs separated by single bases → n_seg = 251,
        // terminal = 256 → overflow.
        let mut input: Vec<u8> = Vec::new();
        for _ in 0..251 {
            input.push(0); // real base
            input.push(5); // spacer
        }
        let err = build_sentinel_transformed_text(&input).unwrap_err();
        match err {
            Error::Index(msg) => assert!(msg.contains("exceeds the byte alphabet")),
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
