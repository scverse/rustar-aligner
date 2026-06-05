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
//! The SA construction path:
//!
//! **Default — segmented arm** (caps-sa 0.4.1+). Spacer-free `u8`
//! text + `caps_sa::SegmentedText` + [`StarSegmentedText`] wrapper
//! that flips the boundary-order to STAR's `spacer-as-largest`
//! convention. Byte-for-byte STAR-compatible (verified by the
//! `segmented_arm_matches_sentinel_arm_byte_for_byte_*` differential
//! tests). Lower peak memory than the sentinel arms — the text
//! shrinks by the spacer-padding overhead, and no per-segment
//! sentinel bytes need allocating in the text. Handles any segment
//! count regardless of alphabet size; future genome+SJ indexes with
//! tens of thousands of junctions "just work" here.
//!
//! **Fallback — sentinel-transform arms**, opt-in via
//! `RUSTAR_USE_SENTINEL_TRANSFORM=1`. When the alphabet fits, this
//! picks the narrowest type:
//!
//! - `5 + n_seg ≤ 255` → **`Vec<u8>` sentinel-transform** (covers
//!   ≤125 chromosomes).
//! - `5 + n_seg ≤ 65535` → **`Vec<u16>` sentinel-transform** (covers
//!   up to ~32 K chromosome-like segments).
//! - Anything larger falls back to the segmented arm anyway.
//!
//! Both arms produce byte-for-byte identical packed SAs; the legacy
//! sentinel-transform path is retained as a differential-test oracle
//! and a runtime safety net. Tests call [`build_impl`] directly with
//! `force_sentinel = true` / `false` to avoid racing on the env var.
//!
//! caps-sa's `Symbol` trait covers both `u8` and `u16` natively
//! (byte-view SIMD LCP at both widths) for the fallback arms, and
//! the byte-level SIMD backs the segmented arm too via
//! [`build_ext_mem_for_positions_with`][caps_sa::build_ext_mem_for_positions_with]
//! / `LimitProvider`.
//!
//! ## WIDE_ALPHABETS — the segmented arm
//!
//! Genome+SJ indexes (where every splice junction wants its own
//! distinct sentinel to avoid cross-junction LCP collisions) can push
//! `n_seg` into the tens of thousands easily. Bumping the text width
//! to encode those sentinels — `u16` (×2 memory), `[u8; 3]` (×3),
//! `u32` (×4) — directly multiplies the genome's resident memory
//! and the LCP wall (the LCP SIMD does 64 bytes per iter regardless
//! of symbol width, so wider symbols mean fewer comparable symbols
//! per iter).
//!
//! The segmented arm avoids that entirely:
//!
//! - **Text stays at `u8`** (just real ACGTN bases — no spacers, no
//!   sentinels). `build_spacer_free` strips spacers from the original
//!   text and produces a sorted list of [`SpacerFreeSegment`]s
//!   recording the orig↔spacer-free coordinate mapping.
//! - **Segment boundaries** are passed to caps-sa as a
//!   `SegmentedText { ends: Vec<u64> }` (cumulative end positions).
//!   For 50 K segments that's 400 KB total — vs ~6 extra GB a `u16`
//!   text or the 18 GB a `[u8; 3]` text would need on the human
//!   genome.
//! - **LCP scans stop at segment boundaries.** caps-sa's `merge` /
//!   `cascade_merge` compute `lim_a = lp.lim_at(p_a)` (a binary
//!   search in the boundary list, amortisable to once per merge
//!   output by caching in the merge state) and cap the SIMD LCP
//!   call's `max_ctx` accordingly. The inner LCP loop is unchanged
//!   from the non-segmented path.
//! - **`pack_one`** wraps the caller's encoder to translate each
//!   spacer-free SA position back to an original coordinate via
//!   `sf_to_orig` — a binary search in the segments list.
//! - **Tie-break matches STAR.** caps-sa's default `SegmentedText`
//!   uses the generalised-SA "shorter-suffix-is-smaller" convention
//!   at boundary ties; STAR uses the opposite ("spacer-as-largest",
//!   equivalently "longer-suffix-is-smaller") with an ascending-
//!   position tie-break when both `lim`s coincide. We override
//!   `LimitProvider::boundary_order` on a [`StarSegmentedText`]
//!   newtype to flip the convention, recovering byte-for-byte
//!   STAR compatibility on the segmented arm.
//!
//! ### Alternatives considered (and not taken)
//!
//! - **`[u8; 3]` (24-bit) text**: would work via caps-sa's blanket
//!   `Symbol for [T; N]` impl, but at 3× the memory of the segmented
//!   arm and ~3× the LCP wall on memory-bound workloads. Cheaper to
//!   implement (no `SegmentedText` plumbing) but a permanent runtime
//!   tax once you use it.
//! - **Bit-packed variable-width text** (e.g. 10 bits/symbol for
//!   1024 sentinels): best memory for the thousand-symbol regime
//!   (1.25 B/sym), but requires a custom sub-byte LCP path in
//!   caps-sa. Not worth it as long as the segmented arm covers the
//!   workload.

use crate::error::Error;
use crate::genome::Genome;
use crate::index::packed_array::PackedArray;
use crate::index::suffix_array::SuffixArray;
use rayon::prelude::*;
use std::path::Path;

/// STAR's spacer byte. Matches `GENOME_spacingChar` in
/// `STAR/source/IncludeDefine.h` and the value `5` used throughout the
/// existing rustar-aligner code.
const SPACER: u8 = 5;

/// Base value for per-run sentinels: run 0 → 5, run 1 → 6, ….
const SENTINEL_BASE: u8 = 5;

/// Symbol type usable for the sentinel-transformed text. Implementors
/// satisfy caps-sa's [`caps_sa::Symbol`] (so the byte-view SIMD LCP
/// path works), encode an input base `0..=4` (ACGT + N), encode a
/// sentinel value `SENTINEL_BASE + run_idx`, and report the largest
/// value they can represent so the dispatch in [`build`] can pick
/// the narrowest width that fits.
trait SaSymbol: caps_sa::Symbol {
    /// Encode an original ACGT / N base (value `0..=4`) as this symbol.
    fn from_base(b: u8) -> Self;
    /// Encode the sentinel for run `idx` (value `SENTINEL_BASE + idx`).
    /// Caller has verified `SENTINEL_BASE as u32 + idx <= Self::MAX_REPRESENTABLE`.
    fn from_sentinel(idx: u32) -> Self;
    /// Largest sentinel value this width can represent. Used by the
    /// dispatch in [`build`] to pick the narrowest `SaSymbol`.
    const MAX_REPRESENTABLE: u32;
}

impl SaSymbol for u8 {
    #[inline]
    fn from_base(b: u8) -> u8 {
        b
    }
    #[inline]
    fn from_sentinel(idx: u32) -> u8 {
        // Caller has bounds-checked against MAX_REPRESENTABLE.
        (SENTINEL_BASE as u32 + idx) as u8
    }
    const MAX_REPRESENTABLE: u32 = u8::MAX as u32;
}

impl SaSymbol for u16 {
    #[inline]
    fn from_base(b: u8) -> u16 {
        b as u16
    }
    #[inline]
    fn from_sentinel(idx: u32) -> u16 {
        // Caller has bounds-checked against MAX_REPRESENTABLE.
        (SENTINEL_BASE as u32 + idx) as u16
    }
    const MAX_REPRESENTABLE: u32 = u16::MAX as u32;
}

/// Build the suffix array for `genome` using the caps-sa sample-sort
/// construction.
///
/// `genome.sequence` must already be of length `2 * genome.n_genome` (forward
/// + reverse complement laid out as `[forward | RC]`). The current call site
///   (`GenomeIndex::build` after `genome.append_sjdb`) already satisfies this.
pub fn build(genome: &Genome) -> Result<SuffixArray, Error> {
    // Production entry: read the `RUSTAR_USE_SENTINEL_TRANSFORM` env
    // var once and delegate. Tests should call [`build_impl`] directly
    // to avoid racing on the env var with parallel tests.
    let force_sentinel = matches!(
        std::env::var("RUSTAR_USE_SENTINEL_TRANSFORM")
            .ok()
            .as_deref(),
        Some("1" | "true" | "yes" | "on")
    );
    build_impl(genome, force_sentinel)
}

/// Inner build that takes `force_sentinel` as an explicit argument
/// rather than reading the env var, so tests can exercise the
/// sentinel-transform fallback without racing on the shared
/// `RUSTAR_USE_SENTINEL_TRANSFORM` env var when cargo runs them in
/// parallel.
///
/// When `force_sentinel` is `true` and the alphabet fits (`u8` or
/// `u16`), the corresponding sentinel-transform arm is used.
/// Otherwise the segmented arm (the default since
/// caps-sa 0.4.1) is used.
pub(crate) fn build_impl(genome: &Genome, force_sentinel: bool) -> Result<SuffixArray, Error> {
    let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
    let gstrand_mask = (1u64 << gstrand_bit) - 1;
    let word_length = gstrand_bit + 1;

    // (1) Pre-count the kept positions so we can pre-size the
    //     in-memory `PackedArray`. The streaming path
    //     ([`build_streaming_impl`]) doesn't need this — it just
    //     pipes each entry into the caller's `emit` closure — but
    //     since the `PackedArray` constructor needs a fixed length,
    //     we count here.
    let n_genome = genome.n_genome as usize;
    let n2 = 2 * n_genome;
    let n_sa_kept: usize = genome.sequence[..n2.min(genome.sequence.len())]
        .par_iter()
        .filter(|&&b| b < 4)
        .count();

    let mut data = PackedArray::new(word_length, n_sa_kept);
    let mut out_idx: usize = 0;

    // Stream construction into the in-memory `PackedArray`. This is
    // exactly the path used by tests and by any non-streaming
    // caller; production [`genomeGenerate`] uses
    // [`build_streaming`] directly and bypasses this allocation.
    build_streaming_impl(genome, force_sentinel, None, |packed_value| {
        data.write(out_idx, packed_value);
        out_idx += 1;
        Ok(())
    })?;

    debug_assert_eq!(out_idx, n_sa_kept);
    Ok(SuffixArray {
        data,
        gstrand_bit,
        gstrand_mask,
    })
}

/// Public streaming entry. Calls `emit(packed_value)` for each SA
/// entry in lexicographic order, packed in STAR's strand-bit
/// encoding (forward `p → p`, reverse `p → (p - n_genome) |
/// (1 << gstrand_bit)`). Returns `(gstrand_bit, gstrand_mask,
/// n_entries)` — `n_entries` matches what
/// [`SuffixArray::data.len()`][PackedArray::len] would have been.
/// The caller is responsible for whatever sink the entries land in
/// (an on-disk SA file, an in-RAM PackedArray, both).
///
/// Reads `RUSTAR_USE_SENTINEL_TRANSFORM` once at entry; tests
/// should go through [`build_streaming_impl`] directly to avoid
/// racing on the env var.
pub fn build_streaming<F>(
    genome: &Genome,
    temp_dir: Option<&Path>,
    emit: F,
) -> Result<(u32, u64, usize), Error>
where
    F: FnMut(u64) -> Result<(), Error>,
{
    let force_sentinel = matches!(
        std::env::var("RUSTAR_USE_SENTINEL_TRANSFORM")
            .ok()
            .as_deref(),
        Some("1" | "true" | "yes" | "on")
    );
    let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
    let gstrand_mask = (1u64 << gstrand_bit) - 1;
    let n_kept = build_streaming_impl(genome, force_sentinel, temp_dir, emit)?;
    Ok((gstrand_bit, gstrand_mask, n_kept))
}

/// Shared implementation behind both [`build_impl`] (in-memory
/// `PackedArray` sink) and [`build_streaming`] (on-disk SA file
/// sink). Counts kept positions, dispatches to the right caps-sa
/// arm, and invokes `emit(packed_value)` for each SA entry. Returns
/// the number of entries emitted.
pub(crate) fn build_streaming_impl<F>(
    genome: &Genome,
    force_sentinel: bool,
    temp_dir: Option<&Path>,
    mut emit: F,
) -> Result<usize, Error>
where
    F: FnMut(u64) -> Result<(), Error>,
{
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
    let n2_bit = 1u64 << gstrand_bit;

    // (1) Count spacer runs so we can pick the narrowest alphabet
    //     width that fits. The build itself is a separate pass through
    //     the genome that emits the typed `Vec<S>` for the chosen S.
    let n_seg = count_spacer_runs(&genome.sequence[..n2]);
    let alphabet_max = SENTINEL_BASE as u32 + n_seg;
    log::info!("sa_build: counted {n_seg} per-segment sentinels (alphabet max = {alphabet_max})");

    // (2) Count the kept positions for the log line + the
    //     debug-only post-condition. We don't materialise the
    //     position list — caps-sa accepts a streaming predicate
    //     (`build_ext_mem_for_filter_*`) and walks it once to build
    //     a 1-bit-per-position bitmap + tiny popcount prefix-sum.
    //     On the human genome that's ~770 MB total vs the 47 GB
    //     `Vec<u64>` the previous `_for_positions` path required
    //     (~60× reduction).
    let sparse_d: u64 = 1;
    debug_assert_eq!(
        sparse_d, 1,
        "non-default sparse_d isn't wired through this path"
    );
    let n_sa_kept: usize = genome.sequence[..n2].par_iter().filter(|&&b| b < 4).count();
    log::info!("sa_build: {n_sa_kept} entries after ACGT + sparse-d={sparse_d} filter");

    let n_genome_u64 = n_genome as u64;
    let mut emit_count: usize = 0;

    // The packer is shared between the in-memory and ext-mem paths.
    // With caps-sa's filter API every emitted position is already
    // a kept ACGT position (in original coordinates after the
    // spacer-free → original translation done inside
    // `dispatch_caps_sa_segmented`), so the only work here is the
    // strand-bit encoding + forwarding to `emit`. caps-sa 0.6's
    // `try_*` APIs propagate this callback error immediately, so
    // SA-file write failures abort the build instead of running to
    // completion and reporting afterwards.
    let mut pack_one = |orig_pos: u64| -> Result<(), Error> {
        let packed_value = if (orig_pos as usize) < n_genome {
            orig_pos
        } else {
            (orig_pos - n_genome_u64) | n2_bit
        };
        emit(packed_value)?;
        emit_count += 1;
        Ok(())
    };

    // The default path is the segmented arm: caps-sa's
    //     `SegmentedText` over a spacer-free `u8` text, wrapped in
    //     `StarSegmentedText` to flip the boundary-order to STAR's
    //     `spacer-as-largest` convention. This is byte-for-byte
    //     STAR-compatible (see `segmented_arm_matches_sentinel_arm_*`
    //     tests), uses less memory than the sentinel-transform arms
    //     (no spacers in the text, no per-segment sentinel bytes),
    //     and handles any segment count regardless of alphabet size —
    //     so genome+SJ indexes with tens of thousands of junctions
    //     "just work" with no width dispatch.
    //
    //     Setting `RUSTAR_USE_SENTINEL_TRANSFORM=1` (or passing
    //     `force_sentinel = true` to [`build_impl`] from tests) opts
    //     into the legacy sentinel-transform arms when the alphabet
    //     fits — useful as a fallback and as a cross-checking oracle
    //     in differential tests.
    if force_sentinel && alphabet_max <= <u8 as SaSymbol>::MAX_REPRESENTABLE {
        log::info!(
            "sa_build: RUSTAR_USE_SENTINEL_TRANSFORM=1, alphabet fits u8 — \
             using sentinel-transform arm"
        );
        let t_prime: Vec<u8> = build_sentinel_transformed_text(&genome.sequence[..n2], n_seg);
        dispatch_caps_sa(t_prime, &genome.sequence[..n2], temp_dir, &mut pack_one)?;
    } else if force_sentinel && alphabet_max <= <u16 as SaSymbol>::MAX_REPRESENTABLE {
        log::info!(
            "sa_build: RUSTAR_USE_SENTINEL_TRANSFORM=1, alphabet fits u16 — \
             using sentinel-transform arm"
        );
        let t_prime: Vec<u16> = build_sentinel_transformed_text(&genome.sequence[..n2], n_seg);
        dispatch_caps_sa(t_prime, &genome.sequence[..n2], temp_dir, &mut pack_one)?;
    } else {
        if force_sentinel {
            log::warn!(
                "sa_build: RUSTAR_USE_SENTINEL_TRANSFORM=1 requested but \
                 alphabet_max={alphabet_max} exceeds the u16 limit ({}) — \
                 falling back to the segmented arm",
                u16::MAX
            );
        } else {
            log::info!(
                "sa_build: using segmented arm (default; \
                 alphabet_max={alphabet_max}, {n_seg} segments)"
            );
        }
        dispatch_caps_sa_segmented(&genome.sequence[..n2], temp_dir, &mut pack_one)?;
    }

    debug_assert_eq!(
        emit_count, n_sa_kept,
        "sa_build: emitted {emit_count} entries but counted {n_sa_kept}",
    );
    Ok(emit_count)
}

/// Record describing one segment in both the original (spacer-bordered)
/// coordinate system and the spacer-free coordinate system caps-sa sees.
/// The spacer-free length is `sf_end - sf_start`; the original-text
/// length is `orig_end - orig_start` and includes any internal N bytes
/// (which are real symbols, not boundaries) but excludes the spacer
/// run that follows the segment.
/// Legacy spacer-free segment record. Production builds no longer
/// use the spacer-free copy ([`dispatch_caps_sa_segmented`] hands
/// `genome.sequence` directly to caps-sa with a segment-end list in
/// original coords); the helpers below are kept `#[cfg(test)]`-only
/// as part of the test-side oracle that the spacer-free → original
/// coordinate mapping is correct.
#[cfg(test)]
#[derive(Clone, Debug)]
struct SpacerFreeSegment {
    orig_start: u64,
    sf_start: u64,
    sf_end: u64,
}

/// Walk the original spacer-bordered text and produce
///
/// 1. a fresh `Vec<u8>` containing only the non-spacer bytes (real
///    ACGTN positions), preserving order, and
/// 2. a sorted list of [`SpacerFreeSegment`]s — one per maximal run
///    of non-spacer bytes — giving the orig↔spacer-free mapping.
///
/// Segments are split exactly where spacer bytes appear: a segment is
/// a maximal run of `text[p] != SPACER`. N bytes (value 4) stay inside
/// the segment they were in; only spacer bytes (value 5) end one.
#[cfg(test)]
fn build_spacer_free(original: &[u8]) -> (Vec<u8>, Vec<SpacerFreeSegment>) {
    let mut text_sf = Vec::with_capacity(original.len());
    let mut segments: Vec<SpacerFreeSegment> = Vec::new();
    let mut cur: Option<(u64, u64)> = None; // (orig_start, sf_start)
    for (i, &b) in original.iter().enumerate() {
        if b == SPACER {
            if let Some((orig_start, sf_start)) = cur.take() {
                segments.push(SpacerFreeSegment {
                    orig_start,
                    sf_start,
                    sf_end: text_sf.len() as u64,
                });
                let _ = orig_start; // silence warning if log feature disabled
                let _ = sf_start;
            }
        } else {
            if cur.is_none() {
                cur = Some((i as u64, text_sf.len() as u64));
            }
            text_sf.push(b);
        }
    }
    if let Some((orig_start, sf_start)) = cur {
        segments.push(SpacerFreeSegment {
            orig_start,
            sf_start,
            sf_end: text_sf.len() as u64,
        });
        let _ = orig_start;
        let _ = sf_start;
    }
    (text_sf, segments)
}

/// Map a position in the spacer-free coordinate system back to its
/// original (spacer-bordered) position. Binary-searches the segments
/// by their `sf_end` to find the containing segment in `O(log n_seg)`.
#[cfg(test)]
fn sf_to_orig(segments: &[SpacerFreeSegment], sf_pos: u64) -> u64 {
    // First segment with `sf_end > sf_pos` contains this position.
    let i = segments.partition_point(|s| s.sf_end <= sf_pos);
    debug_assert!(
        i < segments.len(),
        "sf_to_orig: sf_pos {sf_pos} past end of spacer-free text"
    );
    let seg = &segments[i];
    seg.orig_start + (sf_pos - seg.sf_start)
}

/// caps-sa [`LimitProvider`][caps_sa::LimitProvider] that wraps a
/// [`caps_sa::SegmentedText`] and overrides `boundary_order` to
/// match STAR's `spacer-as-largest` comparator at cross-segment
/// ties: the suffix that hits its boundary first is *larger*,
/// equivalently the longer-`lim` one is smaller, with an
/// ascending-position tie-break when both `lim`s coincide.
///
/// Spacer-free positions are sufficient for the position tie-break:
/// the `orig ↔ spacer-free` mapping ([`sf_to_orig`]) is monotonic
/// on every segment, so two positions' relative order is preserved
/// across the translation. A tie-break by spacer-free `p_a` / `p_b`
/// gives the same outcome as a tie-break by their original
/// counterparts.
struct StarSegmentedText {
    inner: caps_sa::SegmentedText,
}

impl caps_sa::LimitProvider for StarSegmentedText {
    #[inline]
    fn lim_at(&self, p: usize) -> usize {
        self.inner.lim_at(p)
    }

    #[inline]
    fn boundary_order(
        &self,
        p_a: usize,
        lim_a: usize,
        p_b: usize,
        lim_b: usize,
    ) -> std::cmp::Ordering {
        // longer-`lim` is smaller (STAR's spacer-as-largest);
        // position tie-break on equal `lim`s (STAR's ascending-`p`
        // convention).
        lim_b.cmp(&lim_a).then(p_a.cmp(&p_b))
    }
}

/// Drive caps-sa over the **original** spacer-bordered `&[u8]` text
/// with a [`StarSegmentedText`] limit provider whose segment ends
/// are the original-coordinate spacer positions. Byte-for-byte
/// STAR-compatible at cross-segment ties.
///
/// The previous implementation built a spacer-free copy of the text
/// (`text_sf`, a fresh ~6.3 GB `Vec<u8>` on the human genome) so
/// caps-sa's `SegmentedText` could model contiguous segments
/// without sentinel bytes between them. That copy is unnecessary —
/// the predicate `|p| genome.sequence[p] < 4` already filters
/// spacer (5) and N (4) positions out of the sort, so caps-sa's
/// `lim_at` is **only ever called on ACGT positions**. For those
/// positions, returning `next_spacer_pos - p` from the original
/// text's segment-ends list gives exactly the same `max_ctx` the
/// spacer-free path computed, and the LCP scans inside `lim_at`'s
/// bound never touch a spacer byte. The result is byte-identical
/// to the spacer-free path (verified by the
/// `segmented_arm_matches_sentinel_arm_byte_for_byte_*` differential
/// tests) and saves the 6.3 GB copy.
///
/// Each emitted SA position is already in original coordinates;
/// `pack_one` is the only thing needed before bit-packing into the
/// output sink.
fn dispatch_caps_sa_segmented(
    original: &[u8],
    temp_dir: Option<&Path>,
    mut pack_one: impl FnMut(u64) -> Result<(), Error>,
) -> Result<(), Error> {
    let n = original.len();
    if n == 0 {
        return Ok(());
    }

    // Segment ends in original coords: for each maximal non-spacer
    // run, the index of the first spacer byte after it (or `n` for
    // the trailing run). caps-sa's `SegmentedText::from_ends`
    // requires `ends.last() == text_len`, so we always close out
    // with `n` even if the text ended on a spacer (treating the
    // trailing spacer region as a nominal "segment" — never
    // queried because the predicate filters out spacer positions).
    let ends_orig = compute_spacer_ends(original);
    let n_seg_runs = if ends_orig.last() == Some(&(n as u64))
        && (ends_orig.len() < 2 || ends_orig[ends_orig.len() - 2] != n as u64 - 1)
    {
        ends_orig.len()
    } else {
        ends_orig.len().saturating_sub(1)
    };
    let lp = StarSegmentedText {
        inner: caps_sa::SegmentedText::from_ends(n, ends_orig),
    };
    log::info!(
        "sa_build: invoking caps-sa segmented filter path \
         (text len {n}, {n_seg_runs} non-spacer segments, no spacer-free copy)"
    );

    if use_ext_mem(n) {
        let opts = caps_sa_ext_mem_opts(temp_dir);
        // Predicate accepts ACGT only (rejects N at 4, spacer at 5).
        // Borrows `original` via `&[u8]` — `Send + Sync` is satisfied.
        let original_ref: &[u8] = original;
        caps_sa::try_build_ext_mem_for_filter_with(
            original,
            |p| original_ref[p as usize] < 4,
            &lp,
            &opts,
            &mut pack_one,
        )
        .map_err(map_caps_sa_error)?;
    } else {
        // Small-input in-memory path. Materialising the position
        // list at this scale is harmless (≤ 16 MB text → ≤ ~100 K
        // kept positions); keeps the in-mem fast-path simple.
        let positions: Vec<u64> = (0..n as u64)
            .filter(|&p| original[p as usize] < 4)
            .collect();
        let sa: Vec<u64> = caps_sa::build_in_memory_for_positions_with(
            original,
            positions,
            &lp,
            &caps_sa::Opts::default(),
        );
        for &orig_pos in &sa {
            pack_one(orig_pos)?;
        }
    }
    Ok(())
}

/// Compute segment ends in original-coord space — for each maximal
/// non-spacer run, the position of the first spacer byte after the
/// run (or `text.len()` for the trailing run). Always closes out
/// with `text.len()` so [`caps_sa::SegmentedText::from_ends`]
/// accepts the result (its constructor requires
/// `ends.last() == text_len`).
fn compute_spacer_ends(text: &[u8]) -> Vec<u64> {
    let mut ends: Vec<u64> = Vec::new();
    let mut in_seg = false;
    for (i, &b) in text.iter().enumerate() {
        if b == SPACER {
            if in_seg {
                ends.push(i as u64);
                in_seg = false;
            }
        } else {
            in_seg = true;
        }
    }
    let n = text.len() as u64;
    if ends.last() != Some(&n) {
        ends.push(n);
    }
    ends
}

/// Drive caps-sa over a typed `&[S]` sentinel-transformed text.
/// Encapsulates the in-mem vs ext-mem branch so both alphabet widths
/// share one path. Uses the streaming filter API in ext-mem mode —
/// the predicate `|p| p < n2 && original[p] < 4` selects ACGT-only
/// positions and skips the per-segment sentinel bytes + terminal
/// sentinel; caps-sa builds a ~`n / 8`-byte bitmap + popcount
/// prefix-sum internally and never materialises a `Vec<u64>`.
fn dispatch_caps_sa<S>(
    t_prime: Vec<S>,
    original: &[u8],
    temp_dir: Option<&Path>,
    mut pack_one: impl FnMut(u64) -> Result<(), Error>,
) -> Result<(), Error>
where
    S: caps_sa::Symbol,
{
    let n2 = original.len();
    let symbol_width = std::mem::size_of::<S>();
    if use_ext_mem(t_prime.len()) {
        log::info!(
            "sa_build: invoking caps-sa::build_ext_mem_for_filter \
             (text len {}, {symbol_width}-byte alphabet)",
            t_prime.len()
        );
        let opts = caps_sa_ext_mem_opts(temp_dir);
        // Predicate must skip both the per-segment sentinels (bytes
        // ≥ 5 in `t_prime`, encoded at the spacer-run positions in
        // `original`) and the terminal sentinel at index `n2`. The
        // shorthand `p < n2 && original[p] < 4` covers both: the
        // terminal-sentinel index n2 fails the first guard and the
        // spacer indices (where `original[p] == 5`) fail the second.
        caps_sa::try_build_ext_mem_for_filter(
            &t_prime,
            |p| (p as usize) < n2 && original[p as usize] < 4,
            &opts,
            &mut pack_one,
        )
        .map_err(map_caps_sa_error)?;
        drop(t_prime);
    } else {
        log::info!(
            "sa_build: invoking caps-sa::build_in_memory_for_positions \
             (text len {}, {symbol_width}-byte alphabet, small input)",
            t_prime.len()
        );
        // Small-input in-memory path — see the sibling segmented arm
        // for the rationale on keeping `_for_positions` here.
        let positions: Vec<u64> = (0..n2 as u64)
            .filter(|&p| original[p as usize] < 4)
            .collect();
        let sa: Vec<u64> = caps_sa::build_in_memory_for_positions(&t_prime, positions);
        drop(t_prime);
        for &sa_pos in &sa {
            pack_one(sa_pos)?;
        }
    }
    Ok(())
}

fn caps_sa_ext_mem_opts(temp_dir: Option<&Path>) -> caps_sa::ExtMemOpts {
    let mut opts = caps_sa::ExtMemOpts::from_env();
    if let Some(dir) = temp_dir {
        opts = opts.work_dir(dir);
    } else if let Some(dir) =
        std::env::var_os("RUSTAR_TMPDIR").or_else(|| std::env::var_os("RUSTAR_TEMP_DIR"))
    {
        opts = opts.work_dir(dir);
    }
    opts
}

fn map_caps_sa_error(err: caps_sa::BuildError<Error>) -> Error {
    match err {
        caps_sa::BuildError::Io(e) => {
            Error::Index(format!("caps-sa external-memory I/O failed: {e}"))
        }
        caps_sa::BuildError::Emit(e) => e,
    }
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
/// > then we filter them out, but the sort itself does `O(spacer_run_len²)`
/// > work because every spacer-starting suffix shares a near-maximal LCP
/// > with every other).
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
    if let Ok(v) = std::env::var("RUSTAR_USE_IN_MEM")
        && matches!(v.as_str(), "1" | "true" | "yes" | "on")
    {
        return false;
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

/// Count the maximal runs of spacer bytes (value `5`) in `genome`.
/// One sentinel value is reserved per run, so this is the number that
/// drives the alphabet-width decision in [`build`].
fn count_spacer_runs(genome: &[u8]) -> u32 {
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
    n_seg
}

/// Build the per-segment sentinel-transformed text `T'` of length `n + 1`
/// at the chosen symbol width `S`.
///
/// Each maximal run of spacer bytes (value `5`) is numbered in position
/// order; the run at index `i` has every position stamped with the
/// sentinel value `SENTINEL_BASE + i`, encoded via `S::from_sentinel(i)`.
/// After the last byte we append one extra terminal sentinel
/// `SENTINEL_BASE + n_seg` (a value larger than every per-run sentinel),
/// so the final RC suffix has a sentinel-terminator and the `caps-sa`
/// "implicit smallest sentinel at the end" never affects any kept
/// suffix's order.
///
/// `n_seg` must be the value returned by [`count_spacer_runs`] for the
/// same input and must satisfy `SENTINEL_BASE + n_seg <= S::MAX_REPRESENTABLE`
/// — the caller in [`build`] dispatches on exactly that check.
fn build_sentinel_transformed_text<S: SaSymbol>(genome: &[u8], n_seg: u32) -> Vec<S> {
    let n = genome.len();
    debug_assert!(
        SENTINEL_BASE as u32 + n_seg <= S::MAX_REPRESENTABLE,
        "build_sentinel_transformed_text: alphabet width {} cannot represent \
         terminal sentinel {} (n_seg={n_seg})",
        std::mem::size_of::<S>(),
        SENTINEL_BASE as u32 + n_seg,
    );

    let mut out: Vec<S> = Vec::with_capacity(n + 1);
    let mut run_idx: u32 = 0;
    let mut in_run = false;
    for &b in genome {
        if b == SPACER {
            in_run = true;
            out.push(S::from_sentinel(run_idx));
        } else {
            if in_run {
                in_run = false;
                run_idx += 1;
            }
            out.push(S::from_base(b));
        }
    }
    if in_run {
        run_idx += 1;
    }
    debug_assert_eq!(run_idx, n_seg);

    // Terminal sentinel — distinct from every per-run value, larger than
    // every real base, and the unique-maximum symbol of `T'`.
    out.push(S::from_sentinel(n_seg));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sentinel_transform_counts_runs() {
        // [a, 5, 5, c, 5, 5, 5, t]
        let input: Vec<u8> = vec![0, 5, 5, 1, 5, 5, 5, 3];
        let n_seg = count_spacer_runs(&input);
        assert_eq!(n_seg, 2);
        let out: Vec<u8> = build_sentinel_transformed_text(&input, n_seg);
        // Run 0 → sentinel 5; run 1 → sentinel 6; terminal sentinel = 7.
        assert_eq!(out, vec![0, 5, 5, 1, 6, 6, 6, 3, 7]);
    }

    #[test]
    fn sentinel_transform_trailing_run() {
        // Input ending in a spacer run.
        let input: Vec<u8> = vec![0, 1, 5, 5];
        let n_seg = count_spacer_runs(&input);
        assert_eq!(n_seg, 1);
        let out: Vec<u8> = build_sentinel_transformed_text(&input, n_seg);
        // Run 0 → sentinel 5; terminal sentinel = 6.
        assert_eq!(out, vec![0, 1, 5, 5, 6]);
    }

    #[test]
    fn sentinel_transform_no_spacers() {
        let input: Vec<u8> = vec![0, 1, 2, 3];
        let n_seg = count_spacer_runs(&input);
        assert_eq!(n_seg, 0);
        let out: Vec<u8> = build_sentinel_transformed_text(&input, n_seg);
        assert_eq!(out, vec![0, 1, 2, 3, 5]);
    }

    /// Inputs requiring more than `u8::MAX - SENTINEL_BASE = 250`
    /// per-run sentinels (251 spacer runs → terminal sentinel value
    /// 256) overflow the byte alphabet. The dispatch in [`build`]
    /// switches to the `u16` path here; this test exercises the
    /// `Vec<u16>` builder directly.
    #[test]
    fn sentinel_transform_u16_handles_many_runs() {
        // 300 spacer runs separated by single bases → n_seg = 300,
        // terminal sentinel = 305 → needs u16.
        let mut input: Vec<u8> = Vec::new();
        for _ in 0..300 {
            input.push(0); // real base
            input.push(5); // spacer
        }
        let n_seg = count_spacer_runs(&input);
        assert_eq!(n_seg, 300);

        let out: Vec<u16> = build_sentinel_transformed_text(&input, n_seg);
        assert_eq!(out.len(), input.len() + 1);
        // First spacer run got sentinel SENTINEL_BASE + 0 = 5.
        assert_eq!(out[1], SENTINEL_BASE as u16);
        // Last spacer run (idx 299) got SENTINEL_BASE + 299 = 304.
        assert_eq!(out[599], SENTINEL_BASE as u16 + 299);
        // Terminal sentinel at out[2*300] = SENTINEL_BASE + 300 = 305.
        assert_eq!(*out.last().unwrap(), SENTINEL_BASE as u16 + 300);
    }

    /// The `MAX_REPRESENTABLE` constants and the `From*` casts on
    /// `SaSymbol` line up with what the dispatch in [`build`] expects.
    #[test]
    fn sa_symbol_widths() {
        assert_eq!(<u8 as SaSymbol>::MAX_REPRESENTABLE, 255);
        assert_eq!(<u16 as SaSymbol>::MAX_REPRESENTABLE, 65535);
        assert_eq!(<u8 as SaSymbol>::from_base(2), 2u8);
        assert_eq!(<u16 as SaSymbol>::from_base(2), 2u16);
        assert_eq!(<u8 as SaSymbol>::from_sentinel(3), 8u8); // 5+3
        assert_eq!(<u16 as SaSymbol>::from_sentinel(3000), 3005u16);
    }

    /// Spacer-free transform should drop every spacer byte and emit
    /// one [`SpacerFreeSegment`] per maximal non-spacer run, with the
    /// `orig_start`/`sf_start` mapping consistent.
    #[test]
    fn spacer_free_basic() {
        // Original: [a, c, 5, t, 5, 5, g, n] → 3 segments: [a,c], [t], [g,n]
        let original = vec![0u8, 1, 5, 3, 5, 5, 2, 4];
        let (text_sf, segments) = build_spacer_free(&original);
        assert_eq!(text_sf, vec![0u8, 1, 3, 2, 4]);
        assert_eq!(segments.len(), 3);
        assert_eq!(segments[0].orig_start, 0);
        assert_eq!(segments[0].sf_start, 0);
        assert_eq!(segments[0].sf_end, 2);
        assert_eq!(segments[1].orig_start, 3);
        assert_eq!(segments[1].sf_start, 2);
        assert_eq!(segments[1].sf_end, 3);
        assert_eq!(segments[2].orig_start, 6);
        assert_eq!(segments[2].sf_start, 3);
        assert_eq!(segments[2].sf_end, 5);
    }

    /// `sf_to_orig` must invert `build_spacer_free` exactly: every
    /// non-spacer byte's sf-position maps back to its original
    /// position.
    #[test]
    fn sf_to_orig_round_trip() {
        let original = vec![0u8, 1, 5, 3, 5, 5, 2, 4, 5, 1, 3];
        let (text_sf, segments) = build_spacer_free(&original);
        let mut sf_pos: u64 = 0;
        for (orig_idx, &b) in original.iter().enumerate() {
            if b == SPACER {
                continue;
            }
            assert_eq!(text_sf[sf_pos as usize], b);
            assert_eq!(sf_to_orig(&segments, sf_pos), orig_idx as u64);
            sf_pos += 1;
        }
        assert_eq!(sf_pos as usize, text_sf.len());
    }

    /// Build a small multi-chromosome genome from an inline FASTA
    /// string for the differential tests below. Routes through the
    /// same `Genome::from_fasta` path the rest of the test suite
    /// uses, so the result is byte-identical to what production
    /// `build()` would see.
    fn build_genome_from_fasta(fasta: &str, bin_nbits: u32) -> crate::genome::Genome {
        use crate::params::Parameters;
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(fasta.as_bytes()).unwrap();
        let bin_nbits_str = bin_nbits.to_string();
        let args = vec![
            "rustar-aligner",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            &bin_nbits_str,
        ];
        let params = Parameters::parse_from(args);
        crate::genome::Genome::from_fasta(&params).unwrap()
    }

    /// With [`StarSegmentedText`] flipping the boundary-order
    /// convention to STAR's, the segmented arm's packed array must
    /// be **byte-for-byte identical** to the `Vec<u8>` /
    /// `Vec<u16>` sentinel-transform arm — both implement STAR's
    /// `spacer-as-largest` order on the same kept ACGT positions.
    ///
    /// Driven through [`build_impl`] directly to avoid racing with
    /// parallel tests on the env var.
    fn assert_arms_byte_identical(label: &str, fasta: &str, bin_nbits: u32) {
        let genome = build_genome_from_fasta(fasta, bin_nbits);
        // Default (segmented) arm — the new production path.
        let sa_segmented = build_impl(&genome, false).unwrap();
        // Sentinel-transform fallback — the legacy STAR-faithful path.
        let sa_sentinel = build_impl(&genome, true).unwrap();
        assert_eq!(
            sa_segmented.data.data(),
            sa_sentinel.data.data(),
            "segmented vs sentinel-transform packed array differ on `{label}`"
        );
        assert_eq!(sa_segmented.gstrand_bit, sa_sentinel.gstrand_bit);
    }

    /// The streaming entry [`build_streaming`] must produce the same
    /// byte sequence as the in-memory [`build`] when its packed
    /// values are written through [`PackedStreamWriter`]. Drives both
    /// over a small multi-chromosome fixture and compares the byte
    /// output of each.
    #[test]
    fn streaming_build_matches_in_memory_build() {
        let genome =
            build_genome_from_fasta(">chrA\nACGTACGTAC\n>chrB\nGGGGCCCC\n>chrC\nNNACGTNN\n", 4);
        let in_mem = build_impl(&genome, false).unwrap();

        // Streaming: collect packed values into the same byte layout
        // a `PackedArray` would have, via [`PackedStreamWriter`].
        use crate::index::packed_stream::PackedStreamWriter;
        let word_length = in_mem.data.word_length();
        let mut got: Vec<u8> = Vec::new();
        let mut writer = PackedStreamWriter::new(&mut got, word_length);
        let (gbit, gmask, n) = super::build_streaming(&genome, None, |pv| {
            writer.write_one(pv).unwrap();
            Ok(())
        })
        .unwrap();
        let _w = writer.finish().unwrap();
        assert_eq!(gbit, in_mem.gstrand_bit);
        assert_eq!(gmask, in_mem.gstrand_mask);
        assert_eq!(n, in_mem.data.len());
        assert_eq!(
            got,
            in_mem.data.data(),
            "streaming build's byte stream differs from in-memory build"
        );
    }

    #[test]
    fn segmented_arm_matches_sentinel_arm_byte_for_byte_three_chrs() {
        // Three chromosomes including one with internal N's.
        assert_arms_byte_identical(
            "three-chr-with-Ns",
            ">chrA\nACGTACGTAC\n>chrB\nGGGGCCCC\n>chrC\nNNACGTNN\n",
            4,
        );
    }

    #[test]
    fn segmented_arm_matches_sentinel_arm_byte_for_byte_single_chr() {
        // Single chromosome, no inter-segment spacer (within the
        // same chromosome the sentinel byte still appears at the
        // chromosome's bin-padding boundary). Exercises the
        // longer-is-smaller convention on prefix-relationships
        // within one segment.
        assert_arms_byte_identical("single-chr", ">chr1\nACGTACGTACGTNACGT\n", 4);
    }

    #[test]
    fn segmented_arm_matches_sentinel_arm_byte_for_byte_many_short_chrs() {
        // ~30 short chromosomes — exercises many cross-segment
        // boundary tie-breaks. Still well within the u8 alphabet
        // (~30 sentinels) so the sentinel arm picks the u8 path.
        use std::fmt::Write as _;
        let mut fasta = String::new();
        for i in 0..30 {
            write!(fasta, ">chr{i}\nACGT{}\n", "AC".repeat(i % 5 + 1)).unwrap();
        }
        assert_arms_byte_identical("many-short-chrs", &fasta, 4);
    }

    #[test]
    fn segmented_arm_matches_sentinel_arm_byte_for_byte_repeats() {
        // Chromosomes with shared prefixes — exercises the cross-
        // segment LCP-truncation case where one suffix hits its
        // segment boundary while the other still has bytes left.
        assert_arms_byte_identical(
            "shared-prefixes",
            ">chrA\nACGTACGT\n>chrB\nACGTAC\n>chrC\nACGT\n",
            4,
        );
    }

    #[test]
    fn segmented_arm_matches_sentinel_arm_byte_for_byte_n_heavy() {
        // Multiple N runs inside one chromosome (which are *not*
        // segment boundaries — N is a real base, only the spacer
        // byte 5 ends a segment) so the segment list stays small
        // but kept-position filtering does real work.
        assert_arms_byte_identical(
            "n-heavy",
            ">chr1\nACGTNNNACGTNNNACGTNNN\n>chr2\nNNNACGT\n",
            5,
        );
    }
}
