# rustar-aligner Changelog

<!--
Release notes are extracted from this file by the release workflow.
Each released version needs a heading of the form:

    ## [Version X.Y.Z](https://github.com/scverse/rustar-aligner/releases/tag/vX.Y.Z) - YYYY-MM-DD

Sections commonly used: Features, Bug fixes, Other changes.
-->

## [Unreleased]

### Features

- **`genomeGenerate` peak RSS cut from ~113 GB → ~11 GB** on the human
  genome (GRCh38, 32 threads). The construction pipeline no longer
  materialises three large intermediates that were dominating the peak:

  - The ACGT-only **kept-positions `Vec<u64>`** (~47 GB on the human
    genome) — caps-sa 0.5's `build_ext_mem_for_filter` API takes a
    predicate over text positions instead, and internally maintains
    only a ~770 MB bitmap + popcount prefix sum.
  - The **spacer-free copy of `genome.sequence`** (~6.3 GB) — the new
    `dispatch_caps_sa_segmented` hands the **original** spacer-bordered
    `&genome.sequence[..n2]` to caps-sa with a
    `SegmentedText::from_ends(spacer_positions)` limit provider, so
    LCP scans stop at the next spacer without a copy.
  - The **in-RAM SA `PackedArray`** (~25 GB) — the SA streams directly
    to `genome_dir/SA` via the new `PackedStreamWriter` as each caps-sa
    entry is emitted. `SuffixArray::build` (in-memory) is retained for
    tests; `sa_build::build_streaming` is the production path.

- **SAindex now parallelises across all rayon workers**. Previously
  every SAindex k-mer extraction sat on caps-sa's single-threaded
  phase-4 emit loop (~16 min of pure serial work). The new
  `SaIndex::build_parallel` reads the on-disk SA via chunked `pread`
  (so SA pages live in kernel page cache, **not** process RSS) and
  atomic-mins each k-mer's first-occurrence `sa_idx` into a
  `Vec<AtomicU64>` — the final pack into the SAindex's `PackedArray`
  is a single fast sequential pass.

- **SAindex inner loop now matches STAR's `isaStep + binary-search`
  skip algorithm** (`genomeSAindex.cpp::genomeSAindexChunk` /
  `funSAiFindNextIndex`). Consecutive SA entries share monotonically
  non-decreasing k-mer prefixes; rather than visit every entry, each
  rayon worker jumps forward by `isa_step = nSA / 4^nbases` (≈ 22 on
  human) and only stops to record k-mer boundaries — binary-searching
  inside the last `isa_step` window when `(indFull, iL4)` changes.
  Per-worker `ind0_local[il]` tracks the last-written k-mer index at
  each level so we skip writes inside a constant-prefix run; cross-
  chunk merge still uses `fetch_min` on the shared `Vec<AtomicU64>`.
  Drops the SAindex phase from ~10:50 to ~8 s on the human genome
  (~80× speedup), making rustar-aligner's full `genomeGenerate`
  **faster than upstream STAR 2.7.11b** (7:36 vs 11:26 wall, 32
  threads).

- **SAindex absent-slot encoding now matches STAR's** (`next_present
  _sa_idx | absent_mask` for between-present gaps, `n_entries
  | absent_mask` for tail-gaps). A single backward pass per SAindex
  level fills the gaps in place inside `firsts[]` before encoding
  into the output `PackedArray`. `hierarchical_lookup` is unchanged
  because it only consults the absent flag bit, not the slot's
  value; this change makes the on-disk SAindex bytes closer to
  STAR's (N flag bit is still not tracked — that's the only remaining
  divergence and a `hierarchical_lookup`-irrelevant one).

- **rayon thread pool now bound to `--runThreadN` at the
  run() dispatcher** rather than only inside `align_reads`. On a
  256-core machine this drops the rayon pool from 256 (the
  `num_cpus::get()` default) to whatever `--runThreadN` says,
  eliminating ~15 GB of glibc-arena-style allocator slack from
  256 worker-thread heaps.

- **mimalloc as the global allocator**. Lower per-allocation cost than
  glibc's malloc and per-thread heaps that return whole segments to
  the OS when abandoned, so allocator cache size stays bounded.

### Bumps

- `caps-sa` → `0.5` (adds `build_ext_mem_for_filter*`; see the caps-sa
  v0.5.0 release notes).

### Other

- New module `index::packed_stream` — bit-for-bit-compatible streaming
  writer for STAR's `PackedArray` format. Used by the streaming SA
  emit; documented to match `PackedArray::write` for `word_length ≤ 57`
  (the upper bound where the existing `PackedArray::write` is
  truncation-free; STAR's production widths are 32-37).

- New module-level `GenomeIndex::generate_streaming(params)` — the
  full `genomeGenerate` pipeline, used by the `genomeGenerate`
  run-mode dispatcher. The in-memory `GenomeIndex::build` +
  `GenomeIndex::write` flow remains for tests and any caller that
  needs random access to the SA in RAM.

Initial release of Rust rewrite of STAR.
