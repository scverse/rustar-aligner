# rustar-aligner Changelog

<!--
Release notes are extracted from this file by the release workflow.
Each released version needs a heading of the form:

    ## [Version X.Y.Z](https://github.com/scverse/rustar-aligner/releases/tag/vX.Y.Z) - YYYY-MM-DD

Sections commonly used: Features, Bug fixes, Other changes.
-->

## [Unreleased]

### Features

- **`--clipAdapterType CellRanger4` now matches STAR exactly**, both halves of the
  clip. The 5' TSO trim is an overlap alignment against the first 91 bases of the
  read, replicating `ClipCR4`'s Opal call (ACGTN alphabet; match +1, mismatch -2,
  any-vs-N -2, N-vs-N 0; gap open = gap extend = 2; overlap mode with end
  tracking), followed by STAR's acceptance gate
  `S<20 || (S==20 && L>26) || (S==21 && L>30)`. The 3' trim is STAR's scored
  `polyTail3p` scan, replacing a "trailing run of A >= 8" approximation that was
  wrong in both directions (a read ending in 10 A's was trimmed by 10 where STAR
  trims 0; `A*15 + C + A*15` was trimmed by 15 where STAR trims 31). The poly-A
  transcription is taken from Benjamin Demaille's #148 and moved onto the solo
  path, which is where `CellRanger4` is actually reached; #148's `clip_mate`
  wiring is not used and that PR is superseded.

  The overlap alignment comes from [`hyalite`](https://crates.io/crates/hyalite)
  0.2 — a new dependency, pure Rust with no dependencies of its own. Reads are
  scanned a batch at a time via `Database::scan_all`, mirroring STAR's
  `ClipMate::clipChunk`, which aligns the adapter against a chunk of reads in one
  Opal call.

  Both halves are gated against STAR's own C++: `tests/data/cr4_opal_oracle.cpp`
  links STAR's `opal.cpp` and reproduces `ClipCR4` + `ClipMate::clipChunk`
  verbatim, and the committed `cr4_opal_oracle.tsv` is its output over 938 reads
  chosen to straddle both decision boundaries. `cr4_tso_matches_star_opal_oracle`
  checks every one against the scalar path, a forced-SIMD path, and the
  production batch path.

  Measured on 10x mouse chr19 with both tools under `--clipAdapterType
  CellRanger4`: clip-amount differences against STARsolo drop from 154 to **0**,
  and the soft-clipped share moves from 39% (STAR 41%) to 34.9% (STAR 35.0%).
  Default SE/PE alignment is untouched, as `CellRanger4` is opt-in: SE 8788/8926
  and PE 8390 both-mapped / 0 half-mapped, both unchanged.

  Known limitation: `CellRanger4` combined with a non-zero `--clip5pNbases` /
  `--clip3pNbases` misplaces most reads by `clip5pNbases`. That is a pre-existing
  bug, not introduced here, and is tracked separately; `CellRanger4` on its own
  and the fixed clips on their own are both unaffected.

- **STARsolo single-cell quantification (`--soloType`)** — the 10x
  Chromium / plate-based count-matrix pipeline, ported from STAR and
  verified against real STARsolo (#90).

  - Chemistries: `CB_UMI_Simple` (10x 3'/5'), `CB_UMI_Complex`
    (multi-segment CB), and `SmartSeq` (manifest-driven, SE read counts
    / PE fragment counts).
  - Features: `Gene`, `GeneFull` (pre-mRNA), `SJ` (splice-junction
    counts), and `Velocyto` (spliced/unspliced/ambiguous matrices for
    RNA velocity).
  - Barcode correction: all `--soloCBmatchWLtype` modes
    (`Exact`/`1MM`/`1MM_multi`/`1MM_multi_pseudocounts`/
    `1MM_multi_Nbase_pseudocounts`); all `--soloUMIdedup` and
    `--soloUMIfiltering` modes; real-valued `--soloMultiMappers`
    (`Uniform`/`PropUnique`/`Rescue`/`EM`).
  - Output: `raw/` + `filtered/` matrices, `EmptyDrops_CR` cell calling,
    CellRanger-style `Summary.csv`, and `--soloOutGzip`.

- **`--outSAMattributes GX GN`** — per-read gene-id / gene-name SAM tags
  for the solo `Gene` assignment.

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

### Bug fixes

- **STARsolo `Gene` assignment now requires exon concordance**, matching
  STARsolo: a read counts toward a gene only when every aligned block
  lies within the gene's exons, rather than merely overlapping one. This
  removes over-assignment of reads that extend past an exon boundary into
  an intron. On the mouse-chr19 differential test the `Gene` count matrix
  now matches STAR 2.7.11b (raw matrix identical apart from the
  multimapper tie-break tail that even a byte-faithful reimplementation
  cannot reproduce).
- **STARsolo `features.tsv` column 2 now emits the GTF `gene_name`**
  (symbol), with the STAR gene_id fallback, instead of duplicating the
  gene_id.

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
