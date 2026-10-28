# Dependencies

Every crate this project depends on, and every crate it decided not to depend on.

[CONTRIBUTING.md](CONTRIBUTING.md#new-dependencies-need-prior-discussion) requires that a new
dependency be raised in an issue before the PR. That rule covers the *decision*; this file is the
*record*, so a crate that was considered and turned down does not get re-proposed and re-litigated
from scratch six months later.

Licenses and versions below are the resolved ones from `cargo metadata`, and every "used by" entry
names a file that actually imports the crate. Both are checked by hand when this file changes; if
you bump a dependency, update the row.

## 1. Runtime dependencies

These ship. The bar is high: a crate here has to do something the standard library does not, on all
five supported platforms, without a C toolchain unless the row says otherwise.

| Crate | Version | License | Used by | Why it earns its place |
|---|---|---|---|---|
| `clap` | 4 | MIT OR Apache-2.0 | `params/` | STAR's parameter surface is ~200 flags with per-flag arity and defaults. Hand-rolling that parser would be more code than the aligner core. |
| `anyhow` | 1 | MIT OR Apache-2.0 | `main.rs`, `lib.rs`, the run-mode entry points | Top-level error propagation, where the caller only prints. |
| `thiserror` | 2 | MIT OR Apache-2.0 | `error.rs` | The typed `Error` enum the library returns. |
| `log` + `env_logger` | 0.4 / 0.11 | MIT OR Apache-2.0 | everywhere | `Log.out` and the progress logs. |
| `memmap2` | 0.9 | MIT OR Apache-2.0 | `index/packed_array.rs`, `genome/mod.rs`, `index/mod.rs` | The genome, SA and SAindex are read by mapping the files, not by reading them into a `Vec`. A human genome index is tens of gigabytes; mapping is what makes a shared index across processes possible at all. |
| `byteorder` | 1 | Unlicense OR MIT | `io/bam.rs`, `index/io.rs` | Little-endian field access in the BAM and index formats, both externally specified. |
| `noodles` | 0.116 | MIT | `io/`, `align/`, `chimeric/`, `quant/`, `solo/`, `wasp/`, `bam_dedup.rs` | SAM/BAM/BGZF/FASTQ readers and writers, pure Rust, self-contained on all five platforms. The alternative (`rust-htslib`) is in section 3. |
| `noodles-bgzf` | 0.51 | MIT | `io/bgzf_writer.rs` | BGZF block framing for BAM, with the `libdeflate` feature. |
| `libdeflater` | 1.25.2 | Apache-2.0 | `io/bgzf_writer.rs`, `solo/count.rs` | Deflate for the parallel BGZF writer and gzip for the solo matrix files. Faster than `flate2` on those paths, and already in the tree through `noodles-bgzf`'s `libdeflate` feature. |
| `flate2` (`zlib-rs` backend) | 1 | MIT OR Apache-2.0 | `io/fastq.rs`, `io/compression.rs`, `junction/gtf.rs`, `quant/coverage.rs`, `solo/whitelist.rs`, `solo/count.rs`, `bin/emptydrops.rs` | gzip FASTQ, GTF and whitelist input. The `zlib-rs` backend instead of the default `miniz_oxide`: 2-3x faster inflate and deflate on the decode and BGZF paths, still pure Rust, no C toolchain. Backend is a build-time choice; the API is unchanged. |
| `bstr` | 1 | MIT OR Apache-2.0 | `io/sam.rs`, `chimeric/bam.rs` | Byte strings for SAM fields, which are bytes and not guaranteed UTF-8. |
| `rayon` | 1 | MIT OR Apache-2.0 | `lib.rs`, `index/sa_build.rs`, `index/sa_index.rs`, `index/sjdb_otf.rs`, `solo/`, `degnorm/run.rs` | The per-read parallelism. |
| `rustc-hash` | 2 | Apache-2.0 OR MIT | `align/score.rs`, `align/seed.rs`, `align/stitch.rs`, `solo/count.rs` | A fast non-cryptographic hasher for the hot maps. The inputs are internal keys, never attacker-controlled. |
| `dashmap` | 6 | MIT | `junction/sj_output.rs` | The concurrent junction-count map. Its iteration order is not stable, which is why every path that emits an order sorts on a total key first (see [#210](https://github.com/scverse/rustar-aligner/issues/210), and the determinism test that locks it). |
| `chrono` | 0.4 | MIT OR Apache-2.0 | `lib.rs`, `stats.rs`, `io/log.rs` | Timestamps in `Log.out` / `Log.final.out`, whose format STAR fixes. |
| `tempfile` | 3 | MIT OR Apache-2.0 | `io/bam.rs`, `solo/count.rs` | The sorted-BAM spill runs and the solo matrix body, plus test fixtures. |
| `shlex` | 2 | MIT OR Apache-2.0 | `params/mod.rs` | Quoting the command line for the `@PG` header, and splitting `--readFilesCommand`. |
| `caps-sa` | 0.7 | MIT | `index/sa_build.rs` | Suffix-array construction at `genomeGenerate`. Alternatives are in section 4. |
| `mimalloc` + `libmimalloc-sys` | 0.1 | MIT | `main.rs`, `lib.rs` | The global allocator. glibc malloc keeps one arena per worker thread and caches freed allocations in them indefinitely, which added 10-20 GB of slack to peak RSS on genome-scale runs; mimalloc's per-thread heaps release whole segments back to the OS. It is also cheaper per allocation, which matters for the millions of small allocations SA construction makes. |
| `hyalite` | =0.4.0 | MIT | `solo/mod.rs` | The Smith-Waterman alignment behind `--clipAdapterType CellRanger4`'s TSO clip ([#197](https://github.com/scverse/rustar-aligner/issues/197)). Pinned to an exact version because its output is checked against a committed fixture. |

Optional, off by default (features `bz2`, `zstd`, `xz`, or all three as `compressed-input`). They read FASTQ input compressed with formats beyond gzip, which STAR does not ([#218](https://github.com/scverse/rustar-aligner/issues/218), done in [#272](https://github.com/scverse/rustar-aligner/pull/272)):

| Crate | Version | License | Used by | Why it earns its place |
|---|---|---|---|---|
| `bzip2` | 0.6 | MIT OR Apache-2.0 | `io/compression.rs`, `io/fastq.rs` | bzip2 input (feature `bz2`). Its backend, `libbz2-rs-sys`, is a pure-Rust port under the bzip2-1.0.6 licence, so no C toolchain. |
| `ruzstd` | 0.9.1 | MIT | `io/compression.rs`, `io/fastq.rs` | Zstandard input (feature `zstd`), pure Rust. |
| `lzma-rust2` | 0.21 | Apache-2.0 | `io/compression.rs`, `io/fastq.rs` | XZ input (feature `xz`), pure Rust. |

## 2. Development-only dependencies

Nothing here ships, so the bar is lower: it has to make a test clearer than the same test written
by hand.

| Crate | Version | License | Used by |
|---|---|---|---|
| `assert_cmd` | 2 | MIT OR Apache-2.0 | every integration test that runs the binary |
| `predicates` | 3 | MIT OR Apache-2.0 | `tests/phase9_threading.rs`, `tests/transcriptome_sam.rs` |
| `lzma-rust2` (`encoder`) | 0.21 | Apache-2.0 | the XZ input tests in `io/compression.rs` and `io/fastq.rs` |
| `chrono` (build) | 0.4 | MIT OR Apache-2.0 | `build.rs`, for the build timestamp |

## 3. Considered and declined

Declined means: raised, weighed, and turned down. Reopening one of these needs new evidence, not a
new opinion.

| Crate(s) | Decided | Reason |
|---|---|---|
| `block-aligner`, `ksw2rs`, `parasailors` | Aug 2026 survey ([#201](https://github.com/scverse/rustar-aligner/issues/201)) | Faithfulness to STAR's scoring, extension and tie-breaks *is* the product. An aligner that computes alignments differently, however much better, is a divergence, and the diff would show up as reads moving. |
| `simd-minimizers`, `minimizer-iter`, `minimizer-queue`, `seq-hash`, `sourmash`, `nthash` | Aug 2026 survey ([#203](https://github.com/scverse/rustar-aligner/issues/203)) | STAR seeds by maximal-mappable-prefix search in the suffix array. Sketching changes which seeds exist, which changes which alignments exist. |
| `rust-htslib` | Aug 2026 survey | Brings htslib as a system C dependency. `noodles` already covers SAM/BAM/BGZF with a self-contained build on all five supported platforms, and this crate publishes to crates.io. |
| `rust-bio` | Aug 2026 survey | A broad toolkit for a handful of needed pieces, with readers slower than `noodles` (extra allocations, copying, UTF-8 validation). |
| `sucds`, `vers-vecs`, `sux`, `bitm` | Aug 2026 survey | The packed index layout is dictated by STAR's on-disk format, and nothing queries it with rank/select. |
| `rkyv` | Aug 2026 survey | The index layout is externally specified by STAR, not a serialization of our own types. |
| `niffler` | Oct 2026 ([#218](https://github.com/scverse/rustar-aligner/issues/218)) | bz2 / zstd / xz input was added with the format crates directly, behind features (section 1), so the wrapper wasn't needed. |

## 4. Under evaluation

Open questions. Listed so that "still open" is not mistaken for either "accepted" or "declined",
and so a proposal arrives at the existing thread rather than a new one.

| Crate(s) | Issue | Question |
|---|---|---|
| `sufr` / `libsufr` | [#202](https://github.com/scverse/rustar-aligner/issues/202) | Suffix-array construction against the `caps-sa` incumbent. Needs a measured win on `genomeGenerate` wall time or peak RSS to flip. |
| `libsais` (vendored C via `cc`) | [#162](https://github.com/scverse/rustar-aligner/issues/162) | Same question, with a C toolchain in the build. PR [#109](https://github.com/scverse/rustar-aligner/pull/109) selects it by `--limitGenomeGenerateRAM`. |
| `wide` / `pulp` | [#205](https://github.com/scverse/rustar-aligner/issues/205) | A portable SIMD crate against the hand-rolled intrinsics in `align/simd_scan.rs`. |
| `superintervals`, `coitrees`, `rust-lapper` | [#208](https://github.com/scverse/rustar-aligner/issues/208) | Interval overlap against the in-tree max-end segment tree, including their non-standard license flags. |
| `rapidgzip` | [#224](https://github.com/scverse/rustar-aligner/issues/224) | Parallel gzip input behind an optional feature. PR [#225](https://github.com/scverse/rustar-aligner/pull/225). |

## Adding a dependency

1. Open an issue first, per [CONTRIBUTING.md](CONTRIBUTING.md#new-dependencies-need-prior-discussion). A non-Rust dependency (a `-sys` crate, `bindgen`/`libclang`, a system library) is a supply-chain and five-platform decision, not an implementation detail.
2. Say what it replaces and what it costs: build time, binary size, platforms, license, and how much of it you will actually use.
3. When the decision lands, add the row here — to section 1 or 2 if accepted, to section 3 if not. A declined crate with no row is a decision that will have to be made twice.
