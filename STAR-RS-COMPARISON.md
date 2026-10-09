# STAR-rs vs rustar-aligner — Comprehensive Review & Comparison

*Prepared 2026-07-24. STAR-rs reviewed at `/home/jamfer/repos/STAR-rs` (v3.1.0, commit `f475924`); rustar-aligner reviewed at `/home/jamfer/repos/rustar-aligner` (v0.1.0, branch `fix/log-out-files`). STAR-rs was built and smoke-tested on this machine; both codebases were reviewed at source level by a fan-out of specialist agents.*

---

## 1. Executive Summary

**Both projects are Rust reimplementations of Alexander Dobin's STAR RNA-seq aligner, but they are aiming at fundamentally different targets and are at very different stages of maturity.**

**STAR-rs** (Benjamin Demaille, INSERM/IPNP Paris) is a **byte-for-byte** clean-room reproduction of **STAR 2.7.11b**. Its binary is literally named `STAR` (+ `STARlong`), it recognises all **200** STAR parameters, and it reproduces STAR's output byte-identically, **independent of thread count**, across essentially the entire STAR feature surface — including the parts most ports never touch: **STARsolo** (single-cell, all UMI/CB-match modes, EmptyDrops_CR, velocyto), **WASP**, **STARlong**, genome transform (haploid/diploid), liftOver, `inputAlignmentsFromBAM`, SuperTranscriptome, signal/wig tracks, and chimeric detection in all output modes. It is a 10-crate workspace (~33k src LOC), version **3.1.0**, engineered around reproducibility as a first-class invariant: a dedicated provenance crate (SHA-256 input-graph sidecars), a pinned toolchain, single-codegen-unit/fat-LTO builds, and a **474-test differential harness that runs the real native STAR binary as an oracle** and diffs normalised SAM, backed by committed goldens. It reports 100.0000% byte-identical placement on a 1M-read-pair GRCh38 GIAB run. **I built it here (clean, 1m05s) and ran its suite: 464/474 tests pass; all 10 "failures" are environmental oracle artifacts on this Linux host — 9 are the glibc `-nan%` vs Apple-libc `nan%` `printf` difference (its goldens were locked against a macOS STAR oracle) and 1 is the native STAR oracle itself crashing on `--alignWindowsPerReadNmax 1`. None is a STAR-rs correctness defect.**

**rustar-aligner** (James Ferguson, Garvan) is a **faithful-*port*** of STAR that targets **near-exact statistical concordance** rather than byte-identity. It is a single crate (~30.5k LOC), version **0.1.0**, covering the **core aligner** thoroughly: SE/PE alignment, splice junctions, the annotated-junction database (with in-index sjdb insertion), indels, multimapping, two-pass, `GeneCounts`, `TranscriptomeSAM`, chimeric (SE+PE, 4 tiers), full BAM output (sorted/unsorted/stdout), unmapped FASTX, and BySJout. It **deliberately does not** implement STARsolo, WASP, STARlong, genome transform, liftOver or wig tracks (STARsolo is explicitly deferred to a future phase). Its defining algorithmic choice: it uses a **per-read-seeded `StdRng`** for multimapper tie-breaks instead of emulating STAR's `mt19937`, so it is near-exact on uniquely-mapping reads but picks a different repeat copy than STAR on ties. It reports **99.815% (SE) / 99.883% (PE) "tie-adjusted" faithfulness** on 10k yeast reads (i.e. excluding tie-break divergences), with **0 STAR-only and 0 rustar-only reads** — every disagreement is a documented tie. Its own suffix array is byte-identical to STAR's for the yeast genome. Testing is a large in-tree unit suite (~434 `#[test]`) plus integration tests and **Python differential scripts that consume real STAR output as ground truth**. Unlike STAR-rs, it is **backed by two major bioinformatics community organisations — scverse (the single-cell Python ecosystem) and nf-core (Nextflow pipelines)** — and carries a **full CI/CD pipeline** (per-push/PR build+test+fmt+clippy+`cargo audit`, an Astro docs site to GitHub Pages, multi-arch Docker + crates.io release automation, Dependabot). This institutional backing and automated verification are a significant maintainability and longevity advantage over STAR-rs's single-author project with no test CI.

### Head-to-head at a glance

| Dimension | **STAR-rs** | **rustar-aligner** |
|---|---|---|
| Version / maturity | 3.1.0, public releases | 0.1.0, pre-release |
| Author / org | Benjamin Demaille (INSERM/IPNP) | James Ferguson (Garvan) |
| License | MIT **OR** Apache-2.0 | MIT |
| Rust edition / toolchain | 2021, **pinned 1.96.0** | 2024, `rust-version 1.88` |
| Structure | **10-crate workspace** | single crate + binary |
| Source size | ~33k src (~55k incl. tests) | ~30.5k |
| **Fidelity target** | **byte-for-byte** vs STAR 2.7.11b | ~96–99.8% (tie-adjusted) concordance |
| Determinism | **thread-count-invariant, byte-identical** | seeded `StdRng` (not `mt19937`); not byte-identical on ties |
| STAR params recognised | **200** (locked by a test) | ~86 |
| Feature breadth | **entire STAR surface** (solo, WASP, STARlong, transform, liftOver, wig, chimeric-all) | **core aligner** (no solo/WASP/STARlong/transform/wig) |
| Verification | native-STAR **oracle** + goldens, 474 tests | Python diff vs STAR output + ~434 unit tests |
| Verified this session | ✅ builds; **464/474 pass** (10 = oracle env. artifacts) | ✅ builds; **434/434 pass** |
| `unsafe` | 4 blocks (prefetch, all documented) | **0** |
| `unwrap()` in non-test code | ~26 | a subset of 568 total |
| Provenance | **dedicated crate** (SHA-256 graph) | `Log.out`/`@PG` headers |
| CI / CD | release-artifact workflow only (**no test CI, by design**) | **full CI** (build+test+fmt+clippy+`cargo audit`), docs→Pages, multi-arch Docker + crates.io release, Dependabot |
| Backing / governance | single-author INSERM/IPNP project | **scverse + nf-core** community orgs (repo under `github.com/scverse`) |
| Big-code liability | `align_faithful` (2312 lines), `run` (1126) | `io/sam.rs` (4165), `align/stitch.rs` (3522) — files, not single fns |
| STAR C++ source citations | ~118 `.cpp` refs, 100% module docs, **0 TODOs** | ~72 `.cpp` refs, ~1821 STAR mentions, 6 TODO markers |

### Bottom line

- **On completeness and fidelity, STAR-rs is far ahead.** It is a near-feature-complete, byte-exact, production-grade drop-in for STAR 2.7.11b with an exceptional engineering process. rustar-aligner is an earlier-stage, core-only port that is excellent where it is implemented but covers a fraction of the surface and does not aim for byte-identity.
- **On maintainability the two are more evenly matched, and it depends which lens you use.** STAR-rs is better-factored at the *code* level (clean acyclic crate DAG, minimal `unsafe`, ~28% comment density, zero TODOs, a self-tested SAM normaliser). But its verification is fragile: **correctness is only checkable on a machine with a matching native STAR oracle, and there is no test CI** (proven live by the 10 platform-coupled failures above). rustar-aligner is stronger on *process and governance*: **zero `unsafe`**, a simpler one-crate mental model, and — decisively — a **full automated CI gate** (build+test+fmt+clippy+`cargo audit` on every push/PR) plus docs/release/Docker/crates.io automation and Dependabot, all **backed by the scverse and nf-core communities**. Its weaker points are a large `unwrap()` count and a looser (statistical, RNG-divergent) correctness bar. Net: "cleaner code" is not a clean win for either — STAR-rs leads on documentation and determinism discipline and has the stronger correctness *guarantee*, while rustar-aligner has cleaner error handling (`thiserror` vs mixed `String`/`Box<dyn Error>`), **zero `unsafe`**, and no 2,000-line functions; the multi-crate split is *not* itself a cleanliness advantage (see §6.1). Where rustar-aligner clearly leads is the more robust, portable, community-backed *engineering process* and longevity story.
- **On future extension, the two diverge by design.** STAR-rs's strict "no silent divergence / everything oracle-locked / thread-invariant" contract makes it extremely safe to extend *within STAR 2.7.11b* but actively **hostile to features that go beyond STAR** (any new behaviour must still be byte-justified). rustar-aligner's looser contract, single-crate simplicity, and — importantly — its **scverse/nf-core institutional backing and multi-contributor CI workflow** make it a more natural base for **divergent/novel features** (its stated long-term goal) and for sustained community maintenance, at the cost of weaker byte-level regression guarantees.
- **Choose STAR-rs** if you need a trustworthy, complete, deterministic STAR 2.7.11b replacement today and can supply a matching native oracle to verify it. **Choose rustar-aligner** if you want a smaller, hackable, community-governed core with automated CI to build *new* aligner behaviour on, and can accept partial feature coverage and tie-level divergence from STAR.

The rest of this document details each project (§2–§5), the maintainability/extensibility comparison (§6), a **core-algorithm head-to-head with concrete portable correctness fixes for rustar-aligner (§7)**, and methodology (§8).

---

## 2. STAR-rs — Architecture & Infrastructure

**Repository:** `/home/jamfer/repos/STAR-rs`. Clean-room Rust reimplementation emulating **STAR 2.7.11b**. Workspace version `3.1.0`, edition 2021, dual-licensed `MIT OR Apache-2.0`, org `IPNP-BIPN` (Benjamin Demaille). The README is explicit that it is *not* affiliated with STAR and *not* a "STAR 3" release — `3.x` is its own SemVer line.

### 2.1 Workspace / crate structure (10 crates, acyclic DAG)

| Crate | src LOC | Responsibility |
|---|--:|---|
| `star-core` | 916 | Shared vocabulary: 2-bit DNA, `Strand`, `Interval`, `Cigar`, `AlignmentRecord`, `Score`. Zero deps. No I/O, no algorithm. |
| `star-io` | 1,727 | External-format boundary: FASTA/FASTQ/SAM/BAM (BGZF via noodles + libdeflate). |
| `star-index` | 5,950 | Suffix-array genome index (SA build, `suffix_cmp`), prefix lookup table, sjdb insertion, GTF, supertranscript, transform. Versioned on-disk format (v7). |
| `star-seed` | 957 | Maximal-Mappable-Prefix seeding over the two-strand SA. |
| `star-align` | 12,988 | The algorithmic bulk: windowing, stitching, scoring, PE, chimeric, **STARsolo engine**, WASP, wig, quant, velocyto, emptydrops. |
| `star-filter` | 43 | Contract crate — currently an identity pass-through (real filtering lives in `star-align`). |
| `star-solo` | 27 | Architecture-marker crate only (`is_supported() -> true`); real solo code is in `star-align`/`star-cli`. |
| `star-provenance` | 622 | SHA-256 artifact sidecars + run-dependency graph (Graphviz export). |
| `star-diff` | 283 | Differential-testing harness: locate/run native STAR oracle, normalise SAM. |
| `star-cli` | 9,813 src + 21,683 test | The `STAR`/`STARlong` binaries + orchestration. `commands.rs` alone is 6,384 lines. |

The dependency graph (auto-generated to `docs/architecture/crates.dot` and drift-checked in CI script) is clean and acyclic: `star-core` at the base → `star-io`/`star-index`/`star-seed` → `star-align` → `star-cli` at the top. Two crates (`star-solo`, `star-filter`) are essentially aspirational placeholders whose real logic lives elsewhere — a mild dilution of the "one crate, one job" story, acknowledged in the README.

### 2.2 Dependencies — minimal and heavily justified

Root `Cargo.toml` centralises versions in `[workspace.dependencies]`; **every dependency is documented in `DEPENDENCIES.md` with a why-not-std rationale and an argument that it does not affect output bytes.** 11 runtime deps: `serde`/`serde_json` (versioned index + JSON provenance; BTreeMap → deterministic key order), `sha2` (audited hashing), `clap` v4, `thiserror` v2, `noodles-sam`/`noodles-bam` (BAM codec), `bgzf` 0.4 (libdeflate-backed, chosen over noodles-bgzf for speed with a byte-identity locking test), `smallvec`, `rayon` (parallelism that never affects order), `flate2` **rust_backend only** (pure-Rust gzip, no C zlib), `mimalloc` (global allocator, ~10% speedup, argued output-neutral). The only C-backed pieces (libdeflate, mimalloc) are explicitly defended. `Cargo.lock` is committed as "the pin of record."

### 2.3 Build & reproducibility infrastructure

- **`rust-toolchain.toml` pins channel 1.96.0** + clippy/rustfmt.
- **Release profile: `codegen-units = 1` + `lto = "fat"`** — removes parallel-codegen nondeterminism.
- **`build.rs`** captures git short-commit + rustc version into the provenance-quality `--version`.
- **`scripts/check.sh`** is the authoritative local gate: `cargo fmt --check` + `cargo clippy --workspace --all-targets -- -D warnings` + `cargo test --workspace` + an arch-graph staleness diff.
- **CI is deliberately release-artifact-only** (`.github/workflows/release.yml`; builds Linux+macOS binaries with SHA-256 sidecars on `v*` tags). **There is no test/lint CI, by explicit decision** — the differential tests need a native STAR oracle that CI runners lack.
- **`star-provenance`** is a genuine differentiator: per-output `.prov.json` sidecars record SHA-256 of the inputs that actually fed an output, effective params, tool identity (own semver + `emulated_reference: STAR 2.7.11b` + git + rustc), plus a run-level dependency graph. Determinism is baked in (`RunManifest::canonical()` sorts/dedups before serialising).

### 2.4 Documentation

Unusually rigorous: `README.md`, a **43.8 KB `DIVERGENCES.md`** (every intentional difference D1–D25, each with a locking test), a 27 KB `CHANGELOG.md` (Keep-a-Changelog + SemVer), `ROADMAP.md` (vertical-slice model + candid residual-tail disclosure), `DEPENDENCIES.md`, `CONTRIBUTING.md`, `docs/` (a 38.7 KB window-model port log, benchmark, publication-audit, unsupported-parameters), a JOSS-style `paper/`, `tool.yaml` (bio.tools descriptor), `CITATION.cff`, and 63 per-feature `test-data/` fixture dirs.

---

## 3. STAR-rs — Algorithm & Determinism

The port tracks STAR's C++ at the function — and even dead-code — level, citing specific source files throughout (`SuffixArrayFuns.cpp`, `stitchWindowAligns.cpp`, `stitchAlignToTranscript.cpp`, `Genome_genomeGenerate.cpp`, `sjdbPrepare.cpp`, …).

- **Index (`star-index/src/star_gen.rs`).** Builds STAR's doubled genome `D = G ++ revcomp(G)` and sorts non-spacer offsets with `par_sort_unstable_by(suffix_cmp)`. `suffix_cmp` is a direct port of `funCompareSuffixes` (byte order `0<1<2<3<4<5`, ties broken by smaller position) — a *strict total order*, so the rayon parallel sort is provably identical to a sequential one. The result is **byte-identical to STAR's SA**. sjdb insertion (`sjdb.rs`) reproduces `sjdbBuildIndex`: it does **not** re-sort the whole SA but computes insertion ranks and merges junction suffixes in.
- **Seeding (`star-seed/src/pieces.rs` + `star-index/src/search.rs`).** MMP primitives (`max_mappable_length`, `find_mult_range`, `compare_seq_to_genome`) port `SuffixArrayFuns.cpp`; `map_one_read`/`quality_split` mirror `ReadAlign_mapOneRead.cpp`. Seeds are batched to interleave DRAM-miss chains (~3× seeding speedup, output-neutral).
- **Stitch/score (`star-align/src/star_model.rs`, `star_window.rs`).** `stitch_window_aligns` is a direct recursive include/exclude port; `stitch_align_to_transcript` is the scoring heart (SJ motif classification, junction-shift scan, micro-repeat length, sjdb lookup, gap/intron/deletion scoring). `extend_align` is STAR's greedy left/right extension (not a Smith-Waterman substitute — this is what makes soft-clip boundaries byte-identical). The genomic-length-log2 penalty matches STAR.
- **Multimapping / primary selection — the RNG question.** *By default there is no RNG.* STAR's default `Old_2.4` order + `OneBestScore` primary flag do not consume `mt19937`, so STAR-rs reproduces STAR by reproducing its **SA/window discovery order** plus a fixed tie-break (max score → smaller genomic length → earliest). The one RNG path, `--outMultimapperOrder Random`, is **intentionally divergent (DIVERGENCES D16)**: native STAR seeds per-thread (so its Random output changes with thread count), whereas STAR-rs seeds **per read** via a `splitmix64` Fisher-Yates shuffle — byte-identical across every `--runThreadN`. A faithful `mt19937` *does* exist but only for STARsolo EmptyDrops_CR Monte-Carlo.
- **Determinism mechanisms.** Pure per-read alignment collected in **input order** (`par_chunks().map().collect().concat()`); explicit total-key sorts for every observable ordering; deterministic parallel SA sort; per-read RNG seeds; single-codegen-unit/fat-LTO build; pure-Rust gzip.
- **`unsafe`:** only 4 sites, all in `star-index`, all `// SAFETY:`-documented prefetch-hint/bounds-checked pointer ops. **No mmap** — the index is read fully into RAM via `std::fs::read` (matches STAR's in-memory model, keeps loading deterministic).

---

## 4. STAR-rs — Features, Testing & Maturity

### 4.1 Feature inventory (verified against source — all substantial, non-stub)

- **Core:** SE/PE, splice, sjdb, indels, multimapping, 2-pass, PE mate-overlap merge (`--peOverlapNbasesMin`). *Complete.*
- **Output:** SAM + BAM (Unsorted + SortedByCoordinate), full `--outSAMattributes` (incl. derived NM/MD/jM/jI), `--outWigType` signal tracks, `SJ.out.tab`, `Log.final.out`, BySJout, unmapped (SAM `Within` + FASTX). *Complete.*
- **Quant:** GeneCounts, TranscriptomeSAM (SE+PE), Transcript3p. *Complete.*
- **Chimeric:** `Chimeric.out.junction` (both strands), WithinBAM, SeparateSAMold, PE, `--chimMultimapNmax`, PE-merged. *Near-complete* (one documented output-neutral tag limitation).
- **STARsolo (the biggest surface):** `CB_UMI_Simple`/`Complex`/`SmartSeq`; raw+filtered matrices; all six UMI-dedup modes; all five implemented `--soloCBmatchWLtype` modes (incl. CellRanger's pseudocount posterior); features Gene/GeneFull/SJ/Velocyto/Transcript3p; EmptyDrops_CR cell calling; Summary.csv/*.stats; CB/UB BAM tags; **10x chemistry auto-detection** with 7 bundled whitelists (~26 MB). *Broad and deep.*
- **WASP** allele-specific filtering (vW/vA/vG).
- **Run modes:** genomeGenerate, liftOver, soloCellFiltering, `inputAlignmentsFromBAM` (coverage + bamRemoveDuplicates), genome transform (haploid+diploid with back-transform), SuperTranscriptome, on-the-fly insertSequences.
- **STARlong** long-read chaining-DP binary (incl. chimeric).

### 4.2 Parameter coverage — 200, machine-locked

A test (`recognizes_every_star_2711b_parameter`) reads a committed 200-name fixture and asserts clap recognises every one. Unimplemented-but-recognised knobs follow a **loud-reject, accept-default** discipline: accepted at STAR's default (byte-identical by construction) but **fatally reject any non-default value** with a STAR-style error. The residual rejected set is small and enumerated (`docs/unsupported-parameters.md`): a few `--genomeLoad` shared-mem modes, `EditDist_2`, `CB_samTagOut`, some wig modes, on-the-fly sjdb at align, etc. **No silent no-ops.**

### 4.3 Testing — differential against a live oracle

- **474 test functions** (231 unit incl. `proptest` round-trips on DNA/CIGAR; 243 differential across **104 `differential_*.rs` files**, ~1:1 with feature areas).
- **Method:** run native STAR 2.7.11b + STAR-rs on identical inputs, compare **normalised SAM** (drop volatile `@PG`/`@CO`, sort tags, sort records). The normaliser is itself unit-tested. Oracle located via `$STAR_BIN` → Homebrew path → `PATH`; **skips gracefully if absent.** Most tests also carry a committed **golden**, so both "golden == live STAR" and "STAR-rs == golden" are locked.
- **Determinism tests** run at 1/4/16 threads asserting byte-identity.
- **Real-data validation:** 100.0000% byte-identical placement on GIAB ERR356372 (1M PE, full GRCh38, 0/2,283,844 differing).

### 4.4 Divergences & maturity

- **`DIVERGENCES.md`** enumerates D1–D25, each locked by a test. Several are **STAR *bugs* deliberately fixed** (D20 SuperTranscriptome minus-strand RC; D22 `MultiGeneUMI_All` no-op) under a "fix reproducibly, never silently" policy. D1: the reference oracle is the `2.7.11b_alpha_2024-02-09` build.
- **Version 3.1.0** (2026-07-22), rapid 3.0.0→3.1.0 cadence, exemplary CHANGELOG.
- **Residual known tail (documented, no golden fails):** a <0.05% low-complexity poly-A/T enumeration tail; a ~0.003% microexon+indel-at-annotated-junction case; one output-neutral single-thread stitch speed gap. Performance: ~1.4–2.2× faster than native STAR at 8–16 threads, within ~4% single-threaded.

---

## 5. rustar-aligner — Review (same axes)

**Repository:** `/home/jamfer/repos/rustar-aligner`. A faithful-*port* of STAR, MIT-licensed, v0.1.0, Rust 2024. Single crate + one binary (`rustar-aligner`).

### 5.1 Architecture

Clean domain modules mapping to STAR concepts: `genome/`, `index/` (packed_array, suffix_array, sa_index, io), `align/` (seed, stitch, score, transcript, read_align), `io/` (fastq, sam, bam, log), `junction/` (gtf, sj_output, sjdb_insert), `quant/` (mod, transcriptome), `chimeric/` (detect, segment, score, output), plus `params/`, `error.rs`, `mapq.rs`, `stats.rs`, `cpu.rs`. **~30,468 LOC.** Largest files: `io/sam.rs` (4,165), `align/stitch.rs` (3,522), `quant/transcriptome.rs` (2,454), `align/read_align.rs` (2,177), `params/mod.rs` (1,596) — the top two are large enough to be a navigability concern. Mainstream deps (clap, anyhow, thiserror, noodles 0.109, memmap2, rayon, dashmap, rand 0.10, bitflags, shlex); tuned release profile (`lto="fat"`, `codegen-units=1`, `strip`).

### 5.2 Core algorithm

- **Suffix array** (`index/suffix_array.rs`): custom `compare_suffixes` with a `packed_a.cmp(&packed_b)` tie-break that makes the SA **byte-for-byte identical to STAR's** for the yeast genome (Phase G3). Variable-width `PackedArray` backing store; 35-bit k-mer prefix table. Can load a STAR-generated index.
- **Seeding** (`align/seed.rs`, 1,051 LOC): hierarchical SAindex lookup + MMP binary search, STAR-style.
- **Stitch/score** (`align/stitch.rs` 3,522; `align/score.rs` 1,413): seed clustering into windows + recursive DP stitching + `extendAlign`, tracking STAR functions by name in comments.
- **Multimapping / tie-break — the key divergence.** rustar-aligner uses `rand::StdRng` seeded **per read** via `per_read_seed(run_rng_seed, read_name)`, shuffling only the equal-top-score prefix. This **deliberately does not reproduce STAR's `mt19937`** primary choice. Consequence: near-exact on uniquely-mapping reads, but ties pick a different repeat copy → 299 (SE) / 475 (PE) tie-break diffs, excluded from the "tie-adjusted" metric.
- **PE:** combined-read seeding, joint DP, per-mate seeding, half-mapped fallback, split-combined-WT — reports PE both-mapped 8390 (exact match to STAR), 0 half-mapped, 0 NH/MAPQ diffs.

### 5.3 Feature breadth

- **Implemented:** SE+PE; splice; sjdb (with in-index insertion at genomeGenerate); indels; multimapping; two-pass; `GeneCounts`; `TranscriptomeSAM`; chimeric **SE+PE** (4 tiers) incl. `--chimOutType WithinBAM`; BAM unsorted + coordinate-sorted; `--outStd`; `--outReadsUnmapped Fastx`; BySJout; GTF-tag params; `--outBAMcompression` + `--limitBAMsortRAM`; `--outSAMattrRGline`; `--runRNGseed`; Log.final/out/progress; MAPQ lookup table.
- **NOT implemented (verified absent):** STARsolo/single-cell (explicitly **deferred**), WASP, STARlong, genome transform / liftOver (only static header placeholders), signal/wig tracks. This is the single biggest gap vs STAR-rs.

### 5.4 Parameter coverage

**~86** `--camelCase` params (the figure in the retired CLAUDE.md, ~52, was stale). clap's derive parser errors on unknown flags; some recognised-but-unenforced knobs are accepted with a warning (e.g. `--limitGenomeGenerateRAM`). No loud-reject discipline as strict as STAR-rs's, but it does not silently swallow arbitrary flags.

### 5.5 Testing

- **~434 `#[test]`** in `src/` (the retired CLAUDE.md cited 396); **14 integration tests** in `tests/` (`alignment_features.rs`, `phase9_threading.rs`, `transcriptome_sam.rs`) via `assert_cmd`.
- **Differential testing via Python** (`test/`): `compare_sam.py`, `compare_sam_thorough.py`, `compare_pe.py`, `assess_faithfulness.py` (exact FLAG/RNAME/POS/MAPQ/CIGAR/NH/AS/NM + PE fields + SJ.out.tab), `compare_junctions.py`, `compare_chimeric.py`, etc. — these consume **real STAR output as ground truth**.
- **Target is statistical, not byte-identical:** ~96–99.8% on 10k yeast reads with a **"tie-adjusted"** metric that excludes RNG tie-breaks. **0 STAR-only / 0 rustar-only reads** on SE — every disagreement is a documented tie.

### 5.6 Code quality & maturity

- **Error handling:** `thiserror` enum + `anyhow`. Idiomatic.
- **`unsafe`: 0** (notable given memmap2/byteorder) — mmap access is wrapped safely.
- **`.unwrap()`: 568** (concentrated in tests/lock/parse sites, but a portion in non-test paths — the main quality caveat). `.expect()`: 8. TODO/FIXME/panic!/unimplemented!: 6 total — very low.
- **STAR C++ traceability:** ~1,821 STAR/Dobin mentions, **72 explicit `.cpp` citations** with function/line — a standout maintainability strength.
- **Naming:** STAR camelCase via `#![allow(non_snake_case)]`.
- **Maturity:** v0.1.0, ~304 commits, active numbered-PR history. Known issues are honest and narrow: the 299/475 tie divergences, 1 CIGAR-only insertion-placement tie, and **4 PE AS diffs that are actually rustar-aligner *improvements*** over STAR (finds better-scoring alignments STAR's combined-window approach misses).

---

## 6. Maintainability & Future Extension

### 6.1 Modularity & navigation — and why "multi-crate ≠ cleaner"

An important correction to a framing that is easy to get wrong: **STAR-rs's 10-crate workspace is not, in itself, "cleaner" than rustar-aligner's single crate.** The distinction only matters if the crates are *self-contained, reusable libraries* that add value to the ecosystem. STAR-rs's are not — they are internal pipeline stages, and two are effectively empty markers (`star-solo` is 27 lines of `is_supported() -> true`; `star-filter` is 43 lines of identity pass-through). The same encapsulation is achievable with Rust's `mod` system inside one crate, without the version-pinning, contract, and build ceremony a workspace imposes. For a single-purpose binary that ships no libraries, the crate split is at best organisational preference and at worst clutter; publishing internal-only crates would also namespace-squat crates.io (STAR-rs avoids this by publishing only binaries, not the crates). rustar-aligner's single-crate layout, modularised via `mod`, is a legitimate and arguably tidier choice for a tool.

- **What genuinely helps navigation in STAR-rs** is not the crate count but the *discipline*: an auto-generated + drift-checked dependency/mod graph, a near-1:1 test-to-feature mapping, and module-level docs on 100% of files. A contributor can orient quickly. The real complexity still concentrates in one place (`star-align`, 13k LOC).
- **rustar-aligner** has one mental model and no cross-crate overhead. Its navigation cost is local: two oversized files (`io/sam.rs` 4.2k, `align/stitch.rs` 3.5k) would benefit from being split into submodules.

Net: this axis is roughly a wash. Neither the crate split nor the single crate is inherently cleaner; both projects modularise reasonably.

### 6.2 The shared "giant function/file" problem

Both carry heavy hotspots. STAR-rs's are **single functions**: `align_faithful` (2,312 lines, ~13 levels of nesting) and `run` (1,126). rustar-aligner's are **files** (`sam.rs`, `stitch.rs`). STAR-rs's ~500-line stitch functions are defensible (1:1 ports of STAR's monolithic C++ methods, where structural alignment with the oracle has real value); the CLI god-functions are not. This is the top maintainability flag for *both* projects.

### 6.3 Safety & correctness hygiene

- **`unsafe`:** rustar-aligner **0**; STAR-rs 4 (all documented, output-neutral). Both excellent.
- **Panic surface:** STAR-rs has only ~26 non-test `unwrap()` in 33k lines — very low. rustar-aligner's 568 total is higher and a genuine caveat, though heavily test-concentrated.
- **Error types:** rustar-aligner is more consistent (`thiserror` enum + `anyhow` throughout). STAR-rs is looser (`Result<_, String>` / `Box<dyn Error>` mixed) — its one idiom weak spot.
- **Comments/traceability:** both cite STAR C++ heavily. STAR-rs has ~28% comment density, 100% module docs, **0 TODO markers**; rustar-aligner has 72 `.cpp` cites and only 6 TODO markers. Both are well above scientific-software norms.

### 6.4 Verification model — the decisive maintainability difference

- **STAR-rs** has the stronger *guarantee* (byte-identity, oracle-locked, thread-invariant) but the more *fragile* verification: **there is no test CI, and the tests only pass on a machine with a matching native STAR oracle.** This session proved it live — 10 tests "failed" here purely because the local oracle is Linux/glibc (`-nan%`) while the goldens were locked on macOS (`nan%`), plus one native-STAR crash. A contributor without the right oracle cannot verify a change, and a careless PR can merge without the gate ever running.
- **rustar-aligner** has a weaker *guarantee* (statistical, RNG-divergent) but more *portable and automated* verification: its in-tree unit + integration tests run anywhere with `cargo test` (no oracle required), and it enforces a **real CI gate on every push and PR** — `.github/workflows/ci.yml` runs `cargo build --release && cargo test --release`, `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`, `cargo check`, and `rustsec/audit-check`, behind an `alls-green` status check. It further ships `docs.yml` (Astro site → GitHub Pages), a full `release.yml` (multi-target binaries + multi-arch Docker buildx + `cargo publish` to crates.io), and Dependabot. The Python differential scripts against real STAR output are an optional deeper check on top. **This is the more robust maintainability posture of the two:** regressions are caught automatically by CI that any contributor's PR must pass, whereas STAR-rs's stronger byte-identity guarantee is only ever checked manually on a machine with a matching native oracle.

### 6.5 Building forward

- **Extending *within* STAR 2.7.11b:** STAR-rs is safer — the "no silent divergence + oracle-lock + thread-invariance" contract catches regressions precisely, and the crate boundaries localise change. The cost is high ceremony: every new knob needs a golden, a differential test, provenance, and docs (its documented vertical-slice model).
- **Extending *beyond* STAR (novel features, new science):** rustar-aligner is the more natural base. Its looser fidelity bar, single-crate simplicity, zero `unsafe`, and stated long-term goal ("extra features and divergences from the original will come in later releases/forks") make experimentation cheaper. STAR-rs's whole architecture is *designed to resist* divergence, so novel behaviour fights the grain.
- **Onboarding:** STAR-rs's periphery is easy; its alignment core assumes real STAR-C++ familiarity (mitigated by dense comments). rustar-aligner's single-crate layout is easier to hold in one's head end-to-end, but its two mega-files and higher `unwrap()` density raise the local reading cost.

---

## 7. Algorithm Head-to-Head & Portable Correctness Improvements

A focused cross-repo diff of the core algorithms (seeding, windowing, stitching, scoring, multimapper selection), cross-checked against the real STAR 2.7.11b C++ source. **The good news: the two implementations agree on almost everything that matters numerically.** The gaps are concentrated in a few well-defined places, and every one of them is a spot where **STAR-rs is closer to STAR and its approach can be ported into rustar-aligner** to improve faithfulness — fully consistent with rustar-aligner's own "match STAR exactly, never revert a STAR-faithful change" development rule.

### 7.1 Where the two are already identical (no work needed)

- **Scoring constants — all of them:** sjdbScore=2, scoreStitchSJshift=1, scoreGap/Noncan/GCAG/ATAC = 0/-8/-4/-8, del/ins open+base = -2/-2/-2/-2, alignIntronMin/Max = 21/589824, alignSJstitchMismatchNmax=[0,-1,0,0], scoreGenomicLengthLog2scale=-0.25, match=+1.
- **Splice-motif truth tables** map one-to-one (GT/AG, CT/AC, GC/AG, CT/GC, AT/AC, GT/AT, non-canonical).
- **Genomic-length penalty formula** `ceil(log2(gExtent)·scale − 0.5)`, clamped ≥0.
- **`extendAlign`** (greedy local extension) for the default `Local` end type: budget formulas, match-only max recording, mismatch-before-increment, N-skip, spacer stop — line-for-line equivalent.
- **Junction-shift scan, micro-repeat length, deletion/intron block** — faithfully reproduced.
- **Two-strand suffix array** — rustar-aligner's SA is byte-identical to STAR's for yeast.

### 7.2 Where STAR-rs is more STAR-faithful — prioritized port targets

| # | Gap in rustar-aligner | STAR / STAR-rs behaviour | Impact | Effort | rustar files to change |
|---|---|---|---|---|---|
| **1** | **Multimapper primary uses an unconditional `StdRng` shuffle** of the tied top-score prefix | STAR's **default** primary selection consumes **no RNG** (verified in `ReadAlign_multMapSelect.cpp:69-103`). RNG runs *only* under `--outMultimapperOrder Random`. | **HIGH** — dominant cause of the 299 SE / 475 PE tie diffs | Low | `params/mod.rs` (add `--outMultimapperOrder`, default `Old_2.4`); gate the shuffle at `read_align.rs:379-383` & `1094-1098` |
| **2** | **Wrong primary tie-break key** (`score → n_junction → chr_idx → genome_start → is_reverse`) | STAR's `trBest` key is **`score → smaller gLength → earliest-discovered window`** (`ReadAlign_stitchPieces.cpp:358`). `n_junction`/`chr`/`pos` are not in STAR's key. | **HIGH** — latent second cause of tie diffs even after #1 | Low-Med | `read_align.rs:369-376` (SE) & `1077-1091` (PE); add `g_length` to `transcript.rs`; select by discovery-order min-scan, don't coordinate-sort. Copy `STAR-rs star_window.rs:827-909` |
| **3** | **Forces the last anchor's inclusion** (`can_exclude` gating on `last_anchor_idx`) | STAR 2.7.11b's `WlastAnchor` is **dead code** (init `UINT64_MAX`, guard never true) — it **never** forces a seed, always explores the exclude branch. | **HIGH** — flips transcript sets / primaries; plausibly part of the tie tail | **Very low** — delete the gating | `stitch.rs:2333-2341` (make `can_exclude` always true) + drop `2958` |
| **4** | **sjdb annotated junctions: exact-coordinate boolean lookup only** | STAR-rs has (a) the `sjA` read-adjacent **simple-stitch fast path** and (b) **annotated-boundary snapping** that shifts a non-canonical annotated junction to its annotated coords and applies `sjdbScore`. | **HIGH for `--sjdbGTFfile` runs** (invisible on unannotated yeast) — otherwise mis-places junction boundaries in repeats & forfeits sjdbScore | Med-High | `stitch.rs:1367-1381`; add `sj_a` to `WindowAlignment`; port `STAR-rs star_model.rs:384-413, 602-628` |
| **5** | **Seed array not `rStart`-sorted; dedup keeps `search_rc`** (retains duplicate fwd+rev seeds) | STAR keeps `PC[]` sorted by `rStart` and drops exact `(rStart,length)` duplicates **regardless of direction** (`store_aligns`). This makes window-creation order follow `rStart`. | **MED** — window order feeds the "earliest-window" tie-break (#2) and the CIGAR case | Med | `seed.rs::find_seeds` (92-100); mirror `STAR-rs pieces.rs:184-214` |
| **6** | **No `flagDirMap` reverse-walk suppression** | STAR skips the reverse istart=0 walk when the forward walk already mapped the piece to full length. | **MED** — removes spurious reverse seeds that perturb chains | Low | `seed.rs::find_seeds` / `search_direction_sparse`; mirror `pieces.rs:440,588` |
| **7** | **Genomic-length penalty applied only in a separate post-pass** | STAR-rs folds it into the in-recursion dedup/eviction ranking. | **MED** — can keep/drop the wrong transcript among spliced-vs-unspliced ties | Med | `stitch.rs:2208-2254` |
| **8** | **No `alignEndsType` / `extend_to_end` support** | STAR-rs threads `ext[mate][end]` and has a force-to-end branch. | MED (feature gap) — `EndToEnd`/`Extend5p`/`Extend3p` unsupported | Med | `stitch.rs::extend_alignment` (+ param plumbing); port `star_model.rs:173-212` |
| **9** | Minor: `same_structure` dedup guard; stitch mismatch cap uses `p_mm_max·len` not `outFilterMismatchNmaxTotal`; `alignInsertionFlush Right` unimplemented; `seedPerReadNmax` silently truncates (STAR aborts); chain-stop threshold conflated with store-min-length | Various STAR behaviours | LOW each | Low | `stitch.rs:1530, 2224`; `seed.rs:75-77, 260, 414` |

### 7.3 The two known rustar-aligner disagreements — likely fixes

- **The 299 SE / 475 PE tie-break diffs** (rustar-aligner's headline "single biggest correctness gap"): caused by **#1 + #2**. Removing the default shuffle (#1) and adopting STAR's `trBest` key with preserved discovery order (#2) should close the large majority; **#3** and **#5** address the residual (they restore STAR's exact "earliest-window" ordering). Notably, this makes rustar-aligner *more* STAR-faithful than its current `--runRNGseed`/`StdRng` design — STAR simply does not use RNG here by default, so the fix is "stop shuffling," not "reimplement mt19937."
- **The 1 CIGAR-only disagreement (`ERR12389696.13573895`, insertion at read pos 100 vs STAR's 108):** the causal 71-base seed lands 8 bases off because of a different seed-chain path through a homopolymer. The fix levers are **#5 + #6 + #7-adjacent window-assignment order** (process seeds in `rStart` order rather than length-descending pre-dedup) plus possibly **#3** — i.e. reproduce STAR's exact `Lmapped` chain and `rStart`-ordered, direction-agnostic seed handling.

### 7.4 Recommended sequencing

Do the cheap high-impact ones first and re-benchmark after each with the **raw exact-match** metric (see §7.5 for why the tie-adjusted % misleads):

1. **#3 (delete WlastAnchor forcing)** — one-line-ish, isolates cleanly. ✓ done (neutral on yeast, correct).
2. **#1 (RNG-free default primary)** — gate the shuffle behind `--outMultimapperOrder Random`. ✓ done (+480 raw exact, no regression). **This was the real win — see §7.5.**
3. **#5 (seed `rStart` ordering + direction-agnostic dedup)** — STAR's `storeAligns` semantics. ✓ done (+8 raw exact, no regression). **#6 (flagDirMap) turned out redundant with #5's direction-agnostic dedup.**
4. **#2 (STAR's `gLength` tie-break key)** — ~~apply after #5~~ **DROPPED**: measured to regress both with and without #5 (§7.5). The coordinate key matches STAR better for these ties.
5. **#4 (sjdb snapping/fast-path)** — before any serious `--sjdbGTFfile` benchmarking. Best remaining lever.
6. **#7, #8, #9** as follow-ups.

All line numbers above are current as of this review; treat them as starting points, not guarantees.

### 7.5 Empirical results — measured on the yeast 10k benchmark (this session)

Fixes were implemented and measured against real STAR 2.7.11b. **Use the raw exact-match count, not the "tie-adjusted" %** — the tie-adjusted denominator shifts as reads move in/out of the tie bucket, which masks real changes.

**Fix #3 (remove `WlastAnchor` forcing)** — merged on `fix/wlastanchor-no-force`. Behaviorally neutral on yeast 10k (no read's output changed); correct and STAR-faithful, zero regression. Only alters output on datasets where the last anchor is reached with no anchor yet included and excluding it wins.

**Fix #1 (gate the multimapper RNG shuffle behind `--outMultimapperOrder Random`; default `Old_2.4` is deterministic)** — the big win:

| Config | SE raw exact / 8926 | PE raw exact / 16778 |
|---|---|---|
| Baseline (`origin/main`, shuffle always on) | 8606 (96.4%) | 16278 (97.0%) |
| **#1 — no default shuffle** | **8789 (98.5%)** | **16575 (98.8%)** |
| #1 + #2 (add STAR's `gLength` key) | 8633 (96.7%) | 16261 (96.9%) |

- **#1 alone gains +183 SE and +297 PE exact matches (+480 total), zero regression.** rustar-aligner was shuffling the tied top-score prefix unconditionally, randomising the primary among equal-score loci; STAR's default consumes no RNG. Removing the default shuffle makes the primary deterministic and STAR-matching. Residual tie-break diffs collapse (SE ~299→127, PE ~489→177).

**Fix #5 (STAR-faithful seed dedup + `rStart` ordering)** — measured on top of #1 (`fix/seed-rstart-order`, based on the #1 branch):

| Config | SE raw exact / 8926 | PE raw exact / 16780 |
|---|---|---|
| #1 only | 8789 | 16575 |
| **#5 + #1** | **8790** | **16582** |
| #5 + #1 + #2 (gLength key) | 8653 | 16288 |

- **#5 is a small positive: +1 SE, +7 PE exact matches (+2 more PE both-mapped pairs), zero regression.** It replaces rustar-aligner's direction-keyed seed dedup with STAR's `storeAligns` semantics — sort `PC[]` by `rStart` (longer first) and drop `(rStart, Length)` duplicates regardless of search direction — so window-creation order tracks STAR's. The direct gain is small because #1 already made the primary deterministic and rustar's discovery order was already close; the value is correctness (matches STAR's dedup) and it is a clean prerequisite step.
- **Critical negative result: #5 does NOT unlock #2.** Adding STAR's `gLength` tie-break key *on top of* #5 still regresses hard (SE 8790→8653, PE 16582→16288). This falsifies the earlier hypothesis. The reason: most multimapper ties are **same-`gLength` repeat copies**, so `gLength` doesn't discriminate and the decision falls to the *earliest-window* fallback — and even #5-corrected discovery order tracks STAR's chosen copy *worse* than rustar-aligner's existing explicit coordinate key (`chr, pos`). **Empirically, STAR's documented `gLength` key is the wrong move for rustar-aligner here; the coordinate key matches STAR's actual output better.** #2 should be dropped from the plan (or revisited only with a much more careful look at STAR's exact combined-window `gLength` and window order — not worth it at current margins).

**Fix #4 (sjdb annotated junctions)** — investigation on the annotated (`--sjdbGTFfile`) benchmark surfaced a more fundamental bug than the "boundary snapping" originally scoped. Using STAR-rs as the reference for the snap mechanism (it uses a post-scan `sjdb.find` + precomputed `shift_left/shift_right`) led to checking annotation recognition end-to-end, which revealed:

- **The runtime `SpliceJunctionDb` was populated *only* from `params.sjdb_gtf_file`** (`index/io.rs`). In the standard workflow — build the index once with `--sjdbGTFfile`, then align with only `--genomeDir` — no GTF is passed at align time, so the db was **empty** and *every junction was treated as novel*. Measured: rustar recognised **0** annotated junctions vs STAR's **71** on the yeast GTF run.
- **Fix (this branch, `fix/sjdb-annotated-snapping`):** when no GTF is passed but the index carries `sjdbInfo.txt` junctions, populate the runtime db from them (keyed on the stored post-`sjdbPrepare` donor/acceptor coords). Result:

| Annotated-run metric | Before | **After** | STAR |
|---|---|---|---|
| Annotated junctions recognised | 0 | **70** | 71 |
| Total junctions | 310 | **317** | 322 |
| SE raw exact / 8926 | 8779 | **8786** | — |
| CIGAR (same-pos) diffs | 15 | 14 | — |

- Unannotated path unchanged (8790, no regression — indices without `sjdbInfo.txt` skip the fallback). The residual gap (1 annotated junction, −5 total, 14 CIGAR diffs) is what the *actual* boundary-snapping (STAR-rs `star_model.rs:602-628`) would address — a good follow-up now that recognition works.

Net revised sequencing (corrected by measurement): **#3 ✓ → #1 ✓ (+480) → #5 ✓ (+8) → #4a ✓ (annotated-db load; 0→70 annotated junctions) → #4b (boundary snapping, follow-up)**. **#2 (gLength key) is dropped** — it regresses both with and without #5. The recurring lesson: measure with raw exact-match against the real oracle; STAR's documented logic does not always translate to better agreement when another part of rustar-aligner's pipeline diverges, and end-to-end benchmark checks surface root causes (empty align-time db) that source-reading alone missed.

---

## 8. Methodology & Verification Notes

- **Reviewed by fan-out:** five specialist agents covered STAR-rs architecture, STAR-rs algorithm/determinism, STAR-rs features/testing, STAR-rs code quality, and a parallel rustar-aligner characterisation. All findings above are drawn from source reading; file paths and line numbers were cited by the agents and spot-checked.
- **STAR-rs was built and smoke-tested on this machine (Linux/glibc):**
  - `cargo build --release` → **clean, exit 0, ~1m05s**, all 10 crates + `STAR`/`STARlong` binaries. `STAR --version` → `STAR 3.1.0 / emulated_reference: STAR 2.7.11b / index_format_version: 7`.
  - `cargo test --release --workspace --no-fail-fast` → **464/474 passed, 10 failed.** All 10 failures are **environmental oracle artifacts, not defects**: 9 are the glibc `-nan%` vs Apple-libc `nan%` `printf` NaN-formatting difference in `Log.final.out` (STAR-rs's goldens were locked against a macOS STAR oracle; native STAR here is Linux/glibc), and 1 is the **native STAR oracle itself crashing** on `--alignWindowsPerReadNmax 1`. This directly corroborates the project's own D1 caveat and the §6.4 verification-portability weakness.
  - The native oracle present on this host is genuine **STAR 2.7.11b** (`/home/jamfer/Dropbox/.../STAR/STAR`), the exact version STAR-rs targets, so the differential tests ran meaningfully rather than skipping.
- **rustar-aligner was also built and tested on this machine:** `cargo test --release` → **clean, exit 0, 434/434 tests passed, 0 failed.** Its in-tree suite runs to completion with no external oracle required — corroborating the §6.4 point that its verification is more portable, even though its correctness bar (statistical, tie-adjusted) is weaker than STAR-rs's byte-identity.
