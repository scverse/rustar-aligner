# Divergences from STAR

rustar-aligner is a faithful port of [STAR](https://github.com/alexdobin/STAR) 2.7.11b: the goal is to match STAR's algorithms, thresholds, and output byte-for-byte wherever it is reasonable to do so. This file is the complete, authoritative list of the places where the two **do** differ, and why.

Every entry here is a **deliberate, signed-off** decision or a **known, tracked** residual difference — not an accident. Per [CONTRIBUTING.md](CONTRIBUTING.md), a change that diverges from STAR must be recorded here, must never be presented as "faithful", and must not invent a STAR flag or behaviour that does not exist. If you find behaviour that differs from STAR and is *not* listed here, that is a bug — please open an issue.

Divergences are grouped by kind:

1. [Deliberate algorithmic divergences (affect alignment output)](#1-deliberate-algorithmic-divergences)
2. [Cases where rustar-aligner produces a better alignment than STAR](#2-cases-where-rustar-aligner-outperforms-star)
3. [Output-file metadata divergences (not alignments)](#3-output-file-metadata-divergences)
4. [Implementation divergences with no intended output difference](#4-implementation-divergences-no-intended-output-difference)
5. [Known residual single-read differences (tracked, not chosen)](#5-known-residual-single-read-differences)

---

## 1. Deliberate algorithmic divergences

### 1.1 Multimapper tie-breaking / RNG

**What STAR does.** STAR seeds a single `std::mt19937` per read-chunk/thread (`runRNGseed * (iChunk + 1)`) and advances that state sequentially as it processes reads. Under `--outMultimapperOrder Random`, the primary among equal-scoring loci is chosen from that per-thread stream, so the result depends on how reads are partitioned across threads.

**What rustar-aligner does.** rustar-aligner parallelises **per read** via rayon, so a per-thread sequential RNG would make output depend on thread scheduling. Instead it derives a deterministic per-read seed by folding the read name into `--runRNGseed` (`per_read_seed` in `src/align/read_align.rs`), using an in-tree splitmix64 generator (`src/rng.rs`) rather than mt19937.

**Why.** Determinism and thread-count invariance: the same read produces the same primary regardless of `--runThreadN`. STAR's exact mt19937 stream cannot be reproduced under per-read parallelism, and matching it would forfeit reproducibility.

**Impact.** With the default `--outMultimapperOrder Old_2.4`, **no RNG is consulted at all** — the primary is the deterministic best alignment (max score → smaller genomic length → earliest discovered), which is STAR-faithful. The divergence is observable only under `--outMultimapperOrder Random`, and only in *which* equal-scoring locus is marked primary — never in the set of reported alignments.

This is the reason faithfulness is reported **tie-adjusted**. On the 10k yeast benchmark, 299 SE and 475 PE primary-selection differences are all genuine ties: both tools find the identical alignment set, and differ only in which equal-scoring member is primary (from SA-iteration order or the RNG-seed difference above). Excluding those ties, SE is 99.815% and PE 99.883% exact.

**Source.** `src/rng.rs`, `src/align/read_align.rs` (`per_read_seed`, `shuffle_tied_prefix`), `src/params/mod.rs` (`MultimapperOrder`). STAR: `ReadAlign_multMapSelect.cpp`, `ReadAlignChunk` RNG seeding.

---

## 2. Cases where rustar-aligner outperforms STAR

These are not chosen divergences and not bugs: rustar-aligner reports a **higher-scoring, correct** alignment that STAR misses. They are listed here so the differential benchmark's non-exact reads are fully accounted for.

### 2.1 Four PE alignments scored better than STAR

On the 10k yeast PE benchmark, 4 reads differ in alignment score (AS) because STAR's combined-window stitching fails to place the pair at the better location:

- `ERR12389696.844151` — rustar-aligner finds VIII:451791 with 0 mismatches; STAR reports VII:1001391 with 6 mismatches.
- `ERR12389696.4972950` — rustar-aligner finds the correct **spliced** mate 2; STAR reports it unspliced.

**Impact.** rustar-aligner's result is the better alignment in each case. These are counted against exact faithfulness in the raw metric but are improvements, not regressions.

**Source.** See `CLAUDE.md` (PE status) and `STAR-RS-COMPARISON.md`.

---

## 3. Output-file metadata divergences

### 3.1 `genomeParameters.txt` command-line header

**What STAR does.** At `genomeGenerate`, STAR writes a `### <commandLineFull>` header line reproducing the full command line that built the index.

**What rustar-aligner does.** rustar-aligner emits a fixed skeleton containing the parameters it knows at invocation (`--runMode`, `--runThreadN`, `--genomeDir`, `--genomeFastaFiles`, `--genomeSAindexNbases`, `--sjdbGTFfile`, `--sjdbOverhang`). The remaining value lines of `genomeParameters.txt` match STAR's `genomeParametersWrite.cpp` order and tab/space formatting.

**Why.** The header is informational; reproducing an arbitrary STAR invocation's exact argv byte-for-byte serves no functional purpose and the index loads identically either way.

**Impact.** The `###` header line will not byte-match an arbitrary STAR run. No effect on alignment, index loading, or any downstream tool.

**Source.** `src/genome/mod.rs` (`genomeParameters.txt` writer).

---

### 1.3 CellRanger behaviour is the default on 10x geometry

**What STAR does.** STARsolo's defaults are its own (`1MM_multi`,
`1MM_All`, no UMI filtering, `Hamming` clipping, `outFilterScoreMin 0`)
whatever the barcode geometry. Matching CellRanger requires passing five flags,
listed in STAR's `docs/STARsolo.md`.

**What rustar-aligner does.** When the run is unambiguously 10x —
`CB_UMI_Simple`, a whitelist, a 16-base CB and a 10- or 12-base UMI — those
five flags default to their CellRanger values. Any flag given on the command
line wins, and the substitution is logged in full.

**Why.** A user aligning 10x data and comparing against CellRanger otherwise
gets a successful run and different numbers, with nothing pointing at the five
flags that explain it. Measured against CellRanger 10.0.0 on a 20 000-read
fixture, those flags move the count matrix from 8.96% above CellRanger to
0.03% above it, once #165's `cbMinP` posterior threshold is also applied.
STAR 2.7.11b with the same flags is at +0.09%, so all three agree to within a
fraction of a percent.

**Impact.** This is a **change of default output behaviour** and therefore the
largest divergence in this file. It is confined to a geometry nothing else in
common use shares, it is escapable by naming any flag explicitly, and it is
announced at `INFO` on every run it touches. It needs maintainer sign-off.

**Source.** `src/params/mod.rs` (`looks_like_10x`,
`apply_cellranger_defaults_on_10x`). STAR: `docs/STARsolo.md`, "Matching
CellRanger 4.x and 5.x results".

---

### 3.2 `--soloOutRawBarcodes Observed` (opt-in, non-STAR)

**What STAR does.** STARsolo's raw matrix has one column per whitelist
barcode, whether or not any read carried it. For 10x v3 that is 3 686 400
columns and a 62 MB `barcodes.tsv`, nearly all of it zeros.

**What rustar-aligner does.** The same, by default. `--soloOutRawBarcodes
Observed` narrows the raw matrix to the barcodes that actually hold a count,
which is what CellRanger's `raw_feature_bc_matrix` contains. On a 200-cell
fixture that is 200 columns and a 3.4 kB `barcodes.tsv`.

**Why.** Someone comparing our raw matrix against CellRanger's finds no
overlapping keys at all, because the two files mean different things by "raw".
The flag makes the comparison possible without changing what STARsolo users
get.

**Impact.** The counts are identical either way — same entries, same values,
verified on the fixture — only the columns present differ. This is a non-STAR
flag and needs maintainer sign-off; it is off by default so STARsolo parity is
untouched.

**Source.** `src/solo/count.rs` (`observed_barcodes`), `src/params/mod.rs`
(`solo_out_raw_barcodes`). CellRanger: `outs/raw_feature_bc_matrix/` from a
`cellranger count` run, observed directly rather than taken from its source.

---

## 4. Implementation divergences (no intended output difference)

These differ in *how* a result is produced, not *what* is produced. They are documented so a reviewer chasing a discrepancy knows the mechanism differs by design.

### 4.1 Transcriptome tables built on the fly

For `--quantMode TranscriptomeSAM`, rustar-aligner builds the per-transcript exon map directly from the input GTF at run time, instead of loading STAR's persisted `transcriptInfo.tab` / `exonInfo.tab` files. The projection logic mirrors STAR's `Transcriptome_quantAlign.cpp`; the output (`Aligned.toTranscriptome.out.bam`) is intended to be equivalent.

**Source.** `src/quant/transcriptome.rs`.

### 4.2 In-tree RNG generator

rustar-aligner uses an in-tree splitmix64 (`src/rng.rs`) rather than the `rand` crate, avoiding the `getrandom`/`zerocopy`/`ppv-lite86` dependency chain. This is the generator underlying §1.1; it is called out separately because it is a dependency/implementation choice independent of the tie-break policy.

---

## 5. Known residual single-read differences

These are **not** deliberate divergences — they are tracked residual diffs on the 10k yeast benchmark, kept here for completeness. Each is a single read; none is a systematic behaviour difference.

- **1 SE CIGAR-only diff** — `ERR12389696.13573895`: both tools align to XV:218357, MAPQ 255, identical score (AS=133), but place a 1-base insertion differently (`100M1I45M4S` vs STAR's `108M1I37M4S`). The 71-base seed is found at a different position within a long homopolymer run (a seed-level tie); resolving it requires reproducing STAR's exact Lmapped chain path.
- **1 STAR-only PE mate** — `ERR12389696.18919121`: an SA-level difference.
- **1 rustar-aligner-only PE mate** — `ERR12389696.6302610`: a pre-existing false positive.

See `CLAUDE.md` ("Known Issues" / "PE Status") for the current status of these.

---

## Adding a new divergence

When a change deliberately diverges from STAR (including adding a non-STAR flag, or choosing STAR's *documented* behaviour over its *actual binary* behaviour where they differ):

1. Add an entry to the appropriate section above using the **What STAR does / What rustar-aligner does / Why / Impact / Source** format.
2. Cite the STAR C++ source you checked, so the divergence can be re-verified.
3. Get maintainer sign-off in the PR — see [CONTRIBUTING.md](CONTRIBUTING.md#divergence-from-star-is-allowed--but-must-be-deliberate-and-flagged).
