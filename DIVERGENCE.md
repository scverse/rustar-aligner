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

**Impact.** With the default `--outMultimapperOrder Old_2.4`, **no RNG is consulted at all**: the primary is the deterministic best alignment. The first two keys are STAR's own (`ReadAlign_stitchPieces.cpp:340` compares `maxScore`, then `gLength`). Where those two tie, STAR takes the earliest window in its iteration order and rustar-aligner takes the smallest genomic position; that last key is a divergence, and it is the one described below. The RNG divergence proper is observable only under `--outMultimapperOrder Random`, and in both cases only in *which* equal-scoring locus is marked primary, never in the set of reported alignments.

**On the residual ties.** STAR's window order is deterministic per read (anchor pieces in `PC` order, positions in suffix-array order) and rustar-aligner also builds windows per read, so reproducing it would not cost thread-invariance. It has been measured, not assumed: substituting seed discovery order for the positional key gained 21 reads on the annotated-junction tier and lost 142 on the 10k single-end tier, because this codebase's anchor iteration order is not STAR's `PC` order. Closing the gap means matching how MMP results are recorded, which has not been done. Until it is, the positional key stays, because it is total, cheap and thread-invariant.

This is the reason faithfulness is reported **tie-adjusted**. On the 10k yeast benchmark, 130 of the 131 differing single-end records and 197 of the 202 differing paired-end mate records are genuine ties: both tools find the identical alignment set and differ only in which equal-scoring member is primary. On the annotated-junction tier all 42 remaining differences are of this kind, which is what its 100.000% tie-adjusted figure means. The handful that are not ties are in [§5](#5-known-residual-single-read-differences). Raw counts are reported alongside, never the tie-adjusted figure alone.

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

- **1 SE insertion-placement diff**: `ERR12389696.20597455`, XIV:545446, MAPQ 255, identical score (AS=143), with the 1-base insertion four bases to the left of STAR's choice (`25M1I124M` against STAR's `29M1I120M`). The same read shows on the PE tier. Insertion placement inside `stitchAlignToTranscript`, not a tie-break.
- **1 PE pair where rustar-aligner scores lower**: `ERR12389696.11539725`, AS 224 against STAR's 235; mate 1 is soft-clipped to `31S95M459N24M` where STAR reaches back across a second junction to `22S10M468N94M459N24M`.
- **1 PE pair where rustar-aligner scores higher**: `ERR12389696.4972950`, AS 260 against STAR's 248; mate 2 is spliced `1S33M72N50M186N65M1S` where STAR soft-clips 27 bases to `27S122M1S`. Recorded here rather than in [§2](#2-cases-where-rustar-aligner-outperforms-star) because it is a single observed read, not a characterised behaviour.

Counts are from the 10k yeast SE and PE tiers with an index built by the same binary; the annotated-junction tier has no non-tie residual.

See `CLAUDE.md` ("Known Issues" / "PE Status") for the current status of these.

---

## Adding a new divergence

When a change deliberately diverges from STAR (including adding a non-STAR flag, or choosing STAR's *documented* behaviour over its *actual binary* behaviour where they differ):

1. Add an entry to the appropriate section above using the **What STAR does / What rustar-aligner does / Why / Impact / Source** format.
2. Cite the STAR C++ source you checked, so the divergence can be re-verified.
3. Get maintainer sign-off in the PR — see [CONTRIBUTING.md](CONTRIBUTING.md#divergence-from-star-is-allowed--but-must-be-deliberate-and-flagged).
