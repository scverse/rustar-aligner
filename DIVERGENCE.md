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

### 1.2 Homopolymer UMIs are rejected whichever base repeats

**What STAR does.** STAR means to reject a UMI made of one repeated base, and under `CB_UMI_Simple` it does. Under `CB_UMI_Complex` it rejects only poly-A, letting poly-C, poly-G and poly-T through.

The mechanism is an initialisation-order defect, not a rule. The four packed constants are built in the `SoloReadBarcode` constructor (`SoloReadBarcode.cpp:16-21`):

```cpp
for (uint32 jj=0;jj<4;jj++) {
    homoPolymer[jj]=0;
    for (uint32 ii=0; ii<pSolo.umiL;ii++)
        homoPolymer[jj]=(homoPolymer[jj]<<2)+jj;
};
```

but `CB_UMI_Complex` does not know its UMI length until it has extracted one, and sets it at read-processing time (`SoloReadBarcode_getCBandUMI.cpp:353-354`):

```cpp
if ( pSolo.umiL == 0 )
    pSolo.umiL = umiSeq.size();
```

So when `umiL` is still 0 at construction, the loop body never runs and all four constants stay `0`. The check (`SoloReadBarcode_getCBandUMI.cpp:139`) then compares the packed UMI against `0` four times. Poly-A packs to `0`, because `A` is nucleotide 0, so poly-A is still caught; poly-C, poly-G and poly-T pack to non-zero values and match nothing.

**What rustar-aligner does.** Rejects all four homopolymers, in every barcode geometry.

**Why.** A UMI of one repeated base carries no information whichever base it is, so rejecting only poly-A is not a threshold choice to be matched but an accident of when a length becomes known. STAR's own `CB_UMI_Simple` path, where `umiL` is set from the parameters before construction, rejects all four; the complex path is the outlier. Matching the accident would mean counting poly-T UMIs as real molecules.

**Impact.** Under `--soloType CB_UMI_Complex`, reads whose UMI is a non-A homopolymer are dropped here and counted by STAR. On real data these are rare and are sequencing artefacts. Note this is not a new behaviour in this PR: rustar-aligner has always rejected all four. What is new is that it is now deliberate, locked by a test, and recorded here rather than being an undocumented difference.

**Source.** `src/solo/whitelist.rs` (`check_umi`, `umi_valid_rejects_every_homopolymer`). STAR: `SoloReadBarcode.cpp:16-21`, `SoloReadBarcode_getCBandUMI.cpp:139` and `:353-354`.

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
