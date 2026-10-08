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

### 1.2 `--soloUMIfiltering MultiGeneUMI_All` filters, rather than doing nothing

**What STAR does.** Nothing, in effect — but not because the rule is unimplemented. `SoloFeature_collapseUMIall.cpp:79-88` implements exactly the documented behaviour, zeroing every gene for any UMI seen in more than one:

```cpp
if (pSolo.umiFiltering.MultiGeneUMI_All) {
    for (auto &iu : umiGeneMapCount)
        if (iu.second.size()>1)
            for (auto &ig : iu.second) ig.second=0; //kill all genes for this UMI
};
```

The site that acts on those zeroed counts, however, gates on a different flag (`:116`, `if (pSolo.umiFiltering.MultiGeneUMI && umiGeneMapCount[...]==0)`), and `MultiGeneUMI` and `MultiGeneUMI_All` are set in mutually exclusive branches (`ParametersSolo.cpp:457-462`). Selecting `MultiGeneUMI_All` therefore zeroes the counts and then never reads them, and the run reports unfiltered counts.

**What rustar-aligner does.** The documented behaviour: a UMI seen in more than one gene is removed from all of them.

**Why.** This is a one-line wiring bug in STAR, not a design decision: STAR's own code, immediately above, computes the documented result and then discards it. Matching the binary would mean shipping a flag that silently does nothing to anyone who read either the documentation or STAR's source, which is what #144 was raised about. So the divergence is from STAR's behaviour but *not* from its intent. Single-gene UMIs are untouched, which the tests check across every mode.

**Impact.** Confined to `--soloUMIfiltering MultiGeneUMI_All`. The default (`-`) and the other filtering modes produce identical counts. Inverting the choice is a one-line change, since the test asserts the behaviour either way.

**Source.** `src/solo/count.rs` (`UmiFiltering::MultiGeneUmiAll`, `filter_multi_gene_umi`), locked by `multigene_umi_all_drops_the_umi_from_every_gene` and `multigene_umi_all_parses_to_its_own_variant`. STAR: `SoloFeature_collapseUMIall.cpp`, `ParametersSolo.cpp`.

### 1.3 `EmptyDrops_CR` Simple-Good-Turing with fewer than five distinct frequencies

**What STAR does.** The ambient profile for `--soloCellFilter EmptyDrops_CR` is smoothed with Simple Good-Turing (Elworthy's `SimpleGoodTuring/sgt.h`). `analyse()` returns early, doing nothing, when the frequency spectrum has fewer than five distinct counts — Elworthy's `MinInput` guard. `PZero`, the mass reserved for genes unseen in the ambient droplets, is neither assigned in that case nor initialised at construction, so a caller that asks for it reads whatever the stack held.

**What rustar-aligner does.** `PZero` is zero from construction.

**Why.** There is nothing to reproduce: the value STAR reads is not a decision it made. With fewer than five distinct frequencies there is no basis for reserving unseen mass, and zero says so. Reproducing STAR would mean writing code whose correct behaviour is to emit an uninitialised value, and a test asserting it.

**Impact.** Degenerate inputs only — any dataset large enough to reach the significance test has far more than five distinct frequencies. On those inputs an uninitialised read can place arbitrary mass on unseen genes, which makes the multinomial log-probabilities meaningless; zero keeps them defined.

**Source.** `src/solo/sgt.rs`, locked by `solo::sgt::tests::too_few_frequencies_leaves_the_unseen_mass_at_zero` (asserting the exact bit pattern, since the point is that nothing was written). STAR: `SoloFeature_emptyDrops_CR.cpp`, `SimpleGoodTuring/sgt.h`.

---

## 2. Cases where rustar-aligner outperforms STAR

These are not chosen divergences and not bugs: rustar-aligner reports a **higher-scoring, correct** alignment that STAR misses. They are listed here so the differential benchmark's non-exact reads are fully accounted for.

### 2.1 Four PE alignments scored better than STAR

On the 10k yeast PE benchmark, 4 reads differ in alignment score (AS) because STAR's combined-window stitching fails to place the pair at the better location:

- `ERR12389696.844151` — rustar-aligner finds VIII:451791 with 0 mismatches; STAR reports VII:1001391 with 6 mismatches.
- `ERR12389696.4972950` — rustar-aligner finds the correct **spliced** mate 2; STAR reports it unspliced.

**Impact.** rustar-aligner's result is the better alignment in each case. These are counted against exact faithfulness in the raw metric but are improvements, not regressions.

**Source.** See `STAR-RS-COMPARISON.md`, and the PE benchmark section of `CONTRIBUTING.md`.

---

## 3. Output-file metadata divergences

### 3.1 `genomeParameters.txt` command-line header

**What STAR does.** At `genomeGenerate`, STAR writes a `### <commandLineFull>` header line reproducing the full command line that built the index.

**What rustar-aligner does.** rustar-aligner echoes its own actual command line after `### ` (falling back to a parameter skeleton for API callers constructed without one). Every value line of `genomeParameters.txt` matches STAR's `genomeParametersWrite.cpp` order and tab/space formatting, including the *effective* `sjdbOverhang` (0 when the index has no sjdb, mirroring `mapGen.sjdbOverhang`).

**Why.** The header is informational; the binary path and argument spacing can never byte-match an arbitrary STAR invocation, and the index loads identically either way.

**Impact.** The `###` header line will not byte-match a STAR run (different `argv[0]` and spacing). Every other line matches byte-for-byte. No effect on alignment, index loading, or any downstream tool.

**Source.** `src/genome/mod.rs` (`genomeParameters.txt` writer).

---

### 3.1a `SAindex` N-mark bits adjacent to junction-flank k-mers

**What STAR does.** With `--sjdbGTFfile`, STAR builds the base-genome `SAindex` first (`genomeSAindex.cpp`) and then *patches* it while inserting junction-flank suffixes (`sjdbBuildIndex.cpp`). `SAiMarkNbit` marks — "suffixes for this k-mer slot may border an N" — are placed against the **base-genome** k-mer landscape: a mark lands on the last k-mer that was present *before* the junction flanks were inserted, marks are silently dropped when an inserted flank suffix takes over a slot's first-occurrence value (`sjdbBuildIndex.cpp:228-231` overwrites the packed value, flags included), and flank suffixes that touch the inter-junction spacer get marks via a separate T-fill backward-scan rule (`sjdbBuildIndex.cpp:262-284`).

**What rustar-aligner does.** rustar-aligner builds the final genome+flank text in one pass and replicates `genomeSAindex.cpp`'s serial mark semantics over that final text: the mark lands on the last k-mer *present in the final index* before the N-run.

**Why.** On GRCh38 + GENCODE v49 this changes a handful of bits (2 slots out of 357,913,940 on the measured build): exactly the slots where a k-mer became present only via a junction flank. STAR's placement there is an artifact of its incremental patch, not a semantic choice; reproducing it would mean simulating the two-phase build. Both placements are valid conservative markers — the bit only widens seed-search bounds near Ns.

**Impact.** ≤ a few bytes of the ~1.5 GB `SAindex` differ on sjdb builds (indexes built *without* a GTF are byte-identical). STAR loads either file and produces identical alignments (verified on 100k read pairs). No effect on any coordinate, count, or emitted record.

**Source.** `src/index/sa_index.rs` (`build_parallel`, `build`).

---

### 3.1b `Log.out` in the genome directory

**What STAR does.** `genomeGenerate` writes its free-form run log to `<outFileNamePrefix>Log.out` and copies it into the genome directory, so a STAR-built index always contains a `Log.out`.

**What rustar-aligner does.** The same — a STAR-shaped `Log.out` (version header, command-line/parameter sections, phase timestamps, `DONE: Genome generation, EXITING`) is written to the output prefix and copied into the genome directory.

**Impact.** The file's *content* is a run log (timestamps, host-specific paths) and can never byte-match across runs or tools; only its presence and shape are mirrored. Nothing loads it at align time.

**Source.** `src/io/log.rs` (`write_genome_generate_log`), `src/lib.rs` (`genome_generate`).

---

### 3.2 `CellReads.stats` row order

**What STAR does.** `--soloCellReadStats CB` emits its rows by iterating a libc++ `std::unordered_map`, so the order is a hash-table walk rather than a sort. At the map sizes this produces, libc++ chains new entries at the head of their bucket and walks buckets in order, which comes out as the reverse of each barcode's first appearance in read order.

**What rustar-aligner does.** Emits that same reverse-first-appearance order, including across threads: the per-read accumulator merges in read order, so a threaded run writes the same file as a serial one.

**Why.** Reproducing the order where it is reproducible costs nothing and keeps a byte-comparison against STAR usable on the sizes where it can work at all.

**Impact.** Past libc++'s load factor the map rehashes, and the order then depends on the bucket count, which depends on how many distinct barcodes were seen; beyond that size the order diverges. The **values never do** — only which line they appear on. Reading the file by barcode rather than by position is unaffected either way.

**Source.** `src/solo/cell_reads.rs`, locked by `rows_are_emitted_in_reverse_first_appearance_order`. STAR: `SoloFeature_statsOutput.cpp`.
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

### 3.3 `--soloOutLayout CellRanger` (non-STAR; on by default on 10x geometry)

**What STAR does.** STARsolo writes
`Solo.out/<feature>/{raw,filtered}/{matrix.mtx, barcodes.tsv, features.tsv}`,
uncompressed, with bare barcodes.

**What rustar-aligner does.** The same, by default, on non-10x geometry.
`--soloOutLayout CellRanger` writes the same numbers in the shape
`cellranger count` produces: `outs/{raw,filtered}_feature_bc_matrix/`, all three
files gzipped, a `-1` GEM-well suffix on every barcode, one raw column per
observed barcode, and no per-feature subdirectory when a single feature is
requested. It implies `--soloOutGzip yes`, `--soloOutRawBarcodes Observed` and
`--soloOutFileNames outs/ ...`, each still overridable on the command line.

Under §1.3 the same 10x detection turns this on by default, so a bare 10x run
lands in CellRanger's layout.

**Why.** Tools written against CellRanger's `outs/` (scanpy's `read_10x_mtx`,
Seurat's `Read10X`, any in-house loader) key on those directory names and on
the `-1` suffix. Without them the numbers are right and nothing downstream can
read them without a rename step.

**Impact.** No count changes: verified entry by entry on the 20 000-read
fixture, 13 959 entries and 15 439 counts in both layouts. On 10x geometry it
**changes where output files are written**, which needs maintainer sign-off
alongside §1.3. Against a real `cellranger count` run on the same fixture the
raw barcode sets match exactly, 200 of 200.

**Source.** `src/solo/count.rs` (`write_gene_matrix`, `write_one_barcode`),
`src/params/mod.rs` (`solo_out_layout`, `apply_cellranger_layout`). CellRanger:
`outs/` from a `cellranger count` 10.0.0 run, observed directly rather than
taken from its source.

---

### 3.4 `metrics_summary.csv` under the CellRanger layout (non-STAR)

**What STAR does.** STARsolo writes `Summary.csv`, its own metric set with its
own names. It has no `metrics_summary.csv`.

**What rustar-aligner does.** The same, by default. Under
`--soloOutLayout CellRanger` it *additionally* writes
`metrics_summary.csv` with CellRanger 10.0.0's 20 metrics, in CellRanger's
order and value formats. `Summary.csv` is written unchanged alongside it, so
nothing STARsolo-faithful is altered. Under the CellRanger layout the older
`CellRanger.summary.csv` is not written: its four rows are a subset of the
metrics file.

**Why.** The file is how a 10x pipeline reads run quality. Its 20 names are
also the clearest public statement of what CellRanger measures, which makes it
useful as a target even where our value differs.

**Impact.** No count changes. It costs one extra gene-body overlap query per
read, because the exonic/intronic split needs it and a `Gene`-only run does not
otherwise do it.

**Not all 20 are certainties.** Twelve match a real `cellranger count` 10.0.0
run exactly on the 20 000-read fixture. The other eight follow from our own
tallies under a stated interpretation, because CellRanger does not document the
denominators:

* `Reads Mapped Confidently to *` reads "confidently" as MAPQ 255, i.e. our
  uniquely-mapped set.
* `Reads Mapped Confidently to Transcriptome` is the `Gene` feature's own
  uniquely-assigned read tally.
* `Valid UMI Sequences` is measured over reads that reached the UMI check, i.e.
  those with a valid barcode.
* `Sequencing Saturation` is `1 - molecules / reads` over the reads that
  entered the matrix. This is the largest disagreement on the fixture: 12.7%
  against CellRanger's 7.4%. 10x define it as
  `1 - n_deduped_reads / n_reads`, with `n_deduped_reads` the number of unique
  `(barcode, UMI, gene)` combinations among confidently mapped reads. Taking
  that literally — counting distinct triples *before* UMI correction — gives
  **0.0%** on this fixture, because no two reads here share an exact triple, so
  the literal reading is ruled out and their numerator is the corrected
  molecule count, as ours is. The residual is therefore the denominator: at
  7.4% theirs implies about 16 300 reads where ours counts about 17 300. The
  difference is which reads are "confidently mapped", not the formula.
* `Mean Reads per Cell` is total reads over called cells, not reads-in-cells
  over called cells, which is what reproduces CellRanger's value.

The cell-count-dependent metrics (`Median Genes per Cell`,
`Median UMI Counts per Cell`, `Total Genes Detected`) differ by 2-4 on the
fixture, following the count differences recorded in §1.3 rather than a
different definition.

**Source.** `src/solo/count.rs` (`write_metrics_summary`, `metric_int`,
`metric_pct`), `src/solo/mod.rs` (`Q30Stats`). CellRanger: the
`metrics_summary.csv` of a `cellranger count` 10.0.0 run, observed directly
rather than taken from its source.

---

### 3.5 `--soloCellFilter OrdMag`, and EmptyDrops_CR's initial cell set (non-STAR)

**What STAR does.** STARsolo's cell calling is `CellRanger2.2`: take the
`quantile` of the top `nExpectedCells` barcodes by UMI count and call every
barcode holding at least that over `maxMinRatio`, with `nExpectedCells` fixed
at 3 000. `EmptyDrops_CR` uses the same knee for the set of cells it guarantees
before the Monte-Carlo rescue.

**What CellRanger does.** The same rule, but it *searches* for the expected
cell count instead of assuming one: it minimises `(OrdMag(x) - x)^2 / x` over
`x` from 2 to about 45 000, where `OrdMag(x)` is the number of cells the rule
calls when told to expect `x`. The loss is small where the rule predicts
itself. EmptyDrops then rescues barcodes below that cutoff.

**What rustar-aligner does.** `--soloCellFilter OrdMag maxExpectedCells
quantile ratio` (default `45000 0.99 10`) implements the search, and
`EmptyDrops_CR` now uses it for its initial set, which is the order CellRanger
runs the two steps in. `CellRanger2.2` is untouched and remains the default.

**Two things the 10x page does not specify, decided here.** Every integer in
range is evaluated rather than a geometric grid, so the result is the true
minimum of the stated loss rather than a nearby grid point; each evaluation is
a binary search over the sorted totals, so an exhaustive sweep costs nothing
measurable. And ties go to the smaller `x`, so the result does not depend on
iteration order — determinism, as everywhere else in this codebase.

**Impact.** `EmptyDrops_CR`'s initial cell set can differ from what it was.
On the 20 000-read fixture it does not: `CellRanger2.2`, `OrdMag` and
`EmptyDrops_CR` all call 200 cells, as does CellRanger 10.0.0. That fixture has
a clean plateau, where any of these rules works; the two rules part company on
a graded distribution with no plateau, which is what the unit tests cover.

**Source.** `src/solo/count.rs` (`ordmag_at`, `ordmag_threshold`). CellRanger:
the Gene Expression algorithm page, "Cell Calling", read rather than taken from
its source. Coverage of the rest of that page is tracked in #181.

---

## 4. Implementation divergences (no intended output difference)

These differ in *how* a result is produced, not *what* is produced. They are documented so a reviewer chasing a discrepancy knows the mechanism differs by design.

### 4.1 Transcriptome tables built on the fly

For `--quantMode TranscriptomeSAM`, rustar-aligner builds the per-transcript exon map directly from the input GTF at run time, instead of loading STAR's persisted `transcriptInfo.tab` / `exonInfo.tab` files. The projection logic mirrors STAR's `Transcriptome_quantAlign.cpp`; the output (`Aligned.toTranscriptome.out.bam`) is intended to be equivalent.

**Source.** `src/quant/transcriptome.rs`.

### 4.2 In-tree RNG generator

rustar-aligner uses an in-tree splitmix64 (`src/rng.rs`) rather than the `rand` crate, avoiding the `getrandom`/`zerocopy`/`ppv-lite86` dependency chain. This is the generator underlying §1.1; it is called out separately because it is a dependency/implementation choice independent of the tie-break policy. It is not the only in-tree generator: `--soloCellFilter EmptyDrops_CR` samples with a bit-exact libc++ `mt19937` (`src/solo/libcxx_rng.rs`) so its Monte-Carlo null matches STAR's — a convergence with STAR rather than a divergence from it.

### 4.3 Coordinate-sort memory budget and spilling

**What STAR does.** `--outSAMtype BAM SortedByCoordinate` collects alignments into `--outBAMsortingBinsN` coordinate bins on disk under `outFileTmp` (`--outTmpDir`, else `<prefix>_STARtmp/`), then sorts each bin in memory. `--limitBAMsortRAM 0` is replaced by the genome size plus the SA and SAindex sizes (`STAR.cpp:234-236`). If any single bin needs more than that, STAR exits with a fatal error and asks to be re-run with a larger value (`bamSortByCoordinate.cpp:29-34`).

**What rustar-aligner does.** An external merge sort. Records fill a buffer of `--limitBAMsortRAM` bytes (estimated), which is sorted and spilled as a run when full; the runs are k-way merged at the end, at most 64 open at a time. `--limitBAMsortRAM 0` means 512 MiB. The run never fails for lack of sort memory. Spill runs go to `--outTmpDir` when given, otherwise beside the output, and are always removed (`--outTmpKeep` has no effect).

**Why.** Peak memory stays bounded whatever the output size, without STAR's failure mode. A genome-size default would also be unknown to the writer when it is created, and an unbounded buffer was what this replaced.

**Impact.** None on output: records are byte-identical to an unbounded in-memory sort, ties keep input order (as STAR's do), and this holds through multi-pass merges. Only memory use, temporary files and the absence of the out-of-memory error differ.

**Source.** `src/io/bam.rs` (`CoordinateSorter`).

## 5. Known residual single-read differences

These are **not** deliberate divergences — they are tracked residual diffs on the 10k yeast benchmark, kept here for completeness. Each is a single read; none is a systematic behaviour difference.

- **1 SE CIGAR-only diff** — `ERR12389696.13573895`: both tools align to XV:218357, MAPQ 255, identical score (AS=133), but place a 1-base insertion differently (`100M1I45M4S` vs STAR's `108M1I37M4S`). The 71-base seed is found at a different position within a long homopolymer run (a seed-level tie); resolving it requires reproducing STAR's exact Lmapped chain path.
- **1 STAR-only PE mate** — `ERR12389696.18919121`: an SA-level difference.
- **1 rustar-aligner-only PE mate** — `ERR12389696.6302610`: a pre-existing false positive.

See `ROADMAP.md` and `STAR-RS-COMPARISON.md` for the current status of these.

---

## Adding a new divergence

When a change deliberately diverges from STAR (including adding a non-STAR flag, or choosing STAR's *documented* behaviour over its *actual binary* behaviour where they differ):

1. Add an entry to the appropriate section above using the **What STAR does / What rustar-aligner does / Why / Impact / Source** format.
2. Cite the STAR C++ source you checked, so the divergence can be re-verified.
3. Get maintainer sign-off in the PR — see [CONTRIBUTING.md](CONTRIBUTING.md#divergence-from-star-is-allowed--but-must-be-deliberate-and-flagged).
