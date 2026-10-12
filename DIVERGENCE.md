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

**Impact.** With the default `--outMultimapperOrder Old_2.4`, **no RNG is consulted at all**. Alignments are written in STAR's `trMult` order (window order) and the primary is STAR's `trBest`: best score, then shorter genomic length, then the earlier window (`ReadAlign_multMapSelect.cpp:58-91`, `ReadAlign_stitchPieces.cpp:328-329`). On the yeast 10k and 100k benchmarks, SE and PE, every SAM record, including the primary flag and `HI` order, is identical to STAR's. The divergence is observable only under `--outMultimapperOrder Random`, and only in *which* equal-scoring locus is marked primary and in what order the alignments are written — never in the set of reported alignments.

**Transcriptome BAM primary.** With `--quantMode TranscriptomeSAM`, STAR marks one of a read's transcript alignments primary with a draw from the same per-thread stream, `int(u * nAlignT)`, for every mapped read and whatever `--outMultimapperOrder` says (`ReadAlign_quantTranscriptome.cpp:69`). rustar-aligner draws it from STAR's generator (`mt19937` seeded with `--runRNGseed`, `ReadAlign.cpp:11-12`) in the single writer stage, in STAR's read order: the BySJout held reads come after all the others, as in STAR's second stage. That is the stream of STAR's chunk 0, so `Aligned.toTranscriptome.out.bam` is record for record STAR's `--runThreadN 1` output at any `--runThreadN` here. With more threads STAR's own draws depend on which thread maps which chunk. Under `--outMultimapperOrder Random`, STAR's shuffles take numbers from the same stream first, which rustar-aligner does not, so the transcriptome primaries differ there too.

**Source.** `src/rng.rs`, `src/align/read_align.rs` (`per_read_seed`, `shuffle_tied_prefix`, `mult_map_order`), `src/params/mod.rs` (`MultimapperOrder`), `src/lib.rs` (`TrRecords::finalize`, `tr_primary_rng`). STAR: `ReadAlign_multMapSelect.cpp`, `ReadAlign_quantTranscriptome.cpp`, `ReadAlignChunk` RNG seeding.

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

---

### 3.3 SAM/BAM header program lines

**What STAR does.** The header names STAR: `@PG ID:STAR PN:STAR VN:2.7.11b CL:<commandLineFull>` (`samHeaders.cpp:64`), where `commandLineFull` is STAR's re-formatted parameter list, and `@CO user command line: <command line>` (`:84`).

**What rustar-aligner does.** The same two lines, in the same place, name rustar-aligner: `@PG ID:rustar-aligner PN:rustar-aligner VN:<version> CL:<command line>` and `@CO user command line: <command line>`. Every other header line follows STAR's layout and text, including `@HD VN:1.4`, a user `@HD`/`@PG`, the comment file and `@RG` (`samHeaders.cpp:29-100`).

**Why.** The program record says which program wrote the file; claiming to be STAR would be wrong, and the command line differs anyway (`argv[0]`, spacing).

**Impact.** On the yeast 10k runs (SE and PE, SAM, unsorted, sorted and both BAMs, user `@HD`/`@PG`, comment file, read groups) the headers are identical to STAR's apart from these two lines. A tool that detects the aligner from `PN:STAR` will not recognise the file as STAR's.

**Source.** `src/io/sam.rs` (header text), `src/io/bam.rs`.

---

## 4. Implementation divergences (no intended output difference)

These differ in *how* a result is produced, not *what* is produced. They are documented so a reviewer chasing a discrepancy knows the mechanism differs by design.

### 4.1 Transcriptome tables built on the fly (GTF given at mapping time)

With `--sjdbGTFfile` at mapping, STAR writes the GTF's transcript tables to `<prefix>_STARgenome/` and then reads them back (`Transcriptome.cpp:13`). rustar-aligner writes the same tables but builds its in-memory transcript index and gene model directly from the GTF. Without a GTF at mapping, both read the index's `transcriptInfo.tab` / `exonInfo.tab` / `geneInfo.tab`, for `--quantMode TranscriptomeSAM` / `GeneCounts`, STARsolo `Gene` / `GeneFull` / `Velocyto` / SmartSeq and the `GX`/`GN` tags, and exit with STAR's "could not open ... geneInfo.tab" error when neither a GTF nor the tables exist. The output is intended to be the same.

**Source.** `src/index/io.rs` (`GenomeIndex::load`), `src/quant/transcriptome.rs`, `src/quant/mod.rs` (`resolve_gene_annotation`).

### 4.2 In-tree RNG generator

rustar-aligner uses an in-tree splitmix64 (`src/rng.rs`) rather than the `rand` crate, avoiding the `getrandom`/`zerocopy`/`ppv-lite86` dependency chain. This is the generator underlying §1.1; it is called out separately because it is a dependency/implementation choice independent of the tie-break policy. It is not the only in-tree generator: `--soloCellFilter EmptyDrops_CR` and the transcriptome BAM primary (§1.1) draw from a bit-exact `mt19937` (`src/solo/libcxx_rng.rs`) so their results match STAR's — a convergence with STAR rather than a divergence from it.

### 4.3 Coordinate-sort memory budget and spilling

**What STAR does.** `--outSAMtype BAM SortedByCoordinate` collects alignments into `--outBAMsortingBinsN` coordinate bins on disk under `outFileTmp` (`--outTmpDir`, else `<prefix>_STARtmp/`), then sorts each bin in memory. `--limitBAMsortRAM 0` is replaced by the genome size plus the SA and SAindex sizes (`STAR.cpp:234-236`). If any single bin needs more than that, STAR exits with a fatal error and asks to be re-run with a larger value (`bamSortByCoordinate.cpp:29-34`).

**What rustar-aligner does.** An external merge sort. Records fill a buffer of `--limitBAMsortRAM` bytes (estimated), which is sorted and spilled as a run when full; the runs are k-way merged at the end, at most 64 open at a time. `--limitBAMsortRAM 0` means 512 MiB. The run never fails for lack of sort memory. Spill runs go to `--outTmpDir` when given, otherwise beside the output, and are always removed (`--outTmpKeep` has no effect).

**Why.** Peak memory stays bounded whatever the output size, without STAR's failure mode. A genome-size default would also be unknown to the writer when it is created, and an unbounded buffer was what this replaced.

**Impact.** None on output: records are byte-identical to an unbounded in-memory sort, ties keep input order (as STAR's do), and this holds through multi-pass merges. Only memory use, temporary files and the absence of the out-of-memory error differ.

**Source.** `src/io/bam.rs` (`CoordinateSorter`).

### 4.4 Gzipped GTF files

**What STAR does.** `--sjdbGTFfile` is read as plain text; a gzipped GTF fails to parse.

**What rustar-aligner does.** A GTF whose first two bytes are the gzip magic (`1f 8b`) is decompressed while it is read, at `genomeGenerate`, at mapping and in `--runMode liftOver`; the file name does not matter, and multi-member gzip is read to the end. A plain-text GTF is read as before.

**Why.** Annotations are often shipped compressed (10x references: `genes/genes.gtf.gz`), and decompressing them first adds a step and a large file for no change in the result.

**Impact.** None on output: the same GTF gzipped or not gives identical indexes, `_STARgenome/` tables, alignments, `ReadsPerGene.out.tab` and transcriptome BAM. Only a run STAR would reject is accepted.

**Source.** `src/junction/gtf.rs` (`open_annotation`).

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
