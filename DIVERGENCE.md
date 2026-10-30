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

### 1.4 Bulk total-RNA options: unspliced transcriptome targets, `sp` tag, `--quantGeneSplicing`

All three are rustar-aligner extensions with no STAR equivalent, and all are opt-in. Each has its own flag rather than a new value of a STAR parameter (`--quantMode`, `--outSAMattributes`), so a STAR command line never means something different here. With the defaults (`--quantTranscriptomeUnspliced None`, `--outSAMsplicingStatus No`, `--quantGeneSplicing No`) no output changes, which `tests/bulk_unspliced.rs` (`new_options_leave_existing_outputs_unchanged`) checks.

**What STAR does.** `--quantMode TranscriptomeSAM` projects every alignment onto every annotated transcript whose exons contain it (`Transcriptome_quantAlign.cpp`); there are no unspliced targets, so an unspliced read inside an intron that another isoform retains is projected onto the retaining isoform only, and a read in a constitutive intron is not projected at all. STAR classifies reads as spliced / unspliced / ambiguous only for single-cell runs (STARsolo, `Transcriptome_classifyAlign.cpp` + `SoloFeature_countVelocyto.cpp`), and has no per-alignment splicing tag.

**What rustar-aligner does.**

- `--quantTranscriptomeUnspliced Intron|PreMRNA` appends one `<gene_id>-I` target per gene (merged introns of all isoforms plus `--quantTranscriptomeUnsplicedFlank`, or the gene body) after the annotated transcripts, and projects alignments onto them with the unchanged STAR projection, except that fragments crossing a splice junction are kept off the unspliced targets. It also writes `Aligned.toTranscriptome.targets.tsv` and, on request, the unspliced sequences.
- `--outSAMsplicingStatus Yes` and `--quantGeneSplicing Yes` run a port of STARsolo's spliced / unspliced classification (`alignToTranscriptMinOverlap` with `minOverlapMinusOne = 6` and the 1 Mb intron cap, then the per-UMI collapse of `countVelocyto`), one read or pair standing for one UMI. Two details differ from `classifyAlign`: the containment test uses the true leftmost / rightmost aligned base of the pair, where STAR uses the first block's start and the last block's end (a `TODO` next to that line in STAR flags the case where mate 2 ends before mate 1); and blocks are sorted by position before the scan, so STAR's early exit at the last exon also holds for overlapping mates. The `sp` tag collapses over all containing transcripts regardless of gene and strand.

**Why.** In ribo-depleted total RNA a large share of reads is pre-mRNA. Without unspliced targets STAR's projection gives the intronic ones to retained-intron isoforms, so isoform proportions follow the library's pre-mRNA content (see the PR for measurements, and COMBINE-lab/salmon#1229 for the same effect with Salmon decoys). Unspliced targets let the downstream EM share such reads, as splici does for single-cell data. Bulk users have no STAR option that reports how much of a library is unspliced.

**Impact.** Only when the options are given. Unspliced targets add `@SQ` lines and records to `Aligned.toTranscriptome.out.bam` and change `NH` / `HI` / `MAPQ` of the reads they receive; `sp` adds one tag per record.

**Source.** `src/quant/transcriptome.rs` (`with_unspliced_targets`, `unspliced_intervals`, `filter_and_project`), `src/quant/splice_status.rs`, `src/quant/mod.rs` (`QuantContext`), `src/params/mod.rs` (the three flags). STAR: `Transcriptome_quantAlign.cpp`, `Transcriptome_classifyAlign.cpp`, `SoloFeature_countVelocyto.cpp`.

### 1.5 DegNorm degradation normalization (`--quantGeneCoverage`, `--degNormCoverageFiles`)

A rustar-aligner extension with no STAR equivalent, opt-in. Like §1.4 it uses its own flags rather than new values of `--quantMode` or `--runMode`.

**What STAR does.** Nothing: STAR has no transcript-degradation model and no per-gene coverage output.

**What rustar-aligner does.** `--quantGeneCoverage Yes` accumulates per-gene, per-exonic-base coverage during alignment and writes `GeneCoverage.out.bin`. A separate run with `--degNormCoverageFiles` (two or more files, no reads) merges the samples and fits the rank-one NMF over-approximation of [DegNorm](https://nustatbioinfo.github.io/DegNorm/) (Xiong et al., *Genome Biology* 2019), writing Degradation Index scores and degradation-adjusted counts to `DegNorm.out/`.

**Why.** Degradation bias is gene- and sample-specific, so a global size factor cannot remove it. DegNorm normally recovers coverage curves by re-reading sorted, indexed BAMs; the aligner already assigns each unique read to a gene for `ReadsPerGene.out.tab`, so the curves cost almost nothing to produce there.

**Impact.** Only when the options are given; alignment output and `ReadsPerGene.out.tab` are unchanged (`tests/degnorm.rs`). A degradation index is defined across samples, so the merge needs at least two coverage files.

**Relation to the Python DegNorm.** The fit is a port of `degnorm/nmf.py` (upstream commit 1023383): gene filter (`__main__.py:226-232`), rank-one NMF over-approximation, baseline selection, DI, scale factors, outer loop and adjusted counts. Given the same coverage matrices, DI agrees within 7e-6 and adjusted counts within 3e-5 relative (human chr21, four samples, default, `--degNormSkipBaselineSelection` and non-default iteration/bin settings). The remaining differences are float-level: upstream takes the leading singular triplet with `svds`, rustar-aligner by power iteration on the `p x p` Gram matrix. With `--degNormDownsampleRate` above 1, upstream draws each gene's sampling offset from a global random stream; rustar-aligner derives it from `--runRNGseed` and the gene index, so a run does not depend on thread count. The coverage is not DegNorm's: its BAM reader drops reads with any base outside the annotated exons, advances the reference on soft clips and insertions, assigns reads in overlapping genes by containment in one gene's exons, and treats a BAM as paired only when read names end in `.1`/`.2` (`reads.py`), so it reads STAR-style BAMs as single-end. rustar-aligner counts the exonic bases of unique reads (fragments) that overlap the exons of exactly one gene, with mate overlap counted once. On the same alignments the two coverage definitions give different matrices for most genes, so DI scores differ from a Python DegNorm run on the BAMs.

**Source.** `src/quant/coverage.rs`, `src/degnorm/` (ported from DegNorm's `degnorm/nmf.py`, MIT), `src/params/mod.rs` (the flags).

### 1.6 Homopolymer UMIs are rejected whichever base repeats

**What STAR does.** STAR means to reject a UMI made of one repeated base, and under `CB_UMI_Simple` it does. Under `CB_UMI_Complex` it rejects only poly-A, letting poly-C, poly-G and poly-T through.

The mechanism is an initialisation-order defect, not a rule. The four packed constants are built in the `SoloReadBarcode` constructor, after `ParametersSolo.cpp:368` has set `umiL = 0` for this type (`SoloReadBarcode.cpp:16-21`):

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

**Why.** Rejecting only poly-A is an accident of when the UMI length becomes known, not a rule, on three pieces of evidence:

1. **The code builds all four.** The constructor loops over all four bases to build one constant per homopolymer; it only ends up matching poly-A because `umiL` is still 0 (`ParametersSolo.cpp:368`: `umiL = 0; //this will be defined when the first barcode is processed`) and `A` encodes as 0.
2. **STAR's documentation says all homopolymers.** `--soloUMIfiltering -` is "basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)" (`parametersDefault:888`, `docs/STARsolo.md:436`), with no exception by base or by `--soloType`. CellRanger removes all four.
3. **The same check rejects all four under `CB_UMI_Simple`**, where `umiL` comes from `--soloUMIlen` before the constructor runs. Nothing in the source or the commit history gives `CB_UMI_Complex` a rule of its own.

Poly-A being special in RNA does not apply here: the UMI comes from the barcode read, not the cDNA. In oligo-dT chemistries the UMI sits next to the dT, so a UMI that runs into it reads as poly-T, the one homopolymer STAR's complex path lets through. Matching the accident would count poly-T UMIs as molecules.

**Impact.** Under `--soloType CB_UMI_Complex`, reads whose UMI is a non-A homopolymer are dropped here and counted by STAR. On real data these are rare and are sequencing artefacts. The stats files show it too: such reads count as `noUMIhomopolymer` in `Barcodes.stats` here, and under their barcode's match (or `noNoWLmatch`) in STAR's, which also moves `Features.stats` and the "Reads With Valid Barcodes" line of `Summary.csv`. On the 10x v2 mouse chr19 fixture run as `CB_UMI_Complex`, 43 reads; the count matrices are unchanged.

**Source.** `src/solo/whitelist.rs` (`check_umi`, `umi_valid_rejects_every_homopolymer`). STAR: `SoloReadBarcode.cpp:16-21`, `SoloReadBarcode_getCBandUMI.cpp:139` and `:353-354`.

### 1.7 `--clipAdapterType CellRanger4` with `--clip5pNbases`: records that describe the aligned bases

**What STAR does.** With `CellRanger4`, the 5' clip of mate 1 is `ClipMate` type 10. `ClipMate::clip` (`ClipMate_clip.cpp`) runs in three steps:

- **Fixed bases (`:14-24`).** `Lread -= N` and `clippedN += N` always run, but the `memmove` that drops the N bases from the start of the sequence runs only `if (type==0)`. Type 0 is the 5' clip of `--clipAdapterType Hamming`; the guard was not extended to type 10, so the bases are counted as clipped but stay in the buffer.
- **TSO (`:50-54`).** Type 10 then moves the sequence by the TSO length alone, counted from the original start of the read.
- **Bases after the adapter (`:64-70`).** `--clip5pAfterAdapterNbases` has the same `type==0` guard.

`readLoad.cpp:60-61` aligns what is left. So STAR aligns read bases `[TSO, len−N)` (before the 3' clip) where it means to align `[N+TSO, len)`: N bases it was told to discard are aligned, and N real bases at the 3' end are lost. The record is then written with a 5' soft clip of `clippedN` = N+TSO (`ReadAlign_alignBAM.cpp:224-230`, `ReadAlign_calcCIGAR.cpp:16-22`), which places SEQ N bases away from the bases that were aligned at POS. The 3' poly-A search (`:55-57`) also runs on the shifted buffer.

**What rustar-aligner does.** Removes the fixed bases, then the TSO, then the bases after it, and at the 3' end the fixed bases, then the poly-A tail (STAR's amounts and order, `clip::clip_mate_cr4`), and aligns the bases that remain. Each record describes the alignment that was made.

**Why.** STAR's records are not internally consistent, and its own `nM` tag shows it. On the 10x v2 mouse chr19 fixture with `--clipAdapterType CellRanger4 --clip5pNbases 5 --clip3pNbases 3`, recomputing mismatches against the genome from each primary record's POS, CIGAR and SEQ:

| | primary records | recomputed = `nM` as written | recomputed = `nM` with SEQ shifted 5 bases | 0 mismatches as written |
|---|---|---|---|---|
| STAR | 16,880 | 2 | 16,880 | 0 |
| rustar-aligner | 16,936 | 16,936 | 13 | 7,501 |

Every STAR record matches the alignment STAR computed only once SEQ is moved by the N bases that were never removed. rustar-aligner's records match their `nM` as written. Nothing in STAR's documentation describes the shift, and no report of it was found in STAR's issue tracker (searched October 2026).

**Impact.** Only `CellRanger4` with a non-zero `--clip5pNbases` (and, by the same code, `--clip5pAfterAdapterNbases`, which has not been run). Without them the two agree (fixture: all 21,732 records identical). With them, every read differs: of the 16,177 reads both map, 9,353 differ only in POS (−N on forward reads, +N on reverse, relative to rustar-aligner), and 6,824 also differ in CIGAR, locus or strand, because STAR aligned a different stretch of the read. The gene matrix differs too (fixture: 21 raw entries, 448 against 443 UMIs).

**Source.** `src/clip/mod.rs` (`clip_mate_cr4`). STAR: `ClipMate_clip.cpp:14-24, 50-57, 64-70`, `readLoad.cpp:60-61`, `ReadAlign_alignBAM.cpp:224-230`, `ReadAlign_calcCIGAR.cpp:16-22`.

### 1.8 `--soloCellRangerSummary`: a CellRanger-style mapping summary

A rustar-aligner extension with no STAR equivalent, off by default (`No`). With the default no output changes.

**What STAR does.** STARsolo writes one `Summary.csv` per feature (`SoloFeature_statsOutput.cpp`), with reads mapped to the genome and to the feature, but no split of the reads by region.

**What rustar-aligner does.** With `--soloCellRangerSummary Yes` and both `Gene` and `GeneFull` in `--soloFeatures`, it also writes `Solo.out/<feature>/CellRanger.summary.csv`: of all input reads, the fractions that map uniquely and fall in exonic, intronic or intergenic regions, or antisense to a gene, under the names CellRanger's `metrics_summary.csv` uses ("Reads Mapped Confidently to Exonic Regions" and so on). The definitions are rustar-aligner's own (exonic and intronic come from the two features' overlaps) and are not claimed to equal CellRanger's. `Summary.csv` itself is never changed.

**Why.** CellRanger users look at these fractions, and the two overlaps needed are computed anyway when both features run.

**Source.** `src/solo/count.rs` (`write_cellranger_summary`), `src/solo/gene.rs` (`Region`).

---

## 2. Cases where rustar-aligner outperforms STAR

These are not chosen divergences and not bugs: rustar-aligner reports a **higher-scoring, correct** alignment that STAR misses.

There are none at present. The four yeast PE reads once listed here now match STAR exactly, as does every record of the yeast 10k and 100k SE and PE benchmarks.

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

### 3.4 `CB_UMI_Complex`: the `Barcodes.stats` bin of two kinds of rejected read

**What STAR does.** Under `CB_UMI_Complex`, STAR checks the UMI first and then matches each cell-barcode segment to its own whitelist, in order (`SoloReadBarcode_getCBandUMI.cpp:356-420`). The first failure sets `cbMatch`, which picks the read's `Barcodes.stats` bin. Two consequences:
- a read with one mismatch in each of two segments is `cbMatch` -12, counted as `noTooManyMM`;
- a read with a rejected UMI keeps the UMI's code (`noNinUMI` or `noUMIhomopolymer`) even when a segment cannot be extracted.

**What rustar-aligner does.** It matches the concatenated segments against the product of the segment whitelists. A read with mismatches in two segments is two mismatches from every entry, so it is counted as `noNoWLmatch`. A barcode read too short for a segment is counted as `noNoCB` before its UMI is checked.

**Why.** Not chosen: per-segment matching is not ported yet. The product whitelist gives the same match as STAR for a read with at most one mismatched segment, which is every read STAR counts.

**Impact.** Only the bins of reads that both tools reject: the count matrices, `Features.stats` and `Summary.csv` are unchanged. On the 10x v2 mouse chr19 fixture run as `CB_UMI_Complex`, `Barcodes.stats` differs from STAR's only by the reads of §1.6.

**Source.** `src/solo/mod.rs` (`SoloBarcodeLayout::extract`), `src/solo/whitelist.rs` (`CbWhitelist::load_complex`). STAR: `SoloReadBarcode_getCBandUMI.cpp:356-420`.

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

### 4.5 Compressed FASTQ input detected by magic bytes

**What STAR does.** Reads `--readFilesIn` as plain text unless `--readFilesCommand` names a decompressor (`zcat`, `bunzip2 -c`, ...), which STAR runs per file (`Parameters_openReadsFiles.cpp:21-70`). It does no detection of its own.

**What rustar-aligner does.** Without `--readFilesCommand`, it peeks each input's first bytes and decodes gzip (always), and bzip2, Zstandard and XZ when built with the `bz2`, `zstd` and `xz` Cargo features. Anything else is read as plain text, whatever the file name. A recognised format whose feature is not compiled in is an error, not garbage fed to the parser. With `--readFilesCommand`, no detection happens, exactly as in STAR.

**Why.** Extension-based gzip detection was already a local convenience; content detection makes it independent of file names and covers stdin and pipes.

**Impact.** Input STAR would reject (compressed data without `--readFilesCommand`) is accepted. No alignment output changes for input both accept.

**Source.** `src/io/compression.rs`, called from `FastqReader::open` in `src/io/fastq.rs`.

## 5. Known residual single-read differences

These are **not** deliberate divergences. They are the differences `test/parity_suite.sh` finds, kept here so that each one is fixed or explained before v0.2.0 (README, *Comparison with STAR 2.7.11b*). Read names are from that suite's datasets.

There are none at present. The six once listed here (`nM` and mate order, same-span alignments, a pair STAR leaves unmapped, the unmapped-read category, STARsolo `GeneFull` multi-mappers and `SJ` UMI filtering) now match STAR on the suite's yeast, chr21, fusion and STARsolo runs.

---

## Adding a new divergence

When a change deliberately diverges from STAR (including adding a non-STAR flag, or choosing STAR's *documented* behaviour over its *actual binary* behaviour where they differ):

1. Add an entry to the appropriate section above using the **What STAR does / What rustar-aligner does / Why / Impact / Source** format.
2. Cite the STAR C++ source you checked, so the divergence can be re-verified.
3. Get maintainer sign-off in the PR — see [CONTRIBUTING.md](CONTRIBUTING.md#divergence-from-star-is-allowed--but-must-be-deliberate-and-flagged).
