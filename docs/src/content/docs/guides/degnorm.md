---
title: Degradation normalization (DegNorm)
description: Capture per-gene coverage during alignment, then correct read counts for sample- and gene-specific transcript degradation.
---

RNA degrades unevenly. Long transcripts lose more material than short ones, and
the loss differs from sample to sample, so a single global size factor cannot
undo it. [DegNorm](https://nustatbioinfo.github.io/DegNorm/) (Xiong et al.,
*Genome Biology* 2019) models this directly: it fits each gene's coverage matrix
across samples with a rank-one over-approximation, splitting a shared coverage
envelope from per-sample abundances, and reports a Degradation Index (DI) score
per gene per sample plus degradation-adjusted counts.

rustar-aligner implements this in two phases:

1. **During alignment**, `--quantMode GeneCoverage` records per-gene, per-exonic
   base coverage and writes `GeneCoverage.out.bin`.
2. **After alignment**, `--runMode degNorm` merges several samples' coverage
   files and fits the model.

:::note[This is not a STAR feature]
Both flags are rustar-aligner extensions with no STAR counterpart, off by
default, and documented in `DIVERGENCE.md`. Alignment output is unchanged when
they are enabled.
:::

## Why two phases

The DI score is defined *across* samples. With one sample the rank-one fit is
exact and every DI is zero, so the model needs at least two coverage files and
cannot run inside a single alignment. What a single run can do, essentially for
free, is the expensive part DegNorm normally pays for by writing a sorted BAM,
indexing it, and re-reading it: computing the coverage curves themselves.

## Phase 1: capture coverage

```bash
for sample in ctrl1 ctrl2 treat1 treat2; do
  rustar-aligner \
    --genomeDir /path/to/index \
    --readFilesIn ${sample}.fastq.gz \
    --readFilesCommand zcat \
    --sjdbGTFfile annotations.gtf \
    --quantMode GeneCounts GeneCoverage \
    --outFileNamePrefix out/${sample}_
done
```

This writes `out/${sample}_GeneCoverage.out.bin` alongside the usual output.

Coverage is accumulated from the same reads that feed `ReadsPerGene.out.tab`
column 1: uniquely mapped, and assigned to exactly one gene. Multimappers and
gene-ambiguous reads are skipped, matching DegNorm's default. For paired-end
data the fragment is one observation: mate blocks are merged first, so an
overlapping pair contributes one to each base it covers, not two.

| Flag | Default | Meaning |
|---|---|---|
| `--quantMode GeneCoverage` | off | Enable coverage capture (needs `--sjdbGTFfile`) |
| `--degNormSampleId` | basename of `--outFileNamePrefix` | Sample name stored in the file and used as a column header later |

Memory: about 4 bytes per exonic base of the annotation, roughly 280 MB for
human GENCODE. Yeast or a chromosome-scale index costs a few megabytes.

## Phase 2: fit the model

```bash
rustar-aligner \
  --runMode degNorm \
  --degNormCoverageFiles out/ctrl1_GeneCoverage.out.bin \
                         out/ctrl2_GeneCoverage.out.bin \
                         out/treat1_GeneCoverage.out.bin \
                         out/treat2_GeneCoverage.out.bin \
  --outFileNamePrefix out/
```

This phase touches no genome index and no read file: it is pure CPU work over
the coverage matrices, parallelised over genes. All samples must have been
aligned against the same GTF; a mismatching gene model is a fatal error naming
the first gene that differs.

| Flag | Default | DegNorm CLI equivalent |
|---|---|---|
| `--degNormCoverageFiles` | (required, at least 2) | `--bam-files` |
| `--degNormIter` | 5 | `--iter` |
| `--degNormNmfIter` | 100 | `--nmf-iter` |
| `--degNormDownsampleRate` | 1 (no downsampling) | `--downsample-rate` |
| `--degNormMinimaxCoverage` | 0 | `--minimax-coverage` |
| `--degNormSkipBaselineSelection` | off | `--skip-baseline-selection` |
| `--degNormBins` | 20 | (baseline-selection bins) |
| `--degNormMinHighCoverage` | 50 | (minimum high-coverage positions) |

`--runRNGseed` sets the offset used by systematic downsampling, so a run is
reproducible regardless of thread count.

## Output

Everything lands in `<outFileNamePrefix>DegNorm.out/`:

| File | Contents |
|---|---|
| `DegradationIndex.tab` | DI per gene (rows) per sample (columns). 0 = no detected degradation, capped at 0.9 |
| `AdjustedCounts.tab` | Depth-normalised counts divided by `1 - DI` |
| `RawCounts.tab` | The unique counts that entered the model, for provenance |
| `ScaleFactors.tab` | Final per-sample sequencing-depth factors |
| `Summary.txt` | Parameters, gene counts, genes sent through baseline selection per iteration, median DI per sample |

`AdjustedCounts.tab` is the matrix to hand to a downstream differential
expression tool.

## How the fit works

For each gene, the coverage matrix `F` (samples x exonic positions) is scaled by
the current depth factors, then fitted by a rank-one over-approximation: dual
ascent lifts the fit above the observed coverage, so the envelope represents the
undegraded curve rather than the average one. The DI is the area between the
scaled envelope and the observed coverage, divided by the area under the
envelope.

Baseline selection refines this. When a gene shows degradation, the transcript
is split into bins and the bin with the worst relative residual is dropped
repeatedly, until what remains is a region where every sample agrees. The
envelope re-estimated from that baseline region gives a cleaner DI for the whole
transcript. Genes that never enter baseline selection inherit the sample-average
DI, as upstream does.

The outer loop then folds the degradation correction back into the depth
factors, and repeats `--degNormIter` times.

## Limitations

- At least two samples are required.
- Coverage must come from phase 1; third-party BAMs are not read.
- Coverage-curve plots and the estimated coverage matrices that upstream stores
  as `.pkl` files are not produced.
- No MPI or warm-start directories.
