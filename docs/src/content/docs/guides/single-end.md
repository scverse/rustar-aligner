---
title: Single-end alignment
description: Align single-end RNA-seq reads against a rustar-aligner / STAR genome index.
---

Single-end (SE) alignment maps each read in a single FASTQ file independently — there's no mate pair to constrain the alignment.

## Basic command

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outFileNamePrefix sample_
```

## Strand information

For stranded RNA-seq libraries, add the strand field to the output so that downstream tools (Cufflinks, StringTie, etc.) can pick up the strand assignment from junction motifs:

```bash
--outSAMstrandField intronMotif
```

## Gzipped input

Use `--readFilesCommand zcat` to pipe through decompression:

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

On macOS, use `gzcat` instead of `zcat`.

## Threading

Use `--runThreadN` to parallelise the alignment:

```bash
--runThreadN 16
```

The alignment phase scales nearly linearly with thread count up to the IO limit of the storage holding the input/output files.

## Filtering

The defaults match STAR. Two parameters are particularly useful to know:

- `--outFilterMultimapNmax` (default `10`) — max number of alignments per read. A read mapping to more loci is discarded as `MultiMapTooMany`.
- `--outFilterMismatchNoverLmax` (default `0.3`) — max ratio of mismatches to mapped length. Tightening this to `0.04` is a common choice for human RNA-seq.
- `--outFilterScoreMinOverLread` (default `0.66`) — min alignment score relative to read length. Reads scoring below this threshold are unmapped.

See the [CLI parameters reference](/rustar-aligner/reference/cli-parameters/) for the full list.

## Output limits

Multi-mapping reads emit one SAM record per locus by default (up to `--outFilterMultimapNmax`). To cap the number of secondary alignments written:

```bash
--outSAMmultNmax 5
```

## Reproducibility

When two alignments tie on score, rustar-aligner uses a seeded RNG (`StdRng`) to pick the primary. Set the seed for reproducible runs:

```bash
--runRNGseed 42
```

The default seed is `777`.

## What's next

- [Paired-end alignment](/rustar-aligner/guides/paired-end/) for paired libraries
- [Two-pass mode](/rustar-aligner/guides/two-pass/) for novel-junction recovery
- [Output files reference](/rustar-aligner/reference/output-files/) for what each output file contains
