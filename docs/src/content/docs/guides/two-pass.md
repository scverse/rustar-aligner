---
title: Two-pass mode
description: Discover novel splice junctions in pass 1, then re-align with them in pass 2.
---

Two-pass alignment improves recovery of reads spanning novel splice junctions — those not present in the GTF (or when no GTF is supplied). rustar-aligner runs both passes in a single command, matching STAR's `--twopassMode Basic`.

## How it works

1. **Pass 1** aligns all reads (or a subset, see `--twopass1readsN`) against the original index, with output, chimeric detection and quantification switched off. Its junctions go through the `SJ.out.tab` filters and are written to `<prefix>_STARpass1/SJ.out.tab`.
2. Every junction in that file is inserted into the genome, as STAR does: the junction flanks are added to the genome and the suffix array and SAindex are updated, so pass-2 seeds can cross the new junctions.
3. **Pass 2** re-aligns all reads against the extended genome. In pass 2 the new junctions count as annotated, as in STAR (`sjdbScore`, column 6 of `SJ.out.tab`).

## Basic command

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

## Combining with a GTF

You can use two-pass mode together with `--sjdbGTFfile`. The GTF's junctions are inserted before pass 1; the pass-1 junctions are then added to them. When the same junction comes from both, the genome's and the GTF's version is kept, as STAR does.

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile gencode.v45.gtf \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

## Limiting pass 1

For very large datasets you can limit pass 1 to a subset of reads:

```bash
--twopass1readsN 1000000
```

This processes the first million reads in pass 1, then runs pass 2 on the full input. The default (`-1`) uses all reads in both passes.

## Junction filter parameters

The defaults match STAR. The most relevant filters:

- `--outSJfilterOverhangMin` — min overhang per motif category `[noncanon, GT/AG, GC/AG, AT/AC]`. Default `30 12 12 12`.
- `--outSJfilterCountUniqueMin` — min unique-mapping reads per motif. Default `3 1 1 1`.
- `--outSJfilterCountTotalMin` — min total reads per motif. Default `3 1 1 1`.
- `--outSJfilterIntronMaxVsReadN` — junction is filtered when intron length exceeds the threshold for its read count. Default `50000 100000 200000` (1, 2, 3+ supporting reads).

See the [CLI parameters reference](/rustar-aligner/reference/cli-parameters/) for the full list.

## When to use it

Use two-pass mode when:

- Working with a species or sample where the GTF annotation is incomplete
- You expect novel junctions (cancer transcriptomes, non-model organisms, single-cell)
- You want to recover the maximum number of cross-junction reads

Skip it when:

- Speed matters more than novel junction sensitivity
- You have a high-quality, complete GTF and don't need novel discovery

## Output

Two-pass mode writes the same set of files as a single-pass run, plus `<prefix>_STARpass1/` (pass-1 `SJ.out.tab` and `Log.final.out`) and `<prefix>_STARgenome/` (the inserted junctions; see `--sjdbInsertSave`). The main `SJ.out.tab` comes from pass 2.
