---
title: Two-pass mode
description: Discover novel splice junctions in pass 1, then re-align with them in pass 2.
---

Two-pass alignment improves recovery of reads spanning novel splice junctions — those not present in the GTF (or when no GTF is supplied). rustar-aligner runs both passes in a single command, matching STAR's `--twopassMode Basic`.

## How it works

1. **Pass 1** aligns all reads (or a subset, see `--twopass1readsN`) against the original index. Splice junctions discovered during this pass are filtered and collected.
2. The discovered junctions are inserted into a per-run on-the-fly junction database.
3. **Pass 2** re-aligns all reads against the augmented index. Reads that previously soft-clipped over an unannotated junction now stitch through it.

The pass-1 junctions are filtered the same way as STAR: by motif (canonical GT/AG vs non-canonical), unique-mapping support, distance to other junctions, and intron length vs read count. Defaults are STAR-compatible.

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

You can use two-pass mode together with `--sjdbGTFfile`. Annotated junctions from the GTF are loaded into the index alongside the pass-1 discoveries; a junction in both gets the higher-confidence treatment.

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

Two-pass mode writes the same set of files as a single-pass run. The `SJ.out.tab` reflects the final pass-2 set of junctions used during alignment.
