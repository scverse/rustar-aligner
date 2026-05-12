---
title: Paired-end alignment
description: Align paired-end RNA-seq reads with mate rescue.
---

Paired-end (PE) alignment maps both mates of a fragment together, using the expected insert size and pair geometry to constrain the alignment. rustar-aligner handles mate rescue (recovering one mate when only the other anchors), proper-pair flagging, and reports `NH` / `HI` tags consistent with STAR's PE pipeline.

## Basic command

Pass both FASTQ files to `--readFilesIn`. **Mate 1 first, mate 2 second** — order matters.

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq reads_2.fq \
  --outSAMtype SAM \
  --outFileNamePrefix sample_
```

## Gzipped input

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

## Insert size

By default rustar-aligner computes the maximum mate gap from the genome size and `winBinNbits` / `winAnchorDistNbins`. To set an explicit cap:

```bash
--alignMatesGapMax 1000000
```

This is the maximum genomic distance between the leftmost ends of the two mates. The default is suitable for most RNA-seq libraries.

## Spliced mate constraints

Two parameters control the minimum mapped length of a spliced mate. They were both introduced because, for poorly-anchored mates, STAR (and rustar-aligner) can otherwise emit very short, low-confidence spliced alignments:

- `--alignSplicedMateMapLmin` (absolute, default `0` = off) — minimum number of bases that must be mapped on a spliced mate.
- `--alignSplicedMateMapLminOverLmate` (default `0.66`) — same, but expressed as a fraction of read length.

The default config means a spliced mate must map at least 66% of its bases, with no absolute floor.

## Unmapped mates

By default unmapped mates aren't written into the SAM/BAM output. To include them in the alignment file:

```bash
--outSAMunmapped Within
```

Or to also keep pairs together when only one mate is mapped:

```bash
--outSAMunmapped Within KeepPairs
```

Or to write unmapped reads to separate FASTQ files:

```bash
--outReadsUnmapped Fastx
```

This produces `<prefix>Unmapped.out.mate1` and `<prefix>Unmapped.out.mate2`.

## Half-mapped pairs

If only one mate has a confident alignment, rustar-aligner still emits a record for the other mate (flagged unmapped) so the pair stays together in the SAM/BAM. On the yeast benchmark (10,000 paired reads) rustar-aligner produces **8,390 fully-mapped pairs and zero half-mapped pairs**, matching STAR exactly.

## Read groups

If you want `@RG` headers on the output (e.g. for downstream tools that group by sample), supply them with `--outSAMattrRGline`:

```bash
--outSAMattrRGline ID:sample1 SM:sample1 LB:lib1 PL:ILLUMINA
```

For multi-lane runs, separate RG blocks with a literal `,`:

```bash
--outSAMattrRGline ID:lane1 SM:sample1 , ID:lane2 SM:sample1
```

## What's next

- [Two-pass mode](/rustar-aligner/guides/two-pass/) for novel-junction recovery
- [Chimeric detection](/rustar-aligner/guides/chimeric/) for fusion calling on PE data
- [STAR compatibility report](/rustar-aligner/reference/star-compatibility/) for the PE faithfulness numbers
