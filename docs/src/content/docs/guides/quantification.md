---
title: Gene quantification
description: Generate gene-level read counts and transcriptome-coordinate SAM during alignment.
---

rustar-aligner can produce two flavours of quantification output during alignment, controlled by `--quantMode`. Both require a GTF file.

## Gene-level counts (`GeneCounts`)

Produces a `<prefix>ReadsPerGene.out.tab` file with one row per gene and four columns: gene ID, then read counts for unstranded, forward-stranded, and reverse-stranded protocols. Output is identical to STAR's, so any downstream tool that already consumes STAR `ReadsPerGene.out.tab` works unchanged (e.g. `DESeq2`, `edgeR`, `MultiQC`).

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile gencode.v45.gtf \
  --quantMode GeneCounts \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

### Output format

```
N_unmapped              <count>     <count>     <count>
N_multimapping          <count>     <count>     <count>
N_noFeature             <count>     <count>     <count>
N_ambiguous             <count>     <count>     <count>
ENSG00000000003.15      <count>     <count>     <count>
ENSG00000000005.6       <count>     <count>     <count>
...
```

The first four rows are summary categories (matching STAR / `htseq-count` semantics). Subsequent rows are per-gene counts. Pick the column that matches your library protocol:

| Column | Strand assumption |
|--------|-------------------|
| 2 | unstranded |
| 3 | forward / `htseq-count -s yes` |
| 4 | reverse / `htseq-count -s reverse` |

If you don't know your library's strandedness, all three columns help — you can compare them and infer.

## Transcriptome-coordinate SAM (`TranscriptomeSAM`)

Produces a `<prefix>Aligned.toTranscriptome.out.bam` file with reads mapped to *transcriptome* coordinates (one record per matching transcript) instead of genome coordinates. This is the format consumed by [RSEM](https://github.com/deweylab/RSEM) and similar transcript-level quantifiers.

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile gencode.v45.gtf \
  --quantMode TranscriptomeSAM \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

### Output filtering variants

`--quantTranscriptomeSAMoutput` controls what's allowed in the transcriptome BAM:

| Value | Meaning |
|-------|---------|
| `BanSingleEnd_BanIndels_ExtendSoftclip` | Default. RSEM-compatible: drop unpaired records, drop reads with indels, extend soft-clips into matches. |
| `BanSingleEnd_ExtendSoftclip` | Keep indels, still extend soft-clips. |
| `BanSingleEnd` | Keep indels and soft-clips as-is. |

Pick the variant that matches your downstream tool's expectations. RSEM defaults to `BanSingleEnd_BanIndels_ExtendSoftclip`.

## Combining both

You can request both modes in the same run:

```bash
--quantMode GeneCounts TranscriptomeSAM
```

This emits `ReadsPerGene.out.tab` *and* `Aligned.toTranscriptome.out.bam` in addition to the normal alignment output.

## Index-time vs alignment-time

For best speed, supply `--sjdbGTFfile` at `--runMode genomeGenerate` time and the transcript-level data structures get persisted into the genome directory. Then at alignment time you only need `--quantMode TranscriptomeSAM` (or `GeneCounts`); rustar-aligner reuses the persisted annotations.

If you supply `--sjdbGTFfile` only at alignment time, transcript info is rebuilt on the fly each run. This works but adds startup cost.

## Strandedness summary

For `--quantMode GeneCounts`, all three strand columns are written regardless of library protocol — pick the right one downstream.

For `--quantMode TranscriptomeSAM`, the orientation is inferred from the transcript annotation in the GTF; no explicit strand parameter is needed.
