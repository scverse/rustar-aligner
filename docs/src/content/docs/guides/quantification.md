---
title: Gene quantification
description: Generate gene-level read counts and transcriptome-coordinate SAM during alignment.
---

rustar-aligner can produce two flavours of quantification output during alignment, controlled by `--quantMode`. Both need a gene annotation: the tables of an index built with `--sjdbGTFfile`, or `--sjdbGTFfile` at alignment time. With neither, the run stops with STAR's "could not open ... geneInfo.tab" error.

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

## Transcriptome RAD for salmon (`--quantTranscriptomeRAD`)

A rustar-aligner extension. `--quantTranscriptomeRAD Yes` writes the same transcriptome alignments to `Aligned.toTranscriptome.out.rad`, for `salmon quant --rad` (salmon 2.x; tested with 2.8.0). The file has no read names, sequences or qualities: about 8 times smaller than the BAM on paired yeast reads. The BAM is written only if `--quantMode TranscriptomeSAM` is given as well.

```bash
rustar-aligner --genomeDir index --readFilesIn r1.fq.gz r2.fq.gz --quantTranscriptomeRAD Yes
salmon quant --rad Aligned.toTranscriptome.out.rad -l A -o quant
```

`--quantTranscriptomeSAMoutput` and `--quantTranscriptomeUnspliced` apply as they do to the BAM. Each placement carries its alignment score, which the BAM does not, and salmon weights placements by it; `--scoreExp 0` turns that off and gives the same `quant.sf` as `salmon quant -a` on the BAM. For single-end reads give the library type (`-l U`, `-l SF` or `-l SR`), and bias correction needs `-t transcripts.fa`. `--errorModel` has no effect, as the file has no read bases. salmon cannot read this file from a pipe. It is not available with `--soloType`. See [DIVERGENCE.md §1.9](https://github.com/scverse/rustar-aligner/blob/main/DIVERGENCE.md).

## Combining both

You can request both modes in the same run:

```bash
--quantMode GeneCounts TranscriptomeSAM
```

This emits `ReadsPerGene.out.tab` *and* `Aligned.toTranscriptome.out.bam` in addition to the normal alignment output.

## Total RNA (ribo-depleted) libraries

For libraries with a large pre-mRNA content, `--quantTranscriptomeUnspliced` adds one unspliced `<gene_id>-I` target per gene to the transcriptome BAM (so pre-mRNA reads stop landing on retained-intron isoforms), `--outSAMsplicingStatus Yes` tags genomic alignments with their splicing status, and `--quantGeneSplicing Yes` counts spliced / unspliced / ambiguous reads per gene. All three are rustar-aligner extensions; see [Bulk total RNA-seq](/rustar-aligner/guides/bulk-total-rna/).

## Index-time vs alignment-time

For best speed, supply `--sjdbGTFfile` at `--runMode genomeGenerate` time and the transcript-level data structures get persisted into the genome directory. Then at alignment time you only need `--quantMode TranscriptomeSAM` (or `GeneCounts`, or STARsolo's gene features); rustar-aligner reuses the persisted annotations.

If you supply `--sjdbGTFfile` at alignment time, its junctions are inserted into the genome and its transcript and gene info is used, even if the index has its own, as in STAR. The transcripts are put in STAR's order, so the transcriptome BAM's `@SQ` lines are the same as from an index built with that GTF. This is done on every run, so it adds startup time.

## Strandedness summary

For `--quantMode GeneCounts`, all three strand columns are written regardless of library protocol — pick the right one downstream.

For `--quantMode TranscriptomeSAM`, the orientation is inferred from the transcript annotation in the GTF; no explicit strand parameter is needed.
