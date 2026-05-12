---
title: Quick start
description: Generate a genome index, align reads, and inspect the output.
---

This walks through the minimum needed to run rustar-aligner end-to-end. If you already use STAR, the commands will look familiar — they are.

## 1. Generate a genome index

You need to do this once per reference genome.

```bash
rustar-aligner --runMode genomeGenerate \
  --genomeDir /path/to/genome_index \
  --genomeFastaFiles /path/to/genome.fa
```

For a human genome you'll typically also pass a GTF file and tune the SA index parameter:

```bash
rustar-aligner --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir /path/to/genome_index \
  --genomeFastaFiles /path/to/GRCh38.fa \
  --sjdbGTFfile /path/to/gencode.gtf \
  --sjdbOverhang 100
```

See the [genome index guide](/rustar-aligner/guides/genome-index/) for the full set of relevant parameters.

## 2. Align reads (single-end)

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outSAMstrandField intronMotif \
  --outFileNamePrefix sample_
```

This writes `sample_Aligned.out.sam` plus `sample_Log.final.out`, `sample_SJ.out.tab`, and other STAR-style output files into the current directory. Override the destination with `--outFileNamePrefix /path/to/sample_`.

## 3. Align reads (paired-end)

Pass two files to `--readFilesIn`:

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq reads_2.fq \
  --outSAMtype SAM \
  --outFileNamePrefix sample_
```

## 4. Get a sorted BAM directly

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

The output is `sample_Aligned.sortedByCoord.out.bam`. Pass `BAM Unsorted` instead for an unsorted BAM.

## 5. Gzipped FASTQ input

Use `--readFilesCommand zcat` to decompress on the fly:

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

## 6. Inspect the output

After a successful run you get the standard STAR file set:

| File | Contents |
|------|----------|
| `*_Aligned.out.sam` / `*_Aligned.sortedByCoord.out.bam` | The alignments |
| `*_Log.final.out` | Summary statistics (MultiQC-compatible) |
| `*_Log.out` | Verbose run log |
| `*_Log.progress.out` | Per-chunk progress |
| `*_SJ.out.tab` | Splice junctions discovered during alignment |

See the [output files reference](/rustar-aligner/reference/output-files/) for a full breakdown of each file.

## Common next steps

- **Two-pass alignment** for better novel-junction recovery: see the [two-pass guide](/rustar-aligner/guides/two-pass/).
- **Gene-level counts** for differential expression: see the [quantification guide](/rustar-aligner/guides/quantification/).
- **Chimeric detection** for fusion calling: see the [chimeric guide](/rustar-aligner/guides/chimeric/).
