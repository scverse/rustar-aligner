---
title: Output files
description: Every file rustar-aligner writes during alignment, with format and contents.
---

rustar-aligner produces the same set of output files as STAR. All files are written into the directory implied by `--outFileNamePrefix`, with names that match STAR's exactly.

For the examples below, assume `--outFileNamePrefix sample_`. Replace `sample_` with whatever you've configured.

## Alignment output

### `sample_Aligned.out.sam`

Default output. SAM-format alignment file with header and records. Written when `--outSAMtype SAM` (the default).

### `sample_Aligned.out.bam`

Unsorted BAM equivalent of `Aligned.out.sam`. Written when `--outSAMtype BAM Unsorted`.

### `sample_Aligned.sortedByCoord.out.bam`

Coordinate-sorted BAM. Written when `--outSAMtype BAM SortedByCoordinate`. The sort happens in-memory by default; cap the RAM with `--limitBAMsortRAM` if needed.

### `sample_Aligned.toTranscriptome.out.bam`

Transcriptome-coordinate BAM. Written when `--quantMode TranscriptomeSAM` is set. Each record's reference is a transcript ID rather than a chromosome; one record is emitted per transcript that the read aligns within.

## Log files

### `sample_Log.final.out`

The summary statistics file — the most useful output for monitoring runs and feeding into MultiQC. Contains counts and rates for: input reads, uniquely mapped reads, reads mapped to multiple loci, reads mapped to too many loci, unmapped reads (with reasons), splice junctions, mismatch rates, deletion/insertion rates, and run timing.

Format is identical to STAR's, so MultiQC parses it without modification.

### `sample_Log.out`

The verbose run log. Contains parameter values used, genome loading info, per-thread progress, and warnings. Useful when something looks wrong — search for `WARNING` or `ERROR`.

### `sample_Log.progress.out`

Per-chunk progress lines emitted during alignment. Useful for monitoring long runs.

## Splice junctions

### `sample_SJ.out.tab`

Tab-separated table of splice junctions discovered during alignment. STAR-compatible 9-column format:

| # | Column | Meaning |
|---|--------|---------|
| 1 | `chr` | Chromosome |
| 2 | `start` | Intron start (1-based) |
| 3 | `end` | Intron end (1-based) |
| 4 | `strand` | `0` = undefined, `1` = `+`, `2` = `-` |
| 5 | `motif` | `0` = non-canonical, `1` = GT/AG, `2` = CT/AC, `3` = GC/AG, `4` = CT/GC, `5` = AT/AC, `6` = GT/AT |
| 6 | `annotated` | `0` = novel, `1` = present in GTF / sjdb |
| 7 | `unique_reads` | Number of uniquely-mapping reads supporting the junction |
| 8 | `multi_reads` | Number of multi-mapping reads supporting the junction |
| 9 | `max_overhang` | Max overhang of any supporting read |

Filtered using the `--outSJfilter*` parameters (see [CLI reference](/rustar-aligner/reference/cli-parameters/)).

## Chimeric output

Written when `--chimSegmentMin > 0`. Format depends on `--chimOutType`.

### `sample_Chimeric.out.junction`

14-column tab-separated table when `--chimOutType` includes `Junctions` (default). One row per chimeric junction. Format is STAR-compatible. See the [chimeric guide](/rustar-aligner/guides/chimeric/) for the full column breakdown.

### Chimeric records in the primary BAM

When `--chimOutType` includes `WithinBAM`, the chimeric segments are embedded as supplementary alignment records (FLAG `0x800`) in the main BAM output, with `SA` tags linking the donor and acceptor halves. Tools like Arriba and STAR-Fusion know how to read either format.

## Quantification

### `sample_ReadsPerGene.out.tab`

Per-gene read counts. Written when `--quantMode GeneCounts` is set. Four-column format:

```
gene_id    unstranded    forward_stranded    reverse_stranded
```

The first four rows are summary categories: `N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`. Subsequent rows are per-gene counts. Pick the column matching your library's strandedness — see the [quantification guide](/rustar-aligner/guides/quantification/).

## Unmapped reads

### `sample_Unmapped.out.mate1` / `sample_Unmapped.out.mate2`

Written when `--outReadsUnmapped Fastx` is set. FASTQ output of every read that didn't map (and reads that mapped to too many loci). Mate 2 is only present for paired-end input. For PE half-mapped pairs (one mate mapped, one not), both mates are written so the pair can be re-aligned later.

## Two-pass mode

In two-pass mode (`--twopassMode Basic`), pass 1 produces an internal `_STARpass1/SJ.out.tab` that's used to seed pass 2. The final `sample_SJ.out.tab` reflects pass 2's discoveries; pass 1's intermediate file is not retained by default.

## Standard output

When `--outStd` is set to `SAM`, `BAM_Unsorted`, or `BAM_SortedByCoordinate`, the corresponding alignment output is sent to stdout instead of a file. Useful for piping directly into samtools, sambamba, or any other downstream tool:

```bash
rustar-aligner --outStd BAM_Unsorted ... | samtools sort -@ 4 -o sample.sorted.bam
```

The other output files (`Log.final.out`, `SJ.out.tab`, etc.) are still written to disk normally.
