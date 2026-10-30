---
title: Output files
description: Every file rustar-aligner writes during alignment, with format and contents.
---

rustar-aligner produces the same set of output files as STAR. All files are written into the directory implied by `--outFileNamePrefix`, with names that match STAR's exactly.

For the examples below, assume `--outFileNamePrefix sample_`. Replace `sample_` with whatever you've configured.

## Alignment output

### `sample_Aligned.out.sam`

Default output. SAM-format alignment file with header and records. Written when `--outSAMtype SAM` (the default).

The header follows STAR's layout: `@HD VN:1.4` (or the `--outSAMheaderHD` line as given), the `@SQ` lines, the `--outSAMheaderPG` line if any, the program's own `@PG`, the lines of `--outSAMheaderCommentFile` (blank lines skipped), the `@RG` lines, and `@CO user command line: ...`. A coordinate-sorted BAM appends `SO:coordinate` to the `@HD` line. The transcriptome BAM's header has only `@SQ` and `@RG` lines.

### `sample_Aligned.out.bam`

Unsorted BAM equivalent of `Aligned.out.sam`. Written when `--outSAMtype BAM Unsorted` (also together with the sorted BAM: `--outSAMtype BAM Unsorted SortedByCoordinate`).

### `sample_Aligned.sortedByCoord.out.bam`

Coordinate-sorted BAM. Written when `--outSAMtype BAM SortedByCoordinate`. The sort happens in-memory by default; cap the RAM with `--limitBAMsortRAM` if needed.

### `sample_Aligned.toTranscriptome.out.bam`

Transcriptome-coordinate BAM. Written when `--quantMode TranscriptomeSAM` is set. Each record's reference is a transcript ID rather than a chromosome; one record (one per mate for a pair) is emitted per transcript that the read aligns within, in STAR's order. MAPQ comes from the number of transcript alignments, and one of them is marked primary at random, as STAR does (see the migration guide's tie-breaking notes). There is no `AS` tag, as in STAR; `--quantTranscriptomeBAMattrAS Yes` (rustar-aligner extension) adds `AS`. With `--quantTranscriptomeUnspliced Intron|PreMRNA` (rustar-aligner extension) the references also include one `<gene_id>-I` unspliced target per gene, and `sample_Aligned.toTranscriptome.targets.tsv` (target, gene, `spliced`/`unspliced` status, length) is written next to it, plus `sample_Aligned.toTranscriptome.unspliced.fa` with `--quantTranscriptomeUnsplicedFasta Yes`. See [Bulk total RNA-seq](/rustar-aligner/guides/bulk-total-rna/).

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

When `--chimOutType` includes `WithinBAM`, a chimeric read is written into the main BAM in place of its normal alignment, as STAR does: one segment as a normal record and the other as a supplementary record (FLAG `0x800`, hard-clipped unless `SoftClip`), with `SA` tags linking the two. Paired reads whose mates lie on either side of the junction are written as two normal records. See the [chimeric guide](/rustar-aligner/guides/chimeric/).

## Quantification

### `sample_ReadsPerGene.out.tab`

Per-gene read counts. Written when `--quantMode GeneCounts` is set. Four-column format:

```
gene_id    unstranded    forward_stranded    reverse_stranded
```

The first four rows are summary categories: `N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`. Subsequent rows are per-gene counts. Pick the column matching your library's strandedness — see the [quantification guide](/rustar-aligner/guides/quantification/).

### `sample_ReadsPerGeneSplicing.out.tab` / `sample_ReadsPerGeneSplicing.summary.tsv`

Written when `--quantGeneSplicing Yes` is set (rustar-aligner extension, not in STAR). The table has a header line and one line per gene with spliced, unspliced and ambiguous counts for the unstranded, forward and reverse strand conventions (nine count columns). The summary gives read accounting and the spliced / unspliced / ambiguous fractions per strand convention. See [Bulk total RNA-seq](/rustar-aligner/guides/bulk-total-rna/).

## Unmapped reads

### `sample_Unmapped.out.mate1` / `sample_Unmapped.out.mate2`

Written when `--outReadsUnmapped Fastx` is set. FASTQ output of every read that didn't map (and reads that mapped to too many loci). Mate 2 is only present for paired-end input. For PE half-mapped pairs (one mate mapped, one not), both mates are written so the pair can be re-aligned later.

## Two-pass mode

In two-pass mode (`--twopassMode Basic`), pass 1 writes `sample__STARpass1/SJ.out.tab` and `sample__STARpass1/Log.final.out`; the junctions in that `SJ.out.tab` are inserted into the genome for pass 2. The final `sample_SJ.out.tab` comes from pass 2.

Whenever junctions are inserted at mapping (two-pass, or `--sjdbGTFfile` / `--sjdbFileChrStartEnd` at mapping), `sample__STARgenome/` gets `sjdbInfo.txt` and `sjdbList.out.tab`, the transcript tables when a GTF was given, and with `--sjdbInsertSave All` the extended Genome, SA and SAindex.

## Standard output

When `--outStd` is set to `SAM`, `BAM_Unsorted`, or `BAM_SortedByCoordinate`, the corresponding alignment output is sent to stdout instead of a file. Useful for piping directly into samtools, sambamba, or any other downstream tool:

```bash
rustar-aligner --outStd BAM_Unsorted ... | samtools sort -@ 4 -o sample.sorted.bam
```

The other output files (`Log.final.out`, `SJ.out.tab`, etc.) are still written to disk normally.
