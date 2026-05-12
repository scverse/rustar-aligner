---
title: CLI parameters
description: Every command-line parameter rustar-aligner accepts, grouped by category.
---

rustar-aligner accepts STAR's `--camelCase` parameter names. Defaults match STAR. This page lists the currently-supported parameters; if a parameter is missing, the binary errors out at startup rather than silently ignoring it.

Run `rustar-aligner --help` for the full machine-generated listing.

## Run mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--runMode` | `alignReads` | `alignReads` or `genomeGenerate`. |
| `--runThreadN` | `1` | Number of threads. |
| `--runRNGseed` | `777` | RNG seed for tie-breaking among equal-scoring alignments. |

## Genome

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genomeDir` | `./GenomeDir` | Path to the genome index directory. |
| `--genomeFastaFiles` | — | One or more FASTA files (required for `genomeGenerate`). |
| `--genomeSAindexNbases` | `14` | Length of the SA pre-indexing string (log2). Lower for small genomes. |
| `--genomeChrBinNbits` | `18` | Log2 of chromosome bin size. |
| `--genomeSAsparseD` | `1` | SA sparsity (higher = less RAM, slower mapping). |

## Read input

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--readFilesIn` | — | Input FASTQ file(s); second file is mate 2 for paired-end (required for `alignReads`). |
| `--readFilesCommand` | — | Decompression command, e.g. `zcat` for `.gz`. |
| `--readMapNumber` | `-1` | Number of reads to map (`-1` = all). |
| `--clip5pNbases` | `0` | Bases to clip from the 5' end of each mate. |
| `--clip3pNbases` | `0` | Bases to clip from the 3' end of each mate. |

## Output: file naming and format

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outFileNamePrefix` | `./` | Prefix (path + filename stem) for all output files. |
| `--outSAMtype` | `SAM` | `SAM`, `BAM Unsorted`, `BAM SortedByCoordinate`, or `None`. |
| `--outBAMcompression` | `1` | BGZF level. `-1`/`0` = uncompressed; `1`–`8` = flate2 levels; `≥9` = max. |
| `--limitBAMsortRAM` | `0` | Max RAM (bytes) for sorted BAM. `0` = unlimited. |
| `--outStd` | `None` | Route primary output to stdout: `None`, `SAM`, `BAM_Unsorted`, `BAM_SortedByCoordinate`. |

## Output: SAM/BAM records

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outSAMstrandField` | `None` | `None` or `intronMotif` (sets XS tag from junction motifs). |
| `--outSAMattributes` | `Standard` | Tags to include: `Standard`, `All`, `None`, or an explicit list (e.g. `NH HI AS NM nM MD`). |
| `--outSAMattrRGline` | `-` | Read group line(s). Multiple blocks separated by a literal `,`. |
| `--outSAMunmapped` | `None` | Unmapped reads in SAM: `None`, `Within`, or `Within KeepPairs`. |
| `--outSAMmapqUnique` | `255` | MAPQ value for uniquely-mapping reads. |
| `--outSAMmultNmax` | `-1` | Max alignments per read in SAM (`-1` = all up to `outFilterMultimapNmax`). |

## Output: read filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outFilterType` | `Normal` | `Normal` or `BySJout` (re-filter by discovered SJ pass). |
| `--outFilterMultimapNmax` | `10` | Max number of multi-mapping loci. Reads exceeding this are unmapped (`MultiMapTooMany`). |
| `--outFilterMultimapScoreRange` | `1` | Score range for keeping multi-mappers within best score. |
| `--outFilterMismatchNmax` | `10` | Max mismatches per pair. |
| `--outFilterMismatchNoverLmax` | `0.3` | Max ratio of mismatches to mapped length. |
| `--outFilterScoreMin` | `0` | Min absolute alignment score. |
| `--outFilterScoreMinOverLread` | `0.66` | Min alignment score normalized to read length. |
| `--outFilterMatchNmin` | `0` | Min absolute matched bases. |
| `--outFilterMatchNminOverLread` | `0.66` | Min matched bases normalized to read length. |
| `--outFilterIntronMotifs` | `None` | `None`, `RemoveNoncanonical`, or `RemoveNoncanonicalUnannotated`. |
| `--outFilterIntronStrands` | `RemoveInconsistentStrands` | `None` or `RemoveInconsistentStrands`. |

## Output: unmapped reads

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outReadsUnmapped` | `None` | `None` or `Fastx` (write `Unmapped.out.mate1`/`mate2`). |

## Output: splice junctions

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outSJfilterOverhangMin` | `30 12 12 12` | Min overhang per motif `[noncan, GT/AG, GC/AG, AT/AC]`. |
| `--outSJfilterCountUniqueMin` | `3 1 1 1` | Min unique-mapping reads per motif. |
| `--outSJfilterCountTotalMin` | `3 1 1 1` | Min total reads per motif. |
| `--outSJfilterDistToOtherSJmin` | `10 0 5 10` | Min distance to other SJs per motif. |
| `--outSJfilterIntronMaxVsReadN` | `50000 100000 200000` | Max intron length per supporting-read count tier. |

## Alignment scoring

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--alignIntronMin` | `21` | Min intron size (smaller gaps are deletions). |
| `--alignIntronMax` | `0` | Max intron size; `0` = auto from genome / win params. |
| `--alignMatesGapMax` | `0` | Max genomic distance between PE mates; `0` = auto. |
| `--alignSplicedMateMapLmin` | `0` | Min mapped length for spliced PE mate (absolute). |
| `--alignSplicedMateMapLminOverLmate` | `0.66` | Min mapped length for spliced PE mate (fraction). |
| `--alignSJoverhangMin` | `5` | Min overhang for novel splice junctions. |
| `--alignSJDBoverhangMin` | `3` | Min overhang for annotated junctions. |
| `--alignSJstitchMismatchNmax` | `0 -1 0 0` | Max mismatches for SJ stitching `[noncan, GC/AG, AT/AC, noncan]`. |

## Scoring penalties

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--scoreGap` | `0` | Canonical splice junction penalty. |
| `--scoreGapNoncan` | `-8` | Non-canonical junction penalty. |
| `--scoreGapGCAG` | `-4` | GC/AG junction penalty. |
| `--scoreGapATAC` | `-8` | AT/AC junction penalty. |
| `--scoreDelOpen` | `-2` | Deletion open penalty. |
| `--scoreDelBase` | `-2` | Deletion extension penalty per base. |
| `--scoreInsOpen` | `-2` | Insertion open penalty. |
| `--scoreInsBase` | `-2` | Insertion extension penalty per base. |
| `--scoreStitchSJshift` | `1` | Max score reduction for SJ stitching shift. |
| `--scoreGenomicLengthLog2scale` | `-0.25` | Log-scaled bonus per `log2(genomicLength)`. |

## Seeds and windows

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--winReadCoverageRelativeMin` | `0.5` | Min read coverage for an alignment window (fraction). |
| `--winBinNbits` | `16` | Log2 of window bin size for seed clustering. |
| `--winAnchorDistNbins` | `9` | Max bins for seed anchor distance. |
| `--winFlankNbins` | `4` | Bins to extend each window by on each side. |
| `--winAnchorMultimapNmax` | `50` | Max loci an anchor can map to. |
| `--seedMultimapNmax` | `10000` | Max loci a seed can map to. |
| `--seedPerReadNmax` | `1000` | Max seeds per read. |
| `--seedPerWindowNmax` | `50` | Max seeds per window. |
| `--seedSearchStartLmax` | `50` | Max distance between seed search start positions. |
| `--seedSearchStartLmaxOverLread` | `1.0` | `seedSearchStartLmax` normalised by read length. |
| `--seedSearchLmax` | `0` | Max seed length; `0` = unlimited. |
| `--seedMapMin` | `5` | Min mappable length for seed search termination. |
| `--alignWindowsPerReadNmax` | `10000` | Max alignment windows per read. |
| `--alignTranscriptsPerWindowNmax` | `100` | Max transcripts per window. |

## Splice junction database

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sjdbGTFfile` | — | GTF file with exon annotations. |
| `--sjdbGTFchrPrefix` | `""` | Prefix to add to chromosome names from GTF (e.g. `chr`). |
| `--sjdbGTFfeatureExon` | `exon` | GTF feature type to use as exon. |
| `--sjdbGTFtagExonParentTranscript` | `transcript_id` | GTF attribute for transcript ID. |
| `--sjdbGTFtagExonParentGene` | `gene_id` | GTF attribute for gene ID. |
| `--sjdbOverhang` | `100` | Overhang length around junctions in the index. Set to `read_length - 1`. |
| `--sjdbScore` | `2` | Extra score for alignments crossing annotated junctions. |

## Quantification

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quantMode` | — | `GeneCounts` and/or `TranscriptomeSAM`, space-separated. |
| `--quantTranscriptomeSAMoutput` | `BanSingleEnd_BanIndels_ExtendSoftclip` | Variant for transcriptome BAM: `BanSingleEnd`, `BanSingleEnd_ExtendSoftclip`, or the default RSEM-compatible form. |

## Two-pass mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--twopassMode` | `None` | `None` or `Basic`. |
| `--twopass1readsN` | `-1` | Reads to use in pass 1 (`-1` = all). |

## Chimeric detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--chimSegmentMin` | `0` | Min chimeric segment length. `0` disables chimeric detection. |
| `--chimScoreMin` | `0` | Min total chimeric alignment score. |
| `--chimScoreDropMax` | `20` | Max drop in chimeric score vs read length. |
| `--chimScoreSeparation` | `10` | Min score separation for unique chimeric. |
| `--chimMainSegmentMultNmax` | `10` | Max multimapping for main chimeric segment. |
| `--chimSegmentReadGapMax` | `0` | Max read-space gap between chimeric segments. |
| `--chimJunctionOverhangMin` | `20` | Min overhang at chimeric junction. |
| `--chimScoreJunctionNonGTAG` | `-1` | Score penalty for non-GT/AG chimeric junctions. |
| `--chimOutType` | `Junctions` | `Junctions`, `WithinBAM`, or both space-separated. |

## Debug

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--readNameFilter` | `""` | If set, only emit detailed alignment logs for reads with this name. |
