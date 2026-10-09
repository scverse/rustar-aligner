---
title: Chimeric detection
description: Detect reads spanning two distant genomic locations — fusion candidates, structural variants, and circular RNA.
---

A chimeric alignment is a read whose two halves map to different genomic locations — different chromosomes, the same chromosome with an unrealistically large gap, or the same chromosome but on opposite strands. These are the candidate evidence for gene fusions, large-scale structural variants, and circular RNA back-splicing.

rustar-aligner uses STAR's two chimeric detectors, ported as written; which one runs depends on `--chimMultimapNmax`.

## Enabling chimeric detection

Chimeric detection is **off by default** (`--chimSegmentMin 0`). Enable it by setting a minimum chimeric segment length — STAR's recommended starting value is `12`:

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --readFilesCommand zcat \
  --chimSegmentMin 12 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_
```

Higher values (e.g. 20) produce fewer, more confident calls; lower values (e.g. 10) are more sensitive but noisier.

## What gets detected

Both detectors work on the read's window transcripts, the alignments stitching produced before any filtering. For paired-end reads these are whole-fragment transcripts (`mate1 | RC(mate2)`), so one side of a chimera can be a complete mate pair.

- **Default (`--chimMultimapNmax 0`):** STAR's `chimericDetectionOld`. The best linear alignment is one segment, and the detector looks for the single best partner that covers the rest of the read. A chimera is reported only if it is clearly better than the runner-up (`--chimScoreSeparation`), and at most one per read.
- **`--chimMultimapNmax N > 0`:** STAR's `chimericDetectionMult`. Every pair of window transcripts is tried and re-scored after the junction is placed, and all chimeras within `--chimMultimapScoreRange` of the best are reported, up to `N` (more than `N` reports none). The junction file gains six columns: number of chimeras, maximum possible score, best non-chimeric score, this chimera's score, best chimeric score, and whether the mates were merged.

When the junction falls inside a read, its position is chosen by STAR's scan over the overlap, preferring GT/AG or CT/AC motifs; when the two mates of a pair lie on either side of it, the junction type is `-1`.

## Output formats

Set `--chimOutType` to control the output. Multiple values are allowed.

### `Junctions` (default)

```bash
--chimOutType Junctions
```

Writes a `<prefix>Chimeric.out.junction` file with one row per chimeric junction. The 14-column format matches STAR's; tools like [Arriba](https://github.com/suhrig/arriba) and [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) consume it directly.

### `WithinBAM`

```bash
--chimOutType WithinBAM [HardClip|SoftClip]
--outSAMtype BAM Unsorted        # or SortedByCoordinate; WithinBAM needs BAM output
```

Writes the chimera into the main BAM **in place of the read's normal alignment**, as STAR does: one segment as a normal record and the other as a supplementary record (FLAG `0x800`), with `SA` tags pointing at each other. The supplementary segment is hard-clipped by default (`HardClip`) or soft-clipped (`SoftClip`). For paired-end reads a segment that covers both mates is written as a normal pair; when the mates lie on either side of the junction, both are written as normal records with no supplementary. `NM` is added to the output attributes, as STAR does. This is the input [Arriba](https://github.com/suhrig/arriba) recommends.

Because these reads are written as chimeras, STAR does not count them as uniquely or multi-mapped in `Log.final.out`, and they do not contribute to `SJ.out.tab` or gene counts; rustar-aligner does the same.

### Mixed output

```bash
--chimOutType Junctions WithinBAM
```

Writes both the junction file and the chimeric BAM records.

## Tuning parameters

The most useful chimeric parameters:

- `--chimSegmentMin` — minimum chimeric segment length (also enables/disables detection).
- `--chimScoreMin` — minimum total chimeric alignment score. Default `0`.
- `--chimScoreSeparation` — minimum score gap between the chosen chimeric pair and the next-best alternative. Default `10`.
- `--chimJunctionOverhangMin` — minimum bases on each side of the chimeric junction. Default `20`.
- `--chimMainSegmentMultNmax` — main segment can multimap up to this many loci. Default `10`.
- `--chimScoreJunctionNonGTAG` — score penalty for non-canonical chimeric junctions. Default `-1`.

See the [CLI parameters reference](/rustar-aligner/reference/cli-parameters/) for the rest.

## Output columns

The `Chimeric.out.junction` file has 14 tab-separated columns (STAR-compatible):

| # | Column | Meaning |
|---|--------|---------|
| 1 | `chr_donorA` | Donor (left segment) chromosome |
| 2 | `brkpt_donorA` | Donor breakpoint position |
| 3 | `strand_donorA` | Donor strand |
| 4 | `chr_acceptorB` | Acceptor (right segment) chromosome |
| 5 | `brkpt_acceptorB` | Acceptor breakpoint position |
| 6 | `strand_acceptorB` | Acceptor strand |
| 7 | `junction_type` | -1 (mates on either side of the junction) / 0 (non-canonical) / 1 (GT/AG) / 2 (CT/AC) |
| 8 | `repeat_left_lenA` | Length of repeat to the left |
| 9 | `repeat_right_lenB` | Length of repeat to the right |
| 10 | `read_name` | Source read name |
| 11 | `start_alnA` | Donor segment's first aligned genomic position |
| 12 | `cigar_alnA` | Donor CIGAR; a paired segment carries the second mate after a `p` operation (the genomic gap between mates, negative when they overlap) |
| 13 | `start_alnB` | Acceptor segment's first aligned genomic position |
| 14 | `cigar_alnB` | Acceptor CIGAR, same convention |
