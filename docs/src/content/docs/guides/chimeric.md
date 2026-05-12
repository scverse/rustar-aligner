---
title: Chimeric detection
description: Detect reads spanning two distant genomic locations — fusion candidates, structural variants, and circular RNA.
---

A chimeric alignment is a read whose two halves map to different genomic locations — different chromosomes, the same chromosome with an unrealistically large gap, or the same chromosome but on opposite strands. These are the candidate evidence for gene fusions, large-scale structural variants, and circular RNA back-splicing.

rustar-aligner implements STAR's chimeric detection pipeline with four tiers, all of which run automatically when chimeric detection is enabled.

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

The 4-tier pipeline runs in this order, stopping as soon as a chimeric pair is found:

1. **Tier 1 — transcript-pair search.** Searches the read's existing transcript pool for two segments that together cover most of the read but map to incompatible locations.
2. **Tier 2 — multi-cluster.** When the seed pool produces multiple distinct alignment clusters, evaluates each pair as a candidate chimeric.
3. **Tier 1b — soft-clip re-mapping.** Takes the soft-clipped tail of the primary alignment and re-seeds it against the genome. Recovers chimeric pairs where the original aligner only kept one half.
4. **Tier 3 — residual outer re-seeding.** For reads where Tier 1/1b/2 has already found a chimeric pair, re-seeds the *remaining* uncovered regions of the read (before the donor / after the acceptor). Enables 3-way detection — gene fusions involving three loci, e.g. a complex rearrangement where two breakpoints are present in a single read.

For paired-end data, additional inter-mate detection runs: if mate 1 maps confidently to one location and mate 2 to another that's incompatible with a normal proper pair (different chromosomes, same strand, or >1 Mb gap), the pair is reported as chimeric.

## Output formats

Set `--chimOutType` to control the output. Multiple values are allowed.

### `Junctions` (default)

```bash
--chimOutType Junctions
```

Writes a `<prefix>Chimeric.out.junction` file with one row per chimeric junction. The 14-column format matches STAR's; tools like [Arriba](https://github.com/suhrig/arriba) and [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) consume it directly.

### `WithinBAM`

```bash
--chimOutType WithinBAM
```

Embeds the chimeric segments as supplementary alignment records (FLAG `0x800`) in the primary BAM, with `SA` tags linking the donor and acceptor halves. This is the format expected by tools that process chimeric BAM directly (e.g. for fusion calling on already-sorted BAMs).

### Mixed output

```bash
--chimOutType Junctions WithinBAM
```

Writes both the junction file and the supplementary BAM records. Useful when downstream tools have different format requirements.

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
| 7 | `junction_type` | -1 (encompassing PE) / 0 (non-canonical) / 1 (GT/AG) / 2 (CT/AC) |
| 8 | `repeat_left_lenA` | Length of repeat to the left |
| 9 | `repeat_right_lenB` | Length of repeat to the right |
| 10 | `read_name` | Source read name |
| 11 | `start_alnA` | Donor start position on the read |
| 12 | `cigar_alnA` | Donor CIGAR |
| 13 | `start_alnB` | Acceptor start position on the read |
| 14 | `cigar_alnB` | Acceptor CIGAR |
