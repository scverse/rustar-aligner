# ![rustar-aligner](docs/src/assets/rustar-logo.svg)

A Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), the widely-used RNA-seq aligner originally written in C++ by Alexander Dobin.

## Overview

rustar-aligner aims to be a faithful port of STAR, matching the original behavior as closely as possible. It uses the same genome index format, accepts the same `--camelCase` command-line parameters, and produces compatible SAM/BAM output.

**Current status**: End-to-end single-end and paired-end RNA-seq alignment with splice junction detection, two-pass mode, chimeric alignment detection (including multi-junction Tier 3), gene-level quantification, and multi-threaded parallel processing. 396 tests passing (383 unit + 8 integration + others), 0 clippy warnings.

## Quick Start

### Build

```bash
cargo build --release
```

### Generate genome index

```bash
target/release/rustar-aligner --runMode genomeGenerate \
  --genomeDir /path/to/genome_index \
  --genomeFastaFiles /path/to/genome.fa
```

### Align reads

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outSAMstrandField intronMotif \
  --outFileNamePrefix /path/to/output_
```

### Paired-end alignment

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq reads_2.fq \
  --outSAMtype SAM \
  --outFileNamePrefix /path/to/output_
```

### BAM output

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix /path/to/output_
```

### Coordinate-sorted BAM output

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /path/to/output_
```

### Two-pass mode

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --twopassMode Basic \
  --outFileNamePrefix /path/to/output_
```

### Gene-level counts

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq \
  --sjdbGTFfile /path/to/annotation.gtf \
  --quantMode GeneCounts \
  --outFileNamePrefix /path/to/output_
```

## Accuracy Comparison vs STAR

Benchmarked on 10,000 yeast RNA-seq reads (150 bp, ERR12389696), compared to STAR 2.7.x with identical parameters and genome index.

### Single-End (10k reads, 150 bp SE)

| Metric | rustar-aligner | STAR |
|--------|----------------|------|
| Unique mapped | 82.6% | 82.6% |
| Multi-mapped | 7.4% | 7.4% |
| Total mapped | 90.0% | 90.0% |
| Position agreement | 96.5% raw / **99.815% tie-adjusted** | — |
| STAR-only reads | **0** | — |
| rustar-aligner-only reads | **0** | — |
| CIGAR-only diffs | 1 (seed-level tie in homopolymer) | — |

> **Tie-adjusted**: 299 of 313 disagreements are verified genuine ties — both tools find identical alignment sets but select different copies due to SA-order or RNG tie-breaking differences. Excluding these, faithfulness is 99.815% (8,611/8,627 non-tie reads exact).

### Paired-End (10k read pairs, 150 bp)

| Metric | rustar-aligner | STAR |
|--------|----------------|------|
| Both mates mapped | **8,390** | 8,390 |
| Half-mapped pairs | **0** | 0 |
| Unmapped pairs | 0 | 0 |
| PE faithfulness (tie-adjusted) | **99.883%** | — |
| MAPQ inflations | **0** | — |
| MAPQ deflations | **0** | — |
| NH tag diffs | **0** | — |
| Proper-pair diffs | **0** | — |

> **PE faithfulness**: 16,284 / 16,306 mate alignments exactly match STAR (same position, CIGAR, MAPQ, proper-pair flag, and NH tag). 475 diffs excluded as tie-breaking differences (same MAPQ+NH, different repeat copy chosen).

## Supported Features

- Single-end and paired-end alignment with mate rescue
- SAM, unsorted BAM, and coordinate-sorted BAM output (`--outSAMtype SAM`, `BAM Unsorted`, or `BAM SortedByCoordinate`)
- Multi-threaded parallel alignment (`--runThreadN`)
- GTF-based junction annotation with scoring bonus (`--sjdbGTFfile`)
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
- SJDB insertion into genome index at genomeGenerate time
- Chimeric alignment detection — SE and PE, 4-tier pipeline: transcript-pair search, multi-cluster, soft-clip re-seeding, residual outer re-seeding for multi-junction fusions (`--chimSegmentMin`)
- Gene-level read counting (`--quantMode GeneCounts` → `ReadsPerGene.out.tab`)
- Transcriptome-coordinate SAM output (`--quantMode TranscriptomeSAM`)
- Post-alignment read filtering (`--outFilterType BySJout`)
- Splice junction output (`SJ.out.tab`)
- Unmapped read output to FASTQ (`--outReadsUnmapped Fastx` → `Unmapped.out.mate1` / `mate2`)
- Gzip-compressed FASTQ input (`--readFilesCommand zcat`)
- Read group tags (`--outSAMattrRGline`)
- Seeded RNG for reproducible tie-breaking (`--runRNGseed`)
- SAM optional tags: NH, HI, AS, NM, nM, XS, jM, jI, MD
- `--outSAMattributes` control (Standard/All/None/explicit list)
- SECONDARY flag (0x100) on multi-mapper alignments
- Configurable output limits (`--outSAMmultNmax`)
- Bidirectional seed search with `scoreSeedBest` pre-extension
- Junction boundary optimization (jR scanning)
- Log.final.out statistics file (STAR-compatible, MultiQC-parseable)

## Known Limitations

- No STARsolo single-cell features (Phase 14, deferred)

See [ROADMAP.md](ROADMAP.md) for detailed implementation tracking.

## Building from Source

Requires Rust 2024 edition (rustc 1.85+).

```bash
cargo build --release    # Release build
cargo test               # Run tests
cargo clippy             # Lint
cargo fmt                # Format
```

## Development

The majority of rustar-aligner's code was written by [Claude Code](https://claude.ai/code) (Anthropic's AI coding assistant), with technical direction, architecture decisions, and validation by the project maintainer.

## License

MIT (matching the original STAR license)
