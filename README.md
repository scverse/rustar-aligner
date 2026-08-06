# ![rustar-aligner](docs/src/assets/rustar-logo.svg)

A Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), the widely-used RNA-seq aligner originally written in C++ by Alexander Dobin.

## Overview

rustar-aligner aims to be a faithful port of STAR, matching the original behavior as closely as possible. It uses the same genome index format, accepts the same `--camelCase` command-line parameters, and produces compatible SAM/BAM output.

**Current status**: End-to-end single-end and paired-end RNA-seq alignment with splice junction detection, two-pass mode, chimeric alignment detection (including multi-junction Tier 3), gene-level quantification, **single-cell quantification (STARsolo: Gene / GeneFull / SJ / Velocyto features, barcode correction, UMI dedup, EmptyDrops_CR cell calling)**, WASP allele-specific filtering, paired-end mate-overlap merging, coverage-track output, adapter/CellRanger clipping, and multi-threaded parallel processing. Default SAM output is byte-identical to STAR's (`NH HI AS nM` attributes); solo count matrices are byte-identical to STARsolo's. Pure-Rust core with an in-tree deterministic RNG for alignment tie-breaking. 622 tests passing (591 unit + 31 integration), 0 clippy warnings. See [Performance & Benchmarks](#performance--benchmarks) for a native three-way comparison against STARsolo and CellRanger.

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

### Single-cell (STARsolo)

10x Chromium-style single-cell quantification. `readFilesIn` takes the cDNA
read first, then the barcode read (`R2 R1` for 10x). Writes a
`Solo.out/<feature>/{raw,filtered}/` count matrix per requested feature:

```bash
target/release/rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn cDNA_R2.fq.gz barcode_R1.fq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/3M-february-2018.txt \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloCellFilter EmptyDrops_CR \
  --outSAMtype None \
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

## Performance & Benchmarks

Two single-cell benchmarks: a **mouse** three-way against STARsolo and CellRanger
(throughput on a compute-bound run), and a **human** run that is memory-bound,
where the sparse suffix array is the deciding factor.

### Mouse — 5k PBMC, three-way vs STARsolo and CellRanger

Single-cell (STARsolo) throughput and memory, benchmarked native x86_64 against
**STARsolo 2.7.11b** and **CellRanger 10.0.0** on a real 10x dataset
(`5k_Mouse_PBMCs_5p_gem-x_GEX`, 5′ GEM-X, GRCm39-2024-A reference). All three
tools ran on the **same EC2 instance, same NVMe-staged inputs, 10 threads, page
cache dropped between reps, no BAM output** — the only fair substrate, since
CellRanger is x86-only and STAR can't run natively on Apple Silicon. Harness:
[`test/aws/`](test/aws/).

| Tool | Wall time | Peak RSS | Cells called | Genes detected |
|------|-----------|----------|--------------|----------------|
| STARsolo 2.7.11b | **87 s** | 28.3 GB | 4,061 | 15,439 |
| rustar-aligner | 121 s | 25.7 GB | 3,689 | 16,284 |
| rustar-aligner `--genomeSAsparseD 2` | 119 s | **17.7 GB** | 3,692 | 16,311 |
| CellRanger 10.0.0 | 347 s | 13.1 GB | 3,858 | 16,858 |

- **Correctness**: rustar-aligner's `Solo.out/Gene/raw` count matrix is
  **byte-identical to STARsolo's** on the CellRanger-style flag set (verified
  deterministically, 3/3 runs). Differing "cells called" / "genes detected"
  reflect the cell-calling filter (`EmptyDrops_CR` Monte-Carlo), not the counts.
- **Memory**: `--genomeSAsparseD 2` (sparse suffix array, byte-identical to
  STAR's own D=2 SA) trades ~2% wall time for a **31% smaller peak RSS**,
  giving rustar the smallest memory footprint of any suffix-array aligner here.
  rustar's reported RSS includes reclaimable mmap'd index pages, so it overstates
  real memory pressure.
- **Since this run**: a segment-tree gene-assignment optimization (O(log n + k)
  overlap query, replacing STAR's linear scan) landed in `main` and independently
  measured a **further ~14% reduction in rustar's solo wall time (→ ~105 s)** on
  the same dataset/instance — the table above predates it and is a conservative
  floor. Count matrices are unchanged.

### Human — GRCh38, memory-bound (sparse vs dense index)

A full human snRNA-seq sample (GSE163577 `SRR13278444`, ~155M reads, 10x 3′ v3,
GRCh38) hits a different bottleneck. The **32 GB dense index does not stay
resident in 48 GB RAM**, so the run is dominated by page faults rather than
alignment — a production dense run took **9.5 h with 558M page faults** (≈4 h of
it kernel time). The **sparse suffix array (`--genomeSAsparseD 2`)** shrinks the
index to **16 GB** so it fits in RAM. Measured on a 500k-read subset of the same
sample (8 threads, warm cache):

| rustar-aligner (500k subset) | Wall | Page faults | Unique mapping |
|------------------------------|------|-------------|----------------|
| dense index (32 GB)          | 165–460 s\* | 5.1 M | 77.0% |
| `--genomeSAsparseD 2` (16 GB) | **54 s** | **0.32 M** | 77.6% |

<sub>\*dense wall swings with cache pressure because the 32 GB index can't stay
resident; the 16 GB sparse index is stable at ~54 s.</sub>

- **~3× faster resident, up to ~8× when the dense index is paging**, with **16×
  fewer page faults** — purely because 16 GB fits in 48 GB RAM and 32 GB does not.
- **Sensitivity is equal** (sparse 77.6% vs dense 77.0% unique). Sparse output is
  byte-identical to STAR's own **D=2** SA (not to a dense D=1 run), so counts
  differ from a dense run by ~1% of matrix entries — run a whole cohort on one
  index for consistency.
- **Takeaway**: on a RAM-constrained box, index size — not raw aligner speed — is
  the dominant cost for a human genome; `--genomeSAsparseD 2` is the lever.

## Supported Features

- Single-end and paired-end alignment with mate rescue
- Read-end alignment mode (`--alignEndsType Local` (default) / `EndToEnd` / `Extend5pOfRead1` / `Extend5pOfReads12` / `Extend3pOfRead1`)
- SAM, unsorted BAM, and coordinate-sorted BAM output (`--outSAMtype SAM`, `BAM Unsorted`, or `BAM SortedByCoordinate`). The coordinate sort is external: it buffers up to `--limitBAMsortRAM`, spills sorted runs beside the output, and merges them, so peak memory is independent of output size
- Multi-threaded parallel alignment (`--runThreadN`)
- GTF-based junction annotation with scoring bonus (`--sjdbGTFfile`)
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
- SJDB insertion into genome index at genomeGenerate time
- Chimeric alignment detection — SE and PE, 4-tier pipeline: transcript-pair search, multi-cluster, soft-clip re-seeding, residual outer re-seeding for multi-junction fusions (`--chimSegmentMin`)
- Gene-level read counting (`--quantMode GeneCounts` → `ReadsPerGene.out.tab`)
- Transcriptome-coordinate SAM output (`--quantMode TranscriptomeSAM`)
- **Single-cell quantification (STARsolo)** — `--soloType CB_UMI_Simple`, `CB_UMI_Complex` (multi-segment barcodes), and `SmartSeq` (plate-based, SE + PE); features `Gene`, `GeneFull` (pre-mRNA), `SJ`, and `Velocyto` (spliced/unspliced/ambiguous); barcode correction (`--soloCBmatchWLtype` Exact/1MM/1MM_multi/…), UMI dedup (`--soloUMIdedup` 1MM_All/1MM_CR/1MM_Directional/…), multi-gene UMI filtering, multi-mapper resolution (`--soloMultiMappers` Uniform/PropUnique/EM/Rescue), cell calling (`--soloCellFilter` CellRanger2.2/TopCells/EmptyDrops_CR), gzip output, and `Summary.csv` — writes STARsolo-compatible `Solo.out/<feature>/{raw,filtered}/{matrix.mtx, barcodes.tsv, features.tsv}`
- WASP allele-specific-mapping filter (`--waspOutputMode SAMtag`, `--varVCFfile`) — vW/vA/vG tags, SE + PE
- Paired-end mate-overlap merging (`--peOverlapNbasesMin`, `--peOverlapMMp`)
- Coverage-track output (`--outWigType bedGraph` → `Signal.{Unique,UniqueMultiple}.str{1,2}.out.bg`)
- Adapter and fixed clipping (`--clip5pNbases`, `--clip3pNbases`, `--clip3pAdapterSeq`, `--clipAdapterType` including CellRanger4 TSO + poly-A) — clipped bases soft-clipped in SAM, matching STAR's convention
- Sparse suffix array (`--genomeSAsparseD`) for reduced index memory
- Post-alignment read filtering (`--outFilterType BySJout`)
- Splice junction output (`SJ.out.tab`)
- Unmapped read output to FASTQ (`--outReadsUnmapped Fastx` → `Unmapped.out.mate1` / `mate2`)
- Gzip-compressed FASTQ input (`--readFilesCommand zcat`)
- Read group tags (`--outSAMattrRGline`)
- Deterministic input-order output regardless of thread count (`--outSAMorder Paired` / `PairedKeepInputOrder`, both accepted)
- Deterministic in-tree RNG for reproducible alignment tie-breaking (`--runRNGseed`)
- SAM optional tags: NH, HI, AS, nM, NM, XS, jM, jI, MD (default `NH HI AS nM` matches STAR; `NM` opt-in)
- `--outSAMattributes` control (Standard/All/None/explicit list)
- SECONDARY flag (0x100) on multi-mapper alignments
- Configurable output limits (`--outSAMmultNmax`)
- Bidirectional seed search with `scoreSeedBest` pre-extension
- Junction boundary optimization (jR scanning)
- Log.final.out statistics file (STAR-compatible, MultiQC-parseable)

## Known Limitations

- Solo cell-calling counts can differ from STARsolo by a few cells due to the
  Monte-Carlo `EmptyDrops_CR` rescue; the underlying count matrix is identical.
- A handful of SE/PE tie-breaking differences vs STAR remain (verified genuine
  ties; see the accuracy tables above and [ROADMAP.md](ROADMAP.md)).

See [ROADMAP.md](ROADMAP.md) for detailed implementation tracking.

## Building from Source

Requires Rust 2024 edition (rustc 1.85+).

```bash
cargo build --release       # Release build
cargo test                  # Run tests
cargo clippy --all-targets  # Lint
cargo fmt                   # Format
```

## Development

The majority of rustar-aligner's code was written by [Claude Code](https://claude.ai/code) (Anthropic's AI coding assistant), with technical direction, architecture decisions, and validation by the project maintainer.

## License

MIT (matching the original STAR license)
