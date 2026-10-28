# ![rustar-aligner](docs/src/assets/rustar-logo.svg)

A Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), the widely-used RNA-seq aligner originally written in C++ by Alexander Dobin.

## Overview

rustar-aligner aims to be a faithful port of STAR, matching the original behavior as closely as possible. It uses the same genome index format, accepts the same `--camelCase` command-line parameters, and produces compatible SAM/BAM output.

**Current status (v0.1, working toward v0.2.0):** single-end and paired-end alignment, splice junctions, two-pass mode, chimeric detection, `--quantMode GeneCounts` / `TranscriptomeSAM`, STARsolo (Gene / GeneFull / SJ / Velocyto, barcode correction, UMI dedup, multi-mapper resolution, cell calling), WASP filtering, mate-overlap merging, coverage tracks, adapter clipping, sparse suffix arrays and multi-threading. On the [parity suite](#comparison-with-star-2711b) the genome indices are byte-identical to STAR's, and alignments, `SJ.out.tab`, gene counts, chimeric junctions and STARsolo output match STAR 2.7.11b record for record, apart from the [residual differences](#residual-differences) listed below. Every deliberate difference is in [DIVERGENCE.md](DIVERGENCE.md); the release plan is in [ROADMAP.md](ROADMAP.md#release-plan).

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

## Comparison with STAR 2.7.11b

Measured on 2026-10-28 with `main` at e598b5d and the STAR 2.7.11b release binary, 8 threads, on an 8-core Intel i7-11800H. [`test/parity_suite.sh`](test/parity_suite.sh) reproduces every number here. It downloads the public data (ENA runs and pinned nf-core/test-datasets commits), builds both indices, runs both tools with the same flags and compares every output file:

```bash
DATA=path/to/testdata STAR=path/to/STAR-2.7.11b RUSTAR=target/release/rustar-aligner \
  test/parity_suite.sh all
```

Alignment records are compared per read: with more than one thread STAR writes reads in a thread-dependent order, while rustar-aligner always writes them in input order.

| Dataset | Source |
|---------|--------|
| Yeast | *S. cerevisiae* R64-1-1 + Ensembl 110 GTF; [ERR12389696](https://www.ebi.ac.uk/ena/browser/view/ERR12389696), 150 bp PE, 10k and 100k pairs (`seqtk sample -s100`) |
| Human chr21 | GRCh38 chr21 + GENCODE chr21 GTF (nf-core/test-datasets); first 2M pairs of GEUVADIS [ERR188021](https://www.ebi.ac.uk/ena/browser/view/ERR188021), 75 bp PE. Whole-transcriptome reads, so ~90% don't map to chr21 |
| Fusion | nf-core/test-datasets 40 kb chr22 slice + STAR-Fusion test reads (733k pairs) |
| Single-cell | nf-core/test-datasets scrnaseq: mouse chr19 (GRCm38.p6, GENCODE vM19), 10x v2 sample S10, 50k-barcode whitelist |

### Genome indices

`Genome`, `SA`, `SAindex`, the chromosome tables and all the annotation tables (`sjdbList.out.tab`, `exonInfo.tab`, `transcriptInfo.tab`, `geneInfo.tab`, `exonGeTrInfo.tab`, `sjdbInfo.txt`, ...) are byte-identical to STAR's for yeast, yeast with `--genomeSAsparseD 2`, human chr21, the chr22 slice and mouse chr19. Only `genomeParameters.txt` (command line) and `Log.out` differ.

### Alignment

| Dataset | Options | Records identical / STAR's records | Reads differing | `SJ.out.tab` |
|---------|---------|-------------------|-----------------|--------------|
| Yeast 10k | SE | 9,724 / 9,724 | 0 | identical |
| Yeast 10k | PE | 18,600 / 18,600 | 0 | identical |
| Yeast 100k | SE | 97,754 / 97,754 | 0 | identical |
| Yeast 100k | PE | 188,258 / 188,258 | 0 | identical |
| Yeast 100k | PE, `--quantMode TranscriptomeSAM GeneCounts`, sorted BAM | 188,258 / 188,258 | 0 | identical |
| Yeast 100k | PE, `--twopassMode Basic` | 188,262 / 188,262 | 0 | identical |
| Yeast 100k | PE, `--outFilterType BySJout` | 187,646 / 187,646 | 0 | identical |
| Yeast 100k | PE, `--outSAMattributes All --outSAMstrandField intronMotif --outSAMunmapped Within` | 218,994 / 218,996 | 1 (unmapped pair, `uT`) | identical |
| Yeast 100k | PE, `--genomeSAsparseD 2` index | 188,436 / 188,436 | 0 | identical |
| Yeast 100k | SE / PE, chimeric (`--chimSegmentMin 12`, `--chimMultimapNmax` 0 and 10) | all | 0 | identical |
| Human chr21 | SE (mate 1) | 541,935 / 541,935 | 0 | identical |
| Human chr21 | PE | 902,186 / 902,190 | 3 | 3,071 / 3,073 lines shared |
| Fusion | PE, chimeric (`--chimMultimapNmax` 0 and 10) | 80,592 / 80,656 | 24 | identical |

- **Gene counts** (`ReadsPerGene.out.tab`): identical on yeast. On chr21 PE, 274 of 276 lines are identical; the other 2 come from the reads that differ.
- **Transcriptome BAM** (`--quantMode TranscriptomeSAM`): every record is identical on yeast (138,084) and chr21 (28,972) except for *which* alignment of a multi-mapping read carries the primary flag. STAR draws that from a per-thread random stream, so its choice depends on how reads are split across threads. rustar-aligner reproduces STAR's `--runThreadN 1` choice at any thread count, and against STAR run single-threaded the yeast transcriptome BAM is identical record for record. See [DIVERGENCE.md §1.1](DIVERGENCE.md#11-multimapper-tie-breaking--rng).
- **Chimeric junctions** (`Chimeric.out.junction`): the same lines as STAR in every configuration: yeast SE 194, PE 717 / 882 (`--chimMultimapNmax` 0 / 10); chr21 PE 301 / 3,643; fusion 31 / 675. Only the order differs, for the threading reason above.
- **`Log.final.out`**: identical on yeast SE and 10k PE. Elsewhere the mapped counts differ only by the reads in the table, and a few unmapped reads are counted as "too short" rather than "other" (yeast 100k PE 1, chr21 PE 17, fusion 73).

### Single-cell (STARsolo)

Mouse chr19, 10x v2, `--soloType CB_UMI_Simple --soloFeatures Gene GeneFull SJ Velocyto`:

| Options | BAM records identical | `Solo.out` files identical |
|---------|-----------------------|----------------------------|
| defaults, sorted BAM with `CR UR CB UB GX GN sS sQ sM` | 19,955 / 19,955 | 36 / 36 |
| `--soloCellFilter EmptyDrops_CR --soloMultiMappers Uniform EM --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR` | (no BAM) | 34 / 40 |

In the default run every file matches: raw and filtered matrices, `barcodes.tsv` / `features.tsv`, `Summary.csv`, `Features.stats`, `Barcodes.stats` and `UMIperCellSorted.txt`. In the second run the EmptyDrops_CR cell calls match too, because rustar-aligner draws from the same random-number generator as STAR. The 6 files that differ are explained under [residual differences](#residual-differences).

### Residual differences

These aren't deliberate. Each is tracked to be fixed or explained before v0.2.0; [DIVERGENCE.md §5](DIVERGENCE.md#5-known-residual-single-read-differences) has the read names.

- **`nM` and mate order on some multi-mapping pairs** (chr21 2 pairs, fusion 21 pairs). Same loci, CIGARs, `AS` and MAPQ. On some alignments `nM` is 1–2 lower than STAR's, and in about half of these pairs the two mates are written in the opposite order.
- **Same-span alignments collapsed** (fusion, 3 pairs). STAR reports two alignments with the same mate positions but different indel or splice placement inside a mate; rustar-aligner reports one, so `NH` is one lower.
- **A pair STAR leaves unmapped** (chr21, 1 pair). rustar-aligner maps it to 10 loci (`AS` 136, spliced). This adds 2 lines to `SJ.out.tab`.
- **Unmapped-read category** ("too short" versus "other", `uT`), for the counts listed above. The yeast pair: STAR finds no alignment window (`uT:0`); rustar-aligner finds a 31-point alignment and rejects it as too short (`uT:1`).
- **STARsolo `GeneFull` multi-mappers** (`--soloMultiMappers`). STAR counts one multi-gene read whose barcode needed 1MM_multi correction, and rustar-aligner drops it. This changes `UniqueAndMult-*.mtx` and one `Features.stats` line.
- **STARsolo `SJ` with `--soloUMIdedup 1MM_CR`**. STAR collapses 4 more SJ UMIs than rustar-aligner does (SJ `matrix.mtx`, `Features.stats`, `Summary.csv`).

### Speed and memory

Same machine, inputs read from a warm page cache, `--outSAMtype BAM Unsorted`. Yeast is the median of 3 alternating runs; chr21 is a single run.

| Run | Threads | STAR wall | rustar-aligner wall | STAR peak RSS | rustar-aligner peak RSS |
|-----|---------|-----------|---------------------|---------------|-------------------------|
| Yeast, 1M pairs (150 bp) | 8 | 131.9 s | 85.4 s | 1.26 GB | 0.23 GB |
| Yeast, 1M pairs (150 bp) | 16 | 128.2 s | 81.6 s | 2.18 GB | 0.28 GB |
| Human chr21, 2M pairs (75 bp) | 16 | 1,073 s | 594 s | 2.05 GB | 0.52 GB |

Index builds at 8 threads are slower than STAR's (yeast 3.0 s vs 2.5 s, chr21 12.4 s vs 10.5 s, mouse chr19 20.2 s vs 16.0 s) and use much less memory (222 MB vs 341 MB, 388 MB vs 1,084 MB, 459 MB vs 1,645 MB). rustar-aligner memory-maps the index, so its peak RSS counts only the index pages a run touches. These are small genomes; whole-genome numbers on server hardware are still to come.

The yeast reads are `seqtk sample -s100` 1,000,000-pair subsamples of ERR12389696 (the same procedure as the parity data), decompressed. Each run was (and the same with `--runThreadN 16`):

```bash
/usr/bin/time -f "%e s %M KB" $BIN --runThreadN 8 --genomeDir $INDEX \
  --readFilesIn y1m_1.fq y1m_2.fq --outSAMtype BAM Unsorted --outFileNamePrefix out/
```

## Supported Features

- Single-end and paired-end alignment with mate rescue
- Read-end alignment mode (`--alignEndsType Local` (default) / `EndToEnd` / `Extend5pOfRead1` / `Extend5pOfReads12` / `Extend3pOfRead1`)
- SAM, unsorted BAM, and coordinate-sorted BAM output (`--outSAMtype SAM`, `BAM Unsorted`, `BAM SortedByCoordinate`, or both BAM files). The coordinate sort is external: it buffers up to `--limitBAMsortRAM`, spills sorted runs beside the output, and merges them, so peak memory is independent of output size
- Multi-threaded parallel alignment (`--runThreadN`)
- GTF-based junction annotation with scoring bonus (`--sjdbGTFfile`)
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
- SJDB insertion into genome index at genomeGenerate time
- Chimeric alignment detection — SE and PE, 4-tier pipeline: transcript-pair search, multi-cluster, soft-clip re-seeding, residual outer re-seeding for multi-junction fusions (`--chimSegmentMin`)
- Gene-level read counting (`--quantMode GeneCounts` → `ReadsPerGene.out.tab`)
- Transcriptome-coordinate SAM output (`--quantMode TranscriptomeSAM`)
- Transcript-degradation normalization, DegNorm-style (`--quantGeneCoverage Yes` during alignment, then `--degNormCoverageFiles` across samples) — per-gene Degradation Index scores and degradation-adjusted counts; not a STAR feature
- Bulk total RNA-seq: unspliced `<gene_id>-I` targets in the transcriptome BAM (`--quantTranscriptomeUnspliced`), an `sp:A` splicing-status tag (`--outSAMsplicingStatus`) and spliced / unspliced gene counts (`--quantGeneSplicing`); not STAR features
- **Single-cell quantification (STARsolo)** — `--soloType CB_UMI_Simple`, `CB_UMI_Complex` (multi-segment barcodes), and `SmartSeq` (plate-based, SE + PE); features `Gene`, `GeneFull` (pre-mRNA), `SJ`, and `Velocyto` (spliced/unspliced/ambiguous); barcode correction (`--soloCBmatchWLtype` Exact/1MM/1MM_multi/…), UMI dedup (`--soloUMIdedup` 1MM_All/1MM_CR/1MM_Directional/…), multi-gene UMI filtering, multi-mapper resolution (`--soloMultiMappers` Uniform/PropUnique/EM/Rescue), cell calling (`--soloCellFilter` CellRanger2.2/TopCells/EmptyDrops_CR), gzip output, and `Summary.csv` — writes STARsolo-compatible `Solo.out/<feature>/{raw,filtered}/{matrix.mtx, barcodes.tsv, features.tsv}`
- WASP allele-specific-mapping filter (`--waspOutputMode SAMtag`, `--varVCFfile`) — vW/vA/vG tags, SE + PE; BAM output only, as in STAR
- Paired-end mate-overlap merging (`--peOverlapNbasesMin`, `--peOverlapMMp`)
- Coverage-track output (`--outWigType bedGraph` → `Signal.{Unique,UniqueMultiple}.str{1,2}.out.bg`)
- Adapter and fixed clipping (`--clip5pNbases`, `--clip3pNbases`, `--clip3pAdapterSeq`, `--clipAdapterType` including CellRanger4 TSO + poly-A) — clipped bases soft-clipped in SAM, matching STAR's convention
- Sparse suffix array (`--genomeSAsparseD`) for reduced index memory
- Post-alignment read filtering (`--outFilterType BySJout`)
- Splice junction output (`SJ.out.tab`)
- Unmapped read output to FASTQ (`--outReadsUnmapped Fastx` → `Unmapped.out.mate1` / `mate2`)
- Gzip-compressed FASTQ input (`--readFilesCommand zcat`)
- Read group tags (`--outSAMattrRGline`)
- Deterministic input-order output regardless of thread count (`--outSAMorder Paired` / `PairedKeepInputOrder`, both accepted; rustar always preserves FASTQ order)
- Deterministic in-tree RNG for reproducible tie-breaking (`--runRNGseed`; no `rand` dependency)
- SAM optional tags: NH, HI, AS, nM, NM, XS, jM, jI, MD (default `NH HI AS nM` matches STAR; `NM` opt-in)
- `--outSAMattributes` control (Standard/All/None/explicit list)
- SECONDARY flag (0x100) on multi-mapper alignments
- Configurable output limits (`--outSAMmultNmax`)
- Bidirectional seed search with `scoreSeedBest` pre-extension
- Junction boundary optimization (jR scanning)
- Log.final.out statistics file (STAR-compatible, MultiQC-parseable)

## Known Limitations

Beyond the [residual differences](#residual-differences) above:

- **SmartSeq** (`--soloType SmartSeq`) runs, but its output is far from STARsolo's: only `Gene/raw` is written, and its counts differ.
- **`--soloCBwhitelist None`** writes no count matrices yet.
- **`CB_UMI_Complex`** barcode matching differs from STAR's for some rejected reads ([DIVERGENCE.md §3.4](DIVERGENCE.md)), and `EditDist_2` matching isn't implemented.
- **Not implemented** (rejected with an error): `--outWigType wiggle`, `--outWigStrand Unstranded`, `--chimOutType SeparateSAMold`, `--outSAMfilter KeepOnlyAddedReferences` / `KeepAllAddedReferences`, `--readFilesType SAM`, `--soloCBtype String`, `--soloBarcodeMate 2`.
- **`--genomeLoad`** is accepted but has no effect: the index is always memory-mapped, and STAR's shared-memory modes aren't implemented.
- **STARlong** (long-read alignment) isn't ported.

See [ROADMAP.md](ROADMAP.md) for the plan and [DIVERGENCE.md](DIVERGENCE.md) for every known difference.

## Building from Source

Requires Rust 2024 edition (rustc 1.89+).

```bash
cargo build --release       # Release build
cargo test                  # Run tests
cargo clippy --all-targets  # Lint
cargo fmt                   # Format
```

## Development

The majority of rustar-aligner's code was written by [Claude Code](https://claude.ai/code) (Anthropic's AI coding assistant), with technical direction, architecture decisions, and validation by the project maintainer.

## License

MIT (matching the original STAR license). Code ported from STAR, STAR-rs, DegNorm and
David Elworthy's `sgt.h` stays under its own licence, and libc++ is credited for the
random-number behaviour STARsolo's EmptyDrops_CR depends on; all of these notices are in
[`LICENSES/`](LICENSES/). The licences of the crates linked into the binaries ship with
them as `THIRD_PARTY_LICENSES.md` (`scripts/third-party-licenses.sh`).
