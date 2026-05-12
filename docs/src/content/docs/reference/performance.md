---
title: Performance
description: Build time, alignment throughput, and memory usage benchmarks.
---

:::caution[Placeholder]
This page is a placeholder. Comprehensive performance benchmarks vs STAR — across multiple genome sizes, thread counts, and storage tiers — are still in progress and will land here when the data is in.
:::

For now, the headline points:

## What we know today

- **Algorithmic parity.** rustar-aligner implements the same core algorithms as STAR. There's no fundamental reason for it to be faster *or* slower; the constant factors come down to memory layout, allocator behaviour, and how aggressively the Rust compiler vectorises the inner loops.
- **Single-threaded runtime is comparable to STAR** on the yeast benchmark. Detailed numbers pending.
- **Memory footprint** for genome loading and alignment is in the same ballpark as STAR. The on-disk index format is identical, so the in-memory representation is similar by construction.
- **Multi-threaded scaling** uses Rust's `rayon` for the alignment phase. Scaling is approximately linear up to the IO limit of the input/output storage.

## What's coming

The benchmark suite under construction will report:

### Genome index generation

- Build time vs genome size (yeast, *Drosophila*, mouse, human GRCh38)
- Peak RAM during build
- Effect of `--genomeSAindexNbases` and `--runThreadN`

### Alignment

- Reads/second per core, single-end and paired-end
- Multi-thread scaling curve up to 32 threads
- Comparison vs STAR with identical parameters
- Effect of two-pass mode vs single-pass on total runtime
- Effect of chimeric detection on runtime

### Memory

- Peak RSS during alignment
- Effect of `--limitBAMsortRAM` on coordinate-sorted BAM output

### Output

- Stream-to-stdout (`--outStd`) vs file-to-disk
- BGZF compression level vs throughput tradeoff

## Reporting your own numbers

If you've benchmarked rustar-aligner on a real dataset, we'd love to hear about it — please [open an issue](https://github.com/Psy-Fer/rustar-aligner/issues) with the genome, read count, hardware, parameters, and timing. Real-world data shapes the optimisation work.
