---
title: Introduction
description: What rustar-aligner is, why it exists, and where it fits.
---

**rustar-aligner** is a Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference) — the widely-used RNA-seq aligner originally written in C++ by Alexander Dobin.

It aims to be a *faithful port* of STAR. That means:

- The same **genome index format** on disk
- The same `--camelCase` **command-line parameters**
- **SAM/BAM output** that matches STAR's byte-for-byte where the algorithm is deterministic, and is provably equivalent (same alignment set, different tie-break) where it isn't

If you have a STAR-based pipeline, you can swap the binary and the rest of the pipeline doesn't need to know.

## Why a Rust port?

The original STAR is excellent, battle-tested, and unlikely to be replaced. rustar-aligner exists to:

- Bring STAR's behaviour to a memory-safe toolchain with modern dependencies
- Provide a maintainable base for future RNA-seq alignment work that diverges from STAR (new features, different tradeoffs)
- Verify, by reimplementing it, that STAR's algorithms and thresholds are well-understood — every divergence found in the process becomes a documented edge case rather than a mystery

## What's supported today

End-to-end RNA-seq alignment with all of the features most pipelines need:

- Single-end and paired-end alignment with mate rescue
- SAM, unsorted BAM, and coordinate-sorted BAM output
- Multi-threaded parallel alignment (`--runThreadN`)
- GTF-based junction annotation with scoring bonus (`--sjdbGTFfile`)
- Two-pass mode for novel junction discovery (`--twopassMode Basic`)
- Chimeric alignment detection — 4-tier pipeline including multi-junction Tier 3
- Gene-level read counting (`--quantMode GeneCounts`)
- Transcriptome-coordinate SAM output (`--quantMode TranscriptomeSAM`)
- Splice junction output (`SJ.out.tab`)
- Unmapped read output to FASTQ
- Gzip-compressed FASTQ input
- Read group tags (`--outSAMattrRGline`)
- Seeded RNG for reproducible tie-breaking (`--runRNGseed`)
- All the standard SAM optional tags: NH, HI, AS, NM, nM, XS, jM, jI, MD
- `Log.final.out` statistics file (STAR-compatible, MultiQC-parseable)

See the [STAR compatibility report](/rustar-aligner/reference/star-compatibility/) for the detailed alignment-by-alignment comparison.

## What's not yet there

- **STARsolo** single-cell features

## Project status

rustar-aligner is the result of a phased reimplementation of STAR's algorithms in Rust. It's reached behavioural parity for the core alignment pipeline; the test suite (396 tests, 0 clippy warnings) plus a yeast differential-testing benchmark against STAR 2.7.x keeps it honest. Most of the code was written by [Claude Code](https://claude.ai/code) under direction from the project maintainer.

The codebase is licensed MIT (matching STAR's license).

## Where to next

- [Installation](/rustar-aligner/getting-started/installation/) — build from source
- [Quick start](/rustar-aligner/getting-started/quick-start/) — index a genome and align some reads
- [Migrating from STAR](/rustar-aligner/guides/migrating-from-star/) — swap-in checklist
