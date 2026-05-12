---
title: STAR compatibility
description: How identical rustar-aligner's output is to STAR's, with the known divergences listed in detail.
---

rustar-aligner is a faithful port. The goal is byte-for-byte identical output to STAR for every read where the algorithm is deterministic, and provably equivalent output (same alignment set, different tie-break) for every read where it isn't. This page is the long-form scoreboard.

The benchmark below uses **10,000 yeast RNA-seq reads** (150 bp, ERR12389696), aligned with both tools using identical parameters and the same genome index.

## Single-end summary

| Metric | rustar-aligner | STAR |
|--------|----------------|------|
| Unique mapped | 82.6% | 82.6% |
| Multi-mapped | 7.4% | 7.4% |
| Total mapped | 90.0% | 90.0% |
| **Position agreement (raw)** | **96.5%** | — |
| **Position agreement (tie-adjusted)** | **99.815%** | — |
| Reads mapped only by STAR | **0** | — |
| Reads mapped only by rustar-aligner | **0** | — |
| CIGAR-only differences | 1 | — |

**Tie-adjusted** means: of the 313 raw disagreements, 299 are verified genuine ties — both tools find *the same set of alignments*, but pick different copies as primary because of differences in suffix-array iteration order or RNG-based tie-breaking. Excluding those ties, faithfulness is 8,611 / 8,627 = **99.815%**.

## Paired-end summary

| Metric | rustar-aligner | STAR |
|--------|----------------|------|
| Both mates mapped | **8,390** | 8,390 |
| Half-mapped pairs | **0** | 0 |
| Unmapped pairs | 0 | 0 |
| **PE faithfulness (tie-adjusted)** | **99.883%** | — |
| MAPQ inflations | **0** | — |
| MAPQ deflations | **0** | — |
| NH tag differences | **0** | — |
| Proper-pair flag differences | **0** | — |

**Tie-adjusted** for paired-end: 16,284 / 16,306 mate alignments exactly match STAR (same position, CIGAR, MAPQ, proper-pair flag, NH tag). 475 differences are excluded as tie-breaking only (same MAPQ + same NH, different repeat copy chosen).

## Index format

The genome index format is identical. After Phase G3 (the SA tie-breaking fix), the suffix array for the yeast benchmark genome is **byte-for-byte identical** between STAR and rustar-aligner: 10,862 entry differences fixed → **0 remaining**.

This means an index built with one tool is loadable by the other.

## Where the differences come from

There are three categories of remaining difference, in order of size:

### 1. Tie-breaking (the bulk: ~3% of reads)

When two alignments have the same score, both tools pick a primary using a deterministic but different procedure:

- **STAR** uses a Mersenne Twister RNG (`mt19937`) seeded by `--runRNGseed`.
- **rustar-aligner** uses Rust's `StdRng` (ChaCha) seeded by `--runRNGseed`.

Both honour the seed and are reproducible across runs of the same tool — but the produced sequences differ between tools. So for ~3% of multi-mapped reads, primary vs. secondary flips. Total NH counts, AS scores, CIGARs, and the full alignment set are unaffected.

A subset (~100 of 313 SE diffs) involves "diff-chr ties": the read maps to several copies in a multi-copy region (e.g. rDNA), and the two tools pick different copies. Same alignment quality, different one chosen.

This is what the "tie-adjusted" numbers above account for.

### 2. CIGAR placement in homopolymer runs (1 read in 10,000)

`ERR12389696.13573895`: both tools align to `XV:218357 MAPQ=255 AS=133`, but rustar-aligner emits `100M1I45M4S` (insertion at read position 100) while STAR emits `108M1I37M4S` (insertion at 108). The 71-base seed is found at RC pos 29 (rustar-aligner) vs RC pos 37 (STAR) due to a different `Lmapped` chain path through a long homopolymer. Same diagonal, same score, different starting position → different insertion placement.

This is a seed-level tie. Real impact for downstream tools: effectively zero. To match STAR exactly we'd need to replicate STAR's exact `Lmapped` chain, which is a high-effort, low-value fix.

### 3. PE-specific cases where rustar-aligner finds a *better* alignment

Four PE alignments have a higher AS in rustar-aligner than STAR. These are not bugs — they're improvements:

- `ERR12389696.844151`: rustar-aligner finds `VIII:451791` with 0 mismatches; STAR finds `VII:1001391` with 6 mismatches.
- `ERR12389696.4972950`: rustar-aligner correctly emits a spliced mate 2; STAR's combined-window approach fails to stitch the spliced mate at the better location and emits unspliced.

In both, STAR's combined-window seeding fails to reach the higher-scoring alignment. We've decided not to artificially regress these.

## Faithfulness over time

The faithfulness numbers have been moving upwards as more STAR algorithm details have been replicated:

| Phase | PE tie-adjusted faithfulness |
|-------|------------------------------|
| Pre-Phase F1 | 99.755% |
| Phase F1 (`--runRNGseed`) | 99.755%* (RNG change reset baseline) |
| Phase G1 (junction-shift fix) | 99.865% |
| Phase G2 (recursion budget + overflow fix) | **99.883%** |

\* Phase F1 changed PE tie-breaking from SA-order to seeded `StdRng`, which shuffled which reads count as "tie-broken" without changing the underlying alignment quality.

## Other compatibility checks

- **Genome index files** (`Genome`, `SA`, `SAindex`, `chrName.txt`, `chrStart.txt`, `chrNameLength.txt`, `sjdbList.fromGTF.out.tab`, `transcriptInfo.tab`, `exonInfo.tab`): identical or equivalent.
- **`SJ.out.tab`**: matching format and contents for the yeast benchmark.
- **`Log.final.out`**: format matches STAR; MultiQC parses it without modification.
- **`Chimeric.out.junction`**: 14-column STAR-compatible format. Tools like Arriba and STAR-Fusion read it directly.
- **`ReadsPerGene.out.tab`**: 4-column STAR-compatible format.

## What this means in practice

For the overwhelming majority of bioinformatics workflows — read counting, differential expression, splice junction analysis, fusion calling, MultiQC reports — rustar-aligner's output is interchangeable with STAR's. The exceptions are:

- Workflows that depend on the *specific copy* a multi-mapper landed on (rare; in those cases, you should be using the full alignment set anyway, not just the primary).
- Reproducing exact byte-equality of SAM output across STAR and rustar-aligner runs (not achievable today; the RNG difference is the dominant source).

If you find a divergence not described here, please [open an issue](https://github.com/Psy-Fer/rustar-aligner/issues) — the project goal is to keep this list as short as possible.
