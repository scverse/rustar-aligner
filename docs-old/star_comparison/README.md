[← Back to docs/](../)

# STAR vs rustar-aligner: Source Code Comparison

This directory contains detailed comparisons of each pipeline step between the original STAR C++ source and the rustar-aligner Rust port. The goal is to identify any remaining divergences and track their expected impact on alignment output.

All STAR source references are from the STAR master branch at https://github.com/alexdobin/STAR

## Pipeline Overview

```
FASTQ Input
    ↓
Seed Finding         (STAR: ReadAlign.cpp / seed functions)         → rustar-aligner: src/align/seed.rs
    ↓
Seed Clustering      (STAR: ReadAlign_stitchWindowSeeds context)    → rustar-aligner: cluster_seeds() in stitch.rs
    ↓
Pre-DP Seed Eval     (STAR: stitchWindowSeeds.cpp)                  → rustar-aligner: stitch_seeds_with_jdb_debug()
    ↓
Recursive Stitching  (STAR: stitchWindowAligns.cpp)                 → rustar-aligner: stitch_recurse() + stitch_seeds_core()
    ↓
Per-Step Stitching   (STAR: stitchAlignToTranscript.cpp)            → rustar-aligner: stitch_align_to_transcript()
    ↓
Seed Extension       (STAR: extendAlign.cpp)                        → rustar-aligner: extend_alignment() / finalize_transcript()
    ↓
Quality Filter       (STAR: ReadAlign_mappedFilter.cpp)             → rustar-aligner: filter in read_align.rs
    ↓
SAM Output           (STAR: ReadAlign_outputAlignments.cpp)         → rustar-aligner: src/io/sam.rs
```

## Documents

| File | Topic | Status |
|------|-------|--------|
| [01_recursive_stitcher.md](01_recursive_stitcher.md) | `stitchWindowAligns.cpp` vs `stitch_recurse` | Draft |
| [02_stitch_align_to_transcript.md](02_stitch_align_to_transcript.md) | `stitchAlignToTranscript.cpp` vs `stitch_align_to_transcript` | Draft |
| [03_extend_align.md](03_extend_align.md) | `extendAlign.cpp` vs `extend_alignment` | Draft |
| [04_seed_window_dp.md](04_seed_window_dp.md) | `stitchWindowSeeds.cpp` vs pre-DP in stitch.rs | Draft |
| [05_mapped_filter.md](05_mapped_filter.md) | `ReadAlign_mappedFilter.cpp` vs quality filter | Draft |
| [DIFFERENCES.md](DIFFERENCES.md) | All identified differences, impact assessment | **Active** |

## Quick Reference: STAR Source Files

| STAR file | Purpose | rustar-aligner equivalent |
|-----------|---------|-------------------|
| `stitchWindowAligns.cpp` | Recursive include/exclude stitcher | `stitch_recurse` + base case in `stitch_seeds_core` |
| `stitchAlignToTranscript.cpp` | Per-step gap/splice/indel scoring | `stitch_align_to_transcript` |
| `extendAlign.cpp` | Soft-clip extension at read ends | `extend_alignment` (in `finalize_transcript`) |
| `ReadAlign_stitchWindowSeeds.cpp` | Forward-DP seed chain selection (pre-recursive) | Phase 16.7b pre-DP in `stitch_seeds_with_jdb_debug` |
| `ReadAlign_mappedFilter.cpp` | Post-alignment quality filtering | `filter_transcripts` in `read_align.rs` |
| `ReadAlign_outputAlignments.cpp` | Write SAM/BAM, record junctions | `lib.rs` + `sam.rs` + `junction/sj_output.rs` |
| `genomeSAindex.cpp` | SAindex construction | `src/index/sa_index.rs` |

## Key Conventions

- **STAR coordinate**: `rAend` = last base of seed A (0-based, inclusive). `gAend` = same for genome.
- **rustar-aligner coordinate**: `last_exon.read_end` = first base AFTER seed A (exclusive). `last_exon.genome_end` = same.
- Therefore: `STAR rAend = rustar-aligner last_exon.read_end - 1`; `STAR rBstart = rustar-aligner eff_read_pos`.
- **STAR jR**: distance from `rAend` into the gap/beyond; `jR=0` means junction right at `rAend` (no shift relative to seed end).
- **rustar-aligner jr_shift**: distance from end of shared region; `jr_shift=0` means junction at start of seed B.
- **Relationship**: `STAR_jR = jr_shift + shared` (where `shared = rGap = STAR's read gap`).
- **STAR WA_gStart**: stored in FORWARD genome coordinates (even for reverse-strand seeds), converted via `a1 = nGenome - (aLength + a1)`. rustar-aligner follows the same convention after Phase 16.27.
