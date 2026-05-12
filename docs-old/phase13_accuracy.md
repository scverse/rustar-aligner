[← Back to ROADMAP](../ROADMAP.md)

# Phase 13: Performance + Accuracy Optimization

**Status**: Complete (Phases 13.1-13.14)

**Goal**: Optimize alignment performance to approach STAR speeds, fix classification issues, and reduce accuracy gaps.

---

## Phase 13.1: Critical Bug Fixes ✅ (2026-02-07)

- Fixed clustering explosion (11,982 clusters → ~200 with STAR limits)
- Implemented exon building from CIGAR (line 354 TODO completed)
- Added match scoring (+1 per matched base for DP stitching)
- Implemented soft clip logic for partial alignments
- 100 reads: 100% mapped in ~5s (was: >60s, 0% mapped)

**Files**: `src/params.rs`, `src/align/stitch.rs`, `src/align/read_align.rs`
**See**: [ALIGNMENT_FIXES.md](../ALIGNMENT_FIXES.md)

---

## Phase 13.2: Mismatch Counting + Seed Expansion Bugs ✅ (2026-02-07)

**Problem**: 66% unmapped reads (STAR: 11%) due to inflated mismatch counts.

**Root Causes & Fixes**:
1. **Reverse-strand genome offset** — `count_mismatches()` was reading forward genome for reverse reads. Fixed: add `n_genome` offset.
2. **Seed length overestimation** — `extend_match()` verified at `sa_start` only. Fixed: `verify_match_at_position()` per expanded seed.
3. **DP chain start wrong** — Used `expanded_seeds.first()` instead of tracing `prev_seed`. Fixed chain start tracing.
4. **Combined gap CIGAR missing** — When `read_gap > 0` AND `genome_gap > 0`, seeds merged incorrectly. Fixed to emit `Match(shared) + Del/Ins(excess)`.
5. **Removed incorrect read reverse-complement** — Seeds match read-as-is against RC genome region.

| Metric | Before | After | STAR |
|--------|--------|-------|------|
| Unique | 32% | **82.4%** | 82% |
| Multi | 2% | **8.3%** | 7% |
| Unmapped | 66% | **9.3%** | 11% |

**Files**: `src/align/stitch.rs`, `src/align/score.rs`, `src/align/read_align.rs`

---

## Phase 13.3: Performance Optimization ✅ (2026-02-07)

**Profiling** (1000 reads, single-threaded): 29.5% malloc, 18.1% PackedArray::read, 15.9% cluster_seeds, 14.7% stitch_seeds.

**Optimizations**:
1. **Eliminated Vec allocations** — `Seed::genome_positions()` lazy iterator, pre-allocated Vecs
2. **Optimized PackedArray::read** — Fast path: direct 8-byte slice read
3. **Binary search `position_to_chr`** — O(n) → O(log n) via `partition_point()`
4. **Deduplicated expanded seeds** — Cap at 200 per cluster
5. **Deferred CIGAR clone** — Build only once for winner

| Metric | Before | After |
|--------|--------|-------|
| Wall time (1000 reads) | ~3.0s | **~0.77s** (3.9x faster) |

**Files**: `src/align/seed.rs`, `src/index/packed_array.rs`, `src/genome/mod.rs`, `src/align/stitch.rs`

---

## Phase 13.4: CIGAR Integer Overflow + Coordinate Bugs ✅ (2026-02-09)

**Symptoms**: CIGAR `4294953882D` (near 2³²), junction coords beyond chr boundaries.

**Root Causes**:
1. **Integer overflow** — Gap calculations produce negative values for overlapping seeds; `as u32` wraps
2. **CIGAR merging failure** — Consecutive Match ops not merged
3. **Global vs per-chromosome coords** — Missing `pos - genome.chr_start[chr_idx]`

**Files**: `src/align/stitch.rs`, `src/io/sam.rs`, `src/junction/sj_output.rs`

---

## Phase 13.5: Scoring Fix ✅ (2026-02-09)

- `outFilterIntronMotifs` default → `None` (was wrongly `RemoveNoncanonical`)
- `seedMultimapNmax` restored to 10000 (was 100)
- Removed multi-chr anchor check (STAR processes each independently)
- STAR-style seed overlap trimming (advance start, don't skip)
- Gap mismatch scoring: `shared_score = shared_bases - 2*mismatches`
- `alignSJstitchMismatchNmax` filter

**Files**: `src/params.rs`, `src/align/stitch.rs`

---

## Phase 13.6: Alignment Extension (extendAlign) ✅ (2026-02-09)

**Problem**: 65% soft clips vs STAR's 26%, only 42% position agreement.

**Implementation**: `extend_alignment()` walks base-by-base, scoring +1 match / -1 mismatch, stops at mismatch limit. Applied at both ends of seed chain, replacing unconditional soft clips.

| Metric | Before | After | STAR |
|--------|--------|-------|------|
| Soft clip rate | 65% | **26.6%** | 25.8% |
| Position agreement | 42% | **51.0%** | — |

**Files**: `src/align/stitch.rs`, `src/align/score.rs`, `src/lib.rs`

---

## Phase 13.7: Reverse-Strand Splice Motif Fix ✅

- Added `CtAc`, `CtGc`, `GtAt` variants to `SpliceMotif` enum
- `detect_splice_motif()` always reads forward genome, checks all 6 patterns
- 100% motif agreement on 31/31 shared junctions (was 89.7%)

**Files**: `src/align/score.rs`

---

## Phase 13.8: Splice Junction Overhang Minimum ✅

Enforced `alignSJoverhangMin` (5bp) and `alignSJDBoverhangMin` (3bp) during DP stitching. Previously only used in two-pass junction filtering.

**Files**: `src/align/score.rs`, `src/align/stitch.rs`, `src/lib.rs`

---

## Phase 13.8b: Enforce alignIntronMax ✅

When `alignIntronMax=0` (default), STAR computes max as `2^winBinNbits * winAnchorDistNbins = 589,824bp`. Gaps exceeding this are treated as deletions. Added enforcement in both DP scoring and two-pass filtering.

**Files**: `src/align/score.rs`, `src/align/stitch.rs`, `src/lib.rs`, `src/junction/mod.rs`

---

## Phase 13.8c: Reduce False Non-Canonical Splice Junctions ✅

- **P0**: Fixed overhang calculation (was hardcoded `5u32`). Walks CIGAR for `min(left_exon, right_exon)`.
- **P1**: `outFilterIntronStrands = RemoveInconsistentStrands` (STAR default). Rejects conflicting junction strands.
- **P2**: `outSJfilter*` params (OverhangMin, CountUniqueMin, CountTotalMin, DistToOtherSJmin) — 4-element Vecs indexed by motif category.
- **P3**: Two-pass `filter_novel_junctions()` uses motif-specific thresholds.

**Files**: `src/lib.rs`, `src/params.rs`, `src/align/score.rs`, `src/align/read_align.rs`, `src/junction/sj_output.rs`, `src/junction/mod.rs`

---

## Phase 13.9: Fix Position Agreement ✅

**Root Cause**: SA reverse-strand position encoding bug. Positions stored as RC genome offsets `[0, n_genome)` but used directly as forward coords. ~94% of disagreements were different-chromosome.

**Fixes**:
1. `sa_pos_to_forward()` — converts `n_genome - sa_pos - match_length` for reverse strand
2. `cluster_seeds()` — convert before `position_to_chr()`
3. `stitch_seeds()` — convert for chr bounds; keep raw for DP genome access
4. SAM SEQ/QUAL — reverse-complement SEQ, reverse QUAL for FLAG & 16

| Metric | Before | After |
|--------|--------|-------|
| Position agreement | 51.0% | **94.5%** |
| Diff-chr MAPQ=255 | 3497 | **2** |
| STAR-only mapped | 629 | **42** |

**Files**: `src/index/mod.rs`, `src/align/stitch.rs`, `src/io/sam.rs`, `src/io/fastq.rs`

---

## Phase 13.9b: CIGAR Reversal + Splice Motif Fix + Genomic Length Penalty ✅

1. **CIGAR reversal** for reverse strand — DP builds in RC order, SAM requires forward
2. **Strand-aware splice motif** — `score_gap_with_strand()` converts RC donor: `forward_donor = n_genome - rc_donor - intron_len`
3. **`scoreGenomicLengthLog2scale`** (-0.25) — penalizes long-spanning alignments: `ceil(log2(span) * scale - 0.5)`

| Metric | Before | After |
|--------|--------|-------|
| Position agreement | 94.5% | **95.3%** |
| CIGAR agreement | 84.3% | **96.5%** |
| Splice rate | 5.8% | **4.1%** |

**Files**: `src/align/stitch.rs`, `src/align/score.rs`, `src/params.rs`

---

## Phase 13.9c: Deterministic Multi-Mapper Tie-Breaking ✅

Added secondary sort keys for equal-score transcripts: smallest chr_idx → smallest genome_start → forward strand first. Output now deterministic across runs.

**Files**: `src/align/read_align.rs`

---

## Phase 13.10: Accuracy Parity with STAR ✅

### Sub-phases:
- **13.10a**: Extension mismatch boundary — confirmed STAR uses `>` (not `>=`). No change.
- **13.10b**: Terminal exon overhang — 12bp min for novel, 3bp for annotated (`alignSJDBoverhangMin`). `stitch_seeds_with_jdb()`.
- **13.10c**: `outSJfilterIntronMaxVsReadN` — [50000, 100000, 200000] max intron by read count.
- **13.10d**: `winReadCoverageRelativeMin` (0.5) — discard sparse clusters.
- **13.10e**: Annotation bonus during DP — `sjdbScore` (+2) for annotated junctions in DP transitions.
- **13.10f**: Seed/window caps — `seedPerReadNmax` (1000), `seedPerWindowNmax` (50), `alignWindowsPerReadNmax` (10000).

| Metric | Before | After | STAR |
|--------|--------|-------|------|
| Position agreement | 95.3% | **96.3%** | — |
| CIGAR agreement | 96.5% | **97.4%** | — |
| Splice rate | 4.1% | **0.4%** | 2.5% |
| rustar-aligner-only junctions | 33 | **3** | — |

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`, `src/align/seed.rs`, `src/params.rs`, `src/junction/sj_output.rs`

---

## Phase 13.11: Bidirectional R→L Seed Search ✅

**Implementation** (`src/align/seed.rs`):
- `Seed.search_rc: bool` field, `reverse_complement_read()` helper
- `find_seeds()` runs L→R loop then R→L loop on RC read, shared `seedPerReadNmax` cap
- `genome_positions()` converts RC seeds: `(n_genome - pos - len, !is_rev)`

| Metric | Before | After |
|--------|--------|-------|
| Multi-mapped | 3.65% | **4.92%** |
| Shared junctions | 9 | **30** |
| Splice rate | 0.4% | **0.9%** |

**Files**: `src/align/seed.rs`, `src/align/stitch.rs`

---

## Phase 13.12: SJ.out.tab Motif/Strand Fix ✅

**Root Cause**: Junction strand derived from `transcript.is_reverse` (wrong). STAR derives from splice motif dinucleotides.

**Fixes**:
1. Strand from `motif.implied_strand()` (GT/AG→1, CT/AC→2, None→0)
2. Direct `encode_motif()` mapping without strand transformation

**Result**: Motif agreement 80% → **100%** on shared junctions.

**Files**: `src/lib.rs`, `src/junction/sj_output.rs`

---

## Phase 13.13: Relax Terminal Exon Overhang Filter ✅

Removed 12bp terminal exon floor from DP (Phase 13.10b was too aggressive). STAR only applies 12bp at SJ.out.tab write time, not during DP. Simplified to `alignSJoverhangMin` (5bp novel, 3bp annotated).

| Metric | Before | After | STAR |
|--------|--------|-------|------|
| Splice rate | 0.9% | **3.4%** | 2.2% |
| Shared junctions | 30 | **50** | 72 total |

**Files**: `src/align/stitch.rs`

---

## Phase 13.14: Implement `outFilterBySJout` ✅

After all reads aligned, compute surviving junctions via `outSJfilter*` thresholds, then filter reads with non-surviving junctions.

- `compute_surviving_junctions()` in `sj_output.rs`
- BySJout buffering in `lib.rs` — buffers all SAM records, filters, writes survivors
- `undo_mapped_record_bysj()` in `stats.rs` — atomic CAS for stat correction

| Mode | Position | CIGAR | Splice rate |
|------|----------|-------|-------------|
| Normal | 95.7% | 97.3% | 3.4% |
| **BySJout** | **96.7%** | **98.3%** | **1.1%** |

**Key Insight**: BySJout too aggressive without GTF — all junctions novel → strict thresholds. With GTF, annotated junctions bypass filters.

**Files**: `src/junction/sj_output.rs`, `src/junction/mod.rs`, `src/lib.rs`, `src/stats.rs`
