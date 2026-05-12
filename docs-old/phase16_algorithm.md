[← Back to ROADMAP](../ROADMAP.md)

# Phase 16: Accuracy + Algorithm Parity

**Status**: In Progress (Phases 16.1-16.10 + 16.7b/16.7c complete)

**Goal**: Close remaining accuracy gaps vs STAR. Fix over-splicing, rDNA MAPQ, seed parameters, DP junction optimization, PE mate rescue, and multi-transcript DP.

---

## Phase 16.1: max_cluster_dist from winBinNbits ✅ (2026-02-13)

**Problem**: Hardcoded `max_cluster_dist = 100kb`. STAR computes `2^winBinNbits * winAnchorDistNbins` = 589,824bp (6x larger).

**Fix**:
- Added `--winBinNbits` (16) and `--winAnchorDistNbins` (9) params
- `win_bin_window_dist()` helper: `2^winBinNbits * winAnchorDistNbins`
- Replaced hardcoded 100kb and literal 589,824 with computed value

**Impact**: Splice rate **3.4% → 2.2%** (matches STAR exactly). CIGAR agreement 97.3% → 97.8%.

**Files**: `src/params.rs`, `src/align/read_align.rs`, `src/align/score.rs`, `src/junction/mod.rs`

---

## Phase 16.2: RemoveNoncanonicalUnannotated + GTF Testing ✅ (2026-02-13)

**Problem**: `RemoveNoncanonicalUnannotated` fell through to `RemoveNoncanonical`, rejecting ALL non-canonical junctions even when annotated.

**Fix**: `zip(junction_motifs, junction_annotated)` — only reject when `NonCanonical && !annotated`.

**GTF Testing**: Established differential baseline. STAR gains 8 more junctions with GTF because it inserts splice sites into genome index at alignment time (rustar-aligner does not).

**Files**: `src/align/read_align.rs`

---

## Phase 16.3: Junction Position Optimization (jR Scanning) ✅ (2026-02-13)

STAR's 3-phase jR scanning algorithm as post-DP optimization:
1. `find_best_junction_position()` — scan left, scan right (match quality + motif), repeat detection
2. `optimize_junction_positions()` — walks winning chain's CIGAR, applies to each RefSkip
3. `jr_shift` clamped to `[-prev_match_len, next_match_len]` to prevent CIGAR corruption

**Architecture**: Post-DP on winning chain's ~1-3 junctions. Zero performance impact.

**Results**: Neutral on yeast (+1 shared junction, 42 total). Expected benefit on mammalian genomes.

**Files**: `src/align/score.rs`, `src/align/stitch.rs`

---

## Phase 16.4: Seed Search Params + Sparse Infrastructure ✅ (2026-02-13)

- Added `seedSearchStartLmax` (50), `seedSearchStartLmaxOverLread` (1.0), `seedSearchLmax` (0), `seedMapMin` (5)
- Refactored `find_seed_at_position()` → `MmpResult` (always provides MMP advance length)
- `search_direction_sparse()` with STAR-matching Lmapped tracking — kept dormant (DP needs dense seeds)

**Files**: `src/params.rs`, `src/align/seed.rs`

---

## Phase 16.5: MAPQ Formula Fix ✅ (2026-02-16)

Replaced formula with STAR's lookup table: n=1→255, n=2→3, n=3-4→1, n≥5→0.

Threaded `n_for_mapq` through pipeline: `align_read()` → `align_paired_read()` → SAM builders. Currently `n_for_mapq = transcripts.len()` (no inflation).

**Files**: `src/mapq.rs`, `src/align/read_align.rs`, `src/io/sam.rs`, `src/io/bam.rs`, `src/lib.rs`

---

## Phase 16.5b: rDNA Window-Model Fix — DEFERRED

~157 chrXII rDNA reads get MAPQ=255 vs STAR 1-3. Root cause: 589kb clusters merge all tandem repeat copy seeds into ONE cluster → ONE transcript → MAPQ=255.

**Attempted**: Anchor-bin counting (3 iterations — raw, competitive, quality-filtered). All failed: 589kb max_cluster_dist means ALL clusters share the same seeds → identical transcripts.

**Correct fix**: Split clusters into ~65kb sub-windows before DP. Requires significant pipeline refactoring. `SeedCluster.anchor_bin` field kept for future use.

---

## Phase 16.6: Sparse Seed Bug Fixes ✅ (2026-02-17)

3 bugs fixed in `search_direction_sparse()`:
1. **While condition**: exit when `pos + seed_map_min >= read_len`
2. **Nstart**: `read_len.div_ceil(effective_start_lmax)` (not +1)
3. **RC read_pos**: convert `read_pos = original_read_len - read_pos - length` for R→L seeds

**Activation tested and reverted**: 91.1% pos (was 94.5%), 4.3% splice (was 2.1%). Root cause: sparse seeds at mismatch positions produce spurious locations without enough neighbors to vote them down. Bug-fixed function kept dormant.

**Files**: `src/align/seed.rs`

---

## Phase 16.7: Bin-Based Windowing ✅ (2026-02-18)

Replaced proximity-based clustering with STAR's bin-based windowing architecture:
1. `cluster_seeds()` rewritten: identify anchors → create windows from anchor bins → merge nearby windows (±winAnchorDistNbins) → extend by ±winFlankNbins → assign seeds by bin HashMap lookup
2. `WindowAlignment` struct with verified (seed_idx, sa_pos, strand, read_pos, length) — no SA range re-expansion
3. `stitch_seeds_with_jdb()` converts WA entries directly to ExpandedSeeds
4. Added `--winFlankNbins` param (default 4)

**Performance**: 12× faster (10s vs 120s on 10k reads) — bin HashMap lookup vs O(n²) proximity checks.

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`, `src/params.rs`

---

## Phase 16.7b: Pre-DP Seed Extension ✅ (2026-02-20)

**Problem**: Bin-based windowing (16.7) caused splice rate regression: 2.1% → 3.6%. Wide ~589KB windows allow coincidental short seed matches to create false splice junctions. STAR prevents this via seed extension during DP initialization.

**STAR's approach** (ReadAlign_stitchWindowSeeds.cpp):
1. Pass 1: DP base score = `seed_length + left_extend_score` — true positions extend well (+40bp), coincidental extend ~0bp
2. Pass 2: Right extension for chain endpoint selection — `dp_score + right_ext_score`
3. Post-DP: authoritative extensions replace approximate pre-DP scores

**Implementation** (all in `stitch_seeds_with_jdb()`):
1. **Pre-DP left extension**: `left_ext_scores: Vec<i32>` computed via `extend_alignment()` for every expanded seed before DP loop
2. **DP base case**: `score = seed_length + left_ext_scores[i]` (was `seed_length` only)
3. **Right extension for chain selection**: Replaced `max_by_key(score)` with loop computing rightward extension per chain endpoint, selecting by `dp[i].score + right_ext_score`
4. **Post-DP score adjustment**: `adjusted_score = best_state.score - left_ext_scores[chain_start_idx] + left_extend + right_extend`

| Metric | Pre-16.7 | 16.7 (WA-DP) | 16.7b (pre-DP ext) | STAR |
|--------|----------|--------------|---------------------|------|
| Position agree | 94.5% | 94.2% | **96.2%** | — |
| CIGAR agree | 97.8% | 97.1% | **97.6%** | — |
| Splice rate | 2.1% | 3.6% | **2.7%** | 2.2% |
| Shared junctions | 42 | 48 | **49** | 72 |

**Files**: `src/align/stitch.rs`

---

## Phase 16.8: PE Mate Rescue + Half-Mapped Output ✅ (2026-02-18)

**Problem**: 12.9% unmapped pairs because both mates must independently produce alignments.

**Fix**: 3-tier PE alignment with mate rescue:
1. Both map → pair concordantly
2. One fails → `rescue_unmapped_mate()` (genome-wide seeds → filter chr ± `alignMatesGapMax` → cluster → stitch) → pair or HalfMapped
3. Neither maps → unmapped

**`PairedAlignmentResult` enum**: `BothMapped(Box<PairedAlignment>)` | `HalfMapped { mapped_transcript, mate1_is_mapped }`

**Half-mapped output**: mapped mate FLAG 0x8, unmapped mate FLAG 0x4, co-located.

| Metric | With Rescue | Pre-Rescue | STAR |
|--------|-------------|------------|------|
| Both mapped | 8706 (96.6%) | 8714 | 8390 |
| Half-mapped | **311** | 0 (dropped) | 0 |
| Unmapped pairs | **0** | 286 | 0 |
| Per-mate pos agree | **95.6%** | 95.7% | — |
| STAR-only mates | **98** | 184 | — |
| Shared junctions | **76** | 72 | 90 total |

**Files**: `src/align/read_align.rs`, `src/io/sam.rs`, `src/stats.rs`, `src/lib.rs`

---

## Phase 16.7c: Anchor Threshold Fix ✅ (2026-02-20)

**Problem**: `max_loci_for_anchor` was hardcoded to 10. STAR uses `winAnchorMultimapNmax` (default 50).

**Fix**: Changed `max_loci_for_anchor` from hardcoded 10 to `params.win_anchor_multimap_nmax` in `align_read()` and `rescue_unmapped_mate()`.

**Impact**: Position **96.2% → 97.4%**. Shared junctions **49 → 57**. Splice rate 2.7% → 3.0% (slight increase).

**Files**: `src/align/read_align.rs`

---

## Phase 16.9: MMP SA Range Narrowing ✅ (2026-02-20)

**Problem**: `extend_match()` only narrowed `sa_start`, overestimating match lengths for positions at the end of the SA range.

**Fix**: Replaced `extend_match()` with `max_mappable_length()` — binary searches within SA range to narrow both boundaries. Ports STAR's `maxMappableLength` + `compareSeqToGenome` + `findMultRange`.

New functions:
- `compare_seq_to_genome()` — starts from offset l_start, returns (match_len, is_read_greater)
- `max_mappable_length()` — binary search within SA range
- `find_mult_range()` — narrow to exact boundary
- `median_uint2()` — safe median for binary search

Removed anchor fallback in `cluster_seeds()` (no longer needed with accurate SA ranges).

**Impact**: Position **97.4% → 97.9%**. CIGAR **97.6% → 98.3%**. Splice rate **3.0% → 2.2%** (matches STAR!). Shared junctions **57 → 62**. False-positive junctions **3 → 0**.

**Files**: `src/align/seed.rs`, `src/align/stitch.rs`

---

## Phase 16.10: Multi-Transcript DP (MAPQ Fix) ✅ (2026-02-20)

**Problem**: rustar-aligner produced exactly **one transcript per window** from DP stitching. For rDNA reads (chrXII tandem repeats), all ~6 copies landed in a single window. STAR found multiple transcripts → NH=6 → MAPQ=0. rustar-aligner found one → NH=1 → MAPQ=255. 323 reads had inflated MAPQ.

**Fix**: Multi-endpoint DP — collect top-N DP endpoints instead of single best:

1. **`build_transcript_from_endpoint()`** — extracted transcript-building logic (chain traceback → junction optimization → extend → CIGAR → exons → Transcript) into helper function
2. **Top-N endpoint collection** — collect all `(total_score, endpoint_idx)` pairs, sort by score descending
3. **Chain-start dedup** — different endpoints in the same chain produce the same alignment, so skip duplicates via `seen_chain_starts` HashSet
4. **Score-range early termination** — `score < best_score - 1` breaks early (worse endpoints won't survive `outFilterMultimapScoreRange` filtering)
5. **`max_transcripts_per_window`** parameter — `alignTranscriptsPerWindowNmax` (default 100) controls the limit
6. **Chimeric detection** — `stitch_seeds()` convenience wrapper passes `max_transcripts_per_window=1` (chimeric detection doesn't need multi-transcript)

| Metric | Before (16.9) | After (16.10) | Improvement |
|--------|---------------|---------------|-------------|
| MAPQ inflation (rustar-aligner=255, STAR<255) | **323** | **62** | -81% |
| MAPQ agreement | 96.1% | **99.1%** | +3.0pp |
| Position agreement | 97.9% | **97.4%** | see note |
| CIGAR agree (of pos-agree) | 98.3% | **98.5%** | +0.2pp |
| Multi-mapped reads | 355 | **627** | STAR: 661 |
| Splice rate | 2.2% | **1.9%** | 122 false splices fixed |
| Disagreements | 190 | **173** | -17 |
| Shared junctions | 62 | **62** | stable |

**Position agreement note**: The 0.5pp drop is expected — diff-chr disagreements are all MAPQ-tied multi-mappers where both tools pick different locations among equally-valid alternatives. All 100 diff-chr cases show matching MAPQ.

**Splice rate improvement**: 122 reads lost false splices (100kb+ spurious introns from wrong DP endpoint). Their CIGARs now match STAR exactly.

**rDNA MAPQ**: Was 304 reads all MAPQ=255, now 26 at MAPQ=1 + 220 at MAPQ=3 + 65 at MAPQ=255 — closely matches STAR's distribution (29 at 1, 225 at 3, 60 at 255).

**PE impact** (multi-transcript DP improves both SE and PE since each mate uses `align_read()`):

| PE Metric | Before (16.9) | After (16.10) | Improvement |
|-----------|---------------|---------------|-------------|
| Per-mate position agree | 95.6% | **97.8%** | +2.2pp |
| Per-mate CIGAR agree | 93.1% | **96.0%** | +2.9pp |
| Both mates mapped | 8706 (96.6%) | **8761 (97.1%)** | +55 pairs |
| Half-mapped pairs | 311 (3.4%) | **263 (2.9%)** | -48 pairs |
| Shared junctions | 76 | **85** | +9 |
| rustar-aligner-only junctions | — | **3** | — |

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`

---

## Phase 16.11b: Fix extend_alignment() to Match STAR's extendAlign() ✅ (2026-02-20)

**Problem**: rustar-aligner's `extend_alignment()` passed actual cumulative alignment length (e.g., seed length ~20) as `len_prev`, activating the proportional mismatch check (`pMMmax * total_length`). STAR always passes `Lprev=100000`, making only the absolute `nMMmax` limit apply.

**Example**: For a seed of length 20 extending leftward with `pMMmax=0.3`, `nMMmax=10`:
- STAR: allows up to 10 mismatches (proportional check disabled by Lprev=100000)
- rustar-aligner: allows only `min(0.3*(20+i+1), 10)` ≈ 6 mismatches at start

**Fix**:
1. **All 4 call sites**: Changed `len_prev` to `100_000` (pre-DP left, endpoint right, post-DP left, post-DP right)
2. **Post-DP left extension**: Changed `n_mm_prev` from `best_state.n_mismatch` to `0` (matches STAR)
3. **Loop body restructured** to match STAR's `extendAlign()` exactly:
   - Record max score only on match (was: on every base)
   - Break only on mismatch, checking limit before incrementing nMM (was: after)
   - Break condition uses full `max_extend` length (was: current position `i+1`)

**Impact**: Neutral on 10k yeast dataset (all metrics unchanged). The fix ensures correctness for edge cases (very short seeds, high mismatch regions, stricter filter settings) and matches STAR's calling convention exactly.

**Files**: `src/align/stitch.rs`

---

## Phase 16.12: Diagnose and Fix Remaining SE Disagreements ✅ (2026-02-20)

**Problem**: 173 SE disagreements remain. Of these, **100 are diff-chr multi-mapper ties** — reads that map equally well to multiple chromosomes (all MAPQ 1 or 3), where the choice of primary alignment depends on implementation-specific processing order (SA iteration, window sequence). These are unfixable without matching STAR's low-level tie-breaking.

**Adjusted position agreement definition**: Excluding diff-chr multi-mapper ties (both mapped, different chromosome, both MAPQ < 255) from the denominator. These are not meaningful disagreements — both tools report valid alignments, just with different tie-breaking among equally-good locations.

**Infrastructure added**:
1. **`--readNameFilter` parameter** — when set, produces detailed alignment trace on stderr for a specific read: seed counts, cluster details, DP expanded seeds + endpoints + scores, transcript filtering decisions, final output
2. **`compare_sam_thorough.py` updated** — new "Adjusted Agreement" section counting diff-chr ties and reporting adjusted metrics
3. **`extract_disagreement_reads.py`** — parses SAM comparison, categorizes disagreements (false_splice, missed_splice, same_chr_close, same_chr_far, star_only, rustar_aligner_only, diff_chr_tie), extracts representative reads for targeted debugging

**Adjusted metrics**:
- Raw position agreement: 97.4% (8620/8793 both-mapped reads)
- Diff-chr multi-mapper ties: ~100
- **Adjusted position agreement: 99.2%** (excluding ties)
- Actionable disagreements: 73 same-chr + 33 STAR-only + 26 rustar-aligner-only = 132

**Files**: `src/params.rs`, `src/align/read_align.rs`, `src/align/stitch.rs`, `test/compare_sam_thorough.py`, `test/extract_disagreement_reads.py`

---

## Phase 16.26: SA Range Narrowing Fix ✅ (2026-02-XX)

**Problem**: Two bugs in `find_mult_range()` and `max_mappable_length()` in `seed.rs` caused incorrect SA range widths for near-identical tandem repeats (rDNA copies ~9kb apart, chr II/XVI paralogs).

**Bug 1**: `find_mult_range()` returned early when `l1 >= l3`. STAR's `findMultRange` continues searching outward using history variables `(ia, ib)`. Fixed to match STAR: set up search range based on whether `l1 < l3`, `l1a < l1`, or `l1a >= l1`.

**Bug 2**: `max_mappable_length()` unconditionally shifted history (`i1a = i1b; i1b = i1;`). STAR guards with `if (L3>L1)`. Added guard: only shift when `l3 > l1` / `l3 > l2`.

**Impact**: Near-identical tandem repeats get correct SA range width (2 instead of 1) → correct NH → correct MAPQ.
- MAPQ inflation: 23 → 12 (-11)
- STAR-only mapped: 2 → 0
- Actionable disagreements: 46 → 37 (-9)
- Splice rate: 2.1% → 2.2% (matches STAR exactly)

**Files**: `src/align/seed.rs`

---

## Phase 16.27: Reverse-Strand Stitcher Coordinate Fix ✅ (2026-03-02)

**Problem**: For reverse-strand clusters, gap-fill scoring in the recursive stitcher compared read bases against the wrong genome region, inflating false splice scores. E.g., a 563kb false intron scored 126 instead of the correct 124, beating the true unspliced alignment (AS=125).

**Root cause**: STAR stores `WA_gStart` in **forward** genome coordinates (converting via `a1 = nGenome - (aLength + a1)` in `ReadAlign_stitchPieces.cpp`) and uses `Read1[1]` (the RC read) for reverse-strand stitching. Our code used SA positions (RC genome offsets) with the forward read, so gap-fill between seeds was scored against genome adjacent to the wrong seed.

**Fix** (in `stitch_seeds_with_jdb_debug`, `src/align/stitch.rs:~1840`):
1. For reverse-strand clusters: RC the read (`read_seq.iter().rev().map(|&b| 3-b)`)
2. Convert WA entries: `wa.read_pos = read_len - (wa.read_pos + wa.length)` (→ positive-strand read pos) and `wa.sa_pos = wa.genome_pos` (→ forward genome pos)
3. Stitch as if forward strand (`stitch_cluster.is_reverse = false`)
4. Restore `transcript.is_reverse = true` and original `read_seq` in output transcripts

**Impact**:
- All 4 false splice reads fixed (introns of 563kb, 150kb, 139kb, 72bp)
- Actionable disagreements: 37 → 28 (-9, -24%)
- MAPQ inflation: 12 → 5 (-7)
- MAPQ deflation: 1 → 2 (+1, two new cases of extra unspliced secondary)
- Same-chr disagreements: 35 → 26 (-9)
- Adjusted pos agreement: 99.6% → 99.7%
- Splice rate: 2.2% (unchanged)

**Files**: `src/align/stitch.rs`

---

## Phase 16.PE1: Recursive Combinatorial Stitcher ✅

**Problem**: Forward DP stitcher could not produce multiple valid transcripts from a window when seeds had multiple reasonable include/exclude combinations.

**Fix**: Replaced forward DP with a recursive include/exclude stitcher (`stitch_recurse` in `stitch.rs`) modelled on STAR's `stitchWindowAligns`. For each WA entry, explores both INCLUDE branch (stitch onto current transcript via `stitch_align_to_transcript`) and EXCLUDE branch (skip entry). Deduplication via `blocks_overlap`: existing higher-scoring transcripts dominate subsets. Anchor constraint: the last anchor must appear in at least one transcript. MAX_RECURSION=10,000 limit prevents blowup on large windows.

**Files**: `src/align/stitch.rs`

---

## Phase 16.PE2: PE Joint DP — Combined-Read Path ✅

**Problem**: rustar-aligner aligned each PE mate independently then tried to pair. 263 pairs half-mapped because one mate had no seeds — but the mate had seeds when combined with the other in the correct orientation.

**STAR's approach**: Build a combined read `[mate1_fwd | SPACER=11 | RC(mate2_fwd)]`. Align it as a single SE read. The stitcher produces a joint transcript spanning both mates; split it at the mate boundary into two SAM records. RC(combined) = `[mate2_fwd | RC_SPACER | RC(mate1_fwd)]` — the "stitch read" for reverse clusters.

**Key implementation details**:
- `mate_id: u8` field added to `WindowAlignment` and `ExonBlock` (0=mate1, 1=mate2, 2=SE/untagged)
- `assign_seed_mate_ids()` tags seeds by position in combined read
- Mate-boundary gap handling in `stitch_align_to_transcript`: when crossing mates, skip junction scoring, check `alignMatesGapMax` instead
- `split_working_transcript()`: splits joint WT at first mate1 exon, adjusts mate2 read coords by `-(len1 + SPACER_LEN)`
- `find_mate_boundary()`: detects joint WT by requiring both mate0 and mate1 exons present

**Score threshold**: `wt.score < 0.66 * (len1 + SPACER_LEN + len2 - 1)`. Check BEFORE `split_working_transcript` — the split copies `wt.score` to both halves, so checking `wt1.score + wt2.score` would double-count the threshold.

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`

---

## Phase 16.PE3: PE STAR-Faithful Architecture Refactor ✅

**Problem**: rustar-aligner had a hybrid architecture: combined-read DP first, then independent SE alignment of each mate, then cross-product pairing. STAR only uses the combined-read path — no independent SE fallback, no cross-product. This produced ~426 false-positive BothMapped pairs.

**Fix**: Removed Phase 2 (independent SE via `align_read()`) and cross-product entirely. Decision tree:
1. Joint pairs from combined-read path → BothMapped
2. Single-mate clusters → discarded (score < combined threshold)
3. No joint pairs → TooShort / unmapped

Note: `adjust_mate2_coords()` helper in `stitch.rs` shifts mate2-only WT coords from combined-read space `[len1+spacer, ...)` → mate2-local space `[0, len2)`.

**Result**: false-positive BothMapped pairs eliminated. Both-mapped: 8612 → stabilised toward STAR's 8390.

**Files**: `src/align/read_align.rs`, `src/align/stitch.rs`

---

## Phase 16.28: extendAlign EXTEND_ORDER + Float Comparison Fix ✅ (2026-03-11)

**Problem 1** (EXTEND_ORDER): STAR's `extendAlign` uses `EXTEND_ORDER=1` — extends the 5' end of the **read** first. For forward reads, 5' = left (extend left first). For reverse-strand reads, 5' = **right** end (extend the rightmost coordinate first). rustar-aligner always extended left first, so reverse-strand reads extended in the wrong direction.

**STAR source** (`stitchWindowAligns.cpp` lines 23-33):
```cpp
uint vOrder[2];
if (tR->Str==0) { vOrder[0]=0; vOrder[1]=1; } // fwd: extend left first
else            { vOrder[0]=1; vOrder[1]=0; } // rev: extend right first
```

**Fix**: Added `original_is_reverse: bool` to `finalize_transcript`. When true, extend right before left. All callers updated: SE passes `stitch_is_reverse`; PE fwd cluster passes `false`; PE rev cluster passes `cluster_is_reverse`.

**Problem 2** (float comparison): STAR uses `double` for mismatch limit: `min(pMMmax * double(Lprev+L), double(nMMmax))`. rustar-aligner truncated to `u32` first, underestimating the limit when the product < 1.0 → premature break on very short extensions.

**Fix**: Changed `extend_alignment()` to use `f64` throughout.

**Impact**: SE +7 agreements (8793→8800), actionable 28→27, STAR-only 1→0. PE: 8382→8383 both-mapped.

**Files**: `src/align/stitch.rs`

---

## Phase 16.29: STITCH-SJ Extended Right Range Fix ✅ (2026-03-12)

**Problem**: When a junction shift `jr_shift > 0` but `jr_shift ≤ shared` (e.g., jr_shift=1, shared=2), read bases shifted into seed B territory at the donor side were NOT scored against the donor genome. The old condition `jr_shift > shared as i32` (1 > 2 = FALSE) completely missed this case.

**Root cause**: rustar-aligner's `jr_shift` is measured from the **end** of the shared region (at `r_a_end = last_exon.read_end + shared`), while STAR's `jR` is measured from the end of **seed A**. Therefore `STAR_jR = jr_shift + shared`. Extended right triggers when `STAR_jR > rGap` (= `shared`), i.e. `jr_shift > 0`, with `n_extra = STAR_jR - rGap = jr_shift`.

**Old code**:
```rust
if jr_shift > shared as i32 {
    let n_extra = (jr_shift - shared as i32) as usize;
```

**New code**:
```rust
if jr_shift > 0 {
    let n_extra = jr_shift as usize;
```

**Example**: Read `ERR12389696.7850795` (CIGAR `1S48M9138N101M`, NM=2, GT-AG). Score was 144; STAR expects 142. With jr_shift=1, shared=2: read[48] was shifted to donor side but scored as a match (+1) instead of mismatch (−1) vs donor genome, net +2 error. After fix, score is 142 = correct.

**Also reverted**: A previous session applied `transcript.score - 2 * n_junction` as a proxy workaround for the AS tag in `sam.rs`. With the real scoring fix, this workaround over-subtracts. Reverted all 3 occurrences back to `transcript.score`.

**Impact**: MAPQ inflation 5→4, actionable disagreements 27→26. Two pre-existing MAPQ deflation cases now visible (rustar-aligner correctly lowers primary score → extra unspliced secondary now within score-range threshold). MAPQ deflation 2→4.

**Files**: `src/align/stitch.rs`, `src/io/sam.rs`

---

## Phase 16.30: PE Overlap Check Fix ✅ (2026-03-XX)

**Problem**: Forward cluster overlap check 1 used raw `wt1.first.genome_start` which incorrectly rejected valid short-insert FR pairs where mate1 has a left soft clip (raw genome_start > mate2.start, but post-extension they align correctly).

**Fix**: Use `left_start_ext = wt1_first.genome_start - wt1_first.read_start` (post-extension estimate) for check 1. Restored old reverse cluster checks (check 1: left.start > right.start → reject; check 2: left.end > estimated right.end → reject).

**Key insight**: FR false positives come from the FORWARD cluster; RF false positives come from the REVERSE cluster.

**Files**: `src/align/read_align.rs`

---

## Phase 16.31: scoreGenomicLengthLog2scale PE Penalty ✅ (2026-03-25)

**Problem**: 144 rustar-aligner-only false positive pairs (non-palindromic overlapping pairs). STAR rejects these via `outFilterMatchNminOverLread` because the combined WT nMatch < threshold when only one genomic region is covered (zero-insert / very short insert). rustar-aligner's combined WT score was passing the threshold.

**Fix**: Apply `scoreGenomicLengthLog2scale` penalty to the combined WT score before comparing to `combined_score_threshold`. STAR applies this penalty in `stitchWindowAligns.cpp:262-265`: `Score = max(0, Score + penalty)`. For overlapping pairs, the penalty is typically -2, pushing scores from 200 → 198 (below threshold 198).

**Impact**: 132 of 144 FPs fixed. 12 FPs remain (score inflation TBD). 35 new regressions introduced — pairs with raw scores 198-200 that fall below 198 after penalty (latent stitching/scoring discrepancy vs STAR).

**Files**: `src/align/read_align.rs`

---

## Phase 16.32: outFilterMultimapNmax PE Check ✅ (2026-03-26)

**Problem**: STAR's `mappedFilter` rejects PE reads where `nTr > outFilterMultimapNmax` even after score-range filtering. rustar-aligner was missing this check for the PE joint path.

**Fix**: Added check `if joint_pairs.len() > params.out_filter_multimap_nmax → return TooManyLoci` after `filter_paired_transcripts`, mirroring STAR's ordering.

**Status**: Mechanism correct but blocked for all 12 remaining FPs by score inflation (rustar-aligner scores 36-106 pts higher than STAR for those pairs).

**Files**: `src/align/read_align.rs`

---

## Phase 16.D5: PE Junction Consistency Check ✅ (2026-03-16)

**Problem**: When mates overlap in the genome, any splice junction in the overlap region must be present in both mates' alignments. Missing this check produced a small number of false positive overlapping pairs with inconsistent junction placement.

**Fix**: `pe_junctions_consistent(&t1, &t2)` / `pe_junctions_consistent(&t2, &t1)` added to forward and reverse cluster joint paths respectively, after same-chr check and before `check_proper_pair`. Functions were implemented earlier but not wired in.

**Files**: `src/align/read_align.rs`

---

## Phase 16.33: PE Mate2 Left-Extension Suppression ✅ (2026-03-31)

**Problem**: Zero-insert RF pairs (e.g. `ERR12389696.10454315`: mate1 `V:340570 104M46S`, mate2 `V:340570 46S104M`) were not being mapped. Two cascading bugs:

**Bug 1** — `extlen` unsigned arithmetic in `stitch_align_to_transcript`: When DP takes Path B (`{seed86 → mate2}` instead of `{seed0 → seed86 → mate2}`), `wa.sa_pos=3748442 < first_exon.genome_start=3748527`, causing the old unsigned conditional to fall back to `extlen=wa.read_pos=198` (giving a spurious 2-base left extension of mate2 instead of 1 base).

**Fix**: Signed i64 arithmetic matching STAR's `gBstart - EX_G + EX_R` formula:
```rust
let raw = (wa.sa_pos as i64) - (first_exon.genome_start as i64) + (first_exon.read_start as i64);
let extlen = if raw > 0 { (raw as usize).min(wa.read_pos) } else { 0 };
```

**Bug 2** — Spurious left extension in `finalize_transcript` for mate2: After `split_working_transcript`, `wt2.read_start=46`. In per-mate `finalize_transcript`, `alignment_start=46`, so it tries to extend 46 bases left into adapter. The adapter base rc_mate2[45]='T' matched genome[3748440]='T' by chance, giving 1 extra base (rStart=196 instead of 197, gStart=3748440 instead of 3748441), causing the overlap check to reject the pair.

**Fix**: Added `no_left_ext: bool` parameter to `finalize_transcript`. The joint stitcher already establishes mate2's left boundary via `extlen`; per-mate finalize must not re-extend leftward. Pass `true` for mate2 in the forward joint cluster and `rc_mate1` in the reverse joint cluster.

**Impact**:

| Metric | Before | After |
|--------|--------|-------|
| PE both-mapped | 8245 (−145 vs STAR) | **8392 (+2 vs STAR)** |
| Per-mate pos agree | ~98.4% | **98.8%** |
| Per-mate CIGAR agree | ~98.5% | **98.5%** |
| SE results | unchanged | unchanged |
| Tests | 268/268 | 268/268 |

**Files**: `src/align/stitch.rs`, `src/align/read_align.rs`

---

## Remaining SE Issues (Post Phase 16.33)

| Category | Count | Notes |
|----------|-------|-------|
| Diff-chr multi-mapper ties | ~100 | Unavoidable — same score, different repeat copy |
| Same-chr repeat ties | ~19 | Unavoidable — rDNA/segmental dups |
| Wrong intron choice (same chr) | 4 | rustar-aligner picks 170kb intron, STAR picks 84kb |
| MAPQ inflation | 4 | Missed splice/indel secondary |
| MAPQ deflation | 4 | Extra unspliced secondary |
| rustar-aligner-only false splice | 1 | `ERR12389696.18296181` — 279kb intron, adapter contamination |
| STAR-only mapped | 1 | `ERR12389696.13766843` — NM=10, unknown |
| **Total actionable** | **26** | |

## Remaining PE Issues (Post Phase 16.33)

| Category | Count | Notes |
|----------|-------|-------|
| Diff-chr per-mate disagreements | 75/mate | Unavoidable multi-mapper ties |
| Same-chr per-mate disagreements | ~21-24/mate | Some fixable |
| rustar-aligner-only FPs | ~2 pairs | Score inflation (combined WT 36-106 pts above STAR) — root cause TBD |
| rustar-aligner vs STAR total pair gap | +2 rustar-aligner | rustar-aligner maps 2 more pairs than STAR |

**Latent stitching/scoring discrepancy**: 35 pairs regressed from Phase 16.31 (raw combined WT scores 198-200 fall below 198 after penalty). STAR's raw scores for the same pairs are slightly higher. Root cause is a stitching or scoring difference not yet identified.
