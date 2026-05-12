[← Back to comparison index](README.md)

# STAR vs rustar-aligner: Identified Differences

This file is the primary tracking document for differences between STAR's C++ source and rustar-aligner's Rust port. Each difference is assessed for impact on alignment output.

**Status legend**: 🔴 Likely impacting · 🟡 Uncertain / small impact · 🟢 Confirmed equivalent · ✅ Fixed

---

## D1: `stitchWindowAligns` — Score Range Check Uses Per-Mate Score Too 🟢

**STAR** (`stitchWindowAligns.cpp`, finalization block):
```cpp
if ( Score + P.outFilterMultimapScoreRange >= wTr[0]->maxScore
  || ( trA.iFrag >= 0 && Score + P.outFilterMultimapScoreRange >= RA->maxScoreMate[trA.iFrag] )
  || P.pCh.segmentMin > 0) {
    // record transcript
```

The transcript is recorded if the score is within range of the **global** best OR within range of the **per-mate** best (`maxScoreMate[iFrag]`).

**rustar-aligner** (`stitch_seeds_core` in `stitch.rs`):
```rust
transcripts.retain(|t| t.score >= max_score - params.out_filter_multimap_score_range);
```

Only checks against the global best. The per-mate score check (`maxScoreMate`) is absent.

**Investigation result (2026-03-12)**: The per-mate check has a critical constraint: `trA.iFrag >= 0` is true **only for single-mate transcripts**. Joint transcripts (spanning both mates) always receive `iFrag = -1` in STAR's code:
```cpp
if (trA.exons[0][EX_iFrag] != trA.exons[trA.nExons-1][EX_iFrag]) {
    trA.iFrag = -1;  // multi-fragment (joint)
}
```
And `maxScoreMate` is only **updated** by single-mate transcripts:
```cpp
if (trA.iFrag >= 0) {
    if (Score > RA->maxScoreMate[trA.iFrag]) RA->maxScoreMate[trA.iFrag] = Score;
}
```
Therefore:
- The per-mate check **never benefits joint (BothMapped) transcripts**.
- It only helps retain secondary single-mate transcripts (from the combined-read path) when their score falls below the global max set by a high-scoring joint pair.
- In rustar-aligner, single-mate candidates are collected into `mate1_candidates` / `mate2_candidates` but discarded as `TooShort` (STAR-faithful: no rescue path). Since we have 0 half-mapped reads, this is a non-issue.

**Impact**: **No effect on PE both-mapped count or gap.** Would only affect MAPQ for half-mapped reads, of which there are none. **No action needed.**

---

## D2: `stitchWindowAligns` — Overhang Filters Applied at Finalization Time 🟢

**STAR** checks exon overhang minimums INSIDE the recursive function at finalization:
```cpp
// For non-annotated junctions:
if ( trA.exons[isj][EX_L] < P.alignSJoverhangMin + trA.shiftSJ[isj][0]
  || trA.exons[isj+1][EX_L] < P.alignSJoverhangMin + trA.shiftSJ[isj][1] ) return;
// For annotated junctions (sjdb):
if ( ( trA.exons[isj][EX_L] < P.alignSJDBoverhangMin ... ) ...) return;
```

**rustar-aligner**: Implemented in `finalize_transcript` (`stitch.rs` lines 1377–1426). Checks both left and right exon lengths (including extensions) against `align_sj_overhang_min + shift` for non-annotated, and `align_sjdb_overhang_min` for annotated junctions. The last-exon terminal check is covered because `right_exon_len` includes `right_extend.extend_len` when `isj+1 == wt.exons.len()-1`, matching STAR's `EX_L` after extension. **Confirmed equivalent.**

---

## D3: `stitchWindowAligns` — Soft-Clip at Reference Ends Check 🟢

**STAR**:
```cpp
if (!P.alignSoftClipAtReferenceEnds.yes &&
    ( (trA.exons[trA.nExons-1][EX_G] + Lread - trA.exons[trA.nExons-1][EX_R]) >
       (mapGen.chrStart[trA.Chr] + mapGen.chrLength[trA.Chr])
    || trA.exons[0][EX_G] < (mapGen.chrStart[trA.Chr] + trA.exons[0][EX_R]) ) ) {
    return; //no soft clipping past the ends of the chromosome
}
```

Prevents alignments that would require soft-clipping past chromosome boundaries when `alignSoftClipAtReferenceEnds` is disabled.

**rustar-aligner**: The `alignSoftClipAtReferenceEnds` parameter is not implemented. Since the STAR default is `Yes` (soft-clip at reference ends IS allowed), rustar-aligner's behavior matches the STAR default. Only diverges if the user explicitly sets `--alignSoftClipAtReferenceEnds No`.

**Impact**: Low — default STAR behavior is to allow soft-clipping at ends (`Yes`). Only matters if a user explicitly disables this, which is uncommon. For standard RNA-seq benchmarking, this is equivalent.

---

## D4: `stitchWindowAligns` — `roStart` Computation 🟢

**STAR**:
```cpp
trA.roStart = (trA.roStr == 0) ? trA.rStart : Lread - trA.rStart - trA.rLength;
```

`roStart` is the read-oriented start: for forward reads it equals `rStart`; for reverse reads it's measured from the end.

**rustar-aligner**: This is used to compute soft-clip sizes in SAM output (the left soft-clip = `roStart`, right soft-clip = `Lread - roStart - rLength`). rustar-aligner computes this directly in the SAM writer from CIGAR ops. **Likely equivalent.**

---

## D5: `stitchWindowAligns` — Mate Pair Overlap Consistency Check 🟡

**STAR** (`stitchWindowAligns.cpp` lines ~175-220):
```cpp
if (trA.exons[0][EX_iFrag] != trA.exons[trA.nExons-1][EX_iFrag]) {
    // Joint transcript: first and last exon are from different mates

    // 1. Negative insert size check
    if (trA.exons[trA.nExons-1][EX_G] + trA.exons[trA.nExons-1][EX_L] <= trA.exons[0][EX_G]) return;

    // 2. Overlap check: find mate boundary index (iExFrag1 = first mate2 exon)
    // mate1_end_g = exons[iExFrag1-1][EX_G] + exons[iExFrag1-1][EX_L]
    // mate2_start_g = exons[iExFrag1][EX_G]

    if (mate1_end_g > mate2_start_g) {
        // Mates overlap in genome — verify junction consistency in overlap region

        // Check mate1 junctions that fall within the overlap
        for each mate1 junction (iex1 < iExFrag1-1):
            if junction acceptor (exons[iex1+1][EX_G]) >= mate2_start_g:
                // Must find matching junction in mate2 with same canonSJ motif
                if not found in mate2 → return (reject)
                if found but different motif → return (reject)

        // Check mate2 junctions that fall within the overlap
        for each mate2 junction (iex2 >= iExFrag1):
            if junction donor (exons[iex2][EX_G] + exons[iex2][EX_L]) <= mate1_end_g:
                // Must find matching junction in mate1
                if not found in mate1 → return (reject)
    }
}
```

**rustar-aligner** (`read_align.rs`, forward cluster path ~line 742):
- **Check 1 (negative insert size)**: ✅ Implemented — `t2.genome_end <= t1.genome_start → continue`
- **Check 2 (post-extension overlap start/end)**: ✅ Implemented — the Phase 16.30 fix (uses `left_start_ext` estimate)
- **Junction consistency in overlap region**: ❌ **NOT implemented**

rustar-aligner currently applies the geometric overlap checks but does not verify that splice junctions in the overlapping genome region are consistent between the two mates.

**Investigation result (2026-03-12)**:
D5's junction consistency check is a **false-positive filter** (rejects biologically impossible alignments), not a recovery mechanism for missed pairs. Its impact on the PE gap is:
- It would reduce the **144 rustar-aligner-only false positives** (pairs that pass rustar-aligner but STAR rejects)
- It would **not** help recover the **149 STAR-only missed pairs**
- Net effect on both-mapped gap: **widening** (fewer rustar-aligner-only → larger net gap), unless some rejected false positives happen to compete with (and block) the correct alignment

**Scenarios where D5 matters**:
1. Short-insert pairs (mates overlap) where one mate calls a splice junction in the overlap region and the other doesn't
2. Short-insert pairs where both mates call junctions at the same genomic position but with different motifs

These are relatively rare in yeast (most splicing is clean GT-AG). In human data with abundant splicing and short-insert libraries, D5 would matter more.

**Priority**: Medium for correctness (prevents biologically invalid outputs), but not a gap-recovery mechanism. Implementing it may widen the already-small 5-pair gap slightly.

---

## D6: `stitchAlignToTranscript` — Left Shift Scoring Range 🔴

**STAR** (`stitchAlignToTranscript.cpp`): The jR scoring loop runs from `min(1, jR+1)` to `max(rGap, jR)`:

```cpp
for (int ii=min(1,jR+1); ii<=max(rGap,jR); ii++) {
    uint g1 = (ii <= jR) ? (gAend+ii) : (gBstart1+ii);
    if (G[g1]<4 && R[rAend+ii]<4) {
        if ( R[rAend+ii] == G[g1] ) {
            if (ii>=1 && ii<=rGap) { Score+=scoreMatch; nMatch++; }
        } else {
            Score -= scoreMatch; nMM++;
            if (ii<1 || ii>rGap) { Score -= scoreMatch; nMatch--; }
        };
    };
}
```

When `jR < 0` (junction shifted LEFT into seed A territory, `ii` starts at `jR+1 < 1`):
- For `ii < 1`: genome side is `gAend+ii` (still donor), but `ii < 1 || ii > rGap` so any mismatch subtracts an EXTRA `scoreMatch` to cancel the previously-assumed match score.
- This handles cases where the donor exon loses some bases to the intron due to left shift.

**rustar-aligner** (`stitch_align_to_transcript` in `stitch.rs`):
The "extended left range" block handles `jr_shift < -(shared)` — i.e., when the shift goes PAST all shared bases into exon A. But it does NOT handle the case where `0 > jr_shift > -(shared)` (left-shifted within the shared region but not past it).

Specifically, when `jr_shift < 0` but `|jr_shift| <= shared`:
- STAR: `jR = jr_shift + shared` which is still `>= 0` but `< shared`, so `ii` runs from `min(1, jR+1)` to `max(rGap, jR)` = `max(shared, jR)` = `shared`. This correctly scores all shared bases using the appropriate genome side.
- rustar-aligner: The `junction_offset = (shared + jr_shift).max(0).min(shared)` determines how many shared bases go to donor side. The scoring splits shared bases into donor-side and acceptor-side groups. This **should** be equivalent to STAR's loop, but needs verification for edge cases (e.g., `jr_shift = -1`, `shared = 3`).

**Impact**: Unclear. The logic appears correct at first glance for the in-range left-shift case, but any off-by-one in the range boundary could produce wrong scores. Low priority.

---

## D7: `stitchAlignToTranscript` — Step 1 "Move Left" Uses `scoreStitchSJshift` 🟡

**STAR**: Before scanning right in Step 2, Step 1 moves the initial jR1 leftward while:
```cpp
Score1 + P.scoreStitchSJshift >= 0 && int(trA->exons[trA->nExons-1][EX_L]) + jR1 > 1
```

`scoreStitchSJshift` was used to bias the junction scan starting position. STAR removed this from the score calculation itself (set to 0 in newer versions), but the left-scan step still exists.

**rustar-aligner** (`find_best_junction_position` in `score.rs`): Needs verification that the left scan starting point is computed identically. From Phase 16.3 notes, this was implemented.

**Impact**: Low — `scoreStitchSJshift` defaults to 0 in STAR.

---

## D8: `stitchAlignToTranscript` — jR Scan Upper Bound 🟢

**STAR** scans `jR1` from left limit to `jR1 < rBend - rAend` = `jR1 <= rGap + L - 1`.

This means the scan goes INTO seed B, up to one base before its end. The rightmost position tested is `jR1 = rGap + L - 1`, which places the junction `L-1` bases into seed B.

**rustar-aligner** (`find_best_junction_position` in `score.rs` line 407):
```rust
if jr1 >= r_gap as i32 + next_seed_len as i32 {
    break;
}
```

This breaks when `jr1 >= r_gap + L`, i.e., scans while `jr1 <= r_gap + L - 1`. **Confirmed equivalent to STAR.** ✅

---

## D9: `extendAlign` — `Score > maxScore` (strict) vs `Score >= maxScore` 🟢

**STAR** (`extendAlign.cpp`):
```cpp
if (Score > trA->maxScore) {  // STRICT GREATER THAN
    if (nMM+nMMprev <= min(pMMmax*double(Lprev+i+1), double(nMMmax)) ) {
        trA->extendL = i+1;
        trA->maxScore = Score;
    };
};
```

The record threshold is **strictly greater than**. STAR keeps the first (leftmost/earliest) occurrence of the maximum score.

**rustar-aligner** (`extend_alignment` in `stitch.rs` line 204):
```rust
if score > max_score {
    ...
}
```

**Confirmed: uses strict `>`, matching STAR exactly.** ✅

---

## D10: `stitchWindowSeeds` (forward-DP) — Architecture Difference 🟢

**STAR** has TWO stitching passes per window:
1. `ReadAlign::stitchWindowSeeds()` — forward O(N²) DP selecting the single best chain
2. `stitchWindowAligns()` — recursive include/exclude producing multiple transcripts

The forward DP in `stitchWindowSeeds` selects the "winning chain" for the primary single-transcript output. The result is stored in `trAll[iWrec][0]`. Then `stitchWindowAligns` is called to produce additional transcripts for multi-mapper detection.

**rustar-aligner**: Only uses the recursive `stitch_recurse` approach (equivalent to `stitchWindowAligns`). The pre-DP `stitchWindowSeeds` approach is simulated by Phase 16.7b pre-DP seed extension scoring (which computes left-extension scores before recursion to help select good seed endpoints).

**Investigation result (2026-03-12)**: `stitchWindowSeeds` is gated by `#ifdef COMPILE_FOR_LONG_READS` in STAR's Makefile. Standard short-read STAR (compiled without `COMPILE_FOR_LONG_READS`) **never calls** `stitchWindowSeeds`. The only active stitching pass is `stitchWindowAligns`, which rustar-aligner's `stitch_recurse` correctly implements. D10 is **not applicable** for short-read STAR alignment.

**Impact**: None for standard short-read use. **Closed.**

---

## D11: `ReadAlign_mappedFilter` — `outFilterScoreMinOverLread` uses `Lread-1` 🟢

**STAR** (`ReadAlign_mappedFilter.cpp`):
```cpp
else if ( (trBest->maxScore < P.outFilterScoreMin)
        || (trBest->maxScore < (intScore)(P.outFilterScoreMinOverLread * (Lread-1)))
        || (trBest->nMatch < P.outFilterMatchNmin)
        || (trBest->nMatch < (uint)(P.outFilterMatchNminOverLread * (Lread-1))) ) {
    unmapType = 1; // TooShort
```

Uses `(Lread-1)` not `Lread` for the proportional thresholds.

**rustar-aligner** (`read_align.rs`, `filter_transcripts`): Needs verification that `Lread-1` is used. Phase 16.12 docs mention "Lread-1 filter fix" so this may already be fixed.

**Impact**: Off-by-one in filter threshold. Rarely matters (edge case for very short reads or very tight score thresholds).

---

## D12: Insertion jR Scanning — `alignInsertionFlush` 🟡

**STAR** (`stitchAlignToTranscript.cpp`, insertion path):
```cpp
if (P.alignInsertionFlush.flushRight) {
    for (; jR < (int)rBend-(int)rAend-(int)Ins; jR++) {
        if (R[rAend+jR+1] != G[gAend+jR+1] || G[gAend+jR+1]==4) break;
    };
    if (jR == (int)rBend-(int)rAend-(int)Ins) return -1000009;
}
```

When `alignInsertionFlush = flushRight`, STAR slides the insertion as far right as possible.

**rustar-aligner**: The `alignInsertionFlush` parameter and its application in the insertion jR scan may or may not be implemented. Needs verification.

**Impact**: Low — affects insertion placement within reads, which is rare in RNA-seq (short reads mostly have small indels from real variants or sequencing errors).

---

## D13: Deletion Scoring — Short Deletions vs Introns 🟢

**STAR**: Uses `Del * scoreDelBase + scoreDelOpen` for short deletions (`Del < alignIntronMin`) and `scoreGap + jPen` for introns.

**rustar-aligner** (`stitch_align_to_transcript`): Uses `score_del_open + score_del_base * del` for deletions. This matches STAR. For introns, uses `motif_score` (= `scoreGap + jPen`). **Confirmed equivalent.**

---

## D14: `stitchWindowAligns` — `nUnique` and `nAnchor` Tracking 🟡

**STAR**:
```cpp
if ( WA[iA][WA_Nrep]==1 ) trAi.nUnique++;  // unique piece
if ( WA[iA][WA_Anchor]>0 ) trAi.nAnchor++;  // anchor piece
```

`nUnique` counts how many seeds have SA range = 1 (unique in genome). `nAnchor` counts anchor seeds. These are used in per-window scoring and filtering.

**rustar-aligner**: `wt.n_anchor` is tracked. `n_unique` (equivalent of `nUnique`) may not be tracked. This could affect score-range comparisons if STAR uses `nUnique` anywhere.

**Impact**: Low — `nUnique` is mainly used for MAPQ and multi-mapper reporting, which rustar-aligner handles separately.

---

## D15: `stitchWindowAligns` — Variation Adjustment (SNP-aware alignment) 🟢

**STAR** (finalization):
```cpp
Score += trA.variationAdjust(mapGen, R);
```

Adjusts score for known variants (SNPs). Used when variant databases are provided.

**rustar-aligner**: No variant database support. `variationAdjust` always returns 0 in a standard STAR run (no variant DB). **Equivalent for standard use.**

---

## D16: `stitchWindowAligns` — Chimeric Detection Score Pass-Through 🟡

**STAR**:
```cpp
|| P.pCh.segmentMin > 0)  // include even if out of score range when chimeric detection active
```

When `chimSegmentMin > 0`, STAR records ALL transcripts regardless of score (for chimeric detection tier 2).

**rustar-aligner** (`chimeric/detect.rs`): The chimeric stitcher uses `max_transcripts_per_window=1`, not the full multi-transcript list. Tier 2 chimeric detection uses a separate seed-cluster path. This may not perfectly mirror STAR's approach.

**Impact**: Chimeric detection may miss some edge cases. PE chimeric detection is not yet implemented (Phase 17.3).

---

## D17: PE Palindromic False Positives — Cat A 🟡

**Pattern** (investigated 2026-03-12): 121 high-MAPQ rustar-aligner-only pairs have both mates mapping to the same genomic position with large overlap (63-149bp). All have complementary soft clips (mate1 left-clip + mate2 right-clip, or vice versa). STAR does NOT output these.

**STAR behavior**: STAR correctly rejects them. STAR maps valid zero-insert pairs with the same geometry (`150M|150M`, `1S149M|149M1S`), so neither "same-position" nor "soft clips" alone is the criterion. Most likely: STAR's per-mate SE scoring finds a better non-palindromic alignment for each mate individually (adapter bases hit a different genome region), making the palindromic joint alignment sub-optimal. STAR then fails to form a concordant PE pair at those different positions → both unmapped.

**rustar-aligner behavior**: The combined-read approach `[mate1_fwd | SPACER | RC(mate2_fwd)] = [X | SPACER | X]` has all seeds mapping to position P. rustar-aligner finds a valid joint transcript at P → false BothMapped.

**Impact**: 121 extra false-positive pairs in rustar-aligner. Fixing these alone would WORSEN the 8-pair net gap (removes 121 from BothMapped without recovering any STAR-only missed). Must be fixed in conjunction with recovering STAR-only misses.

**Fix needed**: Implement per-mate score comparison (check if individual mate SE scores would produce a better overall result than the joint alignment). This is a significant architectural addition.

---

## D18: PE 151 STAR-Only Missed Pairs 🔴

**Pattern**: 151 pairs that STAR maps as BothMapped but rustar-aligner doesn't find any joint transcript for.

**Root cause**: Unknown. D10 (`stitchWindowSeeds`) is ruled out (STARlong-only). Diagonal dedup mate guard fix had no impact. Candidates:
- Missing seeds in some clusters (seed finding gap)
- Score threshold too strict at mate boundary
- Different seed extension behavior near the mate boundary

**Impact**: Primary contributor to the 8-pair net gap. Recovering these is the highest priority PE improvement.

**Investigation needed**: Sample 10-20 STAR-only missed pairs, examine what alignments STAR produces vs what rustar-aligner sees in its seed clusters. Run rustar-aligner with debug output for specific read names.

---

## Summary Table

| ID | Component | Description | Impact | Status |
|----|-----------|-------------|--------|--------|
| D1 | stitchWindowAligns | Per-mate score check missing | 🟢 | ✅ No effect on BothMapped — only affects half-mapped MAPQ (0 half-mapped) |
| D2 | stitchWindowAligns | Terminal exon overhang check | 🟢 | ✅ Confirmed equivalent |
| D3 | stitchWindowAligns | Soft-clip at chr ends | 🟢 | Default matches STAR |
| D4 | stitchWindowAligns | `roStart` computation | 🟢 | Equivalent |
| D5 | stitchWindowAligns | PE mate overlap junction consistency | 🟢 | ✅ Implemented 2026-03-12 — false-positive filter, may widen net gap slightly |
| D6 | stitchAlignToTranscript | Left-shift scoring range (in-range) | 🟢 | ✅ Confirmed equivalent |
| D7 | stitchAlignToTranscript | Step 1 left-scan with scoreStitchSJshift | 🟡 | Low risk |
| D8 | stitchAlignToTranscript | jR scan upper bound | 🟢 | ✅ Confirmed correct |
| D9 | extendAlign | Strict `>` vs `>=` for maxScore record | 🟢 | ✅ Confirmed correct |
| D10 | stitchWindowSeeds | Forward-DP not used in rustar-aligner | 🟢 | ✅ STARlong-only (#ifdef COMPILE_FOR_LONG_READS) — not applicable for short reads |
| D17 | PE palindromic pairs | 121 false BothMapped from palindromic reads | 🟡 | Root cause identified, fix requires per-mate score comparison |
| D18 | PE missed pairs | 151 STAR-only pairs not found by rustar-aligner | 🔴 | Root cause unknown — highest priority |
| D11 | mappedFilter | `Lread-1` in proportional filter | 🟢 | Fixed (16.12 docs) |
| D12 | stitchAlignToTranscript | `alignInsertionFlush` | 🟡 | Verify |
| D13 | stitchAlignToTranscript | Del vs intron scoring | 🟢 | Equivalent |
| D14 | stitchWindowAligns | `nUnique` tracking | 🟡 | Low risk |
| D15 | stitchWindowAligns | Variation adjustment | 🟢 | Equivalent |
| D16 | Chimeric | Score pass-through for chimeric | 🟡 | Low risk |

---

## Priority Investigation Order

### Remaining — High Priority
1. **D7**: Left-scan starting point in `find_best_junction_position` — low risk since `scoreStitchSJshift = 0`, but worth verifying.

### Remaining — Medium Priority
2. **D5**: PE mate overlap junction consistency — ✅ Implemented 2026-03-12. Reduces rustar-aligner-only false positives; may widen net gap slightly.

### Remaining — Low Priority
4. **D12**: `alignInsertionFlush` for insertion placement.
5. **D14, D16**: Minor issues for nUnique tracking and chimeric pass-through.

### Closed
- **D1**: Investigated 2026-03-12. Per-mate check (`maxScoreMate`) only applies to single-mate transcripts (`iFrag >= 0`). Joint transcripts have `iFrag = -1` and are unaffected. No impact on BothMapped count. **No action needed.**
- **D2, D3, D6, D8, D9**: All confirmed equivalent to STAR.
- **D4, D11, D13, D15**: Always equivalent.
- **D10**: Investigated 2026-03-12. `stitchWindowSeeds` is `#ifdef COMPILE_FOR_LONG_READS` only — standard STAR never calls it. **Not applicable.**
