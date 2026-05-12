[← Back to comparison index](README.md)

# stitchWindowAligns.cpp vs stitch_recurse()

**STAR file**: `source/stitchWindowAligns.cpp`
**rustar-aligner files**: `src/align/stitch.rs` — `stitch_recurse()` (recursion) + finalization in `finalize_transcript()` + score-range check in `stitch_seeds_core()`

---

## Overview

STAR's `stitchWindowAligns` is a single monolithic recursive function that handles:
1. The recursive include/exclude branching over WA entries
2. Finalization of complete transcripts (extension, overhang checks, dedup, score-range check, insertion into the sorted list)

rustar-aligner splits this into:
- `stitch_recurse`: recursion + dedup at base case
- `finalize_transcript`: extension (extendAlign), CIGAR building, genomic length penalty
- Score-range filter in `stitch_seeds_core` after all recursion

---

## Recursive Structure

**STAR**:
```cpp
void stitchWindowAligns(uint iA, uint nA, int Score, ...) {
    if (iA >= nA && tR2 == 0) return; // empty transcript
    if (iA >= nA) { // base case: finalize
        // ... full finalization (see below)
        return;
    }
    // Include branch
    int dScore = stitchAlignToTranscript(...);
    if (dScore > -1000000) {
        WAincl[iA] = true;
        stitchWindowAligns(iA+1, nA, Score+dScore, ...updated trAi...);
    }
    // Exclude branch (anchor constraint)
    if (WA[iA][WA_Anchor] != 2 || trA.nAnchor > 0) {
        WAincl[iA] = false;
        stitchWindowAligns(iA+1, nA, Score, ...unchanged trA...);
    }
}
```

**rustar-aligner** (`stitch_recurse`):
```rust
fn stitch_recurse(i_a, wt, wa_entries, ...) {
    if i_a >= wa_entries.len() {
        // base case: dedup + push to transcripts
        return;
    }
    let wa = &wa_entries[i_a];
    // Include branch
    if wt.exons.is_empty() { // first seed
        stitch_recurse(i_a+1, new_wt, ...);
    } else if let Some(new_wt) = stitch_align_to_transcript(&wt, wa, ...) {
        stitch_recurse(i_a+1, new_wt, ...);
    }
    // Exclude branch (anchor constraint)
    if can_exclude {
        stitch_recurse(i_a+1, wt, ...);
    }
}
```

**Assessment**: Structure is equivalent. Key differences:
1. STAR's `dScore > -1000000` vs rustar-aligner's `Some(new_wt)` — both gates on successful stitching. ✅
2. STAR's anchor constraint: `WA_Anchor != 2 || trA.nAnchor > 0`. rustar-aligner: `WA[i].is_anchor && i_a == last_anchor → wt.n_anchor > 0`. This checks only the LAST anchor entry. STAR uses `WA_Anchor == 2` to mark the "last anchor" (value 2 vs 1). Needs verification.
3. STAR initializes the transcript from the FIRST included seed using special code (`trAi.rStart = WA[iA][WA_rStart]`). rustar-aligner's first-seed case creates a `WorkingTranscript` with a single exon and score = `wa.length`. Equivalent.

---

## Score Accumulation

**STAR**: `Score` is passed as a parameter and accumulated during recursion. Each step adds `dScore` from `stitchAlignToTranscript`. The base-case score includes all contributions.

**rustar-aligner**: `wt.score` is accumulated in `stitch_align_to_transcript` — the score is embedded in the `WorkingTranscript` struct rather than passed as a parameter. Equivalent result.

**Key difference**: STAR adds `WA[iA][WA_Length] * scoreMatch` when initializing the first seed:
```cpp
for (uint ii=0; ii<WA[iA][WA_Length]; ii++) dScore += scoreMatch;
```
rustar-aligner does this in the first-seed branch of `stitch_recurse`:
```rust
new_wt.score = wa.length as i32; // scoreMatch = 1
```
`scoreMatch = 1` so these are equivalent. ✅

---

## Base Case: Finalization

### STAR: Finalization INSIDE Recursive Function

```cpp
if (iA >= nA) {
    // 1. Extend (EXTEND_ORDER)
    for (int iOrd=0; iOrd<2; iOrd++) {
        switch (vOrder[iOrd]) {
            case 0: extendAlign at start (5' of read for fwd, 3' for rev)
            case 1: extendAlign at end   (3' of read for fwd, 5' for rev)
        }
    }
    // 2. Chr boundary check (alignSoftClipAtReferenceEnds)
    // 3. Compute rLength, gLength
    // 4. Exon overhang checks (alignSJoverhangMin, alignSJDBoverhangMin)
    // 5. Motif strand consistency check
    // 6. outFilterIntronMotifs check
    // 7. Mate length filter (alignSplicedMateMapLmin)
    // 8. BySJout stage 2 check
    // 9. Mate overlap consistency check (PE)
    // 10. Genomic length penalty (scoreGenomicLengthLog2scale)
    // 11. Compute roStart
    // 12. iFrag / maxScoreMate update
    // 13. Variation adjustment
    // 14. Score range check (global + per-mate)
    // 15. Dedup (blocksOverlap) + sorted insertion into wTr[]
}
```

### rustar-aligner: Finalization DEFERRED to finalize_transcript()

The base case of `stitch_recurse` only does dedup. The full finalization happens in `finalize_transcript()` called from `stitch_seeds_core()` after all recursion is complete.

**Missing in rustar-aligner's base case** (applied later or differently):
- Overhang checks at finalization time (steps 4, 7) — rustar-aligner applies these in `filter_transcripts`
- Motif strand consistency (step 5) — applied in `read_align.rs`
- `outFilterIntronMotifs` (step 6) — applied in `read_align.rs`
- Mate overlap consistency (step 9) — partially applied
- `roStart` (step 11) — computed in SAM writer
- `maxScoreMate` update (step 12) — **NOT implemented** (D1 in DIFFERENCES.md)
- Score range check with per-mate score (step 14) — **NOT implemented** (D1)

---

## Score Range Check and Dedup

**STAR** (inside base case, per-transcript):
```cpp
if (Score + P.outFilterMultimapScoreRange >= wTr[0]->maxScore
 || (trA.iFrag >= 0 && Score + P.outFilterMultimapScoreRange >= RA->maxScoreMate[trA.iFrag])
 || P.pCh.segmentMin > 0) {
    // do dedup check, insert into sorted wTr[]
}
```

The list `wTr[]` is sorted by score descending (ties broken by gLength ascending), bounded by `alignTranscriptsPerWindowNmax`. Dedup uses `blocksOverlap`:
- New transcript subset of existing → discard new (if lower score)
- Existing subset of new → remove existing, insert new
- Overlapping but neither subset → keep both

**rustar-aligner** (base case of `stitch_recurse`):
```rust
let dominated = overlap >= wt_len && existing.score >= wt.score && same_structure;
let remove = overlap >= ex_len && wt.score >= existing.score && same_structure;
```

Uses `swap_remove` (unordered). Then after all recursion, `stitch_seeds_core` applies:
```rust
transcripts.retain(|t| t.score >= max_score - params.out_filter_multimap_score_range);
```

**Key differences**:
1. rustar-aligner's dedup only removes transcripts with the **same structure** (same number of exons). STAR's `blocksOverlap` dedup has no such restriction — it can remove an unspliced transcript if a spliced one covers the same bases with equal or higher score. This could affect multi-transcript output.
2. rustar-aligner applies score-range filter AFTER all finalization (post-extension scores). STAR applies it BEFORE extension in the base case. The extension can change the score, so these may disagree in edge cases.
3. rustar-aligner doesn't check per-mate score.

---

## Genomic Length Penalty

**STAR**:
```cpp
if (P.scoreGenomicLengthLog2scale != 0) {
    Score += int(ceil(log2((double)(trA.exons[trA.nExons-1][EX_G] + trA.exons[trA.nExons-1][EX_L]
                           - trA.exons[0][EX_G])) * P.scoreGenomicLengthLog2scale - 0.5));
    Score = max(0, Score);
}
```

Note: `score = max(0, score)` — floored at zero.

**rustar-aligner** (`finalize_transcript`):
```rust
let length_penalty = scorer.genomic_length_penalty(genomic_span);
let final_score = (adjusted_score + length_penalty).max(0);
```

**Assessment**: ✅ Equivalent. The `max(0, ...)` floor is present in rustar-aligner.

---

## Motif Strand Consistency

**STAR**:
```cpp
if (trA.intronMotifs[1] > 0 && trA.intronMotifs[2] == 0) trA.sjMotifStrand = 1;
else if (trA.intronMotifs[1] == 0 && trA.intronMotifs[2] > 0) trA.sjMotifStrand = 2;
else trA.sjMotifStrand = 0;

if (trA.intronMotifs[1] > 0 && trA.intronMotifs[2] > 0
 && P.outFilterIntronStrands == "RemoveInconsistentStrands") return;

if (sjN > 0 && trA.sjMotifStrand == 0 && P.outSAMstrandField.type == 1) return;
```

**rustar-aligner**: Motif strand consistency is checked in `read_align.rs` (`filter_inconsistent_strand_junctions`). The check for `sjMotifStrand == 0` with `outSAMstrandField = intronMotif` is applied. Needs verification of exact conditions.

---

## STAR's WA_Anchor Values

STAR uses:
- `WA_Anchor = 0`: non-anchor seed
- `WA_Anchor = 1`: anchor seed (can be excluded if transcript already has one)
- `WA_Anchor = 2`: LAST anchor seed (can only be excluded if transcript already has an anchor)

rustar-aligner uses `wa.is_anchor: bool` and `last_anchor_idx = wa_entries.iter().rposition(|wa| wa.is_anchor)`. The constraint "can only skip the last anchor if `wt.n_anchor > 0`" matches STAR's `WA_Anchor == 2` logic. ✅
