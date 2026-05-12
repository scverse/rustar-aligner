[← Back to comparison index](README.md)

# ReadAlign_stitchWindowSeeds.cpp vs pre-DP in stitch.rs

**STAR file**: `source/ReadAlign_stitchWindowSeeds.cpp`
**rustar-aligner equivalent**: Pre-DP seed extension scoring in `stitch_seeds_with_jdb_debug()` (Phase 16.7b)

---

## Overview

`ReadAlign_stitchWindowSeeds` implements a forward O(N²) DP that selects the single best seed chain from a window. It's called in addition to (or as a precursor to) the full recursive `stitchWindowAligns` pass.

**Important**: rustar-aligner does NOT call an equivalent of `stitchWindowSeeds`. Instead it goes directly to the recursive `stitch_recurse`. Phase 16.7b added left-extension pre-scoring (analogous to `stitchWindowSeeds`' left extension for DP scoring), but the full N² chain DP is absent.

---

## STAR's stitchWindowSeeds Algorithm

### Phase 1: Forward DP with Left Extension

```cpp
for (uint iS1=0; iS1<nWA[iW]; iS1++) {
    // Try every previous seed as predecessor
    for (uint iS2=0; iS2<=iS1; iS2++) {
        if (iS2 < iS1) {
            // DP: stitch iS2 → iS1
            score2 = stitchAlignToTranscript(...iS2...iS1...);
            // Check exon length, then:
            if (exonLongEnough && score2 > 0
                && score2 + scoreSeedBest[iS2] > scoreSeedBest[iS1]) {
                scoreSeedBest[iS1] = score2 + scoreSeedBest[iS2];
                scoreSeedBestInd[iS1] = iS2;  // best predecessor
            }
        } else { // iS2 == iS1: base case (this seed starts a chain)
            score2 = WA_Length;
            // Left extension from this seed:
            extendAlign(R, G, rStart-1, gStart-1, -1, -1, rStart, 100000, 0, ..., &trA1);
            score2 += trA1.maxScore;  // Add left extension score
            if (exonLongEnough && score2 > scoreSeedBest[iS1]) {
                scoreSeedBest[iS1] = score2;
                scoreSeedBestInd[iS1] = iS1;  // self (chain start)
            }
        }
    }
}
```

### Phase 2: Select Best Chain Endpoint via Right Extension

```cpp
intScore scoreBest = 0;
for (uint iS1=0; iS1<nWA[iW]; iS1++) {
    // Right extension from this seed:
    extendAlign(R, G, tR2, tG2, +1, +1, Lread-tR2, 100000, scoreSeedBestMM[iS1], ..., &trA1);
    scoreSeedBest[iS1] += trA1.maxScore;
    if (exonLongEnough && scoreSeedBest[iS1] > scoreBest) {
        scoreBest = scoreSeedBest[iS1];
        scoreBestInd = iS1;  // best chain endpoint
    }
}
```

### Phase 3: Trace Back the Chain

```cpp
uint seedN = 0;
while (true) {
    seedChain[seedN++] = scoreBestInd;
    if (scoreBestInd > scoreSeedBestInd[scoreBestInd]) {
        scoreBestInd = scoreSeedBestInd[scoreBestInd];
    } else {
        break; // self-referencing = chain start
    }
}
```

### Phase 4: Build Transcript from Chain

```cpp
// Process seeds in chain order (seedN-1 to 0)
// then extend left and right
```

### Output

```cpp
*(trAll[iWrec][0]) = trA; nWinTr[iWrec] = 1;
// OR, if WAexcl != NULL:
*(trAll[iWrec][1]) = trA; nWinTr[iWrec] = 2;
```

---

## rustar-aligner's Phase 16.7b Pre-DP Extension

rustar-aligner approximates `stitchWindowSeeds`'s left-extension scoring:

```rust
// Compute left extension for each expanded seed
let left_ext_scores: Vec<i32> = expanded_seeds.iter()
    .map(|seed| {
        if seed.read_pos == 0 { return 0; }
        // extend_alignment with L = seed.read_pos (max leftward extension)
        match extend_alignment(read_seq, seed.read_pos, seed.genome_pos, ..., false) {
            Some(ext) => ext.score,
            None => 0,
        }
    })
    .collect();

// DP base score includes left extension
score = seed_length + left_ext_scores[i];
```

And for endpoint selection, adds right extension:
```rust
for candidate in endpoints {
    let right_ext = extend_alignment(...) ...;
    if score + right_ext > best_score { best = candidate; }
}
```

---

## Key Differences

### 1. N² DP vs Greedy/Sparse
STAR's forward DP tries ALL pairs `(iS2, iS1)` to find the globally optimal chain. rustar-aligner's Phase 16.7b pre-DP only extends individual seeds; the chain selection in `stitch_recurse` is done recursively.

**Impact**: STAR's N² DP can find chains where the optimal path skips certain seeds. rustar-aligner's recursive stitcher also explores all combinations, so for the actual output this likely produces equivalent results. The difference is mainly in efficiency and in the pre-scoring heuristic.

### 2. `scoreSeedBestMM` Threading
STAR threads the mismatch count through the DP (`scoreSeedBestMM[iS1]`) and uses it as `nMMprev` for the right extension. rustar-aligner doesn't have an equivalent per-seed accumulated mismatch count in the pre-DP.

**Impact**: The right-extension mismatch budget may be slightly different between STAR and rustar-aligner, potentially affecting extension length.

### 3. `exonLongEnough` Check in DP
STAR checks `(WA_Length + trA1.extendL) >= P.alignSJoverhangMin` before accepting a chain step. This prevents very short seeds from starting chains unless they extend enough. rustar-aligner doesn't have this check in the pre-DP.

**Impact**: Low — short seeds are usually rejected by the stitcher or by the overhang filter later.

### 4. `WAexcl` Pass (Two-Pass Stitching)
STAR calls `stitchWindowSeeds(iW, iWrec, WAexcl=NULL)` for the primary path, and then optionally `stitchWindowSeeds(iW, iWrec, WAexcl=waExcl)` with certain seeds excluded. The `WAexcl` array excludes seeds that are already used in a "best transcript" so the second pass finds an ALTERNATIVE transcript.

**rustar-aligner**: The recursive stitcher naturally finds multiple transcripts without needing explicit seed exclusion. This is a fundamentally different but potentially equivalent approach.

---

## Is stitchWindowSeeds Always Called?

**Resolution (2026-03-12)**: `stitchWindowSeeds` is gated by `#ifdef COMPILE_FOR_LONG_READS` in STAR's Makefile. Standard short-read STAR is compiled **without** this flag, so `stitchWindowSeeds` is never compiled in or called. The only active stitching pass is `stitchWindowAligns`.

**Conclusion**: D10 is **not applicable** for standard short-read STAR alignment. rustar-aligner's `stitch_recurse` (equivalent to `stitchWindowAligns`) is the correct and complete equivalent. No implementation needed.
