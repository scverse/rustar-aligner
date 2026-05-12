[← Back to comparison index](README.md)

# extendAlign.cpp vs extend_alignment()

**STAR file**: `source/extendAlign.cpp`
**rustar-aligner file**: `src/align/stitch.rs`, function `extend_alignment()`
**Purpose**: Extend a seed alignment into unaligned read bases at one end (left or right), stopping when the score would drop too far below the best seen.

---

## STAR Function Signature

```cpp
bool extendAlign(char* R, char* G,
    uint rStart, uint gStart, int dR, int dG,
    uint L, uint Lprev, uint nMMprev, uint nMMmax,
    double pMMmax, bool extendToEnd,
    Transcript* trA)
```

- `L`: maximum extension length
- `Lprev`: cumulative previous alignment length (always passed as 100000 to disable proportional check)
- `nMMprev`: accumulated mismatches so far (always 0 in new calls after Phase 16.11b fix)
- `nMMmax`: absolute mismatch limit (`outFilterMismatchNmaxTotal`)
- `pMMmax`: proportional mismatch limit (`outFilterMismatchNoverLmax`)
- `extendToEnd`: if true, extend to spacer or end regardless of score
- Returns: true if any extension was made

---

## `extendToEnd` Path (Mate Extension)

**STAR**:
```cpp
if (extendToEnd) {
    for (iExt=0; iExt<(int)L; iExt++) {
        // Stop at genome padding (5) or spacer base
        if (G[iG]==5) { trA->extendL=0; trA->maxScore=-999999999; return true; }
        if (R[iS]==MARK_FRAG_SPACER_BASE) break;
        if (R[iS]>3 || G[iG]>3) continue; // skip N
        if (G[iG]==R[iS]) { nMatch++; Score+=scoreMatch; }
        else { nMM++; Score-=scoreMatch; }
    }
    if (iExt>0) { trA->extendL=iExt; trA->maxScore=Score; return true; }
    else return false;
}
```

**Key**: Doesn't track `maxScore` during extension — extends ALL the way to the spacer or length limit, reporting the total score. No score-based stopping.

**rustar-aligner**: Needs verification that the `extendToEnd` path (used for PE mate gap extension in `stitch_align_to_transcript`) is implemented this way.

---

## Normal Extension Path

**STAR**:
```cpp
for (int i=0; i<(int)L; i++) {
    iS = dR*i; iG = dG*i;
    // Stop at genome padding, spacer, or invalid base
    if (G[iG]==5 || R[iS]==MARK_FRAG_SPACER_BASE) break;
    if (R[iS]>3 || G[iG]>3) continue; // skip N
    if (G[iG] == R[iS]) {
        nMatch++; Score += scoreMatch;
        if (Score > trA->maxScore) {          // STRICT GREATER THAN
            if (nMM + nMMprev <= min(pMMmax*double(Lprev+i+1), double(nMMmax))) {
                trA->extendL = i+1;
                trA->maxScore = Score;
                trA->nMatch = nMatch;
                trA->nMM = nMM;
            }
        }
    } else {
        if (nMM + nMMprev >= min(pMMmax*double(Lprev+L), double(nMMmax))) {
            break; // too many mismatches
        }
        nMM++; Score -= scoreMatch;
    }
}
bool extDone = trA->extendL > 0;
return extDone;
```

**Critical details**:
1. **Strict `>`** for maxScore record: only saves a new maximum when `Score` STRICTLY exceeds previous max. On ties, keeps the earlier (shorter) extension.
2. **Mismatch limit for record**: uses `Lprev + i + 1` (current position, 1-based).
3. **Break condition uses `Lprev + L`** (full length!) not current position. This means the mismatch check for BREAKING is less strict than for recording — allows more mismatches to occur before stopping.
4. **Check BEFORE incrementing nMM**: `if nMM + nMMprev >= limit { break }` then `nMM++`. This means the limit is checked inclusive — breaks when AT the limit.
5. **N bases** (`R[iS]>3 || G[iG]>3`): continue without scoring (neither match nor mismatch counted).
6. **Genome padding** (`G[iG]==5`): hard stop (break).

---

## rustar-aligner: extend_alignment()

The function in `stitch.rs` was updated in Phase 16.11b to match STAR exactly. Key aspects to verify:

**Strict `>`**: Record `maxScore` only when strictly greater:
```rust
if score > tr.max_score {
    if mm_check_ok { tr.extend_l = i+1; tr.max_score = score; ... }
}
```

**Break uses full length `L`**:
```rust
let limit = (p_mm_max * (l_prev + L) as f64).min(n_mm_max as f64);
if (n_mm + n_mm_prev) as f64 >= limit { break; }
```

**Float types**: Phase 16.28 fixed this to use `f64` throughout, matching STAR's `double`.

**Lprev = 100000**: All call sites use `100_000`, matching STAR's convention (Phase 16.11b fix).

**Assessment**: After Phase 16.11b and 16.28 fixes, this should be fully equivalent. ✅

---

## Call Sites

### In `finalize_transcript` (left and right extension)

**STAR** (`stitchWindowAligns.cpp`):
```cpp
// Case 0: extend at start (5' of read)
extendAlign(R, G, trA.rStart-1, trA.gStart-1, -1, -1,
    trA.rStart, 100000, 0, nMMmax, pMMmax, alignEndsType, &trAstep1)

// Case 1: extend at end (3' of read)
extendAlign(R, G, tR2+1, tG2+1, +1, +1,
    Lread-tR2-1, 100000, scoreSeedBestMM, nMMmax, pMMmax, alignEndsType, &trAstep1)
```

**rustar-aligner** (`finalize_transcript`): The extension order is controlled by `original_is_reverse` (Phase 16.28):
- Forward reads: extend left first (case 0), then right (case 1)
- Reverse reads: extend right first (case 1), then left (case 0)

The `nMMprev` for the right extension is `trA.nMM` (accumulated mismatches from left extension). In STAR it's `scoreSeedBestMM[iS1]`. For the left extension, STAR uses `nMMprev = 0`.

**Potential issue**: For the right extension, STAR uses `scoreSeedBestMM[chain endpoint]` as `nMMprev`. This is the mismatch count from the DP chain (not from any prior extension). rustar-aligner uses `wt.n_mismatch_total` or similar. This needs verification.

---

## extendToEnd in Cross-Fragment (PE Mate Gap)

**STAR** (`stitchAlignToTranscript.cpp`, cross-fragment case):
```cpp
// Extend A rightward into mate gap (extendToEnd=true, length=DEF_readSeqLengthMax)
extendAlign(R, G, rAend+1, gAend+1, 1, 1, DEF_readSeqLengthMax,
    trA->nMatch, trA->nMM, outFilterMismatchNmaxTotal, ...,
    P.alignEndsType.ext[trA->exons[trA->nExons-1][EX_iFrag]][1], &trExtend);

// Extend B leftward from mate boundary
uint extlen = P.alignEndsType.ext[iFragB][1] ? DEF_readSeqLengthMax
                                              : gBstart-trA->exons[0][EX_G]+trA->exons[0][EX_R];
extendAlign(R, G, rBstart-1, gBstart-1, -1, -1, extlen,
    trA->nMatch, trA->nMM, outFilterMismatchNmaxTotal, ...,
    P.alignEndsType.ext[iFragB][1], &trExtend);
```

The `extendToEnd` flag (`P.alignEndsType.ext[iFrag][...]`) controls whether extension stops at score-max or forces to the boundary.

**rustar-aligner** (`stitch_align_to_transcript`, mate boundary): Does NOT do these inward extensions when processing the mate boundary. Just appends the mate2 seed with `new_wt.score += wa.length`. The final `finalize_transcript` call extends outward (away from mates) not inward.

**Impact**: 🟡 For PE reads, if there's a gap between the mates (fragment size > 2×read_length), the inter-mate region won't be scored in rustar-aligner's intermediate step. The STAR approach may produce slightly better final scores for these reads by extending inward during stitch processing.
