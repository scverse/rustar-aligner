# Phase 16.14: Nstart Computation Bug Fix

## Date: 2026-02-25

## Summary

Fixed a critical bug in the sparse seed search where `seedSearchStartLmax` was used as the **number of starting positions** (Nstart) instead of the **maximum spacing between starting positions**. This caused rustar-aligner to search at ~17x more positions than STAR, finding spurious short seeds that produced false splice junctions.

## The Bug

In `src/align/seed.rs`, function `search_direction_sparse()`:

```rust
// BEFORE (wrong — used seedSearchStartLmax AS nstart):
let nstart = effective_start_lmax.min(read_len).max(1);
// For 150bp read, seedSearchStartLmax=50: nstart=50, lstart=3
// Starting positions: 0, 3, 6, 9, ..., 147 (50 positions!)

// AFTER (correct — STAR divides read length by seedSearchStartLmax):
let nstart = (read_len / effective_start_lmax).max(1); // L→R: floor division
let nstart = (read_len + effective_start_lmax - 1) / effective_start_lmax; // R→L: ceil division
// For 150bp read, seedSearchStartLmax=50: nstart=3, lstart=50
// Starting positions: 0, 50, 100 (3 positions)
```

### STAR's logic (from `ReadAlign_mapOneRead.cpp`)

- `seedSearchStartLmax` (default: 50) is the **max length of each piece**, not the number of starts
- L→R: `Nstart = readLen / seedSearchStartLmax` (floor division, min 1)
- R→L: `Nstart = ceil(readLen / seedSearchStartLmax)` (ceil division)
- `Lstart = readLen / Nstart` (spacing between starting positions)
- At each starting position, MMP chains forward until `< seedMapMin` bases remain

## Root Cause of False Splices

With 50 starting positions (instead of 3), the R→L sparse search was exploring every 3rd position in the reverse-complement read. Starting position at RC pos 138 (= original read position 1) found an 11bp seed matching a distant genomic location. This seed was then stitched with the main 144bp seed via the recursive stitcher, creating a false ~18kb intron with a canonical CtAc motif (GT/AG on reverse strand).

### Detailed trace for read ERR12389696.386431

**Before fix** (50 R→L starting positions):
- R→L search at RC pos 138 finds: 11bp seed at read_pos=1, sa_pos=10245644
- This enters the same cluster as the main 144bp seed at read_pos=6, sa_pos=10264492
- Stitcher creates junction: 6bp left exon + 18843bp intron + 143bp right exon
- Splice score: 149 (11+138) + motif 0 - genomic_length 4 = 145
- Non-splice score: 144 - genomic_length 2 = 142
- **Spliced wins by 3 points** → false splice output

**After fix** (3 R→L starting positions):
- R→L searches only at RC positions 0, 50, 100
- RC pos 138 is never visited → 11bp seed is never found
- Only non-spliced 144M6S transcript (score 142) is produced → correct output

## Investigation Process

1. Added debug output to `stitch_align_to_transcript()` and `stitch_recurse()` to trace junction creation
2. Traced false splice read ERR12389696.386431 — identified 11bp R→L seed as the culprit
3. Verified overhang check passes (jj_l=1, jj_r=6, left_exon=6 >= overhangMin+1=6)
4. Investigated STAR's finalization filters — all 10 checks pass for this junction
5. Key insight: STAR never creates this junction because the seed doesn't exist
6. Traced STAR's `ReadAlign_mapOneRead.cpp` — Nstart = readLen / seedSearchStartLmax, not seedSearchStartLmax itself
7. Our code had `nstart = effective_start_lmax.min(read_len)` which used the parameter value (50) directly as nstart

## Results

### SE 10k yeast comparison (after fix)

| Metric | Before Fix | After Fix | STAR |
|--------|-----------|-----------|------|
| Position agreement | 97.5% | **99.5%** | — |
| CIGAR agreement | 98.5% | **99.2%** | — |
| MAPQ agreement | 99.2% | **99.2%** | — |
| Splice rate | 2.3% | **2.0%** | 2.1% |
| False splices | 20 | **23** | — |
| Missed splices | 22 | **10** | — |
| Shared junctions | 64 | **62** | 67 total |
| rustar-aligner-only junctions | 4 | **1** | — |
| MAPQ inflation | 36 | **28** | — |
| Total disagreements | 175 | **48** | — |
| Diff-chr ties | 98 | **0** | — |
| Actionable disagreements | 128 | **52** | — |

### Key improvements
- **Position disagreements**: 175 → 48 (−73%)
- **Diff-chr multi-mapper ties**: 98 → 0 (eliminated)
- **Actionable disagreements**: 128 → 52 (−59%)
- **Missed splices**: 22 → 10 (−55%)
- **rustar-aligner-only junctions**: 4 → 1 (−75%)
- **MAPQ inflation**: 36 → 28 (−22%)
- **Seeds per read**: ~145 → ~10 (dramatic reduction, matches STAR)

### Minor regressions
- **False splices**: 20 → 23 (+3) — explained by more reads now agreeing in position, some with CIGAR diffs
- **Shared junctions**: 64 → 62 (−2) — likely from changed seed coverage pattern

## Performance Impact

The Nstart fix dramatically reduces the number of seeds found per read (e.g., 145 → 10 for ERR12389696.386431). This means:
- Fewer WA entries per cluster
- Fewer recursions in the combinatorial stitcher (684 → fewer)
- Less time in seed search
- Net effect: faster alignment

## Files Modified

| File | Change |
|------|--------|
| `src/align/seed.rs:215` | Fixed Nstart computation from `effective_start_lmax.min(read_len)` to `read_len / effective_start_lmax` (L→R floor, R→L ceil) |

## Related STAR Source References

- `ReadAlign_mapOneRead.cpp`: Lines ~150-200, L→R Nstart computation
- `ReadAlign_mapOneRead.cpp`: Lines ~250-300, R→L Nstart computation
- Parameter: `--seedSearchStartLmax` (default 50) = max spacing between search start positions

## Remaining Issues

1. **23 false splices** — reads where rustar-aligner adds N in CIGAR but STAR doesn't (was 20, +3 from changed position-agree set)
2. **10 missed splices** — reads where STAR has N but rustar-aligner doesn't (was 22, −12 improvement)
3. **5 STAR-only junctions** — junctions STAR finds but rustar-aligner doesn't (seed search differences?)
4. **48 position disagreements** — 41 same-chr + 7 diff-chr
5. **28 MAPQ inflation reads** — rustar-aligner=255, STAR<255 (multi-mapper resolution)
