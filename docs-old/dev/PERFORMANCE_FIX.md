# Performance Fix: Clustering Explosion Bug

**Date**: 2026-02-07
**Issue**: Critical performance bottleneck in seed clustering algorithm
**Status**: ✅ FIXED

---

## Problem Description

The seed clustering algorithm in `src/align/stitch.rs` was creating a **combinatorial explosion** of clusters when reads had no good anchor seeds, causing 10,000x slowdowns.

### Root Cause

When no seeds had ≤10 genomic locations (good anchors), the algorithm fell back to using **all seeds as anchors**. For each anchor, it created a **separate cluster for every genomic position** the seed mapped to.

**Failure mode example (read `ERR12389696.4967095`):**
- 143 seeds found (normal)
- All 143 became anchors (no good anchors found)
- Average ~84 genomic positions per anchor
- **Result: 143 × 84 = 11,982 clusters** (pathological)

---

## Performance Impact

### Before Fix (Yeast RNA-seq, 100 reads)

| Metric | Time | Notes |
|--------|------|-------|
| Pathological read | 8.84s | Single read! |
| - Clustering | 1.937s | 11,982 clusters created |
| - Stitching | 6.857s | Processing all clusters |
| **Total (100 reads)** | **>60s** | **Never completed** |

### After Fix

| Metric | Time | Improvement |
|--------|------|-------------|
| Pathological read | 71.2ms | **124x faster** |
| - Clustering | 37.9ms | 193 clusters (98.4% reduction) |
| - Stitching | 32.2ms | **213x faster** |
| **Total (100 reads)** | **0.233s** | **>250x faster** |

---

## Solution: STAR-Compatible Parameter Limits

Implemented three STAR parameters to prevent clustering explosion:

### 1. `--seedMultimapNmax` (default 10000)
- Seeds mapping to >10000 loci are discarded completely
- Filters highly repetitive sequences

### 2. `--winAnchorMultimapNmax` (default 50)
- Anchors can map to maximum 50 genomic loci
- Skips anchors with too many positions

### 3. `--seedNoneLociPerWindow` (default 10)
- Maximum 10 seed positions considered per window
- Limits cluster creation per anchor

### Additional Improvements

When no good anchors found:
- **Before**: Used all 143 seeds as anchors
- **After**: Sort by SA range size, use best 20 seeds only

---

## Files Modified

### Core Implementation
- `src/params.rs` — Added 3 new STAR parameters
- `src/align/stitch.rs` — Updated `cluster_seeds()` with limits
- `src/align/read_align.rs` — Pass parameters to clustering

### Tests
- `src/params.rs` — Updated defaults test
- `src/align/stitch.rs` — Updated test calls

---

## Test Results

✅ **170/170 unit tests passing**
✅ **3 non-critical clippy warnings** (pre-existing)
✅ **Build: success**

### Known Issue: Integration Tests

The Phase 9 threading integration tests fail because the test genome is pathologically repetitive:
- 50 exact copies of the same 20bp pattern
- Every read maps to 50+ locations
- Now correctly filtered by `winAnchorMultimapNmax=50`

**Resolution**: Integration tests need more realistic test genomes (Phase 13)

---

## References

**STAR Parameters:**
- [STAR parametersDefault](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault)
- [STAR winAnchorMultimapNmax discussion](https://biostar.galaxyproject.org/p/27672/)
- [STAR GitHub](https://github.com/alexdobin/STAR)

---

## Next Steps (Phase 13)

1. Fix integration tests with realistic genomes
2. Profile remaining bottlenecks (100K reads still ~2-3min)
3. Implement full suite of STAR seed parameters:
   - `seedPerReadNmax` (max 1000 seeds per read)
   - `seedPerWindowNmax` (max 50 seeds per window)
   - `seedSearchLmax`, `seedSplitMin`, etc.
4. Optimize hot paths identified by profiling
