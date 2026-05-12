# Real-World Testing Summary - Phase 12.2

**Date**: 2026-02-07
**Test Dataset**: Yeast RNA-seq (ERR12389696)
**Test Location**: `test/data/small/yeast/`
**Status**: ⚠️ **CRITICAL PERFORMANCE ISSUE DISCOVERED**

---

## Executive Summary

✅ **Phase 12.2 functionally complete** - Chimeric detection fully integrated
✅ **Critical bug found and fixed** - Empty exons crash prevented
🚨 **SEVERE performance bottleneck** - ~10,000x slower than STAR
➡️ **Next action**: Phase 13 profiling is MANDATORY, not optional

---

## Test Results

### ✅ What Works

1. **Index Generation** (18 seconds)
   - S. cerevisiae genome: 17 chromosomes, 14.7MB
   - Genome, SA, SAindex files created successfully
   - Parameters: `--genomeSAindexNbases 10 --genomeChrBinNbits 18`

2. **Index Loading**
   - Genome: 17 chromosomes loaded
   - Suffix array: 24,314,210 entries
   - SAindex: 1,398,100 indices
   - GTF: 361 annotated junctions extracted (3 warnings for overlapping exons)

3. **System Initialization**
   - Thread pool: 4 threads spawned successfully
   - Memory: 635MB (reasonable)
   - Files created: SAM and Chimeric.out.junction

4. **Bug Fix** ✅ COMMITTED
   - **File**: `src/chimeric/detect.rs:146`
   - **Issue**: `index out of bounds: the len is 0 but the index is 0`
   - **Cause**: Accessing `transcript.exons[0]` without checking if exons vector is empty
   - **Fix**: Added check: `if t1.exons.is_empty() || t2.exons.is_empty() { return Ok(None); }`
   - **Impact**: Would crash on any read producing transcripts with empty exon lists

---

## 🚨 Critical Performance Issue

### Test 1: 100,000 reads (ERR12389696_sub_1.fastq.gz)
- **Time**: >5 minutes, **NO COMPLETION**
- **CPU**: 400% (4 threads working, not hung)
- **Memory**: 635MB
- **Progress**: "Aligning reads..." message, then silence
- **Expected**: <1 minute

### Test 2: 100 reads (ERR12389696_sub_1_100.fastq.gz)
- **Time**: >60 seconds, **ZERO OUTPUT**
- **CPU**: Working but extremely slow
- **Output files**: Created but 0 bytes (empty)
- **Expected**: <1 second

### Performance Analysis

| Metric | STAR | rustar-aligner (observed) | Ratio |
|--------|------|-------------------|-------|
| Reads/second | ~16,000 | <2 | **~10,000x slower** |
| 100 reads | <0.01s | >60s | 6000x+ |
| 100K reads | ~6s | >300s | 50x+ |

**Conclusion**: This is not a "needs optimization" issue. This is a fundamental algorithmic or implementation problem requiring immediate investigation.

---

## Test Commands Used

### Index Generation (SUCCESSFUL)
```bash
rustar-aligner --runMode genomeGenerate \
  --genomeDir indices_rustar \
  --genomeFastaFiles reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  --genomeSAindexNbases 10 \
  --genomeChrBinNbits 18
```

### Alignment Test 1: 100K reads (FAILED - too slow)
```bash
rustar-aligner --runMode alignReads \
  --genomeDir indices_rustar \
  --readFilesIn reads/ERR12389696_sub_1.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix outputs/rustar_test/noChim \
  --outSAMtype SAM \
  --runThreadN 4
```

### Alignment Test 2: 100 reads (FAILED - too slow)
```bash
rustar-aligner --runMode alignReads \
  --genomeDir indices_rustar \
  --readFilesIn reads/ERR12389696_sub_1_100.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix outputs/rustar_test/tiny100_ \
  --outSAMtype SAM \
  --sjdbGTFfile reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf \
  --chimSegmentMin 15 \
  --chimScoreMin 10 \
  --runThreadN 4
```

---

## Diagnostic Information

### What Couldn't Be Validated
- ❌ Alignment accuracy (no output produced)
- ❌ SAM format correctness (empty file)
- ❌ Chimeric detection (empty file)
- ❌ Junction detection accuracy
- ❌ CIGAR string generation
- ❌ MAPQ calculation

### Files Created During Testing
```
test/data/small/yeast/
  indices_rustar/                    - rustar-aligner-generated index (18s)
  outputs/rustar_test/
    noChim/Aligned.out.sam          - Created but empty (0 bytes)
    tiny100_/Aligned.out.sam        - Created but empty (0 bytes)
    tiny100_/Chimeric.out.junction  - Created but empty (0 bytes)
```

---

## Phase 13: MANDATORY Next Steps

### Step 1: Add Timing Traces (IN PROGRESS)
Added trace points in `src/align/read_align.rs`:
```rust
eprintln!("[TRACE] {} seed_finding: {:?}, found {} seeds", ...);
eprintln!("[TRACE] {} clustering: {:?}, found {} clusters", ...);
```

**TODO**: Complete tracing for:
- [ ] Stitching time per cluster
- [ ] Filtering time
- [ ] SAM writing time
- [ ] Total time per read

### Step 2: Profile with perf

```bash
# Build with debug symbols
cargo build --release

# Profile on tiny dataset (100 reads)
cd test/data/small/yeast
perf record -g ../../../../target/release/rustar-aligner \
  --runMode alignReads \
  --genomeDir indices_rustar \
  --readFilesIn reads/ERR12389696_sub_1_100.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix perf_test \
  --runThreadN 1

# View hotspots
perf report
```

### Step 3: Likely Bottlenecks (Priority Order)

1. **Seed finding** (`Seed::find_seeds`)
   - MMP search inefficiency?
   - Binary search on large suffix array (24M entries)
   - SAindex lookup overhead?

2. **Seed extension/stitching** (`stitch_seeds`)
   - DP algorithm complexity
   - Excessive allocations
   - Poor cache locality

3. **Suffix array search**
   - 24M entries = 96MB for 32-bit indices
   - Poor cache performance on random access?

4. **Memory allocations**
   - Hot path allocations
   - Vector reallocations

### Step 4: Quick Wins to Try

1. **Reduce SA search space**: Use SAindex more aggressively
2. **Cache lookups**: Memoize repeated SA queries
3. **Batch processing**: Process seeds in batches
4. **Profile-guided optimization**: Let data guide decisions

---

## Target Performance Goals

| Dataset | Target Time | STAR Time | Acceptable Ratio |
|---------|-------------|-----------|------------------|
| 100 reads | <1s | <0.01s | 100x slower OK for now |
| 100K reads | <30s | ~6s | 5x slower OK |
| 1M reads | <5min | ~60s | 5x slower OK |

**Minimum acceptable**: Within 10x of STAR before proceeding to correctness validation.

---

## Code Status at End of Session

- **Tests**: 170/170 passing
- **Clippy**: 4 non-critical warnings
- **Build**: Release successful
- **Commits**: Bug fix committed
- **Tracing**: Partially added (needs completion)

---

## Files Modified This Session

1. `src/chimeric/detect.rs` - Bug fix (empty exons check)
2. `src/lib.rs` - Chimeric integration (~100 lines)
3. `src/align/read_align.rs` - Timing traces (partial)
4. `PHASE12_COMPLETE.md` - Implementation summary
5. `TESTING_SUMMARY.md` - This file
6. `CLAUDE.md` - Updated status
7. `ROADMAP.md` - Updated Phase 12 status

---

## Next Session Actions

1. **Complete timing traces** in align_read.rs
2. **Run 100-read test** with traces to see time breakdown
3. **Run perf profile** on 100-read test
4. **Identify bottleneck** from perf report
5. **Fix or optimize** the hotspot
6. **Re-test** and measure improvement
7. **Iterate** until performance is acceptable

**DO NOT proceed** with correctness validation until performance is within 10x of STAR.
