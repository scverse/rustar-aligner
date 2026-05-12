[← Back to comparison index](README.md)

# ReadAlign_mappedFilter.cpp vs Quality Filter in read_align.rs

**STAR file**: `source/ReadAlign_mappedFilter.cpp`
**rustar-aligner file**: `src/align/read_align.rs`, filter logic in `align_read()`

---

## STAR's mappedFilter()

```cpp
void ReadAlign::mappedFilter() {
    unmapType = -1; // mark as mapped

    if (nW == 0) {                      // no good windows
        statsRA.unmappedOther++;
        unmapType = 0;

    } else if (
        (trBest->maxScore < P.outFilterScoreMin)
     || (trBest->maxScore < (intScore)(P.outFilterScoreMinOverLread * (Lread-1)))
     || (trBest->nMatch < P.outFilterMatchNmin)
     || (trBest->nMatch < (uint)(P.outFilterMatchNminOverLread * (Lread-1)))
    ) {
        statsRA.unmappedShort++;
        unmapType = 1;  // TooShort

    } else if (
        (trBest->nMM > outFilterMismatchNmaxTotal)
     || (double(trBest->nMM)/double(trBest->rLength) > P.outFilterMismatchNoverLmax)
    ) {
        statsRA.unmappedMismatch++;
        unmapType = 2;  // TooManyMismatches

    } else if (nTr > P.outFilterMultimapNmax) {
        statsRA.unmappedMulti++;
        unmapType = 3;
    };
}
```

---

## Key Parameters

| STAR parameter | Default | Purpose |
|---------------|---------|---------|
| `outFilterScoreMin` | 0 | Absolute minimum score |
| `outFilterScoreMinOverLread` | 0.66 | Proportional min score: `score >= 0.66 * (Lread-1)` |
| `outFilterMatchNmin` | 0 | Absolute minimum matching bases |
| `outFilterMatchNminOverLread` | 0.66 | Proportional min match: `nMatch >= 0.66 * (Lread-1)` |
| `outFilterMismatchNmax` | 10 | Absolute max mismatches (= `outFilterMismatchNmaxTotal`) |
| `outFilterMismatchNoverLmax` | 0.3 | Proportional max: `nMM/rLength <= 0.3` |
| `outFilterMultimapNmax` | 20 | Max multimapper locations before unmapping |

Note: `Lread-1` is used in proportional thresholds (not `Lread`).

---

## rustar-aligner Equivalent

```rust
fn filter_transcripts(transcripts: &mut Vec<Transcript>, ..., params: &Parameters) {
    let lread = read_len as i32;
    let score_min = (params.out_filter_score_min_over_lread * (lread - 1) as f64) as i32
        .max(params.out_filter_score_min);
    let match_min = (params.out_filter_match_n_min_over_lread * (lread - 1) as f64) as usize
        .max(params.out_filter_match_n_min);
    transcripts.retain(|t|
        t.score >= score_min
        && t.mapped_length() >= match_min
        && t.n_mismatch <= n_mm_max
        && ...
    );
    ...
    if transcripts.len() > params.out_filter_multimap_nmax {
        transcripts.clear(); // TooManyLoci
    }
}
```

---

## Comparison

### `Lread - 1` Usage 🟢

**STAR** uses `(Lread-1)` in both score and match proportional thresholds. Phase 16.12 notes include "Lread-1 filter fix". **Assumed fixed in rustar-aligner.**

### `trBest->nMM / trBest->rLength` — Mismatch Rate 🟡

STAR divides mismatches by `rLength` (the mapped portion of the read, excluding soft-clips). rustar-aligner may divide by `read_len` instead. These differ when there are significant soft clips.

**Investigation needed**: Verify rustar-aligner uses `mapped_length` (= sum of M/I CIGAR ops) not `read_len` in the mismatch rate denominator.

### `nTr` = Number of Transcripts After Score-Range Filter 🟢

STAR's `nTr` is the number of transcripts within `outFilterMultimapScoreRange` of the best. rustar-aligner uses `transcripts.len()` after retaining by score range. **Equivalent.**

### `nW == 0` → unmapType=0 (Other) 🟡

STAR separately tracks "no good windows" as an unmapped reason. rustar-aligner may not distinguish this from general no-alignments case in `UnmappedReason`.

### `outFilterMismatchNmaxTotal` 🟢

STAR computes `outFilterMismatchNmaxTotal = min(outFilterMismatchNmax, floor(outFilterMismatchNoverLmax * Lread))`. This pre-computation creates an effective cap. rustar-aligner should implement this.

---

## `outFilterBySJout` Filter

**STAR** (`ReadAlign_outputAlignments.cpp`):
```cpp
ReadAlign::outFilterBySJout();
if (outFilterBySJoutPass) { ... write SAM ... }
else { statsRA.unmappedOther++; }
```

The `outFilterBySJout` pass filters reads where ALL alignments use novel junctions not present in the SJ filter set (when `--outFilterType BySJout`).

**rustar-aligner**: Implemented as Phase 15 `BySJout` filter. Applied in `lib.rs` as a post-processing pass. Should be equivalent.

---

## Summary

Most filter logic matches STAR. Outstanding items to verify:
1. Mismatch rate denominator: `rLength` (mapped bases) vs `read_len`
2. `outFilterMismatchNmaxTotal` pre-computation
3. `UnmappedReason` enum correctness for "no windows" case
