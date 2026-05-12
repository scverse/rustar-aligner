[← Back to comparison index](README.md)

# stitchAlignToTranscript.cpp vs stitch_align_to_transcript()

**STAR file**: `source/stitchAlignToTranscript.cpp`
**rustar-aligner file**: `src/align/stitch.rs`, function `stitch_align_to_transcript()`
**Purpose**: Given an existing partial transcript (seed A), stitch seed B onto it, scoring the gap between them.

---

## Overview

This function handles all gap types between consecutive seeds:
1. Adjacent (no gap)
2. Equal gap (`rGap == gGap`) — base-by-base scoring in the gap
3. Genomic gap > read gap — deletion or splice junction
4. Read gap > genomic gap — insertion
5. Cross-fragment (different mate in PE combined-read)

---

## Annotated SJ Fast Path

**STAR**:
```cpp
if (sjAB != ((uint)-1) && trA->exons[trA->nExons-1][EX_sjA] == sjAB
    && trA->exons[trA->nExons-1][EX_iFrag] == iFragB
    && rBstart == rAend+1 && gAend+1 < gBstart) {
    // Repeat check for non-canonical:
    if (mapGen.sjdbMotif[sjAB]==0 &&
        (L <= mapGen.sjdbShiftRight[sjAB] || trA->exons[...][EX_L] <= mapGen.sjdbShiftLeft[sjAB]))
        return -1000006;
    // Simple append
    trA->exons[trA->nExons][EX_L] = L;
    ...set canonSJ, shiftSJ, sjAnnot, sjStr...
    trA->nExons++;
    Score += scoreMatch * L + P.pGe.sjdbScore;
}
```

**rustar-aligner**: The annotated SJ fast path is **not implemented**. rustar-aligner always goes through the full jR scanning path even for annotated junctions. The junction is looked up in the DB after scanning, but the fast-path score difference (no gap scoring, direct annotation bonus) may not apply.

**STAR difference**: STAR's fast path skips the jR scan and directly assigns the known motif/strand/shift from the DB. rustar-aligner's general path finds these via scanning, which should produce equivalent results but with extra work.

**Potential issue**: The repeat check `L <= sjdbShiftRight[sjAB]` in the fast path prevents using an annotated junction when the exon is too short to be unambiguous. rustar-aligner may not have this check.

---

## Coordinate Mapping

| STAR | rustar-aligner | Notes |
|------|--------|-------|
| `rAend` (inclusive) | `last_exon.read_end - 1` | STAR is 0-based inclusive, rustar-aligner is exclusive end |
| `gAend` (inclusive) | `last_exon.genome_end - 1` | Same |
| `rBstart` | `eff_read_pos` | Start of seed B (after overlap trimming) |
| `gBstart` | `eff_genome_pos` | Genome start of seed B |
| `L` | `eff_length` | Length of seed B |
| `rGap = rBstart - rAend - 1` | `read_gap = eff_read_pos - last_exon.read_end` | Number of unshared read bases |
| `gGap = gBstart - gAend - 1` | `genome_gap = eff_genome_pos - last_exon.genome_end` | Genome gap |
| `Del = gGap - rGap` | `del = (genome_gap - read_gap) as u32` | Intron/deletion size |
| `shared` (rustar-aligner) | `rGap` (STAR) | Shared gap bases (same value) |
| `gBstart1 = gBstart - rGap - 1` | *(computed inline)* | Last intron base when jR=0 |

---

## Overlap Handling

**STAR**:
```cpp
if (rBend <= rAend) return -1000001;           // Full read overlap
if (gBend <= gAend && iFragA==iFragB) return -1000002;  // Full genome overlap
if (rBstart <= rAend) {                        // Partial overlap: trim B start
    gBstart += rAend - rBstart + 1;
    rBstart = rAend + 1;
    L = rBend - rBstart + 1;
}
```

**rustar-aligner**:
```rust
if last_exon.read_end > eff_read_pos {
    let overlap = last_exon.read_end - eff_read_pos;
    if overlap >= eff_length { return None; }   // Full overlap
    eff_read_pos = last_exon.read_end;
    eff_genome_pos += overlap as u64;
    eff_length -= overlap;
}
if last_exon.genome_end > eff_genome_pos && eff_genome_pos > last_exon.genome_start {
    let g_overlap = (last_exon.genome_end - eff_genome_pos) as usize;
    if g_overlap >= eff_length { return None; }
    ...
}
```

**Assessment**: ✅ Equivalent. rustar-aligner handles both read and genome overlap trimming.

---

## Equal Gap (rGap == gGap > 0): Simple Fill

**STAR**:
```cpp
for (int ii=1; ii<=rGap; ii++) {
    if (G[gAend+ii]<4 && R[rAend+ii]<4) {
        if (R[rAend+ii]==G[gAend+ii]) { Score+=scoreMatch; nMatch++; }
        else { Score-=scoreMatch; nMM++; }
    }
}
```

**rustar-aligner**:
```rust
let shared = read_gap as usize;
gap_mm = count_mismatches_in_region(read_seq, last_exon.read_end, last_exon.genome_end, shared, ...);
d_score += shared as i32 - 2 * gap_mm as i32;
```

**Assessment**: ✅ Equivalent. `count_mismatches_in_region` returns the number of mismatches. Score = matches - mismatches = (shared - mm) - mm = shared - 2*mm. ✓

---

## Large Genomic Gap (gGap > rGap): Junction/Deletion

### Step 1: Move Left

**STAR**:
```cpp
int Score1=0;
int jR1=1;
do {
    jR1--;
    if (R[rAend+jR1] != G[gBstart1+jR1] && G[gBstart1+jR1]<4 && R[rAend+jR1]==G[gAend+jR1])
        Score1 -= scoreMatch;
} while (Score1 + P.scoreStitchSJshift >= 0
      && int(trA->exons[trA->nExons-1][EX_L]) + jR1 > 1);
```

Starting at jR1=1, moves LEFT (jR1--) while: (score penalty ≤ shift bonus) AND (still within exon A by >1 base). Computes score difference between "use acceptor side" vs "use donor side" for each shifted position.

**rustar-aligner** (`find_best_junction_position` in `score.rs`): Implements this scan. Needs verification that the exact termination condition matches, especially the exon length guard `int(EX_L) + jR1 > 1`.

### Step 2: Scan Right to Find Best Junction

**STAR**:
```cpp
int maxScore2 = -999999;
Score1 = 0;
do {
    if (R[rAend+jR1]==G[gAend+jR1] && R[rAend+jR1]!=G[gBstart1+jR1]) Score1 += scoreMatch;
    if (R[rAend+jR1]!=G[gAend+jR1] && R[rAend+jR1]==G[gBstart1+jR1]) Score1 -= scoreMatch;
    int jCan1 = -1; int jPen1 = 0; int Score2 = Score1;
    if (Del >= P.alignIntronMin) {
        // Check 6 canonical motifs at position jR1
        // GTAG=1, CTAC=2, GCAG=3(pen), CTGC=4(pen), ATAC=5(pen), GTAT=6(pen), else 0(non-can)
        Score2 += jPen1;
    }
    if (maxScore2 < Score2) { maxScore2=Score2; jR=jR1; jCan=jCan1; jPen=jPen1; }
    jR1++;
} while (jR1 < (int)rBend - (int)rAend);  // jR1 < rGap + L
```

Key: The scan goes up to `jR1 = rGap + L - 1` (one before `rBend - rAend`). This means the scan extends deep into seed B.

**rustar-aligner** (`find_best_junction_position`): The `eff_length` parameter passed as upper bound. If rustar-aligner's upper bound is `eff_length` (= L), it matches STAR's `rGap + L - 1` in rustar-aligner's local coordinate system (where 0 = start of shared region).

### Step 3: Repeat Length Search

**STAR** searches backward and forward from `gAend + jR` to find repeat length `jjL` and `jjR`. **rustar-aligner** implements this equivalently.

### Step 4: Flush Non-Canonical/Deletions Left

**STAR**:
```cpp
if (jCan <= 0) { // non-canonical or deletion
    jR -= jjL;   // flush left
    jjR += jjL; jjL = 0;
    if (int(trA->exons[...][EX_L]) + jR < 1) return -1000005;
}
```

**rustar-aligner**: Needs verification that non-canonical junctions are flushed left by `jjL`.

### Step 5: Score Donor + Acceptor Bases (The Extended Range)

**STAR**:
```cpp
for (int ii=min(1,jR+1); ii<=max(rGap,jR); ii++) {
    uint g1 = (ii <= jR) ? (gAend+ii) : (gBstart1+ii);
    if (G[g1]<4 && R[rAend+ii]<4) {
        if (R[rAend+ii] == G[g1]) {
            if (ii >= 1 && ii <= rGap) { Score += scoreMatch; nMatch++; }
        } else {
            Score -= scoreMatch; nMM++;
            if (ii < 1 || ii > rGap) { Score -= scoreMatch; nMatch--; }
        }
    }
}
```

This loop covers:
- **Normal range** (`1 <= ii <= rGap`): Score shared bases against correct genome side (donor or acceptor depending on jR).
- **Extended right** (`rGap < ii <= jR` when `jR > rGap`): Score seed B bases shifted to donor side. No positive score/match increment (they were already counted as B-seed matches). Mismatches subtract TWICE (cancel the presumed match + penalty).
- **Extended left** (`jR+1 <= ii < 1` when `jR < 0`): Score exon-A bases shifted to intron. Same cancel-and-penalize logic.

**rustar-aligner** has three separate code blocks for these cases:
1. Donor-side shared bases (`junction_offset > 0`)
2. Acceptor-side shared bases (`acceptor_bases > 0`)
3. Extended left range (`jr_shift < -(shared)`)
4. Extended right range (`jr_shift > 0`) — **Fixed in Phase 16.29**

**Potential remaining issue** (D6 in DIFFERENCES.md): When `0 > jr_shift > -(shared)` (junction shifted left within the shared region), rustar-aligner's `junction_offset = (shared + jr_shift)` and `acceptor_bases = shared - junction_offset = -jr_shift`. This means some shared bases are scored against the acceptor genome side. Let's verify this gives the same result as STAR's unified loop for this case:

For `jr_shift = -1`, `shared = 3`:
- STAR: jR = jr_shift + shared = 2. Loop ii: min(1,3)=1 to max(3,2)=3.
  - ii=1: g1 = gAend+1 (donor). ii <= jR=2, so donor.
  - ii=2: g1 = gAend+2 (donor). ii <= jR=2, so donor.
  - ii=3: g1 = gBstart1+3 (acceptor). ii > jR=2, so acceptor.
  - All within [1,3] range → normal score logic
- rustar-aligner: junction_offset = (3-1) = 2, acceptor_bases = 1.
  - Score 2 bases against donor genome, then 1 base against acceptor.
  - Same as STAR's loop ✓

**Conclusion on D6**: For the in-range left-shift case, rustar-aligner's split donor/acceptor scoring is equivalent to STAR's unified loop. **D6 severity downgraded to 🟢.**

---

## Cross-Fragment (Mate Boundary, canonSJ = -3)

**STAR**:
```cpp
} else if (gBstart + trA->exons[0][EX_R] + P.alignEndsProtrude.nBasesMax >= trA->exons[0][EX_G]
        || trA->exons[0][EX_G] < trA->exons[0][EX_R]) {
    // STAR: extend A's 3' end rightward into the gap, then extend B's 5' leftward
    extendAlign(R, G, rAend+1, gAend+1, 1, 1, DEF_readSeqLengthMax, ...); // A extends right
    extendAlign(R, G, rBstart-1, gBstart-1, -1, -1, extlen, ...);          // B extends left
    trA->canonSJ[trA->nExons-1] = -3;
```

**rustar-aligner** (`stitch_align_to_transcript`, mate boundary block): Just scores the new seed directly (`new_wt.score += wa.length as i32`) without extending into the inter-mate gap. The extensions happen later in `finalize_transcript`.

**Assessment**: rustar-aligner skips the inter-mate extension in `stitch_align_to_transcript`. STAR extends both mates inward toward the gap when processing the mate boundary. This could affect how many bases are "soft-clipped" vs. matched in the inter-mate region for PE reads. The final extensions at `finalize_transcript` time may or may not replicate this.

**Impact**: 🟡 Possible PE soft-clip difference for reads where the mates don't overlap.

---

## Mismatch Acceptance Check

**STAR**:
```cpp
if ((trA->nMM + nMM) <= outFilterMismatchNmaxTotal
  && (jCan < 0 || (jCan < 7 && nMM <= (uint)P.alignSJstitchMismatchNmax[(jCan+1)/2]))) {
```

Two conditions:
1. Total mismatches cumulative over whole read must not exceed `outFilterMismatchNmaxTotal`
2. For junctions (jCan ≥ 0), local gap mismatches must not exceed `alignSJstitchMismatchNmax[motif_category]`

**rustar-aligner** (`scorer.stitch_mismatch_allowed(&motif, gap_mm)`): Implements the per-junction stitch mismatch check. Needs verification that the cumulative total mismatch check is also applied in the right place.

**Assessment**: 🟡 Verify that rustar-aligner also checks the cumulative `trA->nMM + nMM` condition.

---

## Score Update and Exon Adjustment

**STAR** (deletion case):
```cpp
trA->exons[trA->nExons-1][EX_L] += jR;           // Extend exon A by jR
trA->exons[trA->nExons][EX_L] = rBend-rAend-jR;  // New exon B length
trA->exons[trA->nExons][EX_R] = rAend+jR+1;      // New exon B read-start
trA->exons[trA->nExons][EX_G] = gBstart1+jR+1;   // New exon B genome-start
```

**rustar-aligner**:
```rust
last.read_end = (last.read_end as i64 + shared as i64 + jr_shift as i64) as usize;
...
let b_read_start = (eff_read_pos as i64 + jr_shift as i64) as usize;
let b_len = (eff_length as i64 - jr_shift as i64).max(0) as usize;
```

**Assessment**: Converting STAR's `rAend+jR+1` to rustar-aligner notation:
- `rAend+jR+1 = (last_exon.read_end-1) + (jr_shift + shared) + 1 = last_exon.read_end + jr_shift + shared`
- rustar-aligner: `b_read_start = eff_read_pos + jr_shift = (last_exon.read_end + shared) + jr_shift` ✓

Both give the same B start position. ✅
