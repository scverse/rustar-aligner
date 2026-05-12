[← Back to ROADMAP](../ROADMAP.md)

# Phase 15: SAM Tags + Output Correctness

**Status**: Complete (Phases 15.1-15.6 + PE alignment fix)

**Goal**: Add all SAM optional tags required by downstream tools (featureCounts, RSEM, StringTie, GATK, Picard, samtools markdup). Fix paired-end output bugs. Implement `--outSAMattributes` enforcement.

---

## Phase 15.1: NH, HI, AS, NM Tags ✅ (2026-02-10)

**Implementation** (`src/io/sam.rs`):
- `compute_edit_distance()` — sums `n_mismatch + Ins(n) + Del(n)` from CIGAR
- Tags inserted via `record.data_mut().insert(Tag, Value::from(i32))`
- Added to both `transcript_to_record()` and `build_paired_mate_record()`

**STAR Comparison** (8484 position-matching reads):

| Tag | Agreement | Notes |
|-----|-----------|-------|
| NH | 98.3% | 144 differ from seeding differences |
| HI | 100% | Perfect |
| AS | 98.7% | 2 same-CIGAR diffs (minor scoring) |
| NM vs nM | 97.7% | 0 unexplained: diffs = indel bases |

**Files**: `src/io/sam.rs`

---

## Phase 15.2: XS Tag + Secondary Flag + outSAMmultNmax ✅ (2026-02-10)

**Implementation**:
1. **SECONDARY flag** — `FLAGS |= SECONDARY` when `hit_index > 1`
2. **XS tag** — `derive_xs_strand()` from junction motifs via `implied_strand()`. Emits `XS:A:+/-` for spliced reads. Gated by `--outSAMstrandField intronMotif`.
3. **outSAMmultNmax** — Caps output records; NH reflects capped count; MAPQ from true n_alignments.

**STAR Comparison** (8869 common reads):

| Metric | Value |
|--------|-------|
| Primary FLAG agree | 99.8% |
| NH agree | 96.2% |
| XS tags (intronMotif) | 221 (109+, 112-) |
| Secondary FLAG correctness | 100% |

**Files**: `src/io/sam.rs`, `src/params.rs`

---

## Phase 15.3: jM, jI, MD Tags ✅ (2026-02-12)

- `junction_annotated: Vec<bool>` added to Transcript/DpState — propagated from DP
- `build_jm_tag()`: motif codes 0-6 per STAR convention, +20 if annotated → `B:c`
- `build_ji_tag()`: walk CIGAR for RefSkip ops → 1-based per-chr intron pairs → `B:i`
- `build_md_tag()`: walk CIGAR, compare read vs genome → match counts / mismatch / `^DEL` → `Z:`

**Verified**: MD on all 9774 records, jM/jI on 448 spliced records.

**Files**: `src/io/sam.rs`, `src/align/transcript.rs`, `src/align/stitch.rs`, `src/align/read_align.rs`, `src/io/bam.rs`, `Cargo.toml`

---

## Phase 15.4: PE FLAG/PNEXT Fixes ✅ (2026-02-12)

**Problems**:
1. FLAG 0x20 assumed opposite strand — should use mate's actual `is_reverse`
2. PNEXT used own global position — should use mate's per-chr position
3. RNEXT used own chr_idx — should use mate's
4. Tags computed from mate1's transcript for both mates

**Fix**: `PairedAlignment` now stores `mate1_transcript` + `mate2_transcript`. Each record gets mate's actual strand/position/chr and own transcript's tags.

**Files**: `src/align/read_align.rs`, `src/io/sam.rs`, `src/lib.rs`

---

## PE Alignment Fix: Independent Mate Alignment + Pairing ✅

**Problem**: Phase 8 pooled seeds from both mates into unified clusters, then split back. Failed because mates map 200-400bp apart — seeds rarely coexist in same cluster. **0% mapped** on real data.

**Fix**: Align each mate independently via SE `align_read()`, then pair by chr + distance:
1. `align_read()` per mate → transcripts
2. Cross-product pairing: same chr, `check_proper_pair()` distance
3. Dedup by location, sort by combined score, deterministic tie-breaking
4. Score range + multimap filtering

**Results** (10k PE):
- **87.1% mapped** (was 0%): 81.3% unique, 5.9% multi, 12.9% unmapped
- **95.7% per-mate position agree** with STAR
- 72 shared junctions, 100% motif agreement

**Files**: `src/align/read_align.rs`, `src/lib.rs`

---

## Phase 15.5: --outSAMattributes Enforcement ✅ (2026-02-12)

**Implementation**:
- `Parameters::sam_attribute_set()` expands to `HashSet<String>`:
  - `Standard` → {NH, HI, AS, NM, nM}
  - `All` → {NH, HI, AS, NM, nM, MD, jM, jI, XS}
  - `None` → {}
  - Explicit list collected as-is
- Each tag insertion gated by `attrs.contains()`
- XS requires both `--outSAMstrandField intronMotif` AND `XS` in attribute set

**Files**: `src/io/sam.rs`, `src/params.rs`

---

## Phase 15.6: nM Tag (Mismatch Count) ✅ (2026-02-13)

STAR outputs `nM:i:N` (mismatches only). rustar-aligner also outputs `NM:i:N` (edit distance = mismatches + indels). Both emitted by default.

- `Tag::new(b'n', b'M')` with value `transcript.n_mismatch`
- nM = NM for reads without indels; nM < NM for reads with indels

**Files**: `src/io/sam.rs`, `src/params.rs`
