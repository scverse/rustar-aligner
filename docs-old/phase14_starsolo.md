[← Back to ROADMAP](../ROADMAP.md)

# Phase 14: STARsolo (Single-Cell)

**Status**: In progress — **MVP complete (14.1–14.4)**

**Goal**: A faithful port of STARsolo — turn the aligner into a single-cell RNA-seq
quantifier that matches STAR's `--soloType` output (count matrices, barcode/UMI
correction, cell calling, SAM tags) as closely as the bulk aligner already
matches STAR.

**Prerequisite (met)**: position agreement >99% — SE 99.815% (tie-adjusted),
PE 99.883%. Phase unblocked 2026-06-10.

---

## Architecture

STARsolo is a **layer around** the existing aligner, not a change to it. The core
alignment is untouched:

```
 readFilesIn[0] = cDNA read  ──► existing SE alignment ──► Transcript(s)
 readFilesIn[1] = barcode read (R1: CB+UMI) ──► parse ──► correct vs whitelist
                                                              │
              Transcript + corrected CB + UMI ──► gene assignment (overlapping_genes)
                                                              │
                              collate per (CB, gene) ──► UMI dedup ──► count
                                                              │
                                            Solo.out/<Feature>/raw/matrix.mtx
```

Key reuse points already in the codebase:
- `Transcript` (`src/align/transcript.rs`) carries `chr_idx`, `genome_start/end`,
  `is_reverse`, `exons` — everything gene assignment needs.
- `GeneAnnotation::overlapping_genes()` (`src/quant/mod.rs`) maps an alignment to
  gene indices and is directly reusable for per-cell counting.
- The SE parallel batch loop (`align_reads_single_end` in `src/lib.rs`) is where
  per-read barcode info threads through to a per-cell accumulator.

**Read-file convention** (matches STAR): `--readFilesIn cDNA_read barcode_read`.
The cDNA read is file 0, the barcode read is file 1. A solo run therefore supplies
two files but is a *single-end alignment* run.

---

## Sub-phase plan

| Sub-phase | Description | Status |
|-----------|-------------|--------|
| 14.1 | `--solo*` params + barcode-read input plumbing | ✅ Complete |
| 14.2 | Whitelist load + CB correction (`--soloCBmatchWLtype`) + UMI checks | ✅ Complete |
| 14.3 | Per-read gene assignment + CB/UMI threaded into the alignment loop | ✅ Complete |
| 14.4 | UMI dedup + raw `matrix.mtx` (**MVP complete**) | ✅ Complete |
| 14.5 | `Summary.csv` / `Barcodes.stats` / `Features.stats` | ⬜ Planned |
| 14.6 | Cell filtering (`filtered/` matrix) | ⬜ Planned |
| 14.7 | `CB`/`UB`/`GX`/`GN` SAM tags + `CB_samTagOut` | ⬜ Planned |
| 14.8 | More features: GeneFull, SJ, Velocyto | ⬜ Planned |
| 14.9 | Multi-gene resolution (`--soloMultiMappers`) | ⬜ Planned |
| 14.10 | Other chemistries: CB_UMI_Complex, SmartSeq | ⬜ Planned |
| 14.11 | Differential test harness vs STARsolo + integration tests | ⬜ Planned |

**MVP = 14.1–14.5**: a working 10x Chromium `Gene` count matrix.

### Faithfulness risk notes
- **Read ordering**: cDNA read is FIRST in `--readFilesIn`, barcode read second.
- **CB correction** posterior math and the **`1MM_Directional`** UMI-graph collapse
  are the two algorithms where byte-parity with STAR is fiddly — budget extra
  differential-testing time there (14.2, 14.4).
- **Matrix conventions**: MatrixMarket coordinate format, features × barcodes,
  1-based indices — must match Cell Ranger / STARsolo layout exactly.

---

## Phase 14.1: Params + barcode-read plumbing ✅ (2026-06-10)

**Goal**: Accept `--soloType` and the barcode geometry on the CLI, read the barcode
read alongside the cDNA read, and extract CB+UMI — without yet counting.

**Implementation**:

1. **`src/params/mod.rs`** — `SoloType` enum (`None`, `CbUmiSimple` [alias
   `Droplet`], `CbUmiComplex`, `CbSamTagOut`, `SmartSeq`) with `FromStr`/`Display`.
   12 new parameters:
   - `--soloType`, `--soloCBwhitelist`, `--soloCBstart` (1), `--soloCBlen` (16),
     `--soloUMIstart` (17), `--soloUMIlen` (10), `--soloFeatures` (`Gene`),
     `--soloUMIdedup` (`1MM_All`), `--soloCBmatchWLtype` (`1MM_multi`),
     `--soloCellFilter`, `--soloOutFileNames`, `--soloStrand` (`Forward`).
   - Helpers: `solo_enabled()`, `cdna_read_file()`, `barcode_read_file()`,
     `solo_cb_whitelist_none()`.
   - Validation: solo needs exactly 2 read files; `Gene`/`GeneFull` need a GTF;
     CB/UMI length > 0 for `CB_UMI_Simple`.

2. **`src/solo/mod.rs`** (new) —
   - `SoloBarcodeLayout` — fixed-position geometry, 1-based starts converted to
     0-based; `from_params`, `min_read_len`, `extract`.
   - `CellBarcode` — encoded CB/UMI seq + raw Phred qualities; `cb_has_n`,
     `umi_has_n`, `cb_string`, `umi_string`.
   - `SoloReadReader` / `SoloRead` — lockstep reader over the cDNA and barcode
     FASTQ files; `read_batch`; errors on length mismatch. `open_reader(params)`
     factory.

3. **`src/lib.rs`** — `mod solo;`; `run_single_pass` + `run_pass1` compute
   `n_align_files = if solo { 1 } else { read_files_in.len() }` so a 2-file solo
   run routes to the SE cDNA path; `is_paired` excludes solo.

**Boundary**: 14.1 makes a solo run *parse and validate* and aligns the cDNA read
(producing `Aligned.out.sam`). Barcodes are extracted by `SoloReadReader` but not
yet threaded into the parallel alignment loop or counted — that begins in 14.2,
where per-read barcode handling pairs naturally with whitelist correction.

**Tests**: 6 new in `src/solo/mod.rs` (layout conversion, v2 extraction, too-short
read, N-detection, reader pairing, length-mismatch error) + CLI validation smoke
tests. 447 lib tests, 0 clippy warnings.

**Files**: `src/params/mod.rs`, `src/solo/mod.rs` (new), `src/lib.rs`

---

## Phase 14.2: Whitelist load + CB correction ✅ (2026-06-11)

**Goal**: Load the cell-barcode whitelist and match each read's CB to it exactly
as STAR's read stage does, plus validate the UMI.

**Reference**: STAR `source/SoloReadBarcode_getCBandUMI.cpp` (read stage). The
multi-match *posterior* resolution lives in the collation stage, not here — see
the boundary note below.

**Implementation** (`src/solo/whitelist.rs`, new):

- **Packing** — `pack_barcode` 2-bit packs an encoded barcode into a `u64` with
  `seq[0]` in the high bits (matching `convertNuclStrToInt64`). N-handling:
  `NoN(u64)` / `OneN{packed,pos}` / `ManyN`. `unpack_barcode` reverses it.
- **`CbMatchType`** — decodes `--soloCBmatchWLtype` into STAR's `mm1` /
  `mm1_multi` / `mm1_multi_nbase` / `pseudocounts` flags (Exact, 1MM, 1MM_multi
  [default], `_pseudocounts`, `_Nbase_pseudocounts`).
- **`CbWhitelist`** — `List` (sorted unique packed `Vec<u64>` + original-order
  index for `barcodes.tsv` + per-index `exact_counts` atomics) or `NoWhitelist`.
  `load()` reads plain or gzip, validates equal lengths, rejects N-containing
  whitelist entries.
- **`match_cb`** follows STAR exactly: exact binary search (→ `Exact`, bumps the
  exact-count prior); else single-N substitution (all 4 bases at the N position)
  or 1MM enumeration (every position × 3 alternate bases). One candidate →
  `Corrected`; >1 → `Multi(candidates)` when the multi flag is set (records WL
  index + mismatch position + quality for later resolution) else
  `MultMatchRejected`. Rejections map to STAR's cbMatch codes (`NoMatch` -1,
  `NinCb` -2, `MultMatchRejected` -3).
- **`check_umi`** — any N → `NinUmi` (-23); exact homopolymer → `Homopolymer`
  (-24); else `Ok(packed)`.
- **`CbMatchStats`** — atomic counters for STAR's cbMatch categories.

**Params** (`src/params/mod.rs`): `--soloCBmatchWLtype` validity check;
`solo_cb_match_type()` and `solo_cb_whitelist_path()` helpers; rules that
`--soloCBwhitelist None` requires `Exact`, and `--soloCBlen ≤ 32`.

**Boundary**: the count + quality **posterior** that resolves `CbMatch::Multi`
into one corrected barcode needs the *global* `exact_counts` table, which is only
complete after all reads are processed — so it is a collation-stage operation
deferred to Phase 14.4. Phase 14.2 records the candidates (exactly as STAR's
`cbMatchString`) and accumulates the prior. The matcher is also not yet wired
into the alignment loop; that happens in 14.3 alongside gene assignment.

**Tests**: 13 new in `src/solo/whitelist.rs` (pack roundtrip, N-detection, exact
match + count, 1MM correction, ambiguous multi vs reject, no-match, single-N
correction, many-N reject, Exact-only mode, UMI checks, length-mismatch error,
gzip load, match-type parsing) + CLI validation smoke tests. 460 lib tests,
0 clippy warnings.

**Files**: `src/solo/whitelist.rs` (new), `src/solo/mod.rs`, `src/params/mod.rs`

---

## Phase 14.3: Gene assignment + barcode threading ✅ (2026-06-11)

**Goal**: Assign each cDNA alignment to a gene and wire CB/UMI through the
alignment loop so per-cell (CB, UMI, gene) records are collected.

**Gene assignment** (`src/solo/gene.rs`, new):
- `SoloStrand` (`--soloStrand`: Forward [default] / Reverse / Unstranded).
- `assign_gene_se(transcripts, gene_ann, strand)` — the read's gene set is the
  UNION of strand-filtered `GeneAnnotation::overlapping_genes` across ALL its
  alignments. Exactly one gene → `Gene(idx)`; zero → `NoFeature`; >1 →
  `Ambiguous`; no transcripts → `Unmapped`. A multi-locus read whose loci all
  fall in one gene is therefore still gene-unique (matching STARsolo's default
  `--soloMultiMappers Unique`, unlike `quantMode GeneCounts` which drops every
  multimapper).

**Context + recorder** (`src/solo/mod.rs`):
- `SoloContext` — `build(params, genome)` loads the whitelist and builds the
  gene model from `--sjdbGTFfile`; bundles layout + whitelist + match type +
  strand + `CbMatchStats` + `SoloRecorder`, shared as an `Arc` across threads.
- `SoloRecorder` — thread-safe sink for `SoloCountRecord{cb, umi, gene}` plus
  deferred `SoloMultiRecord` (unresolved 1MM_multi CBs, resolved in 14.4).
- `SoloContext::process_read` — CB match → UMI check → gene assign, recording
  stats and producing a record only when all three succeed.

**Loop** (`src/lib.rs`): new `align_reads_solo` reads cDNA (file 0) + barcode
(file 1) in lockstep via `SoloReadReader`, aligns the cDNA exactly like the SE
path (`align_read` → `build_alignment_records`), writes SAM/BAM, runs
`process_read` per read, and appends records to the recorder in the sequential
write phase. `run_single_pass` dispatches solo runs here; `run_single_pass` /
`run_two_pass` thread `solo_ctx`. A run-end summary logs the barcode-match stats
and record count.

**Boundary / limitations**: the solo loop is single-pass and does not yet emit
BySJout / chimeric / transcriptome-SAM side outputs (not part of the MVP). The
count matrix (`raw/matrix.mtx` + `barcodes.tsv` + `features.tsv`) and 1MM_multi
posterior resolution are Phase 14.4. `--soloStrand` validated in params.

**Tests**: 7 new gene-assignment unit tests + end-to-end
`test_starsolo_gene_assignment` (synthetic genome + GTF + whitelist: 16 cDNA
reads → 16 exact CB matches → 16 resolved (CB,UMI,gene) records). 467 lib + 10
integration tests, 0 clippy warnings.

**Files**: `src/solo/gene.rs` (new), `src/solo/mod.rs`, `src/params/mod.rs`,
`src/lib.rs`, `tests/alignment_features.rs`

---

## Phase 14.4: UMI dedup + raw matrix — MVP COMPLETE ✅ (2026-06-11)

**Goal**: Collapse UMIs and write the raw per-cell count matrix — the first
usable single-cell output.

**Reference**: STAR `SoloFeature_collapseUMIall.cpp` (dedup),
`SoloReadFeature_inputRecords.cpp` (CB multi-resolution),
`SoloFeature_outputResults.cpp` (matrix format).

**Implementation** (`src/solo/count.rs`, new):

- **`UmiDedup`** (`--soloUMIdedup`): `Exact` (distinct UMIs), `NoDedup` (reads),
  `1MM_All` (default — connected components where any two UMIs within Hamming-1
  merge transitively, via union-find), `1MM_Directional` / `_UMItools`
  (`count_hub ≥ 2·count_leaf + dirCountAdd`, `dirCountAdd` 0 / −1).
- **Deferred 1MM_multi CB resolution** — `resolve_multi_cb` picks the candidate
  maximizing STAR's posterior weight `exactCount[cand] · 10^(−q/10)` (prior =
  `whitelist.exact_count_snapshot()`, `q` = mismatch-position Phred); rejects
  when no candidate has positive weight.
- **`build_matrix`** groups reads by `(cell, gene)` into UMI→multiplicity maps
  (resolved multi-CB records folded in), then dedups each.
- **`write_gene_matrix`** writes `Solo.out/Gene/raw/`:
  - `matrix.mtx` — `%%MatrixMarket matrix coordinate integer general`; dims
    `nFeatures nBarcodes nEntries`; entries `gene+1 cell+1 count` (1-based),
    iterated in cell-column order.
  - `features.tsv` — `gene_id <TAB> gene_id <TAB> Gene Expression` (CellRanger
    v3; no gene names available so id is repeated).
  - `barcodes.tsv` — full whitelist in sorted order (matrix column order).

Wired into `align_reads` after alignment. `--soloUMIdedup` validated in params.

**Known approximations to revisit** (differential testing, 14.11): the
`1MM_Directional` absorption is a greedy hub model (faithful default path is
`1MM_All`, which is exact); the CB-posterior acceptance uses no `cbMinP`
threshold (always takes the argmax); `barcodes.tsv` uses sorted (not 10x-file)
order; `--soloCBwhitelist None` matrix output is not yet supported.

**Tests**: 8 new unit tests in `count.rs` (each dedup method incl. transitive
chains and the directional thresholds; multi-CB posterior) + end-to-end
`test_starsolo_gene_matrix` (8 reads, one cell, two Hamming-distant UMI clouds →
2 deduped molecules → matrix `1 1 2`, validated `features.tsv` / `barcodes.tsv`).
475 lib + 10 integration tests, 0 clippy warnings.

**Files**: `src/solo/count.rs` (new), `src/solo/mod.rs`, `src/params/mod.rs`,
`src/lib.rs`, `tests/alignment_features.rs`

---

## Phase 14.CR: CellRanger 4.x/5.x matching — VERIFIED vs real STARsolo ✅ (2026-06-12)

**Goal**: Support the [STARsolo CellRanger-matching flag set](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#matching-cellranger-4xx-and-5xx-results)
and prove the output matches real STARsolo.

**Flags** (`--clipAdapterType CellRanger4 --outFilterScoreMin 30
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR`), implemented from STAR source:

- **`1MM_CR`** (`src/solo/count.rs::cellranger_1mm`) — port of STAR
  `umiArrayCorrect_CR`: UMIs sorted ascending by `(count, umi)`, each corrected
  to its highest-count 1MM neighbor, **non-transitive** (points to the neighbor's
  raw UMI), count = distinct corrected UMIs.
- **`MultiGeneUMI_CR`** (`filter_multi_gene_umi`) — keep the top-read-count gene
  of a multi-gene UMI. `build_matrix` restructured to per-cell
  `umi → gene → read_count` so filtering precedes dedup.
- **`1MM_multi_Nbase_pseudocounts`** — +1 pseudocount on the CB posterior prior
  (`resolve_multi_cb`).
- **`CellRanger4` clip** (`src/solo/mod.rs::clip_adapter_cr4`) — TSO 5' clip +
  polyA 3' trim, conservative (no-op on adapter-free reads), applied in
  `align_reads_solo` before fixed Nbases clipping.

All four validated in `params.rs`.

**Differential test** (`test/solo_cellranger_diff.py`): generates a synthetic 10x
dataset (two 2-exon genes, whitelist, cDNA + barcode reads with a planted 1MM
UMI pair), runs the full CellRanger flag set on BOTH rustar-aligner and real
STAR, and compares the decoded `{(barcode, gene_id): count}` matrices.

**Result — byte-identical match, 3/3 deterministic:**
```
(AAAACCCCGGGGTTTT, GENEA) = 2   # 1MM_CR collapsed M(x5)+M-1mm(x1) -> 1, +N(x3) -> 2
(AAAACCCCGGGGTTTT, GENEB) = 1
(ACACACACGTGTGTGT, GENEA) = 1
```

**Why a container**: STAR 2.7.11b reads 0 input reads on Apple-Silicon macOS (a
known STAR/macOS bug — `nextChar=-1` immediate EOF — present in both the homebrew
bottle and a from-source build). The reference therefore runs in a Linux
container (`test/Dockerfile.solodiff` — Debian + `rna-star` 2.7.10b + Rust),
driven by `test/solo_diff_docker.sh` via colima (no Docker Desktop needed). On a
host with a working STAR, `python3 test/solo_cellranger_diff.py` runs it directly.

A committed cargo test (`test_starsolo_cellranger_style_matrix`) asserts the same
CellRanger-style matrix (including the 1MM_CR collapse) without needing STAR, and
each CellRanger algorithm has unit tests in `src/solo/count.rs`.

---

## MVP status

Phases 14.1–14.4 deliver a working **10x Chromium `Gene`** quantifier:
`--soloType CB_UMI_Simple --soloCBwhitelist <wl> --soloFeatures Gene
--sjdbGTFfile <gtf> --readFilesIn cDNA.fq barcode.fq` aligns the cDNA reads and
writes `Solo.out/Gene/raw/{matrix.mtx, barcodes.tsv, features.tsv}`. Remaining
phases (14.5–14.11) add stats files, cell filtering, SAM tags, more features,
multi-gene resolution, other chemistries, and the differential-test harness.
