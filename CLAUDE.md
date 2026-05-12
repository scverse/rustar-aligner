# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important: Git Workflow

**DO NOT commit changes automatically.** The user will review, test, and commit changes themselves. Claude should:
- Make code changes as requested
- Suggest what should be committed
- Let the user handle `git add`, `git commit`, and `git push`

## Project Overview

rustar-aligner is a Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), an RNA-seq aligner originally written in C++ by Alexander Dobin. Licensed under MIT to match the original STAR license.

The primary goal is a faithful port — matching the original STAR behavior as closely as possible. Extra features and divergences from the original will come in later releases/forks. When implementing, refer to the [STAR source code](https://github.com/alexdobin/STAR) to ensure correctness and behavioral parity.

## Build Commands

Rust 2024 edition. Standard Cargo commands:

```bash
cargo build            # Debug build
cargo build --release  # Release build
cargo test             # Run all tests
cargo test <name>      # Run a single test by name
cargo clippy           # Lint
cargo fmt              # Format code
```

Always run `cargo clippy`, `cargo fmt --check`, and `cargo test` before considering a phase complete.

## Current Status

**396 tests passing, 0 clippy warnings.** SE: 8613/8926 compare_sam.py (96.5%; note: lower due to seeded-RNG tie-break PR diverging from STAR's mt19937), **99.815% faithfulness (tie-adjusted)** (8611/8627 non-tie reads exact), 299 tie-breaking diffs excluded. 1 CIGAR-only disagree (ERR12389696.13573895, insertion placement, seed-level tie). **0 STAR-only / 0 rustar-aligner-only SE reads**. PE: **8390 both-mapped** (STAR: 8390), **0 half-mapped**, 0 MAPQ inflations / 0 deflations, **99.883% PE exact faithfulness (tie-adjusted)** (16284/16306, 475 tie-breaking diffs excluded), **0 proper-pair diffs**, **0 NH diffs**. Phase 17.A: `scoreSeedBest` pre-extension. Phase 17.B: per-mate seeding. Phase 17.C: STAR-faithful SCORE-GATE + mappedFilter. Phase 17.D: combined-span penalty fix + dedup ordering. Phase 17.8: `--quantMode GeneCounts`. Phase E fix (2026-04-21): mate_id-aware diagonal dedup. Phase E2 (2026-04-22): STAR-faithful combined-read seeding. Phase E3 (2026-04-22): combined-threshold half-mapped fallback. Phase E4 (2026-04-22): PE-CHECK2 unconditional. Phase E5 (2026-04-23): split_combined_wt n_mismatch propagation. Phase E6 (2026-04-24): tie-adjusted faithfulness metric in assess_faithfulness.py. Phase F1: --runRNGseed + seeded primary tie-break (PR #5). Phase F2: --outSAMattrRGline (PR #6). Phase F3: --quantMode TranscriptomeSAM (PR #7). Phase F4: SJDB insertion into Genome+SA at genomeGenerate (PR #8). Phase G1 (2026-04-29): junction_shifts fix in split_combined_wt (rDNA cross-copy false-splice filter). Phase G2 (2026-04-29): MAX_RECURSION 10k→100k + sa_pos_to_forward overflow fix (ERR12389696.7118031 NH=3→9). Phase 17.2 (2026-04-29): coordinate-sorted BAM output (`--outSAMtype BAM SortedByCoordinate` → `Aligned.sortedByCoord.out.bam`). Phase 17.4 (2026-04-29): `--outReadsUnmapped Fastx` → `Unmapped.out.mate1` / `Unmapped.out.mate2`; writes unmapped + TooManyLoci reads; PE writes both mates for fully-unmapped and half-mapped pairs. Phase 17.6 (2026-05-01): `--outStd SAM/BAM_Unsorted/BAM_SortedByCoordinate` — routes primary alignment output to stdout via `Box<dyn AlignmentWriter>` trait dispatch; `SamStdoutWriter`, `BamStdoutWriter`, `SortedBamStdoutWriter` in sam.rs/bam.rs; verified with samtools pipe (967 records). Phase G3 (2026-05-01): SA tie-breaking fix — `compare_suffixes` tie-breaker changed from `pos_b.cmp(&pos_a)` to `packed_a.cmp(&packed_b)` (ascending by packed SA value with strand bit); rustar-aligner SA is now **byte-for-byte identical** to STAR's SA for the yeast genome (10,862 → 0 entry diffs). diff AS: 6→4 cases (4 remaining are rustar-aligner improvements: .844151 VIII 0mm vs STAR VII 6mm, .4972950 spliced vs unspliced mate2). Phase 17.3 (2026-05-01): PE chimeric detection — `detect_inter_mate_chimeric` in `chimeric/detect.rs`; intra-mate multi-cluster chimeric via cluster splitting + mate2 read_pos adjustment; inter-mate chimeric for discordant pairs (diff chr, same strand, or >1Mb); `align_paired_read` returns 4-tuple including `Vec<ChimericAlignment>`; no benchmark regression (8390 both-mapped, 0 half-mapped). Phase 17.11 (2026-05-01): `--chimOutType WithinBAM` — chimeric alignments written as supplementary records (FLAG 0x800) in primary BAM; donor record has full SEQ + SA tag; acceptor has FLAG 0x800 + SA tag + empty SEQ; `build_within_bam_records` in `chimeric/output.rs`; `chim_out_junctions()` / `chim_out_within_bam()` helpers in params.rs; supports mixed `--chimOutType Junctions WithinBAM`. Phase 17.7 (2026-05-01): GTF tag parameters — `--sjdbGTFchrPrefix`, `--sjdbGTFfeatureExon`, `--sjdbGTFtagExonParentTranscript`, `--sjdbGTFtagExonParentGene`; `_configured` variants in `junction/gtf.rs`, `quant/mod.rs`, `quant/transcriptome.rs`, `junction/mod.rs`; all 4 production paths thread params; backward-compat wrappers preserve zero test disruption. Phase 17.9 (2026-05-01): `--outBAMcompression` (BGZF level -1–9, default 1; -1/0=NONE, 1-8=flate2 levels, ≥9=BEST) + `--limitBAMsortRAM` (bytes, 0=unlimited; aborts sort if ~400 bytes/record estimate exceeds limit); `bgzf_compression()` + `make_bgzf_writer()` helpers in `io/bam.rs`; threaded through all 4 BAM writers (unsorted file, sorted file, unsorted stdout, sorted stdout). PE chimericDetectionOld (2026-05-01): per-mate `detect_chimeric_old` called on `all_m1_transcripts` / `all_m2_transcripts` pools after `filter_paired_transcripts` in `read_align.rs`. Phase 17.12 (2026-05-01): BySJout disk buffering — `BySJReadMeta` struct + `NamedTempFile` SAM temp file replaces `Vec<AlignmentBatchResults>`; `create_bysj_writer` / `bysj_write_records` / `bysj_read_n_records` helpers in `io/sam.rs`; `tempfile` moved to `[dependencies]`. Phase 17.13 (2026-05-01): 8 integration tests in `tests/alignment_features.rs` — synthetic 20kb genome with planted GT-AG intron; tests cover BAM output, PE alignment, spliced reads, BySJout, GeneCounts, unmapped output, two-pass mode. Phase 12.2 (2026-05-04): SE chimeric Tier 1b soft-clip re-mapping — `detect_from_soft_clips` in `chimeric/detect.rs` re-seeds the primary alignment's soft-clipped bases when `detect_chimeric_old` finds no partner; `adjust_read_positions` helper shifts sub-seq coords into full-read space for right clips; called as Step 3c in `read_align.rs`. Phase 17.10 (2026-05-04): Chimeric Tier 3 — `detect_from_chimeric_residuals` in `chimeric/detect.rs` re-seeds outer uncovered read regions (before donor / after acceptor) of each found chimeric pair; enables 3-way gene-fusion detection; called as Step 3d in `read_align.rs`. See [ROADMAP.md](ROADMAP.md) for detailed phase tracking and [docs-old/](docs-old/) for per-phase development notes. The published Astro Starlight docs site is in [docs/](docs/).

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  params.rs        -- ~52 STAR CLI params via clap derive, --camelCase long names
  error.rs         -- Error enum with thiserror (Parameter, Io, Fasta, Index, Alignment, Gtf)
  mapq.rs          -- MAPQ calculation (STAR lookup table: n=1→255, n=2→3, n≥5→0)
  stats.rs         -- Alignment statistics, Log.final.out writer, UnmappedReason enum
  genome/
    mod.rs         -- Genome struct, padding logic, reverse complement, file writing
    fasta.rs       -- FASTA parser, base encoding (A=0,C=1,G=2,T=3,N=4)
  index/
    mod.rs         -- GenomeIndex (build + load + write)
    packed_array.rs-- Variable-width bit packing (1-64 bits per element)
    suffix_array.rs-- SA construction, custom comparator, strand encoding
    sa_index.rs    -- K-mer lookup table (35-bit entries with flags)
    io.rs          -- Load index from disk (Genome, SA, SAindex)
  align/
    mod.rs         -- Module definition
    seed.rs        -- Seed finding via hierarchical SAindex lookup + MMP search
    stitch.rs      -- Seed clustering + DP stitching + alignment extension (extendAlign)
    score.rs       -- Scoring functions (gaps, mismatches, splice junctions)
    transcript.rs  -- Transcript struct (exon coords, CIGAR, scores)
    read_align.rs  -- Per-read alignment driver
  io/
    mod.rs         -- Module exports
    fastq.rs       -- FASTQ reader (plain + gzip, noodles wrapper)
    sam.rs         -- SAM writer (header + records, noodles wrapper)
    bam.rs         -- BAM writer (BGZF compression, streaming unsorted + in-memory sorted output)
  junction/
    mod.rs         -- GTF parsing, junction database, motif detection, two-pass filtering
    sj_output.rs   -- SJ.out.tab writer
    gtf.rs         -- GTF parser (internal)
  quant/
    mod.rs         -- Gene-level read counting (--quantMode GeneCounts, ReadsPerGene.out.tab)
  chimeric/
    mod.rs         -- Module exports
    detect.rs      -- Chimeric detection (Tier 1: transcript-pair, Tier 1b: soft-clip re-seed, Tier 2: multi-cluster, Tier 3: residual re-seed)
    segment.rs     -- ChimericSegment and ChimericAlignment data structures
    score.rs       -- Junction type classification, repeat length calculation
    output.rs      -- Chimeric.out.junction writer (14-column format)
```

## Development Philosophy — Match STAR Exactly

**Never backtrack on changes that make the method closer to STAR's logic, even if metrics temporarily regress.** A regression in alignment statistics (position agreement, splice rate, etc.) after making the code more STAR-faithful does NOT mean the change was wrong — it means the implementation is incomplete. More STAR-matching changes are needed to close the gap.

- It is better to go backwards in metrics temporarily to eventually reach 100% matching results, than to keep approximate methods that can never get there.
- Every algorithm, threshold, and code path should match STAR's C++ source exactly. Any deviation — even a small one — should be corrected.
- Do not revert STAR-faithful changes. Instead, investigate what STAR does differently in the surrounding logic and implement that too.
- When a STAR-matching change causes regression, the response is: "What else does STAR do that we don't?" — not "revert and try something different."

## Key Conventions

- **Crate/binary name is `rustar-aligner`** — lib name in Rust is `rustar_aligner`; `#![allow(non_snake_case)]` in lib.rs is retained for STAR-style camelCase identifiers
- **STAR params use `--camelCase` naming** — clap `#[arg(long = "camelCase")]` maps to snake_case Rust fields
- **Multi-value params** (genomeFastaFiles, readFilesIn, outSAMtype, outSAMattributes, chimOutType, alignSJstitchMismatchNmax, outSJfilterIntronMaxVsReadN) need explicit `num_args`
- **Negative defaults** (scoreGapNoncan=-8, readMapNumber=-1, etc.) need `allow_hyphen_values = true`
- **`outSAMtype`** is parsed as raw `Vec<String>` then structured via `Parameters::out_sam_type()` method
- **Validation** beyond clap's type checking is in `Parameters::validate()` (e.g. genomeGenerate requires FASTA files)
- **No async** — CPU-bound work; async adds complexity with zero benefit
- **Error handling** — `thiserror` for `Error` enum, `anyhow` for top-level result propagation

## Dependencies

```toml
[dependencies]
clap = { version = "4", features = ["derive"] }
anyhow = "1"
thiserror = "2"
log = "0.4"
env_logger = "0.11"
memmap2 = "0.9"
byteorder = "1"
noodles = { version = "0.80", features = ["fastq", "sam", "bam", "bgzf"] }
bstr = "1"
flate2 = "1"
rayon = "1"
dashmap = "6"
chrono = "0.4"

[dev-dependencies]
tempfile = "3"
assert_cmd = "2"
predicates = "3"
```

## Testing Pattern

- Unit tests: `#[cfg(test)]` in each module
- Integration tests: `tests/` directory (future phases)
- Every phase uses differential testing against STAR where applicable
- Test data tiers: synthetic micro-genome → chr22 → full human genome

**Current test status**: 364/364 tests passing (359 unit + 5 integration), 0 clippy warnings

## Known Issues — Disagreement Root Causes (10k SE yeast)

**299 total position disagreements — ALL verified as genuine ties** (SA-order ties + seeded-RNG tie-break divergence from STAR's mt19937):

Both tools find identical alignment sets for all 299 disagreements. Primary selection differs either due to SA iteration order or RNG seed divergence (PR #5: `--runRNGseed`, uses `StdRng`, not `mt19937`).

- **100+ diff-chr ties** — same set of alignments, different repeat copy chosen as primary.
- **27+ same-chr ties** — same alignment set, different primary due to tie-breaking.

**1 CIGAR-only disagreement (same position, different CIGAR):**
- `ERR12389696.13573895`: both tools align to XV:218357 MAPQ=255, but rustar-aligner gives `100M1I45M4S` (insertion at read pos 100) while STAR gives `108M1I37M4S` (insertion at 108). Root cause: both alignments score AS=133. The 71-base seed is found at RC pos 29 (rustar-aligner) vs RC pos 37 (STAR) due to different Lmapped chain paths through a long homopolymer region. Same diagonal, different starting position → different insertion placement. Seed-level tie.

**1 truly actionable SE issue:**

1. **CIGAR-only insertion placement (1 read)** — `ERR12389696.13573895` (see above). Requires matching STAR's exact Lmapped chain to fix.

Previously listed issues now resolved:
- `ERR12389696.18296181` (rustar-aligner-only false splice): filtered out by score threshold (score=95/92 < 98 threshold for 150bp read). Both the 1-exon soft-clip AND the 2-exon false splice fail the minimum score gate — read correctly unmapped.
- `ERR12389696.13766843` (STAR-only NM=10 read): rustar-aligner already maps this correctly (VII:24449, 28S121M1S, MAPQ=255) — it was incorrectly listed as "STAR-only".

**Phase 16.38 (STAR-faithful filter ordering):**
- Moved dedup + score-range filter (multMapSelect) to run BEFORE quality filters (mappedFilter). STAR's ordering: `multMapSelect → mappedFilter`. Old code ran quality filters first, which could remove the high-scoring primary leaving a lower-scoring secondary as apparent "best". No benchmark change (the specific false splice read was already filtered by score threshold), but structurally correct for edge cases.

**Phase 16.37 (alignIntronMax=0 fix) resolved:**
- `ERR12389696.5825571`: now aligns as `XV:80779 121M607028N13M16S` (exact match with STAR). Root cause: rustar-aligner computed a finite intron limit of 589824 when `alignIntronMax=0`, blocking the 607kb intron. STAR uses `>0` guard in `stitchAlignToTranscript.cpp` line 100 — alignIntronMax=0 means no limit. Fix: use `u32::MAX` sentinel in score.rs.
- `ERR12389696.16030539`: MAPQ inflation resolved. Both tools now find two alignments (XV:121224 128M925400N10M12S and XV:598336 128M448288N10M12S), both MAPQ=3. Primary selection differs (tie). MAPQ inflate: 1→0.

See [ROADMAP.md](ROADMAP.md) and [docs-old/](docs-old/) for full issue tracking.

## PE Status (Updated 2026-04-29 — Phase G2)

**Current**: **PE both-mapped = 8390** (STAR: 8390, exact match!), **half-mapped = 0**, **99.883% PE exact faithfulness (tie-adjusted)** (16284/16306, 475 tie-breaking diffs excluded). MAPQ inflations: **0**, deflations: 0. NH diffs: **0**. 1 STAR-only mate (`.18919121`, SA-level diff), 1 rustar-aligner-only mate (`.6302610`, pre-existing FP).

**Phase G2** (2026-04-29): `MAX_RECURSION` 10,000→100,000 + `sa_pos_to_forward` overflow fix. `ERR12389696.7118031` was the sole source of both NH diffs and MAPQ inflations (NH=3 vs STAR's NH=9, rustar-aligner MAPQ=1 vs STAR MAPQ=0). Root cause: the 47-WA rDNA cluster (4 copies × multiple seeds per mate) exhausted 10k recursions before exploring the 4th within-copy pair. Fix: raise the per-cluster recursion budget from 10k to 100k. Also fixed: `sa_pos_to_forward` underflow panic for reverse-strand seeds near genome boundary (now `saturating_sub`). Also added guard in `finalize_transcript` to reject WTs where `adjusted_genome_start + ref_len > n_genome`.

**Phase G1** (2026-04-29): `split_combined_wt` junction_idx fix. Reduced `.16980960`'s pairs from 11 to 9, matching STAR's NH=9.

**Note on faithfulness change**: Phase F1 (--runRNGseed PR) changed PE tie-breaking from SA-order to seeded StdRng, increasing tie-breaking diffs. Phase G1 improved faithfulness from 99.755% → 99.865%. Phase G2: 99.865% → 99.883%.

## Remaining Limitations (Top 2)

- No STARsolo single-cell features — Phase 14 (deferred)
- 4 PE AS diffs (rustar-aligner improvements, not bugs): `.844151` finds VIII:451791 0mm vs STAR's VII:1001391 6mm; `.4972950` finds correct spliced mate2 vs STAR's unspliced. Both cases: STAR's combined-window approach fails to stitch a PE pair at the better location.

See [docs-old/phase17_features.md](docs-old/phase17_features.md) for full feature status.
