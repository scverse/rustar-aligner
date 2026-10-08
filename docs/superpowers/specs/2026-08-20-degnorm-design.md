# DegNorm in rustar-aligner: design

Date: 2026-08-20
Status: approved, awaiting implementation plan

## 1. Goal

Reproduce, natively in rustar-aligner, what the
[DegNorm](https://nustatbioinfo.github.io/DegNorm/) pipeline does: correct RNA-seq
read counts for sample- and gene-specific transcript degradation, producing a
Degradation Index (DI) score per gene per sample plus a degradation-adjusted read
count matrix.

The user request was "do during processing what DegNorm does". The half that can
genuinely move into alignment time is coverage computation; the DI math cannot,
because it is defined across samples (see section 2).

## 2. What can and cannot be done online

DegNorm's estimator is a rank-one over-approximation of a gene's coverage matrix
`F` of shape `p x L_i` (p = samples, L_i = exonic length of gene i). With `p = 1`
the rank-one fit is exact and every DI score is 0. A DI score therefore requires
at least two samples and cannot be produced by a single alignment run.

What a single run *can* produce, essentially for free, is the input DegNorm spends
most of its wall clock computing: the per-gene, per-exonic-base coverage vector,
plus the per-gene raw read count. rustar-aligner already assigns each uniquely
mapped read to a gene for `--quantMode GeneCounts`; adding coverage accumulation
on that same path removes the need to write a sorted BAM, index it, and re-read it
with pysam.

Split, therefore:

- **Phase A (alignment time)**: `--quantMode GeneCoverage` emits
  `GeneCoverage.out.bin` per sample.
- **Phase B (merge time)**: `--runMode degNorm` reads N such files and runs the
  NMF-OA pipeline, emitting DI scores and adjusted counts.

Phase B is a pure-CPU pass over coverage matrices; it never touches the genome
index, the SA, or any read file.

## 3. Phase A: coverage capture during alignment

### 3.1 Parameter

`--quantMode GeneCoverage`, composable with the existing values exactly like STAR
composes `GeneCounts TranscriptomeSAM` (`quant_mode_in: Vec<String>`). Requires a
GTF (`--sjdbGTFfile`), same as `GeneCounts`. Enabling it implies building the
`QuantContext`. When absent, zero cost: no allocation, no branch in the hot path
beyond the existing `Option` check.

### 3.2 Data structure

New module `src/quant/coverage.rs`:

```rust
pub struct GeneCoverage {
    /// Prefix sums of merged-exon lengths, len = n_genes + 1 (transcript space).
    offsets: Vec<u64>,
    /// Flat per-base coverage, len = offsets[n_genes].
    cov: Vec<AtomicU32>,
    /// Per-gene raw unique read/fragment count (same rule as GeneCounts col 1).
    counts: Vec<AtomicU32>,
    /// Total reads/fragments counted into any gene (library size).
    n_counted: AtomicU64,
}
```

`offsets` is built from `GeneAnnotation.gene_exons`, which already holds merged,
sorted, absolute-coordinate exon intervals per gene. Transcript-space position of
an absolute coordinate `x` inside gene `g` is
`sum(len(e) for e in exons[..k]) + (x - exons[k].0)` where `k` is the exon
containing `x`, found by `partition_point` (same technique as the existing
`block_is_exonic`).

Memory: human GENCODE, ~20k genes, ~70 Mb merged exonic bases → ~280 MB for `cov`.
Documented in the flag's help text and in the docs page.

### 3.3 Read-to-coverage rule

Reuse the `GeneCounts` unique-hit rule verbatim so counts and coverage always
agree:

- alignment count != 1 → skipped (multimapper; DegNorm's default
  `--non-unique-alignments` off);
- overlapping genes != 1 → skipped (ambiguous);
- otherwise, for each aligned block of the transcript, intersect with the gene's
  merged exons and increment `cov` over the intersection.

Paired-end: union of both mates' blocks, so an overlapping mate pair contributes 1
to each base it covers, not 2. This matches DegNorm's paired-read handling. It is
implemented by collecting both mates' blocks, sorting, merging overlaps, then
incrementing.

Increments are `fetch_add(1, Relaxed)` on the flat array. Contention is negligible
(distinct genes across threads dominate), and correctness does not depend on
ordering.

### 3.4 Exons shared by several genes

DegNorm discards exons that overlap multiple genes when it builds its gene model.
rustar-aligner's `GeneAnnotation` keeps them. The unique-gene rule in 3.3 already
drops any read that hits two genes, which is the same effect at read level.
Recorded as a documented behavioural note, not a divergence in alignment output.

### 3.5 Output file `GeneCoverage.out.bin`

A single gzip stream (`flate2`, already a dependency), containing:

| field | type | note |
|---|---|---|
| magic | 8 bytes | `RSDGNCOV` |
| version | u32 | 1 |
| flags | u32 | bit 0 = paired-end run |
| n_genes | u64 | |
| total_len | u64 | sum of exonic lengths |
| n_counted | u64 | library size, for depth normalisation |
| sample_id | u16 len + bytes | derived from `--outFileNamePrefix`, overridable with `--degNormSampleId` |
| gene table | n_genes x { offset u64, len u32, chr_idx u32, strand u8, count u32 } | |
| gene id block | n_genes x { u16 len, bytes } | `gene_ids[i]` |
| coverage block | total_len x u32 LE | |

All integers little-endian. Reader and writer live in the same module and are
round-trip tested.

## 4. Phase B: `--runMode degNorm`

### 4.1 Parameters

| flag | default | meaning |
|---|---|---|
| `--degNormCoverageFiles` | (required, >= 2) | Phase A files, one per sample |
| `--degNormIter` | 5 | outer DegNorm iterations (DegNorm `--iter`) |
| `--degNormNmfIter` | 100 | inner NMF-OA iterations (`--nmf-iter`) |
| `--degNormDownsampleRate` | 1 | take-every systematic sampling (`--downsample-rate`) |
| `--degNormMinimaxCoverage` | 0 | gene included only if min-over-samples of max coverage >= this (`--minimax-coverage`) |
| `--degNormSkipBaselineSelection` | off | `--skip-baseline-selection` |
| `--degNormBins` | 20 | baseline-selection bins |
| `--degNormMinHighCoverage` | 50 | minimum high-coverage positions for baseline selection |
| `--runRNGseed` | existing | reused for downsampling's random offset |
| `--outFileNamePrefix` | existing | output location |

Names mirror DegNorm's own CLI so the mapping is obvious; the `degNorm` prefix
keeps them clearly non-STAR.

Validation: at least two files; identical gene id vectors and identical exonic
lengths across files (otherwise error naming the first mismatching gene); files
readable and version-compatible.

### 4.2 Algorithm (ported from DegNorm `degnorm/nmf.py`)

Notation follows the source: `x` = raw count matrix (n_genes x p), `F_i` = gene
i's coverage matrix (p x L_i), `rho` = DI matrix (n_genes x p).

**Rank-one approximation.** DegNorm calls `scipy.sparse.linalg.svds(x, k=1)` and
uses `K = u*s` (p x 1), `E = v` (1 x L). Because `p` is tiny (samples), we compute
it by power iteration on the Gram matrix `G = F F^T` (p x p, built in O(p^2 L)):
iterate `u <- G u / ||G u||` to convergence (tol 1e-10, cap 1000 iters), then
`s = ||F^T u||`, `E = (F^T u / s)^T`, `K = u * s`. Sign is fixed so the dominant
component is non-negative. This is numerically the same leading singular triplet
`svds` returns; unit tests assert agreement on planted rank-one matrices.

**NMF-OA (`nmf`).** With `c = 1/sqrt(nmf_iter)` and `lambda = 0`:

```
K, E = rank1(F); est = K E
repeat nmf_iter times:
    res = est - F
    lambda = max(lambda - c*res, 0)
    K, E = rank1(F + lambda)
    est = K E
```

**Initialisation (`ratio_svd` path).** Per gene, `est = max(rank1(F), F)`
elementwise; `rho = 1 - rowsum(F)/(rowsum(est) + 1)`. Then
`low_di = rowmax(rho) < 0.1`; `count_sums = colsum(x[low_di])` (all genes if none
qualify); `norm_factors = count_sums / median(count_sums)`;
`x_weighted = x / norm_factors`; `scale_factors = norm_factors`.

**Baseline selection (per gene, per outer iteration).** Faithful port, including
its quirks:

1. high-coverage columns: `colmax(F) > 0.1 * max(F)`; intersect with the
   systematic downsample when `downsample_rate > 1`.
2. bail out with `rho = 0` if fewer than `min_high_coverage` such columns, or if
   any sample has zero coverage over them.
3. `K, E = nmf(F_bin)`; `rho = 1 - rowsum(F_bin)/(rowsum(K E) + 1)` (the `+1` is
   DegNorm's, kept).
4. bail out if `median(1 - rho) > 1`.
5. run the bin-dropping loop only if
   `n_hi_cov >= max(2, ceil(200/downsample_rate))` and `min(rho) <= 0.2` and
   baseline selection is not skipped.
6. loop while `max(rho) > 0.1`: per column, `res = max_over_samples(((KE - F)/(F+1))^2)`;
   per bin, mean of `res`; drop the argmax bin; shift remaining bins; refit NMF;
   over-approximate `KE = max(KE, F)`; recompute `rho`; stop when bins <=
   `ceil(0.2*bins)` or the surviving length drops below the minimum.
7. on convergence (`max(rho) < 0.2`): clamp tiny `K`, re-derive the envelope on
   the *full* high-coverage matrix as `E = colwise-max(F_start^T / K)`, recompute
   `rho` from it; if `max(rho) > 0.9`, revert to the pre-baseline fit.
8. otherwise revert to the pre-baseline fit with over-approximation.

**Outer loop** (`degnorm_iter` times):

```
F_adj_i = F_i / scale_factors            (row-wise)
rho     = baseline_selection(F_adj)      (clamped to [0, 0.9])
genes that never ran baseline selection get the sample-average DI
          1 - colsum(x_weighted)/colsum(x_adj)
x_adj        = x_weighted / (1 - rho)
norm_factors = colsum(x_adj) / median(colsum(x_adj))
x_weighted   = x_weighted / norm_factors
scale_factors *= norm_factors
```

Parallelised over genes with `rayon` (DegNorm uses joblib threading over gene
chunks); the per-gene fit is independent, so the result is deterministic
regardless of thread count, except for the downsampling draw, which is derived
from `--runRNGseed` and the gene index rather than a shared RNG.

### 4.3 Outputs

Into `<outFileNamePrefix>DegNorm.out/`:

- `DegradationIndex.tab` — genes x samples, tab-separated, header row of sample ids.
- `AdjustedCounts.tab` — `x_adj`, same shape and header.
- `RawCounts.tab` — `x`, for provenance.
- `ScaleFactors.tab` — final per-sample sequencing-depth scale factors.
- `Summary.txt` — parameters, gene counts (total, included, baseline-selected per
  iteration), per-sample median DI.

Estimated coverage matrices (DegNorm's `.pkl` files, used only for plotting) are
out of scope.

## 5. Module layout

```
src/quant/coverage.rs     GeneCoverage accumulator + binary file writer/reader
src/degnorm/mod.rs        run_degnorm(): load, validate, drive, write outputs
src/degnorm/nmf.rs        rank1 / nmf_oa / ratio_svd / baseline_selection (pure math)
src/degnorm/io.rs         output tables
```

`src/degnorm/nmf.rs` takes and returns plain `&[f64]` matrices with explicit
shapes and has no I/O, so it is unit-testable in isolation. `lib.rs` gains one
dispatch arm for `RunMode::DegNorm` and one call site for coverage output next to
the existing `ReadsPerGene.out.tab` write.

## 6. Testing

Unit:

- absolute-coordinate to transcript-coordinate mapping across multi-exon genes,
  including block boundaries and blocks spanning an intron;
- paired-end mate overlap counted once;
- `GeneCoverage.out.bin` round trip (write, read, compare all fields);
- `rank1` recovers a planted rank-one matrix (agreement with a reference singular
  triplet computed by hand for a small case);
- `nmf_oa` produces an over-approximation and converges;
- DI is ~0 for an undegraded synthetic gene, and recovers the expected ordering
  for a planted 3'-biased truncation in one of three samples;
- gene-set mismatch across coverage files is a clean error, not a panic.

Integration (`tests/degnorm.rs`, following `tests/alignment_features.rs`):

- synthetic genome plus GTF, two simulated samples (one with reads truncated
  towards the 3' end of one gene), align both with `--quantMode GeneCounts
  GeneCoverage`, run `--runMode degNorm`, assert the output files exist, the DI of
  the degraded gene in the degraded sample exceeds that of the control, and
  adjusted counts move in the expected direction.

Optional validation script `scripts/compare_degnorm.py`: runs the Python DegNorm
on the same BAMs and correlates DI matrices. Not part of `cargo test`.

## 7. Divergence and documentation

This is a rustar-aligner extension, not STAR behaviour. It gets:

- an entry in `DIVERGENCE.md` under section 4 (implementation divergences with no
  intended alignment-output difference), stating that `GeneCoverage` and
  `degNorm` are additions with no STAR counterpart, are off by default, and do not
  affect alignment output;
- a docs page under `docs/` describing the two-phase workflow;
- `README.md` feature list and `CHANGELOG.md` entries;
- `ROADMAP.md` phase entry.

Alignment output is bit-identical whether or not `GeneCoverage` is enabled; the
integration suite asserts this on the synthetic fixture.

## 8. Out of scope

- Coverage-curve plots and `.pkl` estimated coverage matrices.
- MPI / multi-node execution (DegNorm's `degnorm_mpi`).
- Warm-start directories.
- Reading third-party BAMs as Phase B input; coverage must come from Phase A.
