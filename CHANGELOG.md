# rustar-aligner Changelog

<!--
Release notes are extracted from this file by the release workflow.
Each released version needs a heading of the form:

    ## [Version X.Y.Z](https://github.com/scverse/rustar-aligner/releases/tag/vX.Y.Z) - YYYY-MM-DD

Sections commonly used: Features, Bug fixes, Other changes.
-->

## [Unreleased]

### Other changes

- noodles 0.113 → 0.116 and noodles-bgzf 0.49 → 0.51, bumped together. The
  direct `noodles-bgzf` dependency exists only to enable its `libdeflate`
  feature; bumping one without the other puts two copies in the tree, and
  noodles then writes BAM through the copy without libdeflate. BAM output is
  byte-identical before and after. (0.117 needs Rust 1.91, above the 1.89 MSRV.)

- Production genome indexing now enables caps-sa 0.7's bounded geometric LCP
  memoization through its stable policy API. On the complete ruSTAR-shaped
  GRCh38 plus GENCODE v50 fixture (6.56 billion symbols, 6.18 billion retained
  suffixes, 1.40 million segments, 32 physical cores), the final caps-sa 0.7
  implementation built the SA in 172.953 s versus 267.592 s for its original
  0.7 baseline: 35.4% faster, with peak RSS reduced from 10,512,408 to
  9,169,892 KiB. The complete output hash was unchanged.

- The splice-junction sorts that produce `SJ.out.tab`, the `SJ` solo-feature
  rows and the `BySJout` survivor set now order on the whole key (chromosome,
  start, end, strand, motif) rather than on coordinates alone. The counts come
  from a `DashMap`, whose iteration order varies with hashing and with
  concurrent insertion, so a tie left to that order would have been a file
  that differs between runs or thread counts. `tests/determinism.rs` locks it:
  the same reads at one and at eight threads, and two runs at eight threads,
  produce byte-identical output in single-pass and two-pass mode. Answers #210.

- `cluster_seeds` reuses its window-bin map across reads on a thread instead
  of rebuilding it per read. Merging two windows re-keys every bin in the
  merged span, so the per-read pre-sizing was only a floor and the map
  rehashed; profiling a human 10x run put that rehash at 2.1% of on-CPU time.
  Output-neutral, verified by an empty diff on 20 M read pairs.

### Features

- **CLI and output parity: SAM/SJ/read-input knobs and the STAR limit
  surface** — 30 further STAR 2.7.11b parameters. (`--outSAMorder` came from #145.)

  - Implemented: `--outSAMmode` (`Full`/`NoQS`/`None`), `--outSJtype None`,
    `--outSJfilterReads Unique`, `--outSAMheaderHD`, `--outSAMheaderPG`,
    `--outSAMheaderCommentFile`, `--readFilesPrefix`, `--readNameSeparator`,
    `--readQualityScoreBase`, `--outQSconversionAdd`.
  - Refused loudly rather than accepted and ignored:
    `--outSAMfilter KeepOnlyAddedReferences` / `KeepAllAddedReferences` and
    `--readFilesType SAM`, both of which need machinery this aligner does not
    have.
  - Accepted and documented as output-neutral: the `--limit*` family,
    `--outTmpDir`, `--outTmpKeep`, `--runDirPerm`, `--genomeFileSizes`,
    `--outBAMsorting*`, `--readMatesLengthsIn`.

  `tests/parameter_surface.rs` now checks every STAR 2.7.11b parameter name
  against the CLI, so the coverage figure is machine-checked and surface drift
  fails a test.

### Bug fixes

- **Annotated junctions are stitched like STAR when a seed crosses them, and
  paired-end windows drop covered sub-alignments.** A seed that crosses an
  inserted (Gsj) junction is split into donor and acceptor halves. Window
  dedup merged the acceptor half into a genomic seed on the same diagonal
  (STAR keys it on the junction index, `WA_sjA`); the stitcher rejected a
  splice whose second exon was shorter than `alignSJoverhangMin` before
  checking whether the junction was annotated (`alignSJDBoverhangMin` is 3);
  and the two halves went through the junction scan, which on a junction
  whose boundary base repeats can settle on an equally scored unannotated
  position (STAR stitches them through the annotated junction directly). For
  pairs, a lower-scoring pair whose blocks are all covered by another pair
  in the same window was reported as an extra multimapper (`66M9S` beside
  `66M2I7M`); STAR drops it, as the single-end path already did. On 2M
  GEUVADIS ERR188021 pairs against chr21, 402 reads now get STAR's exact
  alignment and 4 lose it; annotated splices go from 7,957 to 8,197 (STAR:
  8,203). Runtime +2.5%.

- **`--quantMode TranscriptomeSAM` no longer aborts on reads soft-clipped at
  both ends of a single match block** (e.g. `1S97M2S`). With the default
  `--quantTranscriptomeSAMoutput BanSingleEnd_BanIndels_ExtendSoftclip`, the
  right clip was dropped when folding clips into the CIGAR, the record came out
  shorter than its sequence, and BAM encoding failed with "read
  length-sequence length mismatch", ending the run. Seen on human chr21
  paired-end data. From #283.

- **Unknown `--quantMode` values are rejected**, as STAR does
  (`Parameters.cpp:898-936`): only `TranscriptomeSAM`, `GeneCounts` or `-`. A
  typo such as `Genecounts` used to be ignored silently, producing no counts
  and no error. From #283.

- **`mappedFilter` ported as STAR has it** (`ReadAlign_mappedFilter.cpp`).
  The read is judged on its best alignment alone, in STAR's order: too short
  (score or matched bases), then too many mismatches, then too many loci; the
  whole set is kept or dropped together. Single-end filtered each alignment
  separately, so a secondary could be dropped on its own numbers.
  - Matched bases are STAR's `nMatch` (read equals genome), not the aligned
    length, which counted mismatches and `N` as matches and made
    `--outFilterMatchNminOverLread` looser than STAR's.
  - `--outFilterMismatchNoverLmax` divides by the mapped length (`rLength`),
    not the read length, in single-end and paired-end.
  - `--outFilterIntronMotifs` / `--outFilterIntronStrands` apply when a
    transcript is finalized and before the window's dedup, as in
    `stitchWindowAligns`, so a rejected transcript cannot evict one it covers,
    set the score range or become a chimeric segment. Paired-end did not apply
    them at all. `--outSAMstrandField intronMotif` also drops a spliced
    transcript whose strand is undefined, as STAR does.
  - Paired-end reads with more than `--outFilterMultimapNmax` loci are reported
    as "too many loci"; they were cleared inside the filter and fell through to
    "too short". A failing best pair now reports its own reason, and a pair
    with no transcript in any window is "other", as in single-end.

  Yeast 10k: alignments unchanged; PE unmapped types now match STAR read for
  read (`too many loci` 0 → 21, STAR 21). Human chr21 SE: multi / too many
  loci / too short 6051 / 221 / 2008 → 6114 / 234 / 1956 (STAR 6116 / 232 /
  1955). Supersedes #169 and #252.

- **Multi-member gzip input is no longer truncated.** Compressed input was
  decoded with `flate2::read::GzDecoder`, which stops at the end of the first
  gzip member; a `.gz` made of several concatenated members (bcl2fastq output,
  `cat a.fq.gz b.fq.gz`, any BGZF file) was read partially with no error and no
  warning. All four read paths now use `MultiGzDecoder`: FASTQ input, the solo
  barcode whitelist, solo counting, and the `emptydrops` binary.

- `Log.final.out` splits unmapped paired-end reads between `too short` and
  `other` again. STAR calls a read `other` when no good window was found at
  all and `too short` only when a window existed whose best transcript failed
  the score or length thresholds (`ReadAlign_mappedFilter.cpp`); the
  paired-end path reported every unmapped pair as `too short`, so the `other`
  bucket was permanently zero. On the nf-core/rnaseq test data (50 000 pairs,
  yeast chrI + GFP) this moves rustar from `too short 7374 / other 0` to
  `too short 3778 / other 3596`, against STAR's `3766 / 3609`. Closes #48.

- `--soloBarcodeMate 1` no longer runs without clipping. The flag says the
  barcode lives inside mate 1, and nothing else says how many bases that is,
  so STAR refuses the run unless the mate is clipped
  (`ParametersSolo.cpp:145-150`). Without the check the CB+UMI prefix was
  aligned as if it were cDNA (28 bases of it for 10x v3) and the run reported
  nothing. Closes #227.

- Read names are cut at `--readNameSeparator` (default `/`), as STAR does. A
  read named `foo/1` was previously emitted as `foo/1` where STAR emits `foo`.

- **STARsolo single-cell quantification (`--soloType`)** — the 10x
  Chromium / plate-based count-matrix pipeline, ported from STAR and
  verified against real STARsolo (#90).

  - Chemistries: `CB_UMI_Simple` (10x 3'/5'), `CB_UMI_Complex`
    (multi-segment CB), and `SmartSeq` (manifest-driven, SE read counts
    / PE fragment counts).
  - Features: `Gene`, `GeneFull` (pre-mRNA), `SJ` (splice-junction
    counts), and `Velocyto` (spliced/unspliced/ambiguous matrices for
    RNA velocity).
  - Barcode correction: all `--soloCBmatchWLtype` modes
    (`Exact`/`1MM`/`1MM_multi`/`1MM_multi_pseudocounts`/
    `1MM_multi_Nbase_pseudocounts`); all `--soloUMIdedup` and
    `--soloUMIfiltering` modes; real-valued `--soloMultiMappers`
    (`Uniform`/`PropUnique`/`Rescue`/`EM`).
  - Output: `raw/` + `filtered/` matrices, `EmptyDrops_CR` cell calling,
    CellRanger-style `Summary.csv`, and `--soloOutGzip`.

- **`--outSAMattributes GX GN`** — per-read gene-id / gene-name SAM tags
  for the solo `Gene` assignment.

- **`genomeGenerate` peak RSS cut from ~113 GB → ~11 GB** on the human
  genome (GRCh38, 32 threads). The construction pipeline no longer
  materialises three large intermediates that were dominating the peak:

  - The ACGT-only **kept-positions `Vec<u64>`** (~47 GB on the human
    genome) — caps-sa 0.5's `build_ext_mem_for_filter` API takes a
    predicate over text positions instead, and internally maintains
    only a ~770 MB bitmap + popcount prefix sum.
  - The **spacer-free copy of `genome.sequence`** (~6.3 GB) — the new
    `dispatch_caps_sa_segmented` hands the **original** spacer-bordered
    `&genome.sequence[..n2]` to caps-sa with a
    `SegmentedText::from_ends(spacer_positions)` limit provider, so
    LCP scans stop at the next spacer without a copy.
  - The **in-RAM SA `PackedArray`** (~25 GB) — the SA streams directly
    to `genome_dir/SA` via the new `PackedStreamWriter` as each caps-sa
    entry is emitted. `SuffixArray::build` (in-memory) is retained for
    tests; `sa_build::build_streaming` is the production path.

- **SAindex now parallelises across all rayon workers**. Previously
  every SAindex k-mer extraction sat on caps-sa's single-threaded
  phase-4 emit loop (~16 min of pure serial work). The new
  `SaIndex::build_parallel` reads the on-disk SA via chunked `pread`
  (so SA pages live in kernel page cache, **not** process RSS) and
  atomic-mins each k-mer's first-occurrence `sa_idx` into a
  `Vec<AtomicU64>` — the final pack into the SAindex's `PackedArray`
  is a single fast sequential pass.

- **SAindex inner loop now matches STAR's `isaStep + binary-search`
  skip algorithm** (`genomeSAindex.cpp::genomeSAindexChunk` /
  `funSAiFindNextIndex`). Consecutive SA entries share monotonically
  non-decreasing k-mer prefixes; rather than visit every entry, each
  rayon worker jumps forward by `isa_step = nSA / 4^nbases` (≈ 22 on
  human) and only stops to record k-mer boundaries — binary-searching
  inside the last `isa_step` window when `(indFull, iL4)` changes.
  Per-worker `ind0_local[il]` tracks the last-written k-mer index at
  each level so we skip writes inside a constant-prefix run; cross-
  chunk merge still uses `fetch_min` on the shared `Vec<AtomicU64>`.
  Drops the SAindex phase from ~10:50 to ~8 s on the human genome
  (~80× speedup), making rustar-aligner's full `genomeGenerate`
  **faster than upstream STAR 2.7.11b** (7:36 vs 11:26 wall, 32
  threads).

- **SAindex absent-slot encoding now matches STAR's** (`next_present
  _sa_idx | absent_mask` for between-present gaps, `n_entries
  | absent_mask` for tail-gaps). A single backward pass per SAindex
  level fills the gaps in place inside `firsts[]` before encoding
  into the output `PackedArray`. `hierarchical_lookup` is unchanged
  because it only consults the absent flag bit, not the slot's
  value; this change makes the on-disk SAindex bytes closer to
  STAR's (N flag bit is still not tracked — that's the only remaining
  divergence and a `hierarchical_lookup`-irrelevant one).

- **rayon thread pool now bound to `--runThreadN` at the
  run() dispatcher** rather than only inside `align_reads`. On a
  256-core machine this drops the rayon pool from 256 (the
  `num_cpus::get()` default) to whatever `--runThreadN` says,
  eliminating ~15 GB of glibc-arena-style allocator slack from
  256 worker-thread heaps.

- **mimalloc as the global allocator**. Lower per-allocation cost than
  glibc's malloc and per-thread heaps that return whole segments to
  the OS when abandoned, so allocator cache size stays bounded.

- **`--soloCellReadStats CB`** writes `Solo.out/<feature>/CellReads.stats`:
  one row per cell barcode with fifteen counters describing what
  happened to its reads — barcode match quality, unique or multi
  genomic mapping, feature assignment, exonic/intronic and their
  antisense counterparts, mitochondrial, and whether the read reached
  the matrix — plus the per-cell UMI and gene totals. Reads whose
  barcode never resolved are summed into a `CBnotInPasslist` row rather
  than dropped, so the columns account for the whole input.
  `--genomeChrSetMitochondrial` names the chromosomes behind the `mito`
  column.
- **`--runMode soloCellFiltering <raw dir> <output prefix>`** cell-calls
  an existing raw count matrix without aligning anything. Cell calling
  is a decision about a matrix, not about reads: re-calling with
  different `--soloCellFilter` parameters should not mean re-aligning,
  and a matrix produced elsewhere should be callable too. It streams
  the matrix into the same form the align path produces, so the filters
  are the identical code rather than a second implementation.
- **`--soloCellFilter EmptyDrops_CR` now uses CellRanger's actual
  statistics.** The ambient profile is smoothed with Simple Good-Turing,
  as CellRanger and STAR do, instead of an approximation that reserved
  unseen mass from the singleton rate and spread the remainder in
  proportion to raw counts. The Monte-Carlo null is drawn with libc++'s
  `std::mt19937` and `std::discrete_distribution`, seeded
  `19760110 * (isim + 1)` per simulation as STAR seeds it, replacing a
  SplitMix64 stream that could not agree with STAR's over an arbitrary
  number of draws. Cell calls move as a result.

- **`--soloFeatures Transcript3p`** quantifies transcripts rather than
  genes, using how far each read's 3' end sits from each transcript's.
  In a 3'-biased assay that distance discriminates between isoforms: a
  read 200 bases from the end of one and 4000 from the end of another
  is evidence for the first. The distance distribution is estimated
  from the data, then used as the likelihood in an EM over UMIs. Output
  is per *cluster* rather than per cell — `--soloClusterCBfile` (new,
  and required for this feature) says which cell is in which cluster,
  because one cell has too few UMIs to resolve isoforms. Reads sharing
  a UMI contribute the intersection of their transcript sets, not the
  union: they came from one molecule. Writes `matrix.mtx`,
  `features.tsv` and `transcriptEndDistanceDistribution.txt` under
  `Solo.out/Transcript3p/raw/`.
- **Chimeric multimapping detection (`--chimMultimapNmax`)**, STAR's
  newer enumeration path. Instead of pinning the best transcript and
  looking for one partner, it walks every transcript pair, keeps those
  within `--chimMultimapScoreRange` of the best chimeric score, and
  reports them all. A read whose fusion partner maps equally well to two
  places yields two junctions rather than an arbitrary one.
  `--chimNonchimScoreDropMin` gates the search on the linear alignment
  leaving enough of the read unexplained, and a read with more surviving
  loci than the cap reports none at all, as STAR does. The default of 0
  keeps the old single-best path.

- **`--chimFilter banGenomicN`** (STAR's default) drops a chimeric
  junction whose flanking genomic bases are not real bases: sequence
  around an assembly gap yields junctions that look clean by score and
  mean nothing. `--chimFilter None` keeps everything.

### Bug fixes

- `--runThreadN 1` ran on every logical core instead of on one. The
  rayon pool was configured only above 1, and skipping it leaves rayon's
  default of one worker per core. Output is unchanged; the run now uses
  the thread count asked for.
- `--soloUMIfiltering MultiGeneUMI_CR` kept every gene tied at the
  highest read count; CellRanger gives a tied UMI to no gene at all.
  Since one read per gene is the ordinary shape of a multi-gene UMI, the
  flag removed nothing in practice. On a 20k-read 10x fixture the count
  matrix moves from 16 465 to 15 414 against STAR's 15 423.
- `--soloUMIfiltering MultiGeneUMI_All` was aliased to `MultiGeneUMI`,
  which is neither STAR's behaviour nor the documented one: in STAR
  2.7.11b the variant is a no-op. It now removes a UMI seen in two or
  more genes from **all** of them, the behaviour the option name
  describes. Recorded in `DIVERGENCE.md` (closes #144).

- **STARsolo `Gene` assignment now requires exon concordance**, matching
  STARsolo: a read counts toward a gene only when every aligned block
  lies within the gene's exons, rather than merely overlapping one. This
  removes over-assignment of reads that extend past an exon boundary into
  an intron. On the mouse-chr19 differential test the `Gene` count matrix
  now matches STAR 2.7.11b (raw matrix identical apart from the
  multimapper tie-break tail that even a byte-faithful reimplementation
  cannot reproduce).
- **STARsolo `features.tsv` column 2 now emits the GTF `gene_name`**
  (symbol), with the STAR gene_id fallback, instead of duplicating the
  gene_id.
- **`--outSAMtype BAM SortedByCoordinate` no longer buffers the whole
  output in RAM.** The coordinate sort is now external: it fills a
  `--limitBAMsortRAM` buffer, spills sorted runs beside the output, and
  k-way merges them on finish. Peak memory is flat in output size
  (measured: 610 MB at a 64 MiB budget for 300 k through 2.4 M records,
  versus 648 MB → 2,132 MB before, a growth of 723 B/record that
  extrapolated to ~116 GB for a 160 M-record human sample). Output is
  byte-identical to the previous in-memory sort — coordinate ties still
  resolve to input order, since runs merge with the run index as
  tiebreak (verified at 2.4 M records through a 219-run multi-pass
  merge). `--limitBAMsortRAM N` now spills above `N` rather than
  aborting the run, and `--limitBAMsortRAM 0` means 512 MiB instead of
  "unlimited". Runs beyond 64 are merged in balanced passes so a small
  budget on a large run cannot exhaust file descriptors. Spill runs go to
  `--outTmpDir` when given (previously accepted and ignored), otherwise beside
  the output, and are always removed. See DIVERGENCE.md §4.3.

### Bumps

- `caps-sa` → `0.5` (adds `build_ext_mem_for_filter*`; see the caps-sa
  v0.5.0 release notes).

### Other

- New module `index::packed_stream` — bit-for-bit-compatible streaming
  writer for STAR's `PackedArray` format. Used by the streaming SA
  emit; documented to match `PackedArray::write` for `word_length ≤ 57`
  (the upper bound where the existing `PackedArray::write` is
  truncation-free; STAR's production widths are 32-37).

- New module-level `GenomeIndex::generate_streaming(params)` — the
  full `genomeGenerate` pipeline, used by the `genomeGenerate`
  run-mode dispatcher. The in-memory `GenomeIndex::build` +
  `GenomeIndex::write` flow remains for tests and any caller that
  needs random access to the SA in RAM.

- Removed `Transcript::read_seq`, a public field that was filled with a
  copy of the read at every finalised alignment and never read. **API
  removal.** Output is unchanged.

- The splice-motif check in the junction scan carries a sliding window
  and looks the motif up in a table, instead of re-reading four genome
  bases and matching on them at every position. Consecutive junction
  positions share two of the four bases, so the scan does two genome
  reads per position rather than four. Output is unchanged (SAM
  byte-identical on 200k reads); about 2.8% off the wall clock on a
  2M-read run, measured with `test/bench_ab.sh`.

Initial release of Rust rewrite of STAR.
