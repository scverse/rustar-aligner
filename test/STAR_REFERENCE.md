# STAR Source Code Reference

Quick reference for mapping rustar-aligner implementation to STAR's C++ codebase.

**STAR Repository:** https://github.com/alexdobin/STAR

---

## Module Mapping

### Genome Generation

| rustar-aligner | STAR |
|--------|------|
| `src/genome/mod.rs` | `source/Genome_genomeGenerate.cpp` |
| `src/genome/fasta.rs` | `source/SequenceFuns.cpp` |
| `src/index/suffix_array.rs` | `source/SuffixArrayFuns.cpp` |
| `src/index/sa_index.rs` | `source/genomeSAindex.cpp` |

**Key algorithms:**
- Genome padding: `Genome_genomeGenerate.cpp:genomeGenerate()`
- SA construction: `SuffixArrayFuns.cpp:suffixArrayGenome()`
- SAindex building: `genomeSAindex.cpp:genomeSAindexGenome()`

---

### Index Loading

| rustar-aligner | STAR |
|--------|------|
| `src/index/io.rs` | `source/Genome.cpp::genomeLoad()` |
| `src/index/packed_array.rs` | `source/PackedArray.h` |

**File formats:**
- `Genome`: 1 byte per base (values 0-5)
- `SA`: Variable-width packed array (typically 33 bits per entry)
- `SAindex`: 35-bit entries (32-bit position + 3 flag bits)

**Key differences:**
- STAR uses memory-mapped files (`mmap`) for large arrays
- rustar-aligner uses `memmap2` crate with same logic

---

### Seed Finding

| rustar-aligner | STAR |
|--------|------|
| `src/align/seed.rs` | `source/ReadAlign_maxMappableLength2strands.cpp` |
|  | `source/ReadAlign_mappedFilter.cpp` |

**Key algorithms:**
- MMP search: `ReadAlign_maxMappableLength2strands.cpp:maxMappableLength2strands()`
- Seed extension: `ReadAlign_mappedFilter.cpp:mappedFilter()`

**Parameters:**
- `--seedSearchStartLmax` → `seed_search_start_lmax` (default: 50)
- `--seedSearchLmax` → `seed_search_lmax` (default: 0 = automatic)
- `--seedMultimapNmax` → `seed_multimap_nmax` (default: 10000)
- `--winAnchorMultimapNmax` → `win_anchor_multimap_nmax` (default: 50)

**Known equivalence (verified Phase 4):**
- SAindex k-mer lookup → binary search SA → extend MMP
- Both implementations use identical suffix comparison logic

---

### Seed Clustering

| rustar-aligner | STAR |
|--------|------|
| `src/align/stitch.rs::cluster_seeds()` | `source/stitchWindowAligns.cpp::stitchWindowSeeds()` |

**Key algorithms:**
- Anchor seeds: Seeds with ≤ `winAnchorMultimapNmax` genome positions
- Clustering: Group seeds within `winBinNbits` proximity window (default: 100kb)
- Seed deduplication: Remove duplicate (read_pos, genome_pos) pairs

**Parameters:**
- `--winAnchorMultimapNmax` → `win_anchor_multimap_nmax` (default: 50)
- `--winBinNbits` → `win_bin_nbits` (default: 16, gives ~65kb window)
- `--seedNoneLociPerWindow` → `seed_none_loci_per_window` (default: 10)

**Known differences:**
- rustar-aligner caps expanded seeds at 200 per cluster (Phase 13 optimization)
- STAR does not have explicit cap (can cause O(n²) in repetitive regions)

---

### DP Stitching

| rustar-aligner | STAR |
|--------|------|
| `src/align/stitch.rs::stitch_seeds_dp()` | `source/stitchWindowAligns.cpp` |
| `src/align/score.rs` | `source/alignSmithWaterman.cpp` |

**Key algorithms:**
- DP chaining: Connect seeds with gap penalties
- Gap scoring: Distinguish indels vs splice junctions by gap length
- CIGAR generation: Build CIGAR string from DP traceback

**Parameters:**
- `--alignIntronMin` → `align_intron_min` (default: 21bp)
- `--alignIntronMax` → `align_intron_max` (default: 0 = automatic)
- `--alignSJoverhangMin` → `align_sj_overhang_min` (default: 5bp)
- `--scoreGap` → `score_gap` (default: 0)
- `--scoreGapNoncan` → `score_gap_noncan` (default: -8)
- `--scoreDelOpen` → `score_del_open` (default: -2)
- `--scoreDelBase` → `score_del_base` (default: -2)
- `--scoreInsOpen` → `score_ins_open` (default: -2)
- `--scoreInsBase` → `score_ins_base` (default: -2)

**Known equivalence (verified Phase 5):**
- Both use same gap penalty logic
- Splice junction detection: GT/AG (0), GC/AG (-4), AT/AC (-8), non-canonical (-8)

---

### Transcript Filtering

| rustar-aligner | STAR |
|--------|------|
| `src/align/read_align.rs::filter_transcripts()` | `source/ReadAlign_outputAlignments.cpp` |

**Parameters:**
- `--outFilterScoreMin` → `out_filter_score_min` (default: 0)
- `--outFilterScoreMinOverLread` → `out_filter_score_min_over_lread` (default: 0.66)
- `--outFilterMatchNmin` → `out_filter_match_nmin` (default: 0)
- `--outFilterMatchNminOverLread` → `out_filter_match_nmin_over_lread` (default: 0.66)
- `--outFilterMismatchNmax` → `out_filter_mismatch_nmax` (default: 10)
- `--outFilterMismatchNoverLmax` → `out_filter_mismatch_nover_lmax` (default: 0.3)
- `--outFilterMultimapNmax` → `out_filter_multimap_nmax` (default: 10)

**Logic:**
1. Filter by score threshold
2. Filter by match count
3. Filter by mismatch count
4. Select best alignments (up to `outFilterMultimapNmax`)

---

### MAPQ Calculation

| rustar-aligner | STAR |
|--------|------|
| `src/mapq.rs` | `source/ReadAlign_outputAlignments.cpp` |

**Formula:**
- Unique mapping (1 alignment): MAPQ = 255
- Multi-mapping (n alignments): MAPQ = -10 * log10(1 - 1/n)
- Capped at 255

**Known equivalence (verified Phase 6):**
- Both implementations use identical formula

---

### SAM Output

| rustar-aligner | STAR |
|--------|------|
| `src/io/sam.rs` | `source/ReadAlign_outputAlignments.cpp` |
|  | `source/ReadAlign_outputTranscriptSAM.cpp` |

**SAM tags (STAR vs rustar-aligner):**
- `NH:i` (number of hits): STAR ✓, rustar-aligner ✗ (deferred)
- `HI:i` (hit index): STAR ✓, rustar-aligner ✗ (deferred)
- `AS:i` (alignment score): STAR ✓, rustar-aligner ✗ (deferred)
- `NM:i` (edit distance): STAR ✓, rustar-aligner ✗ (deferred)
- `nM:i` (mismatches): STAR ✓, rustar-aligner ✗ (deferred)

**Reason for deferral:** noodles crate lifetime complexity

---

### GTF Parsing

| rustar-aligner | STAR |
|--------|------|
| `src/junction/gtf.rs` | `source/sjdbLoadFromFiles.cpp` |
| `src/junction/mod.rs` | `source/SjdbClass.cpp` |

**Key algorithms:**
- Parse GTF exon features
- Group exons by transcript
- Extract intron coordinates (consecutive exon boundaries)
- Build junction database (HashMap for O(1) lookup)

**Parameters:**
- `--sjdbGTFfile` → `sjdb_gtf_file`
- `--sjdbScore` → `sjdb_score` (default: 2, bonus for annotated junctions)

**Known equivalence (verified Phase 7):**
- Both use 1-based coordinates for junction boundaries
- Both apply `sjdbScore` bonus to annotated junctions

---

### SJ.out.tab Output

| rustar-aligner | STAR |
|--------|------|
| `src/junction/sj_output.rs` | `source/outputSJ.cpp` |

**Format (9 columns):**
1. Chromosome
2. First base of intron (1-based)
3. Last base of intron (1-based)
4. Strand (0=undefined, 1=+, 2=-)
5. Intron motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=N-N)
6. 0=unannotated, 1=annotated
7. Number of uniquely mapping reads
8. Number of multi-mapping reads
9. Maximum overhang

**Known differences:**
- Overhang calculation: rustar-aligner uses placeholder (5bp), STAR computes exact value
- Thread-safe accumulation: rustar-aligner uses DashMap, STAR uses locks

---

### Threading

| rustar-aligner | STAR |
|--------|------|
| `src/lib.rs::run_with_threads()` | `source/STAR.cpp::main()` |
|  | `source/ThreadControl.cpp` |

**Parallelization:**
- Read-level parallelism: Each thread processes subset of reads
- Shared index: Genome, SA, SAindex shared via Arc (read-only)
- Thread-safe stats: DashMap for junction counts, Mutex for alignment stats

**Parameters:**
- `--runThreadN` → `run_thread_n` (default: 1)

**Known differences:**
- STAR uses custom thread pool with load balancing
- rustar-aligner uses rayon for simpler parallelization

---

### Two-Pass Mode

| rustar-aligner | STAR |
|--------|------|
| `src/lib.rs::run_two_pass()` | `source/STAR.cpp::main()` (--twopassMode Basic) |

**Workflow:**
1. **Pass 1:** Align reads, collect junction stats, discard alignments
2. **Filter junctions:** Keep junctions with ≥1 unique OR ≥2 multi reads
3. **Pass 2:** Re-align ALL reads with merged junction DB (GTF + novel)

**Parameters:**
- `--twopassMode` → `twopass_mode` (None | Basic)
- `--twopass1readsN` → `twopass1_reads_n` (default: -1 = all reads)

**Output files:**
- `SJ.pass1.out.tab`: Pass 1 junctions
- `SJ.out.tab`: Final junctions (pass 2)
- `Aligned.out.sam/bam`: Pass 2 alignments only

**Known equivalence (verified Phase 11):**
- Both implementations use identical filtering logic
- Both re-align ALL reads in pass 2 (not just subset)

---

### Chimeric Alignment Detection

| rustar-aligner | STAR |
|--------|------|
| `src/chimeric/detect.rs` | `source/ReadAlign_chimericDetection.cpp` |
| `src/chimeric/output.rs` | `source/ReadAlign_outputTranscriptChimeric.cpp` |

**Detection tiers:**
1. **Tier 1 (soft-clip):** Detect from soft-clipped reads (>20% clipped)
2. **Tier 2 (multi-cluster):** Detect from multi-locus seed clusters
3. **Tier 3 (re-mapping):** Re-map soft-clipped regions (NOT YET IMPLEMENTED in rustar-aligner)

**Parameters:**
- `--chimSegmentMin` → `chim_segment_min` (default: 0 = disabled, typical: 20)
- `--chimJunctionOverhangMin` → `chim_junction_overhang_min` (default: 20)
- `--chimScoreMin` → `chim_score_min` (default: 0)
- `--chimOutType` → `chim_out_type` (Junctions | SeparateSAMold | WithinBAM)

**Output file:** `Chimeric.out.junction` (14 columns)

**Known limitations:**
- rustar-aligner: Single-end only, tiers 1-2
- STAR: Single-end and paired-end, tiers 1-3

---

## Parameter Name Mapping

STAR uses `--camelCase` CLI parameters. rustar-aligner mirrors these exactly.

### Commonly Used Parameters

| STAR | rustar-aligner Field | Default |
|------|--------------|---------|
| `--genomeDir` | `genome_dir` | (required) |
| `--readFilesIn` | `read_files_in` | (required) |
| `--readFilesCommand` | `read_files_command` | `"zcat"` |
| `--outFileNamePrefix` | `out_file_name_prefix` | `"./"` |
| `--outSAMtype` | `out_sam_type` | `["SAM"]` |
| `--runThreadN` | `run_thread_n` | `1` |
| `--sjdbGTFfile` | `sjdb_gtf_file` | `""` |
| `--sjdbScore` | `sjdb_score` | `2` |
| `--twopassMode` | `twopass_mode` | `None` |
| `--chimSegmentMin` | `chim_segment_min` | `0` |

See `src/params.rs` for complete parameter list (~40 parameters).

---

## Known Algorithmic Differences

### Phase 13 Optimizations (rustar-aligner-specific)

1. **Lazy seed position iterator:** Avoids Vec allocation per seed
2. **PackedArray fast-path:** Direct 8-byte slice read for aligned positions
3. **Binary search position_to_chr:** O(log n) vs O(n) linear scan
4. **Expanded seed cap:** 200 seeds per cluster (prevents O(n²) in repetitive regions)
5. **Deferred CIGAR clone:** Build CIGAR once for best alignment only

**Impact:** 3.9x speedup (3s → 0.77s for 1000 reads), no accuracy loss

### SAM Optional Tags

**Not yet implemented in rustar-aligner:**
- `NH:i` (number of hits)
- `HI:i` (hit index)
- `AS:i` (alignment score)
- `NM:i` (edit distance)

**Reason:** noodles crate lifetime complexity

**Workaround:** Can be added post-hoc with samtools or custom script

---

## Key Bug Patterns

### Reverse Strand Genome Offset

**Issue:** Genome has 2 copies: forward [0, n_genome) and reverse-complement [n_genome, 2*n_genome)

**Rule:** ANY code accessing genome bases for reverse-strand reads MUST add `n_genome` offset

**Example:**
```rust
// WRONG
let base = genome.get_base(pos);

// CORRECT
let offset = if strand == Strand::Reverse { genome.n_genome } else { 0 };
let base = genome.get_base(pos + offset);
```

**Verified locations:**
- ✓ `src/align/stitch.rs::count_mismatches()` (fixed Phase 12.2)
- ✓ `src/align/read_align.rs::align_read()` (seed finding)
- ✗ `src/align/score.rs::detect_splice_motif()` (low priority, deferred)

---

## Investigation Workflow

When rustar-aligner output differs from STAR:

1. **Identify discrepancy type:**
   - Different chromosome → Seed clustering
   - Different position → Tie-breaking
   - Different CIGAR → DP stitching
   - Different mapping status → Filtering thresholds

2. **Extract discrepant reads:**
   ```bash
   ./test/investigate.sh test/results/.../comparison/alignment_diff.txt
   ```

3. **Re-run with debug logging:**
   ```bash
   RUST_LOG=debug ./target/release/rustar-aligner --runMode alignReads ...
   ```

4. **Compare STAR source:**
   - Find relevant STAR C++ file (see table above)
   - Compare algorithm logic
   - Check parameter defaults

5. **Test hypothesis:**
   - Adjust parameters
   - Re-run comparison
   - Validate fix

---

## Useful STAR Source Files

**Core alignment:**
- `source/ReadAlign.cpp` - Main alignment driver
- `source/ReadAlign_mappedFilter.cpp` - Seed finding
- `source/stitchWindowAligns.cpp` - Seed clustering and DP stitching
- `source/alignSmithWaterman.cpp` - Local alignment

**Output:**
- `source/ReadAlign_outputAlignments.cpp` - SAM output and filtering
- `source/ReadAlign_outputTranscriptSAM.cpp` - SAM record generation
- `source/outputSJ.cpp` - SJ.out.tab output

**Index:**
- `source/Genome_genomeGenerate.cpp` - Genome generation
- `source/SuffixArrayFuns.cpp` - Suffix array construction
- `source/genomeSAindex.cpp` - SAindex building
- `source/Genome.cpp` - Index loading

**Junctions:**
- `source/sjdbLoadFromFiles.cpp` - GTF parsing
- `source/SjdbClass.cpp` - Junction database

**Chimeric:**
- `source/ReadAlign_chimericDetection.cpp` - Chimeric detection
- `source/ReadAlign_outputTranscriptChimeric.cpp` - Chimeric output

---

## Version Information

**STAR version tested against:** 2.7.11b (latest stable as of 2026-02-07)

**rustar-aligner phases completed:** 1-13 (Phases 1-11 = feature parity, 12 = chimeric, 13 = optimization)

**Next phase:** 14 (STARsolo single-cell features)
