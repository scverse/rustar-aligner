---
title: Migrating from STAR
description: A practical guide to swapping STAR for rustar-aligner in an existing pipeline.
---

rustar-aligner is designed to be a drop-in replacement for STAR. In most cases you can substitute the binary, keep your existing parameters and genome indices, and the rest of the pipeline keeps working. This page walks through what's identical, what differs, and the small handful of things to watch for.

## What's identical

### Genome index format

The on-disk index format is the same. An index built by STAR can be loaded by rustar-aligner and vice-versa. There's no separate "build for rustar" step.

For the yeast benchmark genome, the suffix array is **byte-for-byte identical** to STAR's after Phase G3 (10,862 → 0 entry diffs). Other index files (`Genome`, `SAindex`, `chrName.txt`, etc.) are also identical or equivalent.

### Command-line parameters

Every parameter rustar-aligner supports uses STAR's exact `--camelCase` name. Defaults are STAR's defaults. Examples:

```bash
--runMode alignReads
--genomeDir /path/to/index
--readFilesIn reads_1.fq.gz reads_2.fq.gz
--readFilesCommand zcat
--runThreadN 16
--outSAMtype BAM SortedByCoordinate
--outFileNamePrefix sample_
--sjdbGTFfile gencode.gtf
--twopassMode Basic
--outFilterMultimapNmax 20
--alignIntronMax 1000000
```

These are the exact same flags STAR accepts.

### Output files

The set of files produced (`Log.final.out`, `Log.out`, `Log.progress.out`, `SJ.out.tab`, `Aligned.out.sam` / `Aligned.sortedByCoord.out.bam`, `Chimeric.out.junction`, `ReadsPerGene.out.tab`, `Aligned.toTranscriptome.out.bam`, `Unmapped.out.mate{1,2}`) and their formats match STAR's. Tools that consume STAR output — MultiQC, RSEM, htseq-count, Arriba, STAR-Fusion — work without modification.

### SAM/BAM records

For deterministic alignments the records match STAR exactly: same position, same CIGAR, same MAPQ, same NH/HI/AS/NM/nM/XS/jM/jI/MD tags, same proper-pair flag.

## What's different

### Tie-breaking

With the default `--outMultimapperOrder Old_2.4` no random numbers are involved: the alignments are written in STAR's order and the primary is the one STAR picks, so the records match STAR's.

The transcriptome BAM (`--quantMode TranscriptomeSAM`) picks its primary at random among a read's transcript alignments. rustar-aligner draws it from STAR's generator (Mersenne Twister, seeded with `--runRNGseed`) in STAR's read order, so it matches a STAR run with `--runThreadN 1`, whatever `--runThreadN` you give rustar-aligner. A multi-threaded STAR run seeds one generator per thread, and its picks depend on which thread maps which reads.

Under `--outMultimapperOrder Random`, rustar-aligner shuffles each read's equal-scoring alignments with its own per-read generator (seeded from `--runRNGseed` and the read name), so its order and primary differ from STAR's, and so do the transcriptome BAM primaries. The alignment set is the same: NH counts, AS scores and CIGAR strings are unchanged.

### CIGAR placement in homopolymers

In one read out of 10,000 in the yeast benchmark, rustar-aligner emits `100M1I45M4S` where STAR emits `108M1I37M4S` — same diagonal, same score, different insertion placement. This is a seed-level tie inside a homopolymer region. Real impact for downstream tools: effectively zero.

### Faithfulness numbers

On the yeast benchmark (10,000 reads, ERR12389696):

| Metric | Value |
|--------|-------|
| Single-end faithfulness (tie-adjusted) | **99.815%** (8,611 / 8,627 non-tie reads exact) |
| Paired-end faithfulness (tie-adjusted) | **99.883%** (16,284 / 16,306 mate alignments exact) |
| STAR-only reads | **0** |
| rustar-only reads | **0** |
| MAPQ inflations / deflations | **0 / 0** |
| NH tag diffs | **0** |
| Proper-pair diffs | **0** |

See the [STAR compatibility report](/rustar-aligner/reference/star-compatibility/) for the long-form analysis.

## What's not yet there

- **STARsolo** single-cell features (deferred)
- A small set of less-common parameters that haven't been ported yet — open an issue if you hit one

If a parameter is missing rustar-aligner errors out cleanly at startup rather than silently ignoring it.

## Migration checklist

1. **Build rustar-aligner** ([installation](/rustar-aligner/getting-started/installation/)).
2. **Reuse your existing genome index** — no rebuild needed. Or [generate a new one](/rustar-aligner/guides/genome-index/) if you want.
3. **Replace `STAR` with `rustar-aligner`** in your command lines. Keep all other flags.
4. **Run a dry comparison** on a sample input to verify outputs look as expected. Spot-check `Log.final.out`, the sorted BAM, and any quantification files against STAR's.
5. **Set `--runRNGseed`** explicitly if you need reproducible primary-alignment selection across runs.
6. **Update pipeline configuration** to point downstream tools at rustar-aligner's output (file paths and names match STAR's, so this is usually a no-op).

## Running both side-by-side

If you want STAR and rustar-aligner side-by-side during a transition period, write outputs to different directories with different `--outFileNamePrefix` values:

```bash
STAR --runMode alignReads --outFileNamePrefix star_run/sample_ ...
rustar-aligner --runMode alignReads --outFileNamePrefix rustar_run/sample_ ...
```

Then compare with your favourite SAM/BAM diff tool. The repository's `test/compare_sam.py` and `test/compare_pe.py` scripts are what we use internally.

## Reporting incompatibilities

If you find a divergence from STAR that isn't documented in [STAR compatibility](/rustar-aligner/reference/star-compatibility/), please [open an issue](https://github.com/scverse/rustar-aligner/issues) with:

- The STAR command line you ran
- A small input dataset (or pointer to one) that reproduces the divergence
- The differing output records (STAR's vs rustar-aligner's)

The project goal is to make this list as short as possible.
