---
title: Generating a genome index
description: Build a rustar-aligner / STAR-compatible genome index from a FASTA reference.
---

A genome index is the on-disk data structure rustar-aligner uses to look up where each read might map. You build it once per reference (and per parameter set) and then reuse it for every alignment.

The index format is **identical to STAR's**. An index built with rustar-aligner can be loaded by STAR and vice-versa. (After Phase G3 the suffix array is byte-for-byte identical for the yeast benchmark genome.)

## Minimal command

```bash
rustar-aligner --runMode genomeGenerate \
  --genomeDir /path/to/genome_index \
  --genomeFastaFiles /path/to/genome.fa
```

The `--genomeDir` directory must exist, but it should be empty (or contain a previous index that you want to overwrite).

## Realistic command

For real-world genomes with GTF annotations:

```bash
rustar-aligner --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir /path/to/genome_index \
  --genomeFastaFiles /path/to/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /path/to/gencode.v45.annotation.gtf \
  --sjdbOverhang 100 \
  --genomeSAindexNbases 14
```

## Key parameters

### `--genomeFastaFiles`

One or more FASTA files describing the reference. You can pass multiple files (e.g. one per chromosome):

```bash
--genomeFastaFiles chr1.fa chr2.fa chr3.fa ...
```

Sequences are concatenated in the order given. If you need to keep the order stable across runs (you do), be explicit about file order rather than using shell globs that may sort differently on different systems.

### `--sjdbGTFfile`

Optional GTF file with exon annotations. Including it at index time lets rustar-aligner pre-compute splice junction positions and store them in the index, which speeds up later alignments and improves accuracy on annotated junctions. If you don't have a GTF at index time you can still pass `--sjdbGTFfile` at alignment time and the junctions will be inserted dynamically.

GTF format options:

- `--sjdbGTFchrPrefix` — prepend a string to chromosome names from the GTF. Typically `chr` when the GTF uses bare numbers but the FASTA uses `chr1`, `chr2`, etc.
- `--sjdbGTFfeatureExon` — feature column value to treat as an exon (default: `exon`).
- `--sjdbGTFtagExonParentTranscript` — attribute name for the transcript ID (default: `transcript_id`).
- `--sjdbGTFtagExonParentGene` — attribute name for the gene ID (default: `gene_id`).

### `--sjdbOverhang`

The length of the splice junction database overhang on each side of a junction. Defaults to `100`. **Set this to `read_length - 1`** for best results — e.g. `99` for 100 bp reads, `149` for 150 bp reads, etc.

### `--genomeSAindexNbases`

Length of the SA pre-indexing string in log2 units. Defaults to `14`, which is appropriate for vertebrate-sized genomes. **For small genomes** (bacteria, yeast, *Drosophila*) **lower this** to avoid wasting RAM. STAR's recommendation is `min(14, log2(genomeLength)/2 - 1)`.

| Genome | Suggested `--genomeSAindexNbases` |
|--------|-----------------------------------|
| Human / mouse / vertebrate | 14 (default) |
| *Drosophila* | 12 |
| *Caenorhabditis elegans* | 11 |
| Yeast (*S. cerevisiae*) | 11 |
| Bacteria / viral | 9–10 |

### `--genomeChrBinNbits`

Log2 of the chromosome bin size used in the index. Defaults to `18` (256 kb). For genomes with many short scaffolds (e.g. fragmented assemblies) you can lower this to reduce padding.

### `--genomeSAsparseD`

Suffix array sparsity. `1` (default) is dense — fastest mapping, highest RAM. Higher values save RAM at the cost of some alignment speed.

## Output

After `genomeGenerate` finishes, `--genomeDir` will contain a set of files including `Genome`, `SA`, `SAindex`, `chrNameLength.txt`, `chrName.txt`, `chrStart.txt`, and (if a GTF was provided) `sjdbList.fromGTF.out.tab`, `transcriptInfo.tab`, `exonInfo.tab`, etc. These are the same files STAR produces.

## Build time and RAM

Index generation is the most expensive step. As a rough guide:

- **Yeast (~12 Mb)**: seconds, <1 GB RAM
- **Human GRCh38 (~3 Gb)**: 30–60 minutes, ~30 GB RAM with default parameters

Use `--runThreadN` to parallelise across CPU cores.
