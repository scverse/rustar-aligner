---
title: Bulk total RNA-seq
description: Spliced / unspliced gene counts and pre-mRNA-aware transcriptome projection for ribo-depleted bulk libraries.
---

Ribo-depleted ("total") RNA-seq libraries contain a large share of unspliced
pre-mRNA: intronic reads routinely make up a third to more than half of the
fragments. Two rustar-aligner options, both **opt-in and not part of STAR**,
target this kind of data:

- `--quantMode GeneSplicing` counts every uniquely mapped read or pair as
  spliced, unspliced or ambiguous per gene, with STARsolo's rules;
- `--quantTranscriptomePreMRNA BanRetainedIntron` stops
  `--quantMode TranscriptomeSAM` from handing unspliced pre-mRNA reads to
  retained-intron isoforms.

Without these options every output is the same as STAR's.

## Why pre-mRNA matters for transcript quantification

`--quantMode TranscriptomeSAM` projects each genomic alignment onto every
annotated transcript whose exons contain it. A read that lies inside an
intron is compatible with no fully spliced isoform, but GENCODE annotates many
`retained_intron` isoforms whose single exon covers that intron. The read is
therefore projected onto the retained-intron isoform only, and Salmon or RSEM
assign it there. In total RNA the number of such reads follows the pre-mRNA
content of the library, not the abundance of the retained-intron isoform, so
isoform proportions (and tximport's `avgTxLength` offsets) track library
quality rather than biology. The same effect is reported for Salmon with
genome decoys ([COMBINE-lab/salmon#1229](https://github.com/COMBINE-lab/salmon/issues/1229)).

## Spliced / unspliced / ambiguous gene counts (`GeneSplicing`)

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz --readFilesCommand zcat \
  --quantMode GeneCounts GeneSplicing \
  --sjdbGTFfile gencode.v50.annotation.gtf \
  --outFileNamePrefix sample_
```

It needs the same GTF-aware index as `TranscriptomeSAM` (built with
`--sjdbGTFfile` at `genomeGenerate`) and can be combined with any other
`--quantMode` value.

### Classification

The rules are those STARsolo uses for its single-cell spliced / unspliced
matrices, applied to each read (single-end) or read pair (paired-end) as if it
were one UMI:

1. Every annotated transcript that fully contains the alignment is tested.
   Each aligned block is called exonic, intronic or exon/intron-spanning,
   with a 6-base tolerance at exon boundaries. A spliced alignment
   that touches an intron is incompatible with that transcript, and a block in
   an intron longer than 1 Mb is not called intronic.
2. If the compatible transcripts belong to more than one gene, the read is
   `N_multiGene` and not counted.
3. Otherwise, per gene: only-exonic models give **spliced** (mature mRNA; the
   read need not cross a junction), intronic or spanning models with no
   only-exonic model give **unspliced** (pre-mRNA), and a mix gives
   **ambiguous**. A read inside an intron that another isoform of the same gene
   retains is ambiguous.

Unmapped, too-many-loci and multimapping reads are accounted as in
`GeneCounts`; paired-end reads with a single mapped mate count as unmapped.

### Output

`sample_ReadsPerGeneSplicing.out.tab`: a header line, then one line per gene
(in `geneInfo.tab` order) with nine counts: spliced, unspliced and ambiguous
for each strand convention.

```
gene_id  unstranded_spliced  unstranded_unspliced  unstranded_ambiguous  forward_spliced  ...  reverse_ambiguous
```

`forward` keeps transcripts on the strand of read 1 (`htseq-count -s yes`,
STARsolo `--soloStrand Forward`); `reverse` keeps transcripts on the opposite
strand (dUTP / TruSeq Stranded, `-s reverse`). Pick the columns matching the
library, as with `ReadsPerGene.out.tab`.

`sample_ReadsPerGeneSplicing.summary.tsv`: read accounting per strand
convention (`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_multiGene`,
`N_spliced`, `N_unspliced`, `N_ambiguous`) and the three fractions of the
assigned reads (`fraction_spliced`, `fraction_unspliced`,
`fraction_ambiguous`). The unspliced fraction is a direct per-library measure
of pre-mRNA content, usable as a QC metric or covariate.

## Pre-mRNA-aware transcriptome projection

```bash
rustar-aligner ... \
  --quantMode TranscriptomeSAM \
  --quantTranscriptomePreMRNA BanRetainedIntron
```

A **retained-intron interval** is an intron of one isoform that another
isoform of the same gene covers with a single exon reaching into both
flanking exons. With `BanRetainedIntron`, an alignment is not projected onto
the transcripts of a gene when:

- no mate crosses a splice junction, and
- one of its aligned blocks (after the soft-clip extension of
  `--quantTranscriptomeSAMoutput`) overlaps a retained-intron interval of that
  gene.

For paired-end data the rule applies to the fragment: the pair is dropped from
the gene as soon as one mate overlaps and neither mate is spliced. Cassette
exons and alternative 5'/3' splice sites do not create intervals, so their
reads are projected as before; so are reads that cross a junction. Other
alignments and other genes are untouched. The number of alignments and
projections removed is written to the log.

### What to expect

- Retained-intron isoforms no longer absorb pre-mRNA reads; their estimated
  share drops and stops following the library's unspliced fraction.
- The flip side: a retained-intron isoform that is genuinely expressed loses
  the reads that distinguish it (they are indistinguishable from pre-mRNA),
  and exonic pre-mRNA reads that it used to share now go to the other
  isoforms of the gene, as they do for genes with no retained-intron isoform.
  If retained introns are the object of study, keep the default and model
  intron retention explicitly.
- Fewer fragments reach Salmon / RSEM; use `ReadsPerGeneSplicing.summary.tsv`
  to report how much of the library is pre-mRNA.

The benchmark behind these statements (public whole-blood total RNA, GRCh38 +
GENCODE v50) can be reproduced with the scripts in
`scripts/bench_bulk_unspliced/`.
