---
title: Bulk total RNA-seq
description: Unspliced targets in the transcriptome BAM, a splicing-status tag and spliced / unspliced gene counts for ribo-depleted bulk libraries.
---

Ribo-depleted ("total") RNA-seq libraries contain a large share of unspliced
pre-mRNA: intronic reads routinely make up a third to more than half of the
fragments. rustar-aligner has three **opt-in** options for this kind of data.
None of them exists in STAR; without them every output is the same as STAR's.

| Option | What it adds |
|---|---|
| `--quantTranscriptomeUnspliced Intron` or `PreMRNA` | one unspliced target per gene, `<gene_id>-I`, in `Aligned.toTranscriptome.out.bam` |
| `--outSAMattributes ... sp` | an `sp:A` splicing-status tag on every genomic alignment |
| `--quantMode GeneSplicing` | spliced / unspliced / ambiguous read counts per gene |

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

The fix used for single-cell data by
[splici / alevin-fry](https://doi.org/10.1038/s41592-022-01408-3) is to give
the quantifier unspliced targets next to the spliced transcripts, so that a
read compatible with both is shared by the EM instead of being forced onto the
retained-intron isoform. `--quantTranscriptomeUnspliced` does this for the
transcriptome BAM.

## Unspliced targets in the transcriptome BAM

```bash
rustar-aligner \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --quantMode TranscriptomeSAM \
  --quantTranscriptomeUnspliced PreMRNA \
  --quantTranscriptomeUnsplicedFasta Yes \
  --outFileNamePrefix sample_
```

### Targets

One target per gene, named `<gene_id>-I`, is added after the annotated
transcripts (`@SQ` lines in gene order):

- `Intron` (splici-style): the union of the introns of all the gene's
  isoforms, merged, extended on each side by
  `--quantTranscriptomeUnsplicedFlank` bases (default `-1` = the index's
  `sjdbOverhang`, i.e. read length - 1 by STAR's convention), clipped to the
  gene body, merged again and concatenated in genome order. Introns retained
  by an isoform and introns skipped by a cassette exon are included.
- `PreMRNA`: the whole gene body, from the first to the last annotated base.

A gene is taken on the chromosome and strand of its first transcript. In
`Intron` mode a single-exon gene has no target. Targets are built at
alignment time from the transcript tables of the index, so any GTF-aware
index works and the flank can change between runs.

### Projection

Nothing changes in how an alignment is projected: it goes to every target that
contains all its aligned blocks, spliced and unspliced alike. So:

- a read inside a retained intron is written both to the retained-intron
  isoform and to `<gene_id>-I`; `NH`, `HI` and `MAPQ` count all targets, and
  Salmon or RSEM decide;
- a read inside a constitutive intron is written to `<gene_id>-I` only
  (it used to be dropped);
- a read or pair that **crosses a splice junction** is processed RNA and goes
  to spliced targets only;
- `--quantTranscriptomeSAMoutput` rules (indels, soft-clip extension,
  single-end) apply to every target.

`Intron` or `PreMRNA`? `PreMRNA` also offers a home to the exonic part of
pre-mRNA and keeps pairs with one mate in an exon and one in an intron (in
`Intron` mode such a pair fits no target unless an isoform retains that
intron). `Intron` matches splici and keeps exon-only reads away from the
unspliced targets, at the cost of counting the exonic part of pre-mRNA as
mature.

### Files

- `sample_Aligned.toTranscriptome.targets.tsv`: `target_id`, `gene_id`,
  `gene_name`, `status` (`spliced` / `unspliced`), `length`, in `@SQ` order;
  a ready-made tx2gene table (sum by `gene_id` and `status` for spliced and
  unspliced gene counts).
- `sample_Aligned.toTranscriptome.unspliced.fa` (with
  `--quantTranscriptomeUnsplicedFasta Yes`): the `<gene_id>-I` sequences, in
  transcript orientation. It is of the order of the genome size and depends
  only on the index and the flank, so write it once.

### Salmon

Salmon's alignment mode needs every `@SQ` target in its FASTA, with the same
names. With GENCODE:

```bash
gzip -dc gencode.v50.transcripts.fa.gz | sed 's/|.*//' > transcripts.fa
cat transcripts.fa sample_Aligned.toTranscriptome.unspliced.fa > targets.fa
salmon quant -a sample_Aligned.toTranscriptome.out.bam -t targets.fa -l A -o salmon_out
```

Keep the `-I` targets out of isoform-level analyses, and out of tximport's
`avgTxLength` if a spliced-only gene length is wanted.

## Splicing-status tag (`sp`)

`--outSAMattributes Standard sp` adds `sp:A:S` (spliced: compatible only with
mature mRNA, the read need not cross a junction), `sp:A:U` (unspliced: needs
pre-mRNA) or `sp:A:A` (ambiguous: both) to every genomic alignment record,
from the rules below applied to every annotated transcript that contains the
alignment, on either strand and whatever its gene. There is no tag when no
transcript contains the alignment. Both mates of a pair carry the fragment
status. `sp` is in no preset and is not used by STAR or STARsolo.

## Spliced / unspliced gene counts (`GeneSplicing`)

```bash
rustar-aligner ... --quantMode GeneCounts GeneSplicing
```

A cheap per-library measure of pre-mRNA content, and a gene-level table. Each
uniquely mapped read or pair is classified with the rules STARsolo uses for its
single-cell spliced / unspliced matrices, one read standing for one molecule:

1. Every annotated transcript that fully contains the alignment is tested.
   Each aligned block is called exonic, intronic or exon/intron-spanning, with
   a 6-base tolerance at exon boundaries. A spliced alignment that touches an
   intron is incompatible with that transcript, and a block in an intron longer
   than 1 Mb is not called intronic.
2. If the compatible transcripts belong to more than one gene, the read is
   `N_multiGene` and not counted.
3. Otherwise: only-exonic models give **spliced**, intronic or spanning models
   with no only-exonic model give **unspliced**, and a mix gives
   **ambiguous** (for example a read inside an intron that another isoform
   retains).

Unmapped, too-many-loci and multimapping reads are accounted as in
`GeneCounts`; pairs with a single mapped mate count as unmapped.

- `sample_ReadsPerGeneSplicing.out.tab`: a header line, then one line per gene
  (in `geneInfo.tab` order) with spliced, unspliced and ambiguous counts for
  each strand convention (`unstranded_*`, `forward_*`, `reverse_*`).
  `forward` keeps transcripts on the strand of read 1 (`htseq-count -s yes`);
  `reverse` keeps the opposite strand (dUTP / TruSeq Stranded,
  `-s reverse`).
- `sample_ReadsPerGeneSplicing.summary.tsv`: `N_unmapped`, `N_multimapping`,
  `N_noFeature`, `N_multiGene`, `N_spliced`, `N_unspliced`, `N_ambiguous` and
  the three fractions of the assigned reads, per strand convention.

The benchmark behind these options (public whole-blood total RNA, GRCh38 +
GENCODE v50) can be reproduced with the scripts in
`scripts/bench_bulk_unspliced/`.
