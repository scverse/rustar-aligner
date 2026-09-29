# Bulk total-RNA benchmark

Public-data benchmark for `--quantTranscriptomeUnspliced` and
`--quantMode GeneSplicing` (see the "Bulk total RNA-seq" guide).

- Reference: GRCh38 primary assembly + GENCODE v50 comprehensive annotation
  (reference chromosomes) + GENCODE v50 transcript FASTA.
- Reads: ENA project PRJEB57727 (whole blood in PAXgene, rRNA + globin
  depleted total RNA, reverse-stranded, 2x100 bp), runs ERR10501066,
  ERR10501099, ERR10501182, ERR10501185, ERR10501210, ERR10501288,
  ERR10501299, ERR10501321. From each run the first 4 M read pairs are split
  into two pseudo-replicates of 2 M pairs (`rep1`, `rep2`). Pseudo-replicates
  only share library prep and sequencing, so they measure sampling noise, not
  technical-replicate variation.
- The index is built with `--sjdbOverhang 49` (the ENA metadata suggested
  2x50 bp; the reads are 2x100 bp). This is valid but not STAR's recommended
  read length - 1, and it sets the default `Intron` flank to 49 bases.
- Tools used for the committed results: salmon 2.8.0, samtools 1.24,
  Python 3 with numpy + scipy; macOS aarch64, 16 threads, 128 GB RAM.

```bash
export DATA=/path/to/bench
bash scripts/bench_bulk_unspliced/fetch_data.sh             # ~6 GB download
cargo build --release
BIN=target/release/rustar-aligner THREADS=16 bash scripts/bench_bulk_unspliced/run.sh
python3 scripts/bench_bulk_unspliced/analyze.py "$DATA"     # results_PRJEB57727.txt
BIN_MAIN=/path/to/main/rustar-aligner LIB=ERR10501185_rep1 \
  bash scripts/bench_bulk_unspliced/compare_main.sh         # outputs identical to main
```

`run.sh` aligns every library three times (STAR projection, `Intron`,
`PreMRNA`) and runs `salmon quant -a` on each transcriptome BAM.
`results_PRJEB57727.txt` is the output of `analyze.py` for the run behind the
PR.
