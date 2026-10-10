## What

Under `--soloOutLayout CellRanger` (on by default for 10x geometry), annotation, counting, cell calling and `metrics_summary.csv` now follow `cellranger count` 10.0.0, ported from its source rather than from observed tags:

- `src/solo/cr_annot.rs`: `TranscriptAnnotator::annotate_alignment` and `rescue_alignments_se`. Per alignment RE, TX, AN, GX, GN, fx, mm, xf, MAPQ 255 for confident reads; a read counts for a gene when its primary names exactly one.
- `src/solo/cr_dups.rs`: `mark_dups`. UMI correction (more reads, then larger UMI), the staged low-support test, representative read by `(UMI type, read name)`, duplicate flag 0x400, `umi_type` for `molecule_info.h5`.
- `src/solo/cr_cells.rs`: `ordmag` (NumPy `RandomState(0)` bootstrap) plus the non-ambient rescue (`default_rng(42)`, SGT ambient profile, BH at 0.001). NumPy's MT19937 `randint`, PCG64/SeedSequence and `argsort` are reproduced bit for bit.
- metrics: every field computed from those definitions (`Reads Mapped Confidently to Transcriptome` no longer 0.0%; `Valid UMI Sequences` over all reads; genes detected over called cells only; saturation and reads in cells over counted candidates).
- raw matrix holds every detected barcode, counted or not; `CB` is set on every read whose barcode resolves (multi-candidate barcodes resolved after the run), `UB` on every read with a valid UMI.
- `molecule_info.h5`: `umi_type` from the annotation, `barcodes` fixed width 43.
- 10x runs default to `--soloCellFilter EmptyDrops_CR`.

## Evidence

Full pbmc_1k_v3 run against `cellranger count` 10.0.0, same alignment code:

- `metrics_summary.csv`: 18 of 20 fields identical (the other two, `Total Genes Detected` and `Median UMI Counts per Cell`, follow the few reads that align differently).
- raw matrix 99.98% of entries identical; filtered matrix: the same 1,221 barcodes, 99.98% of entries identical.
- On reads whose alignment is identical (1.55M on chr21/22): MAPQ, RE, TX, AN, GX, GN, fx, CB 100%; xf and the duplicate flag 99.998%; UB 99.9999%.
- cell calling reproduces CellRanger's 1,221 barcodes exactly from CellRanger's own raw matrix.

## Notes

- The `Gene` and STARsolo paths are untouched; everything is gated on the CellRanger layout with a transcriptome in the index.
- CellRanger's `FDR` (0.001) and simulation count (100 000) are fixed for 3' v3.
