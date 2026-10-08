#!/usr/bin/env bash
# Fetch the public data for the bulk total-RNA benchmark
# (--quantTranscriptomeUnspliced, --quantGeneSplicing Yes).
#
#   DATA=/path/to/bench bash scripts/bench_bulk_unspliced/fetch_data.sh
#
# Reference: GRCh38 primary assembly + GENCODE v50 comprehensive annotation
# (reference chromosomes) + matching transcript FASTA.
# Reads: 8 libraries of ENA project PRJEB57727 (whole blood in PAXgene,
# rRNA + globin depleted total RNA, 2x100 bp). Only the first
# 2 x NPAIRS read pairs of each run are streamed, then split into two
# pseudo-replicates of NPAIRS pairs each (first half / second half).
set -euo pipefail

DATA=${DATA:?set DATA to a writable directory}
NPAIRS=${NPAIRS:-2000000}
RUNS=${RUNS:-"ERR10501185 ERR10501288 ERR10501210 ERR10501321 ERR10501066 ERR10501182 ERR10501299 ERR10501099"}
GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50

mkdir -p "$DATA/ref" "$DATA/reads"
cd "$DATA/ref"
[ -n "${SKIP_REF:-}" ] || for f in GRCh38.primary_assembly.genome.fa.gz gencode.v50.annotation.gtf.gz gencode.v50.transcripts.fa.gz; do
  [ -s "$f" ] || curl -sfSLO "$GENCODE/$f"
done
if [ -z "${SKIP_REF:-}" ]; then
  [ -s GRCh38.primary_assembly.genome.fa ] || gunzip -k GRCh38.primary_assembly.genome.fa.gz
  [ -s gencode.v50.annotation.gtf ] || gunzip -k gencode.v50.annotation.gtf.gz
fi

cd "$DATA/reads"
lines=$((NPAIRS * 4))
for run in $RUNS; do
  [ -s "${run}_rep2_2.fq.gz" ] && continue
  urls=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp" | tail -1 | cut -f2)
  mate=1
  for u in ${urls//;/ }; do
    # head closes the pipe early; ignore curl's resulting write error.
    (curl -sfL "https://$u" || true) | gzip -dc 2>/dev/null | head -n $((2 * lines)) > "${run}_${mate}.fq" || true
    head -n $lines "${run}_${mate}.fq" | gzip -1 > "${run}_rep1_${mate}.fq.gz"
    tail -n +$((lines + 1)) "${run}_${mate}.fq" | gzip -1 > "${run}_rep2_${mate}.fq.gz"
    rm "${run}_${mate}.fq"
    mate=$((mate + 1))
  done
done
ls -la "$DATA/reads"
