#!/usr/bin/env bash
# Bulk total-RNA benchmark for --quantTranscriptomeUnspliced and
# --quantMode GeneSplicing. Run fetch_data.sh first.
#
#   DATA=/path/to/bench BIN=target/release/rustar-aligner THREADS=16 \
#     bash scripts/bench_bulk_unspliced/run.sh
#   python3 scripts/bench_bulk_unspliced/analyze.py "$DATA"
#
# For every pseudo-replicate library it runs rustar-aligner three times:
#   star/    : --quantMode TranscriptomeSAM GeneSplicing (STAR projection)
#   intron/  : TranscriptomeSAM + --quantTranscriptomeUnspliced Intron
#   premrna/ : TranscriptomeSAM + --quantTranscriptomeUnspliced PreMRNA
# then quantifies each transcriptome BAM with `salmon quant -a`, against the
# GENCODE transcripts (+ the <gene_id>-I sequences written by rustar-aligner).
# Wall time and peak RSS come from /usr/bin/time (-l on macOS, -v on Linux).
set -euo pipefail

DATA=${DATA:?set DATA}
BIN=$(cd "$(dirname "${BIN:-target/release/rustar-aligner}")" && pwd)/$(basename "${BIN:-target/release/rustar-aligner}")
THREADS=${THREADS:-16}
MODES=${MODES:-"star intron premrna"}
REF=$DATA/ref
IDX=$DATA/index
OUT=$DATA/out
if [ "$(uname)" = Darwin ]; then TIMEFLAG=-l; else TIMEFLAG=-v; fi

if [ ! -s "$IDX/SA" ]; then
  mkdir -p "$IDX"
  /usr/bin/time $TIMEFLAG "$BIN" --runMode genomeGenerate --runThreadN "$THREADS" \
    --genomeDir "$IDX" --genomeFastaFiles "$REF/GRCh38.primary_assembly.genome.fa" \
    --sjdbGTFfile "$REF/gencode.v50.annotation.gtf" --sjdbOverhang 49 \
    --outFileNamePrefix "$IDX/" > "$IDX/stdout.log" 2> "$IDX/time.log"
fi

# Salmon's alignment mode needs FASTA names equal to the BAM @SQ names
# (GENCODE transcript_id); drop the |-separated suffix of GENCODE headers.
TXFA=$REF/gencode.v50.transcripts.ids.fa
[ -s "$TXFA" ] || gzip -dc "$REF/gencode.v50.transcripts.fa.gz" | sed 's/|.*//' > "$TXFA"
mkdir -p "$DATA/targets"

for r1 in "$DATA"/reads/*_rep?_1.fq.gz; do
  lib=$(basename "$r1" _1.fq.gz)
  r2=${r1%_1.fq.gz}_2.fq.gz
  for mode in $MODES; do
    o=$OUT/$lib/$mode
    [ -s "$o/salmon/quant.sf" ] && continue
    mkdir -p "$o"
    fa=$DATA/targets/$mode.fa
    case $mode in
      star) extra=(--quantMode TranscriptomeSAM GeneSplicing) ;;
      intron) extra=(--quantMode TranscriptomeSAM --quantTranscriptomeUnspliced Intron) ;;
      premrna) extra=(--quantMode TranscriptomeSAM --quantTranscriptomeUnspliced PreMRNA) ;;
    esac
    # The unspliced sequences depend only on the index and the flank:
    # write them on the first run of each mode.
    if [ "$mode" != star ] && [ ! -s "$fa" ]; then
      extra+=(--quantTranscriptomeUnsplicedFasta Yes)
    fi
    /usr/bin/time $TIMEFLAG "$BIN" --runThreadN "$THREADS" --genomeDir "$IDX" \
      --readFilesIn "$r1" "$r2" \
      --outSAMtype None "${extra[@]}" \
      --outFileNamePrefix "$o/" > "$o/stdout.log" 2> "$o/time.log"
    if [ "$mode" = star ]; then
      [ -s "$fa" ] || ln -sf "$TXFA" "$fa"
    elif [ ! -s "$fa" ]; then
      cat "$TXFA" "$o/Aligned.toTranscriptome.unspliced.fa" > "$fa"
      rm "$o/Aligned.toTranscriptome.unspliced.fa"
    fi
    salmon quant -a "$o/Aligned.toTranscriptome.out.bam" -t "$fa" -l A -p "$THREADS" \
      -o "$o/salmon" > "$o/salmon.log" 2>&1
  done
done
