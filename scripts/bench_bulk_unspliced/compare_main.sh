#!/usr/bin/env bash
# Check that this branch, run WITHOUT the new options, reproduces a main
# build's outputs on real data (and that adding --quantGeneSplicing /
# --quantTranscriptomeUnspliced None changes none of them).
#
#   DATA=... BIN_MAIN=/path/to/main/rustar-aligner BIN=target/release/rustar-aligner \
#     LIB=ERR10501185_rep1 bash scripts/bench_bulk_unspliced/compare_main.sh
#
# Compares Aligned.out.bam records and header (minus @PG/@CO, which carry the
# command line), transcriptome BAM records, SJ.out.tab and Log.final.out
# (minus wall-clock lines). Needs samtools. The transcriptome BAM uses
# --quantTranscriptomeSAMoutput BanSingleEnd: with the default soft-clip
# extension, main aborts on reads such as 1S97M2S (fixed on this branch).
set -euo pipefail

DATA=${DATA:?}; BIN_MAIN=${BIN_MAIN:?}; BIN=${BIN:-target/release/rustar-aligner}
LIB=${LIB:-ERR10501185_rep1}; THREADS=${THREADS:-16}
OUT=$DATA/compare_main/$LIB
r1=$DATA/reads/${LIB}_1.fq.gz; r2=$DATA/reads/${LIB}_2.fq.gz

run() { # bin outdir extra...
  local bin=$1 o=$2; shift 2
  mkdir -p "$o"
  "$bin" --runThreadN "$THREADS" --genomeDir "$DATA/index" --readFilesIn "$r1" "$r2" \
    --outSAMtype BAM Unsorted --quantTranscriptomeSAMoutput BanSingleEnd \
    --outFileNamePrefix "$o/" "$@" > "$o/stdout.log" 2>&1
}
run "$BIN_MAIN" "$OUT/main" --quantMode TranscriptomeSAM
run "$BIN" "$OUT/branch" --quantMode TranscriptomeSAM
run "$BIN" "$OUT/branch_new" --quantMode TranscriptomeSAM --quantGeneSplicing Yes \
  --quantTranscriptomeUnspliced None

digest() { # dir
  local d=$1
  {
    samtools view -H "$d/Aligned.out.bam" | grep -v -e '^@PG' -e '^@CO'
    samtools view "$d/Aligned.out.bam"
  } | shasum | cut -c1-16
  samtools view "$d/Aligned.toTranscriptome.out.bam" | shasum | cut -c1-16
  shasum < "$d/SJ.out.tab" | cut -c1-16
  grep -v -e Started -e Finished -e 'Mapping speed' "$d/Log.final.out" | shasum | cut -c1-16
}
status=0
for v in branch branch_new; do
  if diff <(digest "$OUT/main") <(digest "$OUT/$v") > /dev/null; then
    echo "$LIB: main == $v (Aligned.out.bam, toTranscriptome, SJ, Log.final.out)"
  else
    echo "$LIB: main != $v"; status=1
  fi
done
echo "records: $(samtools view -c "$OUT/main/Aligned.out.bam") genome, $(samtools view -c "$OUT/main/Aligned.toTranscriptome.out.bam") transcriptome"
exit $status
