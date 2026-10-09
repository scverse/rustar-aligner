#!/usr/bin/env bash
# Compare the two suffix-array builders on one genome: wall clock, peak RSS,
# and whether the SA and SAindex they write are byte-identical.
#
#   test/bench_sa_builders.sh <genome.fa> <outdir> [saIndexNbases] [threads]
#
# Which builder runs is decided by --limitGenomeGenerateRAM, so the two arms
# differ only in that flag: 0 means "no limit" and always picks libsais, and a
# limit small enough that the estimate never fits always picks caps-sa.
set -euo pipefail

FASTA=${1:?usage: bench_sa_builders.sh <genome.fa> <outdir> [saIndexNbases] [threads]}
OUT=${2:?usage: bench_sa_builders.sh <genome.fa> <outdir> [saIndexNbases] [threads]}
NBASES=${3:-14}
THREADS=${4:-8}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/target/release/rustar-aligner"

# A genomeGenerate run is minutes long and uses every core, so an unrelated job
# does not merely add noise: it lands squarely on the arm unlucky enough to
# overlap it, and nothing in the resulting numbers shows that it happened.
LOAD=$(uptime | sed -E 's/.*load averages?: *([0-9.]+).*/\1/')
if awk -v l="$LOAD" 'BEGIN{exit !(l > 2.0)}'; then
  echo "load average is $LOAD; other work is running and the numbers would not mean anything." >&2
  echo "wait for the machine to be idle, or set BENCH_IGNORE_LOAD=1 to override." >&2
  [ "${BENCH_IGNORE_LOAD:-0}" = "1" ] || exit 1
fi

cargo build --release --manifest-path "$ROOT/Cargo.toml" >&2
mkdir -p "$OUT"

arm() { # $1 = label, $2 = --limitGenomeGenerateRAM value
  local dir="$OUT/$1"
  rm -rf "$dir"
  mkdir -p "$dir"
  /usr/bin/time -l "$BIN" --runMode genomeGenerate --runThreadN "$THREADS" \
      --genomeDir "$dir" --genomeFastaFiles "$FASTA" \
      --genomeSAindexNbases "$NBASES" --limitGenomeGenerateRAM "$2" \
      > "$OUT/$1.log" 2> "$OUT/$1.time"
  local wall rss
  wall=$(grep ' real' "$OUT/$1.time" | awk '{print $1}')
  rss=$(grep 'maximum resident set size' "$OUT/$1.time" | awk '{print $1}')
  printf "%-8s wall=%ss peakRSS=%.2f GB  (%s)\n" "$1" "$wall" \
      "$(echo "$rss/1073741824" | bc -l)" \
      "$(grep -o 'Building suffix array with [a-z-]*' "$OUT/$1.time" | head -1)"
}

arm libsais 0
arm capssa 1000000

if cmp -s "$OUT/libsais/SA" "$OUT/capssa/SA"; then echo "  SA identical"; else echo "  SA DIFFERS"; exit 1; fi
if cmp -s "$OUT/libsais/SAindex" "$OUT/capssa/SAindex"; then echo "  SAindex identical"; else echo "  SAindex DIFFERS"; exit 1; fi
