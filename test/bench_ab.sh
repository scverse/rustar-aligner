#!/usr/bin/env bash
# Interleaved A/B of two rustar-aligner binaries.
#
#   test/bench_ab.sh <binA> <binB> <genomeDir> <reads.fastq> [threads] [rounds]
#
# Prints one line per round and a median over the rounds that were measured on
# a quiet machine, plus a count of the ones that were not. Read the whole
# output: a median over four surviving rounds out of six means something
# different from a median over six.
#
#     BENCH_MIN_IDLE=88   percent CPU idle a round must hold to count
#     BENCH_MODE="None"   --outSAMtype words (quote the two-word forms)
#     BENCH_ARGS=""       extra flags passed to both sides
#
# ---------------------------------------------------------------------------
# Why this script exists, and what it refuses to do
#
# Three things went wrong often enough while measuring this aligner that they
# are worth encoding rather than remembering.
#
# 1. A busy machine does not add noise, it inverts results. An unrelated job
#    saturating the cores lands on whichever side happens to run under it, and
#    nothing in the timings shows that it happened: the series just drifts. So
#    idle is sampled *before and after every single run*, not once at the
#    start, and a round where either sample falls below the threshold is
#    dropped and counted rather than quietly folded into the median.
#
#    The check is on CPU idle, not load average. Load average is an
#    exponential average over minutes, so it stays high long after the
#    offending job is gone and refuses to measure on a machine that is now
#    perfectly quiet.
#
# 2. Timing a total hides the part that changed. If the change is in the
#    writer, run one configuration with `--outSAMtype None` alongside the BAM
#    one: only the difference is the work under test. On yeast, BAM writing is
#    1-4% of the run, so a total dominated by alignment cannot resolve a change
#    to it at any reachable scale.
#
# 3. Both sides must see byte-identical argv. The `@PG` `CL:` line records the
#    command line verbatim, so running `./old` against `./new`, or writing to
#    different `--outFileNamePrefix` directories, makes the output differ in
#    the header and the compressed size with it. Each side is therefore copied
#    to `rustar-aligner` inside its own directory and run with prefix `./`.
# ---------------------------------------------------------------------------
set -euo pipefail

BIN_A=${1:?usage: bench_ab.sh <binA> <binB> <genomeDir> <reads.fastq> [threads] [rounds]}
BIN_B=${2:?usage: bench_ab.sh <binA> <binB> <genomeDir> <reads.fastq> [threads] [rounds]}
GENOME_DIR=${3:?usage: bench_ab.sh <binA> <binB> <genomeDir> <reads.fastq> [threads] [rounds]}
READS=${4:?usage: bench_ab.sh <binA> <binB> <genomeDir> <reads.fastq> [threads] [rounds]}
# Each round runs from a private temp directory, so every path handed to the
# aligner has to be absolute or it will not resolve from there.
GENOME_DIR=$(cd "$GENOME_DIR" && pwd)
READS=$(cd "$(dirname "$READS")" && pwd)/$(basename "$READS")

THREADS=${5:-16}
ROUNDS=${6:-6}

MIN_IDLE=${BENCH_MIN_IDLE:-88}
MODE=${BENCH_MODE:-None}
EXTRA=${BENCH_ARGS:-}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

NAME_A=$(basename "$BIN_A")
NAME_B=$(basename "$BIN_B")
[ "$NAME_A" = "$NAME_B" ] && { NAME_A="A:$NAME_A"; NAME_B="B:$NAME_B"; }

for side in a b; do
  mkdir -p "$WORK/$side"
done
cp "$BIN_A" "$WORK/a/rustar-aligner"
cp "$BIN_B" "$WORK/b/rustar-aligner"

cpu_idle() { # percent idle, sampled over one second
  # macOS and Linux disagree on how to ask. `top -l` is macOS-only; on Linux it
  # exits non-zero, and under `set -o pipefail` that aborts the whole script
  # with no message, so this must branch rather than rely on a fallback.
  if [ "$(uname -s)" = "Darwin" ]; then
    top -l 2 -n 0 2>/dev/null | awk '/CPU usage/ {gsub("%","",$(NF-1)); v=$(NF-1)} END {print v+0}'
    return
  fi
  # Linux: two /proc/stat samples a second apart. Field 5 is idle; the total is
  # every field, so the ratio is idle time over elapsed jiffies across all CPUs.
  local s1 s2
  s1=$(awk '/^cpu /{t=0; for (j=2; j<=NF; j++) t+=$j; print $5, t; exit}' /proc/stat)
  sleep 1
  s2=$(awk '/^cpu /{t=0; for (j=2; j<=NF; j++) t+=$j; print $5, t; exit}' /proc/stat)
  awk -v a="$s1" -v b="$s2" 'BEGIN {
    split(a, x, " "); split(b, y, " ");
    di = y[1] - x[1]; dt = y[2] - x[2];
    print (dt > 0) ? int(100 * di / dt + 0.5) : 100
  }'
}

# Returns "<seconds> <idle_before> <idle_after>".
run() { # $1 = side directory
  local dir="$WORK/$1"
  local before after t0 t1
  before=$(cpu_idle)
  rm -rf "$dir/out"
  mkdir -p "$dir/out"
  t0=$(python3 -c 'import time; print(time.time())')
  # shellcheck disable=SC2086 # MODE and EXTRA are deliberately word-split
  (cd "$dir" && ./rustar-aligner --genomeDir "$GENOME_DIR" --readFilesIn "$READS" \
     --runThreadN "$THREADS" --outSAMtype $MODE $EXTRA \
     --outFileNamePrefix ./out/ >/dev/null 2>&1)
  t1=$(python3 -c "import time; print(f'{time.time()-$t0:.2f}')")
  after=$(cpu_idle)
  echo "$t1 $before $after"
}

echo "mode='$MODE' threads=$THREADS rounds=$ROUNDS min_idle=${MIN_IDLE}%" >&2
run a >/dev/null   # discarded: warms the page cache for the reads and the index

kept_a=()
kept_b=()
dropped=0

for i in $(seq 1 "$ROUNDS"); do
  # Flip the order on even rounds so a machine that drifts in one direction
  # cannot favour whichever side runs first.
  if (( i % 2 )); then
    read -r ta ia1 ia2 <<<"$(run a)"
    read -r tb ib1 ib2 <<<"$(run b)"
  else
    read -r tb ib1 ib2 <<<"$(run b)"
    read -r ta ia1 ia2 <<<"$(run a)"
  fi

  worst=$(printf '%s\n%s\n%s\n%s\n' "$ia1" "$ia2" "$ib1" "$ib2" | sort -n | head -1)
  if awk -v i="$worst" -v m="$MIN_IDLE" 'BEGIN{exit !(i < m)}'; then
    printf 'round%-2s %s=%ss %s=%ss  DROPPED (idle fell to %s%%)\n' \
      "$i" "$NAME_A" "$ta" "$NAME_B" "$tb" "$worst"
    dropped=$((dropped + 1))
    continue
  fi
  printf 'round%-2s %s=%ss %s=%ss  idle>=%s%%\n' "$i" "$NAME_A" "$ta" "$NAME_B" "$tb" "$worst"
  kept_a+=("$ta")
  kept_b+=("$tb")
done

median() { printf '%s\n' "$@" | sort -n | awk '{v[NR]=$1} END {print (NR%2) ? v[(NR+1)/2] : (v[NR/2]+v[NR/2+1])/2}'; }

echo
if [ ${#kept_a[@]} -eq 0 ]; then
  echo "every round was dropped: the machine was never quiet enough to measure." >&2
  echo "wait, or lower BENCH_MIN_IDLE and say so when reporting the numbers." >&2
  exit 1
fi

ma=$(median "${kept_a[@]}")
mb=$(median "${kept_b[@]}")
printf 'kept %d of %d rounds (%d dropped)\n' "${#kept_a[@]}" "$ROUNDS" "$dropped"
printf 'median  %s=%ss  %s=%ss  delta=%s%%\n' \
  "$NAME_A" "$ma" "$NAME_B" "$mb" \
  "$(python3 -c "print(f'{($mb-$ma)/$ma*100:+.1f}')")"

# A difference smaller than the spread within either side is not a result.
sa=$(python3 -c "v=sorted([$(IFS=,; echo "${kept_a[*]}")]); print(f'{v[-1]-v[0]:.2f}')")
sb=$(python3 -c "v=sorted([$(IFS=,; echo "${kept_b[*]}")]); print(f'{v[-1]-v[0]:.2f}')")
printf 'spread  %s=%ss  %s=%ss\n' "$NAME_A" "$sa" "$NAME_B" "$sb"
python3 - "$ma" "$mb" "$sa" "$sb" <<'PY'
import sys
ma, mb, sa, sb = (float(x) for x in sys.argv[1:5])
if abs(mb - ma) < max(sa, sb):
    print("\nthe difference between the medians is smaller than the spread within a side.")
    print("that is not a result: report it as unmeasured rather than as a speedup.")
PY
