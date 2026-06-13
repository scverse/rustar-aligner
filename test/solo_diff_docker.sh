#!/usr/bin/env bash
# Run the STARsolo CellRanger differential test (rustar-aligner vs real STAR) in
# a consistent Linux container, so the comparison works regardless of the host
# (the macOS STAR build has a FASTQ-read bug; Linux STAR works).
#
# Requires a Docker-compatible runtime. On macOS without Docker Desktop:
#   brew install colima docker && colima start
#
# Usage:  test/solo_diff_docker.sh [N_RUNS]
set -euo pipefail

cd "$(dirname "$0")/.."
RUNS="${1:-1}"
IMAGE=rustar-solodiff

docker build -f test/Dockerfile.solodiff -t "$IMAGE" . >/dev/null

# Build rustar for Linux into a host-mounted dir (persisted across runs), then
# run the harness against the Linux STAR + Linux rustar binary.
docker run --rm -v "$PWD":/work -w /work -e CARGO_TARGET_DIR=/work/target-linux "$IMAGE" bash -c '
  set -e
  cargo build --release 2>&1 | tail -1
  RUSTAR=/work/target-linux/release/rustar-aligner
  STARBIN=$(which STAR)
  for i in $(seq 1 '"$RUNS"'); do
    echo "===== differential run $i ====="
    python3 test/solo_cellranger_diff.py --star "$STARBIN" --rustar "$RUSTAR"
  done
'
