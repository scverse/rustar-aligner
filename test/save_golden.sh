#!/usr/bin/env bash

# save_golden.sh - Save test outputs as golden references
# Usage: ./save_golden.sh <test_name>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
GOLDEN_DIR="$SCRIPT_DIR/golden"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <test_name>"
    echo ""
    echo "Available test names:"
    echo "  yeast_100"
    echo "  yeast_1k"
    echo "  yeast_1k_paired"
    echo "  yeast_1k_twopass"
    echo "  yeast_10k"
    exit 1
fi

TEST_NAME="$1"

log() {
    echo "[save_golden] $*"
}

error() {
    echo "[save_golden] ERROR: $*" >&2
}

# Find most recent test result for this test
LATEST_RESULT=$(find "$RESULTS_DIR" -maxdepth 1 -name "*_${TEST_NAME}" -type d | sort | tail -1)

if [[ -z "$LATEST_RESULT" ]]; then
    error "No test results found for $TEST_NAME"
    error "Run: ./run_tests.sh $TEST_NAME"
    exit 1
fi

# Check if test passed
if [[ ! -f "$LATEST_RESULT/PASSED" ]]; then
    error "Test $TEST_NAME did not pass"
    error "Cannot save failing test as golden output"
    exit 1
fi

log "Found test result: $LATEST_RESULT"

# Create golden directory
GOLDEN_TEST_DIR="$GOLDEN_DIR/$TEST_NAME"
mkdir -p "$GOLDEN_TEST_DIR"

# Extract statistics from rustar-aligner output
RUSTAR_ALIGNER_DIR="$LATEST_RESULT/rustar-aligner"

if [[ ! -d "$RUSTAR_ALIGNER_DIR" ]]; then
    error "rustar-aligner output directory not found: $RUSTAR_ALIGNER_DIR"
    exit 1
fi

log "Extracting statistics..."

# Use Python to extract stats from SAM file
python3 - <<'EOF' "$RUSTAR_ALIGNER_DIR" "$GOLDEN_TEST_DIR"
import sys
import json
import os
from collections import defaultdict

def extract_sam_stats(sam_file):
    """Extract statistics from SAM file."""
    stats = {
        "total_reads": 0,
        "mapped_reads": 0,
        "unmapped_reads": 0,
        "unique_mapped": 0,
        "multi_mapped": 0,
        "per_chromosome": {},
        "mapq_sum": 0,
        "alignment_length_sum": 0,
    }

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            rname = fields[2]
            mapq = int(fields[4])
            cigar = fields[5]

            stats["total_reads"] += 1

            if (flag & 4) != 0:  # Unmapped
                stats["unmapped_reads"] += 1
            else:
                stats["mapped_reads"] += 1
                stats["per_chromosome"][rname] = stats["per_chromosome"].get(rname, 0) + 1
                stats["mapq_sum"] += mapq

                if mapq == 255:
                    stats["unique_mapped"] += 1
                else:
                    stats["multi_mapped"] += 1

                # Calculate alignment length from CIGAR
                import re
                ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                aligned_length = sum(int(length) for length, op in ops if op in 'MDN=X')
                stats["alignment_length_sum"] += aligned_length

    # Calculate means
    if stats["mapped_reads"] > 0:
        stats["mean_mapq"] = stats["mapq_sum"] / stats["mapped_reads"]
        stats["mean_alignment_length"] = stats["alignment_length_sum"] / stats["mapped_reads"]
    else:
        stats["mean_mapq"] = 0.0
        stats["mean_alignment_length"] = 0.0

    # Remove intermediate sums
    del stats["mapq_sum"]
    del stats["alignment_length_sum"]

    return stats

def extract_junction_stats(sj_file):
    """Extract junction statistics from SJ.out.tab."""
    if not os.path.exists(sj_file):
        return {"total_junctions": 0, "junctions": []}

    junctions = []

    with open(sj_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            junction = {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "strand": int(fields[3]),
                "motif": int(fields[4]),
                "annotated": int(fields[5]),
                "unique_reads": int(fields[6]),
                "multi_reads": int(fields[7]),
            }

            junctions.append(junction)

    return {
        "total_junctions": len(junctions),
        "junctions": junctions
    }

# Main
rustar_aligner_dir = sys.argv[1]
golden_dir = sys.argv[2]

# Find SAM or BAM file
sam_file = None
for filename in ["Aligned.out.sam", "Aligned.out.bam"]:
    path = os.path.join(rustar_aligner_dir, filename)
    if os.path.exists(path):
        sam_file = path
        break

if not sam_file:
    print(f"ERROR: No alignment file found in {rustar_aligner_dir}", file=sys.stderr)
    sys.exit(1)

# Extract stats
print(f"Extracting from: {sam_file}")
stats = extract_sam_stats(sam_file)

# Extract junction stats
sj_file = os.path.join(rustar_aligner_dir, "SJ.out.tab")
junction_stats = extract_junction_stats(sj_file)

stats["junctions_detected"] = junction_stats["total_junctions"]

# Add metadata
import datetime
metadata = {
    "test_name": os.path.basename(golden_dir),
    "timestamp": datetime.datetime.now().isoformat(),
    "source": rustar_aligner_dir,
}

# Save to JSON
stats_file = os.path.join(golden_dir, "stats.json")
with open(stats_file, 'w') as f:
    json.dump(stats, f, indent=2)

print(f"Saved statistics to: {stats_file}")

# Save junctions
if junction_stats["total_junctions"] > 0:
    junctions_file = os.path.join(golden_dir, "junctions.json")
    with open(junctions_file, 'w') as f:
        json.dump(junction_stats, f, indent=2)
    print(f"Saved junctions to: {junctions_file}")

# Save metadata
metadata_file = os.path.join(golden_dir, "metadata.json")
with open(metadata_file, 'w') as f:
    json.dump(metadata, f, indent=2)

print(f"Saved metadata to: {metadata_file}")

EOF

log "Golden outputs saved to: $GOLDEN_TEST_DIR"
log ""
log "Files created:"
ls -lh "$GOLDEN_TEST_DIR"

log ""
log "Golden outputs ready for version control"
log "Consider running: git add $GOLDEN_TEST_DIR"
