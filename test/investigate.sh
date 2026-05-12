#!/usr/bin/env bash

# investigate.sh - Investigation helper for alignment discrepancies
# Usage: ./investigate.sh <comparison_output_file>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEBUG_DIR="$SCRIPT_DIR/debug"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <comparison_output_file>"
    echo ""
    echo "Example:"
    echo "  $0 test/results/20260207_120000_yeast_1k/comparison/alignment_diff.txt"
    exit 1
fi

COMPARISON_FILE="$1"

if [[ ! -f "$COMPARISON_FILE" ]]; then
    echo "ERROR: Comparison file not found: $COMPARISON_FILE"
    exit 1
fi

log() {
    echo "[investigate] $*"
}

mkdir -p "$DEBUG_DIR"

log "=========================================="
log "Discrepancy Investigation"
log "=========================================="
log "Comparison file: $COMPARISON_FILE"
log ""

# ==============================================================================
# Parse discrepancies from comparison file
# ==============================================================================

log "Parsing discrepancies..."

# Extract discrepant read names and issues
python3 - <<'EOF' "$COMPARISON_FILE" "$DEBUG_DIR"
import sys
import re

comparison_file = sys.argv[1]
debug_dir = sys.argv[2]

discrepancies = []

with open(comparison_file, 'r') as f:
    in_discrepancies = False
    for line in f:
        if 'discrepant reads:' in line:
            in_discrepancies = True
            continue

        if in_discrepancies:
            # Look for lines like: "  read_137: Different chromosome (rustar-aligner=chr1, STAR=chr2)"
            match = re.match(r'\s+(\S+):\s+(.+)', line)
            if match:
                read_name = match.group(1)
                issue = match.group(2)
                discrepancies.append((read_name, issue))
            elif line.strip().startswith('Summary:'):
                break

# Save discrepant read names
if discrepancies:
    with open(f"{debug_dir}/discrepant_reads.txt", 'w') as f:
        for read_name, issue in discrepancies:
            f.write(f"{read_name}\t{issue}\n")

    print(f"Found {len(discrepancies)} discrepant reads")
    for read_name, issue in discrepancies[:10]:
        print(f"  {read_name}: {issue}")

    if len(discrepancies) > 10:
        print(f"  ... and {len(discrepancies) - 10} more")
else:
    print("No discrepancies found in comparison file")

EOF

if [[ ! -f "$DEBUG_DIR/discrepant_reads.txt" ]]; then
    log "No discrepancies to investigate"
    exit 0
fi

# ==============================================================================
# Extract discrepant reads from FASTQ
# ==============================================================================

log ""
log "Extracting discrepant reads from FASTQ..."

# Find the test directory
TEST_DIR=$(dirname "$(dirname "$COMPARISON_FILE")")
TEST_NAME=$(basename "$TEST_DIR" | sed 's/^[0-9_]*//')

log "Test directory: $TEST_DIR"
log "Test name: $TEST_NAME"

# Find FASTQ file
DATA_DIR="$SCRIPT_DIR/data/small/yeast"
READS_DIR="$DATA_DIR/reads"
FASTQ_FILE=""

for f in "$READS_DIR"/*.fastq.gz; do
    if [[ -f "$f" ]]; then
        FASTQ_FILE="$f"
        break
    fi
done

if [[ -z "$FASTQ_FILE" ]]; then
    log "WARNING: No FASTQ file found for extraction"
else
    log "Extracting from: $FASTQ_FILE"

    # Extract reads
    python3 - <<'EOF' "$FASTQ_FILE" "$DEBUG_DIR/discrepant_reads.txt" "$DEBUG_DIR/discrepant_reads.fq"
import sys
import gzip

fastq_file = sys.argv[1]
discrepant_file = sys.argv[2]
output_file = sys.argv[3]

# Load discrepant read names
discrepant_names = set()
with open(discrepant_file, 'r') as f:
    for line in f:
        read_name = line.split('\t')[0]
        discrepant_names.add(read_name)

# Extract reads from FASTQ
extracted = 0
with gzip.open(fastq_file, 'rt') as fin, open(output_file, 'w') as fout:
    while True:
        # Read FASTQ record (4 lines)
        header = fin.readline()
        if not header:
            break
        seq = fin.readline()
        plus = fin.readline()
        qual = fin.readline()

        # Extract read name (remove @ and everything after first space)
        read_name = header[1:].split()[0]

        if read_name in discrepant_names:
            fout.write(header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)
            extracted += 1

print(f"Extracted {extracted} discrepant reads to {output_file}")

EOF

    log "Extracted discrepant reads to: $DEBUG_DIR/discrepant_reads.fq"
fi

# ==============================================================================
# Suggest investigation steps
# ==============================================================================

log ""
log "=========================================="
log "Investigation Suggestions"
log "=========================================="

# Categorize discrepancies
python3 - <<'EOF' "$DEBUG_DIR/discrepant_reads.txt"
import sys
from collections import defaultdict

discrepant_file = sys.argv[1]

categories = defaultdict(list)

with open(discrepant_file, 'r') as f:
    for line in f:
        read_name, issue = line.strip().split('\t', 1)

        if 'chromosome' in issue.lower():
            categories['different_chromosome'].append((read_name, issue))
        elif 'mapping status' in issue.lower():
            categories['mapping_status'].append((read_name, issue))
        elif 'position' in issue.lower():
            categories['position_diff'].append((read_name, issue))
        elif 'mapq' in issue.lower():
            categories['mapq_diff'].append((read_name, issue))
        elif 'cigar' in issue.lower():
            categories['cigar_diff'].append((read_name, issue))
        elif 'strand' in issue.lower():
            categories['strand_diff'].append((read_name, issue))
        else:
            categories['other'].append((read_name, issue))

# Print suggestions for each category
if categories['different_chromosome']:
    print(f"\n{len(categories['different_chromosome'])} reads mapped to different chromosomes")
    print("  Likely cause: Seed clustering or multi-mapping selection difference")
    print("  Check rustar-aligner: src/align/stitch.rs (cluster_seeds, select_best_alignment)")
    print("  Check STAR:   source/stitchWindowAligns.cpp (stitchWindowSeeds)")
    print("  Suggested action: Compare seed positions for these reads")

if categories['mapping_status']:
    print(f"\n{len(categories['mapping_status'])} reads have different mapping status")
    print("  Likely cause: Filtering threshold difference")
    print("  Check rustar-aligner: src/align/read_align.rs (transcript filtering)")
    print("  Check STAR:   source/ReadAlign_outputAlignments.cpp (outFilterScoreMin)")
    print("  Suggested action: Check outFilterScoreMin, outFilterMatchNmin parameters")

if categories['position_diff']:
    print(f"\n{len(categories['position_diff'])} reads have different positions")
    print("  Likely cause: Tie-breaking in multi-mapping reads")
    print("  Check rustar-aligner: src/align/stitch.rs (best alignment selection)")
    print("  Check STAR:   source/ReadAlign_outputAlignments.cpp")
    print("  Suggested action: Check if MAPQ < 10 (expected for ties)")

if categories['cigar_diff']:
    print(f"\n{len(categories['cigar_diff'])} reads have different CIGAR strings")
    print("  Likely cause: DP stitching or gap penalty difference")
    print("  Check rustar-aligner: src/align/stitch.rs (stitch_seeds_dp)")
    print("  Check STAR:   source/stitchWindowAligns.cpp")
    print("  Suggested action: Re-run with debug logging (RUST_LOG=debug)")

if categories['strand_diff']:
    print(f"\n{len(categories['strand_diff'])} reads have different strand assignments")
    print("  Likely cause: Reverse complement handling")
    print("  Check rustar-aligner: src/align/read_align.rs (reverse complement logic)")
    print("  Check STAR:   source/ReadAlign.cpp")

if categories['mapq_diff']:
    print(f"\n{len(categories['mapq_diff'])} reads have different MAPQ values")
    print("  Likely cause: Different multi-mapping count or MAPQ calculation")
    print("  Check rustar-aligner: src/mapq.rs")
    print("  Check STAR:   source/ReadAlign_outputAlignments.cpp")
    print("  Note: Small MAPQ differences for multi-mappers are expected")

print("\n" + "=" * 50)
print("STAR Source Code Reference")
print("=" * 50)
print("Repository: https://github.com/alexdobin/STAR")
print("\nKey files:")
print("  source/ReadAlign.cpp               - Main alignment driver")
print("  source/ReadAlign_outputAlignments.cpp - Output filtering")
print("  source/stitchWindowAligns.cpp      - Seed stitching")
print("  source/alignSmithWaterman.cpp      - Local alignment")
print("  source/Genome.cpp                  - Genome loading")
print("  source/SuffixArrayFuns.cpp         - Suffix array search")

EOF

# ==============================================================================
# Re-run rustar-aligner with debug logging
# ==============================================================================

log ""
log "=========================================="
log "Debug Re-run"
log "=========================================="

if [[ -f "$DEBUG_DIR/discrepant_reads.fq" ]]; then
    log "To re-run rustar-aligner with debug logging on discrepant reads:"
    log ""
    log "  RUST_LOG=debug $PROJECT_ROOT/target/release/rustar-aligner \\"
    log "    --runMode alignReads \\"
    log "    --genomeDir $DATA_DIR/genome_index \\"
    log "    --readFilesIn $DEBUG_DIR/discrepant_reads.fq \\"
    log "    --outFileNamePrefix $DEBUG_DIR/ \\"
    log "    --outSAMtype SAM"
    log ""
    log "This will generate detailed logs for seed finding, clustering, and stitching."
fi

log ""
log "Investigation complete. Results in: $DEBUG_DIR"
