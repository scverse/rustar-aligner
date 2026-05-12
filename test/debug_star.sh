#!/usr/bin/env bash
# debug_star.sh - Run STAR with read-name tracing on specific reads
#
# Usage:
#   # Trace reads from a rustar-aligner-only or STAR-only PE comparison:
#   ./debug_star.sh pe <rustar-aligner_sam> <star_sam> [n_reads]
#
#   # Trace specific read names directly:
#   ./debug_star.sh reads <read1,read2,...>
#
# Output: debug_star_output/ directory with STAR debug log and SAM

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="/home/jamfer/Dropbox/Bioinformatics/tools/repos/STAR/source/STAR"
GENOME_DIR="$SCRIPT_DIR/data/small/yeast/star_genome"
READS1="$SCRIPT_DIR/data/small/yeast/reads/ERR12389696_sub_1_10k.fastq.gz"
READS2="$SCRIPT_DIR/data/small/yeast/reads/ERR12389696_sub_2_10k.fastq.gz"
OUT_DIR="$SCRIPT_DIR/debug_star_output"

mkdir -p "$OUT_DIR"

mode="${1:-pe}"

extract_reads_from_comparison() {
    local rustar_aligner_sam="$1"
    local star_sam="$2"
    local n="${3:-20}"
    python3 - "$rustar_aligner_sam" "$star_sam" "$n" <<'PYEOF'
import sys
from collections import defaultdict

rustar_aligner_sam, star_sam, n = sys.argv[1], sys.argv[2], int(sys.argv[3])

def parse_pe_sam(path):
    mapped = {}  # qname -> (chr, pos, mapq, cigar)
    with open(path) as f:
        for line in f:
            if line.startswith('@'): continue
            f2 = line.strip().split('\t')
            if len(f2) < 11: continue
            qname, flag = f2[0], int(f2[1])
            if flag & 0x100: continue   # skip secondary
            if flag & 0x4: continue     # skip unmapped
            if not (flag & 0x40): continue  # only mate1
            mapped[qname] = (f2[2], int(f2[3]), int(f2[4]), f2[5])
    return mapped

r = parse_pe_sam(rustar_aligner_sam)
s = parse_pe_sam(star_sam)

rustar_aligner_only = []
for qname in sorted(r.keys()):
    if qname not in s:
        rustar_aligner_only.append(qname)
    # same name but STAR has it unmapped would also count - check both dicts

# Also STAR-only
star_only = [q for q in sorted(s.keys()) if q not in r]

print(f"# rustar-aligner-only (false positives): {len(rustar_aligner_only)}")
print(f"# STAR-only (missed): {len(star_only)}")
print()

out_reads = rustar_aligner_only[:n] + star_only[:min(5, n)]
print(','.join(out_reads[:n]))
PYEOF
}

if [ "$mode" = "pe" ]; then
    if [ $# -lt 3 ]; then
        echo "Usage: $0 pe <rustar_aligner_sam> <star_sam> [n_reads=20]"
        exit 1
    fi
    rustar_aligner_sam="$2"
    star_sam="$3"
    n_reads="${4:-20}"

    echo "=== Extracting read names from PE comparison ==="
    read_list=$(extract_reads_from_comparison "$rustar_aligner_sam" "$star_sam" "$n_reads" | tail -1)
    echo "Reads to trace: $read_list"

elif [ "$mode" = "reads" ]; then
    if [ $# -lt 2 ]; then
        echo "Usage: $0 reads <read1,read2,...>"
        exit 1
    fi
    read_list="$2"
    echo "=== Tracing specified reads ==="
    echo "Reads: $read_list"
else
    echo "Unknown mode: $mode"
    echo "Usage: $0 pe <rustar_aligner_sam> <star_sam> [n_reads]"
    echo "       $0 reads <read1,read2,...>"
    exit 1
fi

echo ""
echo "=== Running STAR with debug tracing ==="
echo "Output dir: $OUT_DIR"
echo "Genome: $GENOME_DIR"
echo ""

# Run STAR with debug tracing, capturing stderr (debug log) separately
STAR_DEBUG_READS="$read_list" \
    "$STAR_BIN" \
    --runMode alignReads \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READS1" "$READS2" \
    --readFilesCommand zcat \
    --outSAMtype SAM \
    --outFileNamePrefix "$OUT_DIR/" \
    --runThreadN 1 \
    2>"$OUT_DIR/star_debug.log"

echo ""
echo "=== Debug log summary ==="
echo "Total debug lines: $(grep -c '^\[STAR_DBG' "$OUT_DIR/star_debug.log" 2>/dev/null || echo 0)"
echo ""

# Extract per-read summaries
python3 - "$OUT_DIR/star_debug.log" "$read_list" <<'PYEOF'
import sys
from collections import defaultdict

log_file = sys.argv[1]
read_list = sys.argv[2].split(',') if len(sys.argv) > 2 else []

events = defaultdict(list)

with open(log_file) as f:
    for line in f:
        if not line.startswith('[STAR_DBG:'):
            continue
        # Parse [STAR_DBG:readname] event: ...
        rest = line[len('[STAR_DBG:'):]
        bracket = rest.find(']')
        if bracket < 0: continue
        rname = rest[:bracket]
        msg = rest[bracket+2:].strip()
        events[rname].append(msg)

for rname in read_list:
    if rname not in events:
        print(f"\n{'='*60}")
        print(f"READ: {rname}")
        print("  (no debug events - not processed by STAR or not in this dataset)")
        continue

    print(f"\n{'='*60}")
    print(f"READ: {rname} ({len(events[rname])} events)")
    print(f"{'='*60}")

    # Group by event type
    for msg in events[rname]:
        tag = msg.split(':')[0] if ':' in msg else msg[:20]
        print(f"  {msg}")

PYEOF

echo ""
echo "Full debug log: $OUT_DIR/star_debug.log"
echo "STAR SAM:       $OUT_DIR/Aligned.out.sam"
