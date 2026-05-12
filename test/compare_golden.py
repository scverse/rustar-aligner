#!/usr/bin/env python3

"""
compare_golden.py - Compare current test output against golden reference

Used by CI to detect regressions.

Usage:
    python compare_golden.py --golden golden/yeast_100/stats.json --current results/.../rustar_aligner --tolerance 0.01
"""

import argparse
import json
import sys
import os


def load_golden_stats(golden_file):
    """Load golden statistics from JSON file."""
    with open(golden_file, 'r') as f:
        return json.load(f)


def extract_current_stats(output_dir):
    """Extract statistics from current test output."""
    # Find SAM or BAM file
    sam_file = None
    for filename in ["Aligned.out.sam", "Aligned.out.bam"]:
        path = os.path.join(output_dir, filename)
        if os.path.exists(path):
            sam_file = path
            break

    if not sam_file:
        raise FileNotFoundError(f"No alignment file found in {output_dir}")

    # Extract stats (same logic as save_golden.sh)
    stats = {
        "total_reads": 0,
        "mapped_reads": 0,
        "unmapped_reads": 0,
        "unique_mapped": 0,
        "multi_mapped": 0,
        "per_chromosome": {},
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

            stats["total_reads"] += 1

            if (flag & 4) != 0:  # Unmapped
                stats["unmapped_reads"] += 1
            else:
                stats["mapped_reads"] += 1
                stats["per_chromosome"][rname] = stats["per_chromosome"].get(rname, 0) + 1

                if mapq == 255:
                    stats["unique_mapped"] += 1
                else:
                    stats["multi_mapped"] += 1

    return stats


def compare_stats(golden, current, tolerance):
    """Compare golden and current statistics."""
    messages = []
    passed = True

    def check_metric(name, golden_val, current_val):
        """Check if metric is within tolerance."""
        nonlocal passed

        if golden_val == 0 and current_val == 0:
            messages.append(f"  {name:25s} {current_val:>8} vs {golden_val:>8}  ✓")
            return

        if golden_val == 0:
            diff_pct = float('inf') if current_val != 0 else 0
        else:
            diff_pct = abs(current_val - golden_val) / golden_val

        status = "✓" if diff_pct <= tolerance else "✗ REGRESSION"
        diff_str = f"({current_val - golden_val:+.0f})" if abs(current_val - golden_val) > 0 else ""

        messages.append(f"  {name:25s} {current_val:>8.0f} vs {golden_val:>8.0f}  {status} {diff_str}")

        if diff_pct > tolerance:
            passed = False

    messages.append("=== Golden Output Comparison ===")

    check_metric("Total reads", golden["total_reads"], current["total_reads"])
    check_metric("Mapped reads", golden["mapped_reads"], current["mapped_reads"])
    check_metric("Unmapped reads", golden["unmapped_reads"], current["unmapped_reads"])
    check_metric("Unique mapped", golden["unique_mapped"], current["unique_mapped"])
    check_metric("Multi-mapped", golden["multi_mapped"], current["multi_mapped"])

    # Per-chromosome comparison
    all_chroms = set(golden.get("per_chromosome", {}).keys()) | set(current.get("per_chromosome", {}).keys())

    if all_chroms:
        messages.append("\n  Per-chromosome:")
        for chrom in sorted(all_chroms):
            golden_count = golden.get("per_chromosome", {}).get(chrom, 0)
            current_count = current.get("per_chromosome", {}).get(chrom, 0)

            if golden_count > 0:
                diff_pct = abs(current_count - golden_count) / golden_count
                if diff_pct > tolerance:
                    messages.append(f"    {chrom}: {current_count} vs {golden_count} ✗")
                    passed = False

    return passed, messages


def main():
    parser = argparse.ArgumentParser(description='Compare current output against golden reference')
    parser.add_argument('--golden', required=True, help='Golden stats.json file')
    parser.add_argument('--current', required=True, help='Current output directory')
    parser.add_argument('--tolerance', type=float, default=0.01, help='Tolerance for differences (default: 0.01 = 1%%)')

    args = parser.parse_args()

    # Load golden stats
    print(f"Loading golden reference: {args.golden}")
    golden_stats = load_golden_stats(args.golden)

    # Extract current stats
    print(f"Extracting current statistics: {args.current}")
    current_stats = extract_current_stats(args.current)

    # Compare
    passed, messages = compare_stats(golden_stats, current_stats, args.tolerance)

    # Print results
    for msg in messages:
        print(msg)

    print("\n" + "=" * 50)
    if passed:
        print("Status: PASS ✓")
        print("Output matches golden reference")
    else:
        print("Status: FAIL ✗")
        print("REGRESSION DETECTED - output differs from golden reference")

    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
