#!/usr/bin/env python3

"""
compare_sam.py - Compare SAM/BAM files from rustar-aligner and STAR

Performs multi-level comparison:
1. Statistics comparison (read counts, mapping rates, MAPQ distribution)
2. Record-level comparison (detailed per-read analysis)
3. CIGAR comparison (alignment structure validation)

Usage:
    python compare_sam.py --star STAR.sam --rustar-aligner rustar-aligner.sam [options]
"""

import argparse
import sys
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import re


@dataclass
class AlignmentStats:
    """Statistics for a SAM file."""
    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    unique_mapped: int = 0  # MAPQ == 255
    multi_mapped: int = 0   # MAPQ < 255
    per_chromosome: Dict[str, int] = None
    mapq_distribution: Dict[int, int] = None
    mean_mapq: float = 0.0
    mean_alignment_length: float = 0.0
    total_junctions: int = 0

    def __post_init__(self):
        if self.per_chromosome is None:
            self.per_chromosome = defaultdict(int)
        if self.mapq_distribution is None:
            self.mapq_distribution = defaultdict(int)


@dataclass
class SamRecord:
    """Simplified SAM record for comparison."""
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str

    def is_unmapped(self) -> bool:
        return (self.flag & 4) != 0

    def is_reverse(self) -> bool:
        return (self.flag & 16) != 0

    def is_paired(self) -> bool:
        return (self.flag & 1) != 0


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse CIGAR string into list of (length, operation) tuples."""
    if cigar == "*":
        return []
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    return [(int(length), op) for length, op in operations]


def cigar_aligned_length(cigar: str) -> int:
    """Calculate total aligned length from CIGAR (M, D, N, =, X operations)."""
    ops = parse_cigar(cigar)
    return sum(length for length, op in ops if op in 'MDN=X')


def cigar_junctions(cigar: str) -> int:
    """Count splice junctions (N operations) in CIGAR."""
    ops = parse_cigar(cigar)
    return sum(1 for _, op in ops if op == 'N')


def parse_sam_file(filename: str) -> Tuple[Dict[str, SamRecord], AlignmentStats]:
    """Parse SAM file and extract records and statistics."""
    records = {}
    stats = AlignmentStats()

    mapq_sum = 0
    alignment_length_sum = 0
    mapped_count = 0

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            rnext = fields[6]
            pnext = int(fields[7])
            tlen = int(fields[8])
            seq = fields[9]

            record = SamRecord(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq)
            records[qname] = record

            stats.total_reads += 1

            if record.is_unmapped():
                stats.unmapped_reads += 1
            else:
                stats.mapped_reads += 1
                stats.per_chromosome[rname] += 1
                stats.mapq_distribution[mapq] += 1
                mapq_sum += mapq
                mapped_count += 1

                if mapq == 255:
                    stats.unique_mapped += 1
                else:
                    stats.multi_mapped += 1

                aligned_length = cigar_aligned_length(cigar)
                alignment_length_sum += aligned_length

                stats.total_junctions += cigar_junctions(cigar)

    if mapped_count > 0:
        stats.mean_mapq = mapq_sum / mapped_count
        stats.mean_alignment_length = alignment_length_sum / mapped_count

    return records, stats


def compare_statistics(star_stats: AlignmentStats, rustar_aligner_stats: AlignmentStats, tolerance: float) -> Tuple[bool, List[str]]:
    """Compare statistics between STAR and rustar-aligner."""
    messages = []
    passed = True

    def check_metric(name: str, star_val: float, rustar_aligner_val: float, allow_diff: bool = False) -> bool:
        """Check if metric is within tolerance."""
        if star_val == 0 and rustar_aligner_val == 0:
            messages.append(f"  {name:25s} {rustar_aligner_val:>8} vs {star_val:>8}  ✓")
            return True

        if star_val == 0:
            diff_pct = float('inf') if rustar_aligner_val != 0 else 0
        else:
            diff_pct = abs(rustar_aligner_val - star_val) / star_val

        status = "✓" if diff_pct <= tolerance else "⚠"
        diff_str = f"({abs(rustar_aligner_val - star_val):+.0f})" if diff_pct > 0.001 else ""

        messages.append(f"  {name:25s} {rustar_aligner_val:>8.0f} vs {star_val:>8.0f}  {status} {diff_str}")

        if diff_pct > tolerance and not allow_diff:
            return False
        return True

    messages.append("\n=== Statistics Comparison: rustar-aligner vs STAR ===")

    check_metric("Total reads", star_stats.total_reads, rustar_aligner_stats.total_reads)
    passed &= check_metric("Mapped reads", star_stats.mapped_reads, rustar_aligner_stats.mapped_reads, allow_diff=True)
    passed &= check_metric("Unmapped reads", star_stats.unmapped_reads, rustar_aligner_stats.unmapped_reads, allow_diff=True)
    passed &= check_metric("Unique mapped", star_stats.unique_mapped, rustar_aligner_stats.unique_mapped, allow_diff=True)
    passed &= check_metric("Multi-mapped", star_stats.multi_mapped, rustar_aligner_stats.multi_mapped, allow_diff=True)
    check_metric("Mean MAPQ", star_stats.mean_mapq, rustar_aligner_stats.mean_mapq, allow_diff=True)
    check_metric("Mean alignment length", star_stats.mean_alignment_length, rustar_aligner_stats.mean_alignment_length, allow_diff=True)
    check_metric("Total junctions", star_stats.total_junctions, rustar_aligner_stats.total_junctions, allow_diff=True)

    # Mapping rate
    if star_stats.total_reads > 0:
        star_rate = 100.0 * star_stats.mapped_reads / star_stats.total_reads
        rustar_aligner_rate = 100.0 * rustar_aligner_stats.mapped_reads / rustar_aligner_stats.total_reads
        messages.append(f"  {'Mapping rate':25s} {rustar_aligner_rate:>7.2f}% vs {star_rate:>7.2f}%")

    # Unique mapping rate
    if star_stats.mapped_reads > 0:
        star_unique_rate = 100.0 * star_stats.unique_mapped / star_stats.mapped_reads
        rustar_aligner_unique_rate = 100.0 * rustar_aligner_stats.unique_mapped / rustar_aligner_stats.mapped_reads
        messages.append(f"  {'Unique rate (of mapped)':25s} {rustar_aligner_unique_rate:>7.2f}% vs {star_unique_rate:>7.2f}%")

    return passed, messages


def compare_records(star_records: Dict[str, SamRecord], rustar_aligner_records: Dict[str, SamRecord], tolerance: float) -> Tuple[bool, List[str]]:
    """Compare individual alignment records."""
    messages = []
    messages.append("\n=== Record-Level Comparison ===")

    # Find common reads
    star_reads = set(star_records.keys())
    rustar_aligner_reads = set(rustar_aligner_records.keys())

    common_reads = star_reads & rustar_aligner_reads
    star_only = star_reads - rustar_aligner_reads
    rustar_aligner_only = rustar_aligner_reads - star_reads

    messages.append(f"  Common reads: {len(common_reads)}")
    if star_only:
        messages.append(f"  Only in STAR: {len(star_only)}")
    if rustar_aligner_only:
        messages.append(f"  Only in rustar-aligner: {len(rustar_aligner_only)}")

    # Compare common reads
    discrepancies = []
    position_diffs = []
    mapq_diffs = []
    cigar_diffs = []
    flag_diffs = []

    for qname in common_reads:
        star_rec = star_records[qname]
        rustar_aligner_rec = rustar_aligner_records[qname]

        issues = []

        # Compare unmapped status
        if star_rec.is_unmapped() != rustar_aligner_rec.is_unmapped():
            issues.append(f"mapping status (rustar-aligner={'unmapped' if rustar_aligner_rec.is_unmapped() else 'mapped'}, STAR={'unmapped' if star_rec.is_unmapped() else 'mapped'})")

        # Compare mapped reads
        if not star_rec.is_unmapped() and not rustar_aligner_rec.is_unmapped():
            # Chromosome
            if star_rec.rname != rustar_aligner_rec.rname:
                issues.append(f"chromosome (rustar-aligner={rustar_aligner_rec.rname}, STAR={star_rec.rname})")

            # Position (allow small differences for ties)
            pos_diff = abs(star_rec.pos - rustar_aligner_rec.pos)
            if pos_diff > 10:
                issues.append(f"position (rustar-aligner={rustar_aligner_rec.pos}, STAR={star_rec.pos}, diff={pos_diff})")
                position_diffs.append(pos_diff)

            # MAPQ (allow differences for multi-mappers)
            mapq_diff = abs(star_rec.mapq - rustar_aligner_rec.mapq)
            if mapq_diff > 5 and not (star_rec.mapq < 10 and rustar_aligner_rec.mapq < 10):
                issues.append(f"MAPQ (rustar-aligner={rustar_aligner_rec.mapq}, STAR={star_rec.mapq})")
                mapq_diffs.append(mapq_diff)

            # CIGAR
            if star_rec.cigar != rustar_aligner_rec.cigar:
                # Parse and compare structure
                star_ops = parse_cigar(star_rec.cigar)
                rustar_aligner_ops = parse_cigar(rustar_aligner_rec.cigar)

                if star_ops != rustar_aligner_ops:
                    issues.append(f"CIGAR (rustar-aligner={rustar_aligner_rec.cigar}, STAR={star_rec.cigar})")
                    cigar_diffs.append((star_rec.cigar, rustar_aligner_rec.cigar))

            # Flags (check important bits)
            if star_rec.is_reverse() != rustar_aligner_rec.is_reverse():
                issues.append(f"strand (rustar-aligner={'-' if rustar_aligner_rec.is_reverse() else '+'}, STAR={'-' if star_rec.is_reverse() else '+'})")
                flag_diffs.append((star_rec.flag, rustar_aligner_rec.flag))

        if issues:
            discrepancies.append((qname, issues))

    # Report discrepancies
    if discrepancies:
        messages.append(f"\n  Found {len(discrepancies)} discrepant reads:")

        # Show first 10 discrepancies
        for qname, issues in discrepancies[:10]:
            messages.append(f"    {qname}: {', '.join(issues)}")

        if len(discrepancies) > 10:
            messages.append(f"    ... and {len(discrepancies) - 10} more")

    # Summary statistics
    total_compared = len(common_reads)
    matching = total_compared - len(discrepancies)
    match_rate = 100.0 * matching / total_compared if total_compared > 0 else 0

    messages.append(f"\n  Summary: {matching}/{total_compared} reads match ({match_rate:.1f}%)")

    if position_diffs:
        messages.append(f"  Position differences: mean={sum(position_diffs)/len(position_diffs):.1f}, max={max(position_diffs)}")

    if mapq_diffs:
        messages.append(f"  MAPQ differences: mean={sum(mapq_diffs)/len(mapq_diffs):.1f}, max={max(mapq_diffs)}")

    # Pass if match rate is within tolerance
    passed = match_rate >= (100.0 * (1.0 - tolerance))

    return passed, messages


def main():
    parser = argparse.ArgumentParser(description='Compare SAM files from rustar-aligner and STAR')
    parser.add_argument('--star', required=True, help='STAR SAM file')
    parser.add_argument('--rustar-aligner', required=True, help='rustar-aligner SAM file')
    parser.add_argument('--tolerance', type=float, default=0.01, help='Tolerance for differences (default: 0.01 = 1%%)')
    parser.add_argument('--output', help='Output file for comparison report')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Parse both SAM files
    print(f"Parsing STAR output: {args.star}")
    star_records, star_stats = parse_sam_file(args.star)

    print(f"Parsing rustar-aligner output: {args.rustar_aligner}")
    rustar_aligner_records, rustar_aligner_stats = parse_sam_file(args.rustar_aligner)

    # Compare statistics
    stats_passed, stats_messages = compare_statistics(star_stats, rustar_aligner_stats, args.tolerance)

    # Compare records
    records_passed, records_messages = compare_records(star_records, rustar_aligner_records, args.tolerance)

    # Combine messages
    all_messages = stats_messages + records_messages

    # Overall status
    overall_passed = stats_passed and records_passed

    all_messages.append("\n" + "=" * 50)
    if overall_passed:
        all_messages.append("Status: PASS ✓")
        all_messages.append(f"Alignment outputs match within {args.tolerance*100:.1f}% tolerance")
    else:
        all_messages.append("Status: FAIL ✗")
        all_messages.append("Significant differences detected")

    # Print and optionally save
    output_text = "\n".join(all_messages)
    print(output_text)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_text)
        print(f"\nComparison report saved to: {args.output}")

    # Exit with appropriate code
    sys.exit(0 if overall_passed else 1)


if __name__ == '__main__':
    main()
