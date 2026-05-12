#!/usr/bin/env python3

"""
compare_junctions.py - Compare SJ.out.tab files from rustar-aligner and STAR

SJ.out.tab format (9 columns):
1. chromosome
2. first base of intron (1-based)
3. last base of intron (1-based)
4. strand (0=undefined, 1=+, 2=-)
5. intron motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=N-N)
6. 0=unannotated, 1=annotated
7. number of uniquely mapping reads
8. number of multi-mapping reads
9. maximum overhang

Usage:
    python compare_junctions.py --star STAR_SJ.out.tab --rustar-aligner rustar-aligner_SJ.out.tab [options]
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set


@dataclass
class Junction:
    """Splice junction record."""
    chrom: str
    start: int
    end: int
    strand: int
    motif: int
    annotated: int
    unique_reads: int
    multi_reads: int
    max_overhang: int

    def key(self) -> Tuple[str, int, int, int]:
        """Return junction key for comparison."""
        return (self.chrom, self.start, self.end, self.strand)

    def length(self) -> int:
        """Return intron length."""
        return self.end - self.start + 1

    def total_reads(self) -> int:
        """Return total supporting reads."""
        return self.unique_reads + self.multi_reads


MOTIF_NAMES = {
    0: "non-canonical",
    1: "GT/AG",
    2: "CT/AC",
    3: "GC/AG",
    4: "CT/GC",
    5: "AT/AC",
    6: "N-N",
}

STRAND_NAMES = {
    0: "undefined",
    1: "+",
    2: "-",
}


def parse_sj_file(filename: str) -> Dict[Tuple[str, int, int, int], Junction]:
    """Parse SJ.out.tab file."""
    junctions = {}

    with open(filename, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = int(fields[3])
            motif = int(fields[4])
            annotated = int(fields[5])
            unique_reads = int(fields[6])
            multi_reads = int(fields[7])
            max_overhang = int(fields[8])

            junction = Junction(
                chrom, start, end, strand, motif, annotated,
                unique_reads, multi_reads, max_overhang
            )

            junctions[junction.key()] = junction

    return junctions


def compare_junctions(star_junctions: Dict, rustar_aligner_junctions: Dict, tolerance: float) -> Tuple[bool, List[str]]:
    """Compare junction sets."""
    messages = []
    messages.append("\n=== Junction Comparison: rustar-aligner vs STAR ===")

    # Find common and unique junctions
    star_keys = set(star_junctions.keys())
    rustar_aligner_keys = set(rustar_aligner_junctions.keys())

    common_keys = star_keys & rustar_aligner_keys
    star_only = star_keys - rustar_aligner_keys
    rustar_aligner_only = rustar_aligner_keys - star_keys

    # Overall statistics
    messages.append(f"  Total junctions (rustar-aligner): {len(rustar_aligner_junctions)}")
    messages.append(f"  Total junctions (STAR):   {len(star_junctions)}")
    messages.append(f"  Common junctions:         {len(common_keys)}")

    if star_only:
        messages.append(f"  Only in STAR:             {len(star_only)}")
    if rustar_aligner_only:
        messages.append(f"  Only in rustar-aligner:           {len(rustar_aligner_only)}")

    # Calculate overlap rate
    total_junctions = len(star_keys | rustar_aligner_keys)
    if total_junctions > 0:
        overlap_rate = 100.0 * len(common_keys) / total_junctions
        messages.append(f"  Overlap rate:             {overlap_rate:.1f}%")

    # Show unique junctions (up to 5 each)
    if rustar_aligner_only:
        messages.append(f"\n  Unique to rustar-aligner ({len(rustar_aligner_only)} total, showing first 5):")
        for key in sorted(rustar_aligner_only)[:5]:
            junc = rustar_aligner_junctions[key]
            motif_name = MOTIF_NAMES.get(junc.motif, "unknown")
            strand_name = STRAND_NAMES.get(junc.strand, "?")
            messages.append(
                f"    {junc.chrom}:{junc.start}-{junc.end} "
                f"({strand_name}, {motif_name}, "
                f"{junc.total_reads()} reads)"
            )

    if star_only:
        messages.append(f"\n  Unique to STAR ({len(star_only)} total, showing first 5):")
        for key in sorted(star_only)[:5]:
            junc = star_junctions[key]
            motif_name = MOTIF_NAMES.get(junc.motif, "unknown")
            strand_name = STRAND_NAMES.get(junc.strand, "?")
            messages.append(
                f"    {junc.chrom}:{junc.start}-{junc.end} "
                f"({strand_name}, {motif_name}, "
                f"{junc.total_reads()} reads)"
            )

    # Compare common junctions
    motif_mismatches = []
    annotated_mismatches = []
    count_differences = []

    for key in common_keys:
        star_junc = star_junctions[key]
        rustar_aligner_junc = rustar_aligner_junctions[key]

        # Check motif classification
        if star_junc.motif != rustar_aligner_junc.motif:
            motif_mismatches.append((key, star_junc.motif, rustar_aligner_junc.motif))

        # Check annotation status
        if star_junc.annotated != rustar_aligner_junc.annotated:
            annotated_mismatches.append((key, star_junc.annotated, rustar_aligner_junc.annotated))

        # Check read counts (allow some tolerance)
        star_total = star_junc.total_reads()
        rustar_aligner_total = rustar_aligner_junc.total_reads()

        if star_total > 0:
            diff_pct = abs(rustar_aligner_total - star_total) / star_total
            if diff_pct > tolerance:
                count_differences.append((key, star_total, rustar_aligner_total, diff_pct))

    # Report mismatches
    if motif_mismatches:
        messages.append(f"\n  Motif classification differences: {len(motif_mismatches)}")
        for key, star_motif, rustar_aligner_motif in motif_mismatches[:5]:
            chrom, start, end, strand = key
            messages.append(
                f"    {chrom}:{start}-{end}: "
                f"rustar-aligner={MOTIF_NAMES.get(rustar_aligner_motif, '?')}, "
                f"STAR={MOTIF_NAMES.get(star_motif, '?')}"
            )

    if annotated_mismatches:
        messages.append(f"\n  Annotation status differences: {len(annotated_mismatches)}")
        for key, star_ann, rustar_aligner_ann in annotated_mismatches[:5]:
            chrom, start, end, strand = key
            messages.append(
                f"    {chrom}:{start}-{end}: "
                f"rustar-aligner={'annotated' if rustar_aligner_ann else 'novel'}, "
                f"STAR={'annotated' if star_ann else 'novel'}"
            )

    if count_differences:
        messages.append(f"\n  Read count differences (>{tolerance*100:.0f}% threshold): {len(count_differences)}")
        for key, star_count, rustar_aligner_count, diff_pct in sorted(count_differences, key=lambda x: x[3], reverse=True)[:5]:
            chrom, start, end, strand = key
            messages.append(
                f"    {chrom}:{start}-{end}: "
                f"rustar-aligner={rustar_aligner_count}, STAR={star_count} "
                f"({diff_pct*100:+.1f}%)"
            )

    # Calculate statistics by motif type
    messages.append("\n  Junction motif distribution:")
    star_motif_counts = defaultdict(int)
    rustar_aligner_motif_counts = defaultdict(int)

    for junc in star_junctions.values():
        star_motif_counts[junc.motif] += 1

    for junc in rustar_aligner_junctions.values():
        rustar_aligner_motif_counts[junc.motif] += 1

    for motif in sorted(set(star_motif_counts.keys()) | set(rustar_aligner_motif_counts.keys())):
        motif_name = MOTIF_NAMES.get(motif, f"motif_{motif}")
        star_count = star_motif_counts[motif]
        rustar_aligner_count = rustar_aligner_motif_counts[motif]
        messages.append(f"    {motif_name:20s} rustar-aligner={rustar_aligner_count:>5}, STAR={star_count:>5}")

    # Calculate annotated vs novel
    star_annotated = sum(1 for j in star_junctions.values() if j.annotated)
    rustar_aligner_annotated = sum(1 for j in rustar_aligner_junctions.values() if j.annotated)
    star_novel = len(star_junctions) - star_annotated
    rustar_aligner_novel = len(rustar_aligner_junctions) - rustar_aligner_annotated

    messages.append(f"\n  Annotated junctions:   rustar-aligner={rustar_aligner_annotated}, STAR={star_annotated}")
    messages.append(f"  Novel junctions:       rustar-aligner={rustar_aligner_novel}, STAR={star_novel}")

    # Determine pass/fail
    # Pass if:
    # 1. Overlap rate >= 95%
    # 2. Common junctions have similar counts (most within tolerance)
    passed = True

    if total_junctions > 0:
        if overlap_rate < 95.0:
            passed = False

    # Allow up to 5% of common junctions to have count differences
    if len(common_keys) > 0:
        count_diff_rate = len(count_differences) / len(common_keys)
        if count_diff_rate > 0.05:
            passed = False

    return passed, messages


def main():
    parser = argparse.ArgumentParser(description='Compare SJ.out.tab files from rustar-aligner and STAR')
    parser.add_argument('--star', required=True, help='STAR SJ.out.tab file')
    parser.add_argument('--rustar-aligner', required=True, help='rustar-aligner SJ.out.tab file')
    parser.add_argument('--tolerance', type=float, default=0.10, help='Tolerance for read count differences (default: 0.10 = 10%%)')
    parser.add_argument('--output', help='Output file for comparison report')

    args = parser.parse_args()

    # Parse both files
    print(f"Parsing STAR junctions: {args.star}")
    star_junctions = parse_sj_file(args.star)

    print(f"Parsing rustar-aligner junctions: {args.rustar_aligner}")
    rustar_aligner_junctions = parse_sj_file(args.rustar_aligner)

    # Compare
    passed, messages = compare_junctions(star_junctions, rustar_aligner_junctions, args.tolerance)

    # Add overall status
    messages.append("\n" + "=" * 50)
    if passed:
        messages.append("Status: PASS ✓")
        messages.append("Junction outputs are consistent")
    else:
        messages.append("Status: FAIL ✗")
        messages.append("Significant junction differences detected")

    # Print and optionally save
    output_text = "\n".join(messages)
    print(output_text)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_text)
        print(f"\nComparison report saved to: {args.output}")

    # Exit with appropriate code
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
