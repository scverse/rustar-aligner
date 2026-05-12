#!/usr/bin/env python3

"""
compare_chimeric.py - Compare Chimeric.out.junction files from rustar-aligner and STAR

Chimeric.out.junction format (14 columns):
1. chr_donorA
2. brkpt_donorA
3. strand_donorA
4. chr_acceptorB
5. brkpt_acceptorB
6. strand_acceptorB
7. junction_type (-1=chimeric, 0=GT/AG, 1=CT/AC, 2=GC/AG, 3=CT/GC, 4=AT/AC, 5=non-canonical)
8. repeat_length_left
9. repeat_length_right
10. read_name
11. start_alnA
12. cigar_alnA
13. start_alnB
14. cigar_alnB

Usage:
    python compare_chimeric.py --star STAR_Chimeric.out.junction --rustar-aligner rustar-aligner_Chimeric.out.junction [options]
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set


@dataclass
class ChimericJunction:
    """Chimeric junction record."""
    chr_donor: str
    brkpt_donor: int
    strand_donor: str
    chr_acceptor: str
    brkpt_acceptor: int
    strand_acceptor: str
    junction_type: int
    repeat_left: int
    repeat_right: int
    read_name: str
    start_aln_donor: int
    cigar_donor: str
    start_aln_acceptor: int
    cigar_acceptor: str

    def breakpoint_key(self) -> Tuple[str, int, str, str, int, str]:
        """Return breakpoint key for comparison (ignoring read-specific info)."""
        return (
            self.chr_donor, self.brkpt_donor, self.strand_donor,
            self.chr_acceptor, self.brkpt_acceptor, self.strand_acceptor
        )

    def is_inter_chromosomal(self) -> bool:
        """Check if fusion is between different chromosomes."""
        return self.chr_donor != self.chr_acceptor

    def is_strand_break(self) -> bool:
        """Check if fusion has strand switch on same chromosome."""
        return self.chr_donor == self.chr_acceptor and self.strand_donor != self.strand_acceptor

    def distance(self) -> int:
        """Calculate distance for intra-chromosomal same-strand fusions."""
        if self.chr_donor == self.chr_acceptor and self.strand_donor == self.strand_acceptor:
            return abs(self.brkpt_acceptor - self.brkpt_donor)
        return -1


JUNCTION_TYPE_NAMES = {
    -1: "chimeric",
    0: "GT/AG",
    1: "CT/AC",
    2: "GC/AG",
    3: "CT/GC",
    4: "AT/AC",
    5: "non-canonical",
}


def parse_chimeric_file(filename: str) -> List[ChimericJunction]:
    """Parse Chimeric.out.junction file."""
    junctions = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) != 14:
                    continue

                junction = ChimericJunction(
                    chr_donor=fields[0],
                    brkpt_donor=int(fields[1]),
                    strand_donor=fields[2],
                    chr_acceptor=fields[3],
                    brkpt_acceptor=int(fields[4]),
                    strand_acceptor=fields[5],
                    junction_type=int(fields[6]),
                    repeat_left=int(fields[7]),
                    repeat_right=int(fields[8]),
                    read_name=fields[9],
                    start_aln_donor=int(fields[10]),
                    cigar_donor=fields[11],
                    start_aln_acceptor=int(fields[12]),
                    cigar_acceptor=fields[13],
                )

                junctions.append(junction)
    except FileNotFoundError:
        # File may not exist if no chimeric alignments found
        pass

    return junctions


def compare_chimeric(star_junctions: List[ChimericJunction], rustar_aligner_junctions: List[ChimericJunction], tolerance: float) -> Tuple[bool, List[str]]:
    """Compare chimeric junction sets."""
    messages = []
    messages.append("\n=== Chimeric Alignment Comparison: rustar-aligner vs STAR ===")

    # Overall counts
    messages.append(f"  Total chimeric junctions (rustar-aligner): {len(rustar_aligner_junctions)}")
    messages.append(f"  Total chimeric junctions (STAR):   {len(star_junctions)}")

    # If both are empty, that's a match
    if len(star_junctions) == 0 and len(rustar_aligner_junctions) == 0:
        messages.append("  No chimeric alignments detected by either aligner")
        return True, messages

    # Group by breakpoint
    star_breakpoints = defaultdict(list)
    rustar_aligner_breakpoints = defaultdict(list)

    for junc in star_junctions:
        star_breakpoints[junc.breakpoint_key()].append(junc)

    for junc in rustar_aligner_junctions:
        rustar_aligner_breakpoints[junc.breakpoint_key()].append(junc)

    # Find common and unique breakpoints
    star_keys = set(star_breakpoints.keys())
    rustar_aligner_keys = set(rustar_aligner_breakpoints.keys())

    common_keys = star_keys & rustar_aligner_keys
    star_only = star_keys - rustar_aligner_keys
    rustar_aligner_only = rustar_aligner_keys - star_keys

    messages.append(f"\n  Unique breakpoints (rustar-aligner): {len(rustar_aligner_breakpoints)}")
    messages.append(f"  Unique breakpoints (STAR):   {len(star_breakpoints)}")
    messages.append(f"  Common breakpoints:          {len(common_keys)}")

    if star_only:
        messages.append(f"  Only in STAR:                {len(star_only)}")
    if rustar_aligner_only:
        messages.append(f"  Only in rustar-aligner:              {len(rustar_aligner_only)}")

    # Calculate overlap rate
    total_breakpoints = len(star_keys | rustar_aligner_keys)
    if total_breakpoints > 0:
        overlap_rate = 100.0 * len(common_keys) / total_breakpoints
        messages.append(f"  Overlap rate:                {overlap_rate:.1f}%")

    # Categorize by fusion type
    def categorize_junctions(junctions):
        inter_chr = sum(1 for j in junctions if j.is_inter_chromosomal())
        strand_break = sum(1 for j in junctions if j.is_strand_break())
        large_distance = sum(1 for j in junctions if not j.is_inter_chromosomal() and not j.is_strand_break() and j.distance() > 1000000)
        return inter_chr, strand_break, large_distance

    star_inter, star_strand, star_large = categorize_junctions(star_junctions)
    rustar_aligner_inter, rustar_aligner_strand, rustar_aligner_large = categorize_junctions(rustar_aligner_junctions)

    messages.append("\n  Fusion categories:")
    messages.append(f"    Inter-chromosomal:  rustar-aligner={rustar_aligner_inter}, STAR={star_inter}")
    messages.append(f"    Strand breaks:      rustar-aligner={rustar_aligner_strand}, STAR={star_strand}")
    messages.append(f"    Large distance:     rustar-aligner={rustar_aligner_large}, STAR={star_large}")

    # Show unique breakpoints (up to 5 each)
    if rustar_aligner_only:
        messages.append(f"\n  Unique to rustar-aligner ({len(rustar_aligner_only)} total, showing first 5):")
        for key in sorted(rustar_aligner_only)[:5]:
            chr_d, brkpt_d, strand_d, chr_a, brkpt_a, strand_a = key
            junctions = rustar_aligner_breakpoints[key]
            messages.append(
                f"    {chr_d}:{brkpt_d}({strand_d}) -> {chr_a}:{brkpt_a}({strand_a}) "
                f"[{len(junctions)} reads]"
            )

    if star_only:
        messages.append(f"\n  Unique to STAR ({len(star_only)} total, showing first 5):")
        for key in sorted(star_only)[:5]:
            chr_d, brkpt_d, strand_d, chr_a, brkpt_a, strand_a = key
            junctions = star_breakpoints[key]
            messages.append(
                f"    {chr_d}:{brkpt_d}({strand_d}) -> {chr_a}:{brkpt_a}({strand_a}) "
                f"[{len(junctions)} reads]"
            )

    # Compare common breakpoints
    read_count_differences = []

    for key in common_keys:
        star_count = len(star_breakpoints[key])
        rustar_aligner_count = len(rustar_aligner_breakpoints[key])

        if star_count > 0:
            diff_pct = abs(rustar_aligner_count - star_count) / star_count
            if diff_pct > tolerance:
                read_count_differences.append((key, star_count, rustar_aligner_count, diff_pct))

    if read_count_differences:
        messages.append(f"\n  Read count differences (>{tolerance*100:.0f}% threshold): {len(read_count_differences)}")
        for key, star_count, rustar_aligner_count, diff_pct in sorted(read_count_differences, key=lambda x: x[3], reverse=True)[:5]:
            chr_d, brkpt_d, strand_d, chr_a, brkpt_a, strand_a = key
            messages.append(
                f"    {chr_d}:{brkpt_d}({strand_d}) -> {chr_a}:{brkpt_a}({strand_a}): "
                f"rustar-aligner={rustar_aligner_count}, STAR={star_count} ({diff_pct*100:+.1f}%)"
            )

    # Determine pass/fail
    passed = True

    # If one found chimeric alignments but the other didn't, that's a significant difference
    if (len(star_junctions) == 0) != (len(rustar_aligner_junctions) == 0):
        passed = False

    # If both found chimeric alignments, check overlap
    if len(star_junctions) > 0 and len(rustar_aligner_junctions) > 0:
        if total_breakpoints > 0 and overlap_rate < 80.0:
            # More lenient threshold for chimeric alignments (can be harder to detect)
            passed = False

    return passed, messages


def main():
    parser = argparse.ArgumentParser(description='Compare Chimeric.out.junction files from rustar-aligner and STAR')
    parser.add_argument('--star', required=True, help='STAR Chimeric.out.junction file')
    parser.add_argument('--rustar-aligner', required=True, help='rustar-aligner Chimeric.out.junction file')
    parser.add_argument('--tolerance', type=float, default=0.20, help='Tolerance for read count differences (default: 0.20 = 20%%)')
    parser.add_argument('--output', help='Output file for comparison report')

    args = parser.parse_args()

    # Parse both files
    print(f"Parsing STAR chimeric junctions: {args.star}")
    star_junctions = parse_chimeric_file(args.star)

    print(f"Parsing rustar-aligner chimeric junctions: {args.rustar_aligner}")
    rustar_aligner_junctions = parse_chimeric_file(args.rustar_aligner)

    # Compare
    passed, messages = compare_chimeric(star_junctions, rustar_aligner_junctions, args.tolerance)

    # Add overall status
    messages.append("\n" + "=" * 50)
    if passed:
        messages.append("Status: PASS ✓")
        messages.append("Chimeric outputs are consistent")
    else:
        messages.append("Status: FAIL ✗")
        messages.append("Significant chimeric differences detected")

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
