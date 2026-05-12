#!/usr/bin/env python3
"""
Diagnose position disagreements between rustar-aligner and STAR.

For reads where both tools map to different chromosomes with MAPQ=255,
this script checks the genome sequence at both positions to determine:
- Hypothesis A: Read maps perfectly at BOTH loci (multi-mapping not detected)
- Hypothesis B: Read has a clearly better match at one locus (seed finding bug)

Usage:
    python3 diagnose_disagreements.py <rustar_aligner_dir> <star_dir> <genome_fasta>
"""

import sys
import os
import re
from collections import defaultdict, Counter


def parse_fasta(path):
    """Parse FASTA file, return dict of chr_name -> sequence."""
    sequences = {}
    current_name = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_name is not None:
        sequences[current_name] = "".join(current_seq)
    return sequences


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def parse_sam(path):
    """Parse SAM file, return dict of read_name -> list of alignment records."""
    reads = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]
            reads[qname].append({
                "flag": flag,
                "rname": rname,
                "pos": pos,
                "mapq": mapq,
                "cigar": cigar,
                "seq": seq,
            })
    return reads


def get_primary(records):
    """Get primary alignment from records."""
    for r in records:
        if not (r["flag"] & 256) and not (r["flag"] & 4):
            return r
    return None


def cigar_to_ref_len(cigar):
    """Calculate reference-consuming length from CIGAR."""
    ref_len = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in ("M", "D", "N", "=", "X"):
            ref_len += length
    return ref_len


def cigar_to_query_len(cigar):
    """Calculate query-consuming length from CIGAR."""
    query_len = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in ("M", "I", "S", "=", "X"):
            query_len += length
    return query_len


def count_mismatches_at_position(read_seq, genome_seq, chrom, pos, cigar, is_reverse):
    """
    Count mismatches between read and genome at the given alignment position.

    Args:
        read_seq: Read sequence from SAM SEQ field (always forward strand of read)
        genome_seq: Dict of chr -> sequence
        chrom: Chromosome name
        pos: 1-based alignment position
        cigar: CIGAR string
        is_reverse: Whether read is reverse-complemented

    Returns:
        (n_mismatch, n_compared) or None if position is invalid
    """
    if chrom not in genome_seq:
        return None

    genome = genome_seq[chrom]
    genome_pos = pos - 1  # Convert to 0-based
    read_pos = 0
    n_mismatch = 0
    n_compared = 0

    # The SEQ in SAM is on the forward strand of the mapping.
    # If FLAG & 16 (reverse), SAM SEQ is the reverse complement of the original read.
    seq = read_seq

    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in ("M", "=", "X"):
            for i in range(length):
                if genome_pos + i >= len(genome) or read_pos + i >= len(seq):
                    break
                gb = genome[genome_pos + i]
                rb = seq[read_pos + i]
                if rb != "N" and gb != "N":
                    n_compared += 1
                    if rb != gb:
                        n_mismatch += 1
            read_pos += length
            genome_pos += length
        elif op == "I":
            read_pos += length
        elif op == "D":
            genome_pos += length
        elif op == "N":
            genome_pos += length
        elif op == "S":
            read_pos += length
        elif op == "H":
            pass

    return (n_mismatch, n_compared)


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <rustar_aligner_dir> <star_dir> <genome_fasta>")
        sys.exit(1)

    rustar_aligner_dir = sys.argv[1]
    star_dir = sys.argv[2]
    genome_fasta = sys.argv[3]

    rustar_aligner_sam = os.path.join(rustar_aligner_dir, "Aligned.out.sam")
    star_sam = os.path.join(star_dir, "Aligned.out.sam")

    for f in [rustar_aligner_sam, star_sam, genome_fasta]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            sys.exit(1)

    print("=" * 80)
    print("DIAGNOSTIC: Position Disagreement Analysis")
    print("=" * 80)

    # Parse genome
    print("\nLoading genome...")
    genome_seq = parse_fasta(genome_fasta)
    print(f"  Loaded {len(genome_seq)} chromosomes: {', '.join(sorted(genome_seq.keys()))}")
    for chrom, seq in sorted(genome_seq.items()):
        print(f"    {chrom}: {len(seq):,} bp")

    # Parse SAM files
    print("\nParsing SAM files...")
    rustar_aligner_reads = parse_sam(rustar_aligner_sam)
    star_reads = parse_sam(star_sam)
    print(f"  rustar-aligner: {len(rustar_aligner_reads)} reads")
    print(f"  STAR:   {len(star_reads)} reads")

    # Find diff-chr disagreements with both MAPQ=255
    diff_chr_both_unique = []
    diff_chr_other = []

    all_reads = set(rustar_aligner_reads.keys()) & set(star_reads.keys())

    for qname in sorted(all_reads):
        r_pri = get_primary(rustar_aligner_reads[qname])
        s_pri = get_primary(star_reads[qname])

        if r_pri is None or s_pri is None:
            continue

        # Both mapped
        if (r_pri["flag"] & 4) or (s_pri["flag"] & 4):
            continue

        # Different chromosome
        if r_pri["rname"] == s_pri["rname"]:
            continue

        if r_pri["mapq"] == 255 and s_pri["mapq"] == 255:
            diff_chr_both_unique.append((qname, r_pri, s_pri))
        else:
            diff_chr_other.append((qname, r_pri, s_pri))

    print(f"\nDifferent-chromosome disagreements: {len(diff_chr_both_unique) + len(diff_chr_other)}")
    print(f"  Both MAPQ=255: {len(diff_chr_both_unique)}")
    print(f"  Other MAPQ:    {len(diff_chr_other)}")

    # Analyze mismatch counts at both positions
    print("\n" + "=" * 80)
    print("MISMATCH ANALYSIS at both alignment positions")
    print("=" * 80)

    categories = Counter()
    examples = {"both_perfect": [], "rustar_aligner_better": [], "star_better": [], "both_imperfect": []}

    for qname, r_pri, s_pri in diff_chr_both_unique:
        r_is_rev = bool(r_pri["flag"] & 16)
        s_is_rev = bool(s_pri["flag"] & 16)

        r_mm = count_mismatches_at_position(
            r_pri["seq"], genome_seq, r_pri["rname"], r_pri["pos"], r_pri["cigar"], r_is_rev
        )
        s_mm = count_mismatches_at_position(
            s_pri["seq"], genome_seq, s_pri["rname"], s_pri["pos"], s_pri["cigar"], s_is_rev
        )

        if r_mm is None or s_mm is None:
            categories["invalid_position"] += 1
            continue

        r_mismatches, r_compared = r_mm
        s_mismatches, s_compared = s_mm

        if r_mismatches == 0 and s_mismatches == 0:
            cat = "both_perfect"
        elif r_mismatches < s_mismatches:
            cat = "rustar_aligner_better"
        elif s_mismatches < r_mismatches:
            cat = "star_better"
        else:
            cat = "both_imperfect"

        categories[cat] += 1

        if len(examples[cat]) < 5:
            examples[cat].append({
                "qname": qname,
                "rustar_aligner_chr": r_pri["rname"],
                "rustar_aligner_pos": r_pri["pos"],
                "rustar_aligner_strand": "-" if r_is_rev else "+",
                "rustar_aligner_cigar": r_pri["cigar"],
                "rustar_aligner_mm": r_mismatches,
                "rustar_aligner_compared": r_compared,
                "star_chr": s_pri["rname"],
                "star_pos": s_pri["pos"],
                "star_strand": "-" if s_is_rev else "+",
                "star_cigar": s_pri["cigar"],
                "star_mm": s_mismatches,
                "star_compared": s_compared,
            })

    total_analyzed = sum(categories.values())
    print(f"\nOf {len(diff_chr_both_unique)} diff-chr both-MAPQ=255 reads:")
    print(f"\n{'Category':<30} {'Count':>8} {'%':>8}")
    print("-" * 50)
    for cat in ["both_perfect", "rustar_aligner_better", "star_better", "both_imperfect", "invalid_position"]:
        count = categories.get(cat, 0)
        pct = 100.0 * count / max(total_analyzed, 1)
        print(f"  {cat:<28} {count:>8} {pct:>7.1f}%")

    # Print examples for each category
    for cat, label in [
        ("both_perfect", "BOTH PERFECT (0 mismatches at both loci)"),
        ("rustar_aligner_better", "rustar-aligner BETTER (fewer mismatches)"),
        ("star_better", "STAR BETTER (fewer mismatches)"),
        ("both_imperfect", "BOTH IMPERFECT (same # mismatches > 0)"),
    ]:
        if examples.get(cat):
            print(f"\n--- Examples: {label} ---")
            for ex in examples[cat][:5]:
                print(f"  {ex['qname'][:40]}")
                print(f"    rustar-aligner: {ex['rustar_aligner_chr']}:{ex['rustar_aligner_pos']}({ex['rustar_aligner_strand']}) "
                      f"CIGAR={ex['rustar_aligner_cigar']} mm={ex['rustar_aligner_mm']}/{ex['rustar_aligner_compared']}")
                print(f"    STAR:   {ex['star_chr']}:{ex['star_pos']}({ex['star_strand']}) "
                      f"CIGAR={ex['star_cigar']} mm={ex['star_mm']}/{ex['star_compared']}")

    # Hypothesis verdict
    print("\n" + "=" * 80)
    print("HYPOTHESIS VERDICT")
    print("=" * 80)

    both_perfect_count = categories.get("both_perfect", 0)
    both_perfect_pct = 100.0 * both_perfect_count / max(total_analyzed, 1)

    if both_perfect_pct > 50:
        print(f"\n>>> HYPOTHESIS A CONFIRMED: {both_perfect_pct:.1f}% of diff-chr reads have")
        print(f"    perfect matches at BOTH loci.")
        print(f"    These are genuine multi-mappers that should have MAPQ < 255.")
        print(f"    Fix: Ensure cluster_seeds() explores enough positions to find all loci.")
    elif categories.get("star_better", 0) > categories.get("rustar_aligner_better", 0):
        print(f"\n>>> HYPOTHESIS B (STAR better): STAR finds better alignments.")
        print(f"    Suggests rustar-aligner seed finding or scoring has gaps.")
    else:
        print(f"\n>>> MIXED RESULTS: No single dominant pattern.")
        print(f"    Both perfect: {both_perfect_pct:.1f}%")
        print(f"    rustar-aligner better: {100.0 * categories.get('rustar_aligner_better', 0) / max(total_analyzed, 1):.1f}%")
        print(f"    STAR better: {100.0 * categories.get('star_better', 0) / max(total_analyzed, 1):.1f}%")

    # Additional analysis: which chromosome pairs are most common?
    print("\n" + "=" * 80)
    print("CHROMOSOME PAIR ANALYSIS")
    print("=" * 80)

    chr_pairs = Counter()
    for qname, r_pri, s_pri in diff_chr_both_unique:
        pair = tuple(sorted([r_pri["rname"], s_pri["rname"]]))
        chr_pairs[pair] += 1

    print(f"\nTop chromosome pairs for diff-chr disagreements:")
    for pair, count in chr_pairs.most_common(20):
        print(f"  {pair[0]:<15} <-> {pair[1]:<15}  {count:>6}")

    # Analysis for same-chr, far-apart disagreements
    print("\n" + "=" * 80)
    print("SAME-CHR FAR-APART ANALYSIS (>500bp)")
    print("=" * 80)

    same_chr_far = []
    for qname in sorted(all_reads):
        r_pri = get_primary(rustar_aligner_reads[qname])
        s_pri = get_primary(star_reads[qname])
        if r_pri is None or s_pri is None:
            continue
        if (r_pri["flag"] & 4) or (s_pri["flag"] & 4):
            continue
        if r_pri["rname"] != s_pri["rname"]:
            continue
        pos_diff = abs(r_pri["pos"] - s_pri["pos"])
        if pos_diff > 500:
            same_chr_far.append((qname, r_pri, s_pri, pos_diff))

    if same_chr_far:
        offset_counter = Counter()
        for qname, r_pri, s_pri, pos_diff in same_chr_far:
            offset_counter[pos_diff] += 1

        print(f"\n{len(same_chr_far)} reads differ by >500bp on same chromosome")
        print(f"\nTop offset distances:")
        for offset, count in offset_counter.most_common(10):
            print(f"  {offset:>10}bp  {count:>6} reads")

        # Check mismatches for first few
        print(f"\nMismatch analysis for same-chr far-apart (first 10):")
        for qname, r_pri, s_pri, pos_diff in same_chr_far[:10]:
            r_is_rev = bool(r_pri["flag"] & 16)
            s_is_rev = bool(s_pri["flag"] & 16)
            r_mm = count_mismatches_at_position(
                r_pri["seq"], genome_seq, r_pri["rname"], r_pri["pos"], r_pri["cigar"], r_is_rev
            )
            s_mm = count_mismatches_at_position(
                s_pri["seq"], genome_seq, s_pri["rname"], s_pri["pos"], s_pri["cigar"], s_is_rev
            )
            r_str = f"mm={r_mm[0]}/{r_mm[1]}" if r_mm else "invalid"
            s_str = f"mm={s_mm[0]}/{s_mm[1]}" if s_mm else "invalid"
            print(f"  {qname[:40]} diff={pos_diff}bp chr={r_pri['rname']}")
            print(f"    rustar-aligner: pos={r_pri['pos']} {r_str} CIGAR={r_pri['cigar']}")
            print(f"    STAR:   pos={s_pri['pos']} {s_str} CIGAR={s_pri['cigar']}")
    else:
        print("\nNo same-chr far-apart disagreements found.")


if __name__ == "__main__":
    main()
