#!/usr/bin/env python3
"""
Extract disagreement reads from rustar-aligner vs STAR comparison for targeted debugging.

Usage:
    python3 extract_disagreement_reads.py <rustar_aligner_sam> <star_sam> <input_fastq> <output_dir>

Reads both SAM files, identifies reads where rustar-aligner and STAR disagree on alignment,
categorizes them, and extracts representative reads into small FASTQ files for
targeted re-alignment with --readNameFilter.

Categories:
  - false_splice: rustar-aligner splices with large intron, STAR doesn't (or vice versa)
  - missed_splice: STAR splices, rustar-aligner doesn't
  - same_chr_close: same chromosome, 6-50bp position difference
  - same_chr_far: same chromosome, >50bp position difference
  - star_only: STAR maps read, rustar-aligner doesn't
  - rustar_aligner_only: rustar-aligner maps read, STAR doesn't
  - diff_chr_tie: different chromosome, both multi-mapped (tie-breaking)
"""

import os
import re
import sys
from collections import defaultdict


def parse_sam_primary(path):
    """Parse SAM, return dict of read_name -> primary alignment record."""
    reads = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            flag = int(fields[1])
            # Skip secondary/supplementary
            if flag & 0x900:
                continue
            reads[qname] = {
                "flag": flag,
                "rname": fields[2],
                "pos": int(fields[3]),
                "mapq": int(fields[4]),
                "cigar": fields[5],
            }
    return reads


def has_splice(cigar):
    """Check if CIGAR contains a splice junction (N operation)."""
    return "N" in cigar


def max_intron(cigar):
    """Return the largest intron (N) size in a CIGAR, or 0 if no introns."""
    introns = [int(x) for x in re.findall(r'(\d+)N', cigar)]
    return max(introns) if introns else 0


def categorize_disagreements(rustar_aligner_reads, star_reads):
    """Categorize all disagreement reads."""
    categories = defaultdict(list)

    all_names = set(rustar_aligner_reads.keys()) | set(star_reads.keys())

    for qname in sorted(all_names):
        r = rustar_aligner_reads.get(qname)
        s = star_reads.get(qname)

        r_mapped = r is not None and not (r["flag"] & 4)
        s_mapped = s is not None and not (s["flag"] & 4)

        if not r_mapped and not s_mapped:
            continue  # Both unmapped

        if r_mapped and not s_mapped:
            categories["rustar_aligner_only"].append((qname, r, None))
            continue

        if not r_mapped and s_mapped:
            categories["star_only"].append((qname, None, s))
            continue

        # Both mapped — check agreement
        same_chr = r["rname"] == s["rname"]
        pos_diff = abs(r["pos"] - s["pos"]) if same_chr else float("inf")

        if same_chr and pos_diff <= 5:
            continue  # Agree

        if not same_chr:
            # Different chromosome
            if r["mapq"] < 255 and s["mapq"] < 255:
                categories["diff_chr_tie"].append((qname, r, s))
            else:
                categories["diff_chr_unique"].append((qname, r, s))
            continue

        # Same chr, disagree on position
        r_spliced = has_splice(r["cigar"])
        s_spliced = has_splice(s["cigar"])

        if r_spliced and not s_spliced and max_intron(r["cigar"]) > 10000:
            categories["false_splice"].append((qname, r, s))
        elif s_spliced and not r_spliced:
            categories["missed_splice"].append((qname, r, s))
        elif r_spliced and s_spliced and pos_diff > 500:
            # Both spliced but very different positions
            categories["same_chr_far"].append((qname, r, s))
        elif pos_diff <= 50:
            categories["same_chr_close"].append((qname, r, s))
        else:
            categories["same_chr_far"].append((qname, r, s))

    return categories


def extract_reads_from_fastq(fastq_path, read_names, output_path):
    """Extract specific reads from a FASTQ file (supports .gz)."""
    import gzip

    target = set(read_names)
    found = 0

    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fin, open(output_path, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            # Extract read name (strip @ and /1 /2 suffix)
            name = header.strip().lstrip("@").split()[0]
            if name in target:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                found += 1

    return found


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <rustar_aligner_sam> <star_sam> <input_fastq> <output_dir>")
        sys.exit(1)

    rustar_aligner_sam = sys.argv[1]
    star_sam = sys.argv[2]
    input_fastq = sys.argv[3]
    output_dir = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    print("Parsing SAM files...")
    rustar_aligner_reads = parse_sam_primary(rustar_aligner_sam)
    star_reads = parse_sam_primary(star_sam)
    print(f"  rustar-aligner: {len(rustar_aligner_reads)} primary alignments")
    print(f"  STAR:   {len(star_reads)} primary alignments")

    print("\nCategorizing disagreements...")
    categories = categorize_disagreements(rustar_aligner_reads, star_reads)

    # Print summary
    print("\n" + "=" * 70)
    print("DISAGREEMENT CATEGORY SUMMARY")
    print("=" * 70)
    total = 0
    for cat in ["false_splice", "missed_splice", "same_chr_close", "same_chr_far",
                 "diff_chr_tie", "diff_chr_unique", "star_only", "rustar_aligner_only"]:
        count = len(categories.get(cat, []))
        total += count
        print(f"  {cat:<25} {count:>6}")
    print(f"  {'TOTAL':<25} {total:>6}")

    # Print examples for each category
    for cat in ["false_splice", "missed_splice", "same_chr_close", "same_chr_far",
                 "star_only", "rustar_aligner_only"]:
        items = categories.get(cat, [])
        if not items:
            continue
        print(f"\n--- {cat} (first 5 of {len(items)}) ---")
        for qname, r, s in items[:5]:
            if r and s:
                print(f"  {qname}")
                print(f"    rustar-aligner: {r['rname']}:{r['pos']} MAPQ={r['mapq']} CIGAR={r['cigar']}")
                print(f"    STAR:   {s['rname']}:{s['pos']} MAPQ={s['mapq']} CIGAR={s['cigar']}")
            elif r:
                print(f"  {qname}")
                print(f"    rustar-aligner: {r['rname']}:{r['pos']} MAPQ={r['mapq']} CIGAR={r['cigar']}")
                print(f"    STAR:   unmapped")
            elif s:
                print(f"  {qname}")
                print(f"    rustar-aligner: unmapped")
                print(f"    STAR:   {s['rname']}:{s['pos']} MAPQ={s['mapq']} CIGAR={s['cigar']}")

    # Collect representative reads for debugging
    debug_reads = []
    for cat in ["false_splice", "missed_splice", "same_chr_close", "same_chr_far",
                 "star_only", "rustar_aligner_only"]:
        items = categories.get(cat, [])
        # Take up to 3 per category
        for qname, r, s in items[:3]:
            debug_reads.append(qname)

    # Write read names file
    names_path = os.path.join(output_dir, "disagreement_reads.txt")
    with open(names_path, "w") as f:
        for name in debug_reads:
            f.write(name + "\n")
    print(f"\nWrote {len(debug_reads)} debug read names to {names_path}")

    # Write all disagreement read names by category
    for cat, items in categories.items():
        cat_path = os.path.join(output_dir, f"{cat}_reads.txt")
        with open(cat_path, "w") as f:
            for qname, r, s in items:
                f.write(qname + "\n")

    # Extract reads from FASTQ
    if os.path.exists(input_fastq):
        fq_path = os.path.join(output_dir, "disagreement_reads.fq")
        found = extract_reads_from_fastq(input_fastq, debug_reads, fq_path)
        print(f"Extracted {found}/{len(debug_reads)} reads to {fq_path}")
    else:
        print(f"\nWARNING: FASTQ file not found: {input_fastq}")
        print("  Skipping FASTQ extraction. Run manually when FASTQ is available.")

    print(f"\nTo debug a specific read:")
    print(f"  cargo run --release -- --genomeDir <genome_dir> --readFilesIn {input_fastq} \\")
    print(f"    --readNameFilter <READ_NAME> 2>debug.log")


if __name__ == "__main__":
    main()
