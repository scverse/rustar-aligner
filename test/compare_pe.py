#!/usr/bin/env python3
"""
PE-aware SAM comparison: rustar-aligner vs STAR.
Separates mate1 (FLAG 0x40) and mate2 (FLAG 0x80), compares per-mate.
"""

import re
import sys
import os
from collections import defaultdict, Counter


def parse_pe_sam(path):
    """Parse PE SAM, return dict of (qname, mate_num) -> list of records."""
    reads = defaultdict(list)
    unmapped_pairs = 0
    half_mapped = 0
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

            # Determine mate number
            if flag & 0x40:
                mate = 1
            elif flag & 0x80:
                mate = 2
            else:
                mate = 0  # SE or unknown

            is_secondary = bool(flag & 0x100)
            is_unmapped = bool(flag & 0x4)
            is_reverse = bool(flag & 0x10)
            mate_unmapped = bool(flag & 0x8)

            reads[(qname, mate)].append({
                "flag": flag,
                "rname": rname,
                "pos": pos,
                "mapq": mapq,
                "cigar": cigar,
                "is_secondary": is_secondary,
                "is_unmapped": is_unmapped,
                "is_reverse": is_reverse,
                "mate_unmapped": mate_unmapped,
                "mate": mate,
            })
    return reads


def get_primary(records):
    """Get primary alignment from records."""
    for r in records:
        if not r["is_secondary"] and not r["is_unmapped"]:
            return r
    return None


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <rustar_aligner_sam> <star_sam>")
        sys.exit(1)

    rustar_aligner_path = sys.argv[1]
    star_path = sys.argv[2]

    # Auto-append Aligned.out.sam if directory given
    if os.path.isdir(rustar_aligner_path):
        rustar_aligner_path = os.path.join(rustar_aligner_path, "Aligned.out.sam")
    if os.path.isdir(star_path):
        star_path = os.path.join(star_path, "Aligned.out.sam")

    print("=" * 80)
    print("PE-AWARE COMPARISON: rustar-aligner vs STAR")
    print("=" * 80)

    print(f"\nParsing {rustar_aligner_path}...")
    rustar_aligner = parse_pe_sam(rustar_aligner_path)
    print(f"Parsing {star_path}...")
    star = parse_pe_sam(star_path)

    # Count unique read pairs
    rustar_aligner_pairs = set(qname for qname, mate in rustar_aligner.keys())
    star_pairs = set(qname for qname, mate in star.keys())
    print(f"\nrustar-aligner: {len(rustar_aligner_pairs)} read pairs, {sum(len(v) for v in rustar_aligner.values())} total records")
    print(f"STAR:   {len(star_pairs)} read pairs, {sum(len(v) for v in star.values())} total records")

    # === Per-mate comparison ===
    print("\n" + "=" * 80)
    print("PER-MATE COMPARISON")
    print("=" * 80)

    for mate_num in [1, 2]:
        print(f"\n--- Mate {mate_num} ---")

        all_reads = set()
        for (qname, m) in set(rustar_aligner.keys()) | set(star.keys()):
            if m == mate_num:
                all_reads.add(qname)

        agree_pos = 0
        agree_strand = 0
        agree_cigar = 0
        disagree_pos = 0
        rustar_aligner_only = 0
        star_only = 0
        both_unmapped = 0
        disagree_diff_chr = 0
        disagree_same_chr_diff_strand = 0
        disagree_same_chr_same_strand = 0

        for qname in sorted(all_reads):
            r_recs = rustar_aligner.get((qname, mate_num), [])
            s_recs = star.get((qname, mate_num), [])

            r_pri = get_primary(r_recs) if r_recs else None
            s_pri = get_primary(s_recs) if s_recs else None

            r_mapped = r_pri is not None
            s_mapped = s_pri is not None

            # Check if unmapped record exists
            r_unmapped = any(r["is_unmapped"] for r in r_recs)
            s_unmapped = any(r["is_unmapped"] for r in s_recs)

            if r_mapped and s_mapped:
                same_chr = r_pri["rname"] == s_pri["rname"]
                pos_diff = abs(r_pri["pos"] - s_pri["pos"])
                same_pos = same_chr and pos_diff <= 5
                same_strand = r_pri["is_reverse"] == s_pri["is_reverse"]
                same_cigar = r_pri["cigar"] == s_pri["cigar"]

                if same_pos:
                    agree_pos += 1
                    if same_strand:
                        agree_strand += 1
                    if same_cigar:
                        agree_cigar += 1
                else:
                    disagree_pos += 1
                    if not same_chr:
                        disagree_diff_chr += 1
                    elif not same_strand:
                        disagree_same_chr_diff_strand += 1
                    else:
                        disagree_same_chr_same_strand += 1
            elif r_mapped and not s_mapped:
                rustar_aligner_only += 1
            elif not r_mapped and s_mapped:
                star_only += 1
            else:
                both_unmapped += 1

        total = agree_pos + disagree_pos + rustar_aligner_only + star_only + both_unmapped
        both_mapped = agree_pos + disagree_pos

        print(f"  Total mates:         {total}")
        print(f"  Both mapped:         {both_mapped}")
        if both_mapped > 0:
            print(f"    Agree position:    {agree_pos} ({100.0*agree_pos/both_mapped:.1f}%)")
            print(f"    Agree strand:      {agree_strand} ({100.0*agree_strand/both_mapped:.1f}%)")
            print(f"    Agree CIGAR:       {agree_cigar} ({100.0*agree_cigar/both_mapped:.1f}%)")
            print(f"    Disagree position: {disagree_pos} ({100.0*disagree_pos/both_mapped:.1f}%)")
            if disagree_pos > 0:
                print(f"      Diff chr:          {disagree_diff_chr}")
                print(f"      Same chr, diff strand: {disagree_same_chr_diff_strand}")
                print(f"      Same chr, same strand: {disagree_same_chr_same_strand}")
        print(f"  rustar-aligner only mapped:  {rustar_aligner_only}")
        print(f"  STAR only mapped:    {star_only}")
        print(f"  Both unmapped:       {both_unmapped}")

    # === Pair-level comparison ===
    print("\n" + "=" * 80)
    print("PAIR-LEVEL COMPARISON")
    print("=" * 80)

    all_pairs = rustar_aligner_pairs | star_pairs

    pair_both_mapped = 0
    pair_agree = 0
    pair_disagree = 0
    pair_rustar_aligner_only = 0
    pair_star_only = 0
    pair_both_unmapped = 0

    # Count half-mapped in rustar-aligner
    rustar_aligner_half_mapped = 0
    for qname in rustar_aligner_pairs:
        r1 = rustar_aligner.get((qname, 1), [])
        r2 = rustar_aligner.get((qname, 2), [])
        r1_pri = get_primary(r1) if r1 else None
        r2_pri = get_primary(r2) if r2 else None
        r1_unmapped = any(r["is_unmapped"] for r in r1)
        r2_unmapped = any(r["is_unmapped"] for r in r2)

        if (r1_pri and not r2_pri and r2_unmapped) or (r2_pri and not r1_pri and r1_unmapped):
            rustar_aligner_half_mapped += 1

    # Count mapping categories
    rustar_aligner_both_mates_mapped = 0
    rustar_aligner_no_mates_mapped = 0
    star_both_mates_mapped = 0
    star_no_mates_mapped = 0

    for qname in rustar_aligner_pairs:
        r1 = get_primary(rustar_aligner.get((qname, 1), []))
        r2 = get_primary(rustar_aligner.get((qname, 2), []))
        if r1 and r2:
            rustar_aligner_both_mates_mapped += 1
        elif not r1 and not r2:
            rustar_aligner_no_mates_mapped += 1

    for qname in star_pairs:
        s1 = get_primary(star.get((qname, 1), []))
        s2 = get_primary(star.get((qname, 2), []))
        if s1 and s2:
            star_both_mates_mapped += 1
        elif not s1 and not s2:
            star_no_mates_mapped += 1

    print(f"\nrustar-aligner pair-level mapping:")
    print(f"  Both mates mapped:   {rustar_aligner_both_mates_mapped} ({100.0*rustar_aligner_both_mates_mapped/len(rustar_aligner_pairs):.1f}%)")
    print(f"  Half-mapped pairs:   {rustar_aligner_half_mapped} ({100.0*rustar_aligner_half_mapped/len(rustar_aligner_pairs):.1f}%)")
    print(f"  Neither mate mapped: {rustar_aligner_no_mates_mapped} ({100.0*rustar_aligner_no_mates_mapped/len(rustar_aligner_pairs):.1f}%)")

    print(f"\nSTAR pair-level mapping:")
    print(f"  Both mates mapped:   {star_both_mates_mapped} ({100.0*star_both_mates_mapped/len(star_pairs):.1f}%)")
    print(f"  Neither mate mapped: {star_no_mates_mapped} ({100.0*star_no_mates_mapped/len(star_pairs):.1f}%)")

    # === Half-mapped SAM format check ===
    print("\n" + "=" * 80)
    print("HALF-MAPPED PAIR DETAILS (rustar-aligner)")
    print("=" * 80)

    hm_mate1_mapped = 0
    hm_mate2_mapped = 0
    hm_correct_flags = 0
    hm_colocated = 0

    for qname in rustar_aligner_pairs:
        r1 = rustar_aligner.get((qname, 1), [])
        r2 = rustar_aligner.get((qname, 2), [])
        r1_pri = get_primary(r1) if r1 else None
        r2_pri = get_primary(r2) if r2 else None
        r1_unmapped = any(r["is_unmapped"] for r in r1)
        r2_unmapped = any(r["is_unmapped"] for r in r2)

        if r1_pri and not r2_pri and r2_unmapped:
            hm_mate1_mapped += 1
            # Check FLAGS
            unmapped_rec = [r for r in r2 if r["is_unmapped"]]
            if unmapped_rec:
                ur = unmapped_rec[0]
                # Unmapped mate should have FLAG & 0x4 and be co-located
                if ur["is_unmapped"] and r1_pri["mate_unmapped"]:
                    hm_correct_flags += 1
                if ur["rname"] == r1_pri["rname"] and ur["pos"] == r1_pri["pos"]:
                    hm_colocated += 1
        elif r2_pri and not r1_pri and r1_unmapped:
            hm_mate2_mapped += 1
            unmapped_rec = [r for r in r1 if r["is_unmapped"]]
            if unmapped_rec:
                ur = unmapped_rec[0]
                if ur["is_unmapped"] and r2_pri["mate_unmapped"]:
                    hm_correct_flags += 1
                if ur["rname"] == r2_pri["rname"] and ur["pos"] == r2_pri["pos"]:
                    hm_colocated += 1

    total_hm = hm_mate1_mapped + hm_mate2_mapped
    print(f"\n  Half-mapped pairs:    {total_hm}")
    print(f"    Mate 1 mapped:      {hm_mate1_mapped}")
    print(f"    Mate 2 mapped:      {hm_mate2_mapped}")
    if total_hm > 0:
        print(f"    Correct FLAGS:      {hm_correct_flags} ({100.0*hm_correct_flags/total_hm:.1f}%)")
        print(f"    Co-located:         {hm_colocated} ({100.0*hm_colocated/total_hm:.1f}%)")

    # === Summary ===
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Overall per-mate position agreement
    total_mates_compared = 0
    total_mates_agree = 0
    total_mates_cigar_agree = 0

    for mate_num in [1, 2]:
        all_reads = set()
        for (qname, m) in set(rustar_aligner.keys()) | set(star.keys()):
            if m == mate_num:
                all_reads.add(qname)

        for qname in all_reads:
            r_pri = get_primary(rustar_aligner.get((qname, mate_num), []))
            s_pri = get_primary(star.get((qname, mate_num), []))

            if r_pri and s_pri:
                total_mates_compared += 1
                same_chr = r_pri["rname"] == s_pri["rname"]
                pos_diff = abs(r_pri["pos"] - s_pri["pos"])
                if same_chr and pos_diff <= 5:
                    total_mates_agree += 1
                    if r_pri["cigar"] == s_pri["cigar"]:
                        total_mates_cigar_agree += 1

    if total_mates_compared > 0:
        print(f"\n  Per-mate position agreement: {total_mates_agree}/{total_mates_compared} ({100.0*total_mates_agree/total_mates_compared:.1f}%)")
        print(f"  Per-mate CIGAR agreement:    {total_mates_cigar_agree}/{total_mates_compared} ({100.0*total_mates_cigar_agree/total_mates_compared:.1f}%)")

    print(f"\n  rustar-aligner mapped pairs:         {rustar_aligner_both_mates_mapped + rustar_aligner_half_mapped}/{len(rustar_aligner_pairs)} ({100.0*(rustar_aligner_both_mates_mapped + rustar_aligner_half_mapped)/len(rustar_aligner_pairs):.1f}%)")
    print(f"    Both mates mapped:         {rustar_aligner_both_mates_mapped} ({100.0*rustar_aligner_both_mates_mapped/len(rustar_aligner_pairs):.1f}%)")
    print(f"    Half-mapped:               {rustar_aligner_half_mapped} ({100.0*rustar_aligner_half_mapped/len(rustar_aligner_pairs):.1f}%)")
    print(f"  STAR mapped pairs:           {star_both_mates_mapped}/{len(star_pairs)} ({100.0*star_both_mates_mapped/len(star_pairs):.1f}%)")


if __name__ == "__main__":
    main()
