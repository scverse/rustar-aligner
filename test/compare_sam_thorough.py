#!/usr/bin/env python3
"""
Thorough comparison of rustar-aligner vs STAR SAM output.

Usage:
    python3 compare_sam_thorough.py <rustar_aligner_dir> <star_dir>

Both directories should contain Aligned.out.sam and SJ.out.tab.
"""

import re
import sys
import os
from collections import defaultdict, Counter

def parse_sam(path):
    """Parse SAM file, return dict of read_name -> list of alignment records."""
    reads = defaultdict(list)
    header_lines = 0
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                header_lines += 1
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
            reads[qname].append({
                "flag": flag,
                "rname": rname,
                "pos": pos,
                "mapq": mapq,
                "cigar": cigar,
            })
    return reads, header_lines


def classify_read(records):
    """Classify a read as unique/multi/unmapped based on its records."""
    if len(records) == 1 and (records[0]["flag"] & 4):
        return "unmapped"
    for r in records:
        if not (r["flag"] & 256):  # not secondary
            if r["flag"] & 4:
                return "unmapped"
            if r["mapq"] == 255:
                return "unique"
            return "multi"
    return "unmapped"


def get_primary(records):
    """Get primary alignment from records."""
    for r in records:
        if not (r["flag"] & 256) and not (r["flag"] & 4):
            return r
    return None


def cigar_category(cigar):
    """Categorize a CIGAR string."""
    cats = set()
    if cigar == "*":
        return {"unmapped"}
    ops = re.findall(r'\d+([MIDNSHP=X])', cigar)
    if "N" in ops:
        cats.add("spliced")
    if "I" in ops or "D" in ops:
        cats.add("indel")
    if "S" in ops:
        cats.add("softclip")
    if not cats:
        cats.add("pure_match")
    return cats


def parse_sj(path):
    """Parse SJ.out.tab, return dict of (chr, start, end) -> record."""
    junctions = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = int(fields[3])
            motif = int(fields[4])
            annotated = int(fields[5])
            uniq_count = int(fields[6])
            multi_count = int(fields[7])
            max_overhang = int(fields[8])
            key = (chrom, start, end)
            junctions[key] = {
                "strand": strand,
                "motif": motif,
                "annotated": annotated,
                "uniq_count": uniq_count,
                "multi_count": multi_count,
                "max_overhang": max_overhang,
            }
    return junctions


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <rustar_aligner_dir> <star_dir>")
        sys.exit(1)

    rustar_aligner_dir = sys.argv[1]
    star_dir = sys.argv[2]

    rustar_aligner_sam = os.path.join(rustar_aligner_dir, "Aligned.out.sam")
    star_sam = os.path.join(star_dir, "Aligned.out.sam")
    rustar_aligner_sj_path = os.path.join(rustar_aligner_dir, "SJ.out.tab")
    star_sj_path = os.path.join(star_dir, "SJ.out.tab")

    for f in [rustar_aligner_sam, star_sam]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            sys.exit(1)

    has_sj = os.path.exists(rustar_aligner_sj_path) and os.path.exists(star_sj_path)

    # ============================================================
    # PARSE FILES
    # ============================================================
    print("=" * 80)
    print("COMPREHENSIVE COMPARISON: rustar-aligner vs STAR")
    print("=" * 80)
    print(f"\nrustar-aligner dir: {rustar_aligner_dir}")
    print(f"STAR dir:   {star_dir}")

    print("\nParsing SAM files...")
    rustar_aligner_reads, rustar_aligner_headers = parse_sam(rustar_aligner_sam)
    star_reads, star_headers = parse_sam(star_sam)

    print(f"  rustar-aligner: {len(rustar_aligner_reads)} unique read names, {rustar_aligner_headers} header lines")
    print(f"  STAR:   {len(star_reads)} unique read names, {star_headers} header lines")

    # ============================================================
    # 1. MAPPING STATS COMPARISON
    # ============================================================
    print("\n" + "=" * 80)
    print("1. MAPPING STATS COMPARISON")
    print("=" * 80)

    rustar_aligner_class = Counter()
    star_class = Counter()

    for qname, records in rustar_aligner_reads.items():
        rustar_aligner_class[classify_read(records)] += 1

    for qname, records in star_reads.items():
        star_class[classify_read(records)] += 1

    rustar_aligner_total = sum(rustar_aligner_class.values())
    star_total = sum(star_class.values())

    print(f"\n{'Category':<15} {'rustar-aligner':>10} {'%':>8} {'STAR':>10} {'%':>8} {'Diff':>8}")
    print("-" * 65)
    for cat in ["unique", "multi", "unmapped"]:
        rc = rustar_aligner_class[cat]
        sc = star_class[cat]
        rp = 100.0 * rc / rustar_aligner_total if rustar_aligner_total else 0
        sp = 100.0 * sc / star_total if star_total else 0
        print(f"{cat:<15} {rc:>10} {rp:>7.1f}% {sc:>10} {sp:>7.1f}% {rc - sc:>+8}")
    print(f"{'TOTAL':<15} {rustar_aligner_total:>10} {'':>8} {star_total:>10}")

    # ============================================================
    # 2. PER-READ AGREEMENT
    # ============================================================
    print("\n" + "=" * 80)
    print("2. PER-READ AGREEMENT")
    print("=" * 80)

    all_reads = set(rustar_aligner_reads.keys()) | set(star_reads.keys())

    both_mapped_agree_pos = 0
    both_mapped_agree_exact = 0
    both_mapped_disagree_pos = 0
    star_only_mapped = 0
    rustar_aligner_only_mapped = 0
    both_unmapped = 0
    both_mapped_agree_strand = 0
    both_mapped_agree_cigar = 0

    # Categorized position disagreements
    disagree_diff_chr = 0
    disagree_diff_strand_same_chr = 0
    disagree_same_chr_1bp = 0       # 1 bp off
    disagree_same_chr_2_5bp = 0     # 2-5 bp off
    disagree_same_chr_6_50bp = 0    # 6-50 bp off
    disagree_same_chr_51_500bp = 0  # 51-500 bp off
    disagree_same_chr_500plus = 0   # >500 bp off

    disagree_examples = []
    # Collect examples per category for deeper analysis
    disagree_diff_chr_examples = []
    disagree_same_chr_close_examples = []     # 1-50bp
    disagree_same_chr_far_examples = []       # 51+bp

    for qname in sorted(all_reads):
        r_records = rustar_aligner_reads.get(qname, [])
        s_records = star_reads.get(qname, [])

        r_class = classify_read(r_records) if r_records else "missing"
        s_class = classify_read(s_records) if s_records else "missing"

        r_mapped = r_class in ("unique", "multi")
        s_mapped = s_class in ("unique", "multi")

        if r_mapped and s_mapped:
            r_pri = get_primary(r_records)
            s_pri = get_primary(s_records)
            if r_pri and s_pri:
                same_chr = r_pri["rname"] == s_pri["rname"]
                pos_diff = abs(r_pri["pos"] - s_pri["pos"])
                same_pos = pos_diff <= 5
                same_strand = (r_pri["flag"] & 16) == (s_pri["flag"] & 16)
                same_cigar = r_pri["cigar"] == s_pri["cigar"]

                if same_chr and same_pos:
                    both_mapped_agree_pos += 1
                    if same_strand:
                        both_mapped_agree_strand += 1
                    if same_cigar:
                        both_mapped_agree_cigar += 1
                    if same_chr and same_pos and same_strand and same_cigar:
                        both_mapped_agree_exact += 1
                else:
                    both_mapped_disagree_pos += 1

                    # Build example record
                    ex = {
                        "qname": qname,
                        "rustar_aligner_chr": r_pri["rname"],
                        "rustar_aligner_pos": r_pri["pos"],
                        "rustar_aligner_strand": "-" if (r_pri["flag"] & 16) else "+",
                        "rustar_aligner_cigar": r_pri["cigar"],
                        "rustar_aligner_mapq": r_pri["mapq"],
                        "rustar_aligner_class": classify_read(r_records),
                        "star_chr": s_pri["rname"],
                        "star_pos": s_pri["pos"],
                        "star_strand": "-" if (s_pri["flag"] & 16) else "+",
                        "star_cigar": s_pri["cigar"],
                        "star_mapq": s_pri["mapq"],
                        "star_class": classify_read(s_records),
                    }

                    if len(disagree_examples) < 30:
                        disagree_examples.append(ex)

                    # Categorize
                    if not same_chr:
                        disagree_diff_chr += 1
                        if len(disagree_diff_chr_examples) < 10:
                            disagree_diff_chr_examples.append(ex)
                    elif not same_strand:
                        disagree_diff_strand_same_chr += 1
                    else:
                        # Same chr, same strand, different position
                        if pos_diff == 1:
                            disagree_same_chr_1bp += 1
                        elif pos_diff <= 5:
                            disagree_same_chr_2_5bp += 1
                        elif pos_diff <= 50:
                            disagree_same_chr_6_50bp += 1
                            if len(disagree_same_chr_close_examples) < 10:
                                ex["pos_diff"] = pos_diff
                                disagree_same_chr_close_examples.append(ex)
                        elif pos_diff <= 500:
                            disagree_same_chr_51_500bp += 1
                            if len(disagree_same_chr_far_examples) < 10:
                                ex["pos_diff"] = pos_diff
                                disagree_same_chr_far_examples.append(ex)
                        else:
                            disagree_same_chr_500plus += 1
                            if len(disagree_same_chr_far_examples) < 10:
                                ex["pos_diff"] = pos_diff
                                disagree_same_chr_far_examples.append(ex)
            else:
                both_mapped_disagree_pos += 1
                disagree_diff_chr += 1  # Can't compare, count as different
        elif r_mapped and not s_mapped:
            rustar_aligner_only_mapped += 1
        elif not r_mapped and s_mapped:
            star_only_mapped += 1
        else:
            both_unmapped += 1

    total_both_mapped = both_mapped_agree_pos + both_mapped_disagree_pos

    print(f"\n{'Category':<48} {'Count':>8} {'%':>8}")
    print("-" * 68)
    print(f"{'Both mapped, agree position (<=5bp)':<48} {both_mapped_agree_pos:>8} {100.0 * both_mapped_agree_pos / len(all_reads):>7.1f}%")
    print(f"{'  ...of which agree strand':<48} {both_mapped_agree_strand:>8} {100.0 * both_mapped_agree_strand / max(both_mapped_agree_pos, 1):>7.1f}%")
    print(f"{'  ...of which agree CIGAR (exact)':<48} {both_mapped_agree_cigar:>8} {100.0 * both_mapped_agree_cigar / max(both_mapped_agree_pos, 1):>7.1f}%")
    print(f"{'  ...of which agree ALL (chr+pos+strand+cigar)':<48} {both_mapped_agree_exact:>8} {100.0 * both_mapped_agree_exact / max(both_mapped_agree_pos, 1):>7.1f}%")
    print(f"{'Both mapped, DISAGREE position':<48} {both_mapped_disagree_pos:>8} {100.0 * both_mapped_disagree_pos / len(all_reads):>7.1f}%")
    print(f"{'STAR only mapped':<48} {star_only_mapped:>8} {100.0 * star_only_mapped / len(all_reads):>7.1f}%")
    print(f"{'rustar-aligner only mapped':<48} {rustar_aligner_only_mapped:>8} {100.0 * rustar_aligner_only_mapped / len(all_reads):>7.1f}%")
    print(f"{'Both unmapped':<48} {both_unmapped:>8} {100.0 * both_unmapped / len(all_reads):>7.1f}%")
    print("-" * 68)
    print(f"{'Total reads':<48} {len(all_reads):>8}")

    # ============================================================
    # 2b. CATEGORIZED POSITION DISAGREEMENTS
    # ============================================================
    print("\n" + "=" * 80)
    print("2b. POSITION DISAGREEMENT BREAKDOWN")
    print("=" * 80)

    disagree_same_chr_close = disagree_same_chr_1bp + disagree_same_chr_2_5bp + disagree_same_chr_6_50bp
    disagree_same_chr_far = disagree_same_chr_51_500bp + disagree_same_chr_500plus

    print(f"\nOf {both_mapped_disagree_pos} position disagreements:")
    print(f"\n{'Category':<48} {'Count':>8} {'% of disagree':>14}")
    print("-" * 72)
    print(f"{'Different chromosome':<48} {disagree_diff_chr:>8} {100.0 * disagree_diff_chr / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print(f"{'Same chr, different strand':<48} {disagree_diff_strand_same_chr:>8} {100.0 * disagree_diff_strand_same_chr / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print(f"{'Same chr+strand, off by 1bp':<48} {disagree_same_chr_1bp:>8} {100.0 * disagree_same_chr_1bp / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print(f"{'Same chr+strand, off by 2-5bp':<48} {disagree_same_chr_2_5bp:>8} {100.0 * disagree_same_chr_2_5bp / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print(f"{'Same chr+strand, off by 6-50bp':<48} {disagree_same_chr_6_50bp:>8} {100.0 * disagree_same_chr_6_50bp / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print(f"{'Same chr+strand, off by 51-500bp':<48} {disagree_same_chr_51_500bp:>8} {100.0 * disagree_same_chr_51_500bp / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print(f"{'Same chr+strand, off by >500bp':<48} {disagree_same_chr_500plus:>8} {100.0 * disagree_same_chr_500plus / max(both_mapped_disagree_pos, 1):>13.1f}%")
    print("-" * 72)

    # MAPQ breakdown for different-chromosome disagreements
    if disagree_diff_chr_examples:
        diff_chr_both_unique = sum(1 for ex in disagree_diff_chr_examples
                                   if ex["rustar_aligner_mapq"] == 255 and ex["star_mapq"] == 255)
        # Count across ALL diff-chr, not just examples
        mapq_counter = Counter()
        for qname in sorted(all_reads):
            r_records = rustar_aligner_reads.get(qname, [])
            s_records = star_reads.get(qname, [])
            r_class = classify_read(r_records) if r_records else "missing"
            s_class = classify_read(s_records) if s_records else "missing"
            if r_class in ("unique", "multi") and s_class in ("unique", "multi"):
                r_pri = get_primary(r_records)
                s_pri = get_primary(s_records)
                if r_pri and s_pri and r_pri["rname"] != s_pri["rname"]:
                    r_mapq_cat = "unique" if r_pri["mapq"] == 255 else "multi"
                    s_mapq_cat = "unique" if s_pri["mapq"] == 255 else "multi"
                    mapq_counter[f"rustar-aligner={r_mapq_cat}, STAR={s_mapq_cat}"] += 1

        print(f"\nDifferent-chromosome MAPQ breakdown ({disagree_diff_chr} reads):")
        for key, count in mapq_counter.most_common():
            print(f"  {key:<40} {count:>6} ({100.0 * count / max(disagree_diff_chr, 1):.1f}%)")

        # Are these reads that have multiple equally-good alignment loci?
        # Check if the rustar-aligner-chosen locus appears in STAR's multi-alignments or vice versa
        print(f"\nDifferent-chromosome examples (first 10):")
        for ex in disagree_diff_chr_examples[:10]:
            print(f"  {ex['qname'][:30]:<32} rustar-aligner={ex['rustar_aligner_chr']}:{ex['rustar_aligner_pos']}({ex['rustar_aligner_strand']}) MAPQ={ex['rustar_aligner_mapq']} CIGAR={ex['rustar_aligner_cigar']}")
            print(f"  {'':>32} STAR  ={ex['star_chr']}:{ex['star_pos']}({ex['star_strand']}) MAPQ={ex['star_mapq']} CIGAR={ex['star_cigar']}")

    # Same-chr close disagreement examples
    if disagree_same_chr_close_examples:
        print(f"\nSame-chr close disagreement examples (6-50bp, first 10):")
        for ex in disagree_same_chr_close_examples[:10]:
            print(f"  {ex['qname'][:30]:<32} diff={ex['pos_diff']}bp")
            print(f"  {'':>32} rustar-aligner={ex['rustar_aligner_chr']}:{ex['rustar_aligner_pos']}({ex['rustar_aligner_strand']}) CIGAR={ex['rustar_aligner_cigar']}")
            print(f"  {'':>32} STAR  ={ex['star_chr']}:{ex['star_pos']}({ex['star_strand']}) CIGAR={ex['star_cigar']}")

    # Same-chr far disagreement examples
    if disagree_same_chr_far_examples:
        print(f"\nSame-chr far disagreement examples (51+bp, first 10):")
        for ex in disagree_same_chr_far_examples[:10]:
            print(f"  {ex['qname'][:30]:<32} diff={ex['pos_diff']}bp")
            print(f"  {'':>32} rustar-aligner={ex['rustar_aligner_chr']}:{ex['rustar_aligner_pos']}({ex['rustar_aligner_strand']}) CIGAR={ex['rustar_aligner_cigar']}")
            print(f"  {'':>32} STAR  ={ex['star_chr']}:{ex['star_pos']}({ex['star_strand']}) CIGAR={ex['star_cigar']}")

    if total_both_mapped > 0:
        concordance = 100.0 * both_mapped_agree_pos / total_both_mapped
        print(f"\nPosition concordance (among both-mapped): {concordance:.1f}%")

    # ============================================================
    # 3. JUNCTION COMPARISON
    # ============================================================
    if has_sj:
        print("\n" + "=" * 80)
        print("3. JUNCTION COMPARISON (SJ.out.tab)")
        print("=" * 80)

        rustar_aligner_sj = parse_sj(rustar_aligner_sj_path)
        star_sj = parse_sj(star_sj_path)

        rustar_aligner_keys = set(rustar_aligner_sj.keys())
        star_keys = set(star_sj.keys())

        shared = rustar_aligner_keys & star_keys
        rustar_aligner_only_sj = rustar_aligner_keys - star_keys
        star_only_sj = star_keys - rustar_aligner_keys

        print(f"\n{'Category':<30} {'Count':>8}")
        print("-" * 40)
        print(f"{'Shared junctions':<30} {len(shared):>8}")
        print(f"{'STAR-only junctions':<30} {len(star_only_sj):>8}")
        print(f"{'rustar-aligner-only junctions':<30} {len(rustar_aligner_only_sj):>8}")
        print(f"{'Total STAR junctions':<30} {len(star_keys):>8}")
        print(f"{'Total rustar-aligner junctions':<30} {len(rustar_aligner_keys):>8}")

        # Motif comparison for shared junctions
        motif_names = {0: "non-canonical", 1: "GT/AG", 2: "CT/AC", 3: "GC/AG", 4: "CT/GC", 5: "AT/AC", 6: "GT/AT"}
        motif_agree = 0
        motif_disagree = 0
        motif_disagree_examples = []

        for key in shared:
            r_motif = rustar_aligner_sj[key]["motif"]
            s_motif = star_sj[key]["motif"]
            if r_motif == s_motif:
                motif_agree += 1
            else:
                motif_disagree += 1
                if len(motif_disagree_examples) < 5:
                    motif_disagree_examples.append((key, r_motif, s_motif))

        if shared:
            print(f"\nMotif agreement for shared junctions: {motif_agree}/{len(shared)} ({100.0*motif_agree/len(shared):.1f}%)")
            if motif_disagree_examples:
                print("\nMotif disagreements (first 5):")
                for key, rm, sm in motif_disagree_examples:
                    print(f"  {key[0]}:{key[1]}-{key[2]}  rustar-aligner={motif_names.get(rm, str(rm))}  STAR={motif_names.get(sm, str(sm))}")

        # Coverage comparison for shared junctions
        if shared:
            print(f"\nCoverage comparison for shared junctions (top 20 by STAR unique count):")
            print(f"  {'Junction':<30} {'rSTAR_uniq':>11} {'STAR_uniq':>11} {'rSTAR_multi':>12} {'STAR_multi':>11}")
            print("  " + "-" * 78)
            shared_sorted = sorted(shared, key=lambda k: star_sj[k]["uniq_count"], reverse=True)
            for key in shared_sorted[:20]:
                r = rustar_aligner_sj[key]
                s = star_sj[key]
                jstr = f"{key[0]}:{key[1]}-{key[2]}"
                print(f"  {jstr:<30} {r['uniq_count']:>11} {s['uniq_count']:>11} {r['multi_count']:>12} {s['multi_count']:>11}")

        # Show STAR-only junctions
        if star_only_sj:
            print(f"\nSTAR-only junctions ({len(star_only_sj)}):")
            star_only_sorted = sorted(star_only_sj, key=lambda k: star_sj[k]["uniq_count"], reverse=True)
            for key in star_only_sorted[:10]:
                s = star_sj[key]
                jstr = f"{key[0]}:{key[1]}-{key[2]}"
                print(f"  {jstr:<30} motif={motif_names.get(s['motif'], str(s['motif'])):<15} uniq={s['uniq_count']:<5} multi={s['multi_count']:<5} annot={s['annotated']}")

        # Show rustar-aligner-only junctions (first 20 by count)
        if rustar_aligner_only_sj:
            print(f"\nrustar-aligner-only junctions (top 20 of {len(rustar_aligner_only_sj)} by unique count):")
            rustar_aligner_only_sorted = sorted(rustar_aligner_only_sj, key=lambda k: rustar_aligner_sj[k]["uniq_count"], reverse=True)
            for key in rustar_aligner_only_sorted[:20]:
                r = rustar_aligner_sj[key]
                jstr = f"{key[0]}:{key[1]}-{key[2]}"
                print(f"  {jstr:<30} motif={motif_names.get(r['motif'], str(r['motif'])):<15} uniq={r['uniq_count']:<5} multi={r['multi_count']:<5} annot={r['annotated']}")

        # Summarize rustar-aligner-only by motif
        if rustar_aligner_only_sj:
            motif_counter = Counter()
            for key in rustar_aligner_only_sj:
                motif_counter[rustar_aligner_sj[key]["motif"]] += 1
            print(f"\nrustar-aligner-only junctions by motif:")
            for motif_id, count in motif_counter.most_common():
                print(f"  {motif_names.get(motif_id, f'unknown({motif_id})'):<20} {count:>5}")
    else:
        print("\n(Skipping junction comparison - SJ.out.tab not found in both directories)")

    # ============================================================
    # 4. CIGAR PATTERN DISTRIBUTION
    # ============================================================
    print("\n" + "=" * 80)
    print("4. CIGAR PATTERN DISTRIBUTION")
    print("=" * 80)

    for label, reads in [("rustar-aligner", rustar_aligner_reads), ("STAR", star_reads)]:
        cats = Counter()
        total_mapped = 0
        for qname, records in reads.items():
            pri = get_primary(records)
            if pri:
                total_mapped += 1
                for cat in cigar_category(pri["cigar"]):
                    cats[cat] += 1

        print(f"\n{label} (mapped reads: {total_mapped}):")
        print(f"  {'Category':<20} {'Count':>8} {'%':>8}")
        print("  " + "-" * 40)
        for cat in ["pure_match", "spliced", "indel", "softclip", "unmapped"]:
            c = cats.get(cat, 0)
            pct = 100.0 * c / total_mapped if total_mapped else 0
            print(f"  {cat:<20} {c:>8} {pct:>7.1f}%")
        print(f"  (Note: categories overlap - a read can be spliced + softclipped)")

    # ============================================================
    # 5. EXAMPLE DISAGREEMENTS
    # ============================================================
    print("\n" + "=" * 80)
    print("5. EXAMPLE READS WITH POSITION DISAGREEMENT (first 10)")
    print("=" * 80)

    if not disagree_examples:
        print("\nNo position disagreements found!")
    else:
        for i, ex in enumerate(disagree_examples[:10]):
            print(f"\n--- Read {i+1}: {ex['qname']} ---")
            print(f"  rustar-aligner: {ex['rustar_aligner_chr']}:{ex['rustar_aligner_pos']} ({ex['rustar_aligner_strand']}) MAPQ={ex['rustar_aligner_mapq']} CIGAR={ex['rustar_aligner_cigar']}")
            print(f"  STAR:   {ex['star_chr']}:{ex['star_pos']} ({ex['star_strand']}) MAPQ={ex['star_mapq']} CIGAR={ex['star_cigar']}")

            if ex['rustar_aligner_chr'] != ex['star_chr']:
                print(f"  >> Different chromosome!")
            else:
                diff = abs(ex['rustar_aligner_pos'] - ex['star_pos'])
                print(f"  >> Same chr, position difference: {diff}bp")

    # ============================================================
    # 6. MAPQ AGREEMENT
    # ============================================================
    print("\n" + "=" * 80)
    print("6. MAPQ AGREEMENT")
    print("=" * 80)

    mapq_agree = 0
    mapq_disagree = 0
    mapq_inflation = 0  # rustar-aligner=255, STAR<255
    mapq_deflation = 0  # rustar-aligner<255, STAR=255

    for qname in sorted(all_reads):
        r_records = rustar_aligner_reads.get(qname, [])
        s_records = star_reads.get(qname, [])
        r_class = classify_read(r_records) if r_records else "missing"
        s_class = classify_read(s_records) if s_records else "missing"
        if r_class in ("unique", "multi") and s_class in ("unique", "multi"):
            r_pri = get_primary(r_records)
            s_pri = get_primary(s_records)
            if r_pri and s_pri:
                # Only compare MAPQ for reads that agree on position
                same_chr = r_pri["rname"] == s_pri["rname"]
                same_pos = same_chr and abs(r_pri["pos"] - s_pri["pos"]) <= 5
                if same_pos:
                    if r_pri["mapq"] == s_pri["mapq"]:
                        mapq_agree += 1
                    else:
                        mapq_disagree += 1
                        if r_pri["mapq"] == 255 and s_pri["mapq"] < 255:
                            mapq_inflation += 1
                        elif r_pri["mapq"] < 255 and s_pri["mapq"] == 255:
                            mapq_deflation += 1

    total_mapq_compared = mapq_agree + mapq_disagree
    if total_mapq_compared > 0:
        print(f"\nMAPQ agreement (pos-agreeing reads): {mapq_agree}/{total_mapq_compared} ({100.0*mapq_agree/total_mapq_compared:.1f}%)")
        print(f"  MAPQ disagree:     {mapq_disagree}")
        print(f"  Inflation (r=255, S<255): {mapq_inflation}")
        print(f"  Deflation (r<255, S=255): {mapq_deflation}")

    # ============================================================
    # 7. ADJUSTED AGREEMENT (excluding diff-chr multi-mapper ties)
    # ============================================================
    print("\n" + "=" * 80)
    print("7. ADJUSTED AGREEMENT (excluding diff-chr multi-mapper ties)")
    print("=" * 80)

    # Count diff-chr multi-mapper ties: both mapped to different chromosomes,
    # both with MAPQ < 255 (multi-mappers where primary locus is a tie-break)
    diff_chr_multimap_ties = 0
    for qname in sorted(all_reads):
        r_records = rustar_aligner_reads.get(qname, [])
        s_records = star_reads.get(qname, [])
        r_class = classify_read(r_records) if r_records else "missing"
        s_class = classify_read(s_records) if s_records else "missing"
        if r_class in ("unique", "multi") and s_class in ("unique", "multi"):
            r_pri = get_primary(r_records)
            s_pri = get_primary(s_records)
            if r_pri and s_pri and r_pri["rname"] != s_pri["rname"]:
                if r_pri["mapq"] < 255 and s_pri["mapq"] < 255:
                    diff_chr_multimap_ties += 1

    # Adjusted = exclude diff-chr multi-mapper ties from disagreement count
    adjusted_total = total_both_mapped - diff_chr_multimap_ties
    adjusted_agree = both_mapped_agree_pos
    actionable_same_chr = both_mapped_disagree_pos - disagree_diff_chr
    # Some diff-chr disagreements involve unique mappers (not ties) — count those too
    actionable_diff_chr = disagree_diff_chr - diff_chr_multimap_ties

    print(f"\nDiff-chr multi-mapper ties:    {diff_chr_multimap_ties}  (unavoidable tie-breaking)")
    if adjusted_total > 0:
        print(f"Adjusted position agreement:   {adjusted_agree}/{adjusted_total} ({100.0*adjusted_agree/adjusted_total:.1f}%)")
    print(f"Actionable disagreements:      {actionable_same_chr} same-chr + {actionable_diff_chr} diff-chr-unique + {star_only_mapped} STAR-only + {rustar_aligner_only_mapped} rustar-aligner-only = {actionable_same_chr + actionable_diff_chr + star_only_mapped + rustar_aligner_only_mapped}")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    total = len(all_reads)
    print(f"""
Total reads:               {total}
Both mapped, agree:        {both_mapped_agree_pos} ({100.0*both_mapped_agree_pos/total:.1f}%)
Both mapped, disagree:     {both_mapped_disagree_pos} ({100.0*both_mapped_disagree_pos/total:.1f}%)
STAR only mapped:          {star_only_mapped} ({100.0*star_only_mapped/total:.1f}%)
rustar-aligner only mapped:        {rustar_aligner_only_mapped} ({100.0*rustar_aligner_only_mapped/total:.1f}%)
Both unmapped:             {both_unmapped} ({100.0*both_unmapped/total:.1f}%)""")

    if has_sj:
        print(f"""
Junctions shared:          {len(shared)}/{len(star_keys)} STAR junctions ({100.0*len(shared)/max(len(star_keys),1):.1f}%)
rustar-aligner extra junctions:    {len(rustar_aligner_only_sj)} (potential false positives)""")
    print()


if __name__ == "__main__":
    main()
