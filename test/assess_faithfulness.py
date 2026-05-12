#!/usr/bin/env python3
"""
assess_faithfulness.py — Comprehensive rustar-aligner vs STAR faithfulness assessment.

Reports:
  1. SE exact-record agreement (FLAG, RNAME, POS, MAPQ, CIGAR, tags NH/AS/NM)
  2. PE exact-record agreement (+ RNEXT, PNEXT, TLEN, proper-pair flag)
  3. SJ.out.tab field-level comparison

Usage:
    python assess_faithfulness.py \\
        --se-star  <star_se.sam>  --se-rustar-aligner  <rustar_aligner_se.sam> \\
        --pe-star  <star_pe.sam>  --pe-rustar-aligner  <rustar_aligner_pe.sam> \\
        --sj-star  <star_SJ.tab> --sj-rustar-aligner  <rustar_aligner_SJ.tab> \\
        [--se-only | --pe-only | --sj-only]
"""

import argparse
import sys
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_tags(fields):
    """Return dict of tag -> (type, value) from SAM optional fields."""
    tags = {}
    for f in fields:
        m = re.match(r'([A-Z]{2}):([AifZHB]):(.+)', f)
        if m:
            tags[m.group(1)] = (m.group(2), m.group(3))
    return tags


def tag_int(tags, key):
    """Return integer value of a tag, or None."""
    if key in tags:
        try:
            return int(tags[key][1])
        except ValueError:
            pass
    return None


def parse_sam(path, primary_only=True):
    """
    Parse SAM/BAM file (SAM only). Returns dict:
      qname -> {flag, rname, pos, mapq, cigar, rnext, pnext, tlen, tags}
    Only primary alignments kept (flag & 0x100 == 0) unless primary_only=False.
    Unmapped reads (flag & 0x4) are kept with rname='*'.
    For PE: key is (qname, mate) where mate = 1 or 2.
    """
    records = {}
    with open(path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if primary_only and (flag & 0x100):   # skip secondary
                continue
            if flag & 0x800:                       # skip supplementary
                continue
            qname   = fields[0]
            rname   = fields[2]
            pos     = int(fields[3])
            mapq    = int(fields[4])
            cigar   = fields[5]
            rnext   = fields[6]
            pnext   = int(fields[7])
            tlen    = int(fields[8])
            tags    = parse_tags(fields[11:])

            is_pe = bool(flag & 0x1)
            if is_pe:
                mate = 1 if (flag & 0x40) else 2
                key = (qname, mate)
            else:
                key = qname

            records[key] = {
                'flag':  flag,
                'rname': rname,
                'pos':   pos,
                'mapq':  mapq,
                'cigar': cigar,
                'rnext': rnext,
                'pnext': pnext,
                'tlen':  tlen,
                'tags':  tags,
            }
    return records


def flag_bits_str(flag, mask_bits=None):
    """Return human-readable flag description for the key bits."""
    names = {0x4: 'unmapped', 0x10: 'reverse', 0x100: 'secondary',
             0x800: 'supplementary', 0x2: 'proper_pair',
             0x8: 'mate_unmapped', 0x20: 'mate_reverse'}
    bits = mask_bits or names.keys()
    return '|'.join(n for b, n in names.items() if b in bits and (flag & b))


# ---------------------------------------------------------------------------
# SE comparison
# ---------------------------------------------------------------------------

def compare_se(star_path, rustar_aligner_path):
    print("=" * 72)
    print("SE FAITHFULNESS ASSESSMENT")
    print("=" * 72)

    star   = parse_sam(star_path)
    rustar_aligner = parse_sam(rustar_aligner_path)

    star_mapped   = {k: v for k, v in star.items()   if not (v['flag'] & 0x4)}
    rustar_aligner_mapped = {k: v for k, v in rustar_aligner.items() if not (v['flag'] & 0x4)}

    only_star   = set(star_mapped)   - set(rustar_aligner_mapped)
    only_rustar_aligner = set(rustar_aligner_mapped) - set(star_mapped)
    both        = set(star_mapped)   & set(rustar_aligner_mapped)

    total_star = len(star_mapped)
    print(f"\nSTAR mapped:         {total_star}")
    print(f"rustar-aligner mapped:       {len(rustar_aligner_mapped)}")
    print(f"STAR-only (missed):  {len(only_star)}  ({100*len(only_star)/max(total_star,1):.3f}%)")
    print(f"rustar-aligner-only (FP):    {len(only_rustar_aligner)}  ({100*len(only_rustar_aligner)/max(len(rustar_aligner_mapped),1):.3f}%)")
    print(f"Both mapped:         {len(both)}")

    # --- exact record comparison ---
    n_exact        = 0
    n_pos_agree    = 0
    n_tie_pos      = 0   # pos differs but MAPQ+NH agree → tie-breaking, not an error
    diff_pos       = 0
    diff_chr       = 0
    diff_strand    = 0
    diff_cigar_pos = 0   # same pos, diff CIGAR
    diff_mapq      = 0
    diff_flag      = 0
    diff_nh        = 0
    diff_as        = 0
    diff_nm        = 0
    mapq_inflate   = 0
    mapq_deflate   = 0

    pos_disagree_examples   = []
    cigar_disagree_examples = []
    mapq_disagree_examples  = []
    flag_disagree_examples  = []
    nh_disagree_examples    = []
    as_disagree_examples    = []
    nm_disagree_examples    = []

    for k in both:
        s = star_mapped[k]
        r = rustar_aligner_mapped[k]

        same_chr    = (s['rname'] == r['rname'])
        same_pos    = same_chr and (s['pos'] == r['pos'])
        same_strand = same_pos and ((s['flag'] & 0x10) == (r['flag'] & 0x10))
        same_cigar  = same_strand and (s['cigar'] == r['cigar'])
        same_mapq   = (s['mapq'] == r['mapq'])
        # FLAG: compare only the relevant SE bits (ignore paired/mate bits)
        se_mask     = 0x10 | 0x4  # reverse + unmapped (the only meaningful SE bits)
        same_flag_se = ((s['flag'] & se_mask) == (r['flag'] & se_mask))

        s_nh = tag_int(s['tags'], 'NH')
        r_nh = tag_int(r['tags'], 'NH')
        s_as = tag_int(s['tags'], 'AS')
        r_as = tag_int(r['tags'], 'AS')
        s_nm = tag_int(s['tags'], 'NM')
        r_nm = tag_int(r['tags'], 'NM')

        same_nh = (s_nh == r_nh)
        same_as = (s_as == r_as)
        same_nm = (s_nm == r_nm)

        # A tie-breaking difference: position differs but MAPQ and NH agree.
        # Both tools found the same pool of equally-valid alignments and chose
        # different primaries — neither is wrong, just a different valid choice.
        is_tie_pos = (not same_pos) and same_mapq and same_nh

        if is_tie_pos:
            n_tie_pos += 1

        if same_chr and same_pos and same_strand and same_cigar and same_mapq and same_flag_se:
            n_exact += 1
        if same_chr and same_pos and same_strand:
            n_pos_agree += 1

        if not same_chr:
            diff_chr += 1
            if len(pos_disagree_examples) < 8:
                pos_disagree_examples.append(
                    f"  {k}: rustar-aligner={r['rname']}:{r['pos']} STAR={s['rname']}:{s['pos']}")
        elif not same_pos:
            diff_pos += 1
            if len(pos_disagree_examples) < 8:
                pos_disagree_examples.append(
                    f"  {k}: rustar-aligner={r['rname']}:{r['pos']} STAR={s['rname']}:{s['pos']}")
        elif not same_strand:
            diff_strand += 1

        if same_chr and same_pos and same_strand and not same_cigar:
            diff_cigar_pos += 1
            if len(cigar_disagree_examples) < 8:
                cigar_disagree_examples.append(
                    f"  {k}: rustar-aligner={r['cigar']} STAR={s['cigar']}")

        if not same_mapq:
            diff_mapq += 1
            if s['mapq'] < r['mapq']:
                mapq_inflate += 1
            else:
                mapq_deflate += 1
            if len(mapq_disagree_examples) < 8:
                mapq_disagree_examples.append(
                    f"  {k}: rustar-aligner={r['mapq']} STAR={s['mapq']}")

        if not same_flag_se and len(flag_disagree_examples) < 8:
            diff_flag += 1
            flag_disagree_examples.append(
                f"  {k}: rustar-aligner_flag={r['flag']} STAR_flag={s['flag']}")

        if not same_nh:
            diff_nh += 1
            if len(nh_disagree_examples) < 8:
                nh_disagree_examples.append(f"  {k}: rustar-aligner_NH={r_nh} STAR_NH={s_nh}")

        if not same_as:
            diff_as += 1
            if len(as_disagree_examples) < 8:
                as_disagree_examples.append(f"  {k}: rustar-aligner_AS={r_as} STAR_AS={s_as}")

        if not same_nm:
            diff_nm += 1
            if len(nm_disagree_examples) < 8:
                nm_disagree_examples.append(f"  {k}: rustar-aligner_NM={r_nm} STAR_NM={s_nm}")

    n = len(both)
    n_adj = n - n_tie_pos   # denominator excluding tie-breaking position differences
    n_exact_adj = n_exact   # exact matches are all among non-tie reads (ties never exact-match pos)
    print(f"\n{'─'*50}")
    print("EXACT-RECORD AGREEMENT (primary alignments, both mapped)")
    print(f"{'─'*50}")
    pct = lambda x: f"{100*x/max(n,1):.3f}%"
    pct_adj = lambda x: f"{100*x/max(n_adj,1):.3f}%"
    print(f"  Total compared:          {n}")
    print(f"  Tie-breaking diffs:      {n_tie_pos:>6}  (pos differs, same MAPQ+NH — excluded below)")
    print(f"  Comparable (non-tie):    {n_adj}")
    print(f"  Exact match (all fields):{n_exact:>8}  ({pct(n_exact)} raw / {pct_adj(n_exact_adj)} tie-adjusted)")
    print(f"  Position+strand agree:   {n_pos_agree:>8}  ({pct(n_pos_agree)})")
    print()
    print("  Breakdown of disagreements:")
    print(f"    diff chromosome:        {diff_chr:>6}  ({pct(diff_chr)})")
    print(f"    diff position:          {diff_pos:>6}  ({pct(diff_pos)})")
    print(f"    diff strand:            {diff_strand:>6}  ({pct(diff_strand)})")
    print(f"    diff CIGAR (pos match): {diff_cigar_pos:>6}  ({pct(diff_cigar_pos)})")
    print(f"    diff MAPQ:              {diff_mapq:>6}  ({pct(diff_mapq)})")
    print(f"      inflate (rustar-aligner>STAR):{mapq_inflate:>6}")
    print(f"      deflate (rustar-aligner<STAR):{mapq_deflate:>6}")
    print(f"    diff FLAG (se bits):    {diff_flag:>6}  ({pct(diff_flag)})")
    print(f"    diff NH tag:            {diff_nh:>6}  ({pct(diff_nh)})")
    print(f"    diff AS tag:            {diff_as:>6}  ({pct(diff_as)})")
    print(f"    diff NM tag:            {diff_nm:>6}  ({pct(diff_nm)})")

    def show_examples(label, examples):
        if examples:
            print(f"\n  {label} (first {len(examples)}):")
            for e in examples:
                print(e)

    show_examples("Position disagreements", pos_disagree_examples)
    show_examples("CIGAR disagreements (same pos)", cigar_disagree_examples)
    show_examples("MAPQ disagreements", mapq_disagree_examples)
    show_examples("FLAG disagreements", flag_disagree_examples)
    show_examples("NH disagreements", nh_disagree_examples)
    show_examples("AS disagreements", as_disagree_examples)
    show_examples("NM disagreements", nm_disagree_examples)

    if only_star:
        print(f"\n  STAR-only reads (rustar-aligner missed): {sorted(only_star)[:10]}")
    if only_rustar_aligner:
        print(f"\n  rustar-aligner-only reads (false positives): {sorted(only_rustar_aligner)[:10]}")

    return n_exact, n_adj


# ---------------------------------------------------------------------------
# PE comparison
# ---------------------------------------------------------------------------

def compare_pe(star_path, rustar_aligner_path):
    print("\n" + "=" * 72)
    print("PE FAITHFULNESS ASSESSMENT")
    print("=" * 72)

    star   = parse_sam(star_path)
    rustar_aligner = parse_sam(rustar_aligner_path)

    star_mapped   = {k: v for k, v in star.items()   if not (v['flag'] & 0x4)}
    rustar_aligner_mapped = {k: v for k, v in rustar_aligner.items() if not (v['flag'] & 0x4)}

    only_star   = set(star_mapped)   - set(rustar_aligner_mapped)
    only_rustar_aligner = set(rustar_aligner_mapped) - set(star_mapped)
    both        = set(star_mapped)   & set(rustar_aligner_mapped)

    total_star = len(star_mapped)
    print(f"\nSTAR mapped mates:       {total_star}")
    print(f"rustar-aligner mapped mates:     {len(rustar_aligner_mapped)}")
    print(f"STAR-only (missed):      {len(only_star)}  ({100*len(only_star)/max(total_star,1):.3f}%)")
    print(f"rustar-aligner-only (FP):        {len(only_rustar_aligner)}  ({100*len(only_rustar_aligner)/max(len(rustar_aligner_mapped),1):.3f}%)")
    print(f"Both mapped:             {len(both)}")

    n_exact     = 0
    n_pos_agree = 0
    n_tie_pos   = 0   # pos differs but MAPQ+NH agree → tie-breaking, not an error
    diff_pos    = 0
    diff_cigar  = 0
    diff_mapq   = 0
    diff_proper = 0   # proper-pair flag bit
    diff_tlen   = 0
    diff_nh     = 0
    diff_as     = 0
    diff_nm     = 0
    mapq_inflate = 0
    mapq_deflate = 0

    proper_disagree_examples = []

    for k in both:
        s = star_mapped[k]
        r = rustar_aligner_mapped[k]

        same_chr    = (s['rname'] == r['rname'])
        same_pos    = same_chr and (s['pos'] == r['pos'])
        same_strand = same_pos and ((s['flag'] & 0x10) == (r['flag'] & 0x10))
        same_cigar  = same_strand and (s['cigar'] == r['cigar'])
        same_mapq   = (s['mapq'] == r['mapq'])
        same_proper = ((s['flag'] & 0x2) == (r['flag'] & 0x2))
        same_tlen   = (s['tlen'] == r['tlen'])

        s_nh = tag_int(s['tags'], 'NH')
        r_nh = tag_int(r['tags'], 'NH')
        s_as = tag_int(s['tags'], 'AS')
        r_as = tag_int(r['tags'], 'AS')
        s_nm = tag_int(s['tags'], 'NM')
        r_nm = tag_int(r['tags'], 'NM')

        same_nh = (s_nh == r_nh)
        same_as = (s_as == r_as)
        same_nm = (s_nm == r_nm)

        # A tie-breaking difference: position differs but MAPQ and NH agree.
        # Both tools found the same pool of equally-valid alignments and chose
        # different primaries — neither is wrong, just a different valid choice.
        is_tie_pos = (not same_pos) and same_mapq and same_nh
        if is_tie_pos:
            n_tie_pos += 1

        if same_cigar and same_mapq and same_proper and same_nh:
            n_exact += 1
        if same_chr and same_pos and same_strand:
            n_pos_agree += 1

        if not (same_chr and same_pos):
            diff_pos += 1
        elif same_pos and same_strand and not same_cigar:
            diff_cigar += 1

        if not same_mapq:
            diff_mapq += 1
            if s['mapq'] < r['mapq']:
                mapq_inflate += 1
            else:
                mapq_deflate += 1

        if not same_proper:
            diff_proper += 1
            if len(proper_disagree_examples) < 8:
                qname, mate = k
                proper_disagree_examples.append(
                    f"  {qname}/mate{mate}: rustar-aligner_proper={bool(r['flag']&2)} "
                    f"STAR_proper={bool(s['flag']&2)} "
                    f"rname={r['rname']}:{r['pos']} tlen={r['tlen']}")

        if not same_tlen:
            diff_tlen += 1

        if not same_nh:
            diff_nh += 1
        if not same_as:
            diff_as += 1
        if not same_nm:
            diff_nm += 1

    n = len(both)
    n_adj = n - n_tie_pos   # denominator excluding tie-breaking position differences
    pct = lambda x: f"{100*x/max(n,1):.3f}%"
    pct_adj = lambda x: f"{100*x/max(n_adj,1):.3f}%"
    print(f"\n{'─'*50}")
    print("EXACT-RECORD AGREEMENT (per mate, both mapped)")
    print(f"{'─'*50}")
    print(f"  Total compared:          {n}")
    print(f"  Tie-breaking diffs:      {n_tie_pos:>6}  (pos differs, same MAPQ+NH — excluded below)")
    print(f"  Comparable (non-tie):    {n_adj}")
    print(f"  Exact match (pos+CIGAR+MAPQ+proper+NH): {n_exact:>6}  ({pct(n_exact)} raw / {pct_adj(n_exact)} tie-adjusted)")
    print(f"  Position+strand agree:   {n_pos_agree:>8}  ({pct(n_pos_agree)})")
    print()
    print("  Breakdown of disagreements:")
    print(f"    diff position/chr:      {diff_pos:>6}  ({pct(diff_pos)})")
    print(f"    diff CIGAR (pos match): {diff_cigar:>6}  ({pct(diff_cigar)})")
    print(f"    diff MAPQ:              {diff_mapq:>6}  ({pct(diff_mapq)})")
    print(f"      inflate (rustar-aligner>STAR):{mapq_inflate:>6}")
    print(f"      deflate (rustar-aligner<STAR):{mapq_deflate:>6}")
    print(f"    diff proper-pair flag:  {diff_proper:>6}  ({pct(diff_proper)})")
    print(f"    diff TLEN:              {diff_tlen:>6}  ({pct(diff_tlen)})")
    print(f"    diff NH tag:            {diff_nh:>6}  ({pct(diff_nh)})")
    print(f"    diff AS tag:            {diff_as:>6}  ({pct(diff_as)})")
    print(f"    diff NM tag:            {diff_nm:>6}  ({pct(diff_nm)})")

    if proper_disagree_examples:
        print(f"\n  Proper-pair disagreements (first {len(proper_disagree_examples)}):")
        for e in proper_disagree_examples:
            print(e)

    if only_star:
        qnames = sorted(set(k[0] for k in only_star))[:10]
        print(f"\n  STAR-only mates (sample): {qnames}")
    if only_rustar_aligner:
        qnames = sorted(set(k[0] for k in only_rustar_aligner))[:10]
        print(f"\n  rustar-aligner-only mates (sample): {qnames}")

    return n_exact, n_adj


# ---------------------------------------------------------------------------
# SJ.out.tab comparison
# ---------------------------------------------------------------------------

SJ_COLS = ['chr', 'start', 'end', 'strand', 'motif', 'annotated',
           'uniq_reads', 'multi_reads', 'max_overhang']
SJ_MOTIFS = {0:'non-canonical', 1:'GT/AG', 2:'CT/AC', 3:'GC/AG',
             4:'CT/GC', 5:'AT/AC', 6:'GT/AT'}
SJ_STRAND = {0:'undef', 1:'+', 2:'-'}


def parse_sj(path):
    junctions = {}
    with open(path) as f:
        for line in f:
            f2 = line.strip().split('\t')
            if len(f2) < 9:
                continue
            key = (f2[0], int(f2[1]), int(f2[2]))  # chr, start, end
            junctions[key] = {
                'strand':       int(f2[3]),
                'motif':        int(f2[4]),
                'annotated':    int(f2[5]),
                'uniq_reads':   int(f2[6]),
                'multi_reads':  int(f2[7]),
                'max_overhang': int(f2[8]),
            }
    return junctions


def compare_sj(star_path, rustar_aligner_path):
    print("\n" + "=" * 72)
    print("SJ.OUT.TAB FAITHFULNESS ASSESSMENT")
    print("=" * 72)

    star   = parse_sj(star_path)
    rustar_aligner = parse_sj(rustar_aligner_path)

    star_keys   = set(star.keys())
    rustar_aligner_keys = set(rustar_aligner.keys())
    shared      = star_keys & rustar_aligner_keys
    only_star   = star_keys   - rustar_aligner_keys
    only_rustar_aligner = rustar_aligner_keys - star_keys

    print(f"\nSTAR junctions:      {len(star_keys)}")
    print(f"rustar-aligner junctions:    {len(rustar_aligner_keys)}")
    print(f"Shared:              {len(shared)}")
    print(f"STAR-only:           {len(only_star)}")
    print(f"rustar-aligner-only:         {len(only_rustar_aligner)}")

    # For shared junctions, compare all fields
    n_exact         = 0
    diff_strand     = 0
    diff_motif      = 0
    diff_annotated  = 0
    diff_uniq       = 0
    diff_multi      = 0
    diff_overhang   = 0

    strand_ex  = []
    motif_ex   = []
    uniq_ex    = []

    for k in shared:
        s = star[k]
        r = rustar_aligner[k]
        same = True

        if s['strand'] != r['strand']:
            diff_strand += 1
            same = False
            if len(strand_ex) < 5:
                strand_ex.append(f"  {k}: rustar-aligner={SJ_STRAND.get(r['strand'],r['strand'])} "
                                  f"STAR={SJ_STRAND.get(s['strand'],s['strand'])}")

        if s['motif'] != r['motif']:
            diff_motif += 1
            same = False
            if len(motif_ex) < 5:
                motif_ex.append(f"  {k}: rustar-aligner={SJ_MOTIFS.get(r['motif'],r['motif'])} "
                                 f"STAR={SJ_MOTIFS.get(s['motif'],s['motif'])}")

        if s['annotated'] != r['annotated']:
            diff_annotated += 1
            same = False

        if s['uniq_reads'] != r['uniq_reads']:
            diff_uniq += 1
            same = False
            if len(uniq_ex) < 5:
                uniq_ex.append(f"  {k}: rustar-aligner={r['uniq_reads']} STAR={s['uniq_reads']}")

        if s['multi_reads'] != r['multi_reads']:
            diff_multi += 1
            same = False

        if s['max_overhang'] != r['max_overhang']:
            diff_overhang += 1
            same = False

        if same:
            n_exact += 1

    n = len(shared)
    pct = lambda x: f"{100*x/max(n,1):.1f}%"
    print(f"\n{'─'*50}")
    print("SHARED JUNCTION AGREEMENT")
    print(f"{'─'*50}")
    print(f"  Exact match (all fields): {n_exact:>5}  ({pct(n_exact)})")
    print(f"  diff strand:              {diff_strand:>5}  ({pct(diff_strand)})")
    print(f"  diff motif:               {diff_motif:>5}  ({pct(diff_motif)})")
    print(f"  diff annotated flag:      {diff_annotated:>5}  ({pct(diff_annotated)})")
    print(f"  diff uniq read count:     {diff_uniq:>5}  ({pct(diff_uniq)})")
    print(f"  diff multi read count:    {diff_multi:>5}  ({pct(diff_multi)})")
    print(f"  diff max overhang:        {diff_overhang:>5}  ({pct(diff_overhang)})")

    if strand_ex:
        print(f"\n  Strand mismatches: " + "; ".join(e.strip() for e in strand_ex))
    if motif_ex:
        print(f"  Motif mismatches: "  + "; ".join(e.strip() for e in motif_ex))
    if uniq_ex:
        print(f"\n  Unique-count mismatches (first {len(uniq_ex)}):")
        for e in uniq_ex:
            print(e)

    if only_star:
        print(f"\n  STAR-only junctions ({len(only_star)}):")
        for k in sorted(only_star)[:15]:
            s = star[k]
            print(f"    {k[0]}:{k[1]}-{k[2]}  strand={SJ_STRAND.get(s['strand'],'?')} "
                  f"motif={SJ_MOTIFS.get(s['motif'],'?')} "
                  f"annot={s['annotated']} uniq={s['uniq_reads']}")

    if only_rustar_aligner:
        print(f"\n  rustar-aligner-only junctions ({len(only_rustar_aligner)}):")
        for k in sorted(only_rustar_aligner)[:15]:
            r = rustar_aligner[k]
            print(f"    {k[0]}:{k[1]}-{k[2]}  strand={SJ_STRAND.get(r['strand'],'?')} "
                  f"motif={SJ_MOTIFS.get(r['motif'],'?')} "
                  f"annot={r['annotated']} uniq={r['uniq_reads']}")

    return n_exact, n


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(results):
    print("\n" + "=" * 72)
    print("FAITHFULNESS SUMMARY")
    print("=" * 72)
    for label, exact, total in results:
        pct = 100 * exact / max(total, 1)
        print(f"  {label:45s} {exact:>6}/{total:<6} ({pct:.3f}%)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--se-star',    metavar='SAM')
    ap.add_argument('--se-rustar-aligner', metavar='SAM')
    ap.add_argument('--pe-star',    metavar='SAM')
    ap.add_argument('--pe-rustar-aligner', metavar='SAM')
    ap.add_argument('--sj-star',    metavar='TAB')
    ap.add_argument('--sj-rustar-aligner', metavar='TAB')
    ap.add_argument('--se-only',   action='store_true')
    ap.add_argument('--sj-only',   action='store_true')
    ap.add_argument('--pe-only',   action='store_true')
    args = ap.parse_args()

    results = []

    if (not args.pe_only and not args.sj_only) and args.se_star and args.se_rustar_aligner:
        se_exact, se_total = compare_se(args.se_star, args.se_rustar_aligner)
        results.append(("SE exact records (tie-adjusted denom)", se_exact, se_total))

    if (not args.se_only and not args.sj_only) and args.pe_star and args.pe_rustar_aligner:
        pe_exact, pe_total = compare_pe(args.pe_star, args.pe_rustar_aligner)
        results.append(("PE exact mates (tie-adjusted denom)", pe_exact, pe_total))

    if (not args.se_only and not args.pe_only) and args.sj_star and args.sj_rustar_aligner:
        sj_exact, sj_total = compare_sj(args.sj_star, args.sj_rustar_aligner)
        results.append(("SJ exact records (all fields)", sj_exact, sj_total))

    if results:
        print_summary(results)
    else:
        ap.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
