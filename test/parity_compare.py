#!/usr/bin/env python3
"""Compare a STAR output directory with a rustar-aligner one (used by parity_suite.sh).

  parity_compare.py index STAR_DIR RUSTAR_DIR NAME   genome index files, byte for byte
  parity_compare.py align STAR_DIR RUSTAR_DIR NAME   every alignment output

Alignment records are compared as a multiset per read (samtools view), because with more
than one thread STAR writes reads in a thread-dependent order. Tables are compared line by
line: IDENTICAL (same bytes), SAME-SET (same lines, different order) or DIFFER.
"""
import collections
import filecmp
import gzip
import os
import subprocess
import sys

# Differ by design: command line, version strings, paths.
INDEX_SKIP = {"Log.out", "genomeParameters.txt"}
LOG_SKIP = ("Started job on", "Started mapping on", "Finished on", "Mapping speed")
BAMS = ["Aligned.out.sam", "Aligned.out.bam", "Aligned.sortedByCoord.out.bam",
        "Aligned.toTranscriptome.out.bam"]
TABLES = ["SJ.out.tab", "ReadsPerGene.out.tab", "Chimeric.out.junction"]


def timing(path):
    try:
        wall, rss = open(path).read().split()[-2:]
        return f"{float(wall):.1f} s, {int(rss) / 1024:.0f} MB"
    except (OSError, ValueError):
        return "n/a"


def records(path, ignore_primary=False):
    out = subprocess.run(["samtools", "view", path], capture_output=True, text=True,
                         check=True).stdout
    by_read = collections.defaultdict(collections.Counter)
    for line in out.splitlines():
        f = line.split("\t")
        if ignore_primary:
            f[1] = str(int(f[1]) & ~0x100)
        by_read[f[0]]["\t".join(f)] += 1
    return by_read


def compare_records(a, b, ignore_primary=False):
    ra, rb = records(a, ignore_primary), records(b, ignore_primary)
    n_a = sum(sum(c.values()) for c in ra.values())
    n_b = sum(sum(c.values()) for c in rb.values())
    same = sum(sum((ra[q] & rb[q]).values()) for q in ra.keys() & rb.keys())
    reads = ra.keys() | rb.keys()
    differ = sorted(q for q in reads if ra.get(q) != rb.get(q))
    return (f"records {n_a}/{n_b}, identical {same}; reads {len(reads)}, "
            f"differing {len(differ)}" + (f" (e.g. {', '.join(differ[:3])})" if differ else ""))


def lines(path):
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as f:
        return f.read().splitlines()


def compare_lines(a, b):
    la, lb = lines(a), lines(b)
    if la == lb:
        verdict = "IDENTICAL"
    elif collections.Counter(la) == collections.Counter(lb):
        verdict = "SAME-SET"
    else:
        verdict = "DIFFER"
    same = sum((collections.Counter(la) & collections.Counter(lb)).values())
    return f"{verdict} (lines {len(la)}/{len(lb)}, shared {same})"


def walk(d):
    return {os.path.relpath(os.path.join(r, f), d)
            for r, _, fs in os.walk(d) for f in fs}


def index(a, b, name):
    files = sorted(f for f in walk(a) if f not in INDEX_SKIP)
    differ = [f for f in files
              if not os.path.exists(os.path.join(b, f))
              or not filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False)]
    print(f"== index {name}: {len(files) - len(differ)}/{len(files)} files identical"
          + (f", differ: {' '.join(differ)}" if differ else ""))
    print(f"   STAR {timing(a + '.time')}; rustar-aligner {timing(b + '.time')}")


def align(a, b, name):
    print(f"== {name}")
    print(f"   STAR {timing(os.path.join(a, 'time.txt'))}; "
          f"rustar-aligner {timing(os.path.join(b, 'time.txt'))}")
    for f in BAMS:
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if not os.path.exists(pa):
            continue
        if not os.path.exists(pb):
            print(f"   {f}: missing from rustar-aligner")
            continue
        print(f"   {f}: {compare_records(pa, pb)}")
        if f == "Aligned.toTranscriptome.out.bam":
            # STAR draws the primary from a per-thread RNG stream (DIVERGENCE.md 1.1).
            print(f"   {f} ignoring the 0x100 flag: "
                  f"{compare_records(pa, pb, ignore_primary=True)}")
    for f in TABLES:
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if os.path.exists(pa):
            print(f"   {f}: {compare_lines(pa, pb) if os.path.exists(pb) else 'missing'}")
    la = [x for x in lines(os.path.join(a, "Log.final.out")) if not x.strip().startswith(LOG_SKIP)]
    lb = [x for x in lines(os.path.join(b, "Log.final.out")) if not x.strip().startswith(LOG_SKIP)]
    diff = [f"{x.split('|')[0].strip()}: {x.split('|')[1].strip()} vs {y.split('|')[1].strip()}"
            for x, y in zip(la, lb) if x != y and "|" in x and "|" in y]
    print("   Log.final.out: " + ("identical" if la == lb else "; ".join(diff) or "differs"))
    sa = os.path.join(a, "Solo.out")
    if os.path.isdir(sa):
        sb = os.path.join(b, "Solo.out")
        fa, fb = walk(sa), walk(sb) if os.path.isdir(sb) else set()
        differ = sorted(f for f in fa & fb if not filecmp.cmp(
            os.path.join(sa, f), os.path.join(sb, f), shallow=False))
        same = len(fa & fb) - len(differ)
        print(f"   Solo.out: {same}/{len(fa)} files identical"
              + (f"; differ: {' '.join(differ)}" if differ else "")
              + (f"; only STAR: {' '.join(sorted(fa - fb))}" if fa - fb else "")
              + (f"; only rustar-aligner: {' '.join(sorted(fb - fa))}" if fb - fa else ""))


if __name__ == "__main__":
    if len(sys.argv) != 5 or sys.argv[1] not in ("index", "align"):
        sys.exit(__doc__)
    {"index": index, "align": align}[sys.argv[1]](*sys.argv[2:])
