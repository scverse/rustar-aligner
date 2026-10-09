#!/usr/bin/env python3
"""Differential check for --clipAdapterType CellRanger4 + --clip5pNbases.

Builds a small synthetic solo fixture, runs STAR 2.7.11b and rustar-aligner
with the same flags, and reports the per-read POS and leading-soft-clip
differences. This is the measurement behind issue #199, shrunk to something
that runs in seconds instead of needing the 10x mouse chr19 dataset.
"""
import os
import random
import subprocess
import sys
from pathlib import Path

TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
CB_LEN, UMI_LEN = 16, 12
READ_LEN = 90
N_READS = 400

RUSTAR = sys.argv[1] if len(sys.argv) > 1 else "./target/release/rustar-aligner"
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "/tmp/cr4diff")


def lcg(seed, n):
    bases = "ACGT"
    state = seed
    out = []
    for _ in range(n):
        state = (state * 1103515245 + 12345) & 0xFFFFFFFF
        out.append(bases[(state >> 16) & 3])
    return "".join(out)


def main():
    rng = random.Random(20260827)
    OUT.mkdir(parents=True, exist_ok=True)
    genome = lcg(88888, 20000)
    (OUT / "genome.fa").write_text(">chr1\n" + genome + "\n")

    # A minimal annotation: one long exon, so gene assignment never filters.
    (OUT / "genes.gtf").write_text(
        'chr1\tsyn\texon\t1\t20000\t.\t+\t.\tgene_id "G1"; transcript_id "G1_T1";\n'
    )

    cdna, barcode, whitelist = [], [], []
    for i in range(N_READS):
        start = rng.randrange(200, 19000 - READ_LEN)
        body = genome[start : start + READ_LEN]
        # Four kinds of read: clean, full TSO prefix, TSO a few bases in,
        # TSO with mismatches. Each exercises a different branch of the rule.
        kind = i % 4
        if kind == 0:
            seq = body
        elif kind == 1:
            seq = (TSO + body)[:READ_LEN]
        elif kind == 2:
            seq = ("GATC" + TSO + body)[:READ_LEN]
        else:
            tso = list(TSO)
            for p in (3, 11, 19):
                tso[p] = "C" if tso[p] == "A" else "A"
            seq = ("".join(tso) + body)[:READ_LEN]
        cdna.append((f"r{i}", seq))
        cb = lcg(1000 + (i % 8), CB_LEN)
        umi = lcg(7000 + i, UMI_LEN)
        barcode.append((f"r{i}", cb + umi))
        whitelist.append(cb)

    def write_fq(path, records):
        with open(path, "w") as f:
            for name, seq in records:
                f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

    write_fq(OUT / "cdna.fq", cdna)
    write_fq(OUT / "bc.fq", barcode)
    (OUT / "whitelist.txt").write_text("\n".join(sorted(set(whitelist))) + "\n")

    star_idx, rustar_idx = OUT / "star_idx", OUT / "rustar_idx"
    for idx, exe in ((star_idx, "STAR"), (rustar_idx, RUSTAR)):
        idx.mkdir(exist_ok=True)
        subprocess.run(
            [exe, "--runMode", "genomeGenerate", "--genomeDir", str(idx),
             "--genomeFastaFiles", str(OUT / "genome.fa"),
             "--genomeSAindexNbases", "7",
             "--sjdbGTFfile", str(OUT / "genes.gtf"), "--sjdbOverhang", "89",
             "--outFileNamePrefix", str(OUT / f"{idx.name}_")],
            check=True, capture_output=True,
        )

    common = [
        "--readFilesIn", str(OUT / "cdna.fq"), str(OUT / "bc.fq"),
        "--soloType", "CB_UMI_Simple",
        "--soloCBwhitelist", str(OUT / "whitelist.txt"),
        "--soloCBstart", "1", "--soloCBlen", str(CB_LEN),
        "--soloUMIstart", str(CB_LEN + 1), "--soloUMIlen", str(UMI_LEN),
        "--soloFeatures", "Gene",
        "--sjdbGTFfile", str(OUT / "genes.gtf"),
        *(["--clipAdapterType", "CellRanger4"] if os.environ.get("CR4", "1") == "1" else []),
        "--clip5pNbases", "5",
        "--clip3pNbases", "3",
        "--outSAMtype", "SAM",
    ]
    subprocess.run(["STAR", "--genomeDir", str(star_idx), *common,
                    "--outFileNamePrefix", str(OUT / "star_")],
                   check=True, capture_output=True)
    subprocess.run([RUSTAR, "--runMode", "alignReads", "--genomeDir", str(rustar_idx),
                    *common, "--outFileNamePrefix", str(OUT / "rustar_")],
                   check=True, capture_output=True)

    def primary(path):
        rows = {}
        for line in open(path):
            if line.startswith("@"):
                continue
            f = line.split("\t")
            flag = int(f[1])
            if flag & 0x900 or flag & 0x4:
                continue
            rows[f[0]] = (int(f[3]), f[5])
        return rows

    a = primary(OUT / "star_Aligned.out.sam")
    b = primary(OUT / "rustar_Aligned.out.sam")
    shared = sorted(set(a) & set(b))

    def lead_clip(cigar):
        n = ""
        for c in cigar:
            if c.isdigit():
                n += c
            else:
                return int(n) if c == "S" else 0
        return 0

    deltas = {}
    clip_deltas = {}
    for r in shared:
        d = b[r][0] - a[r][0]
        deltas[d] = deltas.get(d, 0) + 1
        cd = lead_clip(b[r][1]) - lead_clip(a[r][1])
        clip_deltas[cd] = clip_deltas.get(cd, 0) + 1

    print(f"STAR primary: {len(a)}  rustar primary: {len(b)}  shared: {len(shared)}")
    print("POS delta (rustar - STAR):", dict(sorted(deltas.items())))
    print("leading soft-clip delta:", dict(sorted(clip_deltas.items())))
    agree = deltas.get(0, 0)
    print(f"identical POS: {agree}/{len(shared)}")
    if len(shared) and agree == len(shared):
        print("CR4 CLIP DIFF: NONE")
    else:
        for r in shared[:5]:
            if b[r][0] != a[r][0]:
                print(f"  {r}: STAR {a[r]}  rustar {b[r]}")


main()
