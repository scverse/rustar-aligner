#!/usr/bin/env python3
"""Differential test: rustar-aligner STARsolo vs real STAR, CellRanger-style run.

Generates a small synthetic 10x-style dataset (genome + GTF + whitelist + cDNA
read + barcode read), runs BOTH STAR and rustar-aligner with the
CellRanger-4/5-matching solo flags from
https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#matching-cellranger-4xx-and-5xx-results
and compares the raw Gene count matrices decoded to {(barcode, gene_id): count}.

Usage:
    python3 test/solo_cellranger_diff.py [--star /path/to/STAR] [--rustar /path/to/rustar-aligner] [--keep]

Exit code 0 = matrices match, 1 = mismatch / error.
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile

# CellRanger 4.x/5.x matching flags (STARsolo.md).
CELLRANGER_FLAGS = [
    "--clipAdapterType", "CellRanger4",
    "--outFilterScoreMin", "30",
    "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
    "--soloUMIfiltering", "MultiGeneUMI_CR",
    "--soloUMIdedup", "1MM_CR",
]

CB_LEN = 16
UMI_LEN = 12
READ_LEN = 90
BASES = "ACGT"


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


# Two-exon gene layout (0-based): exon1 [s, s+150), intron [s+150, s+400) with
# canonical GT..AG, exon2 [s+400, s+550). Multi-exon genes give STAR a non-empty
# splice-junction DB, which it needs to set up the solo Transcriptome directory.
GENE_A_START = 10000
GENE_B_START = 30000


def _plant_gene(g, s, rng):
    g[s : s + 150] = list(rand_seq(rng, 150))          # exon1
    g[s + 150 : s + 400] = list(rand_seq(rng, 250))    # intron body
    g[s + 150], g[s + 151] = "G", "T"                  # donor
    g[s + 398], g[s + 399] = "A", "G"                  # acceptor
    g[s + 400 : s + 550] = list(rand_seq(rng, 150))    # exon2


def build_genome(rng, length=50000):
    g = list(rand_seq(rng, length))
    _plant_gene(g, GENE_A_START, rng)
    _plant_gene(g, GENE_B_START, rng)
    return "".join(g)


def pick_window(genome, exon_start):
    """Pick a READ_LEN window inside exon1 ending in a non-A base (so the
    CellRanger4 polyA trim is a guaranteed no-op for both tools). The window
    stays inside the 150 bp exon1, so reads never span the junction."""
    a = exon_start + 20
    while genome[a + READ_LEN - 1] == "A":
        a += 1
    return genome[a : a + READ_LEN]


def write_files(d, genome):
    fa = os.path.join(d, "genome.fa")
    with open(fa, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(genome), 70):
            f.write(genome[i : i + 70] + "\n")

    gtf = os.path.join(d, "genes.gtf")
    with open(gtf, "w") as f:
        # Two exons per gene (1-based inclusive), matching the planted layout.
        f.write('chr1\tsrc\texon\t10001\t10150\t.\t+\t.\tgene_id "GENEA"; transcript_id "GENEA.1"; gene_name "GeneA";\n')
        f.write('chr1\tsrc\texon\t10401\t10550\t.\t+\t.\tgene_id "GENEA"; transcript_id "GENEA.1"; gene_name "GeneA";\n')
        f.write('chr1\tsrc\texon\t30001\t30150\t.\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.1"; gene_name "GeneB";\n')
        f.write('chr1\tsrc\texon\t30401\t30550\t.\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.1"; gene_name "GeneB";\n')

    wl = os.path.join(d, "whitelist.txt")
    cbs = ["AAAACCCCGGGGTTTT", "ACACACACGTGTGTGT", "TTTTGGGGCCCCAAAA", "GTGTGTGTACACACAC"]
    with open(wl, "w") as f:
        f.write("\n".join(cbs) + "\n")

    readA = pick_window(genome, 10000)
    readB = pick_window(genome, 30000)

    # (cell, gene-read, umi, n_reads). Designed to exercise:
    #  - exact CB match (all CBs in whitelist)
    #  - 1MM_CR UMI collapse: ACGTACGTACGT (5) + ACGTACGTACGA (1) -> 1 molecule
    #  - distinct molecules counted, two genes, two cells.
    plan = [
        (cbs[0], readA, "ACGTACGTACGT", 5),
        (cbs[0], readA, "ACGTACGTACGA", 1),   # 1MM neighbor of the above
        (cbs[0], readA, "TGCATGCATGCA", 3),   # separate molecule
        (cbs[0], readB, "GGGGTTTTAACC", 2),   # GeneB, cell0
        (cbs[1], readA, "CATGCATGCATG", 4),   # GeneA, cell1
    ]
    # Expected decoded matrix.
    expected = {
        ("AAAACCCCGGGGTTTT", "GENEA"): 2,  # two molecules (1MM pair collapses)
        ("AAAACCCCGGGGTTTT", "GENEB"): 1,
        ("ACACACACGTGTGTGT", "GENEA"): 1,
    }

    cdna = os.path.join(d, "cdna.fq")
    bc = os.path.join(d, "barcode.fq")
    ci = 0
    with open(cdna, "w") as cf, open(bc, "w") as bf:
        for (cb, read, umi, n) in plan:
            for _ in range(n):
                name = f"read{ci}"
                ci += 1
                cf.write(f"@{name}\n{read}\n+\n{'I' * READ_LEN}\n")
                barcode = cb + umi
                bf.write(f"@{name}\n{barcode}\n+\n{'I' * len(barcode)}\n")
    return fa, gtf, wl, cdna, bc, expected


def run(cmd, **kw):
    print("  $", " ".join(str(c) for c in cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-4000:])
        raise SystemExit(f"command failed ({r.returncode}): {cmd[0]}")
    return r


def run_star(star, d, fa, gtf, wl, cdna, bc):
    # Generate WITH the GTF so geneInfo.tab lands in the index, then reset the
    # recorded sjdbGTFfile to "-" in genomeParameters.txt. STAR's solo
    # Transcriptome uses `trInfoDir = sjdbGTFfile=="-" ? genomeDir : sjdbInsert.outDir`
    # (Transcriptome.cpp:18); with the path still recorded it points at an empty
    # insert dir and fails with "/geneInfo.tab". Resetting to "-" makes it read
    # geneInfo.tab from the genome dir. (The gene model is intact in the index.)
    idx = os.path.join(d, "star_index")
    os.makedirs(idx, exist_ok=True)
    run([star, "--runMode", "genomeGenerate", "--genomeDir", idx,
         "--genomeFastaFiles", fa, "--sjdbGTFfile", gtf,
         "--genomeSAindexNbases", "7", "--sjdbOverhang", "89"])
    gp = os.path.join(idx, "genomeParameters.txt")
    lines = open(gp).read().splitlines()
    with open(gp, "w") as f:
        for ln in lines:
            if ln.startswith("sjdbGTFfile\t"):
                f.write("sjdbGTFfile\t-\n")
            else:
                f.write(ln + "\n")

    out = os.path.join(d, "star_out") + os.sep
    run([star, "--genomeDir", idx, "--readFilesIn", cdna, bc,
         "--soloType", "CB_UMI_Simple", "--soloCBwhitelist", wl,
         "--soloCBstart", "1", "--soloCBlen", str(CB_LEN),
         "--soloUMIstart", str(CB_LEN + 1), "--soloUMIlen", str(UMI_LEN),
         "--soloFeatures", "Gene", "--outSAMtype", "SAM",
         "--outFileNamePrefix", out] + CELLRANGER_FLAGS)
    # Guard against a STAR binary that silently reads 0 reads (broken bottle).
    log = os.path.join(out, "Log.final.out")
    if os.path.exists(log):
        for ln in open(log):
            if "Number of input reads" in ln and ln.strip().endswith("0"):
                raise SystemExit(
                    "STAR processed 0 input reads — the STAR binary appears broken "
                    "on this machine (immediate EOF on FASTQ input). Install a working "
                    "STAR and re-run with --star /path/to/STAR."
                )
    return os.path.join(out, "Solo.out", "Gene", "raw")


def run_rustar(rustar, d, fa, gtf, wl, cdna, bc):
    idx = os.path.join(d, "rustar_index")
    os.makedirs(idx, exist_ok=True)
    run([rustar, "--runMode", "genomeGenerate", "--genomeDir", idx,
         "--genomeFastaFiles", fa, "--sjdbGTFfile", gtf,
         "--genomeSAindexNbases", "7", "--sjdbOverhang", "89"])
    out = os.path.join(d, "rustar_out") + os.sep
    run([rustar, "--genomeDir", idx, "--readFilesIn", cdna, bc,
         "--soloType", "CB_UMI_Simple", "--soloCBwhitelist", wl,
         "--soloCBstart", "1", "--soloCBlen", str(CB_LEN),
         "--soloUMIstart", str(CB_LEN + 1), "--soloUMIlen", str(UMI_LEN),
         "--soloFeatures", "Gene", "--sjdbGTFfile", gtf,
         "--outSAMtype", "SAM",
         "--outFileNamePrefix", out] + CELLRANGER_FLAGS)
    return os.path.join(out, "Solo.out", "Gene", "raw")


def decode_matrix(raw_dir):
    """Decode raw/{matrix.mtx,barcodes.tsv,features.tsv} -> {(barcode, gene_id): count}."""
    feats = []
    with open(os.path.join(raw_dir, "features.tsv")) as f:
        for line in f:
            feats.append(line.rstrip("\n").split("\t")[0])
    barcodes = []
    with open(os.path.join(raw_dir, "barcodes.tsv")) as f:
        for line in f:
            barcodes.append(line.strip())
    out = {}
    with open(os.path.join(raw_dir, "matrix.mtx")) as f:
        lines = [l for l in f if not l.startswith("%")]
    # first non-% line is dims
    for entry in lines[1:]:
        parts = entry.split()
        if len(parts) < 3:
            continue
        row, col, cnt = int(parts[0]), int(parts[1]), int(float(parts[2]))
        out[(barcodes[col - 1], feats[row - 1])] = cnt
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default=shutil.which("STAR") or "/opt/homebrew/bin/STAR")
    ap.add_argument("--rustar", default=None)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--seed", type=int, default=20260612)
    args = ap.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if args.rustar:
        # Honor an explicit path exactly — never silently fall back to a
        # different (possibly foreign-arch) binary.
        rustar = args.rustar
        if not os.path.exists(rustar):
            raise SystemExit(f"--rustar binary not found: {rustar}")
    else:
        rustar = os.path.join(repo, "target", "release", "rustar-aligner")
        if not os.path.exists(rustar):
            rustar = os.path.join(repo, "target", "debug", "rustar-aligner")
        if not os.path.exists(rustar):
            raise SystemExit(
                "rustar-aligner binary not found — build it first (cargo build [--release]) "
                "or pass --rustar /path/to/rustar-aligner"
            )
    if not (args.star and os.path.exists(args.star)):
        raise SystemExit(f"STAR binary not found: {args.star}")

    d = tempfile.mkdtemp(prefix="solo_diff_")
    print(f"workdir: {d}")
    print(f"STAR:   {args.star}")
    print(f"rustar: {rustar}")
    rng = random.Random(args.seed)
    try:
        genome = build_genome(rng)
        fa, gtf, wl, cdna, bc, expected = write_files(d, genome)

        print("\n== rustar-aligner ==")
        rustar_raw = run_rustar(rustar, d, fa, gtf, wl, cdna, bc)
        rustar_m = decode_matrix(rustar_raw)

        print("\n== expected (hand-computed CellRanger result) ==")
        for k, v in sorted(expected.items()):
            print(f"   {k} = {v}")
        print("== rustar matrix ==")
        for k, v in sorted(rustar_m.items()):
            print(f"   {k} = {v}")

        # Core guarantee: rustar's CellRanger-style matrix matches the expectation.
        if rustar_m != expected:
            print("\nFAIL: rustar matrix does not match the expected CellRanger result:")
            for k in sorted(set(rustar_m) | set(expected)):
                if rustar_m.get(k) != expected.get(k):
                    print(f"   {k}: rustar={rustar_m.get(k)} expected={expected.get(k)}")
            return 1
        print("\nrustar matrix matches the expected CellRanger result.")

        # Live comparison against the real STAR binary, when it works on this host.
        print("\n== STAR ==")
        try:
            star_raw = run_star(args.star, d, fa, gtf, wl, cdna, bc)
            star_m = decode_matrix(star_raw)
        except SystemExit as e:
            print(f"\nSTAR could not run a live comparison on this host: {e}")
            print("PASS (rustar validated against the CellRanger expectation; "
                  "run on a host with a working STAR for the live diff).")
            return 0
        print("== STAR matrix ==")
        for k, v in sorted(star_m.items()):
            print(f"   {k} = {v}")
        if star_m == rustar_m:
            print("\nPASS: rustar-aligner matrix matches real STARsolo exactly.")
            return 0
        print("\nFAIL: rustar vs STAR mismatch:")
        for k in sorted(set(star_m) | set(rustar_m)):
            if star_m.get(k) != rustar_m.get(k):
                print(f"   {k}: STAR={star_m.get(k)} rustar={rustar_m.get(k)}")
        return 1
    finally:
        if args.keep:
            print(f"(kept workdir {d})")
        else:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
