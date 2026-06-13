#!/usr/bin/env python3
"""Runtime + output-stats benchmark: CellRanger vs STARsolo vs rustar-aligner.

Runs inside the amd64 benchmark container (test/Dockerfile.bench) so all three
tools run in one consistent Linux/x86_64 environment. Mouse GRCm39-2024-A
reference (built from the CellRanger refdata fasta+gtf for STAR/rust; CellRanger
uses the refdata directly), 5' GEM-X chemistry.

Each step is wall-clock + peak-RSS timed via /usr/bin/time -v. Output stats are
read from each tool's raw matrix (+ CellRanger metrics_summary.csv).

Usage (inside container):
  python3 test/solo_bench.py \
     --fasta REF/genome.fa --gtf REF/genes.gtf \
     --whitelist WL.txt --r1 R1.fq --r2 R2.fq \
     --cellranger /work/bench/cellranger-10.0.0/cellranger \
     --transcriptome /work/bench/refdata-gex-GRCm39-2024-A \
     --sample 5k_Mouse_PBMCs_5p_gem-x_GEX --fastqdir /work/bench/gex \
     --rustar /work/target-linux/release/rustar-aligner \
     --star $(which STAR) --threads 14 --mem-gb 36 --out /work/bench/results
"""
import argparse
import csv
import gzip
import json
import os
import re
import subprocess
import sys
import time

# CellRanger 4/5-matching solo flags (3' clip omitted; 5' chemistry).
SOLO_COMMON = [
    "--soloType", "CB_UMI_Simple",
    "--soloCBstart", "1", "--soloCBlen", "16",
    "--soloUMIstart", "17", "--soloUMIlen", "12",
    "--soloFeatures", "Gene",
    "--soloStrand", "Reverse",                 # 5' GEX (SC5P-R2 strandedness "-")
    "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
    "--soloUMIfiltering", "MultiGeneUMI_CR",
    "--soloUMIdedup", "1MM_CR",
]

TIME = ["/usr/bin/time", "-v"]


def timed(cmd, logpath, env=None):
    """Run cmd under /usr/bin/time -v; return (seconds, peak_rss_gb, ok)."""
    print("  $", " ".join(str(c) for c in cmd), flush=True)
    t0 = time.time()
    with open(logpath, "w") as lf:
        r = subprocess.run(TIME + list(map(str, cmd)), stdout=lf, stderr=subprocess.STDOUT, env=env)
    wall = time.time() - t0
    peak = None
    with open(logpath) as lf:
        txt = lf.read()
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)
    if m:
        peak = int(m.group(1)) / 1024 / 1024  # KB -> GB (GNU time reports KB)
    if r.returncode != 0:
        print(f"    !! exit {r.returncode}; tail:\n" + "\n".join(txt.splitlines()[-15:]))
    return wall, peak, r.returncode == 0


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def matrix_stats(raw_dir):
    """Read a MatrixMarket raw dir -> {n_barcodes_with_counts, total_umi, n_genes_detected}."""
    mtx = None
    for name in ("matrix.mtx.gz", "matrix.mtx"):
        p = os.path.join(raw_dir, name)
        if os.path.exists(p):
            mtx = p
            break
    if not mtx:
        return None
    cells, genes, total = set(), set(), 0
    with opener(mtx) as f:
        header_done = False
        for line in f:
            if line.startswith("%"):
                continue
            if not header_done:
                header_done = True  # dims line
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            g, c, v = int(parts[0]), int(parts[1]), int(float(parts[2]))
            if v > 0:
                genes.add(g)
                cells.add(c)
                total += v
    return {"n_barcodes_with_counts": len(cells), "n_genes_detected": len(genes), "total_umi": total}


def cellranger_metrics(outs_dir):
    p = os.path.join(outs_dir, "metrics_summary.csv")
    if not os.path.exists(p):
        return {}
    with open(p) as f:
        rows = list(csv.reader(f))
    if len(rows) >= 2:
        return dict(zip(rows[0], rows[1]))
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--whitelist", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--cellranger", required=True)
    ap.add_argument("--transcriptome", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastqdir", required=True)
    ap.add_argument("--rustar", required=True)
    ap.add_argument("--star", default="STAR")
    ap.add_argument("--threads", type=int, default=14)
    ap.add_argument("--mem-gb", type=int, default=36)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sa-nbases", default="14")
    ap.add_argument("--skip", default="", help="comma list: cellranger,star,rustar")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    logs = os.path.join(args.out, "logs")
    os.makedirs(logs, exist_ok=True)
    skip = set(s.strip() for s in args.skip.split(",") if s.strip())
    results = {}

    # ---- STARsolo -------------------------------------------------------
    if "star" not in skip:
        print("\n===== STARsolo =====")
        star_idx = os.path.join(args.out, "star_idx")
        os.makedirs(star_idx, exist_ok=True)
        s_gen, s_gen_rss, ok = timed(
            [args.star, "--runMode", "genomeGenerate", "--genomeDir", star_idx,
             "--genomeFastaFiles", args.fasta, "--sjdbGTFfile", args.gtf,
             "--sjdbOverhang", "89", "--genomeSAindexNbases", args.sa_nbases,
             "--runThreadN", args.threads],
            os.path.join(logs, "star_genomeGenerate.log"))
        star_out = os.path.join(args.out, "star_out") + "/"
        os.makedirs(star_out, exist_ok=True)
        s_run, s_run_rss, ok2 = timed(
            [args.star, "--genomeDir", star_idx, "--readFilesIn", args.r2, args.r1,
             "--runThreadN", args.threads, "--outSAMtype", "None",
             "--soloCBwhitelist", args.whitelist, "--outFileNamePrefix", star_out]
            + SOLO_COMMON,
            os.path.join(logs, "star_solo.log"))
        raw = os.path.join(star_out, "Solo.out", "Gene", "raw")
        results["STARsolo"] = {
            "index_build_s": round(s_gen, 1), "index_build_rss_gb": round(s_gen_rss or 0, 2),
            "count_s": round(s_run, 1), "count_rss_gb": round(s_run_rss or 0, 2),
            "stats": matrix_stats(raw), "ok": ok and ok2,
        }

    # ---- rustar-aligner -------------------------------------------------
    if "rustar" not in skip:
        print("\n===== rustar-aligner =====")
        rust_idx = os.path.join(args.out, "rust_idx")
        os.makedirs(rust_idx, exist_ok=True)
        r_gen, r_gen_rss, ok = timed(
            [args.rustar, "--runMode", "genomeGenerate", "--genomeDir", rust_idx,
             "--genomeFastaFiles", args.fasta, "--sjdbGTFfile", args.gtf,
             "--sjdbOverhang", "89", "--genomeSAindexNbases", args.sa_nbases,
             "--runThreadN", args.threads],
            os.path.join(logs, "rustar_genomeGenerate.log"))
        rust_out = os.path.join(args.out, "rust_out") + "/"
        os.makedirs(rust_out, exist_ok=True)
        r_run, r_run_rss, ok2 = timed(
            [args.rustar, "--genomeDir", rust_idx, "--readFilesIn", args.r2, args.r1,
             "--sjdbGTFfile", args.gtf, "--runThreadN", args.threads,
             "--outSAMtype", "SAM",
             "--soloCBwhitelist", args.whitelist, "--outFileNamePrefix", rust_out]
            + SOLO_COMMON,
            os.path.join(logs, "rustar_solo.log"))
        raw = os.path.join(rust_out, "Solo.out", "Gene", "raw")
        results["rustar-aligner"] = {
            "index_build_s": round(r_gen, 1), "index_build_rss_gb": round(r_gen_rss or 0, 2),
            "count_s": round(r_run, 1), "count_rss_gb": round(r_run_rss or 0, 2),
            "stats": matrix_stats(raw), "ok": ok and ok2,
        }

    # ---- CellRanger -----------------------------------------------------
    if "cellranger" not in skip:
        print("\n===== CellRanger =====")
        cr_dir = os.path.join(args.out, "cr")
        # cellranger count writes to ./<id>; run in args.out
        if os.path.exists(os.path.join(args.out, "cr_run")):
            subprocess.run(["rm", "-rf", os.path.join(args.out, "cr_run")])
        c_run, c_rss, ok = timed(
            [args.cellranger, "count", "--id", "cr_run",
             "--transcriptome", args.transcriptome,
             "--fastqs", args.fastqdir, "--sample", args.sample,
             "--create-bam", "false", "--nosecondary",
             "--localcores", str(args.threads), "--localmem", str(args.mem_gb)],
            os.path.join(logs, "cellranger_count.log"),
            env={**os.environ})
        outs = os.path.join(args.out, "cr_run", "outs")
        raw = os.path.join(outs, "raw_feature_bc_matrix")
        results["CellRanger"] = {
            "count_s": round(c_run, 1), "count_rss_gb": round(c_rss or 0, 2),
            "stats": matrix_stats(raw),
            "metrics": cellranger_metrics(outs), "ok": ok,
        }

    # ---- report ---------------------------------------------------------
    with open(os.path.join(args.out, "benchmark.json"), "w") as f:
        json.dump(results, f, indent=2)

    print("\n================ BENCHMARK SUMMARY ================")
    hdr = f"{'tool':<16}{'idx build(s)':>14}{'count(s)':>11}{'peak RSS(GB)':>14}{'barcodes':>10}{'genes':>8}{'total UMI':>12}"
    print(hdr)
    print("-" * len(hdr))
    for tool, r in results.items():
        st = r.get("stats") or {}
        idx = r.get("index_build_s", "-")
        peak = max(r.get("index_build_rss_gb", 0) or 0, r.get("count_rss_gb", 0) or 0)
        print(f"{tool:<16}{str(idx):>14}{str(r.get('count_s','-')):>11}{peak:>14.2f}"
              f"{str(st.get('n_barcodes_with_counts','-')):>10}"
              f"{str(st.get('n_genes_detected','-')):>8}{str(st.get('total_umi','-')):>12}")
    if "CellRanger" in results and results["CellRanger"].get("metrics"):
        m = results["CellRanger"]["metrics"]
        keys = ["Estimated Number of Cells", "Mean Reads per Cell", "Median Genes per Cell",
                "Median UMI Counts per Cell", "Reads Mapped Confidently to Transcriptome"]
        print("\nCellRanger reported metrics:")
        for k in keys:
            if k in m:
                print(f"  {k}: {m[k]}")
    print(f"\nFull results: {os.path.join(args.out, 'benchmark.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
