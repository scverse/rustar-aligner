#!/usr/bin/env python3
"""Summarise the bulk total-RNA benchmark produced by run.sh.

    python3 scripts/bench_bulk_unspliced/analyze.py "$DATA"

Per library (run x pseudo-replicate) and TranscriptomeSAM mode (star = STAR
projection, intron / premrna = --quantTranscriptomeUnspliced Intron / PreMRNA):
  - share of Salmon-assigned spliced-target fragments on GENCODE
    retained_intron isoforms (RI share);
  - share of Salmon-assigned fragments on the <gene_id>-I targets;
  - unspliced / ambiguous fractions from ReadsPerGeneSplicing.summary.tsv
    (star runs; strand column picked from the data: most assigned reads);
  - fragments in Aligned.toTranscriptome.out.bam per input pair, and the
    Salmon mapping rate;
  - wall time and peak RSS of the rustar-aligner run.
Across libraries: Spearman rho and least-squares slope of RI share against
unspliced fraction.
Between the two pseudo-replicates of each run: Pearson r of within-gene
isoform proportions (annotated isoforms only) and of log tximport-style
average transcript lengths, for genes with >= MIN_READS spliced fragments.
Needs numpy and scipy only.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

MIN_READS = 100
MODES = ("star", "intron", "premrna")

data = Path(sys.argv[1])
out = data / "out"

# transcript_id -> (gene_id, transcript_type)
tx = {}
with open(data / "ref" / "gencode.v50.annotation.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.split("\t", 9)
        if c[2] != "transcript":
            continue
        a = c[8]
        tid = re.search(r'transcript_id "([^"]+)"', a).group(1)
        gid = re.search(r'gene_id "([^"]+)"', a).group(1)
        tt = re.search(r'transcript_type "([^"]+)"', a).group(1)
        tx[tid] = (gid, tt)


def quant(path):
    names, efflen, reads = [], [], []
    with open(path) as f:
        next(f)
        for line in f:
            n, _l, el, _tpm, nr = line.rstrip("\n").split("\t")
            names.append(n)
            efflen.append(float(el))
            reads.append(float(nr))
    return names, np.array(efflen), np.array(reads)


def time_log(path):
    txt = path.read_text()
    m = re.search(r"([\d.]+) real", txt)
    wall = float(m.group(1)) if m else float("nan")
    m = re.search(r"(\d+)\s+maximum resident set size", txt)  # macOS, bytes
    if m:
        rss = int(m.group(1)) / 2**30
    else:
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt)  # GNU
        rss = int(m.group(1)) / 2**20 if m else float("nan")
        m = re.search(r"Elapsed \(wall clock\) time.*: ([\d:.]+)", txt)
        if m:
            parts = [float(p) for p in m.group(1).split(":")]
            wall = sum(p * 60**i for i, p in enumerate(reversed(parts)))
    return wall, rss


def log_final(path):
    d = {}
    for line in path.read_text().splitlines():
        if "|" in line:
            k, v = line.split("|", 1)
            d[k.strip()] = v.strip()
    return d


def splicing(path):
    rows = {}
    for line in path.read_text().splitlines()[1:]:
        k, *v = line.split("\t")
        rows[k] = v
    assigned = [sum(int(rows[k][i]) for k in ("N_spliced", "N_unspliced", "N_ambiguous")) for i in range(3)]
    col = 1 if assigned[1] > assigned[2] else 2
    return {
        "strand": ["unstranded", "forward", "reverse"][col],
        "unspliced": float(rows["fraction_unspliced"][col]),
        "ambiguous": float(rows["fraction_ambiguous"][col]),
    }


rows = []
q = {}
for libdir in sorted(p for p in out.iterdir() if p.is_dir()):
    lib = libdir.name
    for mode in MODES:
        d = libdir / mode
        sf = d / "salmon" / "quant.sf"
        if not sf.exists():
            continue
        names, efflen, reads = quant(sf)
        unspl = np.array([n.endswith("-I") and n not in tx for n in names])
        ri = np.array([tx.get(n, ("", ""))[1] == "retained_intron" for n in names])
        q[(lib, mode)] = (names, efflen, reads, unspl)
        meta = json.loads((d / "salmon" / "aux_info" / "meta_info.json").read_text())
        lf = log_final(d / "Log.final.out")
        n_in = int(lf["Number of input reads"])
        wall, rss = time_log(d / "time.log")
        spliced_total = reads[~unspl].sum()
        r = {
            "lib": lib,
            "mode": mode,
            "input_pairs": n_in,
            "unique_pct": lf["Uniquely mapped reads %"],
            "trsam_frag_per_pair": meta["num_processed"] / n_in,
            "salmon_rate": meta["percent_mapped"] / 100,
            "ri_share": reads[ri].sum() / spliced_total,
            "unspliced_target_share": reads[unspl].sum() / reads.sum(),
            "wall_s": wall,
            "rss_gib": rss,
        }
        vs = d / "ReadsPerGeneSplicing.summary.tsv"
        if vs.exists():
            r.update(splicing(vs))
        rows.append(r)

print("| library | mode | pairs | unique % | trSAM frag/pair | Salmon mapped | RI share | -I share | unspliced (GeneSplicing) | wall s | RSS GiB |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(
        f"| {r['lib']} | {r['mode']} | {r['input_pairs']} | {r['unique_pct']} | "
        f"{r['trsam_frag_per_pair']:.3f} | {r['salmon_rate']:.3f} | {r['ri_share']:.4f} | "
        f"{r['unspliced_target_share']:.4f} | {r.get('unspliced', float('nan')):.4f} | "
        f"{r['wall_s']:.0f} | {r['rss_gib']:.1f} |"
    )

by = {m: {r["lib"]: r for r in rows if r["mode"] == m} for m in MODES}
libs = sorted(set.intersection(*(set(v) for v in by.values())))
unspliced = {l: by["star"][l]["unspliced"] for l in libs}
print()
print("GeneSplicing strand column:", sorted({by["star"][l].get("strand") for l in libs}))
for subset, label in ((libs, "all libraries"), ([l for l in libs if l.endswith("rep1")], "rep1 only")):
    if len(subset) < 3:
        continue
    u = [unspliced[l] for l in subset]
    parts = []
    for m in MODES:
        s = [by[m][l]["ri_share"] for l in subset]
        slope = np.polyfit(u, s, 1)[0]
        parts.append(
            f"{m}: RI share mean {np.mean(s):.4f}, rho vs unspliced {spearmanr(s, u)[0]:.2f}, "
            f"slope {slope:.3f}"
        )
    print(f"{label} (n={len(subset)}), unspliced fraction mean {np.mean(u):.4f}: " + "; ".join(parts))
for m in MODES:
    print(
        f"{m}: mean trSAM frag/pair {np.mean([by[m][l]['trsam_frag_per_pair'] for l in libs]):.3f}, "
        f"Salmon mapped {np.mean([by[m][l]['salmon_rate'] for l in libs]):.3f}, "
        f"-I share {np.mean([by[m][l]['unspliced_target_share'] for l in libs]):.4f}, "
        f"wall {np.mean([by[m][l]['wall_s'] for l in libs]):.0f} s, "
        f"peak RSS {np.max([by[m][l]['rss_gib'] for l in libs]):.1f} GiB"
    )


def rep_stats(mode):
    iso_r, len_r, n_genes = [], [], []
    runs = sorted({l.rsplit("_rep", 1)[0] for l in libs})
    for run in runs:
        a, b = q.get((run + "_rep1", mode)), q.get((run + "_rep2", mode))
        if not a or not b:
            continue
        names, efflen, ra, unspl = a
        names_b, efflen_b, rb, _ = b
        assert names == names_b
        genes = {}
        for i, n in enumerate(names):
            if not unspl[i]:
                genes.setdefault(tx.get(n, (n, "?"))[0], []).append(i)
        pa, pb, la, lb = [], [], [], []
        for idx in genes.values():
            if len(idx) < 2:
                continue
            idx = np.array(idx)
            ga, gb = ra[idx].sum(), rb[idx].sum()
            if ga < MIN_READS or gb < MIN_READS:
                continue
            pa.extend(ra[idx] / ga)
            pb.extend(rb[idx] / gb)
            # tximport avgTxLength: abundance-weighted effective length
            ta, tb = ra[idx] / efflen[idx], rb[idx] / efflen_b[idx]
            la.append(np.log((ta / ta.sum() * efflen[idx]).sum()))
            lb.append(np.log((tb / tb.sum() * efflen_b[idx]).sum()))
        iso_r.append(pearsonr(pa, pb)[0])
        len_r.append(pearsonr(la, lb)[0])
        n_genes.append(len(la))
    return iso_r, len_r, n_genes


print()
for m in MODES:
    iso_r, len_r, n_genes = rep_stats(m)
    if iso_r:
        print(
            f"pseudo-replicates, {m}: isoform-proportion r median {np.median(iso_r):.3f} "
            f"(range {min(iso_r):.3f}-{max(iso_r):.3f}); log avgTxLength r median {np.median(len_r):.3f} "
            f"(range {min(len_r):.3f}-{max(len_r):.3f}); multi-isoform genes with >= {MIN_READS} "
            f"fragments: median {int(np.median(n_genes))}"
        )

idx_time = data / "index" / "time.log"
if idx_time.exists():
    w, m = time_log(idx_time)
    print(f"\ngenomeGenerate: wall {w:.0f} s, peak RSS {m:.1f} GiB")
