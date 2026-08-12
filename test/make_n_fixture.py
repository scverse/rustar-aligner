#!/usr/bin/env python3
"""Build a reference with planted `N` runs, plus chimeric reads that straddle them.

Assembly gaps are the reason `--chimFilter banGenomicN` exists, and no reference
small enough to test with conveniently has any: the yeast R64-1-1 toplevel FASTA
contains zero `N` bases. So plant them, deterministically, and synthesise reads
whose chimeric junctions land on and beside them.

Nothing here is random. Gap positions are computed from the sequence length and
read construction is table-driven, so two runs on the same input produce
byte-identical output and a failure can be reproduced from the manifest alone.

    python3 test/make_n_fixture.py --fasta REF.fa --out DIR [--gap-every 200000]

Writes into DIR:

    genome_N.fa          the reference with `N` runs planted
    gaps.tsv             chr, 0-based start, length -- one row per planted run
    chimeric_reads.fq    synthetic chimeric reads
    reads.tsv            read name, the two loci, strands, and the expected
                         relationship to a gap (`clean` / `donor_gap` /
                         `acceptor_gap`)

The read table deliberately covers both strands on both sides of the junction:
the genomic base STAR examines differs by strand (it mirrors the offset within
the segment and complements the base), so a filter that ignores strand passes a
forward-only fixture and fails here.
"""

import argparse
import os


def read_fasta(path):
    """Return [(name, sequence)] with sequence upper-cased and newline-free."""
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if names:
                    seqs.append("".join(cur))
                    cur = []
                names.append(line[1:].split()[0])
            else:
                cur.append(line.upper())
    if names:
        seqs.append("".join(cur))
    return list(zip(names, seqs))


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def plant_gaps(records, gap_every, gap_len, edge):
    """Plant `N` runs at regular offsets, away from sequence ends.

    Regular spacing keeps the manifest readable and means a read built at a
    known locus has a predictable distance to the nearest gap. `edge` keeps runs
    clear of contig ends so a junction can sit on either side of one.
    """
    out, gaps = [], []
    for name, seq in records:
        s = list(seq)
        pos = gap_every
        while pos + gap_len + edge < len(seq):
            for i in range(pos, pos + gap_len):
                s[i] = "N"
            gaps.append((name, pos, gap_len))
            pos += gap_every
        out.append((name, "".join(s)))
    return out, gaps


COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def build_reads(original, mutated, gaps, seg_len):
    """Chimeric reads: `seg_len` from one locus, `seg_len` from another.

    Each row fixes both segments' strands explicitly rather than sampling, so
    the four strand combinations are all present and named in the manifest.
    """
    # Segments come from the ORIGINAL sequence, so a read carries real bases
    # where the mutated genome carries `N`. That separates the two rules STAR
    # applies at a chimeric junction: `banGenomicN` looks at the genome, while
    # the read-`N` ban (`bR>3`) is unconditional and fires even under
    # `--chimFilter None`. A fixture built from the mutated genome would trip
    # the second rule everywhere and never isolate the first.
    by_name = dict(original)
    mut_by_name = dict(mutated)
    # Two chromosomes with at least one planted gap each, in manifest order.
    chroms = []
    for name, pos, glen in gaps:
        if name not in [c[0] for c in chroms]:
            chroms.append((name, pos, glen))
        if len(chroms) == 2:
            break
    if len(chroms) < 2:
        raise SystemExit("need planted gaps on at least two sequences")
    (cA, gapA, glenA), (cB, gapB, glenB) = chroms

    # A locus that is comfortably clear of any gap, for the `clean` controls.
    cleanA, cleanB = gapA // 2, gapB // 2

    plan = []
    for d_rev in (False, True):
        for a_rev in (False, True):
            tag = f"{'r' if d_rev else 'f'}{'r' if a_rev else 'f'}"
            # Clean: both segments far from any planted run.
            plan.append((f"clean_{tag}", cA, cleanA, cB, cleanB, d_rev, a_rev, "clean"))
            # Donor's junction-side end abuts the gap.
            plan.append(
                (f"donorgap_{tag}", cA, gapA - seg_len, cB, cleanB, d_rev, a_rev, "donor_gap")
            )
            # Acceptor starts inside the gap.
            plan.append(
                (f"acceptorgap_{tag}", cA, cleanA, cB, gapB, d_rev, a_rev, "acceptor_gap")
            )

    # A few reads that DO carry `N`, to exercise the unconditional read ban.
    for d_rev in (False, True):
        plan.append(
            (f"readn_{'r' if d_rev else 'f'}f", cA, gapA, cB, cleanB, d_rev, False, "read_n")
        )

    reads, table = [], []
    for name, c0, p0, c1, p1, d_rev, a_rev, kind in plan:
        src = mut_by_name if kind == "read_n" else by_name
        d = src[c0][p0 : p0 + seg_len]
        a = src[c1][p1 : p1 + seg_len]
        if len(d) < seg_len or len(a) < seg_len:
            continue
        seq = (revcomp(d) if d_rev else d) + (revcomp(a) if a_rev else a)
        reads.append((name, seq))
        table.append((name, c0, p0, "-" if d_rev else "+", c1, p1, "-" if a_rev else "+", kind))
    return reads, table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="reference FASTA to mutate")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--gap-every", type=int, default=200_000)
    ap.add_argument("--gap-len", type=int, default=50)
    ap.add_argument("--edge", type=int, default=50_000, help="keep gaps this far from contig ends")
    ap.add_argument("--seg-len", type=int, default=60, help="per-segment read length")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    records = read_fasta(args.fasta)
    mutated, gaps = plant_gaps(records, args.gap_every, args.gap_len, args.edge)

    write_fasta(os.path.join(args.out, "genome_N.fa"), mutated)
    with open(os.path.join(args.out, "gaps.tsv"), "w") as fh:
        fh.write("chr\tstart0\tlength\n")
        for name, pos, glen in gaps:
            fh.write(f"{name}\t{pos}\t{glen}\n")

    reads, table = build_reads(records, mutated, gaps, args.seg_len)
    with open(os.path.join(args.out, "chimeric_reads.fq"), "w") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    with open(os.path.join(args.out, "reads.tsv"), "w") as fh:
        fh.write("read\tdonor_chr\tdonor_pos0\tdonor_strand\tacceptor_chr\tacceptor_pos0\tacceptor_strand\texpect\n")
        for row in table:
            fh.write("\t".join(str(x) for x in row) + "\n")

    n_bases = sum(s.count("N") for _, s in mutated)
    print(f"planted {len(gaps)} gaps ({n_bases} N bases) across {len(mutated)} sequences")
    print(f"wrote {len(reads)} chimeric reads to {args.out}/chimeric_reads.fq")


if __name__ == "__main__":
    main()
