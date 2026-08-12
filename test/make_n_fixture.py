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
    chimeric_reads.fq    synthetic chimeric reads (single-end)
    pe_intra_{1,2}.fq    paired: mate1 is itself chimeric, mate2 is a normal
                         read from the donor locus. Exercises the intra-mate
                         (Tier 2) PE path, where the junction lies inside one
                         mate and STAR's `N` scan applies.
    pe_inter_{1,2}.fq    paired: donor on mate1, acceptor on mate2, so the mates
                         *bracket* the junction. STAR takes its `EX_iFrag`
                         branch here and applies no `N` scan at all, so these
                         two files together separate the filtered path from the
                         exempt one.
    reads.tsv            read name, the two loci, strands, and the expected
                         relationship to a gap (`clean` / `donor_gap` /
                         `acceptor_gap` / `read_n`)

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


def write_pe(outdir, reads, table, original, seg_len):
    """Two paired layouts from the same chimeric reads.

    `pe_intra`: the chimeric read stays whole as mate1 and gets an ordinary
    mate2 drawn downstream of the donor locus. The junction is inside mate1, so
    detection runs the intra-mate path and the `N` scan applies.

    `pe_inter`: the same read split at the junction, donor to mate1 and acceptor
    to mate2. The mates now bracket the junction, which STAR handles in a branch
    that never reads a genomic base.

    Mate2 is reverse-complemented in both, as an FR pair.
    """
    # table rows are (name, donor_chr, donor_pos0, donor_strand,
    #                 acceptor_chr, acceptor_pos0, acceptor_strand, kind)
    by_name = {t[0]: (t[4], t[5]) for t in table}

    def fq(path, records):
        with open(path, "w") as fh:
            for n, s in records:
                fh.write(f"@{n}\n{s}\n+\n{'I' * len(s)}\n")

    intra1, intra2, inter1, inter2 = [], [], [], []
    for name, seq in reads:
        acc_chr, acc_pos = by_name[name]
        # Mate2 for the intra layout sits on the ACCEPTOR, downstream of the
        # junction — the shape of a real fusion-spanning pair. Drawing it from
        # the donor instead leaves most of the fragment explainable linearly,
        # and STAR then rejects the chimera on `chimScoreDropMax`.
        m2_start = acc_pos + seg_len
        mate2 = original[acc_chr][m2_start : m2_start + seg_len]
        if len(mate2) < seg_len:
            continue
        intra1.append((name, seq))
        intra2.append((name, revcomp(mate2)))

        half = len(seq) // 2
        inter1.append((name, seq[:half]))
        inter2.append((name, revcomp(seq[half:])))

    fq(os.path.join(outdir, "pe_intra_1.fq"), intra1)
    fq(os.path.join(outdir, "pe_intra_2.fq"), intra2)
    fq(os.path.join(outdir, "pe_inter_1.fq"), inter1)
    fq(os.path.join(outdir, "pe_inter_2.fq"), inter2)
    print(f"wrote {len(intra1)} intra-mate and {len(inter1)} inter-mate PE pairs")


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

    write_pe(args.out, reads, table, dict(records), args.seg_len)
    with open(os.path.join(args.out, "reads.tsv"), "w") as fh:
        fh.write("read\tdonor_chr\tdonor_pos0\tdonor_strand\tacceptor_chr\tacceptor_pos0\tacceptor_strand\texpect\n")
        for row in table:
            fh.write("\t".join(str(x) for x in row) + "\n")

    n_bases = sum(s.count("N") for _, s in mutated)
    print(f"planted {len(gaps)} gaps ({n_bases} N bases) across {len(mutated)} sequences")
    print(f"wrote {len(reads)} chimeric reads to {args.out}/chimeric_reads.fq")


if __name__ == "__main__":
    main()
