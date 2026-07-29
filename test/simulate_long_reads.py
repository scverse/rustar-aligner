#!/usr/bin/env python3
"""Draw long reads from spliced transcripts, for the STARlong differential.

The bundled tiers are 150 bp Illumina reads, which never reach the long-read
code path: `--winReadCoverageRelativeMin` and `--winReadCoverageBasesMin` only
filter windows when the long-read stitcher runs, and a 150 bp read produces one
window with near-total coverage. Reads have to be long enough, and span enough
junctions, for the coverage filter to be able to discard anything.

Each read is the full spliced sequence of one transcript, so it crosses every
junction that transcript has. Transcripts shorter than `--min-length` are
skipped; longer ones are emitted whole. Errors are substitutions only, at a
fixed rate, drawn from a seeded RNG so the FASTQ is byte-reproducible.

    simulate_long_reads.py <genome.fa> <annotation.gtf> <out.fastq> [n] [min_len]
"""

import random
import sys
from collections import defaultdict

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_exons(path):
    """transcript_id -> (chrom, strand, [(start, end), ...]) in 1-based inclusive GTF coords."""
    exons = defaultdict(list)
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs = f[8]
            key = 'transcript_id "'
            i = attrs.find(key)
            if i < 0:
                continue
            j = attrs.find('"', i + len(key))
            tid = attrs[i + len(key) : j]
            exons[tid].append((int(f[3]), int(f[4])))
            meta[tid] = (f[0], f[6])
    return exons, meta


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    genome_path, gtf_path, out_path = sys.argv[1:4]
    n_reads = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    min_len = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
    error_rate = 0.02

    rng = random.Random(100)
    genome = read_fasta(genome_path)
    exons, meta = read_exons(gtf_path)

    # Multi-exon transcripts only: a single-exon read crosses no junction and
    # exercises none of the long-read window logic.
    candidates = []
    for tid, blocks in exons.items():
        if len(blocks) < 2:
            continue
        chrom, strand = meta[tid]
        if chrom not in genome:
            continue
        blocks = sorted(blocks)
        length = sum(e - s + 1 for s, e in blocks)
        if length >= min_len:
            candidates.append((tid, chrom, strand, blocks, length))

    candidates.sort()
    if not candidates:
        sys.exit(f"no multi-exon transcript of at least {min_len} bp found")
    rng.shuffle(candidates)

    written = 0
    with open(out_path, "w") as out:
        while written < n_reads:
            tid, chrom, strand, blocks, _ = candidates[written % len(candidates)]
            chrom_seq = genome[chrom]
            spliced = "".join(chrom_seq[s - 1 : e] for s, e in blocks).upper()
            if strand == "-":
                spliced = spliced.translate(COMPLEMENT)[::-1]
            if "N" in spliced:
                # An N run would be clipped rather than aligned and tells us
                # nothing about the coverage filter.
                written += 1
                continue
            bases = list(spliced)
            for i in range(len(bases)):
                if rng.random() < error_rate:
                    bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
            seq = "".join(bases)
            out.write(f"@sim_{written}_{tid}_{len(blocks)}ex\n{seq}\n+\n{'I' * len(seq)}\n")
            written += 1

    print(f"wrote {written} reads to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
