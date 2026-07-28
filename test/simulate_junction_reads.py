#!/usr/bin/env python3
"""Simulate junction-spanning reads from a genome + GTF.

The existing test tiers all use an unannotated index, which makes the
annotated-junction code paths (`sjA` shortcut, annotated-boundary snapping)
invisible: they simply never run. This generates reads that are guaranteed to
cross annotated exon-exon junctions, so those paths are exercised and any
divergence from STAR shows up as a junction count or CIGAR difference.

Reads are drawn from spliced transcript sequences, centred on junctions, with a
fixed seed so the FASTQ is reproducible byte-for-byte.

Usage:
    simulate_junction_reads.py <genome.fa> <annotation.gtf> <out.fastq> [n_reads] [read_len]
"""

import gzip
import random
import sys
from collections import defaultdict


def read_fasta(path):
    """Return {contig_name: sequence}. Names are truncated at first whitespace."""
    opener = gzip.open if path.endswith(".gz") else open
    seqs, name, chunks = {}, None, []
    with opener(path, "rt") as fh:
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
    """Return {transcript_id: (contig, strand, [(start, end), ...])}, 1-based inclusive."""
    opener = gzip.open if path.endswith(".gz") else open
    exons = defaultdict(list)
    meta = {}
    with opener(path, "rt") as fh:
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
    out = {}
    for tid, blocks in exons.items():
        blocks.sort()
        contig, strand = meta[tid]
        out[tid] = (contig, strand, blocks)
    return out


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    genome_path, gtf_path, out_path = sys.argv[1:4]
    n_reads = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    read_len = int(sys.argv[5]) if len(sys.argv) > 5 else 100

    genome = read_fasta(genome_path)
    transcripts = read_exons(gtf_path)

    # Collect every (transcript, junction) pair with enough flank on both sides
    # to host a read that actually crosses the junction.
    candidates = []
    for tid, (contig, strand, blocks) in transcripts.items():
        if contig not in genome or len(blocks) < 2:
            continue
        # Offset of each junction within the spliced transcript.
        offset = 0
        for k in range(len(blocks) - 1):
            offset += blocks[k][1] - blocks[k][0] + 1
            candidates.append((tid, offset))

    if not candidates:
        sys.exit("no multi-exon transcripts found; nothing to simulate")

    rng = random.Random(20260728)
    rng.shuffle(candidates)

    spliced_cache = {}

    def spliced(tid):
        if tid not in spliced_cache:
            contig, strand, blocks = transcripts[tid]
            seq = "".join(genome[contig][s - 1 : e] for s, e in blocks)
            spliced_cache[tid] = (seq, strand)
        return spliced_cache[tid]

    # Yeast has only a few hundred introns, so one read per junction is not
    # enough to be interesting. Cycle the junction list, varying the overhang
    # each pass, until the requested count is reached or a pass places nothing.
    written = 0
    with open(out_path, "w") as out:
        while written < n_reads:
            placed_this_pass = 0
            for tid, junction_offset in candidates:
                if written >= n_reads:
                    break
                seq, _strand = spliced(tid)
                if len(seq) < read_len:
                    continue
                # Put the junction in the middle third of the read so both
                # overhangs comfortably clear alignSJDBoverhangMin.
                overhang = rng.randint(read_len // 3, 2 * read_len // 3)
                start = junction_offset - overhang
                if start < 0 or start + read_len > len(seq):
                    continue
                read = seq[start : start + read_len].upper()
                if "N" in read:
                    continue
                # Half the reads on the opposite strand, so both orientations
                # of the annotated-junction path get exercised.
                if written % 2 == 1:
                    read = revcomp(read)
                out.write(f"@sim_{written}_{tid}_j{junction_offset}\n")
                out.write(read + "\n+\n" + "I" * read_len + "\n")
                written += 1
                placed_this_pass += 1
            if placed_this_pass == 0:
                break

    print(f"wrote {written} junction-spanning reads to {out_path}", file=sys.stderr)
    if written < n_reads:
        print(
            f"note: only {written} of {n_reads} requested reads were placeable",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
