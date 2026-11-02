---
title: Streaming input
description: Align FASTQ read from a pipe or a download, without writing it to disk.
---

`--readFilesIn` accepts anything that can be opened and read from start to end: FIFOs, process substitution and `/dev/stdin`, plain or compressed. The reads are never written to disk and the output is the same as from a file.

```bash
rustar-aligner --runMode alignReads --genomeDir /path/to/genome_index \
  --readFilesIn <(curl -sSf https://example.org/r1.fastq.gz) \
                <(curl -sSf https://example.org/r2.fastq.gz)
```

Comma-separated lists work too. Files of a list are opened in turn, so a stream later in the list waits, with its pipe full, until the earlier ones are done. A download may time out while it waits.

## Two-pass mode needs files

`--twopassMode Basic` reads the input twice, and a stream can only be read once, so it stops at once with an error. Save the reads to a regular file first, or use `--twopassMode None`. See `DIVERGENCE.md` for how this differs from STAR.

## Checking that the whole input arrived

The exit status of a process feeding `<(...)` does not reach the aligner: when `curl` fails, the aligner only sees its input end. What it can detect:

| The stream stops... | Result |
| --- | --- |
| inside a gzip member (or a bzip2, zstd or xz stream) | error |
| in the middle of a plain-text read | error |
| in one mate of a pair, the other going on | error |
| after a `--readFilesCommand` that exits non-zero | error |
| **between two plain-text reads** | **exit 0, fewer reads** |
| **between two gzip members** (BGZF and concatenated gzip files), or between frames or streams of the other formats | **exit 0, fewer reads** |
| **before its first byte** (a failed request, for one, gives an empty stream) | **exit 0, 0 reads** |

The last three leave an input that is valid and shorter, in both mates if they stop together, and the aligner cannot tell it from a smaller file. The only protection is to compare `Number of input reads` in `Log.final.out` with the number you expect, for example from the read count that the archive publishes.
