---
title: Installation
description: Install rustar-aligner from crates.io, the GitHub Container Registry, or build from source.
---

Three ways to get rustar-aligner. Pick whichever fits your workflow:

1. [**`cargo install`** from crates.io](#1-install-from-cratesio) — easiest if you already have Rust.
2. [**Docker image** from the GitHub Container Registry](#2-docker-image-from-ghcr) — good for pipelines and reproducible runs, no Rust toolchain needed.
3. [**Build from source**](#3-build-from-source) — for development or when you need a custom build.

## 1. Install from crates.io

[rustar-aligner is published to crates.io](https://crates.io/crates/rustar-aligner). With Rust installed (via [rustup](https://rustup.rs/)), one command builds and installs the binary:

```bash
cargo install rustar-aligner
```

The binary lands in `~/.cargo/bin/rustar-aligner` (which `rustup` already adds to your `$PATH`).

```bash
rustar-aligner --version
```

To upgrade later, repeat the same command — `cargo install` overwrites in place. To pick a specific version, pass `--version 0.1.0`.

## 2. Docker image from GHCR

Pre-built multi-arch images are published to the [GitHub Container Registry](https://github.com/Psy-Fer/rustar-aligner/pkgs/container/rustar-aligner) for every release. Linux x86_64 and aarch64 are both supported.

```bash
docker pull ghcr.io/psy-fer/rustar-aligner:latest
docker run --rm ghcr.io/psy-fer/rustar-aligner:latest --version
```

### Available tags

| Tag | Meaning |
|-----|---------|
| `latest` | Latest stable release (multi-arch). |
| `0.1.0`, `0.1`, `0` | A specific release, with major / minor aliases. |
| `dev` | Latest commit on `main` (may be unstable). |
| `latest-avx2` / `latest-avx512` / `latest-sve` | SIMD-optimised single-arch builds for hosts that support the relevant instruction set. Use these for top performance on modern x86_64 (AVX2/AVX512) or ARM (SVE) hardware. |

### Aligning with the Docker image

Mount your genome index and read files into the container, and pass paths inside the container to the CLI:

```bash
docker run --rm \
  -v /local/genome_index:/genome \
  -v /local/reads:/reads \
  -v /local/output:/output \
  ghcr.io/psy-fer/rustar-aligner:latest \
    --genomeDir /genome \
    --readFilesIn /reads/sample_1.fq.gz /reads/sample_2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix /output/sample_
```

## 3. Build from source

For development work, or when you want a build tuned for your specific machine, build from source.

### Requirements

- **Rust** 2024 edition (rustc 1.85 or newer). Install via [rustup](https://rustup.rs/) if you don't have it already.
- **Linux**, **macOS**, or **Windows**. CI tests Linux (x86_64, x86-64-v3, aarch64), macOS (aarch64), and Windows (x86_64).
- **Disk space** for the genome index — same requirements as STAR (e.g. ~30 GB for human GRCh38).

### Clone and build

```bash
git clone https://github.com/Psy-Fer/rustar-aligner.git
cd rustar-aligner
cargo build --release
```

The binary lands at `target/release/rustar-aligner`. Add it to your `$PATH` or invoke it directly.

```bash
target/release/rustar-aligner --version
```

### Verify the build

```bash
cargo test            # full test suite (~396 tests)
cargo clippy          # lint (zero warnings expected)
cargo fmt --check     # formatting check
```

### Debug builds

For development you can use `cargo build` (no `--release`). The debug binary is significantly slower but builds faster — useful for iterating on rustar-aligner itself, not for aligning real datasets.

```bash
cargo build
target/debug/rustar-aligner --runMode alignReads ...
```

## What's next

- [Quick start](/rustar-aligner/getting-started/quick-start/) — generate an index and run an alignment
