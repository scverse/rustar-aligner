---
title: Contributing
description: How to build, test, and contribute to rustar-aligner.
---

Contributions are welcome. The repository is on GitHub at [Psy-Fer/rustar-aligner](https://github.com/Psy-Fer/rustar-aligner).

## Building and testing

Rust 2024 edition. Standard Cargo commands:

```bash
cargo build            # debug build
cargo build --release  # release build
cargo test             # run all tests
cargo clippy           # lint (zero warnings expected)
cargo fmt --check      # formatting check
```

CI runs on Linux (x86_64, x86-64-v3, aarch64), macOS (aarch64), and Windows (x86_64). PRs must pass all CI checks before merging.

## Test data

Small synthetic and yeast test data lives in `test/`. Integration tests in `tests/` use the synthetic genome. Differential testing against STAR reference outputs is done via `test/compare_sam.py` and `test/compare_pe.py`.

## Project history

rustar-aligner was written as a faithful port of [STAR](https://github.com/alexdobin/STAR) by Alexander Dobin. Up to the initial release, the goal was behavioural parity with STAR — matching its algorithms, thresholds, and output formats as closely as possible.

Future development is not bound by that constraint. Adding STARsolo, new features, or diverging from STAR behaviour is entirely welcome.

## Working on the website

The website you're reading lives in `docs/`. It's an [Astro Starlight](https://starlight.astro.build/) project.

```bash
cd docs
pnpm install
pnpm dev          # local dev server at http://localhost:4321/rustar-aligner/
pnpm build        # production build into docs/dist/
```

Edit content under `docs/src/content/docs/`. Each page has YAML frontmatter (`title`, `description`) and a markdown body. Sidebar order is configured in `docs/astro.config.mjs`.

## License

MIT, matching the original STAR license.
