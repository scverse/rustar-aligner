# Contributing to rustar-aligner

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

rustar-aligner was written as a faithful port of [STAR](https://github.com/alexdobin/STAR) by Alexander Dobin. Up to the initial release, the goal was behavioral parity with STAR — matching its algorithms, thresholds, and output formats as closely as possible. Notes from that development phase are in `docs-old/` (`docs-old/dev/` and the `phase*.md` files).

Future development is not bound by that constraint. Adding STARsolo, new features, or diverging from STAR behavior is entirely welcome.

## Documentation site

The user-facing docs site is an [Astro Starlight](https://starlight.astro.build/) project under `docs/`:

```bash
cd docs
pnpm install
pnpm dev          # local dev server
pnpm build        # production build into docs/dist/
```

Content lives under `docs/src/content/docs/` as Markdown / MDX files with YAML frontmatter (`title`, `description`). Sidebar order is configured in `docs/astro.config.mjs`. Site-wide design tokens (colours, fonts, graph-paper background, wave dividers) live in `docs/src/styles/custom.css` and can be tuned in one place.

## License

MIT, matching the original STAR license.
