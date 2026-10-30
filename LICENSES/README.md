The original MIT STAR license is included here for attribution to the original author Alexander Dobin

The LICENSE of this repo has been chosen to also be MIT to match that of the STAR repo https://github.com/alexdobin/STAR

`STAR_RS_LICENSE` is the MIT license of [STAR-rs](https://github.com/BenjaminDEMAILLE/STAR-rs)
(Copyright (c) 2026 Benjamin Demaille), included for attribution of modules ported
from STAR-rs into rustar-aligner. STAR-rs is dual-licensed `MIT OR Apache-2.0`; the
ported code is taken under its MIT option, which is compatible with this repo's MIT
license. See `docs-old/dev/porting-from-star-rs.md` for the porting conventions.

Other third-party code ported into rustar-aligner, each under a licence compatible with MIT:

- `DEGNORM_LICENSE`: [DegNorm](https://github.com/NUStatBioinfo/DegNorm) (MIT,
  Copyright (c) 2018 Northwestern Statistical Bioinformatics Lab). `src/degnorm/` is a Rust
  port of its `degnorm/nmf.py` and the gene selection in `degnorm/__main__.py`.
- `SGT_LICENSE`: David Elworthy's Simple Good-Turing estimator `sgt.h` (BSD 2-clause,
  Copyright (c) David Elworthy 2004), which STAR vendors in `source/SimpleGoodTuring/`.
  `src/solo/sgt.rs` is a Rust port of it. Elworthy's code adapts Geoffrey Sampson's
  Simple Good-Turing program (University of Sussex), which Sampson issued for any use;
  it is acknowledged here as he asks.
- `LIBCXX_LICENSE`: libc++ from the [LLVM Project](https://llvm.org) (Apache License 2.0
  with LLVM Exceptions). `src/solo/libcxx_rng.rs` reimplements in Rust the behaviour of
  libc++'s `std::mt19937`, `std::generate_canonical` and `std::discrete_distribution`, which
  STAR's EmptyDrops_CR draws from, so that the Monte-Carlo draws match STAR's bit for bit.
  Included for attribution.
- `SALMON_LICENSE`: [salmon](https://github.com/COMBINE-lab/salmon) and
  [libradicl](https://github.com/COMBINE-lab/libradicl) (BSD 3-clause, Copyright (c) COMBINE
  lab). `src/io/rad.rs` writes the RAD layout of salmon's `crates/salmon-rad`, so that
  `salmon quant --rad` can read it; no code is linked. Included for attribution.


These cover code ported into rustar-aligner's own source. The licences of the Rust crates
linked into the binaries (and the C libraries some of them compile in) are collected at
build time into `THIRD_PARTY_LICENSES.md`, shipped in the release archives and the Docker
image; generate it locally with `scripts/third-party-licenses.sh`.
