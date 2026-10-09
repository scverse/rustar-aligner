# Contributing to rustar-aligner

## Building and testing

Rust 2024 edition. Standard Cargo commands:

```bash
cargo build                 # debug build
cargo build --release       # release build
cargo test                  # run all tests
cargo clippy --all-targets  # lint (zero warnings expected)
cargo fmt --check           # formatting check
```

CI runs on Linux (x86_64, x86-64-v3, aarch64), macOS (aarch64), and Windows (x86_64). PRs must pass all CI checks before merging.

A green `cargo test` is necessary but **not sufficient** — it is the floor, not the gate. See [Validating alignment, counting, and output changes](#validating-alignment-counting-and-output-changes) below.

## Pull request standards

These are hard requirements. A PR that does not meet them will be sent back rather than reviewed in detail.

### You are the author

Use whatever tools you like, including AI assistants — but **you** are the author, and you are accountable for every line and every claim in the PR. Before you open it, you must have read and understood the whole diff and be able to explain and defend it yourself.

When a reviewer asks a question, answer it yourself, in your own words, from your own understanding of the change. You are welcome to go back to a tool to dig up specifics — file locations, exact numbers, a STAR source reference — and include them as supporting detail. What you should not do is make pasted tool output the substance of your reply: the review is a conversation between the contributor and the maintainer, not a relay for an assistant the contributor hasn't digested.

### The description must match the code

The PR description, CHANGELOG entry, and in-code comments are part of the contribution and are held to the same standard as the code:

- **Every test named in the description must exist and pass.** Do not list a test that isn't in the diff.
- **Every parameter, flag, or mechanism the description refers to must be present in the diff.** Do not describe behaviour the code does not implement.
- **Do not describe dead code as if it runs.** If a function is not reached by any production code path, either wire it in or leave it out — don't ship it with a description implying it is active.
- **Every benchmark number must be reproducible from the branch as submitted**, using a command anyone can run. Numbers produced by a locally-patched binary that isn't in the PR are not acceptable.
- **Claims about what STAR does must be checked against the STAR source**, not asserted from memory. If a comment says "this is what STAR does", it must be verifiable in STAR's C++. A comment that misstates STAR's behaviour is a defect even if the code around it is fine.

Inaccurate descriptions waste reviewer time and erode trust in the whole PR. If we find one fabricated claim, we will assume the rest of the description is unverified.

### One theme per PR

- One logical change per PR. Do not bundle unrelated features (a solo counting fix and a new RNG module; an index builder and two new genome types).
- **No hidden changes.** If the diff introduces a new algorithm, a new dependency, or a change to a shared/hot code path, it must be named in the description — reviewers should never discover a 250-line stitcher the summary didn't mention.
- No dead code. A module that changes no behaviour yet does not belong in a behaviour PR; land it with the change that uses it, behind its own validation.
- Keep changes to shared files (e.g. `src/params/mod.rs`) minimal and additive, since most PRs touch them and will otherwise conflict.

**Exception for maintainers landing their own reviewed work.** A maintainer pushing changes they have already reviewed may bundle related work into one PR to `main`, split cleanly by commit — the PR is the push, not the review, because the review already happened. This exception belongs to whoever carries that review responsibility; it is not a general licence to bundle. If you are contributing from a fork, the one-theme rule applies. A maintainer may bundle multiple PRs together when it would make sense to do so during review, and then submit a single bundled PR, closing out the included PRs. These bundled PRs should be linked to ensure tracability within github.

### Divergence from STAR is allowed — but must be deliberate and flagged

Diverging from STAR (including adding non-STAR flags, or choosing STAR's *documented* behaviour over its *actual binary* behaviour where they differ) is welcome, but it must be an explicit, signed-off decision — never an accident and never presented as faithfulness:

- State the divergence plainly in the PR and add an entry to [`DIVERGENCE.md`](DIVERGENCE.md) with the rationale and the STAR behaviour it departs from.
- Do not label a divergence "faithful", and do not invent a STAR flag/behaviour that does not exist and present it as parity.
- A non-STAR flag or a change to default output behaviour needs maintainer sign-off before merge.

### When a STAR-faithful change makes the numbers worse

A change that brings an algorithm closer to STAR's can move the benchmark the
wrong way. That is not evidence the change was wrong; it usually means the
surrounding logic has not caught up yet. The question to ask is "what else does
STAR do here that we do not", not "which approximation scored better".

Keeping an approximation because it currently measures better is how a port
stops converging. Say so in the PR when this happens, with the before and
after, so the regression is a recorded step rather than a surprise.

### Code conventions

- The crate and binary are `rustar-aligner`; the Rust library name is
  `rustar_aligner`. `#![allow(non_snake_case)]` in `lib.rs` is deliberate, so
  STAR's camelCase identifiers can be kept where they help.
- STAR parameters keep their `--camelCase` spelling through clap's
  `#[arg(long = "camelCase")]`, mapping to snake_case Rust fields.
- Multi-value parameters need an explicit `num_args`, and parameters with
  negative defaults need `allow_hyphen_values = true`.
- Validation that clap cannot express lives in `Parameters::validate()`.
- No async. The work is CPU-bound; async would add machinery and buy nothing.
- `thiserror` for the library's `Error` enum, `anyhow` at the top level.

### New dependencies need prior discussion

Adding a dependency — **especially a non-Rust one** (a C library via a `-sys` crate, anything needing `bindgen`/`libclang` or a system library) — must be raised in an issue *before* the PR. This project is published to crates.io and builds on five platforms including Windows; a new C dependency is a maintenance and supply-chain decision, not an implementation detail.

### Accepted-but-inert parameters

If you add a CLI flag that parses but is not yet implemented, mark it as such in the parameter-surface test and document it — do not silently accept a flag that does nothing. A user passing a flag should never be quietly ignored.

## Test data

Integration tests in `tests/` use a bundled synthetic micro-genome and need no downloads. The differential benchmark below uses a small **public** yeast RNA-seq dataset that is not vendored; fetch it once and point `DATA` at wherever you keep it.

- **Reference genome + annotation:** *Saccharomyces cerevisiae* R64-1-1, [Ensembl release 110](https://ftp.ensembl.org/pub/release-110/).
- **Reads:** ENA run [ERR12389696](https://www.ebi.ac.uk/ena/browser/view/ERR12389696) — paired-end 150 bp yeast RNA-seq.

Requires `seqtk` and a **STAR 2.7.11b** binary (built from <https://github.com/alexdobin/STAR>) on `PATH`.

```bash
export DATA=path/to/testdata          # your choice; used by every command below
mkdir -p "$DATA"/{reference,reads}

# 1. Reference genome + GTF
cd "$DATA/reference"
wget https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
gunzip *.gz

# 2. Full read pair from ENA
cd "$DATA/reads"
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12389696/ERR12389696_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12389696/ERR12389696_2.fastq.gz

# 3. Deterministic subsamples. The SAME seed (-s100) on both mates keeps pairs aligned.
#    The benchmark uses the 10k tier; test/run_tests.sh also uses the 100 and 1000 tiers.
for tier in "100:100" "1000:1000" "10k:10000"; do
  name=${tier%%:*}; n=${tier##*:}
  seqtk sample -s100 ERR12389696_1.fastq.gz "$n" | gzip > "ERR12389696_sub_1_${name}.fastq.gz"
  seqtk sample -s100 ERR12389696_2.fastq.gz "$n" | gzip > "ERR12389696_sub_2_${name}.fastq.gz"
done
cd -
```

Build both indices from that reference — STAR's for the reference alignments, and rustar-aligner's own. **Never reuse a STAR index for rustar-aligner**; each tool builds its own.

```bash
STAR --runMode genomeGenerate --genomeDir "$DATA/indices_star" \
  --genomeFastaFiles "$DATA/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa" \
  --sjdbGTFfile "$DATA/reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf" \
  --sjdbOverhang 149 --genomeSAindexNbases 10 --runThreadN 4

./target/release/rustar-aligner --runMode genomeGenerate --genomeDir "$DATA/indices_rustar" \
  --genomeFastaFiles "$DATA/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa" \
  --sjdbGTFfile "$DATA/reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf" \
  --sjdbOverhang 149 --genomeSAindexNbases 10 --runThreadN 4
```

## Validating alignment, counting, and output changes

`cargo test` does not catch faithfulness regressions. Anything that touches the aligner, counting, or record output **must** run the relevant differential harness against a reference produced by **STAR 2.7.11b** (<https://github.com/alexdobin/STAR>) and report before/after numbers in the PR. The dataset and both indices come from [Test data](#test-data) above (`$DATA`).

**SE and PE alignment** (10k yeast reads) — this is the project's defining metric; a core-aligner PR without these numbers is incomplete.

```bash
# STAR reference alignments (the SAM the comparison diffs against)
STAR --genomeDir "$DATA/indices_star" --readFilesCommand zcat --outSAMtype SAM --runThreadN 1 \
  --readFilesIn "$DATA/reads/ERR12389696_sub_1_10k.fastq.gz" \
  --outFileNamePrefix "$DATA/star_10k_/"
STAR --genomeDir "$DATA/indices_star" --readFilesCommand zcat --outSAMtype SAM --runThreadN 1 \
  --readFilesIn "$DATA/reads/ERR12389696_sub_1_10k.fastq.gz" "$DATA/reads/ERR12389696_sub_2_10k.fastq.gz" \
  --outFileNamePrefix "$DATA/star_10k_pe_/"

# rustar-aligner alignments (its own index)
./target/release/rustar-aligner --runMode alignReads --genomeDir "$DATA/indices_rustar" \
  --readFilesCommand zcat --outSAMtype SAM --runThreadN 1 \
  --readFilesIn "$DATA/reads/ERR12389696_sub_1_10k.fastq.gz" \
  --outFileNamePrefix "$DATA/rustar_10k_/"
./target/release/rustar-aligner --runMode alignReads --genomeDir "$DATA/indices_rustar" \
  --readFilesCommand zcat --outSAMtype SAM --runThreadN 1 \
  --readFilesIn "$DATA/reads/ERR12389696_sub_1_10k.fastq.gz" "$DATA/reads/ERR12389696_sub_2_10k.fastq.gz" \
  --outFileNamePrefix "$DATA/rustar_10k_pe_/"

# compare. NOTE the arg styles differ: compare_sam.py is NAMED, compare_pe.py is POSITIONAL.
python3 test/compare_sam.py \
  --rustar-aligner "$DATA/rustar_10k_/Aligned.out.sam" \
  --star           "$DATA/star_10k_/Aligned.out.sam"
python3 test/compare_pe.py "$DATA/rustar_10k_pe_/Aligned.out.sam" "$DATA/star_10k_pe_/Aligned.out.sam"
```

Baseline to beat: SE ~99.8% and PE ~99.9% tie-adjusted faithfulness (see `ROADMAP.md`). Faithfulness must not regress; if a STAR-faithful change regresses a metric, the fix is more STAR-matching work, not a revert (a deliberate, signed-off divergence is the documented exception).

- **STARsolo** counting changes: run `test/solo_diff_docker.sh` against a real STARsolo oracle. A change to default barcode/UMI counting is not validated by unit tests alone.
- **STARlong / long-read** changes: validate against real STARlong output, not just synthetic unit fixtures.
- **genomeGenerate / index** changes: correctness first — a suffix array is unique for a given text, so the builder's output must be shown **byte-identical** to the existing shipped path (the code path that actually ships, not a function only reached by unit tests). Any **speed or memory claim** must be measured across **multiple genome sizes** — at minimum yeast (~12 Mb), a chromosome-scale genome, and a mammalian-scale genome — with a **fair comparison**: same machine, same `--runThreadN`, same input, cold vs warm cache stated, and peak RSS measured the same way for both builders. One number on one genome is not evidence of a general improvement.

Report the harness used and the raw before/after counts in the PR. "Tests pass" is not a validation result.

## Project history

rustar-aligner was written as a faithful port of [STAR](https://github.com/alexdobin/STAR) by Alexander Dobin. Up to the initial release, the goal was behavioral parity with STAR — matching its algorithms, thresholds, and output formats as closely as possible. Notes from that development phase are in `docs-old/` (`docs-old/dev/` and the `phase*.md` files).

Future development is not bound by that constraint. Adding STARsolo, new features, or diverging from STAR behavior is entirely welcome AFTER a parity version has been released. Additions that don't impact parity are welcome before that point.

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
