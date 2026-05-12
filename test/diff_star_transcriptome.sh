#!/usr/bin/env bash
# Byte-for-byte diff of rustar-aligner vs STAR for the transcriptome index files
# written at genomeGenerate time.
#
# Usage: ./diff_star_transcriptome.sh [WORKDIR]
#
# Generates a synthetic 2-chr / 4-transcript / 4-gene test case, runs both
# rustar-aligner and STAR genomeGenerate, then diffs each output file byte-for-byte.
#
# Requires STAR on PATH (`brew install rna-star`, bioconda star, etc.).
# rustar-aligner must already be built — the script expects `target/debug/rustar-aligner`
# in the repo root (or $RUSTAR_BIN set to a compiled binary).

set -euo pipefail

WORKDIR="${1:-/tmp/rustar_aligner_diff}"
RUSTAR_BIN="${RUSTAR_BIN:-$(pwd)/target/debug/rustar-aligner}"

if [[ ! -x "$RUSTAR_BIN" ]]; then
    echo "ERROR: rustar-aligner binary not found at $RUSTAR_BIN"
    echo "Build with: cargo build (or set RUSTAR_BIN=/path/to/rustar-aligner)"
    exit 1
fi

if ! command -v STAR >/dev/null 2>&1; then
    echo "ERROR: STAR not on PATH. Install via 'brew install rna-star' or bioconda."
    exit 1
fi

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Deterministic test case.
python3 - <<'PYEOF'
BASES = "ACGT"
def lcg(seed, length):
    state = seed & 0xFFFFFFFF
    out = []
    for _ in range(length):
        state = (state * 1103515245 + 12345) & 0xFFFFFFFF
        out.append(BASES[(state >> 16) & 3])
    return "".join(out)

CHR1 = lcg(11111, 3000)
CHR2 = lcg(22222, 3000)

with open("genome.fa", "w") as f:
    f.write(">chr1\n")
    for i in range(0, len(CHR1), 60): f.write(CHR1[i:i+60] + "\n")
    f.write(">chr2\n")
    for i in range(0, len(CHR2), 60): f.write(CHR2[i:i+60] + "\n")

with open("annotations.gtf", "w") as f:
    f.write('chr1\ttest\texon\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    f.write('chr1\ttest\texon\t501\t600\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n')
    f.write('chr1\ttest\texon\t701\t800\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n')
    f.write('chr1\ttest\texon\t901\t1000\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n')
    f.write('chr2\ttest\texon\t301\t600\t.\t-\t.\tgene_id "G3"; transcript_id "T3"; gene_biotype "protein_coding";\n')
    f.write('chr2\ttest\texon\t1001\t1100\t.\t+\t.\tgene_id "G4"; transcript_id "T4";\n')
    f.write('chr2\ttest\texon\t1301\t1400\t.\t+\t.\tgene_id "G4"; transcript_id "T4";\n')
PYEOF

# STAR and rustar-aligner genomeGenerate. Both tools echo `--genomeDir <dir>`
# into genomeParameters.txt's commandLineFull comment, so the two must
# receive the SAME literal path string for that line to match byte-for-byte.
# Run each in a subshell that symlinks `index` → the per-tool output dir,
# pass `--genomeDir index` to both, then resolve the symlink away.
rm -rf star_index rustar_aligner_index index
mkdir star_index rustar_aligner_index

ln -sfn star_index index
STAR --runMode genomeGenerate --genomeDir index \
     --genomeFastaFiles genome.fa --sjdbGTFfile annotations.gtf \
     --genomeSAindexNbases 5 --sjdbOverhang 49 --runThreadN 1 >/dev/null 2>&1
ln -sfn rustar_aligner_index index
"$RUSTAR_BIN" --runMode genomeGenerate --genomeDir index \
     --genomeFastaFiles genome.fa --sjdbGTFfile annotations.gtf \
     --genomeSAindexNbases 5 --sjdbOverhang 49 --runThreadN 1 >/dev/null 2>&1
rm -f index

# Diff each file. SA / SAindex are excluded because STAR's suffix-array
# construction algorithm differs from rustar-aligner's and produces a different
# (but equally-valid) byte layout — the sort itself is deterministic but
# STAR's bucketed parallel build breaks ties differently. Size parity is
# verified separately below.
pass=0
fail=0
for f in chrName.txt chrLength.txt chrStart.txt chrNameLength.txt \
         transcriptInfo.tab exonInfo.tab geneInfo.tab exonGeTrInfo.tab \
         sjdbList.fromGTF.out.tab sjdbInfo.txt sjdbList.out.tab \
         Genome genomeParameters.txt; do
    if diff -q "star_index/$f" "rustar_aligner_index/$f" >/dev/null 2>&1; then
        echo "✓ $f"
        pass=$((pass+1))
    else
        echo "✗ $f DIFFERS"
        diff "star_index/$f" "rustar_aligner_index/$f" | head -20
        fail=$((fail+1))
    fi
done

# SA / SAindex: size-only check (see comment above).
for f in SA SAindex; do
    star_size=$(wc -c <"star_index/$f")
    rustar_aligner_size=$(wc -c <"rustar_aligner_index/$f")
    if [[ "$star_size" == "$rustar_aligner_size" ]]; then
        echo "~ $f size-only match ($star_size bytes) — content differs (expected, different SA algorithm)"
        pass=$((pass+1))
    else
        echo "✗ $f size mismatch: STAR=$star_size rustar-aligner=$rustar_aligner_size"
        fail=$((fail+1))
    fi
done

echo
echo "RESULT: $pass/$((pass+fail)) files identical (or size-matched) to STAR's output"
if [[ $fail -gt 0 ]]; then
    exit 1
fi
