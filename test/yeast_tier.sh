#!/usr/bin/env bash
# Per-read agreement against STAR on the project's yeast tier.
#
# The full benchmark in CONTRIBUTING.md wants the whole ERR12389696 run; this
# fetches the first ~120 MB of each mate instead, which is 10 000 pairs and
# enough to see a faithfulness regression in the numbers that matter (position
# agreement, unique/multi counts, NH depth). It runs in about a minute.
#
#   test/yeast_tier.sh                       # build, run, compare
#   DATA=~/yeast-tier test/yeast_tier.sh     # keep the downloads somewhere
#
# Requires STAR 2.7.11b on PATH.
set -euo pipefail

DATA="${DATA:-/tmp/rustar-yeast-tier}"
RUSTAR="${RUSTAR:-./target/release/rustar-aligner}"
PAIRS="${PAIRS:-10000}"
THREADS="${THREADS:-4}"

mkdir -p "$DATA"
cd "$DATA"

if [ ! -f genome.fa ]; then
  echo "fetching the yeast reference"
  curl -sfL "https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" -o genome.fa.gz
  gunzip -kf genome.fa.gz
fi

lines=$((PAIRS * 4))
for m in 1 2; do
  if [ ! -f "r${m}.fq" ]; then
    echo "fetching mate ${m} (partial)"
    # A byte range rather than the whole run: the reads are ordered, so the
    # first N records of each mate still pair up.
    curl -sfL -r 0-120000000 \
      "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12389696/ERR12389696_${m}.fastq.gz" \
      -o "r${m}_part.fastq.gz"
    gunzip -c "r${m}_part.fastq.gz" 2>/dev/null | head -"${lines}" > "r${m}.fq" || true
  fi
done

for pair in "rustar:$OLDPWD/$RUSTAR:ridx" "star:STAR:sidx"; do
  name="${pair%%:*}"; rest="${pair#*:}"; exe="${rest%%:*}"; idx="${rest##*:}"
  if [ ! -f "${idx}/SA" ]; then
    echo "building the ${name} index"
    mkdir -p "$idx"
    "$exe" --runMode genomeGenerate --genomeDir "$idx" \
      --genomeFastaFiles genome.fa --genomeSAindexNbases 11 \
      --outFileNamePrefix "${name}_idx_" > /dev/null
  fi
done

echo "aligning with rustar-aligner"
"$OLDPWD/$RUSTAR" --runMode alignReads --genomeDir ridx --readFilesIn r1.fq r2.fq \
  --outSAMtype SAM --runThreadN "$THREADS" --outFileNamePrefix rustar_ > /dev/null

echo "aligning with STAR"
STAR --genomeDir sidx --readFilesIn r1.fq r2.fq \
  --outSAMtype SAM --runThreadN "$THREADS" --outFileNamePrefix star_ > /dev/null

echo
for f in star rustar; do
  printf '%-7s ' "$f"
  grep -E "Uniquely mapped reads number|mapped to multiple loci \|" "${f}_Log.final.out" \
    | sed 's/^ *//' | tr '\n' ' '
  echo
done
echo
python3 "$OLDPWD/test/sam_agreement.py" star_Aligned.out.sam rustar_Aligned.out.sam
