#!/usr/bin/env bash
# Parity suite: run rustar-aligner and STAR 2.7.11b on the same public data with the same
# flags, and compare every output file. The README's comparison tables come from this.
#
#   DATA=path/to/testdata STAR=path/to/STAR RUSTAR=target/release/rustar-aligner \
#     test/parity_suite.sh fetch|index|run|compare|all
#
# Needs curl, seqtk, samtools and python3. THREADS (default 8) sets --runThreadN. Files already
# in $DATA are reused, so a slow download (the two ERR12389696 files are 1.6 GB each) can be
# replaced by copies checked against ENA's MD5s.
# Data (all public, pinned):
#   yeast   Ensembl R64-1-1 release 110 + ENA ERR12389696 (seqtk -s100 subsamples, as in
#           CONTRIBUTING.md)
#   chr21   GRCh38 chr21 + GENCODE chr21 GTF from nf-core/test-datasets (modules branch),
#           reads: the first 2,000,000 pairs of GEUVADIS ERR188021 (75 bp PE)
#   fusion  nf-core/test-datasets 40 kb chr22 slice + test_starfusion_rnaseq reads
#   solo    nf-core/test-datasets scrnaseq: mouse chr19 (GRCm38.p6, GENCODE vM19),
#           10x v2 sample S10, 50k-barcode whitelist
set -euo pipefail

: "${DATA:?set DATA}" "${STAR:?set STAR to a STAR 2.7.11b binary}"
RUSTAR=${RUSTAR:-target/release/rustar-aligner}
THREADS=${THREADS:-8}
HERE=$(cd "$(dirname "$0")" && pwd)
NFC=https://raw.githubusercontent.com/nf-core/test-datasets
MODULES=39ab78127391a80748b062614735329dffd351ef
SCRNASEQ=d934d6e8367fe2626184496b1889671cf2b02dab
P=$DATA/parity

# make FILE CMD...: run CMD with stdout to FILE.part, then rename, so an interrupted run never
# leaves a partial FILE behind. A FILE that already exists is kept (you can put copies there).
make() {
  local f=$1; shift
  [ -s "$f" ] && return
  "$@" > "$f.part" && mv "$f.part" "$f"
}
get() { make "$2" curl -fsSL "$1"; }
gz() { curl -fsSL "$1" | gunzip; }
# The first N lines of a gzipped FASTQ, streamed. head closes the pipe early, which makes
# curl and gunzip fail with SIGPIPE, so their status is ignored; the caller checks the count.
head_fq() { { curl -fsSL "$1" | gunzip | head -n "$2" || true; } | gzip; }
subsample() { seqtk sample -s100 "$1" "$2" | gzip; }

fetch() {
  mkdir -p "$DATA"/{reference,reads} "$P"/{chr21,fusion,solo}
  local e=https://ftp.ensembl.org/pub/release-110 y=Saccharomyces_cerevisiae.R64-1-1
  make "$DATA/reference/$y.dna.toplevel.fa" gz "$e/fasta/saccharomyces_cerevisiae/dna/$y.dna.toplevel.fa.gz"
  make "$DATA/reference/$y.110.gtf" gz "$e/gtf/saccharomyces_cerevisiae/$y.110.gtf.gz"
  local r=$DATA/reads/ERR12389696
  for m in 1 2; do
    get "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12389696/ERR12389696_$m.fastq.gz" "${r}_$m.fastq.gz"
    make "${r}_sub_${m}_10k.fastq.gz" subsample "${r}_$m.fastq.gz" 10000
    make "${r}_sub_${m}_100k.fastq.gz" subsample "${r}_$m.fastq.gz" 100000
  done
  local g=data/genomics/homo_sapiens/genome
  get "$NFC/$MODULES/$g/chr21/sequence/genome.fasta" "$P/chr21/genome.fasta"
  get "$NFC/$MODULES/$g/chr21/sequence/chr21_gencode.gtf" "$P/chr21/genes.gtf"
  for m in 1 2; do  # the first 2M pairs (8M lines)
    local h=$P/chr21/ERR188021_2M_$m.fastq.gz
    if [ ! -s "$h" ]; then
      head_fq "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188021/ERR188021_$m.fastq.gz" 8000000 > "$h.part"
      [ "$(gunzip -c "$h.part" | wc -l)" = 8000000 ] && mv "$h.part" "$h" \
        || { echo "ERR188021_$m: download incomplete" >&2; exit 1; }
    fi
    get "$NFC/$MODULES/$g/test_starfusion_rnaseq_$m.fastq.gz" "$P/fusion/test_starfusion_rnaseq_$m.fastq.gz"
  done
  get "$NFC/$MODULES/$g/genome.fasta" "$P/fusion/genome.fasta"
  get "$NFC/$MODULES/$g/genome.gtf" "$P/fusion/genes.gtf"
  get "$NFC/$SCRNASEQ/reference/GRCm38.p6.genome.chr19.fa" "$P/solo/genome.fa"
  get "$NFC/$SCRNASEQ/reference/gencode.vM19.annotation.chr19.gtf" "$P/solo/genes.gtf"
  get "$NFC/$SCRNASEQ/reference/10xv2_whitelist_sub50K.txt" "$P/solo/whitelist.txt"
  for r in R1 R2; do get "$NFC/$SCRNASEQ/testdata/S10_L001_${r}_001.fastq.gz" "$P/solo/S10_L001_${r}_001.fastq.gz"; done
}

TOOLS=(star rustar)
bin() { [ "$1" = star ] && echo "$STAR" || echo "$RUSTAR"; }
timed() { /usr/bin/time -f "%e %M" -o "$1" "${@:2}"; }

index() {  # genome fasta gtf sjdbOverhang SAindexNbases [extra args]
  for t in "${TOOLS[@]}"; do
    local d=$P/idx/${t}_$1; mkdir -p "$d"
    timed "$d.time" "$(bin $t)" --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$d" \
      --genomeFastaFiles "$2" --sjdbGTFfile "$3" --sjdbOverhang "$4" --genomeSAindexNbases "$5" \
      --outFileNamePrefix "$d.log_" "${@:6}" > "$d.stdout" 2>&1
  done
}

indices() {
  local y=$DATA/reference/Saccharomyces_cerevisiae.R64-1-1
  index yeast "$y.dna.toplevel.fa" "$y.110.gtf" 149 10
  index yeast_sa2 "$y.dna.toplevel.fa" "$y.110.gtf" 149 10 --genomeSAsparseD 2
  index chr21 "$P/chr21/genome.fasta" "$P/chr21/genes.gtf" 74 11
  index chr22 "$P/fusion/genome.fasta" "$P/fusion/genes.gtf" 72 8
  index mm19 "$P/solo/genome.fa" "$P/solo/genes.gtf" 89 11
}

align() {  # name genome args...
  local name=$1 g=$2; shift 2
  for t in "${TOOLS[@]}"; do
    local o=$P/out/${t}_$name/; mkdir -p "$o"
    timed "$o/time.txt" "$(bin $t)" --runThreadN "$THREADS" --genomeDir "$P/idx/${t}_$g" \
      --readFilesCommand zcat --outFileNamePrefix "$o" "$@" > "$o/stdout.txt" 2>&1 \
      || echo "FAIL $t $name" >&2
  done
}

CONFIGS=(y10k_se y10k_pe y100k_se y100k_pe y100k_pe_quant y100k_pe_2pass y100k_pe_bysj
  y100k_pe_attr y100k_pe_sa2 y100k_se_chim y100k_pe_chim0 y100k_pe_chim10 c21_se c21_pe c21_pe_quant
  c21_pe_chim0 c21_pe_chim10 fus_pe_chim0 fus_pe_chim10 solo_def solo_cr)

run() {
  local r=$DATA/reads/ERR12389696_sub h=$P/chr21/ERR188021_2M f=$P/fusion/test_starfusion_rnaseq
  local chim=(--chimSegmentMin 12 --chimOutType Junctions)
  local ype=(--readFilesIn "${r}_1_100k.fastq.gz" "${r}_2_100k.fastq.gz")
  local hpe=(--readFilesIn "${h}_1.fastq.gz" "${h}_2.fastq.gz")
  local fpe=(--readFilesIn "${f}_1.fastq.gz" "${f}_2.fastq.gz")
  align y10k_se yeast --readFilesIn "${r}_1_10k.fastq.gz"
  align y10k_pe yeast --readFilesIn "${r}_1_10k.fastq.gz" "${r}_2_10k.fastq.gz"
  align y100k_se yeast --readFilesIn "${r}_1_100k.fastq.gz"
  align y100k_pe yeast "${ype[@]}"
  align y100k_pe_quant yeast "${ype[@]}" --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM SortedByCoordinate
  align y100k_pe_2pass yeast "${ype[@]}" --twopassMode Basic
  align y100k_pe_bysj yeast "${ype[@]}" --outFilterType BySJout
  align y100k_pe_attr yeast "${ype[@]}" --outSAMattributes All --outSAMstrandField intronMotif \
    --outSAMunmapped Within --outSAMtype BAM Unsorted
  align y100k_pe_sa2 yeast_sa2 "${ype[@]}"
  align y100k_se_chim yeast --readFilesIn "${r}_1_100k.fastq.gz" "${chim[@]}"
  align y100k_pe_chim0 yeast "${ype[@]}" "${chim[@]}"
  align y100k_pe_chim10 yeast "${ype[@]}" "${chim[@]}" --chimMultimapNmax 10
  align c21_se chr21 --readFilesIn "${h}_1.fastq.gz"
  align c21_pe chr21 "${hpe[@]}"
  align c21_pe_quant chr21 "${hpe[@]}" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted
  align c21_pe_chim0 chr21 "${hpe[@]}" "${chim[@]}" --outSAMtype None
  align c21_pe_chim10 chr21 "${hpe[@]}" "${chim[@]}" --chimMultimapNmax 10 --outSAMtype None
  align fus_pe_chim0 chr22 "${fpe[@]}" "${chim[@]}"
  align fus_pe_chim10 chr22 "${fpe[@]}" "${chim[@]}" --chimMultimapNmax 10
  local solo=(--readFilesIn "$P/solo/S10_L001_R2_001.fastq.gz" "$P/solo/S10_L001_R1_001.fastq.gz"
    --soloType CB_UMI_Simple --soloCBwhitelist "$P/solo/whitelist.txt"
    --soloFeatures Gene GeneFull SJ Velocyto)
  align solo_def mm19 "${solo[@]}" --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
  align solo_cr mm19 "${solo[@]}" --outSAMtype None --soloCellFilter EmptyDrops_CR \
    --soloMultiMappers Uniform EM --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
}

compare() {
  for g in yeast yeast_sa2 chr21 chr22 mm19; do
    python3 "$HERE/parity_compare.py" index "$P/idx/star_$g" "$P/idx/rustar_$g" "$g"
  done
  for c in "${CONFIGS[@]}"; do
    python3 "$HERE/parity_compare.py" align "$P/out/star_$c" "$P/out/rustar_$c" "$c"
  done
}

case ${1:-all} in
  fetch) fetch ;; index) indices ;; run) run ;; compare) compare ;;
  all) fetch; indices; run; compare ;;
  *) echo "usage: $0 fetch|index|run|compare|all" >&2; exit 2 ;;
esac
