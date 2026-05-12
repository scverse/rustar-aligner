# Project structure
mkdir -p data/{small,medium,large}/{reference,reads,indices,outputs/{star,rustar-aligner}}

# small dataset

# 1. Download yeast reference genome and annotation
cd data/small/reference
wget https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
gunzip *.gz
cd ..

# 2. Download ERR12389696 reads
cd reads

# Using ENA's FTP (faster, recommended)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12389696/ERR12389696_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12389696/ERR12389696_2.fastq.gz

# Create subsampled version for quick testing (100K read pairs)
seqtk sample -s100 ERR12389696_1.fastq.gz 100000 | gzip > ERR12389696_sub_1.fastq.gz
seqtk sample -s100 ERR12389696_2.fastq.gz 100000 | gzip > ERR12389696_sub_2.fastq.gz

cd ..

# 3. Build STAR index (shared between STAR and rustar-aligner)
STAR --runMode genomeGenerate \
     --genomeDir indices \
     --genomeFastaFiles reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     --sjdbGTFfile reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf \
     --sjdbOverhang 149 \
     --genomeSAindexNbases 10 \
     --runThreadN 4

# 4. Quick test alignment with STAR (subsampled data)
STAR --genomeDir indices \
     --readFilesIn reads/ERR12389696_sub_1.fastq.gz reads/ERR12389696_sub_2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix outputs/star/quick_test_ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --runThreadN 4

# parameter explanations
# --sjdbOverhang 149: For 150bp reads (readLength - 1)
# --genomeSAindexNbases 10: Reduced from default 14 for yeast's small genome
# --quantMode GeneCounts: Generates gene count tables for validation