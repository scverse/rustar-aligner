pub mod fasta;

use std::path::Path;

use crate::error::Error;
use crate::params::Parameters;

use fasta::parse_fasta_files;

/// STAR's genome spacing character (used for inter-chromosome padding).
const GENOME_SPACING_CHAR: u8 = 5;

/// Packed genome with chromosome metadata.
///
/// The genome sequence is stored as one byte per base:
/// - A=0, C=1, G=2, T=3, N=4, padding=5
/// - Chromosomes are concatenated with padding to bin boundaries
/// - The reverse complement occupies the second half of the `sequence` buffer
#[derive(Clone)]
pub struct Genome {
    /// Forward genome (0..n_genome) + reverse complement (n_genome..2*n_genome).
    /// Initialized to GENOME_SPACING_CHAR (5), then overwritten with actual bases.
    pub sequence: Vec<u8>,

    /// Total length of the forward (padded) genome.
    pub n_genome: u64,

    /// Number of real chromosomes (not including scaffold/contigs if excluded).
    pub n_chr_real: usize,

    /// Chromosome names.
    pub chr_name: Vec<String>,

    /// True (unpadded) chromosome lengths.
    pub chr_length: Vec<u64>,

    /// Padded start positions of each chromosome in the genome.
    /// Length = n_chr_real + 1; the last entry is n_genome (total size).
    pub chr_start: Vec<u64>,
}

impl Genome {
    /// Build a genome from FASTA files, matching STAR's layout.
    ///
    /// # Arguments
    /// - `params`: CLI parameters (genomeFastaFiles, genomeChrBinNbits)
    ///
    /// # Returns
    /// A `Genome` with forward + reverse complement sequences and metadata.
    pub fn from_fasta(params: &Parameters) -> Result<Self, Error> {
        let chromosomes = parse_fasta_files(&params.genome_fasta_files)?;

        // Compute padding bin size
        let bin_nbits = params.genome_chr_bin_nbits;
        let bin_size = 1u64 << bin_nbits;

        // First pass: compute padded positions and total genome size
        let mut chr_name = Vec::new();
        let mut chr_length = Vec::new();
        let mut chr_start = Vec::new();

        let mut n: u64 = 0; // current position in the padded genome

        for chrom in &chromosomes {
            let len = chrom.sequence.len() as u64;

            if len == 0 {
                return Err(Error::Fasta(format!(
                    "chromosome '{}' has zero length",
                    chrom.name
                )));
            }

            // Apply STAR's padding formula before this chromosome (except for the first)
            if n > 0 {
                n = ((n + 1) / bin_size + 1) * bin_size;
            }

            chr_name.push(chrom.name.clone());
            chr_length.push(len);
            chr_start.push(n);

            n += len;
        }

        // Final padding after the last chromosome
        n = ((n + 1) / bin_size + 1) * bin_size;
        let n_genome = n;
        chr_start.push(n_genome); // STAR adds this extra entry

        let n_chr_real = chromosomes.len();

        // Allocate buffer: forward (0..n_genome) + reverse (n_genome..2*n_genome)
        let total_len = (n_genome * 2) as usize;
        let mut sequence = vec![GENOME_SPACING_CHAR; total_len];

        // Second pass: copy actual chromosome sequences into the buffer
        for (i, chrom) in chromosomes.iter().enumerate() {
            let start = chr_start[i] as usize;
            let len = chrom.sequence.len();
            sequence[start..start + len].copy_from_slice(&chrom.sequence);
        }

        // Build reverse complement in the second half
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        Ok(Genome {
            sequence,
            n_genome,
            n_chr_real,
            chr_name,
            chr_length,
            chr_start,
        })
    }

    /// Append a splice-junction flanking-sequence buffer (`Gsj`) to the
    /// forward genome and rebuild the reverse complement over the extended
    /// forward range. Matches STAR's in-memory layout after
    /// `sjdbBuildIndex.cpp:293` (`memcpy(G+chrStart[nChrReal], Gsj, nGsj)`):
    /// the Gsj bytes live immediately after the forward real-genome bytes
    /// (`chr_start[n_chr_real]` stays pinned at the pre-sjdb forward total,
    /// matching STAR's `chrStart.txt` — verified against STAR docker output).
    ///
    /// `n_genome` grows to include Gsj. `chr_start` / `chr_length` /
    /// `chr_name` / `n_chr_real` are NOT updated — Gsj lives outside the
    /// chromosome accounting.
    pub fn append_sjdb(&mut self, gsj: &[u8]) {
        let old_n = self.n_genome;
        let new_n = old_n + gsj.len() as u64;

        let mut new_seq = vec![GENOME_SPACING_CHAR; (new_n * 2) as usize];
        new_seq[..old_n as usize].copy_from_slice(&self.sequence[..old_n as usize]);
        new_seq[old_n as usize..new_n as usize].copy_from_slice(gsj);

        // Rebuild RC over the extended forward range (STAR stores Gsj_RC
        // implicitly — rustar-aligner keeps the explicit `[fwd | RC]` layout).
        for i in 0..new_n as usize {
            let base = new_seq[i];
            let complement = if base < 4 { 3 - base } else { base };
            new_seq[2 * new_n as usize - 1 - i] = complement;
        }

        self.sequence = new_seq;
        self.n_genome = new_n;
    }

    /// Access a base from the genome (forward or reverse strand).
    ///
    /// # Arguments
    /// - `pos`: Position in the genome (0..2*n_genome)
    ///
    /// # Returns
    /// The base value (0-3 for ACGT, 4 for N, 5 for padding), or None if out of bounds.
    pub fn get_base(&self, pos: u64) -> Option<u8> {
        if pos < self.sequence.len() as u64 {
            Some(self.sequence[pos as usize])
        } else {
            None
        }
    }

    /// Get the chromosome containing a given genomic position.
    ///
    /// Uses binary search on chr_start for O(log n) lookup.
    ///
    /// # Returns
    /// `(chr_index, offset_within_chr)` or None if position is in padding.
    pub fn position_to_chr(&self, pos: u64) -> Option<(usize, u64)> {
        // Binary search: find the last chr_start that is <= pos
        let idx = self.chr_start[..self.n_chr_real].partition_point(|&start| start <= pos);
        if idx == 0 {
            return None;
        }
        let i = idx - 1;
        let start = self.chr_start[i];
        if pos < start + self.chr_length[i] {
            Some((i, pos - start))
        } else {
            None // Position is in padding between chromosomes
        }
    }

    /// Write genome index files to the specified directory.
    ///
    /// Creates:
    /// - `Genome` — raw binary file (n_genome bytes, forward strand only)
    /// - `chrName.txt` — chromosome names, one per line
    /// - `chrLength.txt` — chromosome lengths, one per line
    /// - `chrStart.txt` — chromosome start positions + final n_genome entry
    /// - `chrNameLength.txt` — tab-separated name + length
    /// - `genomeParameters.txt` — key-value pairs of genome generation parameters
    pub fn write_index_files(&self, dir: &Path, params: &Parameters) -> Result<(), Error> {
        use std::fs;
        use std::io::Write;

        // Create output directory if needed
        fs::create_dir_all(dir).map_err(|e| Error::io(e, dir))?;

        // Write Genome file (forward strand only, n_genome bytes)
        let genome_path = dir.join("Genome");
        fs::write(&genome_path, &self.sequence[..self.n_genome as usize])
            .map_err(|e| Error::io(e, &genome_path))?;

        // Write chrName.txt
        let chr_name_path = dir.join("chrName.txt");
        let mut f = fs::File::create(&chr_name_path).map_err(|e| Error::io(e, &chr_name_path))?;
        for name in &self.chr_name {
            writeln!(f, "{}", name).map_err(|e| Error::io(e, &chr_name_path))?;
        }

        // Write chrLength.txt
        let chr_length_path = dir.join("chrLength.txt");
        let mut f =
            fs::File::create(&chr_length_path).map_err(|e| Error::io(e, &chr_length_path))?;
        for &len in &self.chr_length {
            writeln!(f, "{}", len).map_err(|e| Error::io(e, &chr_length_path))?;
        }

        // Write chrStart.txt (includes the extra n_genome entry)
        let chr_start_path = dir.join("chrStart.txt");
        let mut f = fs::File::create(&chr_start_path).map_err(|e| Error::io(e, &chr_start_path))?;
        for &start in &self.chr_start {
            writeln!(f, "{}", start).map_err(|e| Error::io(e, &chr_start_path))?;
        }

        // Write chrNameLength.txt (tab-separated)
        let chr_name_length_path = dir.join("chrNameLength.txt");
        let mut f = fs::File::create(&chr_name_length_path)
            .map_err(|e| Error::io(e, &chr_name_length_path))?;
        for (name, &len) in self.chr_name.iter().zip(&self.chr_length) {
            writeln!(f, "{}\t{}", name, len).map_err(|e| Error::io(e, &chr_name_length_path))?;
        }

        // Write genomeParameters.txt — byte-for-byte matching STAR's
        // `genomeParametersWrite.cpp` layout (order, tab/space separators,
        // trailing whitespace on vector values). STAR's loader reads these
        // keys via `<<` streaming; the leading `###` comment lines are
        // skipped.
        self.write_genome_parameters_txt(dir, params)?;

        Ok(())
    }

    fn write_genome_parameters_txt(&self, dir: &Path, params: &Parameters) -> Result<(), Error> {
        use std::fs;
        use std::io::Write;

        let path = dir.join("genomeParameters.txt");
        let mut f = fs::File::create(&path).map_err(|e| Error::io(e, &path))?;

        // STAR writes: `### <commandLineFull>\n` where commandLineFull is
        // "<argv[0]>   --<name1> <val1>   --<name2> <val2> ...".  We emit
        // the same skeleton using our known-at-invocation parameters.
        // Not exposed for retrospective exact-byte match against an arbitrary
        // STAR run's commandLineFull — see `docs/genome_params_divergence.md`
        // for the short list of parameters we echo.
        let fasta_list = params
            .genome_fasta_files
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(" ");
        let gtf = params
            .sjdb_gtf_file
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "-".to_string());
        writeln!(
            f,
            "### STAR   --runMode genomeGenerate      --runThreadN {thr}   --genomeDir {dir}   --genomeFastaFiles {fa}      --genomeSAindexNbases {sai}   --sjdbGTFfile {gtf}   --sjdbOverhang {ov}",
            thr = params.run_thread_n,
            dir = dir.display(),
            fa = fasta_list,
            sai = params.genome_sa_index_nbases,
            gtf = gtf,
            ov = params.sjdb_overhang,
        )
        .map_err(|e| Error::io(e, &path))?;

        // GstrandBit: floor(log2(nGenome + limitSjdbInsertNsj*sjdbLength))+1,
        // clamped at a minimum of 32. STAR's default limitSjdbInsertNsj is
        // 1e6, so the log-derived value is always >= 32 for sane genomes;
        // use the clamped minimum directly.
        writeln!(f, "### GstrandBit 32").map_err(|e| Error::io(e, &path))?;

        // Value lines — exact order + tab/space separators match STAR's
        // genomeParametersWrite.cpp:11-42.
        writeln!(f, "versionGenome\t2.7.4a").map_err(|e| Error::io(e, &path))?;
        writeln!(f, "genomeType\tFull").map_err(|e| Error::io(e, &path))?;

        // Vectors: `key\t<v1> <v2> ... \n` with trailing space after each
        // value (STAR emits one ostream `<<` per element).
        write!(f, "genomeFastaFiles\t").map_err(|e| Error::io(e, &path))?;
        for p in &params.genome_fasta_files {
            write!(f, "{} ", p.display()).map_err(|e| Error::io(e, &path))?;
        }
        writeln!(f).map_err(|e| Error::io(e, &path))?;

        writeln!(f, "genomeSAindexNbases\t{}", params.genome_sa_index_nbases)
            .map_err(|e| Error::io(e, &path))?;
        writeln!(f, "genomeChrBinNbits\t{}", params.genome_chr_bin_nbits)
            .map_err(|e| Error::io(e, &path))?;
        writeln!(f, "genomeSAsparseD\t{}", params.genome_sa_sparse_d)
            .map_err(|e| Error::io(e, &path))?;

        writeln!(f, "genomeTransformType\tNone").map_err(|e| Error::io(e, &path))?;
        writeln!(f, "genomeTransformVCF\t-").map_err(|e| Error::io(e, &path))?;

        writeln!(f, "sjdbOverhang\t{}", params.sjdb_overhang).map_err(|e| Error::io(e, &path))?;

        // sjdbFileChrStartEnd: empty vector → `-` plus STAR's trailing space.
        writeln!(f, "sjdbFileChrStartEnd\t- ").map_err(|e| Error::io(e, &path))?;

        let gtf_str = params
            .sjdb_gtf_file
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "-".to_string());
        writeln!(f, "sjdbGTFfile\t{}", gtf_str).map_err(|e| Error::io(e, &path))?;
        writeln!(f, "sjdbGTFchrPrefix\t-").map_err(|e| Error::io(e, &path))?;
        writeln!(f, "sjdbGTFfeatureExon\texon").map_err(|e| Error::io(e, &path))?;
        writeln!(f, "sjdbGTFtagExonParentTranscript\ttranscript_id")
            .map_err(|e| Error::io(e, &path))?;
        writeln!(f, "sjdbGTFtagExonParentGene\tgene_id").map_err(|e| Error::io(e, &path))?;

        writeln!(f, "sjdbInsertSave\tBasic").map_err(|e| Error::io(e, &path))?;

        // genomeFileSizes: tab before the FIRST value, spaces between
        // subsequent values (STAR's pattern of `<< "\t"` then `<< " "`).
        // SA is 0 at this point because `Genome::write_index_files` runs
        // before SA is known; `GenomeIndex::write` patches it after.
        writeln!(f, "genomeFileSizes\t{} 0", self.n_genome).map_err(|e| Error::io(e, &path))?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_params(fasta_paths: Vec<std::path::PathBuf>, bin_nbits: u32) -> Parameters {
        use clap::Parser;
        let mut args = vec!["rustar-aligner", "--runMode", "genomeGenerate"];

        for path in &fasta_paths {
            args.push("--genomeFastaFiles");
            args.push(path.to_str().unwrap());
        }

        let bin_str = bin_nbits.to_string();
        args.extend(["--genomeChrBinNbits", &bin_str]);

        Parameters::parse_from(args)
    }

    #[test]
    fn single_chromosome_padding() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap(); // 4 bases

        let params = make_params(vec![file.path().to_path_buf()], 3); // bin_size = 8
        let genome = Genome::from_fasta(&params).unwrap();

        // Padding formula: n=4, then ((4+1)/8 + 1)*8 = (0+1)*8 = 8
        assert_eq!(genome.n_genome, 8);
        assert_eq!(genome.n_chr_real, 1);
        assert_eq!(genome.chr_start, vec![0, 8]);
        assert_eq!(genome.chr_length, vec![4]);

        // Check bases
        assert_eq!(genome.get_base(0), Some(0)); // A
        assert_eq!(genome.get_base(1), Some(1)); // C
        assert_eq!(genome.get_base(2), Some(2)); // G
        assert_eq!(genome.get_base(3), Some(3)); // T
        assert_eq!(genome.get_base(4), Some(5)); // padding
        assert_eq!(genome.get_base(7), Some(5)); // padding

        // Check reverse complement
        // Formula: pos(2*n - 1 - i) = complement(pos i)
        // For n=8: pos 15 = complement(pos 0), pos 12 = complement(pos 3), etc.
        assert_eq!(genome.get_base(15), Some(3)); // T (complement of A at pos 0)
        assert_eq!(genome.get_base(14), Some(2)); // G (complement of C at pos 1)
        assert_eq!(genome.get_base(13), Some(1)); // C (complement of G at pos 2)
        assert_eq!(genome.get_base(12), Some(0)); // A (complement of T at pos 3)
        assert_eq!(genome.get_base(11), Some(5)); // padding (complement of padding at pos 4)
        assert_eq!(genome.get_base(8), Some(5)); // padding (complement of padding at pos 7)
    }

    #[test]
    fn two_chromosomes_padding() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "AA").unwrap(); // 2 bases
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TT").unwrap(); // 2 bases

        let params = make_params(vec![file.path().to_path_buf()], 2); // bin_size = 4
        let genome = Genome::from_fasta(&params).unwrap();

        // chr1 starts at 0, length 2
        // After chr1: n=2, padding ((2+1)/4 + 1)*4 = (0+1)*4 = 4
        // chr2 starts at 4, length 2
        // After chr2: n=6, padding ((6+1)/4 + 1)*4 = (1+1)*4 = 8
        assert_eq!(genome.n_genome, 8);
        assert_eq!(genome.chr_start, vec![0, 4, 8]);
        assert_eq!(genome.chr_length, vec![2, 2]);

        // chr1 bases
        assert_eq!(genome.get_base(0), Some(0)); // A
        assert_eq!(genome.get_base(1), Some(0)); // A
        assert_eq!(genome.get_base(2), Some(5)); // padding
        assert_eq!(genome.get_base(3), Some(5)); // padding

        // chr2 bases
        assert_eq!(genome.get_base(4), Some(3)); // T
        assert_eq!(genome.get_base(5), Some(3)); // T
        assert_eq!(genome.get_base(6), Some(5)); // padding
        assert_eq!(genome.get_base(7), Some(5)); // padding
    }

    #[test]
    fn reverse_complement_correctness() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "ACGTN").unwrap(); // A=0, C=1, G=2, T=3, N=4

        let params = make_params(vec![file.path().to_path_buf()], 3); // bin_size = 8
        let genome = Genome::from_fasta(&params).unwrap();

        let n = genome.n_genome as usize;

        // Forward: A C G T N (then padding)
        assert_eq!(genome.sequence[0], 0); // A
        assert_eq!(genome.sequence[1], 1); // C
        assert_eq!(genome.sequence[2], 2); // G
        assert_eq!(genome.sequence[3], 3); // T
        assert_eq!(genome.sequence[4], 4); // N

        // Reverse complement should be at positions [2n-1, 2n-2, 2n-3, 2n-4, 2n-5]
        // which maps to the reverse of [0,1,2,3,4]
        assert_eq!(genome.sequence[2 * n - 1], 3); // T (complement of A at pos 0)
        assert_eq!(genome.sequence[2 * n - 1 - 1], 2); // G (complement of C at pos 1)
        assert_eq!(genome.sequence[2 * n - 1 - 2], 1); // C (complement of G at pos 2)
        assert_eq!(genome.sequence[2 * n - 1 - 3], 0); // A (complement of T at pos 3)
        assert_eq!(genome.sequence[2 * n - 1 - 4], 4); // N (complement of N at pos 4)
    }

    #[test]
    fn position_to_chr_mapping() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "AAA").unwrap();
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TTT").unwrap();

        let params = make_params(vec![file.path().to_path_buf()], 2); // bin_size = 4
        let genome = Genome::from_fasta(&params).unwrap();

        // chr1: positions 0-2 (starts at 0, length 3)
        // After chr1: n=3, padding ((3+1)/4 + 1)*4 = (1+1)*4 = 8
        // chr2: positions 8-10 (starts at 8, length 3)
        // After chr2: n=11, padding ((11+1)/4 + 1)*4 = (3+1)*4 = 16

        assert_eq!(genome.n_genome, 16);
        assert_eq!(genome.chr_start, vec![0, 8, 16]);

        assert_eq!(genome.position_to_chr(0), Some((0, 0)));
        assert_eq!(genome.position_to_chr(1), Some((0, 1)));
        assert_eq!(genome.position_to_chr(2), Some((0, 2)));
        assert_eq!(genome.position_to_chr(3), None); // padding
        assert_eq!(genome.position_to_chr(7), None); // padding
        assert_eq!(genome.position_to_chr(8), Some((1, 0)));
        assert_eq!(genome.position_to_chr(9), Some((1, 1)));
        assert_eq!(genome.position_to_chr(10), Some((1, 2)));
        assert_eq!(genome.position_to_chr(11), None); // padding
    }

    #[test]
    fn append_sjdb_extends_forward_and_rebuilds_rc() {
        // Forward "ACGT" + spacer padding (bin 3 → 8 bytes padded).
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        let params = make_params(vec![file.path().to_path_buf()], 3);
        let mut genome = Genome::from_fasta(&params).unwrap();
        assert_eq!(genome.n_genome, 8);

        // Gsj bytes: donor "AC" + acceptor "GT" + spacer 5 (five total).
        let gsj: Vec<u8> = vec![0, 1, 2, 3, 5];
        let gsj_len = gsj.len();
        genome.append_sjdb(&gsj);

        // n_genome grows by gsj_len; chr_start pinned at pre-sjdb boundary.
        assert_eq!(genome.n_genome, 8 + gsj_len as u64);
        assert_eq!(genome.chr_start, vec![0, 8]);
        assert_eq!(genome.n_chr_real, 1);

        // Forward is [real 0..8 | gsj 8..13].
        assert_eq!(&genome.sequence[..4], &[0, 1, 2, 3]);
        assert_eq!(&genome.sequence[8..13], gsj.as_slice());

        // RC over the extended forward range. sequence[2n-1-i] = complement(sequence[i]).
        let new_n = genome.n_genome as usize;
        assert_eq!(genome.sequence[2 * new_n - 1 - 8], 3); // complement of A at fwd[8]=0
        assert_eq!(genome.sequence[2 * new_n - 1 - 12], 5); // spacer stays 5
        assert_eq!(genome.sequence.len(), 2 * new_n);
    }

    #[test]
    fn append_sjdb_with_empty_gsj_is_noop() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        let params = make_params(vec![file.path().to_path_buf()], 3);
        let mut genome = Genome::from_fasta(&params).unwrap();
        let before = genome.sequence.clone();
        let before_n = genome.n_genome;
        genome.append_sjdb(&[]);
        assert_eq!(genome.n_genome, before_n);
        assert_eq!(genome.sequence, before);
    }
}
