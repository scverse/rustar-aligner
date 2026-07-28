use std::num::NonZeroUsize;
use std::path::PathBuf;

use clap::{CommandFactory, Parser};

/// Parse a memory string into bytes. Accepts plain integers or a suffix:
/// K/k = ×1024, M/m = ×1024², G/g = ×1024³, T/t = ×1024⁴.
/// Examples: `31000000000`, `64G`, `512M`, `1T`.
fn parse_mem_bytes(s: &str) -> Result<u64, String> {
    let s = s.trim();
    let (digits, shift) = match s.chars().last() {
        Some('K' | 'k') => (&s[..s.len() - 1], 10),
        Some('M' | 'm') => (&s[..s.len() - 1], 20),
        Some('G' | 'g') => (&s[..s.len() - 1], 30),
        Some('T' | 't') => (&s[..s.len() - 1], 40),
        _ => (s, 0),
    };
    let base: u64 = digits.trim().parse().map_err(|_| {
        format!("invalid memory value '{s}' — expected a number with optional K/M/G/T suffix")
    })?;
    base.checked_shl(shift)
        .ok_or_else(|| format!("memory value '{s}' overflows u64"))
}

mod sam;

pub use sam::{OutSamFormat, OutSamSortOrder, OutSamType, OutSamUnmapped, SamAttributes};

// ---------------------------------------------------------------------------
// Run mode enum
// ---------------------------------------------------------------------------

/// STAR's `--runMode` values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunMode {
    AlignReads,
    GenomeGenerate,
    InputAlignmentsFromBAM,
    LiftOver,
}

impl std::str::FromStr for RunMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "alignReads" => Ok(Self::AlignReads),
            "genomeGenerate" => Ok(Self::GenomeGenerate),
            "inputAlignmentsFromBAM" => Ok(Self::InputAlignmentsFromBAM),
            "liftOver" => Ok(Self::LiftOver),
            _ => Err(format!(
                "unknown runMode '{s}'; expected 'alignReads', 'genomeGenerate', \
                 'inputAlignmentsFromBAM', or 'liftOver'"
            )),
        }
    }
}

impl std::fmt::Display for RunMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AlignReads => write!(f, "alignReads"),
            Self::GenomeGenerate => write!(f, "genomeGenerate"),
            Self::InputAlignmentsFromBAM => write!(f, "inputAlignmentsFromBAM"),
            Self::LiftOver => write!(f, "liftOver"),
        }
    }
}

// ---------------------------------------------------------------------------
// Junction motif filter enum
// ---------------------------------------------------------------------------

/// Filter mode for splice junction motifs (outFilterIntronMotifs)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntronMotifFilter {
    /// Accept all junction motifs (no filtering)
    None,
    /// Remove alignments with non-canonical junctions (STAR default for RNA-seq)
    RemoveNoncanonical,
    /// Remove non-canonical junctions only if unannotated
    RemoveNoncanonicalUnannotated,
}

impl std::str::FromStr for IntronMotifFilter {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "RemoveNoncanonical" => Ok(Self::RemoveNoncanonical),
            "RemoveNoncanonicalUnannotated" => Ok(Self::RemoveNoncanonicalUnannotated),
            _ => Err(format!(
                "unknown outFilterIntronMotifs '{s}'; expected 'None', 'RemoveNoncanonical', or 'RemoveNoncanonicalUnannotated'"
            )),
        }
    }
}

impl std::fmt::Display for IntronMotifFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, "None"),
            Self::RemoveNoncanonical => write!(f, "RemoveNoncanonical"),
            Self::RemoveNoncanonicalUnannotated => write!(f, "RemoveNoncanonicalUnannotated"),
        }
    }
}

// ---------------------------------------------------------------------------
// Read-end alignment type (alignEndsType)
// ---------------------------------------------------------------------------

/// Read-end extension policy (`--alignEndsType`). Mirrors STAR's
/// `alignEndsType.ext[iMate][iEnd]` boolean matrix, where `iEnd == 0` is the
/// read 5' end and `iEnd == 1` the 3' end. `true` forces full end-to-end
/// extension of that end (no terminal soft-clip); `false` is local alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignEndsType {
    /// `ext[iMate][iEnd]`: force full extension of mate `iMate`'s `iEnd` end.
    pub ext: [[bool; 2]; 2],
}

impl Default for AlignEndsType {
    /// `Local` — no forced extension on any end.
    fn default() -> Self {
        Self {
            ext: [[false; 2]; 2],
        }
    }
}

impl std::str::FromStr for AlignEndsType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // STAR Parameters.cpp: ext[iMate][iEnd], iEnd 0 = 5', 1 = 3'.
        let mut ext = [[false; 2]; 2];
        match s {
            "Local" => {}
            "EndToEnd" => ext = [[true, true], [true, true]],
            "Extend5pOfRead1" => ext[0][0] = true,
            "Extend5pOfReads12" => {
                ext[0][0] = true;
                ext[1][0] = true;
            }
            "Extend3pOfRead1" => ext[0][1] = true,
            other => {
                return Err(format!(
                    "unknown/unimplemented value for --alignEndsType: '{other}'; expected \
                     'Local', 'EndToEnd', 'Extend5pOfRead1', 'Extend5pOfReads12', or 'Extend3pOfRead1'"
                ));
            }
        }
        Ok(Self { ext })
    }
}

// ---------------------------------------------------------------------------
// Intron strand filter enum
// ---------------------------------------------------------------------------

/// Filter mode for intron strand consistency (outFilterIntronStrands)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntronStrandFilter {
    /// Accept all alignments regardless of strand consistency
    None,
    /// Remove alignments where junction motifs imply conflicting transcript strands
    RemoveInconsistentStrands,
}

impl std::str::FromStr for IntronStrandFilter {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "RemoveInconsistentStrands" => Ok(Self::RemoveInconsistentStrands),
            _ => Err(format!(
                "unknown outFilterIntronStrands '{s}'; expected 'None' or 'RemoveInconsistentStrands'"
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// Multimapper output order / primary selection
// ---------------------------------------------------------------------------

/// STAR's `--outMultimapperOrder` — the order in which multi-mapping alignments
/// are reported and how the primary is chosen among equal-scoring loci.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum MultimapperOrder {
    /// STAR default. Alignments are reported in a deterministic order and the
    /// primary is `trBest` (max score → smaller genomic length → earliest
    /// discovered). No RNG is used for primary selection.
    #[default]
    Old24,
    /// Multimappers are shuffled with a per-read RNG and a best-scoring
    /// alignment is marked primary. Deterministic per read (thread-count
    /// invariant), unlike STAR's per-thread RNG.
    Random,
}

impl std::str::FromStr for MultimapperOrder {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Old_2.4" => Ok(Self::Old24),
            "Random" => Ok(Self::Random),
            _ => Err(format!(
                "unknown outMultimapperOrder value: '{s}'; expected 'Old_2.4' or 'Random'"
            )),
        }
    }
}

impl std::fmt::Display for MultimapperOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Old24 => write!(f, "Old_2.4"),
            Self::Random => write!(f, "Random"),
        }
    }
}

// ---------------------------------------------------------------------------
// Standard output streaming
// ---------------------------------------------------------------------------

/// STAR's `--outStd` — route primary alignment output to stdout.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum OutStd {
    #[default]
    None,
    Sam,
    BamUnsorted,
    BamSortedByCoordinate,
}

impl std::str::FromStr for OutStd {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "SAM" => Ok(Self::Sam),
            "BAM_Unsorted" => Ok(Self::BamUnsorted),
            "BAM_SortedByCoordinate" => Ok(Self::BamSortedByCoordinate),
            _ => Err(format!(
                "unknown outStd value: '{s}'; expected None, SAM, BAM_Unsorted, or BAM_SortedByCoordinate"
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// Unmapped reads FASTQ output
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum OutReadsUnmapped {
    #[default]
    None,
    Fastx,
}

impl std::str::FromStr for OutReadsUnmapped {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "Fastx" => Ok(Self::Fastx),
            _ => Err(format!(
                "unknown outReadsUnmapped value: '{s}'; expected 'None' or 'Fastx'"
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// Output filter type
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum OutFilterType {
    #[default]
    Normal,
    BySJout,
}

impl std::str::FromStr for OutFilterType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Normal" => Ok(Self::Normal),
            "BySJout" => Ok(Self::BySJout),
            _ => Err(format!("unknown outFilterType value: '{s}'")),
        }
    }
}

// ---------------------------------------------------------------------------
// SAM primary-flag assignment
// ---------------------------------------------------------------------------

/// Which alignment(s) get the SAM primary flag (`--outSAMprimaryFlag`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutSamPrimaryFlag {
    /// Only the single best-tie-broken alignment is primary (STAR default).
    #[default]
    OneBestScore,
    /// Every alignment tied for the best score is marked primary.
    AllBestScore,
}

impl std::str::FromStr for OutSamPrimaryFlag {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "OneBestScore" => Ok(Self::OneBestScore),
            "AllBestScore" => Ok(Self::AllBestScore),
            _ => Err(format!(
                "unknown outSAMprimaryFlag value: '{s}'; expected 'OneBestScore' or 'AllBestScore'"
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// SAM read-ID naming
// ---------------------------------------------------------------------------

/// QNAME source for SAM/FASTX output (`--outSAMreadID`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutSamReadId {
    /// Use the FASTQ read name as-is (default).
    #[default]
    Standard,
    /// Replace the QNAME with the read's 1-based input index, uniformly across every record
    /// (mapped, unmapped, FASTX) the read emits.
    Number,
}

impl std::str::FromStr for OutSamReadId {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Standard" => Ok(Self::Standard),
            "Number" => Ok(Self::Number),
            _ => Err(format!(
                "unknown outSAMreadID value: '{s}'; expected 'Standard' or 'Number'"
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// Two-pass mode
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum TwopassMode {
    #[default]
    None,
    Basic,
}

impl std::str::FromStr for TwopassMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "Basic" => Ok(Self::Basic),
            _ => Err(format!("unknown twopassMode value: '{s}'")),
        }
    }
}

/// STAR's `--waspOutputMode`: WASP allele-specific-mapping filter output.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum WaspOutputMode {
    #[default]
    None,
    /// Emit the `vW` SAM tag (and `vA`/`vG` when requested in `--outSAMattributes`).
    SAMtag,
}

impl std::str::FromStr for WaspOutputMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "SAMtag" => Ok(Self::SAMtag),
            _ => Err(format!("unknown waspOutputMode value: '{s}'")),
        }
    }
}

// ---------------------------------------------------------------------------
// STARsolo (single-cell) type
// ---------------------------------------------------------------------------

/// STAR's `--soloType` — selects the single-cell barcode geometry.
///
/// Mirrors STAR's `ParametersSolo::typeStr` values. Only `None` and
/// `CB_UMI_Simple` (droplet 10x-style) are functional in Phase 14.1; the
/// remaining variants are parsed so the CLI accepts them and later sub-phases
/// can fill in behavior.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum SoloType {
    /// Not a single-cell run (default).
    #[default]
    None,
    /// One cell barcode + one UMI at fixed positions in the barcode read
    /// (10x Chromium, Drop-seq, inDrops-simple, etc.). STAR alias: `Droplet`.
    CbUmiSimple,
    /// Multi-segment cell barcode and/or UMI, optionally adapter-anchored.
    CbUmiComplex,
    /// Barcodes passed through as SAM tags only (no collapsing).
    CbSamTagOut,
    /// Plate-based Smart-seq: one cell per read-group, no UMI.
    SmartSeq,
}

impl std::str::FromStr for SoloType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            // STAR accepts both the descriptive name and the `Droplet` alias.
            "CB_UMI_Simple" | "Droplet" => Ok(Self::CbUmiSimple),
            "CB_UMI_Complex" => Ok(Self::CbUmiComplex),
            "CB_samTagOut" => Ok(Self::CbSamTagOut),
            "SmartSeq" => Ok(Self::SmartSeq),
            _ => Err(format!(
                "unknown soloType '{s}'; expected None, CB_UMI_Simple, CB_UMI_Complex, CB_samTagOut, or SmartSeq"
            )),
        }
    }
}

impl std::fmt::Display for SoloType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::None => "None",
            Self::CbUmiSimple => "CB_UMI_Simple",
            Self::CbUmiComplex => "CB_UMI_Complex",
            Self::CbSamTagOut => "CB_samTagOut",
            Self::SmartSeq => "SmartSeq",
        };
        write!(f, "{s}")
    }
}

// ---------------------------------------------------------------------------
// Parameters struct
// ---------------------------------------------------------------------------

/// rustar-aligner command-line parameters, matching STAR's `--camelCase` argument names.
///
/// Only the ~40 most important parameters are included; more will be added
/// incrementally as later phases need them.
#[derive(Debug, Clone, Parser)]
#[command(
    name = "rustar-aligner",
    about = "RNA-seq aligner (Rust reimplementation of STAR)",
    version,
    long_version = concat!(env!("CARGO_PKG_VERSION"), "\n", env!("VERSION_BODY")),
)]
pub struct Parameters {
    // ── Run ─────────────────────────────────────────────────────────────
    /// Run mode: alignReads or genomeGenerate
    #[arg(long = "runMode", default_value = "alignReads")]
    pub run_mode: RunMode,

    /// Number of threads
    #[arg(long = "runThreadN", default_value_t = NonZeroUsize::new(1).unwrap())]
    pub run_thread_n: NonZeroUsize,

    /// Random number generator seed for tie-breaking among equal-scoring alignments
    #[arg(long = "runRNGseed", default_value_t = 777)]
    pub run_rng_seed: u64,

    /// Temporary directory for large external-memory index-construction files
    #[arg(long = "tempDir")]
    pub temp_dir: Option<PathBuf>,

    // ── Genome ──────────────────────────────────────────────────────────
    /// Path to genome index directory
    #[arg(long = "genomeDir", default_value = "./GenomeDir")]
    pub genome_dir: PathBuf,

    /// FASTA file(s) with genome reference sequences (for genomeGenerate)
    #[arg(long = "genomeFastaFiles", num_args = 1..)]
    pub genome_fasta_files: Vec<PathBuf>,

    /// UCSC chain file(s) for `--runMode liftOver`. Only the first is ever
    /// used (matches STAR's `Chain` dispatch, whose loop over chain files
    /// unconditionally exits after the first iteration).
    #[arg(long = "genomeChainFiles", num_args = 1..)]
    pub genome_chain_files: Vec<PathBuf>,

    /// Length of SA pre-indexing string (log2-based)
    #[arg(long = "genomeSAindexNbases", default_value_t = 14)]
    pub genome_sa_index_nbases: u32,

    /// Log2(chromosome bin size) for genome storage
    #[arg(long = "genomeChrBinNbits", default_value_t = 18)]
    pub genome_chr_bin_nbits: u32,

    /// Suffix array sparsity (larger = less RAM, slower mapping)
    #[arg(long = "genomeSAsparseD", default_value_t = 1)]
    pub genome_sa_sparse_d: u32,

    /// Substitute VCF alleles into the genome at genomeGenerate (`None`,
    /// `Haploid`, or `Diploid`). Requires `--genomeTransformVCF`; incompatible
    /// with `--sjdbGTFfile`. `Diploid` is genotype-aware and duplicates the
    /// genome into `_h1`/`_h2` haplotype chromosomes (see `Parameters::validate`)
    #[arg(long = "genomeTransformType", default_value = "None")]
    pub genome_transform_type: String,

    /// VCF of variants for `--genomeTransformType`
    #[arg(long = "genomeTransformVCF")]
    pub genome_transform_vcf: Option<PathBuf>,

    // ── Read files ──────────────────────────────────────────────────────
    /// Input read file(s); second file is mate 2 for paired-end
    #[arg(long = "readFilesIn", num_args = 1..=2)]
    pub read_files_in: Vec<PathBuf>,

    /// Command to decompress input files (e.g. "zcat" for .gz)
    #[arg(long = "readFilesCommand")]
    pub read_files_command: Option<String>,

    /// `--soloType SmartSeq` manifest: a TSV with `read1 <TAB> read2 <TAB> cellID`
    /// per line (`read2` = `-` for single-end). Each line is one plate-well cell;
    /// reads are counted per gene with no UMI.
    #[arg(long = "readFilesManifest")]
    pub read_files_manifest: Option<PathBuf>,

    /// Number of reads to map; -1 = all
    #[arg(long = "readMapNumber", default_value_t = -1, allow_hyphen_values = true)]
    pub read_map_number: i64,

    /// Bases to clip from the 5' end of each mate. One value applies to both
    /// mates; two values are per-mate (mate 1, then mate 2).
    #[arg(long = "clip5pNbases", num_args = 1..=2, default_values_t = vec![0u32])]
    pub clip5p_nbases: Vec<u32>,

    /// Bases to clip from the 3' end of each mate. One value applies to both
    /// mates; two values are per-mate (mate 1, then mate 2).
    #[arg(long = "clip3pNbases", num_args = 1..=2, default_values_t = vec![0u32])]
    pub clip3p_nbases: Vec<u32>,

    /// Adapter clipping type applied to the cDNA read: `Hamming` (default,
    /// adapter-sequence based, no-op when no adapter is configured) or
    /// `CellRanger4` (clip the 10x TSO from the 5' end and trim the 3' polyA
    /// tail, to match CellRanger ≥ 4.0).
    #[arg(long = "clipAdapterType", default_value = "Hamming")]
    pub clip_adapter_type: String,

    /// 3' adapter sequence to clip (Hamming scan), `-` = none
    #[arg(long = "clip3pAdapterSeq", default_value = "-")]
    pub clip3p_adapter_seq: String,

    /// Max mismatch proportion for the 3' adapter clip
    #[arg(long = "clip3pAdapterMMp", default_value_t = 0.1)]
    pub clip3p_adapter_mmp: f64,

    /// Extra bases to clip from the 3' end after the adapter
    #[arg(long = "clip3pAfterAdapterNbases", default_value_t = 0)]
    pub clip3p_after_adapter_nbases: u32,

    /// Extra bases to clip from the 5' end after `clip5pNbases`
    #[arg(long = "clip5pAfterAdapterNbases", default_value_t = 0)]
    pub clip5p_after_adapter_nbases: u32,

    /// Coverage signal output: `None` (default) or `bedGraph`. `wiggle` and a 2nd
    /// word (`read1_5p`/`read2`) are rejected (see `Parameters::validate`)
    #[arg(long = "outWigType", num_args = 1.., default_values_t = vec![String::from("None")])]
    pub out_wig_type: Vec<String>,

    /// Coverage-signal normalization: `RPM` (default) or `None` (raw counts)
    #[arg(long = "outWigNorm", default_value = "RPM")]
    pub out_wig_norm: String,

    /// Coverage-signal strandedness: `Stranded` (default). `Unstranded` is
    /// rejected (see `Parameters::validate`)
    #[arg(long = "outWigStrand", default_value = "Stranded")]
    pub out_wig_strand: String,

    // ── Output ──────────────────────────────────────────────────────────
    /// Output file name prefix (including path)
    #[arg(long = "outFileNamePrefix", default_value = "./")]
    pub out_file_name_prefix: String,

    /// (`--runMode inputAlignmentsFromBAM`) the input BAM to re-process
    #[arg(long = "inputBAMfile", default_value = "-")]
    pub input_bam_file: String,

    /// (`--runMode inputAlignmentsFromBAM`) mark PCR duplicates in the input BAM: `-` = off,
    /// `UniqueIdentical` (marks multimappers too) or `UniqueIdenticalNotMulti` (leaves them unmarked)
    #[arg(long = "bamRemoveDuplicatesType", default_value = "-")]
    pub bam_remove_duplicates_type: String,

    /// Compare this many mate-2 SEQ bases when deduplicating (RAMPAGE); 0 = don't compare SEQ
    #[arg(
        long = "bamRemoveDuplicatesMate2basesN",
        default_value_t = 0,
        allow_hyphen_values = true
    )]
    pub bam_remove_duplicates_mate2_bases_n: i64,

    #[command(flatten)]
    pub out_sam_type: OutSamType,

    /// BAM compression level: -1 = uncompressed, 1 (default) to 9 (maximum)
    #[arg(
        long = "outBAMcompression",
        default_value_t = 1,
        allow_hyphen_values = true
    )]
    pub out_bam_compression: i32,

    /// Maximum RAM for coordinate-sorted BAM sorting. Accepts bytes or a suffix: 8G, 512M, 1T. 0 = unlimited.
    #[arg(long = "limitBAMsortRAM", default_value = "0", value_parser = parse_mem_bytes)]
    pub limit_bam_sort_ram: u64,

    /// Maximum RAM for genome generation. Accepts bytes or a suffix: 64G, 512M, 1T.
    #[arg(long = "limitGenomeGenerateRAM", default_value = "31G", value_parser = parse_mem_bytes)]
    pub limit_genome_generate_ram: u64,

    /// Route primary alignment output to stdout instead of a file.
    /// Values: None (default), SAM, BAM_Unsorted, BAM_SortedByCoordinate.
    #[arg(long = "outStd", default_value = "None")]
    pub out_std: OutStd,

    /// Strand field: None or intronMotif
    #[arg(long = "outSAMstrandField", default_value = "None")]
    pub out_sam_strand_field: String,

    /// SAM output order: `Paired` (default) or `PairedKeepInputOrder`. rustar-aligner
    /// always emits records in input (FASTQ) order regardless of thread count, so both
    /// values behave identically; the flag is accepted for STAR/pipeline compatibility.
    #[arg(long = "outSAMorder", default_value = "Paired")]
    pub out_sam_order: String,

    /// SAM optional-tag set (`Standard`, `All`, `None`, or any combination
    /// of NH HI AS NM nM MD jM jI XS RG). Parsed into a bitflags struct.
    #[command(flatten)]
    pub out_sam_attributes: SamAttributes,

    /// Read group line(s) for the `@RG` SAM header. Space-separated fields;
    /// a bare `,` separates multiple RG blocks. Each block must start with `ID:`.
    /// Default `-` means no `@RG` line (matches STAR).
    #[arg(long = "outSAMattrRGline", num_args = 1.., default_values_t = vec!["-".to_string()])]
    pub out_sam_attr_rg_line: Vec<String>,

    #[command(flatten)]
    pub out_sam_unmapped: OutSamUnmapped,

    /// Output unmapped reads to FASTQ file(s): None or Fastx
    #[arg(long = "outReadsUnmapped", default_value = "None")]
    pub out_reads_unmapped: OutReadsUnmapped,

    /// MAPQ value for unique mappers
    #[arg(long = "outSAMmapqUnique", default_value_t = 255)]
    pub out_sam_mapq_unique: u8,

    /// Max number of multiple alignments per read in SAM output (-1 = all)
    #[arg(long = "outSAMmultNmax", default_value_t = -1, allow_hyphen_values = true)]
    pub out_sam_mult_nmax: i32,

    /// Bits OR-ed into every mapped record's FLAG (`--outSAMflagOR`, default 0)
    #[arg(long = "outSAMflagOR", default_value_t = 0)]
    pub out_sam_flag_or: u32,

    /// Bits AND-ed into every mapped record's FLAG (`--outSAMflagAND`, default 65535 = keep all)
    #[arg(long = "outSAMflagAND", default_value_t = 65535)]
    pub out_sam_flag_and: u32,

    /// Which alignment(s) get the SAM primary flag: `OneBestScore` (default) or `AllBestScore`
    #[arg(long = "outSAMprimaryFlag", default_value = "OneBestScore")]
    pub out_sam_primary_flag: OutSamPrimaryFlag,

    /// Start value of the `HI` SAM attribute (STAR default 1; CellRanger convention uses 0)
    #[arg(long = "outSAMattrIHstart", default_value_t = 1)]
    pub out_sam_attr_ih_start: u32,

    /// QNAME source for SAM/FASTX output: `Standard` (default, the FASTQ read name) or
    /// `Number` (the read's 1-based input index)
    #[arg(long = "outSAMreadID", default_value = "Standard")]
    pub out_sam_read_id: OutSamReadId,

    /// TLEN calculation: 1 (default, whole combined-transcript span) or 2 (per-mate span,
    /// signed by whichever mate is genomically leftmost)
    #[arg(long = "outSAMtlen", default_value_t = 1)]
    pub out_sam_tlen: u8,

    /// Output filter type: Normal or BySJout
    #[arg(long = "outFilterType", default_value = "Normal")]
    pub out_filter_type: OutFilterType,

    /// Max multimap loci (reads mapping to more loci are unmapped)
    #[arg(long = "outFilterMultimapNmax", default_value_t = 10)]
    pub out_filter_multimap_nmax: u32,

    /// Score range for multi-mapping (keep alignments within this range of best score)
    #[arg(long = "outFilterMultimapScoreRange", default_value_t = 1)]
    pub out_filter_multimap_score_range: i32,

    /// Max mismatches per pair
    #[arg(long = "outFilterMismatchNmax", default_value_t = 10)]
    pub out_filter_mismatch_nmax: u32,

    /// Max ratio of mismatches to mapped length
    #[arg(long = "outFilterMismatchNoverLmax", default_value_t = 0.3)]
    pub out_filter_mismatch_nover_lmax: f64,

    /// Min alignment score (absolute)
    #[arg(long = "outFilterScoreMin", default_value_t = 0)]
    pub out_filter_score_min: i32,

    /// Min alignment score normalized to read length
    #[arg(long = "outFilterScoreMinOverLread", default_value_t = 0.66)]
    pub out_filter_score_min_over_lread: f64,

    /// Min matched bases (absolute)
    #[arg(long = "outFilterMatchNmin", default_value_t = 0)]
    pub out_filter_match_nmin: u32,

    /// Min matched bases normalized to read length
    #[arg(long = "outFilterMatchNminOverLread", default_value_t = 0.66)]
    pub out_filter_match_nmin_over_lread: f64,

    /// Order of multi-mapping alignments / primary selection
    #[arg(long = "outMultimapperOrder", default_value = "Old_2.4")]
    pub out_multimapper_order: MultimapperOrder,

    /// Filter alignments based on junction motifs
    #[arg(long = "outFilterIntronMotifs", default_value = "None")]
    pub out_filter_intron_motifs: IntronMotifFilter,

    /// Filter alignments with inconsistent intron strand motifs
    #[arg(
        long = "outFilterIntronStrands",
        default_value = "RemoveInconsistentStrands"
    )]
    pub out_filter_intron_strands: IntronStrandFilter,

    /// SJ filter: min overhang per motif category [noncanon, GT/AG, GC/AG, AT/AC]
    #[arg(long = "outSJfilterOverhangMin", num_args = 4,
          default_values_t = vec![30, 12, 12, 12])]
    pub out_sj_filter_overhang_min: Vec<i32>,

    /// SJ filter: min unique-mapping reads per motif category [noncanon, GT/AG, GC/AG, AT/AC]
    #[arg(long = "outSJfilterCountUniqueMin", num_args = 4,
          default_values_t = vec![3, 1, 1, 1])]
    pub out_sj_filter_count_unique_min: Vec<i32>,

    /// SJ filter: min total (unique+multi) reads per motif category [noncanon, GT/AG, GC/AG, AT/AC]
    #[arg(long = "outSJfilterCountTotalMin", num_args = 4,
          default_values_t = vec![3, 1, 1, 1])]
    pub out_sj_filter_count_total_min: Vec<i32>,

    /// SJ filter: min distance to other SJs per motif category [noncanon, GT/AG, GC/AG, AT/AC]
    #[arg(long = "outSJfilterDistToOtherSJmin", num_args = 4,
          default_values_t = vec![10, 0, 5, 10])]
    pub out_sj_filter_dist_to_other_sjmin: Vec<i32>,

    /// SJ filter: max intron length vs supporting read count
    /// [1_read, 2_reads, 3+_reads] — junctions with intron > threshold for their read count are filtered
    #[arg(long = "outSJfilterIntronMaxVsReadN", num_args = 3,
          default_values_t = vec![50000, 100000, 200000])]
    pub out_sj_filter_intron_max_vs_read_n: Vec<i64>,

    // ── Alignment scoring ───────────────────────────────────────────────
    /// Min intron size (smaller gaps are deletions)
    #[arg(long = "alignIntronMin", default_value_t = 21)]
    pub align_intron_min: u32,

    /// Max intron size; 0 = auto
    #[arg(long = "alignIntronMax", default_value_t = 0)]
    pub align_intron_max: u32,

    /// Max genomic distance between mates; 0 = auto
    #[arg(long = "alignMatesGapMax", default_value_t = 0)]
    pub align_mates_gap_max: u32,

    /// Read-end alignment mode: `Local` (default, soft-clip allowed),
    /// `EndToEnd` (force full extension, no soft-clip), `Extend5pOfRead1`,
    /// `Extend5pOfReads12`, or `Extend3pOfRead1`.
    #[arg(long = "alignEndsType", default_value = "Local")]
    pub align_ends_type: String,

    /// Min overlap (bases) between mates required to trigger merge-and-realign; 0 = off
    #[arg(long = "peOverlapNbasesMin", default_value_t = 0)]
    pub pe_overlap_nbases_min: u64,

    /// Max proportion of mismatches allowed in the mate-overlap region
    #[arg(long = "peOverlapMMp", default_value_t = 0.01)]
    pub pe_overlap_mmp: f64,

    /// Min mapped length of spliced mates (absolute, default 0 = off)
    #[arg(long = "alignSplicedMateMapLmin", default_value_t = 0)]
    pub align_spliced_mate_map_lmin: u32,

    /// Min mapped length of spliced mates as fraction of read length (default 0.66)
    #[arg(long = "alignSplicedMateMapLminOverLmate", default_value_t = 0.66)]
    pub align_spliced_mate_map_lmin_over_lmate: f64,

    /// Min overhang for novel spliced alignments
    #[arg(long = "alignSJoverhangMin", default_value_t = 5)]
    pub align_sj_overhang_min: u32,

    /// Min overhang for annotated splice junctions
    #[arg(long = "alignSJDBoverhangMin", default_value_t = 3)]
    pub align_sjdb_overhang_min: u32,

    /// Max mismatches for stitching SJs (4 ints: noncanonical, GC/AG, AT/AC, noncanonical)
    #[arg(long = "alignSJstitchMismatchNmax", num_args = 4,
          default_values_t = vec![0, -1, 0, 0], allow_hyphen_values = true)]
    pub align_sj_stitch_mismatch_nmax: Vec<i32>,

    /// Splice junction penalty (canonical)
    #[arg(long = "scoreGap", default_value_t = 0)]
    pub score_gap: i32,

    /// Non-canonical junction penalty
    #[arg(long = "scoreGapNoncan", default_value_t = -8, allow_hyphen_values = true)]
    pub score_gap_noncan: i32,

    /// GC/AG junction penalty
    #[arg(long = "scoreGapGCAG", default_value_t = -4, allow_hyphen_values = true)]
    pub score_gap_gcag: i32,

    /// AT/AC junction penalty
    #[arg(long = "scoreGapATAC", default_value_t = -8, allow_hyphen_values = true)]
    pub score_gap_atac: i32,

    /// Deletion open penalty
    #[arg(long = "scoreDelOpen", default_value_t = -2, allow_hyphen_values = true)]
    pub score_del_open: i32,

    /// Deletion extension penalty per base
    #[arg(long = "scoreDelBase", default_value_t = -2, allow_hyphen_values = true)]
    pub score_del_base: i32,

    /// Insertion open penalty
    #[arg(long = "scoreInsOpen", default_value_t = -2, allow_hyphen_values = true)]
    pub score_ins_open: i32,

    /// Insertion extension penalty per base
    #[arg(long = "scoreInsBase", default_value_t = -2, allow_hyphen_values = true)]
    pub score_ins_base: i32,

    /// Max score reduction for SJ stitching shift
    #[arg(long = "scoreStitchSJshift", default_value_t = 1)]
    pub score_stitch_sj_shift: i32,

    /// Extra score log-scaled with genomic length: scoreGenomicLengthLog2scale*log2(genomicLength)
    #[arg(long = "scoreGenomicLengthLog2scale", default_value_t = -0.25, allow_hyphen_values = true)]
    pub score_genomic_length_log2_scale: f64,

    // ── Seed and anchor parameters ──────────────────────────────────────
    /// Min read coverage for a window (relative to read length)
    #[arg(long = "winReadCoverageRelativeMin", default_value_t = 0.5)]
    pub win_read_coverage_relative_min: f64,

    /// Log2 of window bin size for seed clustering
    #[arg(long = "winBinNbits", default_value_t = 16)]
    pub win_bin_nbits: u32,

    /// Max number of bins for seed anchor distance
    #[arg(long = "winAnchorDistNbins", default_value_t = 9)]
    pub win_anchor_dist_nbins: u32,

    /// Number of bins to extend each alignment window by on each side
    #[arg(long = "winFlankNbins", default_value_t = 4)]
    pub win_flank_nbins: u32,

    /// Max number of loci a seed can map to (seeds with more loci are discarded)
    #[arg(long = "seedMultimapNmax", default_value_t = 10000)]
    pub seed_multimap_nmax: usize,

    /// Max number of seeds per read
    #[arg(long = "seedPerReadNmax", default_value_t = 1000)]
    pub seed_per_read_nmax: usize,

    /// Max number of seeds per window
    #[arg(long = "seedPerWindowNmax", default_value_t = 50)]
    pub seed_per_window_nmax: usize,

    /// Max distance between seed search start positions (defines Nstart = readLen/seedSearchStartLmax + 1)
    #[arg(long = "seedSearchStartLmax", default_value_t = 50)]
    pub seed_search_start_lmax: usize,

    /// seedSearchStartLmax normalized by read length (effective = min(seedSearchStartLmax, this * (readLen-1)))
    #[arg(long = "seedSearchStartLmaxOverLread", default_value_t = 1.0)]
    pub seed_search_start_lmax_over_lread: f64,

    /// Max seed length; 0 = unlimited (default)
    #[arg(long = "seedSearchLmax", default_value_t = 0)]
    pub seed_search_lmax: usize,

    /// Min mappable length for seed search while-loop termination (STAR default: 5)
    #[arg(long = "seedMapMin", default_value_t = 5)]
    pub seed_map_min: usize,

    /// Max number of loci anchors are allowed to map to
    #[arg(long = "winAnchorMultimapNmax", default_value_t = 50)]
    pub win_anchor_multimap_nmax: usize,

    /// Max number of alignment windows per read
    #[arg(long = "alignWindowsPerReadNmax", default_value_t = 10000)]
    pub align_windows_per_read_nmax: usize,

    /// Max number of transcripts per window
    #[arg(long = "alignTranscriptsPerWindowNmax", default_value_t = 100)]
    pub align_transcripts_per_window_nmax: usize,

    // ── Splice junction database ────────────────────────────────────────
    /// GTF file for splice junction annotations
    #[arg(long = "sjdbGTFfile")]
    pub sjdb_gtf_file: Option<PathBuf>,

    /// TSV file(s) of `chr start end [strand]` junctions (1-based intron first/last base) to
    /// insert into the sjdb, unioned with any `--sjdbGTFfile` junctions
    #[arg(long = "sjdbFileChrStartEnd", num_args = 1..)]
    pub sjdb_file_chr_start_end: Vec<PathBuf>,

    /// Prefix to add to chromosome names from GTF file (e.g. "chr" when GTF uses bare numbers)
    #[arg(long = "sjdbGTFchrPrefix", default_value = "")]
    pub sjdb_gtf_chr_prefix: String,

    /// Feature type in GTF file to be used as exons for transcript annotation
    #[arg(long = "sjdbGTFfeatureExon", default_value = "exon")]
    pub sjdb_gtf_feature_exon: String,

    /// GTF attribute name for parent transcript ID of exon features
    #[arg(
        long = "sjdbGTFtagExonParentTranscript",
        default_value = "transcript_id"
    )]
    pub sjdb_gtf_tag_exon_parent_transcript: String,

    /// GTF attribute name for parent gene ID of exon features
    #[arg(long = "sjdbGTFtagExonParentGene", default_value = "gene_id")]
    pub sjdb_gtf_tag_exon_parent_gene: String,

    /// GTF attribute name(s) for parent gene name of exon features; when several are given and
    /// several match, the last one in the list wins
    #[arg(long = "sjdbGTFtagExonParentGeneName", default_values_t = vec!["gene_name".to_string()], num_args = 1..)]
    pub sjdb_gtf_tag_exon_parent_gene_name: Vec<String>,

    /// GTF attribute name(s) for parent gene type of exon features; when several are given and
    /// several match, the last one in the list wins
    #[arg(
        long = "sjdbGTFtagExonParentGeneType",
        default_values_t = vec!["gene_type".to_string(), "gene_biotype".to_string()],
        num_args = 1..
    )]
    pub sjdb_gtf_tag_exon_parent_gene_type: Vec<String>,

    /// Overhang length for splice junction database
    #[arg(long = "sjdbOverhang", default_value_t = 100)]
    pub sjdb_overhang: u32,

    /// Extra score for alignments crossing annotated junctions
    #[arg(long = "sjdbScore", default_value_t = 2)]
    pub sjdb_score: i32,

    // ── Quantification ──────────────────────────────────────────────────
    /// Quantification mode(s): GeneCounts, TranscriptomeSAM, or empty for none.
    /// Space-separated, e.g. `--quantMode GeneCounts`.
    #[arg(long = "quantMode", num_args = 0..)]
    pub quant_mode: Vec<String>,

    /// Output format variants for `--quantMode TranscriptomeSAM`:
    ///   * `BanSingleEnd_BanIndels_ExtendSoftclip` (default, RSEM-compatible)
    ///   * `BanSingleEnd` — keep indels and soft-clips
    ///   * `BanSingleEnd_ExtendSoftclip` — keep indels, extend soft-clips
    #[arg(
        long = "quantTranscriptomeSAMoutput",
        default_value = "BanSingleEnd_BanIndels_ExtendSoftclip"
    )]
    pub quant_transcriptome_sam_output: crate::quant::transcriptome::QuantTranscriptomeSAMoutput,

    // ── Two-pass ────────────────────────────────────────────────────────
    /// Two-pass mode: None or Basic
    #[arg(long = "twopassMode", default_value = "None")]
    pub twopass_mode: TwopassMode,

    /// Reads to process in first pass; -1 = all
    #[arg(long = "twopass1readsN", default_value_t = -1, allow_hyphen_values = true)]
    pub twopass1_reads_n: i64,

    // ── WASP allele-specific filtering ──────────────────────────────────
    /// WASP output mode: None or SAMtag (emit the vW tag)
    #[arg(long = "waspOutputMode", default_value = "None")]
    pub wasp_output_mode: WaspOutputMode,

    /// VCF of heterozygous SNVs for WASP filtering (required with --waspOutputMode SAMtag)
    #[arg(long = "varVCFfile")]
    pub var_vcf_file: Option<PathBuf>,

    // ── Chimeric ────────────────────────────────────────────────────────
    // ── Debug ───────────────────────────────────────────────────────
    /// Filter for debug logging: only log detailed alignment info for reads matching this name
    #[arg(long = "readNameFilter", default_value = "")]
    pub read_name_filter: String,

    // ── Chimeric ────────────────────────────────────────────────────────
    /// Min chimeric segment length; 0 = disable chimeric detection
    #[arg(long = "chimSegmentMin", default_value_t = 0)]
    pub chim_segment_min: u32,

    /// Min total chimeric score
    #[arg(long = "chimScoreMin", default_value_t = 0)]
    pub chim_score_min: i32,

    /// Max drop in chimeric score vs read length (chimericDetectionOld)
    #[arg(long = "chimScoreDropMax", default_value_t = 20)]
    pub chim_score_drop_max: i32,

    /// Min score separation for unique chimeric alignment
    #[arg(long = "chimScoreSeparation", default_value_t = 10)]
    pub chim_score_separation: i32,

    /// Max multimapping of main chimeric segment
    #[arg(long = "chimMainSegmentMultNmax", default_value_t = 10)]
    pub chim_main_segment_mult_nmax: u32,

    /// Max read gap between chimeric segments
    #[arg(long = "chimSegmentReadGapMax", default_value_t = 0)]
    pub chim_segment_read_gap_max: u32,

    /// Min overhang at chimeric junction
    #[arg(long = "chimJunctionOverhangMin", default_value_t = 20)]
    pub chim_junction_overhang_min: u32,

    /// Score penalty for non-GT/AG chimeric junction
    #[arg(long = "chimScoreJunctionNonGTAG", default_value_t = -1, allow_hyphen_values = true)]
    pub chim_score_junction_non_gtag: i32,

    /// Chimeric output type
    #[arg(long = "chimOutType", num_args = 1..=2, default_values_t = vec!["Junctions".to_string()])]
    pub chim_out_type: Vec<String>,

    // ── STARsolo (single-cell) ──────────────────────────────────────────
    /// Single-cell barcode geometry; `None` disables solo processing.
    #[arg(long = "soloType", default_value = "None")]
    pub solo_type: SoloType,

    /// Cell-barcode whitelist file (one barcode per line, plain or gzipped).
    /// The literal `None` means "no whitelist" (all observed barcodes kept).
    /// Multiple files are allowed for `CB_UMI_Complex` (one per CB segment).
    #[arg(long = "soloCBwhitelist", num_args = 1.., default_values_t = vec!["None".to_string()])]
    pub solo_cb_whitelist: Vec<String>,

    /// 10x chemistry preset. Sets the CB/UMI geometry so it does not have to be
    /// spelled out with `--soloCBstart` and friends. `-` (the default) leaves
    /// the geometry to those flags.
    #[arg(long = "soloChemistry", default_value = "-")]
    pub solo_chemistry: String,

    /// How cell barcodes are represented: `Sequence` (2-bit packed ACGT, the
    /// default).
    #[arg(long = "soloCBtype", default_value = "Sequence")]
    pub solo_cb_type: String,

    /// 1-based start position of the cell barcode in the barcode read.
    #[arg(long = "soloCBstart", default_value_t = 1)]
    pub solo_cb_start: u32,

    /// Length of the cell barcode in bases.
    #[arg(long = "soloCBlen", default_value_t = 16)]
    pub solo_cb_len: u32,

    /// 1-based start position of the UMI in the barcode read.
    #[arg(long = "soloUMIstart", default_value_t = 17)]
    pub solo_umi_start: u32,

    /// Length of the UMI in bases (10x v2 = 10, v3 = 12).
    #[arg(long = "soloUMIlen", default_value_t = 10)]
    pub solo_umi_len: u32,

    /// Which mate carries the barcode (CB+UMI): 0 = a separate barcode read
    /// (default, 3' 10x); 1 = the barcode is a prefix of mate 1, which also
    /// carries cDNA (5' 10x, paired-end — both mates are aligned). 2 is not yet
    /// supported.
    #[arg(long = "soloBarcodeMate", default_value_t = 0)]
    pub solo_barcode_mate: u32,

    /// Barcode-read length check: 1 = require the barcode read length to equal
    /// soloCBlen + soloUMIlen; 0 = do not check (needed when the barcode read is
    /// longer, e.g. the 5' mate-1 read that continues into cDNA).
    #[arg(long = "soloBarcodeReadLength", default_value_t = 1)]
    pub solo_barcode_read_length: i64,

    /// Accepted for STAR-command compatibility; rustar always mmaps the index, so
    /// STAR's shared-memory genome-loading modes are a no-op here.
    #[arg(long = "genomeLoad", default_value = "NoSharedMemory")]
    pub genome_load: String,

    /// `CB_UMI_Complex` cell-barcode segment positions, one per segment, as
    /// `startAnchor_startDist_endAnchor_endDist`. Only read-start anchoring
    /// (`anchor = 0`, fixed positions) is supported, e.g. `0_0_0_7 0_8_0_15`.
    #[arg(long = "soloCBposition", num_args = 0..)]
    pub solo_cb_position: Vec<String>,

    /// `CB_UMI_Complex` UMI position as `startAnchor_startDist_endAnchor_endDist`
    /// (read-start anchoring only), e.g. `0_16_0_25`.
    #[arg(long = "soloUMIposition", default_value = "")]
    pub solo_umi_position: String,

    /// Genomic features to quantify per cell: Gene, GeneFull, SJ, Velocyto, …
    #[arg(long = "soloFeatures", num_args = 1.., default_values_t = vec!["Gene".to_string()])]
    pub solo_features: Vec<String>,

    /// UMI collapsing strategy: 1MM_All, 1MM_Directional, 1MM_Directional_UMItools,
    /// Exact, or NoDedup.
    #[arg(long = "soloUMIdedup", num_args = 1.., default_values_t = vec!["1MM_All".to_string()])]
    pub solo_umi_dedup: Vec<String>,

    /// Cell-barcode-to-whitelist matching: Exact, 1MM, 1MM_multi,
    /// 1MM_multi_pseudocounts, 1MM_multi_Nbase_pseudocounts.
    #[arg(long = "soloCBmatchWLtype", default_value = "1MM_multi")]
    pub solo_cb_match_wl_type: String,

    /// Cell-calling / matrix filtering: None, CellRanger2.2, EmptyDrops_CR, TopCells.
    #[arg(long = "soloCellFilter", num_args = 1.., default_values_t = vec!["CellRanger2.2".to_string(), "3000".to_string(), "0.99".to_string(), "10".to_string()])]
    pub solo_cell_filter: Vec<String>,

    /// Counting method for reads mapping to multiple genes: Unique (default,
    /// drop), Uniform, Rescue, PropUnique, EM. Non-Unique methods additionally
    /// write `UniqueAndMult-<method>.mtx` (real-valued) per Gene/GeneFull feature.
    #[arg(long = "soloMultiMappers", num_args = 1.., default_values_t = vec!["Unique".to_string()])]
    pub solo_multi_mappers: Vec<String>,

    /// Third column of `features.tsv`. STAR's default is `Gene Expression`;
    /// the sentinel `-` suppresses the column entirely.
    #[arg(
        long = "soloOutFormatFeaturesGeneField3",
        default_value = "Gene Expression"
    )]
    pub solo_out_format_features_gene_field3: String,

    /// Output directory name for solo matrices (relative to `--outFileNamePrefix`).
    #[arg(long = "soloOutFileNames", num_args = 1.., default_values_t = vec!["Solo.out/".to_string(), "features.tsv".to_string(), "barcodes.tsv".to_string(), "matrix.mtx".to_string()])]
    pub solo_out_file_names: Vec<String>,

    /// Gzip the solo `matrix.mtx` / `barcodes.tsv` / `features.tsv` and append a
    /// `.gz` suffix (CellRanger-style output). Default `no` keeps the plain files
    /// that STARsolo writes (so the byte-for-byte STARsolo comparison still holds).
    #[arg(long = "soloOutGzip", default_value = "no")]
    pub solo_out_gzip: String,

    /// Velocyto ambiguous-molecule handling (rustar extension beyond STARsolo).
    /// `yes` (default) writes the three `spliced`/`unspliced`/`ambiguous` matrices
    /// like STARsolo — exon-only molecules with no junction/intron evidence stay in
    /// `ambiguous`. `no` resolves those molecules to `spliced` (an exon-only read is
    /// most likely mature mRNA; cf. He, Soneson & Patro 2023) and writes only
    /// `spliced`/`unspliced`, with no `ambiguous.mtx`.
    #[arg(long = "soloVelocytoAmbiguous", default_value = "yes")]
    pub solo_velocyto_ambiguous: String,

    /// Strand of the read relative to the gene for counting: Forward, Reverse, Unstranded.
    #[arg(long = "soloStrand", default_value = "Forward")]
    pub solo_strand: String,

    /// UMI filtering of multi-gene UMIs: `-`/`None` (default, no filtering),
    /// `MultiGeneUMI`, `MultiGeneUMI_CR`, or `MultiGeneUMI_All`. The `_CR`
    /// variant matches CellRanger > 3.0.
    #[arg(long = "soloUMIfiltering", num_args = 1.., default_values_t = vec!["-".to_string()])]
    pub solo_umi_filtering: Vec<String>,
    /// `Chimeric.out.junction` format: `0` (plain, default) or `1` (append a STAR-Fusion-style
    /// comment header with the command line and read counts)
    #[arg(long = "chimOutJunctionFormat", default_value_t = 0)]
    pub chim_out_junction_format: u8,

    /// Full command line as invoked, embedded in the BAM `@PG` `CL:` field.
    #[arg(skip)]
    pub command_line: Option<String>,
}

impl Parameters {
    /// Build an output path by concatenating `suffix` onto `out_file_name_prefix`.
    pub fn output_path(&self, suffix: &str) -> PathBuf {
        PathBuf::from(format!("{}{suffix}", self.out_file_name_prefix))
    }

    /// Whether the run produces per-read alignment records (SAM/BAM). False only
    /// for `--outSAMtype None` written to a file (no `--outStd`): the alignment
    /// loops then skip building SAM records entirely, which is a large saving for
    /// solo / quant-only runs that only need the count matrix.
    pub fn emits_alignments(&self) -> bool {
        !matches!(self.out_std, OutStd::None) || self.out_sam_type.format != OutSamFormat::None
    }

    /// Whether `--chimOutType` includes `Junctions` (write Chimeric.out.junction).
    pub fn chim_out_junctions(&self) -> bool {
        self.chim_out_type.iter().any(|s| s == "Junctions")
    }

    /// Whether `--chimOutType` includes `WithinBAM` (write supplementary BAM records).
    pub fn chim_out_within_bam(&self) -> bool {
        self.chim_out_type.iter().any(|s| s == "WithinBAM")
    }

    /// True if the user provided a non-default `--outSAMattrRGline`.
    pub fn rg_line_set(&self) -> bool {
        !self.out_sam_attr_rg_line.is_empty() && self.out_sam_attr_rg_line[0] != "-"
    }

    /// Parse `--outSAMattrRGline` into one tab-joined body per `@RG` block.
    ///
    /// Mirrors `Parameters_readFilesInit.cpp:65-82`: tokens are split on bare
    /// `,` separators, and each resulting block's first token must begin with
    /// `ID:`. An empty block (adjacent commas or a trailing comma) is an error.
    pub fn parsed_rg_lines(&self) -> Result<Vec<String>, crate::error::Error> {
        if !self.rg_line_set() {
            return Ok(Vec::new());
        }
        self.out_sam_attr_rg_line
            .split(|tok| tok == ",")
            .map(|block| {
                let first = block.first().ok_or_else(|| {
                    crate::error::Error::Parameter(
                        "--outSAMattrRGline: empty RG block".into(),
                    )
                })?;
                if !first.starts_with("ID:") {
                    return Err(crate::error::Error::Parameter(format!(
                        "--outSAMattrRGline: first field of each RG line must start with 'ID:', got '{first}'"
                    )));
                }
                Ok(block.join("\t"))
            })
            .collect()
    }

    /// Read group ID emitted on SAM records (the first RG line's `ID:` value).
    /// Returns `None` when no RG line is configured.
    pub fn primary_rg_id(&self) -> Result<Option<String>, crate::error::Error> {
        Ok(self.parsed_rg_lines()?.first().and_then(|body| {
            body.split('\t')
                .next()?
                .strip_prefix("ID:")
                .map(str::to_owned)
        }))
    }

    /// Per-file read group ID, replicated from a single RG line if needed.
    /// Returns empty vec when no RG line is set.
    pub fn rg_ids(&self) -> Result<Vec<String>, crate::error::Error> {
        let lines = self.parsed_rg_lines()?;
        if lines.is_empty() {
            return Ok(Vec::new());
        }
        let ids: Vec<String> = lines
            .iter()
            .map(|body| {
                let first = body.split('\t').next().unwrap_or("");
                first.trim_start_matches("ID:").to_string()
            })
            .collect();
        let n_files = self.read_files_in.len().max(1);
        if ids.len() > 1 && ids.len() != n_files {
            return Err(crate::error::Error::Parameter(format!(
                "--outSAMattrRGline: {} RG entries does not match --readFilesIn count {} (must be 1 or N)",
                ids.len(),
                n_files
            )));
        }
        if ids.len() == 1 {
            Ok(vec![ids[0].clone(); n_files])
        } else {
            Ok(ids)
        }
    }

    /// Compute the default window distance: 2^winBinNbits * winAnchorDistNbins.
    /// Used for max_cluster_dist and as default alignIntronMax (when 0).
    pub fn win_bin_window_dist(&self) -> u64 {
        (1u64 << self.win_bin_nbits) * self.win_anchor_dist_nbins as u64
    }

    /// Redefine window parameters based on genome size and intron/gap limits.
    /// Ports STAR's Genome_genomeLoad.cpp logic that recomputes winBinNbits,
    /// winFlankNbins, and winAnchorDistNbins after loading the genome.
    ///
    /// IMPORTANT: winBinNbits is only redefined when alignIntronMax > 0 OR
    /// alignMatesGapMax > 0. When both are 0, winBinNbits stays at its default (16).
    pub fn redefine_window_params(&mut self, n_genome: u64) {
        let intron_max = self.align_intron_max as u64;
        let gap_max = self.align_mates_gap_max as u64;

        if intron_max == 0 && gap_max == 0 {
            // STAR: no redefinition when both are 0. Log effective max intron.
            let max_intron = (1u64 << self.win_bin_nbits) * self.win_anchor_dist_nbins as u64;
            log::info!(
                "alignIntronMax=alignMatesGapMax=0, max intron ~= (2^winBinNbits)*winAnchorDistNbins={max_intron}"
            );
            return;
        }

        // STAR: max(max(4, alignIntronMax), alignMatesGapMax==0 ? 1000 : alignMatesGapMax)
        let max_span = std::cmp::max(
            std::cmp::max(4u64, intron_max),
            if gap_max == 0 { 1000 } else { gap_max },
        );

        // winBinNbits = floor(log2(max_span / 4) + 0.5)
        self.win_bin_nbits = ((max_span as f64 / 4.0).log2() + 0.5).floor() as u32;

        // max with genome-based value: floor(log2(nGenome/40000 + 1) + 0.5)
        let genome_based = ((n_genome as f64 / 40000.0 + 1.0).log2() + 0.5).floor() as u32;
        self.win_bin_nbits = self.win_bin_nbits.max(genome_based);

        // Cap at genomeChrBinNbits
        if self.win_bin_nbits > self.genome_chr_bin_nbits {
            self.win_bin_nbits = self.genome_chr_bin_nbits;
        }

        // Redefine winFlankNbins and winAnchorDistNbins
        let max_gap = std::cmp::max(intron_max, gap_max);
        self.win_flank_nbins = (max_gap / (1u64 << self.win_bin_nbits) + 1) as u32;
        self.win_anchor_dist_nbins = 2 * self.win_flank_nbins;

        log::info!(
            "Redefined window params: winBinNbits={}, winFlankNbins={}, winAnchorDistNbins={}",
            self.win_bin_nbits,
            self.win_flank_nbins,
            self.win_anchor_dist_nbins
        );
    }

    /// Parse and validate parameter combinations.
    pub fn parse() -> Self {
        Self::parse_from(std::env::args_os())
    }

    /// Parse and validate parameter combinations from args.
    pub fn parse_from<T: Into<std::ffi::OsString> + Clone>(
        args: impl IntoIterator<Item = T>,
    ) -> Self {
        Self::try_parse_from(args).unwrap_or_else(|e| {
            if cfg!(test) {
                panic!("{e}")
            } else {
                e.format(&mut <Self as CommandFactory>::command()).exit()
            }
        })
    }

    /// Parse and validate parameter combinations.
    pub fn try_parse() -> Result<Self, clap::Error> {
        Self::try_parse_from(std::env::args_os())
    }

    /// Parse and validate parameter combinations from args.
    pub fn try_parse_from<T: Into<std::ffi::OsString> + Clone>(
        args: impl IntoIterator<Item = T>,
    ) -> Result<Self, clap::Error> {
        use clap::error::ErrorKind;

        let args: Vec<_> = args.into_iter().map(Into::into).collect();

        let mut command = <Self as clap::CommandFactory>::command();
        let matches = command.clone().get_matches_from(args.iter());
        let mut params = <Self as clap::FromArgMatches>::from_arg_matches(&matches)?;

        params.command_line = {
            let args: Vec<_> = args.iter().map(|s| s.to_string_lossy()).collect();
            shlex::try_join(args.iter().map(AsRef::as_ref)).ok()
        };

        // genomeGenerate requires FASTA files
        if params.run_mode == RunMode::GenomeGenerate && params.genome_fasta_files.is_empty() {
            return Err(command.error(
                ErrorKind::MissingRequiredArgument,
                "--genomeFastaFiles is required when --runMode genomeGenerate",
            ));
        }

        // Sparse suffix array stride must be >= 1 (1 = dense, STAR default).
        if params.genome_sa_sparse_d == 0 {
            return Err(command.error(
                ErrorKind::ValueValidation,
                "--genomeSAsparseD must be >= 1 (1 = dense suffix array)",
            ));
        }

        // alignReads requires read files — except SmartSeq, which gets its reads
        // from --readFilesManifest instead.
        if params.run_mode == RunMode::AlignReads
            && params.read_files_in.is_empty()
            && params.solo_type != SoloType::SmartSeq
        {
            return Err(command.error(
                ErrorKind::MissingRequiredArgument,
                "--readFilesIn is required when --runMode alignReads",
            ));
        }

        // --genomeTransformType: Haploid and Diploid are implemented. Both require
        // a VCF, and are incompatible with a GTF (STAR itself doesn't combine
        // genomeTransform with sjdb annotation at genomeGenerate).
        if !params.genome_transform_type.eq_ignore_ascii_case("None") {
            if !params.genome_transform_type.eq_ignore_ascii_case("Haploid")
                && !params.genome_transform_type.eq_ignore_ascii_case("Diploid")
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "--genomeTransformType must be None, Haploid, or Diploid",
                ));
            }
            if params.genome_transform_vcf.is_none() {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--genomeTransformType Haploid/Diploid requires --genomeTransformVCF",
                ));
            }
            if params.sjdb_gtf_file.is_some() {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "--genomeTransformType is incompatible with --sjdbGTFfile",
                ));
            }
        }

        // inputAlignmentsFromBAM: only --bamRemoveDuplicatesType is implemented so far
        if params.run_mode == RunMode::InputAlignmentsFromBAM {
            let dedup = params.bam_remove_duplicates_type.as_str();
            if dedup == "-" {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--runMode inputAlignmentsFromBAM requires --bamRemoveDuplicatesType \
                     (UniqueIdentical or UniqueIdenticalNotMulti)",
                ));
            }
            if !dedup.eq_ignore_ascii_case("UniqueIdentical")
                && !dedup.eq_ignore_ascii_case("UniqueIdenticalNotMulti")
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "unknown --bamRemoveDuplicatesType {dedup}; expected UniqueIdentical or \
                         UniqueIdenticalNotMulti"
                    ),
                ));
            }
        }

        // liftOver requires a chain file and a GTF to lift
        if params.run_mode == RunMode::LiftOver {
            if params.genome_chain_files.is_empty() {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--genomeChainFiles is required when --runMode liftOver",
                ));
            }
            if params.sjdb_gtf_file.is_none() {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--sjdbGTFfile is required when --runMode liftOver",
                ));
            }
        }

        // --outWigType at alignReads: only `bedGraph` (stranded, full-length) is
        // implemented. `wiggle`, `--outWigStrand Unstranded`, and the 2nd word
        // (`read1_5p`/`read2`) are STAR features of `--runMode
        // inputAlignmentsFromBAM`, which rustar-aligner doesn't have; reject them
        // loudly rather than silently emitting bedGraph / stranded / full-length
        // tracks instead.
        if params
            .out_wig_type
            .iter()
            .any(|t| !t.eq_ignore_ascii_case("None"))
        {
            if params
                .out_wig_type
                .iter()
                .any(|t| t.eq_ignore_ascii_case("wiggle"))
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "--outWigType wiggle is not implemented; use --outWigType bedGraph",
                ));
            }
            if params.out_wig_strand.eq_ignore_ascii_case("Unstranded") {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "--outWigStrand Unstranded is not implemented; use --outWigStrand Stranded",
                ));
            }
            if params.out_wig_type.len() > 1 {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "the --outWigType 2nd word (read1_5p / read2) is not implemented; omit it",
                ));
            }
        }

        // quantMode GeneCounts requires a GTF file
        if params.quant_gene_counts() && params.sjdb_gtf_file.is_none() {
            return Err(command.error(
                ErrorKind::MissingRequiredArgument,
                "--quantMode GeneCounts requires --sjdbGTFfile",
            ));
        }

        // Read group: `RG` in outSAMattributes without an RG line is a fatal
        // error (STAR: Parameters_samAttributes.cpp:206). STAR's "All" preset
        // does NOT include RG, so only match a literal user-supplied RG flag.
        if params.out_sam_attributes.contains(SamAttributes::RG) && !params.rg_line_set() {
            return Err(command.error(
                ErrorKind::MissingRequiredArgument,
                "--outSAMattributes contains RG tag, but --outSAMattrRGline is not set",
            ));
        }
        params
            .rg_ids()
            .map_err(|e| command.error(ErrorKind::InvalidValue, e))?;

        // Fold runtime-derived bits into out_sam_attributes so writers can read
        // it directly without re-computing. Mirrors STAR Parameters_samAttributes.cpp.
        if params.rg_line_set() {
            params.out_sam_attributes |= SamAttributes::RG;
        }
        // XS is only emitted in intronMotif mode (Parameters_samAttributes.cpp:172-179).
        // intronMotif forces XS on; anything else strips XS even if explicitly listed.
        if params.out_sam_strand_field == "intronMotif" {
            log::info!(
                "--outSAMstrandField=intronMotif, therefore rustar-aligner will output XS attribute"
            );
            params.out_sam_attributes |= SamAttributes::XS;
        } else {
            params.out_sam_attributes.remove(SamAttributes::XS);
        }

        // --alignEndsType: reject unknown/unimplemented values up front.
        if let Err(msg) = params.align_ends_type.parse::<AlignEndsType>() {
            return Err(command.error(ErrorKind::InvalidValue, msg));
        }

        // --outSAMorder: rustar-aligner always preserves input order, so both STAR
        // values are accepted and behave identically. Reject anything else.
        if params.out_sam_order != "Paired" && params.out_sam_order != "PairedKeepInputOrder" {
            return Err(command.error(
                ErrorKind::InvalidValue,
                format!(
                    "unknown --outSAMorder '{}'; expected 'Paired' or 'PairedKeepInputOrder'",
                    params.out_sam_order
                ),
            ));
        }

        // quantMode TranscriptomeSAM requires transcript annotations —
        // either via --sjdbGTFfile or pre-generated transcriptInfo.tab
        // et al in --genomeDir (persisted at genomeGenerate time). At
        // validation time we can only enforce the genomeGenerate rule;
        // for alignReads, GenomeIndex::load checks for the on-disk files
        // and surfaces a clear error if neither source is available.
        if params.run_mode == RunMode::GenomeGenerate
            && params.quant_transcriptome_sam()
            && params.sjdb_gtf_file.is_none()
        {
            return Err(command.error(
                ErrorKind::MissingRequiredArgument,
                "--quantMode TranscriptomeSAM requires --sjdbGTFfile at genomeGenerate",
            ));
        }

        // --soloChemistry presets. STAR applies these before validating the
        // geometry, so a preset and an explicit --soloCBstart cannot disagree:
        // the preset wins and says so, rather than half-applying.
        if params.solo_chemistry != "-" {
            let geometry = match params.solo_chemistry.as_str() {
                // (cb_start, cb_len, umi_start, umi_len), 1-based as STAR takes them.
                "CR_2" | "CR_3" | "CR_3.1" | "CR_4" | "SC3Pv1" => Some((1, 14, 15, 10)),
                "SC3Pv2" => Some((1, 16, 17, 10)),
                "SC3Pv3" | "SC3Pv4" | "SC5P" => Some((1, 16, 17, 12)),
                _ => None,
            };
            let Some((cb_start, cb_len, umi_start, umi_len)) = geometry else {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "unknown --soloChemistry '{}'; expected SC3Pv1, SC3Pv2, SC3Pv3, \
                         SC3Pv4, SC5P, or -",
                        params.solo_chemistry
                    ),
                ));
            };
            for (flag, given) in [
                ("--soloCBstart", params.solo_cb_start),
                ("--soloCBlen", params.solo_cb_len),
                ("--soloUMIstart", params.solo_umi_start),
                ("--soloUMIlen", params.solo_umi_len),
            ] {
                if given != 0 {
                    log::warn!(
                        "--soloChemistry {} overrides the geometry given by {flag}",
                        params.solo_chemistry
                    );
                }
            }
            params.solo_cb_start = cb_start;
            params.solo_cb_len = cb_len;
            params.solo_umi_start = umi_start;
            params.solo_umi_len = umi_len;
        }

        // --soloCBtype: only 2-bit packed ACGT barcodes are supported. `String`
        // needs a different whitelist representation entirely, so refuse rather
        // than silently pack a non-ACGT barcode into nonsense.
        if params.solo_cb_type != "Sequence" {
            return Err(command.error(
                ErrorKind::InvalidValue,
                format!(
                    "--soloCBtype {} is not supported; expected Sequence",
                    params.solo_cb_type
                ),
            ));
        }

        // ── STARsolo validation ─────────────────────────────────────────
        if params.run_mode == RunMode::AlignReads && params.solo_enabled() {
            // CB_UMI_Complex needs one CB position + whitelist per segment.
            if params.solo_type == SoloType::CbUmiComplex {
                if params.solo_cb_position.is_empty() {
                    return Err(command.error(
                        ErrorKind::MissingRequiredArgument,
                        "--soloType CB_UMI_Complex requires --soloCBposition (one per CB segment)",
                    ));
                }
                if params.solo_cb_whitelist.len() != params.solo_cb_position.len() {
                    return Err(command.error(
                        ErrorKind::InvalidValue,
                        format!(
                            "--soloType CB_UMI_Complex: {} --soloCBposition segments but {} --soloCBwhitelist files (must match)",
                            params.solo_cb_position.len(),
                            params.solo_cb_whitelist.len()
                        ),
                    ));
                }
            }
            // SmartSeq is plate-based (one library per manifest cell, no barcodes).
            if params.solo_type == SoloType::SmartSeq && params.read_files_manifest.is_none() {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--soloType SmartSeq requires --readFilesManifest (a TSV of read1<TAB>read2<TAB>cellID per cell)",
                ));
            }
            // CB_UMI_Simple needs exactly two read files: cDNA + barcode read.
            if matches!(
                params.solo_type,
                SoloType::CbUmiSimple | SoloType::CbUmiComplex | SoloType::CbSamTagOut
            ) && params.read_files_in.len() != 2
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "--soloType {} requires exactly two --readFilesIn files (cDNA read then barcode read); got {}",
                        params.solo_type,
                        params.read_files_in.len()
                    ),
                ));
            }
            // soloBarcodeMate: 0 (separate barcode read) or 1 (barcode on mate 1,
            // 5' paired-end). Mate 2 is not yet supported.
            match params.solo_barcode_mate {
                0 => {}
                1 => {
                    if params.solo_type != SoloType::CbUmiSimple {
                        return Err(command.error(
                            ErrorKind::InvalidValue,
                            "--soloBarcodeMate 1 is only supported with --soloType CB_UMI_Simple",
                        ));
                    }
                }
                other => {
                    return Err(command.error(
                        ErrorKind::InvalidValue,
                        format!(
                            "--soloBarcodeMate {other} not supported (only 0 = separate barcode read, or 1 = barcode on mate 1)"
                        ),
                    ));
                }
            }
            // Gene / GeneFull / SJ / Velocyto are implemented.
            for f in &params.solo_features {
                if !matches!(f.as_str(), "SJ" | "Velocyto")
                    && f.parse::<crate::solo::SoloFeature>().is_err()
                {
                    return Err(command.error(
                        ErrorKind::InvalidValue,
                        format!(
                            "unsupported --soloFeatures '{f}'; supported: Gene, GeneFull, SJ, Velocyto"
                        ),
                    ));
                }
            }
            // soloMultiMappers values.
            for m in &params.solo_multi_mappers {
                if !matches!(
                    m.as_str(),
                    "Unique" | "Uniform" | "Rescue" | "PropUnique" | "EM"
                ) {
                    return Err(command.error(
                        ErrorKind::InvalidValue,
                        format!(
                            "unsupported --soloMultiMappers '{m}'; expected Unique, Uniform, Rescue, PropUnique, or EM"
                        ),
                    ));
                }
            }
            // Gene-level features need a gene model (SJ does not — junctions come
            // from the alignments).
            let needs_gtf = params
                .solo_features
                .iter()
                .any(|f| f == "Gene" || f == "GeneFull" || f == "Velocyto");
            if needs_gtf && params.sjdb_gtf_file.is_none() {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--soloFeatures Gene/GeneFull requires --sjdbGTFfile (a gene model)",
                ));
            }
            // CB length / UMI length sanity.
            if params.solo_type == SoloType::CbUmiSimple
                && (params.solo_cb_len == 0 || params.solo_umi_len == 0)
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "--soloCBlen and --soloUMIlen must be > 0 for soloType CB_UMI_Simple",
                ));
            }
            // Cell barcode cannot exceed a u64 packing (32 bases).
            if params.solo_cb_len as usize > crate::solo::whitelist::CB_LEN_MAX {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "--soloCBlen {} exceeds the maximum of {}",
                        params.solo_cb_len,
                        crate::solo::whitelist::CB_LEN_MAX
                    ),
                ));
            }
            // Validate --soloCBmatchWLtype.
            if params
                .solo_cb_match_wl_type
                .parse::<crate::solo::whitelist::CbMatchType>()
                .is_err()
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "unknown --soloCBmatchWLtype '{}'; expected Exact, 1MM, 1MM_multi, 1MM_multi_pseudocounts, or 1MM_multi_Nbase_pseudocounts",
                        params.solo_cb_match_wl_type
                    ),
                ));
            }
            // Validate --soloUMIdedup (each method string).
            for m in &params.solo_umi_dedup {
                if m.parse::<crate::solo::UmiDedup>().is_err() {
                    return Err(command.error(
                        ErrorKind::InvalidValue,
                        format!(
                            "unknown --soloUMIdedup '{m}'; expected Exact, NoDedup, 1MM_All, 1MM_Directional, or 1MM_Directional_UMItools"
                        ),
                    ));
                }
            }
            // Validate --soloUMIfiltering (each method string).
            for f in &params.solo_umi_filtering {
                if f.parse::<crate::solo::UmiFiltering>().is_err() {
                    return Err(command.error(
                        ErrorKind::InvalidValue,
                        format!(
                            "unknown --soloUMIfiltering '{f}'; expected -, None, MultiGeneUMI, MultiGeneUMI_CR, or MultiGeneUMI_All"
                        ),
                    ));
                }
            }
            // Validate --clipAdapterType.
            if !matches!(
                params.clip_adapter_type.as_str(),
                "Hamming" | "CellRanger4" | "None"
            ) {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "unknown --clipAdapterType '{}'; expected Hamming, CellRanger4, or None",
                        params.clip_adapter_type
                    ),
                ));
            }
            // Validate --soloStrand.
            if params
                .solo_strand
                .parse::<crate::solo::SoloStrand>()
                .is_err()
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    format!(
                        "unknown --soloStrand '{}'; expected Forward, Reverse, or Unstranded",
                        params.solo_strand
                    ),
                ));
            }
            // A whitelist is required for any correction beyond None (SmartSeq
            // has no cell barcodes at all, so the rule does not apply).
            if params.solo_type != SoloType::SmartSeq
                && params.solo_cb_whitelist_none()
                && params.solo_cb_match_wl_type != "Exact"
            {
                return Err(command.error(
                    ErrorKind::InvalidValue,
                    "--soloCBwhitelist None requires --soloCBmatchWLtype Exact (no correction possible without a whitelist)",
                ));
            }
        }

        // WASP SAMtag mode requires a VCF of heterozygous SNVs; fold the vW bit
        // into out_sam_attributes so the writer emits it (vA/vG stay opt-in).
        if params.wasp_output_mode == WaspOutputMode::SAMtag {
            if params.var_vcf_file.is_none() {
                return Err(command.error(
                    ErrorKind::MissingRequiredArgument,
                    "--waspOutputMode SAMtag requires --varVCFfile",
                ));
            }
            params.out_sam_attributes |= SamAttributes::VW;
        }

        Ok(params)
    }

    /// Returns true if `--quantMode GeneCounts` was requested.
    pub fn quant_gene_counts(&self) -> bool {
        self.quant_mode.iter().any(|m| m == "GeneCounts")
    }

    /// Returns true if `--quantMode TranscriptomeSAM` was requested.
    pub fn quant_transcriptome_sam(&self) -> bool {
        self.quant_mode.iter().any(|m| m == "TranscriptomeSAM")
    }

    /// True when a single-cell run is requested (`--soloType` != None).
    pub fn solo_enabled(&self) -> bool {
        self.solo_type != SoloType::None
    }

    /// Path to the cDNA (transcript) read file. For solo runs this is the
    /// FIRST `--readFilesIn` file (STAR convention: `cDNA_read barcode_read`).
    /// Returns `None` if no read files are configured.
    pub fn cdna_read_file(&self) -> Option<&PathBuf> {
        self.read_files_in.first()
    }

    /// Path to the barcode (CB+UMI) read file — the SECOND `--readFilesIn`
    /// file when solo is enabled. `None` if absent.
    pub fn barcode_read_file(&self) -> Option<&PathBuf> {
        if self.solo_enabled() {
            self.read_files_in.get(1)
        } else {
            None
        }
    }

    /// True for a 5' paired-end solo run (`--soloBarcodeMate 1`): the barcode is a
    /// prefix of mate 1 and both `--readFilesIn` files are cDNA mates.
    pub fn solo_barcode_on_mate1(&self) -> bool {
        self.solo_enabled() && self.solo_barcode_mate == 1
    }

    /// The two cDNA mate files (mate 1, mate 2) for a `--soloBarcodeMate 1` run.
    pub fn solo_cdna_mate_files(&self) -> Option<(&PathBuf, &PathBuf)> {
        match (self.read_files_in.first(), self.read_files_in.get(1)) {
            (Some(m1), Some(m2)) => Some((m1, m2)),
            _ => None,
        }
    }

    /// Bases to clip from the 5' end of `mate` (0 or 1). A single configured
    /// value applies to both mates.
    pub fn clip5p(&self, mate: usize) -> usize {
        let v = &self.clip5p_nbases;
        v[mate.min(v.len().saturating_sub(1))] as usize
    }

    /// Bases to clip from the 3' end of `mate` (0 or 1). A single configured
    /// value applies to both mates.
    pub fn clip3p(&self, mate: usize) -> usize {
        let v = &self.clip3p_nbases;
        v[mate.min(v.len().saturating_sub(1))] as usize
    }

    /// Returns true if `--outWigType bedGraph` was requested.
    pub fn out_wig_bedgraph(&self) -> bool {
        self.out_wig_type
            .iter()
            .any(|t| t.eq_ignore_ascii_case("bedGraph"))
    }

    /// Returns true unless `--outWigNorm None` was requested (RPM is the default).
    pub fn out_wig_rpm(&self) -> bool {
        !self.out_wig_norm.eq_ignore_ascii_case("None")
    }

    /// True when the literal `None` whitelist was given (keep all barcodes).
    pub fn solo_cb_whitelist_none(&self) -> bool {
        self.solo_cb_whitelist.len() == 1 && self.solo_cb_whitelist[0] == "None"
    }

    /// Path to the (first) cell-barcode whitelist file, or `None` for the
    /// literal `None` whitelist.
    pub fn solo_cb_whitelist_path(&self) -> Option<PathBuf> {
        if self.solo_cb_whitelist_none() {
            None
        } else {
            self.solo_cb_whitelist.first().map(PathBuf::from)
        }
    }

    /// Parsed `--soloCBmatchWLtype` flags. Falls back to the `1MM_multi`
    /// default if somehow unset (validation rejects invalid strings).
    pub fn solo_cb_match_type(&self) -> crate::solo::whitelist::CbMatchType {
        self.solo_cb_match_wl_type
            .parse()
            .unwrap_or(crate::solo::whitelist::CbMatchType {
                mm1: true,
                mm1_multi: true,
                mm1_multi_nbase: false,
                pseudocounts: false,
            })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: parse a STAR-style command line (without program name).
    fn try_parse(args: &[&str]) -> Result<Parameters, clap::Error> {
        let mut full = vec!["rustar-aligner"];
        full.extend_from_slice(args);
        Parameters::try_parse_from(&full)
    }

    #[test]
    fn defaults() {
        let p = try_parse(&["--readFilesIn", "reads.fq"]).unwrap();
        assert_eq!(p.run_mode, RunMode::AlignReads);
        assert_eq!(p.run_thread_n, NonZeroUsize::new(1).unwrap());
        assert_eq!(p.run_rng_seed, 777);
        assert_eq!(p.genome_dir, PathBuf::from("./GenomeDir"));
        assert_eq!(p.genome_sa_index_nbases, 14);
        assert_eq!(p.genome_chr_bin_nbits, 18);
        assert_eq!(p.genome_sa_sparse_d, 1);
        assert_eq!(p.read_map_number, -1);
        assert_eq!(p.clip5p(0), 0);
        assert_eq!(p.clip3p(0), 0);
        assert_eq!(p.out_file_name_prefix, "./");
        assert_eq!(p.out_sam_type, OutSamType::default());
        assert_eq!(p.out_sam_strand_field, "None");
        assert_eq!(p.out_sam_attributes, SamAttributes::STANDARD);
        assert_eq!(p.out_sam_unmapped, OutSamUnmapped::None);
        assert_eq!(p.out_sam_mapq_unique, 255);
        assert_eq!(p.out_sam_mult_nmax, -1);
        assert_eq!(p.out_filter_type, OutFilterType::Normal);
        assert_eq!(p.out_filter_multimap_nmax, 10);
        assert_eq!(p.out_filter_mismatch_nmax, 10);
        assert!((p.out_filter_mismatch_nover_lmax - 0.3).abs() < f64::EPSILON);
        assert_eq!(p.out_filter_score_min, 0);
        assert!((p.out_filter_score_min_over_lread - 0.66).abs() < f64::EPSILON);
        assert_eq!(p.out_filter_match_nmin, 0);
        assert!((p.out_filter_match_nmin_over_lread - 0.66).abs() < f64::EPSILON);
        assert_eq!(p.align_intron_min, 21);
        assert_eq!(p.align_intron_max, 0);
        assert_eq!(p.align_mates_gap_max, 0);
        assert_eq!(p.align_sj_overhang_min, 5);
        assert_eq!(p.align_sjdb_overhang_min, 3);
        assert_eq!(p.align_sj_stitch_mismatch_nmax, vec![0, -1, 0, 0]);
        assert_eq!(p.score_gap, 0);
        assert_eq!(p.score_gap_noncan, -8);
        assert_eq!(p.score_gap_gcag, -4);
        assert_eq!(p.score_gap_atac, -8);
        assert_eq!(p.score_del_open, -2);
        assert_eq!(p.score_del_base, -2);
        assert_eq!(p.score_ins_open, -2);
        assert_eq!(p.score_ins_base, -2);
        assert_eq!(p.score_stitch_sj_shift, 1);
        assert_eq!(p.seed_multimap_nmax, 10000);
        assert_eq!(p.seed_per_read_nmax, 1000);
        assert_eq!(p.seed_per_window_nmax, 50);
        assert_eq!(p.seed_search_start_lmax, 50);
        assert!((p.seed_search_start_lmax_over_lread - 1.0).abs() < f64::EPSILON);
        assert_eq!(p.seed_search_lmax, 0);
        assert_eq!(p.seed_map_min, 5);
        assert_eq!(p.win_anchor_multimap_nmax, 50);
        assert_eq!(p.align_windows_per_read_nmax, 10000);
        assert_eq!(p.align_transcripts_per_window_nmax, 100);
        assert!((p.win_read_coverage_relative_min - 0.5).abs() < f64::EPSILON);
        assert_eq!(p.win_bin_nbits, 16);
        assert_eq!(p.win_anchor_dist_nbins, 9);
        assert_eq!(p.win_flank_nbins, 4);
        assert!(p.sjdb_gtf_file.is_none());
        assert_eq!(p.sjdb_overhang, 100);
        assert_eq!(p.sjdb_score, 2);
        assert_eq!(p.twopass_mode, TwopassMode::None);
        assert_eq!(p.twopass1_reads_n, -1);
        assert_eq!(p.chim_segment_min, 0);
        assert_eq!(p.chim_score_min, 0);
        assert_eq!(p.chim_out_type, vec!["Junctions".to_string()]);
        assert_eq!(
            p.out_filter_intron_strands,
            IntronStrandFilter::RemoveInconsistentStrands
        );
        assert_eq!(p.out_sj_filter_overhang_min, vec![30, 12, 12, 12]);
        assert_eq!(p.out_sj_filter_count_unique_min, vec![3, 1, 1, 1]);
        assert_eq!(p.out_sj_filter_count_total_min, vec![3, 1, 1, 1]);
        assert_eq!(p.out_sj_filter_dist_to_other_sjmin, vec![10, 0, 5, 10]);
        assert_eq!(
            p.out_sj_filter_intron_max_vs_read_n,
            vec![50_000, 100_000, 200_000]
        );
    }

    #[test]
    fn genome_generate_mode() {
        let p = try_parse(&[
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            "/data/genome",
            "--genomeFastaFiles",
            "chr1.fa",
            "chr2.fa",
            "--runThreadN",
            "8",
            "--tempDir",
            "/scratch2/tmp",
            "--genomeSAindexNbases",
            "11",
        ])
        .unwrap();
        assert_eq!(p.run_mode, RunMode::GenomeGenerate);
        assert_eq!(p.genome_dir, PathBuf::from("/data/genome"));
        assert_eq!(
            p.genome_fasta_files,
            vec![PathBuf::from("chr1.fa"), PathBuf::from("chr2.fa")]
        );
        assert_eq!(p.run_thread_n, NonZeroUsize::new(8).unwrap());
        assert_eq!(p.temp_dir, Some(PathBuf::from("/scratch2/tmp")));
        assert_eq!(p.genome_sa_index_nbases, 11);
    }

    #[test]
    fn typical_align_command() {
        let p = try_parse(&[
            "--runMode",
            "alignReads",
            "--genomeDir",
            "/idx/hg38",
            "--readFilesIn",
            "R1.fq.gz",
            "R2.fq.gz",
            "--readFilesCommand",
            "zcat",
            "--runThreadN",
            "16",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outFileNamePrefix",
            "/out/sample1_",
            "--outFilterMultimapNmax",
            "20",
            "--alignIntronMax",
            "1000000",
            "--sjdbGTFfile",
            "gencode.gtf",
            "--twopassMode",
            "Basic",
        ])
        .unwrap();
        assert_eq!(p.run_mode, RunMode::AlignReads);
        assert_eq!(p.genome_dir, PathBuf::from("/idx/hg38"));
        assert_eq!(
            p.read_files_in,
            vec![PathBuf::from("R1.fq.gz"), PathBuf::from("R2.fq.gz")]
        );
        assert_eq!(p.read_files_command, Some("zcat".to_string()));
        assert_eq!(p.run_thread_n, NonZeroUsize::new(16).unwrap());
        assert_eq!(p.out_sam_type.format, OutSamFormat::Bam);
        assert_eq!(
            p.out_sam_type.sort_order,
            Some(OutSamSortOrder::SortedByCoordinate)
        );
        assert_eq!(p.out_file_name_prefix, "/out/sample1_");
        assert_eq!(p.out_filter_multimap_nmax, 20);
        assert_eq!(p.align_intron_max, 1_000_000);
        assert_eq!(p.sjdb_gtf_file, Some(PathBuf::from("gencode.gtf")));
        assert_eq!(p.twopass_mode, TwopassMode::Basic);
    }

    #[test]
    fn clip_nbases_per_mate() {
        // A single value applies to both mates.
        let p = try_parse(&["--readFilesIn", "r1.fq", "r2.fq", "--clip5pNbases", "7"]).unwrap();
        assert_eq!(p.clip5p(0), 7);
        assert_eq!(p.clip5p(1), 7);
        // Two values are per-mate (mate 1, mate 2).
        let p = try_parse(&[
            "--readFilesIn",
            "r1.fq",
            "r2.fq",
            "--clip5pNbases",
            "39",
            "0",
            "--clip3pNbases",
            "1",
            "2",
        ])
        .unwrap();
        assert_eq!((p.clip5p(0), p.clip5p(1)), (39, 0));
        assert_eq!((p.clip3p(0), p.clip3p(1)), (1, 2));
    }

    #[test]
    fn solo_barcode_mate_validation() {
        let with_mate = |mate: &str| {
            try_parse(&[
                "--readFilesIn",
                "R1.fq",
                "R2.fq",
                "--soloType",
                "CB_UMI_Simple",
                "--soloCBwhitelist",
                "None",
                "--soloCBmatchWLtype",
                "Exact",
                "--sjdbGTFfile",
                "g.gtf",
                "--soloFeatures",
                "Gene",
                "--soloBarcodeMate",
                mate,
            ])
        };
        // Mate 1 (5' paired-end) is accepted; the helper reports it.
        let p = with_mate("1").unwrap();
        assert!(p.solo_barcode_on_mate1());
        // Mate 0 (default) is the standard SE-solo path.
        assert!(!with_mate("0").unwrap().solo_barcode_on_mate1());
        // Mate 2 is rejected with a clear message.
        let err = with_mate("2").unwrap_err().to_string();
        assert!(err.contains("soloBarcodeMate"), "unexpected error: {err}");
    }

    #[test]
    fn scoring_overrides() {
        let p = try_parse(&[
            "--readFilesIn",
            "reads.fq",
            "--scoreGap",
            "0",
            "--scoreGapNoncan",
            "-12",
            "--scoreGapGCAG",
            "-6",
            "--scoreGapATAC",
            "-10",
            "--scoreDelOpen",
            "-3",
            "--scoreDelBase",
            "-1",
            "--scoreInsOpen",
            "-3",
            "--scoreInsBase",
            "-1",
        ])
        .unwrap();
        assert_eq!(p.score_gap, 0);
        assert_eq!(p.score_gap_noncan, -12);
        assert_eq!(p.score_gap_gcag, -6);
        assert_eq!(p.score_gap_atac, -10);
        assert_eq!(p.score_del_open, -3);
        assert_eq!(p.score_del_base, -1);
        assert_eq!(p.score_ins_open, -3);
        assert_eq!(p.score_ins_base, -1);
    }

    #[test]
    fn validate_genome_generate_needs_fasta() {
        let err = try_parse(&["--runMode", "genomeGenerate"]).unwrap_err();
        assert!(err.to_string().contains("genomeFastaFiles"));
    }

    #[test]
    fn validate_align_needs_reads() {
        let err = try_parse(&["--runMode", "alignReads"]).unwrap_err();
        assert!(err.to_string().contains("readFilesIn"));
    }

    #[test]
    fn out_sam_type_parsing() {
        let p = try_parse(&["--readFilesIn", "r.fq", "--outSAMtype", "SAM"]).unwrap();
        assert_eq!(p.out_sam_type.format, OutSamFormat::Sam);
        assert_eq!(p.out_sam_type.sort_order, None);

        let p = try_parse(&["--readFilesIn", "r.fq", "--outSAMtype", "BAM", "Unsorted"]).unwrap();
        assert_eq!(p.out_sam_type.format, OutSamFormat::Bam);
        assert_eq!(p.out_sam_type.sort_order, Some(OutSamSortOrder::Unsorted));

        let p = try_parse(&["--readFilesIn", "r.fq", "--outSAMtype", "None"]).unwrap();
        assert_eq!(p.out_sam_type.format, OutSamFormat::None);

        assert!(try_parse(&["--readFilesIn", "r.fq", "--outSAMtype", "BOGUS"]).is_err());
        assert!(try_parse(&["--readFilesIn", "r.fq", "--outSAMtype", "BAM"]).is_err());
    }

    #[test]
    fn out_sam_unmapped_parsing() {
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert_eq!(p.out_sam_unmapped, OutSamUnmapped::None);

        let p = try_parse(&["--readFilesIn", "r.fq", "--outSAMunmapped", "None"]).unwrap();
        assert_eq!(p.out_sam_unmapped, OutSamUnmapped::None);

        let p = try_parse(&["--readFilesIn", "r.fq", "--outSAMunmapped", "Within"]).unwrap();
        assert_eq!(p.out_sam_unmapped, OutSamUnmapped::Within);

        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--outSAMunmapped",
            "Within",
            "KeepPairs",
        ])
        .unwrap();
        assert_eq!(p.out_sam_unmapped, OutSamUnmapped::WithinKeepPairs);

        assert!(try_parse(&["--readFilesIn", "r.fq", "--outSAMunmapped", "Bogus"]).is_err());
        assert!(
            try_parse(&[
                "--readFilesIn",
                "r.fq",
                "--outSAMunmapped",
                "Within",
                "Bogus"
            ])
            .is_err()
        );
    }

    #[test]
    fn chimeric_params() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "20",
            "--chimScoreMin",
            "10",
            "--chimOutType",
            "WithinBAM",
            "SoftClip",
        ])
        .unwrap();
        assert_eq!(p.chim_segment_min, 20);
        assert_eq!(p.chim_score_min, 10);
        assert_eq!(
            p.chim_out_type,
            vec!["WithinBAM".to_string(), "SoftClip".to_string()]
        );
    }

    #[test]
    fn chimeric_params_extended() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "20",
            "--chimScoreDropMax",
            "30",
            "--chimScoreSeparation",
            "15",
            "--chimMainSegmentMultNmax",
            "5",
            "--chimSegmentReadGapMax",
            "3",
            "--chimJunctionOverhangMin",
            "12",
            "--chimScoreJunctionNonGTAG",
            "-2",
        ])
        .unwrap();
        assert_eq!(p.chim_score_drop_max, 30);
        assert_eq!(p.chim_score_separation, 15);
        assert_eq!(p.chim_main_segment_mult_nmax, 5);
        assert_eq!(p.chim_segment_read_gap_max, 3);
        assert_eq!(p.chim_junction_overhang_min, 12);
        assert_eq!(p.chim_score_junction_non_gtag, -2);
    }

    #[test]
    fn chimeric_params_defaults() {
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert_eq!(p.chim_score_drop_max, 20);
        assert_eq!(p.chim_score_separation, 10);
        assert_eq!(p.chim_main_segment_mult_nmax, 10);
        assert_eq!(p.chim_segment_read_gap_max, 0);
        assert_eq!(p.chim_junction_overhang_min, 20);
        assert_eq!(p.chim_score_junction_non_gtag, -1);
    }

    #[test]
    fn win_bin_window_dist_default() {
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert_eq!(p.win_bin_window_dist(), 589_824); // 2^16 * 9
    }

    #[test]
    fn win_bin_window_dist_custom() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--winBinNbits",
            "14",
            "--winAnchorDistNbins",
            "5",
        ])
        .unwrap();
        assert_eq!(p.win_bin_window_dist(), 81_920); // 2^14 * 5
    }

    #[test]
    fn rg_line_default_unset() {
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert!(!p.rg_line_set());
        assert_eq!(p.parsed_rg_lines().unwrap(), Vec::<String>::new());
        assert_eq!(p.rg_ids().unwrap(), Vec::<String>::new());
        assert!(!p.out_sam_attributes.contains(SamAttributes::RG));
    }

    #[test]
    fn rg_line_single() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--outSAMattrRGline",
            "ID:foo",
            "SM:bar",
            "LB:lib1",
        ])
        .unwrap();
        assert!(p.rg_line_set());
        assert_eq!(
            p.parsed_rg_lines().unwrap(),
            vec!["ID:foo\tSM:bar\tLB:lib1".to_string()]
        );
        assert_eq!(p.rg_ids().unwrap(), vec!["foo".to_string()]);
        assert!(p.out_sam_attributes.contains(SamAttributes::RG));
    }

    #[test]
    fn rg_line_multi() {
        let p = try_parse(&[
            "--readFilesIn",
            "r1.fq",
            "r2.fq",
            "--outSAMattrRGline",
            "ID:a",
            "SM:a",
            ",",
            "ID:b",
            "LB:x",
        ])
        .unwrap();
        let lines = p.parsed_rg_lines().unwrap();
        assert_eq!(
            lines,
            vec!["ID:a\tSM:a".to_string(), "ID:b\tLB:x".to_string()]
        );
        assert_eq!(p.rg_ids().unwrap(), vec!["a".to_string(), "b".to_string()]);
    }

    #[test]
    fn rg_line_single_replicates_for_multi_file() {
        let p = try_parse(&[
            "--readFilesIn",
            "r1.fq",
            "r2.fq",
            "--outSAMattrRGline",
            "ID:foo",
        ])
        .unwrap();
        assert_eq!(
            p.rg_ids().unwrap(),
            vec!["foo".to_string(), "foo".to_string()]
        );
    }

    #[test]
    fn rg_line_missing_id_prefix_errors() {
        let err =
            try_parse(&["--readFilesIn", "r.fq", "--outSAMattrRGline", "SM:oops"]).unwrap_err();
        assert!(err.to_string().contains("ID:"));
    }

    #[test]
    fn rg_line_count_mismatch_errors() {
        // 1 input file, 2 RG entries — mismatch (ids.len()>1 && != n_files).
        let err = try_parse(&[
            "--readFilesIn",
            "r1.fq",
            "--outSAMattrRGline",
            "ID:a",
            ",",
            "ID:b",
        ])
        .unwrap_err();
        assert!(err.to_string().contains("does not match"));
    }

    #[test]
    fn validate_rg_attr_without_line_errors() {
        let err =
            try_parse(&["--readFilesIn", "r.fq", "--outSAMattributes", "NH", "RG"]).unwrap_err();
        assert!(err.to_string().contains("RG"));
    }

    #[test]
    fn run_rng_seed_override() {
        let p = try_parse(&["--readFilesIn", "r.fq", "--runRNGseed", "42"]).unwrap();
        assert_eq!(p.run_rng_seed, 42);
    }

    #[test]
    fn quant_transcriptome_sam_default() {
        use crate::quant::transcriptome::QuantTranscriptomeSAMoutput;
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert!(!p.quant_transcriptome_sam());
        assert_eq!(
            p.quant_transcriptome_sam_output,
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip
        );
    }

    #[test]
    fn quant_transcriptome_sam_enabled() {
        let p = try_parse(&["--readFilesIn", "r.fq", "--quantMode", "TranscriptomeSAM"]).unwrap();
        assert!(p.quant_transcriptome_sam());
        assert!(!p.quant_gene_counts());
    }

    #[test]
    fn quant_transcriptome_sam_output_override() {
        use crate::quant::transcriptome::QuantTranscriptomeSAMoutput;
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--quantTranscriptomeSAMoutput",
            "BanSingleEnd",
        ])
        .unwrap();
        assert_eq!(
            p.quant_transcriptome_sam_output,
            QuantTranscriptomeSAMoutput::BanSingleEnd
        );

        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--quantTranscriptomeSAMoutput",
            "BanSingleEnd_ExtendSoftclip",
        ])
        .unwrap();
        assert_eq!(
            p.quant_transcriptome_sam_output,
            QuantTranscriptomeSAMoutput::BanSingleEndExtendSoftclip
        );
    }

    #[test]
    fn validate_transcriptome_sam_at_genome_generate_needs_gtf() {
        let err = try_parse(&[
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            "g.fa",
            "--quantMode",
            "TranscriptomeSAM",
        ])
        .unwrap_err();
        assert!(err.to_string().contains("TranscriptomeSAM"));
        assert!(err.to_string().contains("sjdbGTFfile"));
    }

    #[test]
    fn validate_transcriptome_sam_at_align_reads_tolerates_no_gtf() {
        // alignReads: if --sjdbGTFfile is absent, the check is deferred to
        // GenomeIndex::load which will either find transcriptInfo.tab in
        // --genomeDir or surface a clear error at load time.
        try_parse(&["--readFilesIn", "r.fq", "--quantMode", "TranscriptomeSAM"]).unwrap();
    }

    #[test]
    fn validate_transcriptome_sam_with_gtf_ok() {
        try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--quantMode",
            "TranscriptomeSAM",
            "--sjdbGTFfile",
            "genes.gtf",
        ])
        .unwrap();
    }

    #[test]
    fn sj_stitch_mismatch() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--alignSJstitchMismatchNmax",
            "1",
            "-1",
            "2",
            "3",
        ])
        .unwrap();
        assert_eq!(p.align_sj_stitch_mismatch_nmax, vec![1, -1, 2, 3]);
    }

    #[test]
    fn output_path_bare_dot_prefix() {
        let p = try_parse(&["--readFilesIn", "r.fq", "--outFileNamePrefix", "SAMPLE."]).unwrap();
        assert_eq!(p.out_file_name_prefix, "SAMPLE.");
        assert_eq!(
            p.output_path("Aligned.out.bam"),
            PathBuf::from("SAMPLE.Aligned.out.bam")
        );
        assert_eq!(
            p.output_path("Log.final.out"),
            PathBuf::from("SAMPLE.Log.final.out")
        );
    }

    #[test]
    fn output_path_trailing_slash_prefix() {
        let p = try_parse(&["--readFilesIn", "r.fq", "--outFileNamePrefix", "out/"]).unwrap();
        assert_eq!(
            p.output_path("Aligned.out.bam"),
            PathBuf::from("out/Aligned.out.bam")
        );
    }

    #[test]
    fn output_path_default_prefix() {
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert_eq!(
            p.output_path("Aligned.out.bam"),
            PathBuf::from("./Aligned.out.bam")
        );
    }

    #[test]
    fn output_path_path_with_underscore() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--outFileNamePrefix",
            "/out/sample1_",
        ])
        .unwrap();
        assert_eq!(
            p.output_path("Aligned.out.bam"),
            PathBuf::from("/out/sample1_Aligned.out.bam")
        );
    }

    #[test]
    fn test_parse_mem_bytes_raw() {
        assert_eq!(parse_mem_bytes("0").unwrap(), 0);
        assert_eq!(parse_mem_bytes("1024").unwrap(), 1024);
        assert_eq!(parse_mem_bytes("31000000000").unwrap(), 31_000_000_000);
    }

    #[test]
    fn test_parse_mem_bytes_suffixes() {
        assert_eq!(parse_mem_bytes("1K").unwrap(), 1024);
        assert_eq!(parse_mem_bytes("1k").unwrap(), 1024);
        assert_eq!(parse_mem_bytes("1M").unwrap(), 1024 * 1024);
        assert_eq!(parse_mem_bytes("1m").unwrap(), 1024 * 1024);
        assert_eq!(parse_mem_bytes("1G").unwrap(), 1024 * 1024 * 1024);
        assert_eq!(parse_mem_bytes("1g").unwrap(), 1024 * 1024 * 1024);
        assert_eq!(parse_mem_bytes("1T").unwrap(), 1024_u64.pow(4));
        assert_eq!(parse_mem_bytes("1t").unwrap(), 1024_u64.pow(4));
        assert_eq!(parse_mem_bytes("64G").unwrap(), 64 * 1024 * 1024 * 1024);
        assert_eq!(parse_mem_bytes("31G").unwrap(), 31 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_parse_mem_bytes_errors() {
        assert!(parse_mem_bytes("abc").is_err());
        assert!(parse_mem_bytes("1X").is_err());
        assert!(parse_mem_bytes("").is_err());
        assert!(parse_mem_bytes("-1G").is_err());
    }

    #[test]
    fn align_ends_type_parses_ext_matrix() {
        use std::str::FromStr;
        assert_eq!(
            AlignEndsType::from_str("Local").unwrap().ext,
            [[false; 2]; 2]
        );
        assert_eq!(
            AlignEndsType::from_str("EndToEnd").unwrap().ext,
            [[true, true], [true, true]]
        );
        assert_eq!(
            AlignEndsType::from_str("Extend5pOfRead1").unwrap().ext,
            [[true, false], [false, false]]
        );
        assert_eq!(
            AlignEndsType::from_str("Extend5pOfReads12").unwrap().ext,
            [[true, false], [true, false]]
        );
        assert_eq!(
            AlignEndsType::from_str("Extend3pOfRead1").unwrap().ext,
            [[false, true], [false, false]]
        );
        assert!(AlignEndsType::from_str("Bogus").is_err());
    }

    #[test]
    fn out_sam_order_accepts_star_values_rejects_others() {
        assert!(try_parse(&["--readFilesIn", "r.fq", "--outSAMorder", "Paired"]).is_ok());
        assert!(
            try_parse(&[
                "--readFilesIn",
                "r.fq",
                "--outSAMorder",
                "PairedKeepInputOrder"
            ])
            .is_ok()
        );
        assert!(try_parse(&["--readFilesIn", "r.fq", "--outSAMorder", "Nope"]).is_err());
    }

    #[test]
    fn align_ends_type_rejected_at_cli() {
        assert!(try_parse(&["--readFilesIn", "r.fq", "--alignEndsType", "EndToEnd"]).is_ok());
        assert!(try_parse(&["--readFilesIn", "r.fq", "--alignEndsType", "Bogus"]).is_err());
    }

    #[test]
    fn xs_strand_field_intron_motif_adds_xs_to_attrs() {
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--outSAMstrandField",
            "intronMotif",
        ])
        .unwrap();
        assert!(p.out_sam_attributes.contains(SamAttributes::XS));
        assert_eq!(p.out_sam_strand_field, "intronMotif");
    }

    #[test]
    fn xs_without_intron_motif_is_stripped() {
        // Explicit XS in attrs without --outSAMstrandField intronMotif: XS gets stripped.
        // Users must set both to get XS output.
        let p = try_parse(&[
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "NH",
            "HI",
            "XS",
        ])
        .unwrap();
        assert_eq!(p.out_sam_strand_field, "None");
        assert!(!p.out_sam_attributes.contains(SamAttributes::XS));
    }

    #[test]
    fn xs_absent_by_default() {
        let p = try_parse(&["--readFilesIn", "r.fq"]).unwrap();
        assert!(!p.out_sam_attributes.contains(SamAttributes::XS));
        assert_eq!(p.out_sam_strand_field, "None");
    }

    #[test]
    fn xs_stripped_from_all_preset_without_intron_motif() {
        // "All" includes XS but it is stripped unless intronMotif is also set.
        let p = try_parse(&["--readFilesIn", "r.fq", "--outSAMattributes", "All"]).unwrap();
        assert_eq!(p.out_sam_strand_field, "None");
        assert!(!p.out_sam_attributes.contains(SamAttributes::XS));
    }

    #[test]
    fn solo_chemistry_presets_set_the_geometry() {
        let base = [
            "--readFilesIn",
            "cdna.fq",
            "bc.fq",
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            "wl.txt",
            "--sjdbGTFfile",
            "genes.gtf",
        ];
        let with = |extra: &[&str]| {
            let mut a = base.to_vec();
            a.extend_from_slice(extra);
            try_parse(&a)
        };

        // v2: 16-base CB, 10-base UMI.
        let p = with(&["--soloChemistry", "SC3Pv2"]).unwrap();
        assert_eq!(
            (
                p.solo_cb_start,
                p.solo_cb_len,
                p.solo_umi_start,
                p.solo_umi_len
            ),
            (1, 16, 17, 10)
        );

        // v3 lengthened the UMI to 12.
        let p = with(&["--soloChemistry", "SC3Pv3"]).unwrap();
        assert_eq!(
            (
                p.solo_cb_start,
                p.solo_cb_len,
                p.solo_umi_start,
                p.solo_umi_len
            ),
            (1, 16, 17, 12)
        );

        // v1 used a 14-base CB.
        let p = with(&["--soloChemistry", "SC3Pv1"]).unwrap();
        assert_eq!((p.solo_cb_len, p.solo_umi_len), (14, 10));

        // The preset wins over an explicit geometry rather than half-applying.
        let p = with(&["--soloChemistry", "SC3Pv3", "--soloCBlen", "14"]).unwrap();
        assert_eq!(p.solo_cb_len, 16);

        // Unknown presets are refused.
        assert!(with(&["--soloChemistry", "SC9Pv9"]).is_err());
    }

    #[test]
    fn solo_cb_type_string_is_refused_rather_than_mispacked() {
        let base = [
            "--readFilesIn",
            "cdna.fq",
            "bc.fq",
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            "wl.txt",
            "--sjdbGTFfile",
            "genes.gtf",
            "--soloCBstart",
            "1",
            "--soloCBlen",
            "16",
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "10",
        ];
        let with = |extra: &[&str]| {
            let mut a = base.to_vec();
            a.extend_from_slice(extra);
            try_parse(&a)
        };
        assert!(with(&["--soloCBtype", "Sequence"]).is_ok());
        // `String` needs a whitelist representation this does not have; packing
        // a non-ACGT barcode into 2 bits per base would produce nonsense.
        assert!(with(&["--soloCBtype", "String"]).is_err());
    }
}
