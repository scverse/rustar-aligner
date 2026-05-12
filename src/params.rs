use std::collections::HashSet;
use std::path::PathBuf;

use clap::Parser;

// ---------------------------------------------------------------------------
// Run mode enum
// ---------------------------------------------------------------------------

/// STAR's `--runMode` values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunMode {
    AlignReads,
    GenomeGenerate,
}

impl std::str::FromStr for RunMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "alignReads" => Ok(Self::AlignReads),
            "genomeGenerate" => Ok(Self::GenomeGenerate),
            _ => Err(format!(
                "unknown runMode '{s}'; expected 'alignReads' or 'genomeGenerate'"
            )),
        }
    }
}

impl std::fmt::Display for RunMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AlignReads => write!(f, "alignReads"),
            Self::GenomeGenerate => write!(f, "genomeGenerate"),
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
// SAM output type enums
// ---------------------------------------------------------------------------

/// STAR's `--outSAMtype` format component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutSamFormat {
    Sam,
    Bam,
    None,
}

/// STAR's `--outSAMtype` sort order component (only applies to BAM).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutSamSortOrder {
    Unsorted,
    SortedByCoordinate,
}

/// Combined `--outSAMtype` value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutSamType {
    pub format: OutSamFormat,
    pub sort_order: Option<OutSamSortOrder>,
}

impl Default for OutSamType {
    fn default() -> Self {
        Self {
            format: OutSamFormat::Sam,
            sort_order: None,
        }
    }
}

impl std::fmt::Display for OutSamType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (&self.format, &self.sort_order) {
            (OutSamFormat::Sam, _) => write!(f, "SAM"),
            (OutSamFormat::None, _) => write!(f, "None"),
            (OutSamFormat::Bam, Some(OutSamSortOrder::SortedByCoordinate)) => {
                write!(f, "BAM SortedByCoordinate")
            }
            (OutSamFormat::Bam, _) => write!(f, "BAM Unsorted"),
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
// SAM unmapped output
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum OutSamUnmapped {
    #[default]
    None,
    Within,
    WithinKeepPairs,
}

impl std::str::FromStr for OutSamUnmapped {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(Self::None),
            "Within" => Ok(Self::Within),
            "Within KeepPairs" => Ok(Self::WithinKeepPairs),
            _ => Err(format!("unknown outSAMunmapped value: '{s}'")),
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
    #[arg(long = "runThreadN", default_value_t = 1)]
    pub run_thread_n: usize,

    /// Random number generator seed for tie-breaking among equal-scoring alignments
    #[arg(long = "runRNGseed", default_value_t = 777)]
    pub run_rng_seed: u64,

    // ── Genome ──────────────────────────────────────────────────────────
    /// Path to genome index directory
    #[arg(long = "genomeDir", default_value = "./GenomeDir")]
    pub genome_dir: PathBuf,

    /// FASTA file(s) with genome reference sequences (for genomeGenerate)
    #[arg(long = "genomeFastaFiles", num_args = 1..)]
    pub genome_fasta_files: Vec<PathBuf>,

    /// Length of SA pre-indexing string (log2-based)
    #[arg(long = "genomeSAindexNbases", default_value_t = 14)]
    pub genome_sa_index_nbases: u32,

    /// Log2(chromosome bin size) for genome storage
    #[arg(long = "genomeChrBinNbits", default_value_t = 18)]
    pub genome_chr_bin_nbits: u32,

    /// Suffix array sparsity (larger = less RAM, slower mapping)
    #[arg(long = "genomeSAsparseD", default_value_t = 1)]
    pub genome_sa_sparse_d: u32,

    // ── Read files ──────────────────────────────────────────────────────
    /// Input read file(s); second file is mate 2 for paired-end
    #[arg(long = "readFilesIn", num_args = 1..=2)]
    pub read_files_in: Vec<PathBuf>,

    /// Command to decompress input files (e.g. "zcat" for .gz)
    #[arg(long = "readFilesCommand")]
    pub read_files_command: Option<String>,

    /// Number of reads to map; -1 = all
    #[arg(long = "readMapNumber", default_value_t = -1, allow_hyphen_values = true)]
    pub read_map_number: i64,

    /// Bases to clip from 5' end of each mate
    #[arg(long = "clip5pNbases", default_value_t = 0)]
    pub clip5p_nbases: u32,

    /// Bases to clip from 3' end of each mate
    #[arg(long = "clip3pNbases", default_value_t = 0)]
    pub clip3p_nbases: u32,

    // ── Output ──────────────────────────────────────────────────────────
    /// Output file name prefix (including path)
    #[arg(long = "outFileNamePrefix", default_value = "./")]
    pub out_file_name_prefix: PathBuf,

    /// Output type: SAM, BAM Unsorted, BAM SortedByCoordinate, None.
    /// Provide as space-separated tokens, e.g. "BAM SortedByCoordinate".
    #[arg(long = "outSAMtype", num_args = 1..=2, default_values_t = vec!["SAM".to_string()])]
    pub out_sam_type_raw: Vec<String>,

    /// BAM compression level: -1 = uncompressed, 1 (default) to 9 (maximum)
    #[arg(
        long = "outBAMcompression",
        default_value_t = 1,
        allow_hyphen_values = true
    )]
    pub out_bam_compression: i32,

    /// Maximum RAM (bytes) for coordinate-sorted BAM. 0 = unlimited.
    #[arg(long = "limitBAMsortRAM", default_value_t = 0)]
    pub limit_bam_sort_ram: u64,

    /// Route primary alignment output to stdout instead of a file.
    /// Values: None (default), SAM, BAM_Unsorted, BAM_SortedByCoordinate.
    #[arg(long = "outStd", default_value = "None")]
    pub out_std: OutStd,

    /// Strand field: None or intronMotif
    #[arg(long = "outSAMstrandField", default_value = "None")]
    pub out_sam_strand_field: String,

    /// SAM attributes to include (Standard, All, None, or explicit list)
    #[arg(long = "outSAMattributes", num_args = 1.., default_values_t = vec!["Standard".to_string()])]
    pub out_sam_attributes: Vec<String>,

    /// Read group line(s) for the `@RG` SAM header. Space-separated fields;
    /// a bare `,` separates multiple RG blocks. Each block must start with `ID:`.
    /// Default `-` means no `@RG` line (matches STAR).
    #[arg(long = "outSAMattrRGline", num_args = 1.., default_values_t = vec!["-".to_string()])]
    pub out_sam_attr_rg_line: Vec<String>,

    /// Unmapped reads in SAM output: None or Within
    #[arg(long = "outSAMunmapped", default_value = "None")]
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
}

impl Parameters {
    /// Parse the raw `--outSAMtype` tokens into a structured `OutSamType`.
    pub fn out_sam_type(&self) -> Result<OutSamType, String> {
        match self
            .out_sam_type_raw
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .as_slice()
        {
            ["SAM"] => Ok(OutSamType {
                format: OutSamFormat::Sam,
                sort_order: None,
            }),
            ["None"] => Ok(OutSamType {
                format: OutSamFormat::None,
                sort_order: None,
            }),
            ["BAM", "Unsorted"] => Ok(OutSamType {
                format: OutSamFormat::Bam,
                sort_order: Some(OutSamSortOrder::Unsorted),
            }),
            ["BAM", "SortedByCoordinate"] => Ok(OutSamType {
                format: OutSamFormat::Bam,
                sort_order: Some(OutSamSortOrder::SortedByCoordinate),
            }),
            other => Err(format!("unknown outSAMtype: {:?}", other)),
        }
    }

    /// Whether `--chimOutType` includes `Junctions` (write Chimeric.out.junction).
    pub fn chim_out_junctions(&self) -> bool {
        self.chim_out_type.iter().any(|s| s == "Junctions")
    }

    /// Whether `--chimOutType` includes `WithinBAM` (write supplementary BAM records).
    pub fn chim_out_within_bam(&self) -> bool {
        self.chim_out_type.iter().any(|s| s == "WithinBAM")
    }

    /// Expand `--outSAMattributes` into a set of individual tag names.
    ///
    /// - `"Standard"` → {NH, HI, AS, NM, nM}
    /// - `"All"`      → {NH, HI, AS, NM, nM, MD, jM, jI, XS}
    /// - `"None"`     → {} (empty)
    /// - Explicit list (e.g. ["NH", "AS"]) → collected as-is
    ///
    /// `RG` is auto-appended when `--outSAMattrRGline` is set (STAR behavior,
    /// `Parameters_samAttributes.cpp:201`).
    pub fn sam_attribute_set(&self) -> HashSet<String> {
        let mut attrs: HashSet<String> = match self
            .out_sam_attributes
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .as_slice()
        {
            ["Standard"] => ["NH", "HI", "AS", "NM", "nM"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            ["All"] => ["NH", "HI", "AS", "NM", "nM", "MD", "jM", "jI", "XS"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            ["None"] => HashSet::new(),
            tags => tags.iter().map(|s| s.to_string()).collect(),
        };
        if self.rg_line_set() {
            attrs.insert("RG".to_string());
        }
        attrs
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
                        "--outSAMattrRGline: first field of each RG line must start with 'ID:', got '{}'",
                        first
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
                "alignIntronMax=alignMatesGapMax=0, max intron ~= (2^winBinNbits)*winAnchorDistNbins={}",
                max_intron
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

    /// Validate parameter combinations that clap alone cannot enforce.
    pub fn validate(&self) -> Result<(), crate::error::Error> {
        // genomeGenerate requires FASTA files
        if self.run_mode == RunMode::GenomeGenerate && self.genome_fasta_files.is_empty() {
            return Err(crate::error::Error::Parameter(
                "--genomeFastaFiles is required when --runMode genomeGenerate".into(),
            ));
        }

        // alignReads requires read files
        if self.run_mode == RunMode::AlignReads && self.read_files_in.is_empty() {
            return Err(crate::error::Error::Parameter(
                "--readFilesIn is required when --runMode alignReads".into(),
            ));
        }

        // Validate outSAMtype
        self.out_sam_type()
            .map_err(crate::error::Error::Parameter)?;

        // Thread count must be at least 1
        if self.run_thread_n == 0 {
            return Err(crate::error::Error::Parameter(
                "--runThreadN must be >= 1".into(),
            ));
        }

        // quantMode GeneCounts requires a GTF file
        if self.quant_gene_counts() && self.sjdb_gtf_file.is_none() {
            return Err(crate::error::Error::Parameter(
                "--quantMode GeneCounts requires --sjdbGTFfile".into(),
            ));
        }

        // Read group: `RG` in outSAMattributes without an RG line is a fatal
        // error (STAR: Parameters_samAttributes.cpp:206). STAR's "All" preset
        // does NOT include RG, so only match a literal "RG" token here. Also
        // parse the RG line to validate its ID: prefix and per-file RG count.
        let user_wants_rg_attr = self.out_sam_attributes.iter().any(|a| a == "RG");
        if !self.rg_line_set() && user_wants_rg_attr {
            return Err(crate::error::Error::Parameter(
                "--outSAMattributes contains RG tag, but --outSAMattrRGline is not set".into(),
            ));
        }
        self.rg_ids()?;

        // quantMode TranscriptomeSAM requires transcript annotations —
        // either via --sjdbGTFfile or pre-generated transcriptInfo.tab
        // et al in --genomeDir (persisted at genomeGenerate time). At
        // validation time we can only enforce the genomeGenerate rule;
        // for alignReads, GenomeIndex::load checks for the on-disk files
        // and surfaces a clear error if neither source is available.
        if self.run_mode == RunMode::GenomeGenerate
            && self.quant_transcriptome_sam()
            && self.sjdb_gtf_file.is_none()
        {
            return Err(crate::error::Error::Parameter(
                "--quantMode TranscriptomeSAM requires --sjdbGTFfile at genomeGenerate".into(),
            ));
        }

        Ok(())
    }

    /// Returns true if `--quantMode GeneCounts` was requested.
    pub fn quant_gene_counts(&self) -> bool {
        self.quant_mode.iter().any(|m| m == "GeneCounts")
    }

    /// Returns true if `--quantMode TranscriptomeSAM` was requested.
    pub fn quant_transcriptome_sam(&self) -> bool {
        self.quant_mode.iter().any(|m| m == "TranscriptomeSAM")
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: parse a STAR-style command line (without program name).
    fn parse(args: &[&str]) -> Parameters {
        let mut full = vec!["rustar-aligner"];
        full.extend_from_slice(args);
        Parameters::parse_from(full)
    }

    #[test]
    fn defaults() {
        let p = parse(&["--readFilesIn", "reads.fq"]);
        assert_eq!(p.run_mode, RunMode::AlignReads);
        assert_eq!(p.run_thread_n, 1);
        assert_eq!(p.run_rng_seed, 777);
        assert_eq!(p.genome_dir, PathBuf::from("./GenomeDir"));
        assert_eq!(p.genome_sa_index_nbases, 14);
        assert_eq!(p.genome_chr_bin_nbits, 18);
        assert_eq!(p.genome_sa_sparse_d, 1);
        assert_eq!(p.read_map_number, -1);
        assert_eq!(p.clip5p_nbases, 0);
        assert_eq!(p.clip3p_nbases, 0);
        assert_eq!(p.out_file_name_prefix, PathBuf::from("./"));
        assert_eq!(p.out_sam_type_raw, vec!["SAM".to_string()]);
        assert_eq!(p.out_sam_strand_field, "None");
        assert_eq!(p.out_sam_attributes, vec!["Standard".to_string()]);
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
            vec![50000, 100000, 200000]
        );
    }

    #[test]
    fn genome_generate_mode() {
        let p = parse(&[
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            "/data/genome",
            "--genomeFastaFiles",
            "chr1.fa",
            "chr2.fa",
            "--runThreadN",
            "8",
            "--genomeSAindexNbases",
            "11",
        ]);
        assert_eq!(p.run_mode, RunMode::GenomeGenerate);
        assert_eq!(p.genome_dir, PathBuf::from("/data/genome"));
        assert_eq!(
            p.genome_fasta_files,
            vec![PathBuf::from("chr1.fa"), PathBuf::from("chr2.fa")]
        );
        assert_eq!(p.run_thread_n, 8);
        assert_eq!(p.genome_sa_index_nbases, 11);
    }

    #[test]
    fn typical_align_command() {
        let p = parse(&[
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
        ]);
        assert_eq!(p.run_mode, RunMode::AlignReads);
        assert_eq!(p.genome_dir, PathBuf::from("/idx/hg38"));
        assert_eq!(
            p.read_files_in,
            vec![PathBuf::from("R1.fq.gz"), PathBuf::from("R2.fq.gz")]
        );
        assert_eq!(p.read_files_command, Some("zcat".to_string()));
        assert_eq!(p.run_thread_n, 16);
        assert_eq!(
            p.out_sam_type_raw,
            vec!["BAM".to_string(), "SortedByCoordinate".to_string()]
        );
        let sam_type = p.out_sam_type().unwrap();
        assert_eq!(sam_type.format, OutSamFormat::Bam);
        assert_eq!(
            sam_type.sort_order,
            Some(OutSamSortOrder::SortedByCoordinate)
        );
        assert_eq!(p.out_file_name_prefix, PathBuf::from("/out/sample1_"));
        assert_eq!(p.out_filter_multimap_nmax, 20);
        assert_eq!(p.align_intron_max, 1_000_000);
        assert_eq!(p.sjdb_gtf_file, Some(PathBuf::from("gencode.gtf")));
        assert_eq!(p.twopass_mode, TwopassMode::Basic);
    }

    #[test]
    fn scoring_overrides() {
        let p = parse(&[
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
        ]);
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
        let p = parse(&["--runMode", "genomeGenerate"]);
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("genomeFastaFiles"));
    }

    #[test]
    fn validate_align_needs_reads() {
        let p = parse(&["--runMode", "alignReads"]);
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("readFilesIn"));
    }

    #[test]
    fn out_sam_type_parsing() {
        let p = parse(&["--readFilesIn", "r.fq", "--outSAMtype", "SAM"]);
        let t = p.out_sam_type().unwrap();
        assert_eq!(t.format, OutSamFormat::Sam);
        assert_eq!(t.sort_order, None);

        let p = parse(&["--readFilesIn", "r.fq", "--outSAMtype", "BAM", "Unsorted"]);
        let t = p.out_sam_type().unwrap();
        assert_eq!(t.format, OutSamFormat::Bam);
        assert_eq!(t.sort_order, Some(OutSamSortOrder::Unsorted));

        let p = parse(&["--readFilesIn", "r.fq", "--outSAMtype", "None"]);
        let t = p.out_sam_type().unwrap();
        assert_eq!(t.format, OutSamFormat::None);
    }

    #[test]
    fn chimeric_params() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "20",
            "--chimScoreMin",
            "10",
            "--chimOutType",
            "WithinBAM",
            "SoftClip",
        ]);
        assert_eq!(p.chim_segment_min, 20);
        assert_eq!(p.chim_score_min, 10);
        assert_eq!(
            p.chim_out_type,
            vec!["WithinBAM".to_string(), "SoftClip".to_string()]
        );
    }

    #[test]
    fn chimeric_params_extended() {
        let p = parse(&[
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
        ]);
        assert_eq!(p.chim_score_drop_max, 30);
        assert_eq!(p.chim_score_separation, 15);
        assert_eq!(p.chim_main_segment_mult_nmax, 5);
        assert_eq!(p.chim_segment_read_gap_max, 3);
        assert_eq!(p.chim_junction_overhang_min, 12);
        assert_eq!(p.chim_score_junction_non_gtag, -2);
    }

    #[test]
    fn chimeric_params_defaults() {
        let p = parse(&["--readFilesIn", "r.fq"]);
        assert_eq!(p.chim_score_drop_max, 20);
        assert_eq!(p.chim_score_separation, 10);
        assert_eq!(p.chim_main_segment_mult_nmax, 10);
        assert_eq!(p.chim_segment_read_gap_max, 0);
        assert_eq!(p.chim_junction_overhang_min, 20);
        assert_eq!(p.chim_score_junction_non_gtag, -1);
    }

    #[test]
    fn win_bin_window_dist_default() {
        let p = parse(&["--readFilesIn", "r.fq"]);
        assert_eq!(p.win_bin_window_dist(), 589_824); // 2^16 * 9
    }

    #[test]
    fn win_bin_window_dist_custom() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--winBinNbits",
            "14",
            "--winAnchorDistNbins",
            "5",
        ]);
        assert_eq!(p.win_bin_window_dist(), 81_920); // 2^14 * 5
    }

    #[test]
    fn rg_line_default_unset() {
        let p = parse(&["--readFilesIn", "r.fq"]);
        assert!(!p.rg_line_set());
        assert_eq!(p.parsed_rg_lines().unwrap(), Vec::<String>::new());
        assert_eq!(p.rg_ids().unwrap(), Vec::<String>::new());
        assert!(!p.sam_attribute_set().contains("RG"));
    }

    #[test]
    fn rg_line_single() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--outSAMattrRGline",
            "ID:foo",
            "SM:bar",
            "LB:lib1",
        ]);
        assert!(p.rg_line_set());
        assert_eq!(
            p.parsed_rg_lines().unwrap(),
            vec!["ID:foo\tSM:bar\tLB:lib1".to_string()]
        );
        assert_eq!(p.rg_ids().unwrap(), vec!["foo".to_string()]);
        assert!(p.sam_attribute_set().contains("RG"));
    }

    #[test]
    fn rg_line_multi() {
        let p = parse(&[
            "--readFilesIn",
            "r1.fq",
            "r2.fq",
            "--outSAMattrRGline",
            "ID:a",
            "SM:a",
            ",",
            "ID:b",
            "LB:x",
        ]);
        let lines = p.parsed_rg_lines().unwrap();
        assert_eq!(
            lines,
            vec!["ID:a\tSM:a".to_string(), "ID:b\tLB:x".to_string()]
        );
        assert_eq!(p.rg_ids().unwrap(), vec!["a".to_string(), "b".to_string()]);
    }

    #[test]
    fn rg_line_single_replicates_for_multi_file() {
        let p = parse(&[
            "--readFilesIn",
            "r1.fq",
            "r2.fq",
            "--outSAMattrRGline",
            "ID:foo",
        ]);
        assert_eq!(
            p.rg_ids().unwrap(),
            vec!["foo".to_string(), "foo".to_string()]
        );
    }

    #[test]
    fn rg_line_missing_id_prefix_errors() {
        let p = parse(&["--readFilesIn", "r.fq", "--outSAMattrRGline", "SM:oops"]);
        let err = p.parsed_rg_lines().unwrap_err();
        assert!(err.to_string().contains("ID:"));
    }

    #[test]
    fn rg_line_count_mismatch_errors() {
        // 1 input file, 2 RG entries — mismatch (ids.len()>1 && != n_files).
        let p = parse(&[
            "--readFilesIn",
            "r1.fq",
            "--outSAMattrRGline",
            "ID:a",
            ",",
            "ID:b",
        ]);
        let err = p.rg_ids().unwrap_err();
        assert!(err.to_string().contains("does not match"));
    }

    #[test]
    fn validate_rg_attr_without_line_errors() {
        let p = parse(&["--readFilesIn", "r.fq", "--outSAMattributes", "NH", "RG"]);
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("RG"));
    }

    fn run_rng_seed_override() {
        let p = parse(&["--readFilesIn", "r.fq", "--runRNGseed", "42"]);
        assert_eq!(p.run_rng_seed, 42);
    }

    #[test]
    fn quant_transcriptome_sam_default() {
        use crate::quant::transcriptome::QuantTranscriptomeSAMoutput;
        let p = parse(&["--readFilesIn", "r.fq"]);
        assert!(!p.quant_transcriptome_sam());
        assert_eq!(
            p.quant_transcriptome_sam_output,
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip
        );
    }

    #[test]
    fn quant_transcriptome_sam_enabled() {
        let p = parse(&["--readFilesIn", "r.fq", "--quantMode", "TranscriptomeSAM"]);
        assert!(p.quant_transcriptome_sam());
        assert!(!p.quant_gene_counts());
    }

    #[test]
    fn quant_transcriptome_sam_output_override() {
        use crate::quant::transcriptome::QuantTranscriptomeSAMoutput;
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--quantTranscriptomeSAMoutput",
            "BanSingleEnd",
        ]);
        assert_eq!(
            p.quant_transcriptome_sam_output,
            QuantTranscriptomeSAMoutput::BanSingleEnd
        );

        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--quantTranscriptomeSAMoutput",
            "BanSingleEnd_ExtendSoftclip",
        ]);
        assert_eq!(
            p.quant_transcriptome_sam_output,
            QuantTranscriptomeSAMoutput::BanSingleEndExtendSoftclip
        );
    }

    #[test]
    fn validate_transcriptome_sam_at_genome_generate_needs_gtf() {
        let p = parse(&[
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            "g.fa",
            "--quantMode",
            "TranscriptomeSAM",
        ]);
        let err = p.validate().unwrap_err();
        assert!(err.to_string().contains("TranscriptomeSAM"));
        assert!(err.to_string().contains("sjdbGTFfile"));
    }

    #[test]
    fn validate_transcriptome_sam_at_align_reads_tolerates_no_gtf() {
        // alignReads: if --sjdbGTFfile is absent, the check is deferred to
        // GenomeIndex::load which will either find transcriptInfo.tab in
        // --genomeDir or surface a clear error at load time.
        let p = parse(&["--readFilesIn", "r.fq", "--quantMode", "TranscriptomeSAM"]);
        assert!(p.validate().is_ok());
    }

    #[test]
    fn validate_transcriptome_sam_with_gtf_ok() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--quantMode",
            "TranscriptomeSAM",
            "--sjdbGTFfile",
            "genes.gtf",
        ]);
        assert!(p.validate().is_ok());
    }

    #[test]
    fn sj_stitch_mismatch() {
        let p = parse(&[
            "--readFilesIn",
            "r.fq",
            "--alignSJstitchMismatchNmax",
            "1",
            "-1",
            "2",
            "3",
        ]);
        assert_eq!(p.align_sj_stitch_mismatch_nmax, vec![1, -1, 2, 3]);
    }
}
