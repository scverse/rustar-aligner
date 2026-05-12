/// Transcript data structures for storing alignment results
use std::fmt;

/// A complete alignment of a read to the genome
#[derive(Debug, Clone)]
pub struct Transcript {
    /// Chromosome index
    pub chr_idx: usize,
    /// Leftmost genomic position (0-based)
    pub genome_start: u64,
    /// Rightmost genomic position (exclusive)
    pub genome_end: u64,
    /// Strand (false = forward, true = reverse)
    pub is_reverse: bool,
    /// Exon segments
    pub exons: Vec<Exon>,
    /// CIGAR operations
    pub cigar: Vec<CigarOp>,
    /// Alignment score
    pub score: i32,
    /// Number of mismatches
    pub n_mismatch: u32,
    /// Number of indels (not splice junctions)
    pub n_gap: u32,
    /// Number of splice junctions
    pub n_junction: u32,
    /// Splice junction motifs (one per junction in CIGAR)
    pub junction_motifs: Vec<crate::align::score::SpliceMotif>,
    /// Whether each junction is annotated in the GTF (for jM +20 offset)
    pub junction_annotated: Vec<bool>,
    /// Original read sequence
    pub read_seq: Vec<u8>,
}

/// An exon segment in a transcript.
///
/// `i_frag` carries the mate index (0 for mate1 / single-end, 1 for mate2),
/// matching STAR's `Transcript::exons[i][EX_iFrag]`
/// (`source/IncludeDefine.h:209`). PE alignments span both mates; this
/// field marks which mate each exon block belongs to and is the input for
/// STAR's single-end filter in transcriptome projection
/// (`ReadAlign_quantTranscriptome.cpp:17`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Exon {
    /// Genomic start position (0-based, inclusive)
    pub genome_start: u64,
    /// Genomic end position (0-based, exclusive)
    pub genome_end: u64,
    /// Read start position (0-based, inclusive)
    pub read_start: usize,
    /// Read end position (0-based, exclusive)
    pub read_end: usize,
    /// Mate index (0 for mate1/SE, 1 for mate2). Matches STAR's `EX_iFrag`.
    pub i_frag: u8,
}

/// CIGAR operation
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CigarOp {
    /// M: match/mismatch (default mode)
    Match(u32),
    /// =: exact match (optional)
    Equal(u32),
    /// X: mismatch (optional)
    Diff(u32),
    /// I: insertion to reference
    Ins(u32),
    /// D: deletion from reference
    Del(u32),
    /// N: splice junction (skipped reference region)
    RefSkip(u32),
    /// S: soft clip (clipped sequence present in read)
    SoftClip(u32),
    /// H: hard clip (clipped sequence not present)
    HardClip(u32),
}

impl CigarOp {
    /// Get the operation character
    pub fn op_char(&self) -> char {
        match self {
            CigarOp::Match(_) => 'M',
            CigarOp::Equal(_) => '=',
            CigarOp::Diff(_) => 'X',
            CigarOp::Ins(_) => 'I',
            CigarOp::Del(_) => 'D',
            CigarOp::RefSkip(_) => 'N',
            CigarOp::SoftClip(_) => 'S',
            CigarOp::HardClip(_) => 'H',
        }
    }

    /// Get the operation length
    pub fn len(&self) -> u32 {
        match self {
            CigarOp::Match(n)
            | CigarOp::Equal(n)
            | CigarOp::Diff(n)
            | CigarOp::Ins(n)
            | CigarOp::Del(n)
            | CigarOp::RefSkip(n)
            | CigarOp::SoftClip(n)
            | CigarOp::HardClip(n) => *n,
        }
    }

    /// Check if operation is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Check if operation consumes query bases
    pub fn consumes_query(&self) -> bool {
        matches!(
            self,
            CigarOp::Match(_)
                | CigarOp::Equal(_)
                | CigarOp::Diff(_)
                | CigarOp::Ins(_)
                | CigarOp::SoftClip(_)
        )
    }

    /// Check if operation consumes reference bases
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            CigarOp::Match(_)
                | CigarOp::Equal(_)
                | CigarOp::Diff(_)
                | CigarOp::Del(_)
                | CigarOp::RefSkip(_)
        )
    }
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.op_char())
    }
}

impl Transcript {
    /// Format CIGAR string
    pub fn cigar_string(&self) -> String {
        self.cigar.iter().map(|op| op.to_string()).collect()
    }

    /// Calculate number of matched bases (for filtering)
    pub fn n_matched(&self) -> u32 {
        self.cigar
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(n) | CigarOp::Equal(n) => Some(*n),
                _ => None,
            })
            .sum()
    }

    /// Calculate read length from CIGAR
    pub fn read_length(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|op| op.consumes_query())
            .map(|op| op.len())
            .sum()
    }

    /// Calculate reference length from CIGAR
    pub fn reference_length(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op.len())
            .sum()
    }

    /// Count soft-clipped bases on left and right ends
    ///
    /// Returns (left_clip, right_clip) in bases
    pub fn count_soft_clips(&self) -> (u32, u32) {
        let mut left_clip = 0u32;
        let mut right_clip = 0u32;

        // Check first operation for left clip
        if let Some(CigarOp::SoftClip(n)) = self.cigar.first() {
            left_clip = *n;
        }

        // Check last operation for right clip
        if let Some(CigarOp::SoftClip(n)) = self.cigar.last() {
            right_clip = *n;
        }

        (left_clip, right_clip)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_op_display() {
        assert_eq!(CigarOp::Match(50).to_string(), "50M");
        assert_eq!(CigarOp::Ins(3).to_string(), "3I");
        assert_eq!(CigarOp::Del(2).to_string(), "2D");
        assert_eq!(CigarOp::RefSkip(1000).to_string(), "1000N");
        assert_eq!(CigarOp::SoftClip(5).to_string(), "5S");
    }

    #[test]
    fn test_cigar_string() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 250,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(100),
                CigarOp::Match(50),
            ],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        assert_eq!(transcript.cigar_string(), "50M100N50M");
    }

    #[test]
    fn test_cigar_consumes() {
        assert!(CigarOp::Match(10).consumes_query());
        assert!(CigarOp::Match(10).consumes_reference());

        assert!(CigarOp::Ins(5).consumes_query());
        assert!(!CigarOp::Ins(5).consumes_reference());

        assert!(!CigarOp::Del(3).consumes_query());
        assert!(CigarOp::Del(3).consumes_reference());

        assert!(!CigarOp::RefSkip(1000).consumes_query());
        assert!(CigarOp::RefSkip(1000).consumes_reference());

        assert!(CigarOp::SoftClip(5).consumes_query());
        assert!(!CigarOp::SoftClip(5).consumes_reference());
    }

    #[test]
    fn test_transcript_lengths() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 250,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::Match(45),
                CigarOp::Ins(3),
                CigarOp::Match(2),
                CigarOp::RefSkip(100),
                CigarOp::Match(50),
            ],
            score: 100,
            n_mismatch: 0,
            n_gap: 1,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0; 100],
        };

        // Query: 45 + 3 + 2 + 50 = 100
        assert_eq!(transcript.read_length(), 100);

        // Reference: 45 + 2 + 100 + 50 = 197 (no insertion)
        assert_eq!(transcript.reference_length(), 197);

        // Matched bases: 45 + 2 + 50 = 97
        assert_eq!(transcript.n_matched(), 97);
    }

    #[test]
    fn test_exon() {
        let exon = Exon {
            genome_start: 1000,
            genome_end: 1050,
            read_start: 0,
            read_end: 50,
            i_frag: 0,
        };

        assert_eq!(exon.genome_start, 1000);
        assert_eq!(exon.genome_end, 1050);
        assert_eq!(exon.read_start, 0);
        assert_eq!(exon.read_end, 50);
    }

    #[test]
    fn test_count_soft_clips_both_ends() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 150,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                CigarOp::SoftClip(10),
                CigarOp::Match(50),
                CigarOp::SoftClip(15),
            ],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let (left, right) = transcript.count_soft_clips();
        assert_eq!(left, 10);
        assert_eq!(right, 15);
    }

    #[test]
    fn test_count_soft_clips_left_only() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 150,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::SoftClip(20), CigarOp::Match(50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let (left, right) = transcript.count_soft_clips();
        assert_eq!(left, 20);
        assert_eq!(right, 0);
    }

    #[test]
    fn test_count_soft_clips_none() {
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 150,
            is_reverse: false,
            exons: vec![],
            cigar: vec![CigarOp::Match(50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let (left, right) = transcript.count_soft_clips();
        assert_eq!(left, 0);
        assert_eq!(right, 0);
    }
}
