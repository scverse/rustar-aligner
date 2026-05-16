//! Transcript data structures for storing alignment results
use noodles::sam::alignment::record::cigar;
use std::fmt::Write as _;

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
    pub cigar: Vec<cigar::Op>,
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

pub(crate) trait CigarOpExt {
    fn add_len(&self, len: usize) -> Self;
}
impl CigarOpExt for cigar::Op {
    fn add_len(&self, len: usize) -> Self {
        cigar::Op::new(self.kind(), self.len() + len)
    }
}

pub(crate) trait KindExt {
    fn char(self) -> char;
}
impl KindExt for cigar::op::Kind {
    fn char(self) -> char {
        use cigar::op::Kind;
        match self {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        }
    }
}
impl KindExt for cigar::Op {
    fn char(self) -> char {
        self.kind().char()
    }
}

/// Convert CIGAR operations to CIGAR string
pub(crate) fn cigar_to_string(cigar: &[cigar::Op]) -> String {
    cigar.iter().fold(String::new(), |mut c, op| {
        let _ = write!(c, "{}{}", op.len(), op.char()); // infallible
        c
    })
}

impl Transcript {
    /// Format CIGAR string
    pub fn cigar_string(&self) -> String {
        cigar_to_string(&self.cigar)
    }

    /// Calculate number of matched bases (for filtering)
    pub fn n_matched(&self) -> usize {
        use cigar::op::Kind;
        self.cigar
            .iter()
            .filter_map(|op| match op.kind() {
                Kind::Match | Kind::SequenceMatch => Some(op.len()),
                _ => None,
            })
            .sum()
    }

    /// Calculate read length from CIGAR
    pub fn read_length(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| op.kind().consumes_read())
            .map(|op| op.len())
            .sum()
    }

    /// Calculate reference length from CIGAR
    pub fn reference_length(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum()
    }

    /// Count soft-clipped bases on left and right ends
    ///
    /// Returns `[left_clip, right_clip]` in bases
    pub fn count_soft_clips(&self) -> [usize; 2] {
        [self.cigar.first(), self.cigar.last()].map(|c| {
            c.filter(|c| c.kind() == cigar::op::Kind::SoftClip)
                .map_or(0, |c| c.len())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_to_string() {
        use cigar::op::{Kind, Op};
        let cigar = vec![
            Op::new(Kind::Match, 50),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Deletion, 3),
            Op::new(Kind::Skip, 1000),
            Op::new(Kind::SoftClip, 5),
        ];

        assert_eq!(cigar_to_string(&cigar), "50M2I3D1000N5S");
    }

    #[test]
    fn test_cigar_string() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 250,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 50),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 50),
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
        use cigar::op::Kind;

        assert!(Kind::Match.consumes_read());
        assert!(Kind::Match.consumes_reference());

        assert!(Kind::Insertion.consumes_read());
        assert!(!Kind::Insertion.consumes_reference());

        assert!(!Kind::Deletion.consumes_read());
        assert!(Kind::Deletion.consumes_reference());

        assert!(!Kind::Skip.consumes_read());
        assert!(Kind::Skip.consumes_reference());

        assert!(Kind::SoftClip.consumes_read());
        assert!(!Kind::SoftClip.consumes_reference());
    }

    #[test]
    fn test_transcript_lengths() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 250,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 45),
                Op::new(Kind::Insertion, 3),
                Op::new(Kind::Match, 2),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 50),
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
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 150,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::SoftClip, 10),
                Op::new(Kind::Match, 50),
                Op::new(Kind::SoftClip, 15),
            ],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let [left, right] = transcript.count_soft_clips();
        assert_eq!(left, 10);
        assert_eq!(right, 15);
    }

    #[test]
    fn test_count_soft_clips_left_only() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 150,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::SoftClip, 20), Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let [left, right] = transcript.count_soft_clips();
        assert_eq!(left, 20);
        assert_eq!(right, 0);
    }

    #[test]
    fn test_count_soft_clips_none() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 150,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let [left, right] = transcript.count_soft_clips();
        assert_eq!(left, 0);
        assert_eq!(right, 0);
    }
}
