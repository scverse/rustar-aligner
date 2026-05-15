//! ChimericSegment and ChimericAlignment data structures

use crate::align::transcript::cigar_to_string;
use noodles::sam::alignment::record::cigar;

/// A single segment of a chimeric alignment
#[derive(Debug, Clone)]
pub struct ChimericSegment {
    pub chr_idx: usize,
    pub genome_start: u64,
    pub genome_end: u64,
    pub is_reverse: bool,
    pub read_start: usize,
    pub read_end: usize,
    pub cigar: Vec<cigar::Op>,
    pub score: i32,
    pub n_mismatch: u32,
}

impl ChimericSegment {
    /// Get segment length in read coordinates
    pub fn read_length(&self) -> usize {
        self.read_end - self.read_start
    }

    /// Get segment length in genome coordinates
    pub fn genome_length(&self) -> u64 {
        self.genome_end - self.genome_start
    }

    /// Check if segment meets minimum length requirement
    pub fn meets_min_length(&self, min_len: u32) -> bool {
        self.read_length() >= min_len as usize
    }

    /// Format CIGAR string
    pub fn cigar_string(&self) -> String {
        cigar_to_string(&self.cigar)
    }
}

/// A chimeric alignment consisting of two segments
#[derive(Debug, Clone)]
pub struct ChimericAlignment {
    pub donor: ChimericSegment,
    pub acceptor: ChimericSegment,
    pub junction_type: i32,
    pub repeat_len_donor: u32,
    pub repeat_len_acceptor: u32,
    pub total_score: i32,
    pub read_seq: Vec<u8>,
    pub read_name: String,
}

impl ChimericAlignment {
    /// Create a new chimeric alignment
    pub fn new(
        donor: ChimericSegment,
        acceptor: ChimericSegment,
        junction_type: i32,
        repeat_len_donor: u32,
        repeat_len_acceptor: u32,
        read_seq: Vec<u8>,
        read_name: String,
    ) -> Self {
        let total_score = donor.score + acceptor.score;
        Self {
            donor,
            acceptor,
            junction_type,
            repeat_len_donor,
            repeat_len_acceptor,
            total_score,
            read_seq,
            read_name,
        }
    }

    /// Check if both segments meet minimum length requirement
    pub fn meets_min_segment_length(&self, min_len: u32) -> bool {
        self.donor.meets_min_length(min_len) && self.acceptor.meets_min_length(min_len)
    }

    /// Check if total score meets minimum threshold
    pub fn meets_min_score(&self, min_score: i32) -> bool {
        self.total_score >= min_score
    }

    /// Get the breakpoint position on the donor chromosome (1-based)
    pub fn donor_breakpoint(&self) -> u64 {
        if self.donor.is_reverse {
            self.donor.genome_start + 1
        } else {
            self.donor.genome_end
        }
    }

    /// Get the breakpoint position on the acceptor chromosome (1-based)
    pub fn acceptor_breakpoint(&self) -> u64 {
        if self.acceptor.is_reverse {
            self.acceptor.genome_end
        } else {
            self.acceptor.genome_start + 1
        }
    }

    /// Get strand symbol for donor
    pub fn donor_strand(&self) -> char {
        if self.donor.is_reverse { '-' } else { '+' }
    }

    /// Get strand symbol for acceptor
    pub fn acceptor_strand(&self) -> char {
        if self.acceptor.is_reverse { '-' } else { '+' }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cigar::op::{Kind, Op};

    fn mock_segment(
        chr_idx: usize,
        genome_start: u64,
        genome_end: u64,
        read_start: usize,
        read_end: usize,
    ) -> ChimericSegment {
        ChimericSegment {
            chr_idx,
            genome_start,
            genome_end,
            is_reverse: false,
            read_start,
            read_end,
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 2,
        }
    }

    #[test]
    fn test_segment_lengths() {
        let seg = mock_segment(0, 1000, 1050, 0, 50);
        assert_eq!(seg.read_length(), 50);
        assert_eq!(seg.genome_length(), 50);
    }

    #[test]
    fn test_segment_min_length() {
        let seg = mock_segment(0, 1000, 1050, 0, 50);
        assert!(seg.meets_min_length(20));
        assert!(seg.meets_min_length(50));
        assert!(!seg.meets_min_length(51));
    }

    #[test]
    fn test_chimeric_alignment_creation() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            1,
            0,
            0,
            vec![0; 80],
            "READ_001".to_string(),
        );

        assert_eq!(chim.total_score, 200);
        assert_eq!(chim.junction_type, 1);
        assert_eq!(chim.read_name, "READ_001");
    }

    #[test]
    fn test_chimeric_meets_min_segment_length() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            1,
            0,
            0,
            vec![0; 80],
            "READ_001".to_string(),
        );

        assert!(chim.meets_min_segment_length(20));
        assert!(chim.meets_min_segment_length(30));
        assert!(!chim.meets_min_segment_length(51));
    }

    #[test]
    fn test_chimeric_meets_min_score() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            1,
            0,
            0,
            vec![0; 80],
            "READ_001".to_string(),
        );

        assert!(chim.meets_min_score(100));
        assert!(chim.meets_min_score(200));
        assert!(!chim.meets_min_score(201));
    }

    #[test]
    fn test_breakpoint_positions_forward() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            1,
            0,
            0,
            vec![0; 80],
            "READ_001".to_string(),
        );

        assert_eq!(chim.donor_breakpoint(), 1050); // end position (1-based)
        assert_eq!(chim.acceptor_breakpoint(), 2001); // start position + 1 (1-based)
    }

    #[test]
    fn test_strand_symbols() {
        let mut donor = mock_segment(0, 1000, 1050, 0, 50);
        let mut acceptor = mock_segment(1, 2000, 2030, 50, 80);

        donor.is_reverse = false;
        acceptor.is_reverse = true;

        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            1,
            0,
            0,
            vec![0; 80],
            "READ_001".to_string(),
        );

        assert_eq!(chim.donor_strand(), '+');
        assert_eq!(chim.acceptor_strand(), '-');
    }
}
