// Chimeric junction scoring and classification

use crate::genome::Genome;

/// Classify junction type based on donor/acceptor splice motifs
///
/// Returns junction type encoding:
/// - 0 = non-canonical
/// - 1 = GT/AG (canonical, + strand)
/// - 2 = CT/AC (reverse of GT/AG, - strand)
/// - 3 = GC/AG
/// - 4 = CT/GC (reverse of GC/AG)
/// - 5 = AT/AC
/// - 6 = GT/AT (rare)
pub fn classify_junction_type(
    genome: &Genome,
    donor_chr: usize,
    donor_pos: u64,
    donor_strand: bool,
    acceptor_chr: usize,
    acceptor_pos: u64,
    acceptor_strand: bool,
) -> i32 {
    // For inter-chromosomal or different strand breaks, always non-canonical
    if donor_chr != acceptor_chr || donor_strand != acceptor_strand {
        return 0;
    }

    // Extract 2 bases at donor junction and 2 bases at acceptor junction
    let donor_motif = extract_motif(genome, donor_chr, donor_pos, donor_strand, true);
    let acceptor_motif = extract_motif(genome, acceptor_chr, acceptor_pos, acceptor_strand, false);

    match (donor_motif.as_str(), acceptor_motif.as_str()) {
        ("GT", "AG") => 1,
        ("CT", "AC") => 2,
        ("GC", "AG") => 3,
        ("CT", "GC") => 4,
        ("AT", "AC") => 5,
        ("GT", "AT") => 6,
        _ => 0,
    }
}

/// Extract 2-base motif at junction site
///
/// For donor (is_donor=true): extract 2 bases after the junction
/// For acceptor (is_donor=false): extract 2 bases before the junction
fn extract_motif(
    genome: &Genome,
    chr_idx: usize,
    pos: u64,
    is_reverse: bool,
    is_donor: bool,
) -> String {
    let chr_start = genome.chr_start[chr_idx];
    let chr_len = genome.chr_length[chr_idx];

    // Calculate extraction position
    let extract_pos = if is_donor {
        pos // donor: bases after junction
    } else {
        if pos < 2 {
            return "NN".to_string();
        }
        pos - 2 // acceptor: 2 bases before junction
    };

    // Bounds check
    if extract_pos + 2 > chr_len {
        return "NN".to_string();
    }

    let genome_idx = (chr_start + extract_pos) as usize;
    let b1 = genome.sequence.get(genome_idx).copied().unwrap_or(4);
    let b2 = genome.sequence.get(genome_idx + 1).copied().unwrap_or(4);

    // Convert to bases
    let mut motif = vec![base_to_char(b1), base_to_char(b2)];

    // Reverse complement if on reverse strand
    if is_reverse {
        motif.reverse();
        motif = motif.iter().map(|&c| complement(c)).collect();
    }

    motif.into_iter().collect()
}

/// Calculate repeat length at junction
///
/// STAR definition: Number of bases at junction that are identical
/// on both sides (donor and acceptor)
pub fn calculate_repeat_length(
    genome: &Genome,
    donor_chr: usize,
    donor_pos: u64,
    acceptor_chr: usize,
    acceptor_pos: u64,
    max_check: usize,
) -> (u32, u32) {
    // If different chromosomes, no repeat
    if donor_chr != acceptor_chr {
        return (0, 0);
    }

    let chr_start = genome.chr_start[donor_chr];
    let chr_len = genome.chr_length[donor_chr];

    // Calculate how far we can check
    let donor_idx = chr_start + donor_pos;
    //let acceptor_idx = chr_start + acceptor_pos;

    let mut repeat_len_donor = 0u32;
    let mut repeat_len_acceptor = 0u32;

    // Check forward from donor and backward from acceptor
    for i in 0..max_check {
        let d_pos = donor_idx + i as u64;
        let a_pos = if acceptor_pos < i as u64 {
            break;
        } else {
            chr_start + (acceptor_pos - i as u64 - 1)
        };

        // Bounds check
        if d_pos >= chr_start + chr_len || a_pos < chr_start {
            break;
        }

        let d_base = genome.sequence.get(d_pos as usize).copied().unwrap_or(4);
        let a_base = genome.sequence.get(a_pos as usize).copied().unwrap_or(4);

        if d_base == a_base && d_base < 4 {
            // Only count ACGT, not N
            repeat_len_donor += 1;
            repeat_len_acceptor += 1;
        } else {
            break;
        }
    }

    (repeat_len_donor, repeat_len_acceptor)
}

/// Convert base encoding to character
fn base_to_char(base: u8) -> char {
    match base {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

/// Get complement base
fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => 'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;

    fn mock_genome_with_sequence(seq: Vec<u8>) -> Genome {
        Genome {
            sequence: seq,
            n_genome: 100,
            n_genome_real: 100,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_start: vec![0],
            chr_length: vec![100],
        }
    }

    #[test]
    fn test_base_to_char() {
        assert_eq!(base_to_char(0), 'A');
        assert_eq!(base_to_char(1), 'C');
        assert_eq!(base_to_char(2), 'G');
        assert_eq!(base_to_char(3), 'T');
        assert_eq!(base_to_char(4), 'N');
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement('A'), 'T');
        assert_eq!(complement('T'), 'A');
        assert_eq!(complement('C'), 'G');
        assert_eq!(complement('G'), 'C');
        assert_eq!(complement('N'), 'N');
    }

    #[test]
    fn test_extract_motif_donor_forward() {
        // Sequence: ...ACGTAG...
        //                 ^^ donor at pos 4, extract GT
        let seq = vec![0, 1, 2, 3, 0, 2]; // ACGTAG
        let genome = mock_genome_with_sequence(seq);

        let motif = extract_motif(&genome, 0, 4, false, true);
        assert_eq!(motif, "AG");
    }

    #[test]
    fn test_extract_motif_acceptor_forward() {
        // Sequence: ...GTACAG...
        //              ^^ acceptor at pos 4, extract GT (2 bases before)
        let seq = vec![2, 3, 0, 1, 0, 2]; // GTACAG
        let genome = mock_genome_with_sequence(seq);

        let motif = extract_motif(&genome, 0, 4, false, false);
        assert_eq!(motif, "AC");
    }

    #[test]
    fn test_classify_junction_type_canonical() {
        // GT...AG canonical junction
        let seq = vec![2, 3, 0, 0, 0, 0, 2]; // GT....AG
        let genome = mock_genome_with_sequence(seq);

        let jtype = classify_junction_type(&genome, 0, 0, false, 0, 7, false);
        assert_eq!(jtype, 1); // GT/AG
    }

    #[test]
    fn test_classify_junction_type_inter_chromosomal() {
        let genome = mock_genome_with_sequence(vec![0; 10]);
        let jtype = classify_junction_type(&genome, 0, 5, false, 1, 10, false);
        assert_eq!(jtype, 0); // Inter-chromosomal always non-canonical
    }

    #[test]
    fn test_classify_junction_type_strand_break() {
        let genome = mock_genome_with_sequence(vec![0; 10]);
        let jtype = classify_junction_type(&genome, 0, 5, false, 0, 10, true);
        assert_eq!(jtype, 0); // Strand break always non-canonical
    }

    #[test]
    fn test_calculate_repeat_length_no_repeat() {
        // No repeating bases at junction
        let seq = vec![0, 1, 2, 3, 0, 1]; // ACGTAC
        let genome = mock_genome_with_sequence(seq);

        let (rep_donor, rep_acceptor) = calculate_repeat_length(&genome, 0, 2, 0, 4, 10);
        assert_eq!(rep_donor, 0);
        assert_eq!(rep_acceptor, 0);
    }

    #[test]
    fn test_calculate_repeat_length_with_repeat() {
        // AAA at junction: ...CAAA|AAAG...
        // Positions:       0   1234 5678
        let seq = vec![1, 0, 0, 0, 0, 0, 0, 2]; // CAAAAAAG
        let genome = mock_genome_with_sequence(seq);

        // Donor breakpoint at 4, acceptor at 7
        // Repeating region: positions 4,5,6 (3 As)
        let (rep_donor, rep_acceptor) = calculate_repeat_length(&genome, 0, 4, 0, 7, 10);
        assert_eq!(rep_donor, 3);
        assert_eq!(rep_acceptor, 3);
    }

    #[test]
    fn test_calculate_repeat_length_inter_chromosomal() {
        let genome = mock_genome_with_sequence(vec![0; 10]);
        let (rep_donor, rep_acceptor) = calculate_repeat_length(&genome, 0, 5, 1, 10, 10);
        assert_eq!(rep_donor, 0);
        assert_eq!(rep_acceptor, 0);
    }
}
