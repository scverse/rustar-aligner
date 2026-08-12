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
    let b1 = genome.sequence.get(genome_idx).unwrap_or(4);
    let b2 = genome.sequence.get(genome_idx + 1).unwrap_or(4);

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
/// Repeat lengths either side of a chimeric junction, as STAR computes them.
///
/// STAR walks outward from the two junction positions and counts how far the
/// sequences agree (`ReadAlign_chimericDetectionOld.cpp:260-294`). A repeat
/// there means the junction could be placed anywhere within it, which is what
/// columns 8 and 9 of `Chimeric.out.junction` report.
///
/// Three things the earlier version got wrong:
///
/// - It returned `(0, 0)` whenever the two segments were on different
///   chromosomes. STAR does no such check — the positions are genome-absolute,
///   so the comparison is well defined across chromosomes, and inter-chromosomal
///   fusions are precisely the interesting case. This alone zeroed both columns
///   for every chimera in an inter-chromosomal test set.
/// - It added `chr_start` to positions that were already absolute, so it read
///   the wrong bases whenever it read any.
/// - It ignored strand. STAR walks a reverse segment in the opposite direction
///   and complements the base.
///
/// `chimJ0` is the first base past the leading segment's block, `chimJ1` the
/// base before the trailing segment's; both follow STAR's `EX_G`/`EX_L`
/// arithmetic at `:244-257`. Returns `(chimRepeat0, chimRepeat1)`, the order the
/// junction file prints them in.
pub fn chim_repeat_lengths(
    genome: &Genome,
    leading: &crate::chimeric::ChimericSegment,
    trailing: &crate::chimeric::ChimericSegment,
) -> (u32, u32) {
    const MAX_REPEAT: u64 = 100; // STAR's loop bound

    let ex0 = leading.junction_exon(true);
    let ex1 = trailing.junction_exon(false);
    let (rev0, rev1) = (leading.is_reverse, trailing.is_reverse);

    // First base past the leading block / base before the trailing block.
    let chim_j0 = if rev0 {
        ex0.genome_start.checked_sub(1)
    } else {
        Some(ex0.genome_end)
    };
    let chim_j1 = if rev1 {
        Some(ex1.genome_end)
    } else {
        ex1.genome_start.checked_sub(1)
    };
    let (Some(j0), Some(j1)) = (chim_j0, chim_j1) else {
        return (0, 0);
    };

    // A base read in the segment's own orientation: reverse segments are walked
    // the other way and complemented, exactly as STAR does.
    let base = |pos: Option<u64>, rev: bool| -> Option<u8> {
        let b = genome.get_base(pos?)?;
        Some(if rev && b < 4 { 3 - b } else { b })
    };
    let step = |p: u64, delta: i64| -> Option<u64> {
        if delta < 0 {
            p.checked_sub(delta.unsigned_abs())
        } else {
            p.checked_add(delta as u64)
        }
    };

    // Forward: chimRepeat1.
    let mut repeat1 = 0u32;
    for jr in 0..MAX_REPEAT {
        let jr = jr as i64;
        let b0 = base(step(j0, if rev0 { -jr } else { jr }), rev0);
        let b1 = base(step(j1, if rev1 { -(jr + 1) } else { jr + 1 }), rev1);
        match (b0, b1) {
            (Some(a), Some(b)) if a == b => repeat1 += 1,
            _ => break,
        }
    }

    // Reverse: chimRepeat0.
    let mut repeat0 = 0u32;
    for jr in 0..MAX_REPEAT {
        let jr = jr as i64;
        let b0 = base(step(j0, if rev0 { jr + 1 } else { -(jr + 1) }), rev0);
        let b1 = base(step(j1, if rev1 { jr } else { -jr }), rev1);
        match (b0, b1) {
            (Some(a), Some(b)) if a == b => repeat0 += 1,
            _ => break,
        }
    }

    (repeat0, repeat1)
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
            transform_blocks: None,
            sequence: seq.into(),
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

    /// Build a forward, single-exon segment spanning `[gs, ge)`.
    fn rep_seg(chr_idx: usize, gs: u64, ge: u64) -> crate::chimeric::ChimericSegment {
        use crate::chimeric::{ChimericSegment, ExonSpan};
        ChimericSegment {
            chr_idx,
            genome_start: gs,
            genome_end: ge,
            is_reverse: false,
            read_start: 0,
            read_end: (ge - gs) as usize,
            cigar: Vec::new(),
            score: 0,
            n_mismatch: 0,
            first_exon: ExonSpan {
                genome_start: gs,
                genome_end: ge,
                read_start: 0,
                read_end: (ge - gs) as usize,
            },
            last_exon: ExonSpan {
                genome_start: gs,
                genome_end: ge,
                read_start: 0,
                read_end: (ge - gs) as usize,
            },
        }
    }

    #[test]
    fn repeat_length_is_zero_when_the_flanks_disagree() {
        // ACGTAC. Leading ends at 2 (so chimJ0 = 2), trailing starts at 4
        // (chimJ1 = 3). Forward compares G[2]=G vs G[4]=A, reverse G[1]=C vs
        // G[3]=T — neither agrees, so there is no repeat either side.
        let genome = mock_genome_with_sequence(vec![0, 1, 2, 3, 0, 1]);
        let (r0, r1) = chim_repeat_lengths(&genome, &rep_seg(0, 0, 2), &rep_seg(0, 4, 6));
        assert_eq!((r0, r1), (0, 0));
    }

    #[test]
    fn repeat_length_counts_a_homopolymer_either_side_of_the_junction() {
        // CAAAAAAG: the junction sits inside a run of six As, so it could be
        // placed anywhere within it — which is what the repeat columns report.
        let genome = mock_genome_with_sequence(vec![1, 0, 0, 0, 0, 0, 0, 2]);
        // chimJ0 = 4, chimJ1 = 6.
        let (r0, r1) = chim_repeat_lengths(&genome, &rep_seg(0, 0, 4), &rep_seg(0, 7, 8));
        assert_eq!(r0, 3, "three bases of agreement walking back");
        assert_eq!(r1, 0, "G[4]=A vs G[7]=G stops the forward walk immediately");
    }

    #[test]
    fn repeat_length_is_computed_across_chromosomes_too() {
        // STAR applies no same-chromosome check: the positions are
        // genome-absolute, so the comparison is well defined either way, and an
        // inter-chromosomal fusion is exactly the case worth reporting. This
        // previously returned (0, 0) for any cross-chromosome pair, which
        // zeroed both columns for every chimera in an inter-chromosomal set.
        let genome = mock_genome_with_sequence(vec![1, 0, 0, 0, 0, 0, 0, 2]);
        let same = chim_repeat_lengths(&genome, &rep_seg(0, 0, 4), &rep_seg(0, 7, 8));
        let cross = chim_repeat_lengths(&genome, &rep_seg(0, 0, 4), &rep_seg(1, 7, 8));
        assert_eq!(
            same, cross,
            "the chromosome index must not change the result"
        );
        assert_ne!(cross, (0, 0), "and it must not be short-circuited to zero");
    }
}
