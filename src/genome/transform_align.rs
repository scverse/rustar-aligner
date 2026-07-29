//! Align-time back-transform (`--genomeTransformOutput SAM`).
//!
//! STAR `Transcript_transformGenome.cpp` + `ReadAlign_transformGenome.cpp`.
//!
//! `--genomeTransformType Haploid` bakes a VCF's variants into the genome
//! sequence, so reads carrying those alleles align without mismatches. The
//! price is that every coordinate in the output then refers to a genome nobody
//! else has. `--genomeTransformOutput SAM` pays it back: each alignment is
//! mapped through the conversion blocks written at build time
//! (`transformGenomeBlocks.tsv`) onto the original genome, so an indel baked
//! into the transformed sequence reappears as an `I`/`D` CIGAR operation at the
//! original coordinates.
//!
//! Three things happen, in STAR's order:
//!
//! 1. **Remap.** Each exon is looked up in the block map and split wherever it
//!    straddles a block boundary, since the two halves land at unrelated
//!    original coordinates.
//! 2. **Merge.** Splits that turn out to be seamless in the original genome are
//!    collapsed again, and a block gap that is shorter on one side than the
//!    other has its shared part folded back into the neighbouring exon. What
//!    survives as a gap is a real indel.
//! 3. **Reclassify.** Junction motifs and annotation flags are recomputed
//!    against the *original* genome. A junction that was canonical in the
//!    transformed genome need not be canonical in the original one, and STAR
//!    reports what the original genome says.
//!
//! Only the SAM path is implemented here. `--genomeTransformOutput SJ` and
//! `Quant` are still rejected by parameter validation.

use noodles::sam::alignment::record::cigar::{self, op::Kind};

use crate::align::score::SpliceMotif;
use crate::align::transcript::{Exon, Transcript};
use crate::genome::Genome;
use crate::junction::SpliceJunctionDb;

/// Map a transcript from transformed-genome coordinates back to the original
/// genome.
///
/// `blocks` is the conversion map in `[transformed_start, length,
/// original_start]` order, ascending by `[0]`, terminated by the sentinel
/// [`super::transform::BLOCK_SENTINEL`] — the walk over trailing blocks relies
/// on a stop entry rather than a bounds test, as STAR's does.
///
/// Returns `None` when the transcript cannot be converted: no block covers its
/// start, or it ends up empty.
pub fn transform_transcript(
    orig: &Genome,
    orig_junctions: &SpliceJunctionDb,
    blocks: &[[u64; 3]],
    tr: &Transcript,
    align_intron_min: u64,
    align_intron_max: u64,
) -> Option<Transcript> {
    if tr.exons.is_empty() || blocks.is_empty() {
        return None;
    }

    let exons = remap_exons(blocks, &tr.exons)?;
    let exons = merge_adjacent(exons);

    let (motifs, annotated) = reclassify_junctions(
        orig,
        orig_junctions,
        &exons,
        tr.is_reverse,
        align_intron_min,
    );

    let (chr_idx, _) = orig.position_to_chr(exons[0].genome_start)?;
    let cigar = rebuild_cigar(tr, &exons, align_intron_min, align_intron_max);

    let genome_start = exons[0].genome_start;
    let genome_end = exons[exons.len() - 1].genome_end;
    let n_junction = motifs.len() as u32;
    // Gaps that are not junctions are indels. Counting them off the rebuilt
    // CIGAR rather than the exon list keeps the two consistent: the CIGAR is
    // what decides which gap is `N` and which is `D`.
    let n_gap = cigar
        .iter()
        .filter(|op| matches!(op.kind(), Kind::Insertion | Kind::Deletion))
        .count() as u32;

    Some(Transcript {
        chr_idx,
        genome_start,
        genome_end,
        is_reverse: tr.is_reverse,
        exons,
        cigar,
        score: tr.score,
        n_mismatch: tr.n_mismatch,
        n_gap,
        n_junction,
        junction_motifs: motifs,
        junction_annotated: annotated,
        read_seq: tr.read_seq.clone(),
    })
}

/// Step 1: per-exon remap, splitting at block boundaries.
fn remap_exons(blocks: &[[u64; 3]], exons: &[Exon]) -> Option<Vec<Exon>> {
    let mut out: Vec<Exon> = Vec::with_capacity(exons.len());
    for exon in exons {
        let len = exon.genome_end - exon.genome_start;
        if len == 0 {
            continue;
        }
        let g1 = exon.genome_start;
        let g2 = exon.genome_end - 1;
        let read_start = exon.read_start;

        // The last block starting at or before this exon. STAR reaches it with
        // `upper_bound` then `--`; with nothing at or before it, the exon lies
        // outside the map entirely and the alignment does not convert.
        let idx = blocks.partition_point(|b| b[0] <= g1);
        if idx == 0 {
            return None;
        }
        let mut ci = idx - 1;

        let b_start = blocks[ci][0];
        let b_end = blocks[ci][0] + blocks[ci][1] - 1;
        if g1 <= b_end {
            let piece = if g2 <= b_end { len } else { b_end - g1 + 1 };
            out.push(Exon {
                genome_start: blocks[ci][2] + g1 - b_start,
                genome_end: blocks[ci][2] + g1 - b_start + piece,
                read_start,
                read_end: read_start + piece as usize,
                i_frag: exon.i_frag,
            });
        }

        // Any further blocks this exon reaches into, each a separate piece.
        ci += 1;
        while ci < blocks.len() && g2 >= blocks[ci][0] {
            let piece = if g2 < blocks[ci][0] + blocks[ci][1] {
                g2 - blocks[ci][0] + 1
            } else {
                blocks[ci][1]
            };
            let r = read_start + (blocks[ci][0] - g1) as usize;
            out.push(Exon {
                genome_start: blocks[ci][2],
                genome_end: blocks[ci][2] + piece,
                read_start: r,
                read_end: r + piece as usize,
                i_frag: exon.i_frag,
            });
            ci += 1;
        }
    }
    if out.is_empty() { None } else { Some(out) }
}

/// Step 2: collapse the splits that are seamless in the original genome, and
/// fold the shared part of an uneven gap back into the exons around it.
fn merge_adjacent(exons: Vec<Exon>) -> Vec<Exon> {
    let mut out: Vec<Exon> = Vec::with_capacity(exons.len());
    for exon in exons {
        let Some(prev) = out.last_mut() else {
            out.push(exon);
            continue;
        };
        if prev.i_frag != exon.i_frag {
            out.push(exon);
            continue;
        }
        let gap_r = (exon.read_start - prev.read_end) as u64;
        let gap_g = exon.genome_start - prev.genome_end;
        if gap_r == gap_g {
            // Same gap on both sides: not an indel, so the split was an
            // artefact of the block boundary. Absorb it, gap included.
            prev.genome_end = exon.genome_end;
            prev.read_end = exon.read_end;
        } else {
            // Uneven gap: the smaller side is aligned sequence, not part of the
            // indel. Give it to the following exon so the remaining gap is the
            // indel alone.
            let shared = gap_r.min(gap_g);
            let mut exon = exon;
            if shared > 0 {
                exon.genome_start -= shared;
                exon.read_start -= shared as usize;
            }
            out.push(exon);
        }
    }
    out
}

/// Step 3: junction motifs and annotation flags, read off the original genome.
fn reclassify_junctions(
    orig: &Genome,
    orig_junctions: &SpliceJunctionDb,
    exons: &[Exon],
    is_reverse: bool,
    align_intron_min: u64,
) -> (Vec<SpliceMotif>, Vec<bool>) {
    let mut motifs = Vec::new();
    let mut annotated = Vec::new();
    for pair in exons.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        if a.i_frag != b.i_frag {
            continue; // the mate gap is not a junction
        }
        let gap_g = b.genome_start.saturating_sub(a.genome_end);
        let gap_r = b.read_start - a.read_end;
        if gap_r > 0 || gap_g < align_intron_min {
            continue; // insertion or deletion, not a junction
        }
        let motif = junction_motif(orig, a.genome_end, b.genome_start - 1);
        let is_annot = orig
            .position_to_chr(a.genome_end)
            .is_some_and(|(chr, offset)| {
                let end = b.genome_start - 1 - orig.chr_start[chr];
                let strand = if is_reverse { 2 } else { 1 };
                orig_junctions.is_annotated(chr, offset, end, strand)
                    || orig_junctions.is_annotated(chr, offset, end, 0)
            });
        motifs.push(motif);
        annotated.push(is_annot);
    }
    (motifs, annotated)
}

/// The donor/acceptor dinucleotides of intron `[j_s, j_e]` in the original
/// genome (STAR `Transcript_transformGenome.cpp:144-156`).
fn junction_motif(orig: &Genome, j_s: u64, j_e: u64) -> SpliceMotif {
    let at = |i: u64| orig.get_base(i).unwrap_or(5);
    match (at(j_s), at(j_s + 1), at(j_e.wrapping_sub(1)), at(j_e)) {
        (2, 3, 0, 2) => SpliceMotif::GtAg,
        (1, 3, 0, 1) => SpliceMotif::CtAc,
        (2, 1, 0, 2) => SpliceMotif::GcAg,
        (1, 3, 2, 1) => SpliceMotif::CtGc,
        (0, 3, 0, 1) => SpliceMotif::AtAc,
        (2, 3, 0, 3) => SpliceMotif::GtAt,
        _ => SpliceMotif::NonCanonical,
    }
}

/// Rebuild the CIGAR from the remapped exons.
///
/// The soft clips are the original transcript's: back-transforming moves where
/// the read aligns, never how much of it aligns.
fn rebuild_cigar(
    tr: &Transcript,
    exons: &[Exon],
    align_intron_min: u64,
    align_intron_max: u64,
) -> Vec<cigar::Op> {
    let mut ops: Vec<cigar::Op> = Vec::new();
    let push_match = |ops: &mut Vec<cigar::Op>, len: usize| {
        if len == 0 {
            return;
        }
        match ops.last_mut() {
            Some(op) if op.kind() == Kind::Match => {
                *op = cigar::Op::new(Kind::Match, op.len() + len);
            }
            _ => ops.push(cigar::Op::new(Kind::Match, len)),
        }
    };

    let [left_clip, right_clip] = tr.count_soft_clips();
    if left_clip > 0 {
        ops.push(cigar::Op::new(Kind::SoftClip, left_clip));
    }
    for (i, exon) in exons.iter().enumerate() {
        if i > 0 {
            let prev = &exons[i - 1];
            let gap_r = exon.read_start - prev.read_end;
            let gap_g = exon.genome_start - prev.genome_end;
            if gap_r > 0 {
                ops.push(cigar::Op::new(Kind::Insertion, gap_r));
            }
            if gap_g > 0 {
                let kind = if gap_g >= align_intron_min && gap_g <= align_intron_max {
                    Kind::Skip
                } else {
                    Kind::Deletion
                };
                ops.push(cigar::Op::new(kind, gap_g as usize));
            }
        }
        push_match(&mut ops, exon.read_end - exon.read_start);
    }
    if right_clip > 0 {
        ops.push(cigar::Op::new(Kind::SoftClip, right_clip));
    }
    ops
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::transform::BLOCK_SENTINEL;

    /// A genome of one chromosome whose bases are given as codes 0..3.
    fn genome_from_codes(codes: &[u8]) -> Genome {
        let n = codes.len() as u64;
        let mut seq = vec![0u8; 2 * n as usize];
        seq[..codes.len()].copy_from_slice(codes);
        Genome {
            transform_blocks: None,
            sequence: seq.into(),
            n_genome: n,
            n_genome_real: n,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![n],
            chr_start: vec![0, n],
        }
    }

    fn exon(genome_start: u64, len: u64, read_start: usize) -> Exon {
        Exon {
            genome_start,
            genome_end: genome_start + len,
            read_start,
            read_end: read_start + len as usize,
            i_frag: 0,
        }
    }

    fn transcript(exons: Vec<Exon>, cigar: Vec<cigar::Op>) -> Transcript {
        let genome_start = exons[0].genome_start;
        let genome_end = exons[exons.len() - 1].genome_end;
        Transcript {
            chr_idx: 0,
            genome_start,
            genome_end,
            is_reverse: false,
            exons,
            cigar,
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        }
    }

    /// Blocks for a genome where 5 original bases were deleted at original
    /// position 50: transformed [0,50) is original [0,50), and transformed
    /// [50,..) is original [55,..).
    fn deletion_blocks() -> Vec<[u64; 3]> {
        vec![[0, 50, 0], [50, 100, 55], BLOCK_SENTINEL]
    }

    #[test]
    fn an_exon_inside_one_block_is_shifted_by_that_block() {
        let orig = genome_from_codes(&[0u8; 200]);
        let jdb = SpliceJunctionDb::empty();
        let tr = transcript(vec![exon(60, 20, 0)], vec![cigar::Op::new(Kind::Match, 20)]);
        let out =
            transform_transcript(&orig, &jdb, &deletion_blocks(), &tr, 21, 1_000_000).unwrap();
        // transformed 60 is 10 into the second block, which starts at original 55.
        assert_eq!(out.genome_start, 65);
        assert_eq!(out.exons.len(), 1);
        assert_eq!(out.cigar_string(), "20M");
    }

    /// The point of the whole module: a variant baked into the genome comes
    /// back out as a CIGAR operation.
    #[test]
    fn an_exon_spanning_a_block_boundary_becomes_a_deletion() {
        let orig = genome_from_codes(&[0u8; 200]);
        let jdb = SpliceJunctionDb::empty();
        // One 20-base exon straddling the boundary at transformed 50.
        let tr = transcript(vec![exon(40, 20, 0)], vec![cigar::Op::new(Kind::Match, 20)]);
        let out =
            transform_transcript(&orig, &jdb, &deletion_blocks(), &tr, 21, 1_000_000).unwrap();
        assert_eq!(out.genome_start, 40);
        assert_eq!(out.exons.len(), 2, "the split survives as a real gap");
        // 10 bases, the 5 deleted bases, then 10 more.
        assert_eq!(out.cigar_string(), "10M5D10M");
        assert_eq!(out.n_gap, 1);
    }

    /// A gap the transform did not create must not be turned into one: when the
    /// two pieces are contiguous in the original genome too, the split
    /// disappears again.
    #[test]
    fn a_seamless_split_is_merged_back() {
        let orig = genome_from_codes(&[0u8; 200]);
        let jdb = SpliceJunctionDb::empty();
        // Two blocks that are adjacent in both coordinate systems: the split is
        // pure bookkeeping.
        let blocks = vec![[0, 50, 0], [50, 100, 50], BLOCK_SENTINEL];
        let tr = transcript(vec![exon(40, 20, 0)], vec![cigar::Op::new(Kind::Match, 20)]);
        let out = transform_transcript(&orig, &jdb, &blocks, &tr, 21, 1_000_000).unwrap();
        assert_eq!(out.exons.len(), 1);
        assert_eq!(out.cigar_string(), "20M");
        assert_eq!(out.n_gap, 0);
    }

    /// Soft clips describe the read, not the genome, so they survive unchanged.
    #[test]
    fn soft_clips_are_carried_through() {
        let orig = genome_from_codes(&[0u8; 200]);
        let jdb = SpliceJunctionDb::empty();
        let tr = transcript(
            vec![exon(60, 20, 5)],
            vec![
                cigar::Op::new(Kind::SoftClip, 5),
                cigar::Op::new(Kind::Match, 20),
                cigar::Op::new(Kind::SoftClip, 3),
            ],
        );
        let out =
            transform_transcript(&orig, &jdb, &deletion_blocks(), &tr, 21, 1_000_000).unwrap();
        assert_eq!(out.cigar_string(), "5S20M3S");
    }

    /// A junction is reclassified against the original genome, not the
    /// transformed one. Here the original bases spell GT..AG.
    #[test]
    fn a_junction_motif_is_read_from_the_original_genome() {
        let mut codes = vec![0u8; 200];
        // intron [80, 119]: GT at the donor, AG at the acceptor.
        codes[80] = 2;
        codes[81] = 3;
        codes[118] = 0;
        codes[119] = 2;
        let orig = genome_from_codes(&codes);
        let jdb = SpliceJunctionDb::empty();
        // Both exons sit in the identity block, so the coordinates survive.
        let blocks = vec![[0, 200, 0], BLOCK_SENTINEL];
        let tr = transcript(
            vec![exon(60, 20, 0), exon(120, 20, 20)],
            vec![
                cigar::Op::new(Kind::Match, 20),
                cigar::Op::new(Kind::Skip, 40),
                cigar::Op::new(Kind::Match, 20),
            ],
        );
        let out = transform_transcript(&orig, &jdb, &blocks, &tr, 21, 1_000_000).unwrap();
        assert_eq!(out.junction_motifs, vec![SpliceMotif::GtAg]);
        assert_eq!(out.cigar_string(), "20M40N20M");
        assert_eq!(out.n_junction, 1);
    }

    /// The same geometry over a genome that does not spell a canonical motif
    /// there. STAR reports what the original genome says, even when the
    /// transformed genome said otherwise.
    #[test]
    fn a_junction_canonical_only_in_the_transformed_genome_becomes_noncanonical() {
        let orig = genome_from_codes(&[0u8; 200]); // all A: no motif anywhere
        let jdb = SpliceJunctionDb::empty();
        let blocks = vec![[0, 200, 0], BLOCK_SENTINEL];
        let tr = transcript(
            vec![exon(60, 20, 0), exon(120, 20, 20)],
            vec![
                cigar::Op::new(Kind::Match, 20),
                cigar::Op::new(Kind::Skip, 40),
                cigar::Op::new(Kind::Match, 20),
            ],
        );
        let out = transform_transcript(&orig, &jdb, &blocks, &tr, 21, 1_000_000).unwrap();
        assert_eq!(out.junction_motifs, vec![SpliceMotif::NonCanonical]);
    }

    /// A gap shorter than `--alignIntronMin` is a deletion, not a junction, and
    /// carries no motif.
    #[test]
    fn a_short_gap_is_a_deletion_not_a_junction() {
        let orig = genome_from_codes(&[0u8; 200]);
        let jdb = SpliceJunctionDb::empty();
        let blocks = vec![[0, 200, 0], BLOCK_SENTINEL];
        let tr = transcript(
            vec![exon(60, 20, 0), exon(85, 20, 20)],
            vec![
                cigar::Op::new(Kind::Match, 20),
                cigar::Op::new(Kind::Deletion, 5),
                cigar::Op::new(Kind::Match, 20),
            ],
        );
        let out = transform_transcript(&orig, &jdb, &blocks, &tr, 21, 1_000_000).unwrap();
        assert!(out.junction_motifs.is_empty());
        assert_eq!(out.n_junction, 0);
        assert_eq!(out.cigar_string(), "20M5D20M");
    }

    /// An alignment starting before every block cannot be converted, and is
    /// dropped rather than guessed at.
    #[test]
    fn an_alignment_outside_the_block_map_does_not_convert() {
        let orig = genome_from_codes(&[0u8; 200]);
        let jdb = SpliceJunctionDb::empty();
        let blocks = vec![[100, 50, 100], BLOCK_SENTINEL];
        let tr = transcript(vec![exon(10, 20, 0)], vec![cigar::Op::new(Kind::Match, 20)]);
        assert!(transform_transcript(&orig, &jdb, &blocks, &tr, 21, 1_000_000).is_none());
    }
}
