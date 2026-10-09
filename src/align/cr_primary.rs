//! `--soloOutLayout CellRanger`: which alignment of a multimapper `cellranger count`
//! reports, and with what MAPQ.
//!
//! `cellranger count` keeps one BAM record per read. Its annotation stage
//! (`tx_annotation`, "MAPQ adjustment" of the cr-gex algorithm) looks at every
//! alignment STAR found: an alignment is *transcriptomic* when it is compatible with
//! an annotated transcript on the sense strand (every aligned block inside an exon,
//! every splice matching the transcript's junction). When all transcriptomic
//! alignments of the read belong to a single gene, the read is rescued: the
//! transcriptomic alignment with the best score becomes the primary record, whatever
//! the other alignments scored, and its MAPQ is 255. Otherwise STAR's own choice
//! stands, with STAR's MAPQ (3 for two loci, 1 for three or four, 0 above).

use crate::align::transcript::Transcript;
use crate::quant::transcriptome::{TranscriptomeIndex, align_to_transcripts};

/// MAPQ of a rescued alignment (CellRanger's "confidently mapped").
pub const RESCUED_MAPQ: u8 = 255;

/// The alignment `cellranger count` rescues as the primary, if any
/// (`tx_annotation::read::rescue_alignments_se`): the first alignment in STAR's output
/// order (window order, `star_order`) that is transcriptomic with a single gene, provided
/// every transcriptomic alignment of the read belongs to that one gene.
pub fn rescued_primary(transcripts: &[Transcript], tx: &TranscriptomeIndex) -> Option<usize> {
    if transcripts.len() < 2 {
        return None;
    }
    let mut order: Vec<usize> = (0..transcripts.len()).collect();
    order.sort_by_key(|&i| transcripts[i].star_order);
    let mut seen: Vec<u32> = Vec::new();
    let mut promote: Option<usize> = None;
    for &i in &order {
        let t = &transcripts[i];
        let mut genes: Vec<u32> = Vec::new();
        for projected in align_to_transcripts(t, tx, 0) {
            let tr = projected.chr_idx;
            // Sense: a `+` transcript takes a forward alignment, a `-` one a reverse.
            if (tx.tr_strand[tr] == 2) != t.is_reverse {
                continue;
            }
            let g = tx.tr_gene_idx[tr];
            if !genes.contains(&g) {
                genes.push(g);
            }
        }
        if genes.is_empty() {
            continue;
        }
        if genes.len() == 1 {
            promote = promote.or(Some(i));
        }
        for g in genes {
            if !seen.contains(&g) {
                seen.push(g);
            }
        }
    }
    if seen.len() > 1 { None } else { promote }
}

/// Reduce a read's records to the single one `cellranger count` writes: the rescued
/// alignment (primary flag, MAPQ 255) or, without a rescue, STAR's primary. NH and HI
/// keep describing the full set of alignments.
pub fn keep_primary_only(
    records: &mut Vec<noodles::sam::alignment::RecordBuf>,
    transcripts: &[Transcript],
    rescued: Option<usize>,
) {
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    let keep = match rescued {
        Some(i) if i < records.len() => Some(i),
        _ => records
            .iter()
            .position(|r| !r.flags().contains(Flags::SECONDARY)),
    };
    let Some(keep) = keep else { return };
    let mut record = records.swap_remove(keep);
    // HI is the alignment's place in STAR's output order (window order).
    if let Some(t) = transcripts.get(keep) {
        let hi = 1 + transcripts
            .iter()
            .filter(|o| o.star_order < t.star_order)
            .count();
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        record
            .data_mut()
            .insert(Tag::new(b'H', b'I'), Value::Int32(hi as i32));
    }
    if rescued.is_some() {
        let mut flags = record.flags();
        flags.remove(Flags::SECONDARY);
        *record.flags_mut() = flags;
        *record.mapping_quality_mut() = MappingQuality::new(RESCUED_MAPQ);
    }
    records.clear();
    records.push(record);
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::{Flags, MappingQuality};

    fn record(flags: Flags, mapq: u8) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.flags_mut() = flags;
        *r.mapping_quality_mut() = MappingQuality::new(mapq);
        r
    }

    fn transcript(star_order: u32) -> Transcript {
        Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 10,
            is_reverse: false,
            exons: vec![],
            cigar: vec![],
            score: 10,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            star_order,
        }
    }

    #[test]
    fn keeps_stars_primary_without_a_rescue() {
        let mut records = vec![record(Flags::empty(), 3), record(Flags::SECONDARY, 3)];
        let transcripts = [transcript(1), transcript(0)];
        keep_primary_only(&mut records, &transcripts, None);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].mapping_quality().map(u8::from), Some(3));
        assert!(!records[0].flags().contains(Flags::SECONDARY));
    }

    #[test]
    fn rescued_alignment_becomes_primary_with_mapq_255() {
        let mut records = vec![record(Flags::empty(), 3), record(Flags::SECONDARY, 3)];
        let transcripts = [transcript(1), transcript(0)];
        keep_primary_only(&mut records, &transcripts, Some(1));
        assert_eq!(records.len(), 1);
        assert!(!records[0].flags().contains(Flags::SECONDARY));
        // 255 is "missing" in BAM and is what CellRanger writes for a confident read.
        assert_eq!(records[0].mapping_quality(), None);
    }
}
