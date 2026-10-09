//! CellRanger's read annotation against the transcriptome (`--soloOutLayout
//! CellRanger`).
//!
//! CellRanger does not count the way STARsolo does. It annotates every
//! alignment of a read against the reference's transcripts and genes, then
//! decides from the whole alignment set which alignment is the confident one.
//! The rules below were recovered by comparing against the `RE`, `TX`, `AN`,
//! `GX`, `xf` and MAPQ values of a real `cellranger count` 10.0.0 BAM
//! (`pbmc_1k_v3`, `refdata-gex-GRCh38-2024-A`); on 395k uniquely mapped reads
//! the per-alignment rules reproduce those tags for >99.9% of alignments.
//!
//! Per alignment, with *blocks* the aligned stretches separated by `N`
//! (deletions stay inside a block):
//!
//! * A transcript is **compatible** when the first block starts inside one of
//!   its exons, every block lies within one exon, and every junction joins the
//!   end of an exon to the start of the next one. Soft clips are ignored. The
//!   alignment is then *transcriptomic*: `RE:E`, and the compatible transcripts
//!   on the read's strand go to `TX`, those on the other strand to `AN`.
//! * Otherwise a gene is **hit** when every block has at least half its bases
//!   inside the gene's span. `TX`/`AN` then carry `gene,+` / `gene,-` (sense /
//!   antisense) and the region is `E` when some transcript of a hit gene has at
//!   least half of the aligned bases in its exons and every block touches an
//!   exon, `N` otherwise. No hit gene: `I`.
//!
//! Per read:
//!
//! * One alignment: it is the confident one (MAPQ 255).
//! * Several: the read is rescued, to MAPQ 255 on its first sense-transcriptomic
//!   alignment, when the sense-transcriptomic alignments all belong to a single
//!   gene. Otherwise the read is multi-mapped and not confident.
//! * A confident read counts for a gene when its `TX` names exactly one gene
//!   (intronic reads included, as `--include-introns`, CellRanger's default).

use crate::align::transcript::Transcript;
use crate::quant::transcriptome::TranscriptomeIndex;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};

/// CellRanger's `RE` tag: where an alignment falls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrRegion {
    Exonic,
    Intronic,
    Intergenic,
}

impl CrRegion {
    /// The `RE:A` character.
    pub fn tag_char(self) -> char {
        match self {
            Self::Exonic => 'E',
            Self::Intronic => 'N',
            Self::Intergenic => 'I',
        }
    }
}

struct CrTx {
    start: i64,
    end: i64,
    plus: bool,
    gene: u32,
    /// `(start, end inclusive, cumulative length before this exon)`, ascending.
    exons: Vec<(i64, i64, u32)>,
    len: u32,
}

/// The reference as the annotation needs it: transcripts per chromosome, and
/// gene spans.
pub struct CrModel {
    tx: Vec<CrTx>,
    /// Per chromosome: transcript indices ascending by start, and the running
    /// maximum of their ends.
    by_chr: Vec<(Vec<u32>, Vec<i64>)>,
    /// Transcript ids, gene ids and names (for the tags).
    pub tx_ids: Vec<String>,
    pub gene_ids: Vec<String>,
    pub gene_names: Vec<String>,
    /// Rank of each transcript / gene when ordered by id, since CellRanger lists
    /// ids in that order.
    tx_rank: Vec<u32>,
    gene_rank: Vec<u32>,
}

/// A transcript-level hit of one alignment.
#[derive(Debug, Clone, Copy)]
pub struct TxHit {
    pub tx: u32,
    pub gene: u32,
    /// 0-based position of the alignment start in transcript space (for a `-`
    /// transcript, measured from its 5' end).
    pub pos: u32,
}

/// What one alignment was annotated as.
#[derive(Debug, Clone)]
pub struct AlnAnnot {
    pub region: CrRegion,
    /// Compatible transcripts on the read's strand, ascending by id.
    pub tx_sense: Vec<TxHit>,
    /// Compatible transcripts on the opposite strand.
    pub tx_anti: Vec<TxHit>,
    /// Genes hit at gene level (only when no transcript is compatible).
    pub gene_sense: Vec<u32>,
    pub gene_anti: Vec<u32>,
}

impl AlnAnnot {
    /// Genes named by the `TX` tag (and so by `GX`), ascending by id.
    fn tx_genes(&self, model: &CrModel) -> Vec<u32> {
        let mut g: Vec<u32> = if self.tx_sense.is_empty() {
            self.gene_sense.clone()
        } else {
            self.tx_sense.iter().map(|h| h.gene).collect()
        };
        g.sort_unstable_by_key(|&x| model.gene_rank[x as usize]);
        g.dedup();
        g
    }
}

/// What a whole read was annotated as.
#[derive(Debug, Clone)]
pub struct ReadAnnot {
    pub alns: Vec<AlnAnnot>,
    /// The alignment CellRanger keeps as the read's primary record.
    pub primary: usize,
    /// MAPQ 255: unique, or rescued.
    pub conf: bool,
    /// Genes of the primary's `TX`, ascending by id.
    pub genes: Vec<u32>,
    /// The gene the read counts for: confident and exactly one gene.
    pub gene: Option<u32>,
}

impl ReadAnnot {
    /// `Reads Mapped Antisense to Gene`: confident, with an antisense hit (a
    /// transcript or, for an intronic read, a gene) on the primary and nothing
    /// on the sense strand.
    pub fn antisense(&self) -> bool {
        let a = &self.alns[self.primary];
        self.conf
            && (!a.tx_anti.is_empty() || !a.gene_anti.is_empty())
            && a.tx_sense.is_empty()
            && a.gene_sense.is_empty()
    }

    /// CellRanger's `is_conf_mapped_unique_txomic`: confidently mapped to one
    /// gene with a transcript-level sense hit.
    pub fn txomic(&self) -> bool {
        self.conf && self.gene.is_some() && !self.alns[self.primary].tx_sense.is_empty()
    }

    pub fn primary_region(&self) -> CrRegion {
        self.alns[self.primary].region
    }
}

/// A run of aligned bases between two `N`s: `[start, end)` on the genome.
struct Segment {
    start: i64,
    end: i64,
}

/// `get_cigar_segments`: the aligned stretches of a CIGAR, split at `N`.
fn segments_of(cigar: &[Op], start: i64) -> Vec<Segment> {
    let mut out = Vec::new();
    let mut cur = Segment { start, end: start };
    for op in cigar {
        let n = op.len() as i64;
        match op.kind() {
            Kind::Skip => {
                let next = cur.end + n;
                out.push(std::mem::replace(
                    &mut cur,
                    Segment {
                        start: next,
                        end: next,
                    },
                ));
            }
            Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                cur.end += n;
            }
            _ => {}
        }
    }
    out.push(cur);
    out
}

/// `get_overlap`: CellRanger's fraction of the read interval covered by the
/// reference one. Its end coordinates are mixed inclusive and exclusive, and
/// that is kept.
fn overlap_frac(read_start: i64, read_end: i64, ref_start: i64, ref_end: i64) -> f64 {
    let bases = (ref_end.min(read_end) - ref_start.max(read_start)).max(0);
    bases as f64 / (read_end - read_start) as f64
}

/// `is_read_exonic` with `region_min_overlap` 0.5: every segment overlaps, by at
/// least half, the first exon that ends right of its start.
fn is_read_exonic(segments: &[Segment], exons: &[(i64, i64, u32)]) -> bool {
    segments.iter().all(|seg| {
        let idx = exons.partition_point(|e| e.1 <= seg.start);
        let Some(e) = exons.get(idx) else {
            return false;
        };
        overlap_frac(seg.start, seg.end, e.0, e.1) >= 0.5
    })
}

impl CrModel {
    pub fn from_transcriptome(tr: &TranscriptomeIndex, n_chr: usize) -> Self {
        let n_tx = tr.n_transcripts();
        let mut tx = Vec::with_capacity(n_tx);
        for i in 0..n_tx {
            let mut cum = 0u32;
            let exons: Vec<(i64, i64, u32)> = tr.tr_exons[i]
                .iter()
                .map(|e| {
                    let r = (e.genome_start as i64, e.genome_end as i64 - 1, cum);
                    cum += (e.genome_end - e.genome_start) as u32;
                    r
                })
                .collect();
            let g = tr.tr_gene_idx[i];
            let (s, e) = (tr.tr_start[i] as i64, tr.tr_end[i] as i64 - 1);
            let plus = tr.tr_strand[i] != 2;
            tx.push(CrTx {
                start: s,
                end: e,
                plus,
                gene: g,
                exons,
                len: cum,
            });
        }
        let mut by_chr: Vec<(Vec<u32>, Vec<i64>)> = vec![(Vec::new(), Vec::new()); n_chr];
        let mut per: Vec<Vec<u32>> = vec![Vec::new(); n_chr];
        for i in 0..n_tx {
            if let Some(v) = per.get_mut(tr.tr_chr_idx[i]) {
                v.push(i as u32);
            }
        }
        for (c, mut v) in per.into_iter().enumerate() {
            v.sort_by_key(|&i| (tx[i as usize].start, tx[i as usize].end));
            let mut mx = Vec::with_capacity(v.len());
            let mut m = i64::MIN;
            for &i in &v {
                m = m.max(tx[i as usize].end);
                mx.push(m);
            }
            by_chr[c] = (v, mx);
        }
        let rank = |ids: &[String]| -> Vec<u32> {
            let mut order: Vec<usize> = (0..ids.len()).collect();
            order.sort_by(|&a, &b| ids[a].cmp(&ids[b]));
            let mut r = vec![0u32; ids.len()];
            for (k, &i) in order.iter().enumerate() {
                r[i] = k as u32;
            }
            r
        };
        Self {
            tx,
            by_chr,
            tx_rank: rank(&tr.tr_ids),
            gene_rank: rank(&tr.gene_ids),
            tx_ids: tr.tr_ids.clone(),
            gene_ids: tr.gene_ids.clone(),
            gene_names: tr.gene_names.clone(),
        }
    }

    /// Transcripts overlapping the alignment's extent, as `annotate_alignment`
    /// walks them: start at or before the read's end, end at or after its start.
    fn overlapping(&self, chr: usize, start: i64, end: i64) -> Vec<u32> {
        let Some((order, mx)) = self.by_chr.get(chr) else {
            return Vec::new();
        };
        let hi = order.partition_point(|&i| self.tx[i as usize].start <= end);
        let mut out = Vec::new();
        for k in (0..hi).rev() {
            if mx[k] < start {
                break;
            }
            let t = &self.tx[order[k] as usize];
            if t.end >= start {
                out.push(order[k]);
            }
        }
        out
    }

    /// `TranscriptAnnotator::annotate_alignment` (tx_annotation/src/transcript.rs).
    pub fn annotate(&self, t: &Transcript) -> AlnAnnot {
        let mut out = AlnAnnot {
            region: CrRegion::Intergenic,
            tx_sense: Vec::new(),
            tx_anti: Vec::new(),
            gene_sense: Vec::new(),
            gene_anti: Vec::new(),
        };
        let start = t.genome_start as i64;
        let alen: i64 = t
            .cigar
            .iter()
            .filter(|o| {
                matches!(
                    o.kind(),
                    Kind::Match
                        | Kind::Deletion
                        | Kind::Skip
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                )
            })
            .map(|o| o.len() as i64)
            .sum();
        if alen == 0 {
            return out;
        }
        let end = start + alen - 1;
        let segments = segments_of(&t.cigar, start);
        let read_reverse = t.is_reverse;

        // (transcript, region, tx alignment position if transcript-compatible)
        struct Aln {
            gene: u32,
            antisense: bool,
            tx: Option<TxHit>,
        }
        let mut any_exonic = false;
        let mut any_intronic = false;
        let mut alns: Vec<(u32, Aln)> = Vec::new();
        for i in self.overlapping(t.chr_idx, start, end) {
            let tx = &self.tx[i as usize];
            let is_exonic = is_read_exonic(&segments, &tx.exons);
            let is_intronic = !is_exonic && overlap_frac(start, end, tx.start, tx.end) >= 1.0;
            if !is_exonic && !is_intronic {
                continue;
            }
            let antisense = tx.plus == read_reverse;
            let tx_align = if is_exonic {
                Self::align_to_exons(tx, i, &segments, start, end)
            } else {
                None
            };
            if is_exonic {
                any_exonic = true;
            } else {
                any_intronic = true;
            }
            alns.push((
                i,
                Aln {
                    gene: tx.gene,
                    antisense,
                    tx: tx_align,
                },
            ));
        }
        if any_exonic {
            out.region = CrRegion::Exonic;
            let txome_exists = alns.iter().any(|(_, a)| a.tx.is_some());
            let count_gene_level = !txome_exists;
            for (_, a) in &alns {
                match (&a.tx, a.antisense, count_gene_level) {
                    (Some(h), false, _) => out.tx_sense.push(*h),
                    (Some(h), true, _) => out.tx_anti.push(*h),
                    (None, false, true) => out.gene_sense.push(a.gene),
                    (None, true, true) => out.gene_anti.push(a.gene),
                    _ => {}
                }
            }
        } else if any_intronic {
            out.region = CrRegion::Intronic;
            for (_, a) in &alns {
                if a.antisense {
                    out.gene_anti.push(a.gene);
                } else {
                    out.gene_sense.push(a.gene);
                }
            }
        }
        out.tx_sense
            .sort_unstable_by_key(|h| self.tx_rank[h.tx as usize]);
        out.tx_anti
            .sort_unstable_by_key(|h| self.tx_rank[h.tx as usize]);
        out.gene_sense
            .sort_unstable_by_key(|&g| self.gene_rank[g as usize]);
        out.gene_sense.dedup();
        out.gene_anti
            .sort_unstable_by_key(|&g| self.gene_rank[g as usize]);
        out.gene_anti.dedup();
        out
    }

    /// `find_exons` + `align_junctions` with CellRanger's zero tolerances:
    /// the transcript position of the alignment when its blocks are exactly the
    /// exons (inner ends flush with exon ends, outer ends inside exons).
    fn align_to_exons(
        tx: &CrTx,
        tx_index: u32,
        segments: &[Segment],
        read_start: i64,
        read_end: i64,
    ) -> Option<TxHit> {
        let ex = &tx.exons;
        // First exon whose end lies right of the read start; last exon that
        // starts left of the read end (strictly).
        let ex_start = ex.partition_point(|e| e.1 <= read_start);
        let ex_end = ex.partition_point(|e| e.0 < read_end) as i64 - 1;
        if ex_start >= ex.len() || ex_end < 0 {
            return None;
        }
        let ex_end = ex_end as usize;
        if ex_end < ex_start {
            return None;
        }
        if read_start < ex[ex_start].0 || read_end > ex[ex_end].1 {
            return None;
        }
        let exons = &ex[ex_start..=ex_end];
        if exons.len() != segments.len() {
            return None;
        }
        let mut aligned = 0i64;
        for (k, (seg, e)) in segments.iter().zip(exons).enumerate() {
            aligned += e.1 - e.0 + 1;
            let start_diff = e.0 - seg.start;
            let end_diff = seg.end - e.1 - 1;
            if k == 0 {
                if start_diff > 0 {
                    return None;
                }
                aligned -= start_diff.abs();
            } else if start_diff != 0 {
                return None;
            }
            if k == segments.len() - 1 {
                if end_diff > 0 {
                    return None;
                }
                aligned -= end_diff.abs();
            } else if end_diff != 0 {
                return None;
            }
        }
        let ex_offset = (read_start - ex[ex_start].0).max(0);
        let mut offset = i64::from(ex[ex_start].2) + ex_offset;
        if !tx.plus {
            offset = i64::from(tx.len) - (offset + aligned);
        }
        Some(TxHit {
            tx: tx_index,
            gene: tx.gene,
            pos: offset as u32,
        })
    }

    /// Annotate a read from all its alignments (`transcripts` in the aligner's
    /// order; the first is STAR's primary).
    pub fn annotate_read(&self, transcripts: &[Transcript]) -> Option<ReadAnnot> {
        if transcripts.is_empty() {
            return None;
        }
        let alns: Vec<AlnAnnot> = transcripts.iter().map(|t| self.annotate(t)).collect();
        let (primary, conf) = if alns.len() == 1 {
            (0, true)
        } else {
            let mut genes: Vec<u32> = alns
                .iter()
                .flat_map(|a| a.tx_sense.iter().map(|h| h.gene))
                .collect();
            genes.sort_unstable();
            genes.dedup();
            if genes.len() == 1 {
                let first = alns
                    .iter()
                    .position(|a| !a.tx_sense.is_empty())
                    .unwrap_or(0);
                (first, true)
            } else {
                (0, false)
            }
        };
        let genes = alns[primary].tx_genes(self);
        let gene = (conf && genes.len() == 1).then(|| genes[0]);
        Some(ReadAnnot {
            alns,
            primary,
            conf,
            genes,
            gene,
        })
    }

    /// The `TX` tag value of an alignment (empty when it has none).
    pub fn tx_tag(&self, cigar: &[Op], a: &AlnAnnot) -> String {
        self.tx_like(cigar, &a.tx_sense, &a.gene_sense, '+')
    }

    /// The `AN` tag value of an alignment (empty when it has none).
    pub fn an_tag(&self, cigar: &[Op], a: &AlnAnnot) -> String {
        self.tx_like(cigar, &a.tx_anti, &a.gene_anti, '-')
    }

    fn tx_like(&self, cigar: &[Op], hits: &[TxHit], genes: &[u32], sym: char) -> String {
        use std::fmt::Write as _;
        let mut s = String::new();
        if hits.is_empty() {
            for &g in genes {
                if !s.is_empty() {
                    s.push(';');
                }
                let _ = write!(s, "{},{sym}", self.gene_ids[g as usize]);
            }
            return s;
        }
        for h in hits {
            if !s.is_empty() {
                s.push(';');
            }
            let tx = &self.tx[h.tx as usize];
            let _ = write!(s, "{},{sym}{},", self.tx_ids[h.tx as usize], h.pos);
            // The CIGAR in transcript space: no `N`, reversed on `-` transcripts.
            let mut ops: Vec<(usize, char)> = Vec::new();
            for op in cigar {
                let c = match op.kind() {
                    Kind::Skip => continue,
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                };
                match ops.last_mut() {
                    Some(l) if l.1 == c => l.0 += op.len(),
                    _ => ops.push((op.len(), c)),
                }
            }
            if !tx.plus {
                ops.reverse();
            }
            for (n, c) in ops {
                let _ = write!(s, "{n}{c}");
            }
        }
        s
    }

    /// The `GX` and `GN` tag values of a read's primary.
    pub fn gx_gn(&self, genes: &[u32]) -> (String, String) {
        let gx: Vec<&str> = genes
            .iter()
            .map(|&g| self.gene_ids[g as usize].as_str())
            .collect();
        let gn: Vec<&str> = genes
            .iter()
            .map(|&g| self.gene_names[g as usize].as_str())
            .collect();
        (gx.join(";"), gn.join(";"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::Exon;
    use noodles::sam::alignment::record::cigar::Op;

    /// Two genes on chromosome 0: `GA` (+) with transcripts `TA1` (exons
    /// 100..=199, 300..=399) and `TA2` (100..=199), and `GB` (-) with `TB`
    /// (exon 1000..=1099).
    fn model() -> CrModel {
        let mk = |plus, gene, exons: Vec<(i64, i64)>| {
            let mut cum = 0;
            let ex: Vec<(i64, i64, u32)> = exons
                .iter()
                .map(|&(s, e)| {
                    let r = (s, e, cum);
                    cum += (e - s + 1) as u32;
                    r
                })
                .collect();
            CrTx {
                start: ex[0].0,
                end: ex[ex.len() - 1].1,
                plus,
                gene,
                exons: ex,
                len: cum,
            }
        };
        let tx = vec![
            mk(true, 0, vec![(100, 199), (300, 399)]),
            mk(true, 0, vec![(100, 199)]),
            mk(false, 1, vec![(1000, 1099)]),
        ];
        CrModel {
            by_chr: vec![(vec![0, 1, 2], vec![399, 399, 1099])],
            tx,
            tx_ids: vec!["TA1".into(), "TA2".into(), "TB".into()],
            gene_ids: vec!["GA".into(), "GB".into()],
            gene_names: vec!["gA".into(), "gB".into()],
            tx_rank: vec![0, 1, 2],
            gene_rank: vec![0, 1],
        }
    }

    fn aln(start: u64, cigar: &[(Kind, usize)], reverse: bool) -> Transcript {
        let cigar: Vec<Op> = cigar.iter().map(|&(k, n)| Op::new(k, n)).collect();
        let end = start
            + cigar
                .iter()
                .filter(|o| o.kind().consumes_reference())
                .map(|o| o.len() as u64)
                .sum::<u64>();
        Transcript {
            chr_idx: 0,
            genome_start: start,
            genome_end: end,
            is_reverse: reverse,
            exons: vec![Exon {
                genome_start: start,
                genome_end: end,
                read_start: 0,
                read_end: (end - start) as usize,
                i_frag: 0,
            }],
            cigar,
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            star_order: 0,
        }
    }

    #[test]
    fn a_spliced_read_matching_a_junction_is_transcriptomic() {
        let m = model();
        // 20 bases ending exon 1, 30 starting exon 2: only TA1 has that junction.
        let t = aln(
            180,
            &[(Kind::Match, 20), (Kind::Skip, 100), (Kind::Match, 30)],
            false,
        );
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Exonic);
        assert_eq!(a.tx_sense.len(), 1);
        assert_eq!(m.tx_ids[a.tx_sense[0].tx as usize], "TA1");
        assert_eq!(a.tx_sense[0].pos, 80);
        assert!(a.tx_anti.is_empty());
        assert_eq!(m.tx_tag(&t.cigar, &a), "TA1,+80,50M");
    }

    #[test]
    fn an_exonic_read_on_the_wrong_strand_is_antisense() {
        let m = model();
        let t = aln(1010, &[(Kind::Match, 50)], false);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Exonic);
        assert!(a.tx_sense.is_empty());
        assert_eq!(m.an_tag(&t.cigar, &a), "TB,-40,50M");
        // A `-` transcript counts position from its 5' end.
        let r = m.annotate_read(&[t]).unwrap();
        assert!(r.conf && r.gene.is_none() && r.antisense());
    }

    #[test]
    fn an_intronic_read_names_the_gene_not_a_transcript() {
        let m = model();
        let t = aln(210, &[(Kind::Match, 50)], false);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Intronic);
        assert_eq!(m.tx_tag(&t.cigar, &a), "GA,+");
        let r = m.annotate_read(&[t]).unwrap();
        assert_eq!(r.gene, Some(0));
    }

    #[test]
    fn a_read_beside_every_gene_is_intergenic() {
        let m = model();
        let a = m.annotate(&aln(600, &[(Kind::Match, 50)], false));
        assert_eq!(a.region, CrRegion::Intergenic);
        assert!(a.tx_sense.is_empty() && a.gene_sense.is_empty() && a.gene_anti.is_empty());
    }

    #[test]
    fn exonic_needs_half_of_each_segment_in_its_exon_and_stays_gene_level_when_it_overhangs() {
        let m = model();
        // Mostly outside TB: intergenic.
        let a = m.annotate(&aln(950, &[(Kind::Match, 60)], false));
        assert_eq!(a.region, CrRegion::Intergenic);
        // 80 bases from 1050: 49 inside the exon, so exonic, but the end
        // overhangs the exon, so the hit is the gene's, not a transcript's.
        let t = aln(1050, &[(Kind::Match, 80)], true);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Exonic);
        assert!(a.tx_sense.is_empty() && a.tx_anti.is_empty());
        assert_eq!(a.gene_sense, vec![1]);
    }

    #[test]
    fn a_multimapper_with_one_gene_among_its_transcriptomic_alignments_is_rescued() {
        let m = model();
        let intergenic = aln(600, &[(Kind::Match, 50)], false);
        let txomic = aln(110, &[(Kind::Match, 50)], false);
        let r = m
            .annotate_read(&[intergenic.clone(), txomic.clone()])
            .unwrap();
        assert!(r.conf);
        assert_eq!(r.primary, 1);
        assert_eq!(r.gene, Some(0));
        // Two transcriptomic alignments in different genes: stays multi-mapped.
        let other = aln(1010, &[(Kind::Match, 50)], true);
        let r = m.annotate_read(&[txomic, other]).unwrap();
        assert!(!r.conf);
        assert_eq!(r.primary, 0);
        assert_eq!(r.gene, None);
    }
}
