//! Read annotation in the style of CellRanger (region, transcript and gene
//! hits, confident alignment, MAPQ, and the TX/AN/GX/GN/RE/mm/xf tags).
//!
//! Independent implementation from a behavioural specification.
//!
//! Coordinates are 0-based genome positions (global, as carried by
//! [`Transcript`]). Exon and alignment ends are handled as inclusive last
//! bases unless stated otherwise.

use crate::align::transcript::Transcript;
use crate::quant::transcriptome::TranscriptomeIndex;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use std::fmt::Write as _;

/// Region class of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrRegion {
    Exonic,
    Intronic,
    Intergenic,
}

impl CrRegion {
    /// The `RE` tag character.
    #[must_use]
    pub fn tag_char(self) -> char {
        match self {
            CrRegion::Exonic => 'E',
            CrRegion::Intronic => 'N',
            CrRegion::Intergenic => 'I',
        }
    }
}

/// A transcript hit: transcript, gene and 0-based position from the
/// transcript 5' end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TxHit {
    pub tx: u32,
    pub gene: u32,
    pub pos: u32,
}

/// Annotation of one alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlnAnnot {
    pub region: CrRegion,
    pub tx_sense: Vec<TxHit>,
    pub tx_anti: Vec<TxHit>,
    pub gene_sense: Vec<u32>,
    pub gene_anti: Vec<u32>,
}

impl AlnAnnot {
    fn intergenic() -> Self {
        AlnAnnot {
            region: CrRegion::Intergenic,
            tx_sense: Vec::new(),
            tx_anti: Vec::new(),
            gene_sense: Vec::new(),
            gene_anti: Vec::new(),
        }
    }
}

/// Annotation of one read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadAnnot {
    /// Per-alignment results, in aligner order.
    pub alns: Vec<AlnAnnot>,
    /// Index into `alns` of the chosen alignment.
    pub primary: usize,
    /// Whether the choice is confident.
    pub conf: bool,
    /// Genes named by the primary alignment (GX/GN), ascending by gene id.
    pub genes: Vec<u32>,
    /// The single counted gene, if any.
    pub gene: Option<u32>,
}

impl ReadAnnot {
    /// Confident, antisense hits only on the primary alignment.
    #[must_use]
    pub fn antisense(&self) -> bool {
        if !self.conf {
            return false;
        }
        let a = &self.alns[self.primary];
        let anti = !a.tx_anti.is_empty() || !a.gene_anti.is_empty();
        let sense = !a.tx_sense.is_empty() || !a.gene_sense.is_empty();
        anti && !sense
    }

    /// Confident, counted gene, and transcript-level sense hits.
    #[must_use]
    pub fn txomic(&self) -> bool {
        self.conf && self.gene.is_some() && !self.alns[self.primary].tx_sense.is_empty()
    }

    /// Region of the primary alignment.
    #[must_use]
    pub fn primary_region(&self) -> CrRegion {
        self.alns[self.primary].region
    }
}

/// Gene model used for annotation.
#[derive(Debug, Clone, Default)]
pub struct CrModel {
    pub tx_ids: Vec<String>,
    pub gene_ids: Vec<String>,
    pub gene_names: Vec<String>,
    tx_minus: Vec<bool>,
    tx_gene: Vec<u32>,
    tx_start: Vec<i64>,
    tx_end: Vec<i64>,
    tx_len: Vec<i64>,
    /// Exons as (first base, last base, cumulative length before).
    tx_exons: Vec<Vec<(i64, i64, i64)>>,
    /// Per chromosome: transcripts sorted by start, with running max end.
    by_chr: Vec<ChrIndex>,
}

#[derive(Debug, Clone, Default)]
struct ChrIndex {
    order: Vec<u32>,
    starts: Vec<i64>,
    max_end: Vec<i64>,
}

/// An alignment reduced to what annotation needs.
struct Aln {
    reverse: bool,
    start: i64,
    end: i64,
    /// Half-open reference segments between N operations.
    segs: Vec<(i64, i64)>,
}

fn aln_of(t: &Transcript) -> Option<Aln> {
    let mut segs = Vec::new();
    let mut pos = i64::try_from(t.genome_start).unwrap_or(0);
    let start = pos;
    let mut seg_start = pos;
    for op in &t.cigar {
        let n = i64::try_from(op.len()).unwrap_or(0);
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                pos += n;
            }
            Kind::Skip => {
                if pos > seg_start {
                    segs.push((seg_start, pos));
                }
                pos += n;
                seg_start = pos;
            }
            _ => {}
        }
    }
    if pos > seg_start {
        segs.push((seg_start, pos));
    }
    let first = segs.first()?.0;
    let last = segs.last()?.1;
    debug_assert!(first >= start);
    Some(Aln {
        reverse: t.is_reverse,
        start: first,
        end: last - 1,
        segs,
    })
}

impl CrModel {
    /// Build from the transcriptome index.
    #[must_use]
    pub fn from_transcriptome(t: &TranscriptomeIndex, n_chr: usize) -> CrModel {
        let n = t.tr_ids.len();
        let mut m = CrModel {
            tx_ids: t.tr_ids.clone(),
            gene_ids: t.gene_ids.clone(),
            gene_names: t.gene_names.clone(),
            ..CrModel::default()
        };
        let mut n_chr = n_chr;
        if let Some(mx) = t.tr_chr_idx.iter().max() {
            n_chr = n_chr.max(mx + 1);
        }
        let mut per_chr: Vec<Vec<u32>> = vec![Vec::new(); n_chr];
        for i in 0..n {
            m.tx_minus.push(t.tr_strand[i] == 2);
            m.tx_gene.push(t.tr_gene_idx[i]);
            m.tx_start.push(t.tr_start[i] as i64);
            m.tx_end.push(t.tr_end[i] as i64 - 1);
            m.tx_len.push(i64::from(t.tr_length[i]));
            m.tx_exons.push(
                t.tr_exons[i]
                    .iter()
                    .map(|e| {
                        (
                            e.genome_start as i64,
                            e.genome_end as i64 - 1,
                            i64::from(e.ex_len_cum),
                        )
                    })
                    .collect(),
            );
            per_chr[t.tr_chr_idx[i]].push(i as u32);
        }
        for mut order in per_chr {
            order.sort_by_key(|&i| (m.tx_start[i as usize], i));
            let starts: Vec<i64> = order.iter().map(|&i| m.tx_start[i as usize]).collect();
            let mut max_end = Vec::with_capacity(order.len());
            let mut cur = i64::MIN;
            for &i in &order {
                cur = cur.max(m.tx_end[i as usize]);
                max_end.push(cur);
            }
            m.by_chr.push(ChrIndex {
                order,
                starts,
                max_end,
            });
        }
        m
    }

    /// Transcripts on `chr` with `tx.start <= end` and `tx.end >= start`.
    fn candidates(&self, chr: usize, start: i64, end: i64) -> Vec<u32> {
        let Some(ci) = self.by_chr.get(chr) else {
            return Vec::new();
        };
        let hi = ci.starts.partition_point(|&s| s <= end);
        let mut out = Vec::new();
        let mut i = hi;
        while i > 0 {
            i -= 1;
            if ci.max_end[i] < start {
                break;
            }
            let tx = ci.order[i];
            if self.tx_end[tx as usize] >= start {
                out.push(tx);
            }
        }
        out
    }

    fn exonic(&self, tx: u32, a: &Aln) -> bool {
        let exons = &self.tx_exons[tx as usize];
        a.segs.iter().all(|&(sa, sb)| {
            let Some(&(xs, xe, _)) = exons.iter().find(|e| e.1 > sa) else {
                return false;
            };
            let den = sb - sa;
            let num = (xe.min(sb) - xs.max(sa)).max(0);
            den > 0 && 2 * num >= den
        })
    }

    fn intronic(&self, tx: u32, a: &Aln) -> bool {
        let ts = self.tx_start[tx as usize];
        let te = self.tx_end[tx as usize];
        a.end > a.start && te.min(a.end) - ts.max(a.start) >= a.end - a.start
    }

    /// Transcript-space position if the alignment is compatible with `tx`.
    fn compatible(&self, tx: u32, a: &Aln) -> Option<u32> {
        let exons = &self.tx_exons[tx as usize];
        let first = exons.iter().position(|e| e.1 > a.start)?;
        let last = exons.iter().rposition(|e| e.0 < a.end)?;
        if last < first {
            return None;
        }
        if a.start < exons[first].0 || a.end > exons[last].1 {
            return None;
        }
        let nseg = a.segs.len();
        if last - first + 1 != nseg {
            return None;
        }
        for (k, &(sa, sb)) in a.segs.iter().enumerate() {
            let (xs, xe, _) = exons[first + k];
            if k > 0 && sa != xs {
                return None;
            }
            if k + 1 < nseg && sb - 1 != xe {
                return None;
            }
            if k == 0 && sa < xs {
                return None;
            }
            if k + 1 == nseg && sb - 1 > xe {
                return None;
            }
        }
        let full: i64 = exons[first..=last].iter().map(|e| e.1 - e.0 + 1).sum();
        let aligned = full - (a.start - exons[first].0) - (exons[last].1 - a.end);
        let offset = exons[first].2 + (a.start - exons[first].0);
        let pos = if self.tx_minus[tx as usize] {
            self.tx_len[tx as usize] - (offset + aligned)
        } else {
            offset
        };
        u32::try_from(pos.max(0)).ok()
    }

    fn sort_dedup_genes(&self, v: &mut Vec<u32>) {
        v.sort_by(|&a, &b| self.gene_ids[a as usize].cmp(&self.gene_ids[b as usize]));
        v.dedup();
    }

    /// Annotate one alignment.
    #[must_use]
    pub fn annotate(&self, t: &Transcript) -> AlnAnnot {
        let Some(a) = aln_of(t) else {
            return AlnAnnot::intergenic();
        };
        let mut exonic = Vec::new();
        let mut intronic = Vec::new();
        for tx in self.candidates(t.chr_idx, a.start, a.end) {
            if self.exonic(tx, &a) {
                exonic.push(tx);
            } else if self.intronic(tx, &a) {
                intronic.push(tx);
            }
        }
        let is_anti = |tx: u32| self.tx_minus[tx as usize] != a.reverse;
        let mut out = AlnAnnot::intergenic();
        if !exonic.is_empty() {
            out.region = CrRegion::Exonic;
            let mut any = false;
            for &tx in &exonic {
                if let Some(pos) = self.compatible(tx, &a) {
                    any = true;
                    let hit = TxHit {
                        tx,
                        gene: self.tx_gene[tx as usize],
                        pos,
                    };
                    if is_anti(tx) {
                        out.tx_anti.push(hit);
                    } else {
                        out.tx_sense.push(hit);
                    }
                }
            }
            if !any {
                for &tx in exonic.iter().chain(intronic.iter()) {
                    let g = self.tx_gene[tx as usize];
                    if is_anti(tx) {
                        out.gene_anti.push(g);
                    } else {
                        out.gene_sense.push(g);
                    }
                }
            }
        } else if !intronic.is_empty() {
            out.region = CrRegion::Intronic;
            for &tx in &intronic {
                let g = self.tx_gene[tx as usize];
                if is_anti(tx) {
                    out.gene_anti.push(g);
                } else {
                    out.gene_sense.push(g);
                }
            }
        }
        let key = |h: &TxHit| &self.tx_ids[h.tx as usize];
        out.tx_sense.sort_by(|x, y| key(x).cmp(key(y)));
        out.tx_anti.sort_by(|x, y| key(x).cmp(key(y)));
        self.sort_dedup_genes(&mut out.gene_sense);
        self.sort_dedup_genes(&mut out.gene_anti);
        out
    }

    /// Genes named by an alignment (sense transcript hits, else sense genes).
    fn aln_genes(&self, a: &AlnAnnot) -> Vec<u32> {
        let mut g: Vec<u32> = if a.tx_sense.is_empty() {
            a.gene_sense.clone()
        } else {
            a.tx_sense.iter().map(|h| h.gene).collect()
        };
        self.sort_dedup_genes(&mut g);
        g
    }

    /// Annotate a read from its alignments (aligner order given by
    /// `star_order`, then index). `None` when there are no alignments.
    #[must_use]
    pub fn annotate_read(&self, alns: &[Transcript]) -> Option<ReadAnnot> {
        if alns.is_empty() {
            return None;
        }
        let mut order: Vec<usize> = (0..alns.len()).collect();
        order.sort_by_key(|&i| (alns[i].star_order, i));
        let ann: Vec<AlnAnnot> = order.iter().map(|&i| self.annotate(&alns[i])).collect();
        let (primary, conf) = if ann.len() == 1 {
            (0, true)
        } else {
            let mut genes: Vec<u32> = ann
                .iter()
                .flat_map(|a| a.tx_sense.iter().map(|h| h.gene))
                .collect();
            genes.sort_unstable();
            genes.dedup();
            if genes.len() == 1 {
                let p = ann.iter().position(|a| !a.tx_sense.is_empty()).unwrap_or(0);
                (p, true)
            } else {
                (order.iter().position(|&i| i == 0).unwrap_or(0), false)
            }
        };
        let genes = self.aln_genes(&ann[primary]);
        let gene = (conf && genes.len() == 1).then(|| genes[0]);
        Some(ReadAnnot {
            alns: ann,
            primary,
            conf,
            genes,
            gene,
        })
    }

    /// CIGAR text for a TX/AN entry: N removed, `=`/`X` as `M`, equal
    /// neighbours merged, reversed for a minus-strand transcript.
    fn hit_cigar(cigar: &[Op], minus: bool) -> String {
        let mut ops: Vec<(char, usize)> = Vec::new();
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
                Some((lc, n)) if *lc == c => *n += op.len(),
                _ => ops.push((c, op.len())),
            }
        }
        if minus {
            ops.reverse();
        }
        let mut s = String::new();
        for (c, n) in ops {
            let _ = write!(s, "{n}{c}");
        }
        s
    }

    fn tag(&self, cigar: &[Op], hits: &[TxHit], genes: &[u32], sign: char) -> String {
        if !hits.is_empty() {
            hits.iter()
                .map(|h| {
                    format!(
                        "{},{}{},{}",
                        self.tx_ids[h.tx as usize],
                        sign,
                        h.pos,
                        Self::hit_cigar(cigar, self.tx_minus[h.tx as usize])
                    )
                })
                .collect::<Vec<_>>()
                .join(";")
        } else {
            genes
                .iter()
                .map(|&g| format!("{},{}", self.gene_ids[g as usize], sign))
                .collect::<Vec<_>>()
                .join(";")
        }
    }

    /// The `TX` tag value (empty when there is nothing to write).
    #[must_use]
    pub fn tx_tag(&self, cigar: &[Op], a: &AlnAnnot) -> String {
        self.tag(cigar, &a.tx_sense, &a.gene_sense, '+')
    }

    /// The `AN` tag value (empty when there is nothing to write).
    #[must_use]
    pub fn an_tag(&self, cigar: &[Op], a: &AlnAnnot) -> String {
        self.tag(cigar, &a.tx_anti, &a.gene_anti, '-')
    }

    /// `GX` and `GN` values for a gene list.
    #[must_use]
    pub fn gx_gn(&self, genes: &[u32]) -> (String, String) {
        let gx = genes
            .iter()
            .map(|&g| self.gene_ids[g as usize].as_str())
            .collect::<Vec<_>>()
            .join(";");
        let gn = genes
            .iter()
            .map(|&g| self.gene_names[g as usize].as_str())
            .collect::<Vec<_>>()
            .join(";");
        (gx, gn)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ops(v: &[(Kind, usize)]) -> Vec<Op> {
        v.iter().map(|&(k, n)| Op::new(k, n)).collect()
    }

    /// Plus gene GA/TA: exons [100,149] [200,249]; minus gene GB/TB:
    /// exons [1000,1049] [1100,1149]; plus gene GC/TC: exon [2000,2099]
    /// and TD (gene GD, plus) exon [2000,2099].
    fn model() -> CrModel {
        let tx_ids: Vec<String> = ["TA", "TB", "TC", "TD"].map(String::from).to_vec();
        let ex = |v: &[(i64, i64)]| {
            let mut cum = 0;
            v.iter()
                .map(|&(s, e)| {
                    let r = (s, e, cum);
                    cum += e - s + 1;
                    r
                })
                .collect::<Vec<_>>()
        };
        let exons = vec![
            ex(&[(100, 149), (200, 249)]),
            ex(&[(1000, 1049), (1100, 1149)]),
            ex(&[(2000, 2099)]),
            ex(&[(2000, 2099)]),
        ];
        let tx_start: Vec<i64> = exons.iter().map(|e| e[0].0).collect();
        let tx_end: Vec<i64> = exons.iter().map(|e| e.last().unwrap().1).collect();
        let mut m = CrModel {
            tx_ids,
            gene_ids: ["GA", "GB", "GC", "GD"].map(String::from).to_vec(),
            gene_names: ["na", "nb", "nc", "nd"].map(String::from).to_vec(),
            tx_minus: vec![false, true, false, false],
            tx_gene: vec![0, 1, 2, 3],
            tx_len: vec![100, 100, 100, 100],
            tx_start,
            tx_end,
            tx_exons: exons,
            by_chr: Vec::new(),
        };
        let order: Vec<u32> = (0..4).collect();
        let starts = order.iter().map(|&i| m.tx_start[i as usize]).collect();
        let mut max_end = Vec::new();
        let mut cur = i64::MIN;
        for &i in &order {
            cur = cur.max(m.tx_end[i as usize]);
            max_end.push(cur);
        }
        m.by_chr.push(ChrIndex {
            order,
            starts,
            max_end,
        });
        m
    }

    fn aln(start: u64, rev: bool, cigar: &[(Kind, usize)], order: u32) -> Transcript {
        Transcript {
            chr_idx: 0,
            genome_start: start,
            genome_end: start,
            is_reverse: rev,
            exons: Vec::new(),
            cigar: ops(cigar),
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            star_order: order,
        }
    }

    const M: Kind = Kind::Match;
    const N: Kind = Kind::Skip;

    #[test]
    fn spliced_match() {
        let m = model();
        // 30 bases ending exon 1, 20 bases starting exon 2.
        let t = aln(120, false, &[(M, 30), (N, 50), (M, 20)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Exonic);
        assert_eq!(a.tx_sense.len(), 1);
        assert_eq!(a.tx_sense[0].pos, 20);
        assert_eq!(a.tx_anti, Vec::<TxHit>::new());
        assert_eq!(m.tx_tag(&t.cigar, &a), "TA,+20,50M");
        assert_eq!(m.an_tag(&t.cigar, &a), "");
    }

    #[test]
    fn exact_match_pos() {
        let m = model();
        let t = aln(150 + 50, false, &[(M, 50)], 0);
        let a = m.annotate(&t);
        assert_eq!(m.tx_tag(&t.cigar, &a), "TA,+50,50M");
    }

    #[test]
    fn antisense_minus() {
        let m = model();
        // Forward read on a minus transcript is sense; reverse read is
        // antisense here? plus==reverse => antisense; minus transcript with a
        // forward read: (false) == (false) => antisense.
        let t = aln(1110, false, &[(M, 30)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Exonic);
        assert_eq!(a.tx_sense, Vec::<TxHit>::new());
        assert_eq!(a.tx_anti.len(), 1);
        assert_eq!(m.an_tag(&t.cigar, &a), "TB,-10,30M");
        assert_eq!(m.tx_tag(&t.cigar, &a), "");
        let r = m.annotate_read(&[t]).unwrap();
        assert!(r.antisense());
        assert!(!r.txomic());
        assert_eq!(r.gene, None);
    }

    #[test]
    fn intronic() {
        let m = model();
        let t = aln(160, false, &[(M, 30)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Intronic);
        assert_eq!(a.tx_sense, Vec::<TxHit>::new());
        assert_eq!(m.tx_tag(&t.cigar, &a), "GA,+");
        let r = m.annotate_read(&[t]).unwrap();
        assert_eq!(r.gene, Some(0));
        assert_eq!(m.gx_gn(&r.genes), ("GA".into(), "na".into()));
        assert!(!r.txomic());
    }

    #[test]
    fn one_base_never_intronic() {
        let m = model();
        let t = aln(160, false, &[(M, 1)], 0);
        assert_eq!(m.annotate(&t).region, CrRegion::Intergenic);
    }

    #[test]
    fn intergenic() {
        let m = model();
        let t = aln(5000, false, &[(M, 50)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Intergenic);
        assert!(a.tx_sense.is_empty() && a.gene_sense.is_empty());
        let r = m.annotate_read(&[t]).unwrap();
        assert!(r.conf);
        assert_eq!(r.gene, None);
        assert_eq!(r.genes, Vec::<u32>::new());
    }

    #[test]
    fn zero_ref_length_is_intergenic() {
        let m = model();
        let t = aln(120, false, &[(Kind::SoftClip, 20)], 0);
        assert_eq!(m.annotate(&t).region, CrRegion::Intergenic);
    }

    #[test]
    fn overhang_gene_level_exonic() {
        let m = model();
        // Overhangs the end of exon 1 into the intron: exonic by the loose
        // test (>= half inside) but not transcript compatible.
        let t = aln(130, false, &[(M, 30)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.region, CrRegion::Exonic);
        assert_eq!(a.tx_sense, Vec::<TxHit>::new());
        assert_eq!(a.gene_sense, vec![0]);
        assert_eq!(m.tx_tag(&t.cigar, &a), "GA,+");
        let r = m.annotate_read(&[t]).unwrap();
        assert_eq!(r.gene, Some(0));
        assert_eq!(r.primary_region(), CrRegion::Exonic);
    }

    #[test]
    fn wrong_junction_not_compatible() {
        let m = model();
        // Junction one base off the exon end.
        let t = aln(120, false, &[(M, 29), (N, 51), (M, 20)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.tx_sense, Vec::<TxHit>::new());
        assert_eq!(a.gene_sense, vec![0]);
    }

    #[test]
    fn antisense_reverse_on_plus() {
        let m = model();
        let t = aln(120, true, &[(M, 30)], 0);
        let a = m.annotate(&t);
        // 120..149 inside exon 1 (30 bases), compatible, antisense.
        assert_eq!(a.tx_anti.len(), 1);
        assert_eq!(m.an_tag(&t.cigar, &a), "TA,-20,30M");
    }

    #[test]
    fn minus_strand_cigar_reversed_and_merged() {
        let m = model();
        let c = ops(&[
            (Kind::SoftClip, 3),
            (M, 10),
            (N, 5),
            (M, 10),
            (Kind::Insertion, 2),
            (Kind::SequenceMatch, 4),
        ]);
        assert_eq!(CrModel::hit_cigar(&c, false), "3S20M2I4M");
        assert_eq!(CrModel::hit_cigar(&c, true), "4M2I20M3S");
        let _ = m;
    }

    #[test]
    fn gene_and_tx_order_by_string() {
        let m = model();
        // TC and TD overlap: both compatible, genes GC and GD, two genes.
        let t = aln(2010, false, &[(M, 40)], 0);
        let a = m.annotate(&t);
        assert_eq!(a.tx_sense.iter().map(|h| h.tx).collect::<Vec<_>>(), [2, 3]);
        assert_eq!(m.tx_tag(&t.cigar, &a), "TC,+10,40M;TD,+10,40M");
        let r = m.annotate_read(&[t]).unwrap();
        assert_eq!(r.genes, vec![2, 3]);
        assert_eq!(r.gene, None);
        assert_eq!(m.gx_gn(&r.genes), ("GC;GD".into(), "nc;nd".into()));
    }

    #[test]
    fn rescue_single_gene() {
        let m = model();
        // Primary (order 0) intergenic, second aligns in GA.
        let a0 = aln(5000, false, &[(M, 50)], 0);
        let a1 = aln(100, false, &[(M, 50)], 1);
        let r = m.annotate_read(&[a0, a1]).unwrap();
        assert!(r.conf);
        assert_eq!(r.primary, 1);
        assert_eq!(r.gene, Some(0));
        assert!(r.txomic());
    }

    #[test]
    fn two_gene_multimapper_unrescued() {
        let m = model();
        let a0 = aln(100, false, &[(M, 50)], 0);
        let a1 = aln(2000, false, &[(M, 50)], 1);
        // Input index 0 is the aligner's primary even though star_order of
        // index 0 is larger.
        let a0b = Transcript {
            star_order: 1,
            ..a0.clone()
        };
        let a1b = Transcript {
            star_order: 0,
            ..a1.clone()
        };
        let r = m.annotate_read(&[a0b, a1b]).unwrap();
        assert!(!r.conf);
        assert_eq!(r.alns.len(), 2);
        assert_eq!(r.primary, 1);
        assert_eq!(r.gene, None);
        let r = m.annotate_read(&[a0, a1]).unwrap();
        assert!(!r.conf);
        assert_eq!(r.primary, 0);
    }

    #[test]
    fn same_gene_multimapper_rescues_first_sense() {
        let m = model();
        let a0 = aln(5000, false, &[(M, 50)], 0);
        let a1 = aln(100, false, &[(M, 50)], 1);
        let a2 = aln(200, false, &[(M, 50)], 2);
        let r = m.annotate_read(&[a0, a1, a2]).unwrap();
        assert_eq!(r.primary, 1);
        assert!(r.conf);
    }

    #[test]
    fn empty_read() {
        assert!(model().annotate_read(&[]).is_none());
    }
}
