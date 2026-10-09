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
    start: u64,
    end: u64,
    plus: bool,
    gene: u32,
    /// `(start, end_exclusive, cumulative length before this exon)`, ascending.
    exons: Vec<(u64, u64, u32)>,
    len: u32,
}

struct CrGene {
    start: u64,
    end: u64,
    plus: bool,
}

/// The reference as the annotation needs it: transcripts per chromosome, and
/// gene spans.
pub struct CrModel {
    tx: Vec<CrTx>,
    genes: Vec<CrGene>,
    /// Per chromosome: transcript indices ascending by start, and the running
    /// maximum of their ends.
    by_chr: Vec<(Vec<u32>, Vec<u64>)>,
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

    pub fn primary_region(&self) -> CrRegion {
        self.alns[self.primary].region
    }
}

fn blocks_of(t: &Transcript) -> Vec<(u64, u64)> {
    let mut out = Vec::with_capacity(t.n_junction as usize + 1);
    let mut pos = t.genome_start;
    let mut cur: Option<(u64, u64)> = None;
    for op in &t.cigar {
        let n = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                cur = Some(match cur {
                    Some((s, _)) => (s, pos + n),
                    None => (pos, pos + n),
                });
                pos += n;
            }
            Kind::Skip => {
                if let Some(b) = cur.take() {
                    out.push(b);
                }
                pos += n;
            }
            _ => {}
        }
    }
    if let Some(b) = cur {
        out.push(b);
    }
    out
}

fn overlap(a: (u64, u64), b: (u64, u64)) -> u64 {
    a.1.min(b.1).saturating_sub(a.0.max(b.0))
}

impl CrModel {
    pub fn from_transcriptome(tr: &TranscriptomeIndex, n_chr: usize) -> Self {
        let n_tx = tr.n_transcripts();
        let mut tx = Vec::with_capacity(n_tx);
        let mut gspan: Vec<Option<(u64, u64)>> = vec![None; tr.gene_ids.len()];
        let mut gplus = vec![true; tr.gene_ids.len()];
        for i in 0..n_tx {
            let mut cum = 0u32;
            let exons: Vec<(u64, u64, u32)> = tr.tr_exons[i]
                .iter()
                .map(|e| {
                    let r = (e.genome_start, e.genome_end, cum);
                    cum += (e.genome_end - e.genome_start) as u32;
                    r
                })
                .collect();
            let g = tr.tr_gene_idx[i];
            let (s, e) = (tr.tr_start[i], tr.tr_end[i]);
            let sp = &mut gspan[g as usize];
            *sp = Some(match *sp {
                Some((lo, hi)) => (lo.min(s), hi.max(e)),
                None => (s, e),
            });
            let plus = tr.tr_strand[i] != 2;
            gplus[g as usize] = plus;
            tx.push(CrTx {
                start: s,
                end: e,
                plus,
                gene: g,
                exons,
                len: cum,
            });
        }
        let genes = gspan
            .iter()
            .zip(&gplus)
            .map(|(sp, &plus)| {
                let (start, end) = sp.unwrap_or((0, 0));
                CrGene { start, end, plus }
            })
            .collect();
        let mut by_chr: Vec<(Vec<u32>, Vec<u64>)> = vec![(Vec::new(), Vec::new()); n_chr];
        let mut per: Vec<Vec<u32>> = vec![Vec::new(); n_chr];
        for i in 0..n_tx {
            if let Some(v) = per.get_mut(tr.tr_chr_idx[i]) {
                v.push(i as u32);
            }
        }
        for (c, mut v) in per.into_iter().enumerate() {
            v.sort_by_key(|&i| (tx[i as usize].start, tx[i as usize].end));
            let mut mx = Vec::with_capacity(v.len());
            let mut m = 0;
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
            genes,
            by_chr,
            tx_rank: rank(&tr.tr_ids),
            gene_rank: rank(&tr.gene_ids),
            tx_ids: tr.tr_ids.clone(),
            gene_ids: tr.gene_ids.clone(),
            gene_names: tr.gene_names.clone(),
        }
    }

    /// Transcripts with a block overlapping them.
    fn overlapping(&self, chr: usize, blocks: &[(u64, u64)]) -> Vec<u32> {
        let Some((order, mx)) = self.by_chr.get(chr) else {
            return Vec::new();
        };
        let (qs, qe) = (blocks[0].0, blocks[blocks.len() - 1].1);
        let hi = order.partition_point(|&i| self.tx[i as usize].start < qe);
        let mut out = Vec::new();
        for k in (0..hi).rev() {
            if mx[k] <= qs {
                break;
            }
            let t = &self.tx[order[k] as usize];
            if t.end > qs && blocks.iter().any(|&b| overlap(b, (t.start, t.end)) > 0) {
                out.push(order[k]);
            }
        }
        out
    }

    fn compatible(t: &CrTx, blocks: &[(u64, u64)]) -> Option<u32> {
        let ex = &t.exons;
        let mut k = ex
            .iter()
            .position(|e| e.0 <= blocks[0].0 && blocks[0].0 < e.1)?;
        let start_off = ex[k].2 as u64 + (blocks[0].0 - ex[k].0);
        for (bi, &(s, e)) in blocks.iter().enumerate() {
            let x = ex.get(k)?;
            if !(x.0 <= s && e <= x.1) {
                return None;
            }
            if bi + 1 < blocks.len() {
                if e != x.1 {
                    return None;
                }
                k += 1;
                if ex.get(k)?.0 != blocks[bi + 1].0 {
                    return None;
                }
            }
        }
        Some(start_off as u32)
    }

    /// Annotate one alignment.
    pub fn annotate(&self, t: &Transcript) -> AlnAnnot {
        let blocks = blocks_of(t);
        let mut out = AlnAnnot {
            region: CrRegion::Intergenic,
            tx_sense: Vec::new(),
            tx_anti: Vec::new(),
            gene_sense: Vec::new(),
            gene_anti: Vec::new(),
        };
        if blocks.is_empty() {
            return out;
        }
        let ovl = self.overlapping(t.chr_idx, &blocks);
        if ovl.is_empty() {
            return out;
        }
        let ref_len: u64 = blocks.iter().map(|b| b.1 - b.0).sum();
        let read_plus = !t.is_reverse;
        for &i in &ovl {
            let tx = &self.tx[i as usize];
            if let Some(off) = Self::compatible(tx, &blocks) {
                let pos = if tx.plus {
                    off
                } else {
                    tx.len - (off + ref_len as u32)
                };
                let hit = TxHit {
                    tx: i,
                    gene: tx.gene,
                    pos,
                };
                if tx.plus == read_plus {
                    out.tx_sense.push(hit);
                } else {
                    out.tx_anti.push(hit);
                }
            }
        }
        if !out.tx_sense.is_empty() || !out.tx_anti.is_empty() {
            out.region = CrRegion::Exonic;
            out.tx_sense
                .sort_unstable_by_key(|h| self.tx_rank[h.tx as usize]);
            out.tx_anti
                .sort_unstable_by_key(|h| self.tx_rank[h.tx as usize]);
            return out;
        }
        // Gene level.
        let mut genes: Vec<u32> = ovl.iter().map(|&i| self.tx[i as usize].gene).collect();
        genes.sort_unstable();
        genes.dedup();
        genes.retain(|&g| {
            let gs = &self.genes[g as usize];
            blocks
                .iter()
                .all(|&b| overlap(b, (gs.start, gs.end)) * 2 >= b.1 - b.0)
        });
        if genes.is_empty() {
            return out;
        }
        let tot: u64 = ref_len;
        let exonic = ovl.iter().any(|&i| {
            let tx = &self.tx[i as usize];
            if !genes.contains(&tx.gene) {
                return false;
            }
            let mut sum = 0u64;
            for &b in &blocks {
                let c: u64 = tx.exons.iter().map(|e| overlap(b, (e.0, e.1))).sum();
                if c == 0 {
                    return false;
                }
                sum += c;
            }
            sum * 2 >= tot
        });
        out.region = if exonic {
            CrRegion::Exonic
        } else {
            CrRegion::Intronic
        };
        for g in genes {
            if self.genes[g as usize].plus == read_plus {
                out.gene_sense.push(g);
            } else {
                out.gene_anti.push(g);
            }
        }
        out.gene_sense
            .sort_unstable_by_key(|&g| self.gene_rank[g as usize]);
        out.gene_anti
            .sort_unstable_by_key(|&g| self.gene_rank[g as usize]);
        out
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
    /// 100..200, 300..400) and `TA2` (100..200 only), and `GB` (-) with `TB`
    /// (exons 1000..1100).
    fn model() -> CrModel {
        let mk = |start, end, plus, gene, exons: Vec<(u64, u64)>| {
            let mut cum = 0;
            let ex: Vec<(u64, u64, u32)> = exons
                .iter()
                .map(|&(s, e)| {
                    let r = (s, e, cum);
                    cum += (e - s) as u32;
                    r
                })
                .collect();
            CrTx {
                start,
                end,
                plus,
                gene,
                exons: ex,
                len: cum,
            }
        };
        let tx = vec![
            mk(100, 400, true, 0, vec![(100, 200), (300, 400)]),
            mk(100, 200, true, 0, vec![(100, 200)]),
            mk(1000, 1100, false, 1, vec![(1000, 1100)]),
        ];
        let genes = vec![
            CrGene {
                start: 100,
                end: 400,
                plus: true,
            },
            CrGene {
                start: 1000,
                end: 1100,
                plus: false,
            },
        ];
        CrModel {
            by_chr: vec![(vec![0, 1, 2], vec![400, 400, 1100])],
            tx,
            genes,
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
    fn a_read_overlapping_a_gene_by_under_half_is_intergenic() {
        let m = model();
        // 60 bases, only the first 10 in GB's span.
        let a = m.annotate(&aln(950, &[(Kind::Match, 60)], false));
        assert_eq!(a.region, CrRegion::Intergenic);
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
