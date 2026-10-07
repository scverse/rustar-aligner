//! STAR's chimeric detection, ported as written.
//!
//! STAR has exactly two chimeric detectors, and both work on the same input:
//! the window transcripts `trAll[iW][iTr]` that stitching produced. With
//! `--chimMultimapNmax 0` it runs `ReadAlign::chimericDetectionOld`, which pins
//! the best transcript and looks for one partner; otherwise
//! `ChimericDetection::chimericDetectionMult`, which tries every pair and
//! re-stitches each candidate (`ChimericAlign::chimericStitching`).
//!
//! For a paired read those window transcripts are *combined-read* transcripts:
//! the read is `mate1 | spacer | RC(mate2)`, a transcript may cover both mates,
//! and its exons carry the mate they came from (`EX_iFrag`). Everything below
//! works in that frame, so a chimeric segment can be a whole stitched pair and
//! its CIGAR carries the second mate after a `p` operation, as STAR's does.
//!
//! [`WinTr`] is STAR's `Transcript` reduced to the fields these functions read.

use crate::align::score::{AlignmentScorer, SpliceMotif};
use crate::align::transcript::Transcript;
use crate::chimeric::segment::{
    ChimericAlignment, ChimericSegment, ExonSpan, JunctionLine, MultimapInfo,
};
use crate::genome::Genome;
use crate::params::Parameters;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};

/// One row of STAR's `exons[][]`: genome start (`EX_G`), read start in the
/// transcript's own strand frame (`EX_R`), length (`EX_L`) and mate (`EX_iFrag`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChimExon {
    pub g: u64,
    pub r: usize,
    pub l: usize,
    pub frag: u8,
}

/// STAR's `canonSJ` value for the gap between two consecutive mates.
const SJ_MATE_GAP: i32 = -3;
const SJ_INSERTION: i32 = -2;
const SJ_DELETION: i32 = -1;

/// The parts of STAR's `Transcript` that chimeric detection reads.
#[derive(Debug, Clone)]
pub struct WinTr {
    pub chr: usize,
    /// `Str`: 0 forward, 1 reverse. Read coordinates index `Read1[0]` for 0 and
    /// its reverse complement `Read1[2]` for 1.
    pub str_: u8,
    pub exons: Vec<ChimExon>,
    /// `canonSJ` per gap: -3 between mates, -2 insertion, -1 deletion, 0..6 a
    /// junction motif (0 non-canonical, odd `+`, even `-`).
    pub canon_sj: Vec<i32>,
    pub sj_annot: Vec<bool>,
    pub max_score: i32,
    /// `intronMotifs[sjStr]` counts: unstranded, `+`, `-` junctions.
    pub intron_motifs: [u32; 3],
    /// Sum of exon lengths, as STAR's `rLength` (not the read span).
    pub r_length: usize,
    pub ro_start: usize,
}

/// STAR's `canonSJ` code for a junction motif.
fn motif_code(m: SpliceMotif) -> i32 {
    match m {
        SpliceMotif::NonCanonical => 0,
        SpliceMotif::GtAg => 1,
        SpliceMotif::CtAc => 2,
        SpliceMotif::GcAg => 3,
        SpliceMotif::CtGc => 4,
        SpliceMotif::AtAc => 5,
        SpliceMotif::GtAt => 6,
    }
}

/// Exons and gap codes of one mate's alignment, read from its CIGAR.
///
/// STAR ends an exon at every junction, insertion and deletion, which is
/// exactly where a CIGAR switches away from `M`. `r_offset` places the mate in
/// the combined read; the leading soft clip places the first exon within it.
fn mate_blocks(t: &Transcript, r_offset: usize, frag: u8) -> (Vec<ChimExon>, Vec<i32>, Vec<bool>) {
    let mut exons = Vec::new();
    let mut gaps = Vec::new();
    let mut annot = Vec::new();
    let (mut gpos, mut rpos) = (t.genome_start, r_offset);
    let mut junction_idx = 0usize;
    // Pending gap kind since the last exon: junction beats deletion beats
    // insertion, as one STAR gap holds one `canonSJ`.
    let mut pending: Option<(i32, bool)> = None;
    for op in &t.cigar {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if !exons.is_empty() {
                    let (code, ann) = pending.take().unwrap_or((SJ_DELETION, false));
                    gaps.push(code);
                    annot.push(ann);
                }
                exons.push(ChimExon {
                    g: gpos,
                    r: rpos,
                    l: len,
                    frag,
                });
                gpos += len as u64;
                rpos += len;
            }
            Kind::SoftClip => {
                if exons.is_empty() {
                    rpos += len;
                }
            }
            Kind::Insertion => {
                rpos += len;
                if pending.is_none() {
                    pending = Some((SJ_INSERTION, false));
                }
            }
            Kind::Deletion => {
                gpos += len as u64;
                if !matches!(pending, Some((c, _)) if c >= 0) {
                    pending = Some((SJ_DELETION, false));
                }
            }
            Kind::Skip => {
                gpos += len as u64;
                let m = t
                    .junction_motifs
                    .get(junction_idx)
                    .copied()
                    .unwrap_or(SpliceMotif::NonCanonical);
                let ann = t
                    .junction_annotated
                    .get(junction_idx)
                    .copied()
                    .unwrap_or(false);
                junction_idx += 1;
                pending = Some((motif_code(m), ann));
            }
            _ => {}
        }
    }
    (exons, gaps, annot)
}

impl WinTr {
    /// A transcript covering one mate (or a single-end read).
    pub fn single(
        t: &Transcript,
        str_: u8,
        r_offset: usize,
        frag: u8,
        lread: usize,
    ) -> Option<Self> {
        let (exons, canon_sj, sj_annot) = mate_blocks(t, r_offset, frag);
        Self::build(t.chr_idx, str_, exons, canon_sj, sj_annot, t.score, lread)
    }

    /// A combined-read transcript covering both mates: `first` occupies the
    /// start of the strand frame, `second` follows the spacer.
    #[allow(clippy::too_many_arguments)]
    pub fn pair(
        first: &Transcript,
        first_frag: u8,
        second: &Transcript,
        second_frag: u8,
        second_offset: usize,
        str_: u8,
        score: i32,
        lread: usize,
    ) -> Option<Self> {
        let (mut exons, mut canon_sj, mut sj_annot) = mate_blocks(first, 0, first_frag);
        let (e2, c2, a2) = mate_blocks(second, second_offset, second_frag);
        if exons.is_empty() || e2.is_empty() {
            return None;
        }
        canon_sj.push(SJ_MATE_GAP);
        sj_annot.push(false);
        exons.extend(e2);
        canon_sj.extend(c2);
        sj_annot.extend(a2);
        Self::build(first.chr_idx, str_, exons, canon_sj, sj_annot, score, lread)
    }

    fn build(
        chr: usize,
        str_: u8,
        exons: Vec<ChimExon>,
        canon_sj: Vec<i32>,
        sj_annot: Vec<bool>,
        max_score: i32,
        lread: usize,
    ) -> Option<Self> {
        if exons.is_empty() {
            return None;
        }
        // `intronMotifs[trA.sjStr[iex]]++` over junctions (`canonSJ >= 0`). An
        // unannotated junction's `sjStr` follows its motif; an annotated one
        // takes the annotation's strand, which for a canonical motif agrees.
        // An annotated non-canonical junction's strand is not carried here, so
        // it is left uncounted rather than miscounted as unstranded.
        let mut intron_motifs = [0u32; 3];
        for (&c, &a) in canon_sj.iter().zip(&sj_annot) {
            if c < 0 {
                continue;
            }
            let s = if c == 0 {
                if a {
                    continue;
                }
                0
            } else if c % 2 == 1 {
                1
            } else {
                2
            };
            intron_motifs[s] += 1;
        }
        let r_length = exons.iter().map(|e| e.l).sum();
        let r_start = exons[0].r;
        // `stitchWindowAligns.cpp:299`.
        let ro_start = if str_ == 0 {
            r_start
        } else {
            lread.wrapping_sub(r_start).wrapping_sub(r_length)
        };
        Some(Self {
            chr,
            str_,
            exons,
            canon_sj,
            sj_annot,
            max_score: max_score.max(0),
            intron_motifs,
            r_length,
            ro_start,
        })
    }

    /// `gLength`: genomic span of the whole transcript.
    pub fn g_length(&self) -> u64 {
        let last = self.exons[self.exons.len() - 1];
        last.g + last.l as u64 - self.exons[0].g
    }

    fn last(&self) -> ChimExon {
        self.exons[self.exons.len() - 1]
    }
}

/// The read as chimeric detection sees it.
pub struct ChimRead<'a> {
    /// `Read1[0]`: the read, or `mate1 | spacer | RC(mate2)` for a pair.
    pub fwd: &'a [u8],
    /// `Read1[2]`: reverse complement of `fwd`.
    pub rev: Vec<u8>,
    /// `readLength[0..2]`; the second is 0 for single-end.
    pub read_length: [usize; 2],
    pub paired: bool,
    pub name: &'a str,
    /// The mates as sequenced, for the per-mate segments WithinBAM writes.
    pub mates: [&'a [u8]; 2],
}

impl<'a> ChimRead<'a> {
    pub fn new(
        fwd: &'a [u8],
        read_length: [usize; 2],
        name: &'a str,
        mates: [&'a [u8]; 2],
    ) -> Self {
        let rev = fwd
            .iter()
            .rev()
            .map(|&b| if b < 4 { 3 - b } else { b })
            .collect();
        Self {
            fwd,
            rev,
            read_length,
            paired: read_length[1] > 0,
            name,
            mates,
        }
    }

    /// `Lread`.
    pub fn lread(&self) -> usize {
        self.fwd.len()
    }

    /// `readLengthPairOriginal`.
    fn pair_len(&self) -> usize {
        if self.paired {
            self.read_length[0] + self.read_length[1] + 1
        } else {
            self.read_length[0]
        }
    }
}

/// STAR's `trAll`, plus what `multMapSelect` derived from it.
pub struct Windows {
    pub tr: Vec<Vec<WinTr>>,
}

impl Windows {
    /// Build each window the way `stitchWindowAligns` records into `wTr`
    /// (`stitchWindowAligns.cpp:337-381`), taking transcripts in discovery
    /// order: a transcript whose blocks another already covers is dropped if it
    /// scores lower, one that covers an earlier transcript evicts it, and the
    /// rest are inserted by score, shorter genomic span first on ties.
    ///
    /// With chimeric detection on STAR records every transcript regardless of
    /// score (`|| P.pCh.segmentMin>0`), so there is no score gate here.
    pub fn new(tr: Vec<Vec<WinTr>>, max_per_window: usize) -> Self {
        let tr = tr
            .into_iter()
            .map(|win| {
                let mut w: Vec<WinTr> = Vec::new();
                for t in win {
                    let mapped = t.r_length;
                    let mut i = 0;
                    let mut dropped = false;
                    while i < w.len() {
                        let n = blocks_overlap(&t, &w[i]);
                        let (u_new, u_old) = (mapped - n, w[i].r_length - n);
                        if u_new == 0 && t.max_score < w[i].max_score {
                            dropped = true;
                            break;
                        } else if u_old == 0 {
                            w.remove(i);
                        } else {
                            // `uOld>0 && (uNew>0 || Score>=old)` always holds here.
                            i += 1;
                        }
                    }
                    if dropped {
                        continue;
                    }
                    let at = w
                        .iter()
                        .position(|o| {
                            t.max_score > o.max_score
                                || (t.max_score == o.max_score && t.g_length() < o.g_length())
                        })
                        .unwrap_or(w.len());
                    w.insert(at, t);
                    // STAR overwrites the slot past the cap rather than growing.
                    w.truncate(max_per_window.max(1));
                }
                w
            })
            .filter(|w| !w.is_empty())
            .collect();
        Self { tr }
    }

    /// `trBest` (`ReadAlign_stitchPieces.cpp:358`): the best window's top
    /// transcript, a later window winning only on a higher score or an equal
    /// score with a shorter span.
    fn best(&self) -> Option<usize> {
        let mut best: Option<usize> = None;
        for (iw, w) in self.tr.iter().enumerate() {
            match best {
                None => best = Some(iw),
                Some(b) => {
                    let (cur, top) = (&w[0], &self.tr[b][0]);
                    if cur.max_score > top.max_score
                        || (cur.max_score == top.max_score && cur.g_length() < top.g_length())
                    {
                        best = Some(iw);
                    }
                }
            }
        }
        best
    }

    /// `multMapSelect`: `nTr` and the first two `trMult` entries, as
    /// (window, transcript) indices.
    fn mult(&self, score_range: i32) -> (usize, [Option<(usize, usize)>; 2]) {
        let max_score = self.tr.iter().map(|w| w[0].max_score).max().unwrap_or(0);
        let mut n = 0;
        let mut first = [None, None];
        for (iw, w) in self.tr.iter().enumerate() {
            for (it, t) in w.iter().enumerate() {
                if t.max_score + score_range >= max_score {
                    if n < 2 {
                        first[n] = Some((iw, it));
                    }
                    n += 1;
                }
            }
        }
        (n, first)
    }
}

/// A genomic base as STAR's `G[]` holds it; positions outside the array read
/// as padding, which like `N` is `> 3`.
fn gbase(genome: &Genome, pos: i64) -> u8 {
    if pos < 0 {
        return 5;
    }
    genome.get_base(pos as u64).unwrap_or(5)
}

fn comp(b: u8) -> u8 {
    if b < 4 { 3 - b } else { b }
}

/// `Transcript::alignScore` (`Transcript_alignScore.cpp`): re-score a
/// transcript from its exons.
fn align_score(tr: &WinTr, read: &ChimRead, genome: &Genome, scorer: &AlignmentScorer) -> i32 {
    let r: &[u8] = if tr.str_ == 0 { read.fwd } else { &read.rev };
    let mut score = 0i32;
    for e in &tr.exons {
        for ii in 0..e.l {
            let r1 = r.get(e.r + ii).copied().unwrap_or(5);
            let g1 = gbase(genome, (e.g + ii as u64) as i64);
            if r1 > 3 || g1 > 3 {
            } else if r1 == g1 {
                score += 1;
            } else {
                score -= 1;
            }
        }
    }
    for iex in 0..tr.exons.len().saturating_sub(1) {
        let (a, b) = (tr.exons[iex], tr.exons[iex + 1]);
        if tr.sj_annot[iex] {
            score += scorer.sjdb_score;
            continue;
        }
        score += match tr.canon_sj[iex] {
            SJ_MATE_GAP => 0,
            SJ_INSERTION => {
                (b.r - a.r - a.l) as i32 * scorer.score_ins_base + scorer.score_ins_open
            }
            SJ_DELETION => {
                (b.g - a.g - a.l as u64) as i32 * scorer.score_del_base + scorer.score_del_open
            }
            0 => scorer.score_gap_noncan + scorer.score_gap,
            1 | 2 => scorer.score_gap,
            3 | 4 => scorer.score_gap_gcag + scorer.score_gap,
            5 | 6 => scorer.score_gap_atac + scorer.score_gap,
            _ => 0,
        };
    }
    if scorer.score_genomic_length_log2_scale != 0.0 {
        let span = tr.g_length().max(1) as f64;
        score += (span.log2() * scorer.score_genomic_length_log2_scale - 0.5).ceil() as i32;
    }
    score
}

/// `Transcript::generateCigarP` / `ReadAlign::outputTranscriptCIGARp`: the
/// CIGAR with a `p` operation for the gap between mates, negative when the
/// mates overlap.
fn cigar_p(tr: &WinTr, read: &ChimRead) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    let left = if read.paired { tr.str_ as usize } else { 0 };
    let rl = read.read_length[left];
    let e0 = tr.exons[0];
    let trim_l = e0.r - if e0.r < rl { 0 } else { rl + 1 };
    if trim_l > 0 {
        let _ = write!(s, "{trim_l}S");
    }
    for ii in 0..tr.exons.len() {
        if ii > 0 {
            let (p, c) = (tr.exons[ii - 1], tr.exons[ii]);
            let p_end = p.g + p.l as u64;
            if c.g >= p_end {
                let gap_g = c.g - p_end;
                if tr.canon_sj[ii - 1] == SJ_MATE_GAP {
                    let s1 = rl - (p.r + p.l);
                    let s2 = c.r - (rl + 1);
                    if s1 > 0 {
                        let _ = write!(s, "{s1}S");
                    }
                    let _ = write!(s, "{gap_g}p");
                    if s2 > 0 {
                        let _ = write!(s, "{s2}S");
                    }
                } else {
                    let gap_r = c.r - p.r - p.l;
                    if gap_r > 0 {
                        let _ = write!(s, "{gap_r}I");
                    }
                    if tr.canon_sj[ii - 1] >= 0 || tr.sj_annot[ii - 1] {
                        let _ = write!(s, "{gap_g}N");
                    } else if gap_g > 0 {
                        let _ = write!(s, "{gap_g}D");
                    }
                }
            } else {
                // Overlapping mates: STAR writes only the overlap, no clips.
                let _ = write!(s, "-{}p", p_end - c.g);
            }
        }
        let _ = write!(s, "{}M", tr.exons[ii].l);
    }
    let last = tr.last();
    let end = if last.r < rl { rl } else { read.pair_len() };
    let trim_r = end as i64 - (last.r + last.l) as i64;
    if trim_r > 0 {
        let _ = write!(s, "{trim_r}S");
    }
    s
}

/// STAR's chimeric strand code for a transcript in `chimericDetectionOld`:
/// 0 undefined, 1 same as the RNA, 2 opposite.
fn chim_str_old(t: &WinTr) -> u8 {
    if t.intron_motifs[1] == 0 && t.intron_motifs[2] == 0 {
        0
    } else if (t.str_ == 0) == (t.intron_motifs[1] > 0) {
        1
    } else {
        2
    }
}

/// `ChimericSegment::str`, which also calls a transcript with both motif
/// strands undefined.
fn chim_str_seg(t: &WinTr) -> u8 {
    let (p, m) = (t.intron_motifs[1], t.intron_motifs[2]);
    if (p == 0 && m == 0) || (p > 0 && m > 0) {
        0
    } else if (t.str_ == 0) == (p > 0) {
        1
    } else {
        2
    }
}

/// `roS`/`roE` of a segment, with the spacer removed for the second mate.
fn ro_se(t: &WinTr, read: &ChimRead) -> (usize, usize) {
    let lread = read.lread();
    let (first, last) = (t.exons[0], t.last());
    let mut start = if t.str_ == 0 {
        first.r
    } else {
        lread - last.r - last.l
    };
    let mut end = if t.str_ == 0 {
        last.r + last.l - 1
    } else {
        lread - first.r - 1
    };
    if start > read.read_length[0] {
        start -= 1;
    }
    if end > read.read_length[0] {
        end -= 1;
    }
    (start, end)
}

fn overlap_ro(s1: usize, e1: usize, s2: usize, e2: usize) -> usize {
    if s2 > s1 {
        if s2 > e1 { 0 } else { e1 - s2 + 1 }
    } else if e2 < s1 {
        0
    } else {
        e2 - s1 + 1
    }
}

/// `blocksOverlap` (`blocksOverlap.cpp`): read bases two transcripts place on
/// the same diagonal. Like STAR's, it does not look at chromosome or strand.
fn blocks_overlap(t1: &WinTr, t2: &WinTr) -> usize {
    let (mut i1, mut i2, mut n_overlap) = (0usize, 0usize, 0usize);
    while i1 < t1.exons.len() && i2 < t2.exons.len() {
        let (x, y) = (t1.exons[i1], t2.exons[i2]);
        let (rs1, rs2) = (x.r, y.r);
        let (re1, re2) = (x.r + x.l, y.r + y.l);
        if rs1 >= re2 {
            i2 += 1;
        } else if rs2 >= re1 {
            i1 += 1;
        } else {
            if x.g.wrapping_sub(rs1 as u64) == y.g.wrapping_sub(rs2 as u64) {
                n_overlap += re1.min(re2) - rs1.max(rs2);
            }
            if re1 >= re2 {
                i2 += 1;
            }
            if re2 >= re1 {
                i1 += 1;
            }
        }
    }
    n_overlap
}

/// A placed chimeric junction: STAR's `trChim[0..2]` (or `al1`/`al2`) after
/// the junction shift, with `chimJ*`, `chimMotif` and `chimRepeat*`.
struct Placed {
    t0: WinTr,
    t1: WinTr,
    j0: u64,
    j1: u64,
    motif: i32,
    rep0: u32,
    rep1: u32,
}

/// Junction placement shared by both detectors: the bracketing branch, or the
/// scan for the best junction position within a mate followed by the shift
/// and the repeat measurement (`ReadAlign_chimericDetectionOld.cpp:121-307`,
/// `ChimericAlign_chimericStitching.cpp:21-168`). The two copies in STAR
/// differ only in how they report failure, which the caller handles.
///
/// Returns `None` when an `N` in the read, or in the genome under
/// `banGenomicN`, falls inside the scanned span.
fn place(
    mut t0: WinTr,
    mut t1: WinTr,
    chim_str: u8,
    read: &ChimRead,
    genome: &Genome,
    params: &Parameters,
) -> Option<Placed> {
    let e0 = if t0.str_ == 1 { 0 } else { t0.exons.len() - 1 };
    let e1 = if t1.str_ == 0 { 0 } else { t1.exons.len() - 1 };
    if t0.exons[e0].frag < t1.exons[e1].frag {
        // Mates bracket the junction.
        let x0 = t0.exons[e0];
        let x1 = t1.exons[e1];
        let j0 = if t0.str_ == 1 {
            x0.g.wrapping_sub(1)
        } else {
            x0.g + x0.l as u64
        };
        let j1 = if t1.str_ == 0 {
            x1.g.wrapping_sub(1)
        } else {
            x1.g + x1.l as u64
        };
        return Some(Placed {
            t0,
            t1,
            j0,
            j1,
            motif: -1,
            rep0: 0,
            rep1: 0,
        });
    }

    let ban_n = params.chim_filter.iter().any(|f| f == "banGenomicN");
    let lread = read.lread() as i64;
    let x0 = t0.exons[e0];
    let x1 = t1.exons[e1];
    let ro_start0 = if t0.str_ == 0 {
        x0.r as i64
    } else {
        lread - x0.r as i64 - x0.l as i64
    };
    let ro_start1 = if t1.str_ == 0 {
        x1.r as i64
    } else {
        lread - x1.r as i64 - x1.l as i64
    };
    let (g0, l0) = (x0.g as i64, x0.l as i64);
    let (g1, l1) = (x1.g as i64, x1.l as i64);
    let base0 = |off: i64| -> u8 {
        if t0.str_ == 0 {
            gbase(genome, g0 + off)
        } else {
            comp(gbase(genome, g0 + l0 - 1 - off))
        }
    };
    // `off` is the position in the leading segment's read frame; the trailing
    // segment is read in that same frame, shifted by `roStart0 - roStart1`.
    let base1 = |off: i64| -> u8 {
        if t1.str_ == 0 {
            gbase(genome, g1 - ro_start1 + ro_start0 + off)
        } else {
            comp(gbase(genome, g1 + l1 - 1 + ro_start1 - ro_start0 - off))
        }
    };

    let jr_max = {
        let m = ro_start1 + l1;
        if m > ro_start0 { m - ro_start0 - 1 } else { 0 }
    };
    let (mut motif, mut jr_best) = (0i32, 0i64);
    let (mut j_score, mut j_score_best) = (0i32, -999_999i32);
    let mut jr = 0i64;
    while jr < jr_max {
        // STAR compares the scan offset, not the read position, with the mate
        // boundary when it skips the spacer. Kept as written.
        if jr == read.read_length[0] as i64 {
            jr += 1;
        }
        let br = read
            .fwd
            .get((ro_start0 + jr) as usize)
            .copied()
            .unwrap_or(5);
        let b0 = base0(jr);
        let b1 = base1(jr);
        if (ban_n && (b0 > 3 || b1 > 3)) || br > 3 {
            return None;
        }
        let (b01, b02) = (base0(jr + 1), base0(jr + 2));
        let (b11, b12) = (base1(jr - 1), base1(jr));
        let mut j_motif = 0;
        if b01 == 2 && b02 == 3 && b11 == 0 && b12 == 2 {
            if chim_str != 2 {
                j_motif = 1;
            }
        } else if b01 == 1 && b02 == 3 && b11 == 0 && b12 == 1 && chim_str != 1 {
            j_motif = 2;
        }
        if br == b0 && br != b1 {
            j_score += 1;
        } else if br != b0 && br == b1 {
            j_score -= 1;
        }
        let j_score_j = if j_motif == 0 {
            j_score + params.chim_score_junction_non_gtag
        } else {
            j_score
        };
        if j_score_j > j_score_best || (j_score_j == j_score_best && j_motif > 0) {
            motif = j_motif;
            jr_best = jr;
            j_score_best = j_score_j;
        }
        jr += 1;
    }

    // Shift the junction.
    let jr_best = jr_best as usize;
    let (j0, j1);
    {
        let x = &mut t0.exons[e0];
        if t0.str_ == 1 {
            let d = x.l - jr_best - 1;
            x.r += d;
            x.g += d as u64;
            x.l = jr_best + 1;
            j0 = x.g.wrapping_sub(1);
        } else {
            x.l = jr_best + 1;
            j0 = x.g + x.l as u64;
        }
    }
    {
        let x = &mut t1.exons[e1];
        let new_l = (ro_start1 + x.l as i64 - ro_start0 - jr_best as i64 - 1) as usize;
        if t1.str_ == 0 {
            let d = (ro_start0 + jr_best as i64 + 1 - ro_start1) as usize;
            x.r += d;
            x.g += d as u64;
            x.l = new_l;
            j1 = x.g.wrapping_sub(1);
        } else {
            x.l = new_l;
            j1 = x.g + x.l as u64;
        }
    }

    // Repeats around the junction.
    let side0 = |off: i64| -> u8 {
        if t0.str_ == 0 {
            gbase(genome, j0 as i64 + off)
        } else {
            comp(gbase(genome, j0 as i64 - off))
        }
    };
    let side1 = |off: i64| -> u8 {
        if t1.str_ == 0 {
            gbase(genome, j1 as i64 + off)
        } else {
            comp(gbase(genome, j1 as i64 - off))
        }
    };
    let mut rep1 = 0u32;
    while rep1 < 100 && side0(rep1 as i64) == side1(rep1 as i64 + 1) {
        rep1 += 1;
    }
    let mut rep0 = 0u32;
    while rep0 < 100 && side0(-1 - rep0 as i64) == side1(-(rep0 as i64)) {
        rep0 += 1;
    }

    Some(Placed {
        t0,
        t1,
        j0,
        j1,
        motif,
        rep0,
        rep1,
    })
}

/// `ReadAlign::chimericDetectionOld`.
fn detect_old(
    w: &Windows,
    read: &ChimRead,
    genome: &Genome,
    params: &Parameters,
) -> Option<Placed> {
    let lread = read.lread();
    let seg_min = params.chim_segment_min as usize;
    let ib = w.best()?;
    let tr_best = &w.tr[ib][0];
    let (n_tr, tr_mult) = w.mult(params.out_filter_multimap_score_range);
    let main_mult = params.chim_main_segment_mult_nmax as usize;

    if n_tr > main_mult && n_tr != 2 {
        return None;
    }
    let lb = tr_best.last();
    if !(seg_min > 0
        && tr_best.r_length >= seg_min
        && (lb.r + lb.l + seg_min <= lread || tr_best.exons[0].r >= seg_min)
        && tr_best.intron_motifs[0] == 0
        && (tr_best.intron_motifs[1] == 0 || tr_best.intron_motifs[2] == 0))
    {
        return None;
    }

    let (mut chim_score_best, mut chim_score_next) = (0i32, 0i32);
    let mut tr_chim1: Option<(usize, usize)> = None;
    let (ro_start1, ro_end1) = ro_se(tr_best, read);
    let mut chim_str = chim_str_old(tr_best);
    let mut chim_str_best = 0u8;
    let gap_max = params.chim_segment_read_gap_max as usize;

    for (iw, win) in w.tr.iter().enumerate() {
        for (iwt, t) in win.iter().enumerate() {
            if iw != ib && iwt > 0 {
                break;
            }
            if iw == ib && iwt == 0 {
                continue;
            }
            if t.intron_motifs[0] > 0 {
                continue;
            }
            let chim_str1 = chim_str_old(t);
            if chim_str != 0 && chim_str1 != 0 && chim_str != chim_str1 {
                continue;
            }
            let (ro_start2, ro_end2) = ro_se(t, read);
            let overlap = overlap_ro(ro_start1, ro_end1, ro_start2, ro_end2);
            let rl0 = read.read_length[0];
            let diff_mates =
                (ro_end1 < rl0 && ro_start2 >= rl0) || (ro_end2 < rl0 && ro_start1 >= rl0);
            if ro_end1 > seg_min + ro_start1 + overlap
                && ro_end2 > seg_min + ro_start2 + overlap
                && (diff_mates
                    || (ro_end1 + gap_max + 1 >= ro_start2 && ro_end2 + gap_max + 1 >= ro_start1))
            {
                let chim_score = tr_best.max_score + t.max_score - overlap as i32;
                let mut overlap1 = 0;
                if iwt > 0
                    && chim_score_best > 0
                    && let Some((a, b)) = tr_chim1
                {
                    overlap1 = blocks_overlap(&w.tr[a][b], t);
                }
                if chim_score > chim_score_best {
                    tr_chim1 = Some((iw, iwt));
                    if overlap1 == 0 {
                        chim_score_next = chim_score_best;
                    }
                    chim_score_best = chim_score;
                    chim_str_best = chim_str1;
                } else if chim_score > chim_score_next && overlap1 == 0 {
                    chim_score_next = chim_score;
                }
            }
        }
    }

    let frag_len = (read.read_length[0] + read.read_length[1]) as i32;
    if !(chim_score_best >= params.chim_score_min
        && chim_score_best + params.chim_score_drop_max >= frag_len)
    {
        return None;
    }
    let tc1 = tr_chim1?;
    if n_tr > main_mult && Some(tc1) != tr_mult[0] && Some(tc1) != tr_mult[1] {
        return None;
    }
    if chim_str == 0 {
        chim_str = chim_str_best;
    }
    if chim_score_next + params.chim_score_separation >= chim_score_best {
        return None;
    }

    let (mut c0, mut c1) = (tr_best.clone(), w.tr[tc1.0][tc1.1].clone());
    if c0.ro_start > c1.ro_start {
        std::mem::swap(&mut c0, &mut c1);
    }
    let e0 = if c0.str_ == 1 { 0 } else { c0.exons.len() - 1 };
    let e1 = if c1.str_ == 0 { 0 } else { c1.exons.len() - 1 };
    if c0.exons[e0].frag > c1.exons[e1].frag {
        return None;
    }
    let overhang_min = params.chim_junction_overhang_min as usize;
    if c0.exons[e0].frag == c1.exons[e1].frag
        && !(c0.exons[e0].l >= overhang_min && c1.exons[e1].l >= overhang_min)
    {
        return None;
    }

    let p = place(c0, c1, chim_str, read, genome, params)?;
    if p.motif == 0 {
        chim_score_best += 1 + params.chim_score_junction_non_gtag;
        if !(chim_score_best >= params.chim_score_min
            && chim_score_best + params.chim_score_drop_max >= frag_len)
        {
            return None;
        }
    }

    // Final check: different chromosome or strand, or farther apart than a
    // linear alignment may reach. STAR's distance is unsigned and wraps when
    // the junction runs backwards, which then always counts as far.
    let dist = if p.t0.str_ == 0 {
        p.j1.wrapping_sub(p.j0).wrapping_add(1)
    } else {
        p.j0.wrapping_sub(p.j1).wrapping_add(1)
    };
    let limit = if p.motif >= 0 {
        params.align_intron_max as u64
    } else {
        params.align_mates_gap_max as u64
    };
    if p.t0.str_ != p.t1.str_ || p.t0.chr != p.t1.chr || dist > limit {
        let (x0, x1) = (p.t0.exons[e0], p.t1.exons[e1]);
        if p.motif >= 0
            && (x0.l < overhang_min + p.rep0 as usize || x1.l < overhang_min + p.rep1 as usize)
        {
            return None;
        }
        return Some(p);
    }
    None
}

/// A chimera from `chimericDetectionMult`, with its re-stitched score.
struct MultChim {
    p: Placed,
    score: i32,
}

/// `ChimericDetection::chimericDetectionMult` with `ChimericAlign`'s checks
/// and stitching.
fn detect_mult(
    w: &Windows,
    max_non_chim: i32,
    read: &ChimRead,
    genome: &Genome,
    scorer: &AlignmentScorer,
    params: &Parameters,
) -> (Vec<MultChim>, i32, i32) {
    let seg_min = params.chim_segment_min as usize;
    let gap_max = params.chim_segment_read_gap_max as usize;
    let overhang_min = params.chim_junction_overhang_min as usize;
    let max_possible = (read.read_length[0] + read.read_length[1]) as i32;
    let mut min_score = params.chim_score_min;
    if max_non_chim >= min_score {
        min_score = max_non_chim + 1;
    }
    if max_possible - params.chim_score_drop_max > min_score {
        min_score = max_possible - params.chim_score_drop_max;
    }
    let rl0 = read.read_length[0];
    let seg_ok = |t: &WinTr| t.r_length >= seg_min && t.intron_motifs[0] == 0;

    let mut out: Vec<MultChim> = Vec::new();
    let mut best = 0i32;
    let flat: Vec<(usize, usize)> =
        w.tr.iter()
            .enumerate()
            .flat_map(|(iw, win)| (0..win.len()).map(move |it| (iw, it)))
            .collect();
    for (i, &(iw1, ia1)) in flat.iter().enumerate() {
        let s1 = &w.tr[iw1][ia1];
        if !seg_ok(s1) {
            continue;
        }
        let (ros1, roe1) = ro_se(s1, read);
        let str1 = chim_str_seg(s1);
        // Same order as STAR's nested window/align loops, with the second
        // segment always after the first.
        for &(iw2, ia2) in &flat[i + 1..] {
            let s2 = &w.tr[iw2][ia2];
            if !seg_ok(s2) {
                continue;
            }
            let str2 = chim_str_seg(s2);
            if str1 != 0 && str2 != 0 && str2 != str1 {
                continue;
            }
            // `chimericAlignScore`.
            let (ros2, roe2) = ro_se(s2, read);
            let overlap = if ros2 > ros1 {
                if ros2 > roe1 { 0 } else { roe1 - ros2 + 1 }
            } else if roe2 < ros1 {
                0
            } else {
                roe2 - ros1 + 1
            };
            let diff_mates = (roe1 < rl0 && ros2 >= rl0) || (roe2 < rl0 && ros1 >= rl0);
            let mut chim_score = 0;
            if roe1 > seg_min + ros1 + overlap
                && roe2 > seg_min + ros2 + overlap
                && (diff_mates || (roe1 + gap_max + 1 >= ros2 && roe2 + gap_max + 1 >= ros1))
            {
                chim_score = s1.max_score + s2.max_score - overlap as i32;
            }
            if chim_score < min_score {
                continue;
            }
            // `ChimericAlign` orders its pair by `roStart`.
            let (a1, a2) = if s1.ro_start > s2.ro_start {
                (s2, s1)
            } else {
                (s1, s2)
            };
            let ex1 = if a1.str_ == 1 { 0 } else { a1.exons.len() - 1 };
            let ex2 = if a2.str_ == 0 { 0 } else { a2.exons.len() - 1 };
            // `chimericCheck`.
            let (f1, f2) = (a1.exons[ex1].frag, a2.exons[ex2].frag);
            if f1 > f2
                || !(f1 < f2
                    || (a1.exons[ex1].l >= overhang_min && a2.exons[ex2].l >= overhang_min))
            {
                continue;
            }
            // `chimericStitching`.
            let chim_str = str1.max(str2);
            let Some(p) = place(a1.clone(), a2.clone(), chim_str, read, genome, params) else {
                continue;
            };
            // Unlike the old path, the repeat lengths are not added here
            // (`ChimericAlign_chimericStitching.cpp:172`).
            let score = if p.motif >= 0
                && (p.t0.exons[ex1].l < overhang_min || p.t1.exons[ex2].l < overhang_min)
            {
                0
            } else {
                align_score(&p.t0, read, genome, scorer)
                    + align_score(&p.t1, read, genome, scorer)
                    + if p.motif == 0 {
                        params.chim_score_junction_non_gtag
                    } else {
                        0
                    }
            };
            if score >= min_score {
                if score > best {
                    best = score;
                    if best - params.chim_multimap_score_range > min_score {
                        min_score = best - params.chim_multimap_score_range;
                    }
                }
                out.push(MultChim { p, score });
            }
        }
    }
    (out, best, min_score)
}

/// The per-mate view of a chimeric transcript that WithinBAM records use: the
/// exons of the mate holding the junction, rebased into that mate's own read.
fn mate_segment(tr: &WinTr, junction_exon: usize, read: &ChimRead, score: i32) -> ChimericSegment {
    let frag = tr.exons[junction_exon].frag;
    let ex: Vec<ChimExon> = tr
        .exons
        .iter()
        .copied()
        .filter(|e| e.frag == frag)
        .collect();
    // The mate at the start of the strand frame begins at 0, the other after
    // the spacer. Forward frames hold mate1 first, reverse frames mate2.
    let first_frag = if read.paired { tr.str_ } else { 0 };
    let mate_len = read.read_length[frag as usize];
    let offset = if frag == first_frag {
        0
    } else {
        read.read_length[first_frag as usize] + 1
    };
    // A mate's own orientation: mate1 is reversed in a reverse transcript,
    // mate2 in a forward one (it enters the combined read complemented).
    let is_reverse = if frag == 0 {
        tr.str_ == 1
    } else {
        tr.str_ == 0
    };
    let span = |e: &ChimExon| ExonSpan {
        genome_start: e.g,
        genome_end: e.g + e.l as u64,
        read_start: e.r - offset,
        read_end: e.r - offset + e.l,
    };
    let mut cigar = Vec::new();
    let lead = ex[0].r - offset;
    if lead > 0 {
        cigar.push(Op::new(Kind::SoftClip, lead));
    }
    for (i, e) in ex.iter().enumerate() {
        if i > 0 {
            let p = ex[i - 1];
            let gap_r = e.r - p.r - p.l;
            let gap_g = e.g.saturating_sub(p.g + p.l as u64) as usize;
            if gap_r > 0 {
                cigar.push(Op::new(Kind::Insertion, gap_r));
            }
            if gap_g > 0 {
                let idx = tr.exons.iter().position(|x| *x == p).unwrap_or(0);
                let kind = if tr.canon_sj.get(idx).is_some_and(|&c| c >= 0)
                    || tr.sj_annot.get(idx).copied().unwrap_or(false)
                {
                    Kind::Skip
                } else {
                    Kind::Deletion
                };
                cigar.push(Op::new(kind, gap_g));
            }
        }
        cigar.push(Op::new(Kind::Match, e.l));
    }
    let last = ex[ex.len() - 1];
    let trail = mate_len.saturating_sub(last.r - offset + last.l);
    if trail > 0 {
        cigar.push(Op::new(Kind::SoftClip, trail));
    }
    ChimericSegment {
        chr_idx: tr.chr,
        genome_start: ex[0].g,
        genome_end: last.g + last.l as u64,
        is_reverse,
        read_start: ex[0].r - offset,
        read_end: last.r - offset + last.l,
        cigar,
        score,
        n_mismatch: 0,
        first_exon: span(&ex[0]),
        last_exon: span(&last),
    }
}

fn to_alignment(
    p: &Placed,
    score: i32,
    read: &ChimRead,
    genome: &Genome,
    scorer: &AlignmentScorer,
) -> ChimericAlignment {
    let e0 = if p.t0.str_ == 1 {
        0
    } else {
        p.t0.exons.len() - 1
    };
    let e1 = if p.t1.str_ == 0 {
        0
    } else {
        p.t1.exons.len() - 1
    };
    let donor = mate_segment(&p.t0, e0, read, align_score(&p.t0, read, genome, scorer));
    let acceptor = mate_segment(&p.t1, e1, read, align_score(&p.t1, read, genome, scorer));
    let mate = p.t0.exons[e0].frag as usize;
    let mut chim = ChimericAlignment::new(
        donor,
        acceptor,
        p.motif,
        p.rep0,
        p.rep1,
        read.mates[mate].to_vec(),
        read.name.to_string(),
    );
    chim.total_score = score;
    chim.junction_line = Some(JunctionLine {
        donor_chr: p.t0.chr,
        donor_break: p.j0,
        donor_reverse: p.t0.str_ == 1,
        acceptor_chr: p.t1.chr,
        acceptor_break: p.j1,
        acceptor_reverse: p.t1.str_ == 1,
        donor_start: p.t0.exons[0].g,
        donor_cigar: cigar_p(&p.t0, read),
        acceptor_start: p.t1.exons[0].g,
        acceptor_cigar: cigar_p(&p.t1, read),
    });
    chim
}

/// `ReadAlign::chimericDetection`: run whichever STAR detector the parameters
/// select over a read's window transcripts.
pub fn chimeric_detection(
    windows: Vec<Vec<WinTr>>,
    read: &ChimRead,
    genome: &Genome,
    scorer: &AlignmentScorer,
    params: &Parameters,
) -> Vec<ChimericAlignment> {
    if params.chim_segment_min == 0 {
        return Vec::new();
    }
    let w = Windows::new(windows, params.align_transcripts_per_window_nmax);
    if params.chim_multimap_nmax == 0 {
        let Some(p) = detect_old(&w, read, genome, params) else {
            return Vec::new();
        };
        // The old path's score is not written anywhere; WithinBAM and the
        // single-best bookkeeping use the re-scored segments.
        let score =
            align_score(&p.t0, read, genome, scorer) + align_score(&p.t1, read, genome, scorer);
        return vec![to_alignment(&p, score, read, genome, scorer)];
    }
    let Some(ib) = w.best() else {
        return Vec::new();
    };
    let max_non_chim = w.tr[ib][0].max_score;
    let frag_len = (read.read_length[0] + read.read_length[1]) as i32;
    if max_non_chim > frag_len - params.chim_nonchim_score_drop_min {
        return Vec::new();
    }
    let (chims, best, min_score) = detect_mult(&w, max_non_chim, read, genome, scorer, params);
    if best == 0 {
        return Vec::new();
    }
    let kept: Vec<&MultChim> = chims.iter().filter(|c| c.score >= min_score).collect();
    if kept.len() > params.chim_multimap_nmax {
        return Vec::new();
    }
    let chim_n = kept.len();
    kept.into_iter()
        .map(|c| {
            to_alignment(&c.p, c.score, read, genome, scorer).with_multimap(MultimapInfo {
                chim_n,
                max_possible_score: frag_len,
                max_non_chim_score: max_non_chim,
                chim_score: c.score,
                best_chim_score: best,
                pe_merged: false,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::Exon;
    use crate::genome::Genome;

    const CHR_PAD: u64 = 1024;

    /// Two 1000-base chromosomes of pseudo-random sequence, so that a read
    /// taken from one place matches nowhere else by accident.
    fn genome() -> Genome {
        let n_genome = CHR_PAD * 2;
        let mut seq = vec![5u8; 2 * n_genome as usize];
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        for chr in 0..2u64 {
            for i in 0..1000u64 {
                x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                seq[(chr * CHR_PAD + i) as usize] = ((x >> 33) % 4) as u8;
            }
        }
        Genome {
            transform_blocks: None,
            sequence: seq.into(),
            n_genome,
            n_genome_real: n_genome,
            n_chr_real: 2,
            chr_name: vec!["chr0".into(), "chr1".into()],
            chr_length: vec![1000, 1000],
            chr_start: vec![0, CHR_PAD, n_genome],
        }
    }

    fn slice(gn: &Genome, start: u64, len: usize) -> Vec<u8> {
        (0..len as u64)
            .map(|i| gn.get_base(start + i).unwrap())
            .collect()
    }

    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| comp(b)).collect()
    }

    fn params(extra: &[&str]) -> Parameters {
        let mut args = vec![
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--chimSegmentMin",
            "12",
        ];
        args.extend_from_slice(extra);
        Parameters::parse_from(args)
    }

    /// A one-block transcript at `g_start` with the given CIGAR.
    fn tx(chr: usize, g_start: u64, cigar: &[(Kind, usize)], score: i32) -> Transcript {
        let ops: Vec<Op> = cigar.iter().map(|&(k, l)| Op::new(k, l)).collect();
        let g_len: usize = cigar
            .iter()
            .filter(|(k, _)| k.consumes_reference())
            .map(|&(_, l)| l)
            .sum();
        Transcript {
            chr_idx: chr,
            genome_start: g_start,
            genome_end: g_start + g_len as u64,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: g_start,
                genome_end: g_start + g_len as u64,
                read_start: 0,
                read_end: g_len,
                i_frag: 0,
            }],
            cigar: ops,
            score,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
        }
    }

    fn exon(gn: u64, r: usize, l: usize, frag: u8) -> ChimExon {
        ChimExon { g: gn, r, l, frag }
    }

    fn wintr(exons: Vec<ChimExon>, canon_sj: Vec<i32>, score: i32) -> WinTr {
        let n = canon_sj.len();
        WinTr::build(0, 0, exons, canon_sj, vec![false; n], score, 1000).unwrap()
    }

    #[test]
    fn cigar_blocks_become_star_exons_with_gap_codes() {
        // 5S10M2I10M3D10M100N10M5S: an insertion, a deletion and a GT/AG
        // junction each end an exon, and each gap gets its own `canonSJ`.
        let t = tx(
            0,
            1000,
            &[
                (Kind::SoftClip, 5),
                (Kind::Match, 10),
                (Kind::Insertion, 2),
                (Kind::Match, 10),
                (Kind::Deletion, 3),
                (Kind::Match, 10),
                (Kind::Skip, 100),
                (Kind::Match, 10),
                (Kind::SoftClip, 5),
            ],
            40,
        );
        let (exons, gaps, _) = mate_blocks(&t, 0, 0);
        assert_eq!(
            exons,
            vec![
                exon(1000, 5, 10, 0),
                exon(1010, 17, 10, 0),
                exon(1023, 27, 10, 0),
                exon(1133, 37, 10, 0),
            ]
        );
        assert_eq!(gaps, vec![SJ_INSERTION, SJ_DELETION, 1]);
    }

    #[test]
    fn cigar_p_writes_the_mate_gap_and_its_clips() {
        let fwd = vec![0u8; 21];
        let read = ChimRead::new(&fwd, [10, 10], "r", [&[], &[]]);
        // Mate1 clipped by 2 at the start, mate2 by 3 at the end, 22 bases apart.
        let t = wintr(vec![exon(100, 2, 8, 0), exon(130, 11, 7, 1)], vec![-3], 10);
        assert_eq!(cigar_p(&t, &read), "2S8M22p7M3S");
        // A trailing clip on mate1 is `s1`, written before the `prm`.
        let t = wintr(vec![exon(100, 2, 6, 0), exon(130, 11, 7, 1)], vec![-3], 10);
        assert_eq!(cigar_p(&t, &read), "2S6M2S24p7M3S");
        // Overlapping mates: the overlap only, and STAR drops both clips there.
        let t = wintr(vec![exon(100, 2, 6, 0), exon(103, 13, 7, 1)], vec![-3], 10);
        assert_eq!(cigar_p(&t, &read), "2S6M-3p7M1S");
    }

    #[test]
    fn align_score_is_star_rescoring() {
        let gn = genome();
        let mut fwd = slice(&gn, 100, 60);
        fwd[30] = comp(fwd[30]); // one mismatch
        let read = ChimRead::new(&fwd, [60, 0], "r", [&fwd, &[]]);
        let scorer = AlignmentScorer::from_params_minimal();
        let t = wintr(vec![exon(100, 0, 60, 0)], vec![], 0);
        // 59 matches - 1 mismatch, then ceil(log2(60) * -0.25 - 0.5) = -1.
        assert_eq!(align_score(&t, &read, &gn, &scorer), 57);
        // Split around a 40-base deletion: -2 open -2*40 per base (STAR
        // defaults), and the genomic span grows to 100.
        let t = wintr(
            vec![exon(100, 0, 30, 0), exon(170, 30, 30, 0)],
            vec![SJ_DELETION],
            0,
        );
        let expected = (30 - 2 * 40 - 2) // first half matches, deletion
            + slice(&gn, 170, 30)
                .iter()
                .zip(&fwd[30..])
                .map(|(a, b)| if a == b { 1 } else { -1 })
                .sum::<i32>()
            + ((100f64).log2() * -0.25 - 0.5).ceil() as i32;
        assert_eq!(align_score(&t, &read, &gn, &scorer), expected);
    }

    #[test]
    fn window_dedup_follows_stitch_window_aligns() {
        let w = |e: Vec<ChimExon>, s: i32| wintr(e, vec![], s);
        let a = w(vec![exon(100, 0, 50, 0)], 50);
        // Covered by `a` and scoring lower: never recorded.
        let sub = w(vec![exon(110, 10, 30, 0)], 30);
        // Unrelated diagonal: kept and ordered by score.
        let other = w(vec![exon(500, 0, 40, 0)], 40);
        let win = Windows::new(vec![vec![a.clone(), sub, other]], 100);
        let scores: Vec<i32> = win.tr[0].iter().map(|t| t.max_score).collect();
        assert_eq!(scores, vec![50, 40]);

        // A later transcript that covers an earlier one evicts it, even when
        // the earlier one scored higher (`uOld==0` has no score test).
        let cover = w(vec![exon(100, 0, 60, 0)], 45);
        let win = Windows::new(vec![vec![a, cover]], 100);
        assert_eq!(win.tr[0].len(), 1);
        assert_eq!(win.tr[0][0].max_score, 45);

        // Equal scores: the shorter genomic span goes first.
        let long = w(vec![exon(100, 0, 20, 0), exon(300, 20, 20, 0)], 40);
        let short = w(vec![exon(600, 0, 40, 0)], 40);
        let mut long = long;
        long.canon_sj = vec![1];
        long.sj_annot = vec![false];
        let win = Windows::new(vec![vec![long, short.clone()]], 100);
        assert_eq!(win.tr[0][0].exons, short.exons);
    }

    /// A single-end read, half chr0 and half chr1: the old detector pins the
    /// better half and finds the other as its partner.
    #[test]
    fn old_path_reports_an_inter_chromosomal_se_chimera() {
        let gn = genome();
        let mut fwd = slice(&gn, 100, 60);
        fwd.extend(slice(&gn, CHR_PAD + 500, 60));
        let read = ChimRead::new(&fwd, [120, 0], "r", [&fwd, &[]]);
        let prm = params(&[]);
        let scorer = AlignmentScorer::from_params(&prm);
        let a = tx(0, 100, &[(Kind::Match, 60), (Kind::SoftClip, 60)], 59);
        let b = tx(
            1,
            CHR_PAD + 500,
            &[(Kind::SoftClip, 60), (Kind::Match, 60)],
            58,
        );
        let windows = vec![
            vec![WinTr::single(&a, 0, 0, 0, 120).unwrap()],
            vec![WinTr::single(&b, 0, 0, 0, 120).unwrap()],
        ];
        let out = chimeric_detection(windows, &read, &gn, &scorer, &prm);
        assert_eq!(out.len(), 1);
        let j = out[0].junction_line.as_ref().unwrap();
        assert_eq!((j.donor_chr, j.acceptor_chr), (0, 1));
        // chimJ0 is the first base past the donor, chimJ1 the base before the
        // acceptor.
        assert_eq!((j.donor_break, j.acceptor_break), (160, CHR_PAD + 499));
        assert_eq!(j.donor_cigar, "60M60S");
        assert_eq!(j.acceptor_cigar, "60S60M");
        assert!(out[0].multimap.is_none());
    }

    /// Two equally good partners: `chimScoreNext + chimScoreSeparation >=
    /// chimScoreBest`, so the old detector reports nothing.
    #[test]
    fn old_path_drops_an_ambiguous_partner() {
        let gn = genome();
        let mut fwd = slice(&gn, 100, 60);
        fwd.extend(slice(&gn, CHR_PAD + 500, 60));
        let read = ChimRead::new(&fwd, [120, 0], "r", [&fwd, &[]]);
        let prm = params(&[]);
        let scorer = AlignmentScorer::from_params(&prm);
        let a = tx(0, 100, &[(Kind::Match, 60), (Kind::SoftClip, 60)], 59);
        let b = tx(
            1,
            CHR_PAD + 500,
            &[(Kind::SoftClip, 60), (Kind::Match, 60)],
            58,
        );
        let b2 = tx(
            1,
            CHR_PAD + 800,
            &[(Kind::SoftClip, 60), (Kind::Match, 60)],
            58,
        );
        let windows = vec![
            vec![WinTr::single(&a, 0, 0, 0, 120).unwrap()],
            vec![WinTr::single(&b, 0, 0, 0, 120).unwrap()],
            vec![WinTr::single(&b2, 0, 0, 0, 120).unwrap()],
        ];
        assert!(chimeric_detection(windows, &read, &gn, &scorer, &prm).is_empty());
    }

    /// Mates on two chromosomes bracket the junction: STAR's `chimMotif = -1`
    /// branch, with no scan, and each segment printed as its own mate.
    #[test]
    fn pe_mates_bracketing_the_junction_are_type_minus_one() {
        let gn = genome();
        let m1 = slice(&gn, 100, 60);
        let m2_rc = slice(&gn, CHR_PAD + 500, 60); // RC(mate2) maps forward
        let m2 = rc(&m2_rc);
        let mut fwd = m1.clone();
        fwd.push(11);
        fwd.extend(&m2_rc);
        let read = ChimRead::new(&fwd, [60, 60], "r", [&m1, &m2]);
        let prm = params(&[]);
        let scorer = AlignmentScorer::from_params(&prm);
        let a = tx(0, 100, &[(Kind::Match, 60)], 59);
        let b = tx(1, CHR_PAD + 500, &[(Kind::Match, 60)], 59);
        let windows = vec![
            vec![WinTr::single(&a, 0, 0, 0, 121).unwrap()],
            vec![WinTr::single(&b, 0, 61, 1, 121).unwrap()],
        ];
        let out = chimeric_detection(windows, &read, &gn, &scorer, &prm);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].junction_type, -1);
        let j = out[0].junction_line.as_ref().unwrap();
        assert_eq!((j.donor_break, j.acceptor_break), (160, CHR_PAD + 499));
        assert_eq!(
            (j.donor_cigar.as_str(), j.acceptor_cigar.as_str()),
            ("60M", "60M")
        );
    }

    /// `--chimMultimapNmax` re-scores each stitched chimera with `alignScore`
    /// and reports the run-level columns.
    #[test]
    fn mult_path_rescores_and_reports_context() {
        let gn = genome();
        let mut fwd = slice(&gn, 100, 60);
        fwd.extend(slice(&gn, CHR_PAD + 500, 60));
        let read = ChimRead::new(&fwd, [120, 0], "r", [&fwd, &[]]);
        let prm = params(&["--chimMultimapNmax", "10"]);
        let scorer = AlignmentScorer::from_params(&prm);
        let a = tx(0, 100, &[(Kind::Match, 60), (Kind::SoftClip, 60)], 59);
        let b = tx(
            1,
            CHR_PAD + 500,
            &[(Kind::SoftClip, 60), (Kind::Match, 60)],
            59,
        );
        let windows = vec![
            vec![WinTr::single(&a, 0, 0, 0, 120).unwrap()],
            vec![WinTr::single(&b, 0, 0, 0, 120).unwrap()],
        ];
        let out = chimeric_detection(windows, &read, &gn, &scorer, &prm);
        assert_eq!(out.len(), 1);
        let m = out[0].multimap.unwrap();
        // Two exact 60-base blocks score 60 - 1 each after the span penalty;
        // a non-GT/AG junction costs `chimScoreJunctionNonGTAG` (-1) more.
        let expected = if out[0].junction_type == 0 { 117 } else { 118 };
        assert_eq!(m.chim_score, expected);
        assert_eq!(
            (m.chim_n, m.max_possible_score, m.max_non_chim_score),
            (1, 120, 59)
        );
    }
}
