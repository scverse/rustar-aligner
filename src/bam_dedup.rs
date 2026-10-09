//! `--runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType`: STAR's `bamRemoveDuplicates`.
//!
//! Marks PCR duplicates in a coordinate-sorted BAM. Every unique mapper (and every multimapper
//! when `mark_multi`, i.e. `UniqueIdentical`) is pre-marked with the SAM duplicate flag `0x400`;
//! reads are grouped into maximal overlapping clusters; within a cluster the unique mappers are
//! paired by name and, per run of pairs with an identical (start, flag, S-extended CIGAR) key, the
//! highest-`AS` pair is un-marked (kept). All records are re-emitted in the input order with the
//! flag adjusted. `NH`/`AS` are required (a missing tag is a fatal error).
//!
//! Near-verbatim port of the sister project STAR-rs's `crates/star-io/src/bam.rs`
//! (`bam_remove_duplicates` and its helpers), which cites the exact upstream STAR routines
//! (`funCompareCoordFlagCigarSeq`, `funCigarExtendS`, `funStartExtendS`).

use std::io::{self, Read, Write};
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::Value;

use crate::params::Parameters;

/// Per-record fields for dedup (STAR's dedup key + grouping coordinates).
struct DedupRec {
    tid: i64,
    pos: i64,
    flag: u16, // pre-marked/un-marked in place
    nh: u32,
    as_score: Option<i64>,
    name: Vec<u8>,
    mate_pos: i64,
    /// `funStartExtendS`: alignment start with a leading soft-clip unclipped.
    unclip: i64,
    /// `funCigarExtendS`: CIGAR with leading/trailing `S` folded into the adjacent op, packed
    /// `(len << 4) | op` (STAR op codes M=0,I=1,D=2,N=3,S=4,...).
    s_cigar: Vec<u32>,
    /// BAM-packed SEQ (two 4-bit `seq_nt16` codes per byte, first base in the high nibble) and its
    /// base count; only populated when `--bamRemoveDuplicatesMate2basesN > 0` (mate-2 SEQ compare).
    seq_packed: Vec<u8>,
    seq_len: u32,
}

/// Map an ASCII base to its BAM `seq_nt16` 4-bit code (`=ACMGRSVTWYHKDBN`).
fn nt16(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15, // N and anything else
    }
}

/// Pack ASCII bases into BAM nibble bytes (first base in the high nibble), as STAR reads them.
fn pack_seq(bases: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; bases.len().div_ceil(2)];
    for (i, &b) in bases.iter().enumerate() {
        let code = nt16(b);
        if i.is_multiple_of(2) {
            out[i / 2] |= code << 4;
        } else {
            out[i / 2] |= code;
        }
    }
    out
}

/// STAR `funCompareCoordFlagCigarSeq`'s mate-2 SEQ comparison (`bamRemoveDuplicatesMate2basesN`):
/// compare the first `n` bases of mate 2's stored SEQ (from the 5' end for a forward mate, from the
/// 3' end for a reverse-complemented mate), byte by byte over the packed nibbles, exactly as STAR.
fn cmp_mate2_seq(a2: &DedupRec, b2: &DedupRec, n: u32) -> std::cmp::Ordering {
    use std::cmp::Ordering::Equal;
    if n == 0 {
        return Equal;
    }
    let (sa, sb) = (&a2.seq_packed, &b2.seq_packed);
    let byte = |s: &[u8], i: u32| s.get((i / 2) as usize).copied().unwrap_or(0);
    if a2.flag & 0x10 == 0 {
        // forward: bytes 0..(n-1)/2, then the high nibble of the last byte when n is odd
        let mut ii = 1u32;
        while ii < n {
            let c = byte(sa, ii).cmp(&byte(sb, ii));
            if c != Equal {
                return c;
            }
            ii += 2;
        }
        if !n.is_multiple_of(2) {
            return (byte(sa, ii) >> 4).cmp(&(byte(sb, ii) >> 4));
        }
        Equal
    } else {
        // reverse: the last n bases, starting at seq_len - n
        let len = a2.seq_len;
        let mut ii = len.saturating_sub(n);
        if !ii.is_multiple_of(2) {
            let c = (byte(sa, ii) & 15).cmp(&(byte(sb, ii) & 15));
            if c != Equal {
                return c;
            }
            ii += 1;
        }
        while ii < len {
            let c = byte(sa, ii).cmp(&byte(sb, ii));
            if c != Equal {
                return c;
            }
            ii += 2;
        }
        Equal
    }
}

/// STAR's `funCigarExtendS`: fold a leading/trailing soft-clip into the adjacent op.
fn s_ext_cigar(cigar: &[u32]) -> Vec<u32> {
    let n = cigar.len();
    let mut out: Vec<u32> = if n > 0 && (cigar[0] & 0xf) == 4 {
        let mut o = cigar[1..].to_vec();
        if !o.is_empty() {
            o[0] += (cigar[0] >> 4) << 4; // add the leading-S length to the next op
        }
        o
    } else {
        cigar.to_vec()
    };
    if n > 0 && (cigar[n - 1] & 0xf) == 4 {
        out.pop(); // drop the trailing S ...
        if let Some(last) = out.last_mut() {
            *last += (cigar[n - 1] >> 4) << 4; // ... and add its length to the previous op
        }
    }
    out
}

/// STAR's `funCompareCigarsExtendS`: compare by op count, then op by op.
fn cmp_cigar(a: &[u32], b: &[u32]) -> std::cmp::Ordering {
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

/// Compare a mate pair `(a1, a2)` vs `(b1, b2)` (STAR's `funCompareCoordFlagCigarSeq`). With
/// `mate2_bases_n > 0` the first `n` bases of mate 2's SEQ are also compared (RAMPAGE); at `0` the
/// SEQ is not compared.
fn cmp_pair(
    a1: &DedupRec,
    a2: &DedupRec,
    b1: &DedupRec,
    b2: &DedupRec,
    mate2_bases_n: u32,
) -> std::cmp::Ordering {
    use std::cmp::Ordering::Equal;
    a1.unclip
        .cmp(&b1.unclip)
        .then_with(|| a2.unclip.cmp(&b2.unclip))
        .then_with(|| a1.flag.cmp(&b1.flag))
        .then_with(|| a2.flag.cmp(&b2.flag))
        .then_with(|| cmp_cigar(&a1.s_cigar, &b1.s_cigar))
        .then_with(|| cmp_cigar(&a2.s_cigar, &b2.s_cigar))
        .then_with(|| cmp_mate2_seq(a2, b2, mate2_bases_n))
        .then(Equal)
}

/// STAR `bamRemoveDuplicates`: mark PCR duplicates in a coordinate-sorted BAM.
pub fn bam_remove_duplicates<R: Read, W: Write>(
    input: R,
    output: W,
    mark_multi: bool,
    mate2_bases_n: u32,
) -> io::Result<()> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut reader = bam::io::Reader::new(input);
    let header = reader.read_header()?;

    let mut records = Vec::new();
    let mut recs = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result?;
        let tid = rec.reference_sequence_id().map_or(-1, |i| i as i64);
        let pos = rec
            .alignment_start()
            .map_or(-1, |p| usize::from(p) as i64 - 1);
        let mate_pos = rec
            .mate_alignment_start()
            .map_or(-1, |p| usize::from(p) as i64 - 1);
        let mut flag = u16::from(rec.flags());
        let data = rec.data();
        let nh = data
            .get(b"NH")
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "SAM tag NH is required for --bamRemoveDuplicatesType (re-generate with NH+AS)",
                )
            })?
            .as_int()
            .unwrap_or(1) as u32;
        let as_score = data.get(b"AS").and_then(Value::as_int);
        // Packed CIGAR (len<<4 | STAR op code).
        let cigar: Vec<u32> = rec
            .cigar()
            .as_ref()
            .iter()
            .map(|op| {
                let k = match op.kind() {
                    Kind::Match => 0,
                    Kind::Insertion => 1,
                    Kind::Deletion => 2,
                    Kind::Skip => 3,
                    Kind::SoftClip => 4,
                    Kind::HardClip => 5,
                    Kind::Pad => 6,
                    Kind::SequenceMatch => 7,
                    Kind::SequenceMismatch => 8,
                };
                ((op.len() as u32) << 4) | k
            })
            .collect();
        let unclip = if !cigar.is_empty() && (cigar[0] & 0xf) == 4 {
            pos - (cigar[0] >> 4) as i64
        } else {
            pos
        };
        // Pre-mark as duplicate (unique always; multi only when mark_multi).
        if nh == 1 || (nh > 1 && mark_multi) {
            flag |= 0x400;
        }
        // Mate-2 SEQ, packed as BAM nibbles, only when the RAMPAGE mate-2 compare is on.
        let (seq_packed, seq_len) = if mate2_bases_n > 0 {
            let bases: Vec<u8> = rec.sequence().as_ref().to_vec();
            let len = bases.len() as u32;
            (pack_seq(&bases), len)
        } else {
            (Vec::new(), 0)
        };
        recs.push(DedupRec {
            tid,
            pos,
            flag,
            nh,
            as_score,
            name: rec.name().unwrap_or_default().to_vec(),
            mate_pos,
            unclip,
            s_cigar: s_ext_cigar(&cigar),
            seq_packed,
            seq_len,
        });
        records.push(rec);
    }

    // Group into maximal overlapping clusters; collapse each (unique mappers only).
    let n = recs.len();
    let mut group_start = 0usize;
    let mut right_max: i64 = 0;
    let mut group: Vec<usize> = Vec::new();
    let mut e = 0usize;
    while e <= n {
        let boundary = e == n
            || recs[e].tid != recs[group_start].tid
            || (right_max > 0 && recs[e].pos > right_max);
        if boundary {
            collapse_group(&mut recs, &group, mate2_bases_n)?;
            if e == n {
                break;
            }
            right_max = 0;
            group_start = e;
            group.clear();
        }
        if recs[e].nh == 1 {
            group.push(e);
            if recs[e].mate_pos > recs[e].pos {
                right_max = right_max.max(recs[e].mate_pos);
            }
        }
        e += 1;
    }

    // Re-emit every record with its (possibly adjusted) flag.
    let mut writer = bam::io::Writer::new(output);
    writer.write_header(&header)?;
    for (rec, d) in records.iter_mut().zip(&recs) {
        *rec.flags_mut() = Flags::from(d.flag);
        writer.write_alignment_record(&header, rec)?;
    }
    writer.try_finish()?;
    Ok(())
}

/// Collapse one overlapping cluster of unique mappers (`group` = record indices): pair mates by name,
/// sort the pairs by the dedup key, and un-mark (clear `0x400` on both mates of) the highest-`AS` pair
/// in each run of identical pairs.
fn collapse_group(recs: &mut [DedupRec], group: &[usize], mate2_bases_n: u32) -> io::Result<()> {
    // Pair mates: sort by (name length, name bytes, mate bit 0x80) so each read's two mates adjoin.
    let mut sorted = group.to_vec();
    sorted.sort_by(|&a, &b| {
        let (ra, rb) = (&recs[a], &recs[b]);
        ra.name
            .len()
            .cmp(&rb.name.len())
            .then_with(|| ra.name.cmp(&rb.name))
            .then_with(|| (ra.flag & 0x80).cmp(&(rb.flag & 0x80)))
    });
    // Pairs of adjacent (mate1, mate2) indices.
    let mut pairs: Vec<(usize, usize)> = sorted
        .as_chunks::<2>()
        .0
        .iter()
        .map(|&[m1, m2]| (m1, m2))
        .collect();
    // Sort pairs by the coordinate/flag/CIGAR key.
    pairs.sort_by(|&(a1, a2), &(b1, b2)| {
        cmp_pair(&recs[a1], &recs[a2], &recs[b1], &recs[b2], mate2_bases_n)
    });

    let np = pairs.len();
    let mut best_score = i64::MIN;
    let mut best = 0usize;
    for pp in 0..np {
        let (m1, _m2) = pairs[pp];
        let score = recs[m1].as_score.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "SAM tag AS is required for --bamRemoveDuplicatesType (re-generate with NH+AS)",
            )
        })?;
        if score > best_score {
            best_score = score;
            best = pp;
        }
        let run_end = pp == np - 1 || {
            let (a1, a2) = pairs[pp];
            let (b1, b2) = pairs[pp + 1];
            cmp_pair(&recs[a1], &recs[a2], &recs[b1], &recs[b2], mate2_bases_n)
                != std::cmp::Ordering::Equal
        };
        if run_end {
            // Un-mark (clear the pre-set 0x400) the best pair's two mates; the rest stay marked.
            let (k1, k2) = pairs[best];
            recs[k1].flag ^= 0x400;
            recs[k2].flag ^= 0x400;
            best_score = i64::MIN;
        }
    }
    Ok(())
}

/// `--runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType`: mark PCR duplicates in
/// `--inputBAMfile` and write `<prefix>Processed.out.bam`.
pub fn run(params: &Parameters) -> anyhow::Result<()> {
    let dedup = params.bam_remove_duplicates_type.as_str();
    let mark_multi = dedup.eq_ignore_ascii_case("UniqueIdentical");
    let mate2_bases_n = params.bam_remove_duplicates_mate2_bases_n.max(0) as u32;

    let out = params.output_path("Processed.out.bam");
    if let Some(parent) = out.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)?;
    }

    if params.input_bam_file == "-" {
        let stdin = io::stdin();
        let output = std::fs::File::create(&out)?;
        bam_remove_duplicates(stdin.lock(), output, mark_multi, mate2_bases_n)?;
    } else {
        let input = std::fs::File::open(Path::new(&params.input_bam_file))?;
        let output = std::fs::File::create(&out)?;
        bam_remove_duplicates(input, output, mark_multi, mate2_bases_n)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;

    fn mate2(seq: &[u8], flag: u16) -> DedupRec {
        DedupRec {
            tid: 0,
            pos: 0,
            flag,
            nh: 1,
            as_score: Some(0),
            name: Vec::new(),
            mate_pos: -1,
            unclip: 0,
            s_cigar: Vec::new(),
            seq_packed: pack_seq(seq),
            seq_len: seq.len() as u32,
        }
    }

    #[test]
    fn pack_seq_nibbles() {
        // A=1, C=2, G=4, T=8: two bases per byte, first in the high nibble.
        assert_eq!(pack_seq(b"ACGT"), vec![0x12, 0x48]);
        assert_eq!(pack_seq(b"ACG"), vec![0x12, 0x40]);
    }

    #[test]
    fn cmp_mate2_seq_forward_matches() {
        let a = mate2(b"ACGTACGT", 0);
        let b = mate2(b"ACGTTTTT", 0);
        // First 4 bases identical -> Equal.
        assert_eq!(cmp_mate2_seq(&a, &b, 4), Ordering::Equal);
        // Full length differs at base 5 (T vs T is same at 5..) -> check divergence further out.
        let c = mate2(b"ACGAACGT", 0);
        assert_eq!(cmp_mate2_seq(&a, &c, 4), b'T'.cmp(&b'A'));
    }

    #[test]
    fn cmp_mate2_seq_zero_bases_is_equal() {
        let a = mate2(b"AAAA", 0);
        let b = mate2(b"TTTT", 0x10);
        assert_eq!(cmp_mate2_seq(&a, &b, 0), Ordering::Equal);
    }

    #[test]
    fn cmp_mate2_seq_reverse_reads_from_the_end() {
        // Reverse-complemented mate (flag 0x10): compare the *last* n bases.
        let a = mate2(b"TTTTACGT", 0x10);
        let b = mate2(b"GGGGACGT", 0x10);
        assert_eq!(cmp_mate2_seq(&a, &b, 4), Ordering::Equal);
        let c = mate2(b"GGGGACGA", 0x10);
        assert_eq!(cmp_mate2_seq(&a, &c, 4), b'T'.cmp(&b'A'));
    }

    #[test]
    fn s_ext_cigar_folds_leading_and_trailing_softclip() {
        // 5S10M3S -> packed [ (5<<4|4), (10<<4|0), (3<<4|4) ] -> folds to single op [ (18<<4|0) ]
        let cigar = vec![(5 << 4) | 4, 10 << 4, (3 << 4) | 4];
        assert_eq!(s_ext_cigar(&cigar), vec![18u32 << 4]);
    }

    #[test]
    fn s_ext_cigar_no_softclip_is_unchanged() {
        let cigar = vec![20u32 << 4];
        assert_eq!(s_ext_cigar(&cigar), cigar);
    }

    #[test]
    fn cmp_cigar_by_len_then_ops() {
        let a = vec![10u32 << 4];
        let b = vec![5u32 << 4, 5u32 << 4];
        assert_eq!(cmp_cigar(&a, &b), Ordering::Less); // fewer ops first
        assert_eq!(cmp_cigar(&a, &a), Ordering::Equal);
    }
}
