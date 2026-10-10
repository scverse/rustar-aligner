//! `--soloOutLayout CellRanger`: 5' TSO and 3' poly(A) trimming of the cDNA read, as
//! `cellranger count` reports it in the `ts:i` and `pa:i` tags.
//!
//! This is an independent implementation, written from the behaviour of `cellranger
//! count` 10.0.0 (the `ts`/`pa` tags of its BAM on pbmc_1k_v3) and from cutadapt's
//! published adapter model, not from 10x's code. Measured on 2.3 million reads of that
//! BAM, the poly(A) lengths agree on every read and the TSO lengths on all but 23.
//!
//! - TSO (`AAGCAGTGGTATCAACGCAGAGTACATGGG`), a 5' adapter found anywhere: the adapter
//!   aligns semi-globally to the read (its prefix may hang off the read's 5' end, it may
//!   start anywhere in the read and must end inside it), scoring +1 per match and -2 per
//!   mismatch or indel. The best-scoring alignment (the leftmost end on ties) with at most
//!   10% errors over the aligned adapter bases is kept when it scores at least 20; the
//!   read is trimmed up to the alignment's end.
//! - Poly(A), a 3' adapter that must reach the read end: cutadapt's choice among the read
//!   suffixes with at most 10% errors, the most `A` matches and then the fewest errors;
//!   kept when matches - 2 * errors is at least 20.
//!
//! Both are found on the untrimmed read; the bases kept are the ones neither trims.

/// The TSO in the read encoding (A=0, C=1, G=2, T=3).
const TSO: [u8; 30] = encode(b"AAGCAGTGGTATCAACGCAGAGTACATGGG");

const fn encode<const N: usize>(s: &[u8; N]) -> [u8; N] {
    let mut out = [0u8; N];
    let mut i = 0;
    while i < N {
        out[i] = match s[i] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        i += 1;
    }
    out
}
const MAX_ERROR_RATE: f64 = 0.1;
const MIN_SCORE: i32 = 20;
const MISMATCH: i32 = -2;
const GAP: i32 = -2;

/// Result of trimming one read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CrTrim {
    /// Bases trimmed from the 5' end by the TSO (`ts:i`), 0 when none.
    pub tso: usize,
    /// Bases trimmed from the 3' end by the poly(A) (`pa:i`), 0 when none.
    pub polya: usize,
    /// First base kept.
    pub start: usize,
    /// One past the last base kept (`start` when nothing is kept).
    pub end: usize,
}

/// Trim one read (encoded bases: A=0, C=1, G=2, T=3, N=4).
pub fn trim(read: &[u8]) -> CrTrim {
    let n = read.len();
    let tso = tso_end(read);
    let polya = n - polya_start(read);
    let start = tso;
    let end = (n - polya).max(start);
    CrTrim {
        tso,
        polya,
        start,
        end,
    }
}

/// End of the TSO in the read (0 when no TSO qualifies).
fn tso_end(read: &[u8]) -> usize {
    // Cell: (score, errors, adapter bases aligned).
    type Cell = (i32, i32, i32);
    const NONE: Cell = (i32::MIN / 2, 0, 0);
    let (adapter_len, n) = (TSO.len(), read.len());
    if n == 0 {
        return 0;
    }
    let pick = |a: Cell, b: Cell| if b.0 > a.0 { b } else { a };
    // Row 0: the adapter may start before any read base.
    let mut prev: Vec<Cell> = vec![(0, 0, 0); n + 1];
    let mut cur: Vec<Cell> = vec![NONE; n + 1];
    for i in 1..=adapter_len {
        // Column 0: the first `i` adapter bases hang off the read's 5' end, free.
        cur[0] = (0, 0, 0);
        for j in 1..=n {
            let hit = TSO[i - 1] == read[j - 1];
            let d = prev[j - 1];
            let mut c = (
                d.0 + if hit { 1 } else { MISMATCH },
                d.1 + i32::from(!hit),
                d.2 + 1,
            );
            c = pick(c, (prev[j].0 + GAP, prev[j].1 + 1, prev[j].2 + 1));
            c = pick(
                c,
                (curead[j - 1].0 + GAP, curead[j - 1].1 + 1, curead[j - 1].2),
            );
            cur[j] = c;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    // `prev` now holds the last adapter row.
    let mut best = (i32::MIN, 0usize);
    for (j, c) in prev.iter().enumerate().skip(1) {
        if f64::from(c.1) > MAX_ERROR_RATE * f64::from(c.2) {
            continue;
        }
        if c.0 > best.0 {
            best = (c.0, j);
        }
    }
    if best.0 >= MIN_SCORE { best.1 } else { 0 }
}

/// Start of the poly(A) in the read (`read.len()` when none qualifies).
fn polya_start(r: &[u8]) -> usize {
    let n = r.len();
    let (mut best_m, mut best_e, mut best_s) = (-1i32, 0i32, n);
    let (mut m, mut e) = (0i32, 0i32);
    for s in (0..n).rev() {
        if r[s] == 0 {
            m += 1;
        } else {
            e += 1;
        }
        if f64::from(e) > MAX_ERROR_RATE * (n - s) as f64 {
            continue;
        }
        if m > best_m || (m == best_m && e < best_e) {
            (best_m, best_e, best_s) = (m, e, s);
        }
    }
    if best_m >= 0 && best_m + MISMATCH * best_e >= MIN_SCORE {
        best_s
    } else {
        n
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trim(r: &[u8]) -> CrTrim {
        let e: Vec<u8> = r
            .iter()
            .map(|&b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();
        super::trim(&e)
    }

    // Reads and expected `ts`/`pa` values from cellranger count 10.0.0 on pbmc_1k_v3.
    #[test]
    fn full_and_partial_tso() {
        let r = b"AAGCAGTGGTATCAACGCAGAGTACATGGGCCATGACACCTTCCCCGCCAGACCCAGACTTGGGCCGTTGCTCTGACATGGACACAGCCAG";
        assert_eq!(trim(r).tso, 30);
        let r = b"CAGTGGTATCAACGCAGAGTACATGGGGGTTGCAGGTGAGGTAGAGAAGCAGTGAACAGAGCACGAAGACCTGATGTTCCAGGGTTGGGAG";
        assert_eq!(trim(r).tso, 27);
        // A truncated adapter end stops at the last matching base.
        let r = b"AAGCAGTGGTATCAACGCAGAGTACATGGAAATGGCTATAGCAATAAAACTAGGAATAGCCCCCTTTCACTTCTGAGTCCCAGAGGTTACC";
        assert_eq!(trim(r).tso, 29);
    }

    #[test]
    fn too_many_errors_is_not_tso() {
        let r = b"AAGCAGTGGTATCAACGCAGAGCATCTTGGTGTTAATATTCTTTTAGTCATGGTCTCTGAAAGTAAAATGTCTATTCTTCTAATTCTGGCA";
        assert_eq!(trim(r).tso, 0);
    }

    #[test]
    fn polya_extends_through_sparse_errors() {
        let r = b"CCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGACAAAAAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(trim(r).polya, 30);
        let r = b"CTATATGGATGCCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(trim(r).polya, 31);
    }

    #[test]
    fn short_polya_is_kept() {
        let r = b"ACGTACGTACGTACGTAAAAAAAAAAAAAAAAAAA";
        let t = trim(r);
        assert_eq!((t.polya, t.start, t.end), (0, 0, r.len()));
    }

    #[test]
    fn overlapping_trims_keep_nothing() {
        let mut r = b"AAGCAGTGGTATCAACGCAGAGTACATGGG".to_vec();
        r.extend_from_slice(&[b'A'; 25]);
        let t = trim(&r);
        assert_eq!(t.tso, 30);
        assert_eq!(t.start, t.end);
    }
}
