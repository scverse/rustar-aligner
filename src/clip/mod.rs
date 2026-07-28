//! Adapter-aware read clipping (STAR `ClipMate` / `ParametersClip`, Hamming mode).
//!
//! Ported from STAR-rs `crates/star-align/src/star_clip.rs`. Extends the plain
//! `--clip5pNbases`/`--clip3pNbases` fixed clip with a 3' adapter Hamming scan
//! (`--clip3pAdapterSeq`/`--clip3pAdapterMMp`) and the after-adapter trim
//! (`--clip{5,3}pAfterAdapterNbases`).
//!
//! STAR maps the *clipped* read (the clipped bases don't participate in alignment)
//! but keeps them in the SAM record as **soft-clips**: the output SEQ is the full
//! original read and the CIGAR gains leading/trailing `S` ops for the clipped ends.
//! rustar-aligner matches this — it aligns the clipped read, then `io::sam`'s
//! `apply_read_clips` restores the full read and wraps the CIGAR with soft-clips
//! (strand-aware: a read-3' clip is a leading `S` on a reverse-strand record).
//! Unmapped reads keep the full read too (CIGAR `*`, no soft-clip).
//!
//! STAR clips a mate as two passes in order: 5' first (`clip5pNbases`, then
//! `clip5pAfterAdapterNbases`), then 3' on the 5'-clipped read (`clip3pNbases`,
//! then the 3' adapter Hamming scan, then `clip3pAfterAdapterNbases`).
//!
//! `--clipAdapterType Hamming` (STAR's default) uses the 3' Hamming scan above.
//! `CellRanger4` instead trims a 3' poly-A tail and a 5' TSO; see
//! [`cellranger4`]. A 5' *Hamming* adapter is not a thing STAR supports either.

pub mod cellranger4;

use crate::io::fastq::encode_base;
use crate::params::Parameters;

/// One end's clipping parameters (STAR `ClipMate`). `n` = fixed clip, `adapter` =
/// the adapter sequence (base codes A=0..T=3; empty = none), `ad_mmp` = max
/// mismatch fraction, `n_after` = extra clip after the adapter.
#[derive(Debug, Clone, Default)]
pub struct ClipEnd {
    /// `--clip{5,3}pNbases` for this end.
    pub n: usize,
    /// `--clip3pAdapterSeq` as base codes (A=0..T=3), empty when none. Only
    /// meaningful for the 3' end. Encoded once in [`clip_params_from`] so the
    /// scan matches the already-numeric read without re-encoding.
    pub adapter: Vec<u8>,
    /// `--clip3pAdapterMMp` (default 0.1).
    pub ad_mmp: f64,
    /// `--clip{5,3}pAfterAdapterNbases`.
    pub n_after: usize,
}

/// One mate's clipping parameters: the 5' and 3' ends.
#[derive(Debug, Clone, Default)]
pub struct ClipParams {
    /// 5' end (STAR `ClipMate` type 0).
    pub five: ClipEnd,
    /// 3' end (STAR `ClipMate` type 1).
    pub three: ClipEnd,
}

/// Build [`ClipParams`] for `mate` (0 or 1) from the run's `--clip{5,3}pNbases`
/// / `--clip3pAdapterSeq` / `--clip3pAdapterMMp` / `--clip{5,3}pAfterAdapterNbases`.
/// `-` (STAR's sentinel) means no adapter. The fixed `--clip{5,3}pNbases` are
/// per-mate (`clip5p(mate)`/`clip3p(mate)`); the adapter / mmp / after-adapter
/// clips are single-valued and apply to both mates. Cheap to build per batch.
pub fn clip_params_from(params: &Parameters, mate: usize) -> ClipParams {
    // Encode the adapter to base codes (A=0..T=3) ONCE here. The read reaching
    // clip_mate is already numeric (io::fastq encodes at read time), so the 3'
    // Hamming scan compares numeric-vs-numeric — re-encoding the read there turned
    // every base into N (code 4) and silently disabled all adapter clipping.
    let adapter = if params.clip3p_adapter_seq == "-" {
        Vec::new()
    } else {
        params.clip3p_adapter_seq.bytes().map(encode_base).collect()
    };
    ClipParams {
        five: ClipEnd {
            n: params.clip5p(mate),
            adapter: Vec::new(),
            ad_mmp: 0.0,
            n_after: params.clip5p_after_adapter_nbases as usize,
        },
        three: ClipEnd {
            n: params.clip3p(mate),
            adapter,
            ad_mmp: params.clip3p_adapter_mmp,
            n_after: params.clip3p_after_adapter_nbases as usize,
        },
    }
}

/// STAR `localSearch`: the offset `ixBest` in `x` where `y` best matches (max
/// matches, then min mismatches, subject to `nMM/nMatch <= p_mm`); `nx` when no
/// acceptable match (so nothing clips).
fn local_search(x: &[u8], y: &[u8], p_mm: f64) -> usize {
    let nx = x.len();
    let ny = y.len();
    let (mut best_ix, mut best_match, mut best_mm) = (nx, 0usize, 0usize);
    for ix in 0..nx {
        let (mut n_match, mut n_mm) = (0usize, 0usize);
        for iy in 0..ny.min(nx - ix) {
            if x[ix + iy] > 3 {
                continue;
            }
            if x[ix + iy] == y[iy] {
                n_match += 1;
            } else {
                n_mm += 1;
            }
        }
        if (n_match > best_match || (n_match == best_match && n_mm < best_mm))
            && (n_match == 0 || n_mm as f64 / n_match as f64 <= p_mm)
        {
            best_ix = ix;
            best_match = n_match;
            best_mm = n_mm;
        }
    }
    best_ix
}

/// Clip one mate (STAR's two `ClipMate::clip` passes: 5' then 3'). Returns
/// `(clip5p_total, clip3p_total)`; the clipped region is
/// `read[clip5p_total .. len-clip3p_total]`.
///
/// STAR marks a `ClipMate` inactive (no clip at all) when it has no fixed clip
/// AND no adapter: an end's `n_after` alone is a no-op, matching the oracle.
pub fn clip_mate(read: &[u8], p: &ClipParams) -> (usize, usize) {
    let len = read.len();

    // ---- 5' end (STAR ClipMate type 0) ----
    let five_active = p.five.n > 0;
    let mut c5 = 0;
    if five_active {
        c5 = p.five.n.min(len);
        if p.five.n_after > 0 {
            c5 += p.five.n_after.min(len - c5);
        }
    }

    // ---- 3' end (STAR ClipMate type 1), on the 5'-clipped read ----
    let three_active = p.three.n > 0 || !p.three.adapter.is_empty();
    let s = &read[c5..];
    let sl = s.len();
    let mut c3 = 0;
    if three_active {
        c3 = p.three.n.min(sl);
        if !p.three.adapter.is_empty() {
            // Hamming 3' adapter scan on the (already-numeric) read after the fixed
            // 3' clip. Both `x` and the adapter are base codes; local_search skips
            // N (code > 3) internally.
            let x_num = &s[..sl - c3];
            let nx = x_num.len();
            let ix_best = local_search(x_num, &p.three.adapter, p.three.ad_mmp);
            c3 += nx - ix_best;
        }
        if p.three.n_after > 0 {
            let remaining = sl.saturating_sub(c3);
            c3 += p.three.n_after.min(remaining);
        }
    }

    (c5, c3)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encode an ASCII sequence to base codes (A=0..T=3), matching how `io::fastq`
    /// stores reads. `clip_mate` receives numeric reads in production, so every
    /// test feeds numeric input — passing raw ASCII would exercise a path the real
    /// pipeline never hits (and previously masked the double-encode dead-code bug).
    fn enc(ascii: &[u8]) -> Vec<u8> {
        ascii.iter().map(|&b| encode_base(b)).collect()
    }

    /// 3'-adapter-only params; the adapter is stored as base codes exactly as
    /// [`clip_params_from`] produces it.
    fn hamming_3p(adapter: &[u8], mmp: f64) -> ClipParams {
        ClipParams {
            three: ClipEnd {
                adapter: enc(adapter),
                ad_mmp: mmp,
                ..Default::default()
            },
            ..Default::default()
        }
    }

    #[test]
    fn fixed_5p_3p() {
        let p = ClipParams {
            five: ClipEnd {
                n: 3,
                ..Default::default()
            },
            three: ClipEnd {
                n: 2,
                ..Default::default()
            },
        };
        assert_eq!(clip_mate(&enc(b"AAACCCGGGTT"), &p), (3, 2));
    }

    #[test]
    fn adapter_3p_hamming() {
        // Adapter "AGATCGGAAGAGC"; read = insert + adapter tail -> the 13-base tail is clipped.
        let p = hamming_3p(b"AGATCGGAAGAGC", 0.1);
        let (c5, c3) = clip_mate(&enc(b"CCCCGGGGAGATCGGAAGAGC"), &p);
        assert_eq!((c5, c3), (0, 13));
    }

    #[test]
    fn adapter_3p_mmp_rejects_when_strict() {
        // One mismatch in a 13-base adapter: allowed at 0.1 (1/12<=0.1), rejected at 0.0.
        let read = enc(b"CCCCGGGGAGATCGGATGAGC"); // 'T' at pos 16 -> mismatch vs adapter 'G'
        let (_, c3_loose) = clip_mate(&read, &hamming_3p(b"AGATCGGAAGAGC", 0.1));
        assert_eq!(c3_loose, 13);
        let (_, c3_strict) = clip_mate(&read, &hamming_3p(b"AGATCGGAAGAGC", 0.0));
        assert_eq!(c3_strict, 0);
    }

    #[test]
    fn after_adapter_alone_is_noop() {
        // STAR's inactive-end short-circuit: n_after with no fixed clip and no adapter clips nothing.
        let p = ClipParams {
            five: ClipEnd {
                n_after: 4,
                ..Default::default()
            },
            three: ClipEnd {
                n_after: 4,
                ..Default::default()
            },
        };
        assert_eq!(clip_mate(&enc(b"ACGTACGTACGTACGTACGT"), &p), (0, 0));
    }

    #[test]
    fn after_adapter_nbases_3p() {
        let mut p = hamming_3p(b"AGATCGGAAGAGC", 0.1);
        p.three.n_after = 3;
        let (_, c3) = clip_mate(&enc(b"CCCCGGGGAGATCGGAAGAGC"), &p);
        assert_eq!(c3, 16); // 13 adapter + 3 after
    }

    #[test]
    fn no_adapter_configured_only_fixed_clips() {
        let p = ClipParams {
            five: ClipEnd {
                n: 2,
                ..Default::default()
            },
            three: ClipEnd {
                n: 1,
                ..Default::default()
            },
        };
        assert_eq!(clip_mate(&enc(b"AAACCCGGGTT"), &p), (2, 1));
    }

    /// Regression guard for the double-encode dead-code bug: a numeric read (as the
    /// pipeline delivers) whose 3' tail is the adapter MUST clip. Before the fix,
    /// clip_mate re-encoded the numeric read to all-N and never clipped.
    #[test]
    fn numeric_read_adapter_actually_clips() {
        let p = hamming_3p(b"AGATCGGAAGAGC", 0.1);
        // 8-base insert + 13-base adapter, delivered as base codes.
        let read = enc(b"CCCCGGGGAGATCGGAAGAGC");
        assert!(read.iter().all(|&b| b <= 3), "test read must be numeric");
        let (_, c3) = clip_mate(&read, &p);
        assert_eq!(c3, 13, "adapter tail must be clipped on a numeric read");
    }
}
