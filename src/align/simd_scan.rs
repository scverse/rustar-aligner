//! SIMD-accelerated "find first stop" scan for the SA seed-extension hot path
//! (`compare_seq_to_genome` in `seed.rs`).
//!
//! The scan STAR performs there is: walk read/genome bytes in lockstep,
//! stopping at the first position where the genome byte is padding (`>= 5`)
//! or the two bytes differ. This is exactly a vectorizable "find first
//! difference" scan (like a SIMD memcmp) — most of a seed's extension length
//! matches, so the common case is scanning many equal bytes before the first
//! (or no) stop.
//!
//! Design choice: **only** the cheap "does this 16-byte chunk contain any
//! stop condition at all" check is arch-specific (SSE2 / NEON intrinsics,
//! a handful of lines each). The moment a chunk reports a stop, the exact
//! position is found by a plain, portable, easy-to-audit scalar re-scan of
//! that one 16-byte chunk — this happens at most once per call, so its cost
//! is negligible, and it means the only unsafe/arch-specific surface is a
//! trivial boolean reduction, not position-extraction logic.

/// Find the first index in `read`/`genome` (equal-length slices) where
/// `genome[i] >= 5` (padding) or `read[i] != genome[i]`. Returns `None` if no
/// such position exists (the whole range matches).
///
/// # Panics
/// Debug-asserts `read.len() == genome.len()`; callers must pass equal-length
/// slices (guaranteed by `compare_seq_to_genome`, which bounds both from the
/// same loop counter).
pub fn find_stop(read: &[u8], genome: &[u8]) -> Option<usize> {
    debug_assert_eq!(read.len(), genome.len());

    let mut base = 0usize;
    let (r_chunks, r_rem) = read.as_chunks::<16>();
    let (g_chunks, g_rem) = genome.as_chunks::<16>();
    for (rc, gc) in r_chunks.iter().zip(g_chunks) {
        if !chunk_all_match(rc, gc) {
            for k in 0..16 {
                if gc[k] >= 5 || rc[k] != gc[k] {
                    return Some(base + k);
                }
            }
            unreachable!("chunk_all_match reported a stop but the scalar re-scan found none");
        }
        base += 16;
    }

    // Tail shorter than 16 bytes: plain scalar scan.
    for (k, (&r, &g)) in r_rem.iter().zip(g_rem.iter()).enumerate() {
        if g >= 5 || r != g {
            return Some(base + k);
        }
    }
    None
}

/// `true` iff every lane matches (`read[k] == genome[k]`) and no lane is
/// genome padding (`genome[k] >= 5`, i.e. `> 4` since values are 0..=5).
#[cfg(target_arch = "x86_64")]
#[inline]
fn chunk_all_match(read: &[u8; 16], genome: &[u8; 16]) -> bool {
    use std::arch::x86_64::{
        _mm_cmpeq_epi8, _mm_cmpgt_epi8, _mm_loadu_si128, _mm_movemask_epi8, _mm_set1_epi8,
    };
    // SAFETY: `_mm_loadu_si128` reads exactly 16 bytes from each pointer, and
    // `read`/`genome` are `&[u8; 16]` — always exactly 16 valid bytes. SSE2 is
    // part of the mandatory x86_64 baseline ABI, so no runtime feature check
    // is needed (unlike AVX2/AVX-512, which are opt-in via `target-cpu`).
    unsafe {
        let r = _mm_loadu_si128(read.as_ptr().cast());
        let g = _mm_loadu_si128(genome.as_ptr().cast());
        let eq = _mm_cmpeq_epi8(r, g); // 0xFF per lane where read == genome
        let pad = _mm_cmpgt_epi8(g, _mm_set1_epi8(4)); // 0xFF per lane where genome > 4 (i.e. == 5)
        (_mm_movemask_epi8(eq) as u32 & 0xFFFF) == 0xFFFF && _mm_movemask_epi8(pad) == 0
    }
}

/// `true` iff every lane matches and no lane is genome padding. See the
/// x86_64 variant's doc for the shared semantics.
#[cfg(target_arch = "aarch64")]
#[inline]
fn chunk_all_match(read: &[u8; 16], genome: &[u8; 16]) -> bool {
    use std::arch::aarch64::{vceqq_u8, vcgtq_u8, vdupq_n_u8, vld1q_u8, vmaxvq_u8, vminvq_u8};
    // SAFETY: `vld1q_u8` reads exactly 16 bytes from each pointer, and
    // `read`/`genome` are `&[u8; 16]` — always exactly 16 valid bytes. NEON is
    // part of the mandatory aarch64 baseline (unlike SVE, which is opt-in).
    unsafe {
        let r = vld1q_u8(read.as_ptr());
        let g = vld1q_u8(genome.as_ptr());
        let eq = vceqq_u8(r, g); // 0xFF per lane where read == genome, else 0x00
        let pad = vcgtq_u8(g, vdupq_n_u8(4)); // 0xFF per lane where genome > 4
        vminvq_u8(eq) == 0xFF && vmaxvq_u8(pad) == 0
    }
}

/// Portable fallback for any architecture without a dedicated SIMD kernel
/// above (the caller in `seed.rs` only calls into this module on hot paths;
/// this fallback keeps `find_stop` fully functional everywhere, just without
/// the chunked speedup).
#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
#[inline]
fn chunk_all_match(read: &[u8; 16], genome: &[u8; 16]) -> bool {
    read.iter()
        .zip(genome.iter())
        .all(|(&r, &g)| g < 5 && r == g)
}

#[cfg(test)]
mod tests {
    use super::find_stop;

    /// Reference implementation: the exact scalar semantics `find_stop` must
    /// reproduce, ported directly from the original `compare_seq_to_genome`
    /// loop body.
    fn find_stop_scalar_reference(read: &[u8], genome: &[u8]) -> Option<usize> {
        for (i, (&r, &g)) in read.iter().zip(genome.iter()).enumerate() {
            if g >= 5 || r != g {
                return Some(i);
            }
        }
        None
    }

    /// Small xorshift PRNG (deterministic, no OS entropy) for the fuzz loop.
    struct Xorshift(u64);
    impl Xorshift {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn next_range(&mut self, n: usize) -> usize {
            (self.next_u64() % n as u64) as usize
        }
    }

    #[test]
    fn hand_crafted_edge_cases() {
        // Empty.
        assert_eq!(find_stop(&[], &[]), None);
        // Fully matching, various lengths straddling the 16-byte chunk size.
        for len in [1usize, 15, 16, 17, 31, 32, 33, 63, 64, 100] {
            let read = vec![1u8; len];
            let genome = vec![1u8; len];
            assert_eq!(find_stop(&read, &genome), None, "len={len}");
        }
        // Mismatch at every position for a 33-byte buffer (crosses 2 full
        // chunks + a tail), one position at a time.
        for pos in 0..33 {
            let mut genome = vec![1u8; 33];
            genome[pos] = 2; // mismatch (read stays 1)
            let read = vec![1u8; 33];
            assert_eq!(find_stop(&read, &genome), Some(pos), "mismatch pos={pos}");
        }
        // Padding (value 5) at every position for a 33-byte buffer.
        for pos in 0..33 {
            let mut genome = vec![1u8; 33];
            genome[pos] = 5;
            let read = vec![1u8; 33];
            assert_eq!(find_stop(&read, &genome), Some(pos), "padding pos={pos}");
        }
        // N (value 4) in genome is NOT padding — must NOT stop (only >=5 stops).
        let read = vec![1u8; 20];
        let mut genome = vec![1u8; 20];
        genome[10] = 4;
        assert_eq!(
            find_stop(&read, &genome),
            Some(10),
            "N (4) in genome must still stop on the read/genome mismatch, not be treated as padding-only"
        );
        // But if read ALSO has 4 there (equal), and genome's 4 is not >=5, it's a match:
        let mut read2 = vec![1u8; 20];
        read2[10] = 4;
        assert_eq!(
            find_stop(&read2, &genome),
            None,
            "read==genome==4 (N==N) must match; N is not padding"
        );
        // Earliest of multiple stop conditions must win (first mismatch before
        // a later padding byte).
        let read3 = vec![1u8; 20];
        let mut genome3 = vec![1u8; 20];
        genome3[3] = 2; // mismatch first
        genome3[10] = 5; // padding later
        assert_eq!(find_stop(&read3, &genome3), Some(3));
    }

    #[test]
    fn fuzz_against_scalar_reference() {
        let mut rng = Xorshift(0xD1B5_4A32_D192_ED03);
        for _ in 0..20_000 {
            let len = rng.next_range(200) + 1; // 1..=200
            let read: Vec<u8> = (0..len).map(|_| (rng.next_range(6)) as u8).collect();
            let mut genome: Vec<u8> = (0..len).map(|_| (rng.next_range(6)) as u8).collect();
            // Bias roughly a third of cases toward a genuine long match prefix
            // by copying the read into genome up to a random cut point, then
            // resuming random bytes — otherwise almost every fuzz case stops
            // at index 0 (both bytes are drawn from 0..6 independently) and
            // the SIMD "all match" fast path barely gets exercised.
            if rng.next_range(3) != 0 {
                let cut = rng.next_range(len + 1);
                genome[..cut].copy_from_slice(&read[..cut]);
            }
            let expected = find_stop_scalar_reference(&read, &genome);
            let actual = find_stop(&read, &genome);
            assert_eq!(
                actual, expected,
                "mismatch for read={read:?} genome={genome:?}"
            );
        }
    }
}
