//! Streaming writer for STAR's [`PackedArray`] bit format.
//!
//! [`PackedArray`] stores `length` entries of `word_length` bits each,
//! LSB-first within each byte. The full `data_byte_len_for(word_length,
//! length)` block is `((length - 1) * word_length / 8) + 8` bytes (an
//! extra 8 bytes pad the end so [`PackedArray::read`] can do an
//! unaligned 8-byte load on the last entry without going out of
//! bounds).
//!
//! The genome-scale suffix array is `length = ~6 × 10⁹` and
//! `word_length = 33` bits on the human assembly, giving a ~25 GB
//! packed block. Materialising that block in RAM before writing it to
//! disk is wasteful: every entry is written exactly once, in
//! increasing order, and after construction the SA file is opened by
//! `--runMode alignReads` from disk anyway. [`PackedStreamWriter`]
//! lets the SA construction pipeline pack each entry directly into
//! the output file as caps-sa emits it, without ever holding the full
//! packed block in memory.
//!
//! ## Bit layout invariant
//!
//! Entry `i` occupies bits `[i*word_length, (i+1)*word_length)` of the
//! stream, with bit 0 being the LSB of byte 0. The streaming writer
//! maintains a `u128` accumulator (`pending`) holding the next
//! up-to-128 bits not yet flushed; once `pending` holds 8 or more
//! bits, the low byte is flushed and `pending >>= 8`. With
//! `word_length ≤ 64` and at most 7 bits already in `pending` from
//! the previous entry, accumulating one more entry leaves at most
//! `7 + 64 = 71` bits — comfortably under the `u128` cap, with
//! ample room to flush the freed bytes after.
//!
//! ## Correspondence with [`PackedArray::write`]
//!
//! `PackedArray::write` does read-modify-write on an 8-byte window
//! straddling each entry. The streaming writer never reads —
//! sequential writes mean each bit is set exactly once. As long as
//! we emit bytes in the same little-endian order, the resulting byte
//! sequence is **byte-for-byte identical** to
//! `PackedArray::data()`. The unit tests assert this for several
//! `word_length` values and random inputs.

use std::io::{self, Write};

use crate::index::packed_array::PackedArray;

/// Writes [`PackedArray`]-format entries directly to a `Write` sink
/// without holding the full bit-packed block in RAM.
///
/// See the module-level doc for the bit-layout invariant and the
/// rationale for the streaming approach.
pub struct PackedStreamWriter<W: Write> {
    writer: W,
    word_length: u32,
    /// Bits accumulated, LSB-first; flushed byte-by-byte from the
    /// low end. With `word_length ≤ 64` and at most 7 carry-over
    /// bits from the previous entry, this holds at most 71 bits
    /// before each flush — well under `u128::BITS`.
    pending: u128,
    /// Number of valid low bits in `pending`. Always `< 8` between
    /// `write_one` calls.
    bits: u32,
    /// Total entries committed so far. Drives the `finish` padding
    /// computation.
    n_written: usize,
}

impl<W: Write> PackedStreamWriter<W> {
    /// Create a writer that bit-packs entries of `word_length` bits
    /// each into `writer`.
    pub fn new(writer: W, word_length: u32) -> Self {
        assert!(
            word_length > 0 && word_length <= 64,
            "PackedStreamWriter: word_length must be 1..=64, got {word_length}"
        );
        Self {
            writer,
            word_length,
            pending: 0,
            bits: 0,
            n_written: 0,
        }
    }

    /// Append one entry. The low `word_length` bits of `value` are
    /// packed; higher bits are masked off.
    pub fn write_one(&mut self, value: u64) -> io::Result<()> {
        let mask: u128 = if self.word_length == 64 {
            u64::MAX as u128
        } else {
            (1u128 << self.word_length) - 1
        };
        let masked = (value as u128) & mask;
        self.pending |= masked << self.bits;
        self.bits += self.word_length;

        // Flush whole bytes from the low end. After each flush, the
        // low byte of `pending` is consumed and `bits` decreases by 8.
        while self.bits >= 8 {
            let byte = (self.pending & 0xFF) as u8;
            self.writer.write_all(std::slice::from_ref(&byte))?;
            self.pending >>= 8;
            self.bits -= 8;
        }
        self.n_written += 1;
        Ok(())
    }

    /// Total entries written so far.
    pub fn n_written(&self) -> usize {
        self.n_written
    }

    /// Flush any remaining partial byte, then pad with zero bytes
    /// until the total byte count equals
    /// [`PackedArray::data_byte_len_for`]`(word_length, n_written)`.
    /// Returns the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        // (1) The last partial byte (`bits` in `[1, 7]`) — emit it.
        // The high bits of that byte are implicitly zero, matching
        // the `(value & mask) << bit_shift` semantics of
        // `PackedArray::write` for the last entry's tail.
        if self.bits > 0 {
            let byte = (self.pending & 0xFF) as u8;
            self.writer.write_all(std::slice::from_ref(&byte))?;
            self.pending >>= 8;
            self.bits = 0;
        }

        // (2) Pad with zeros to STAR's `length_byte`. The `+ 8` in
        // STAR's formula reserves the 8-byte read window the `read`
        // path uses on the last entry; with sequential writes we
        // emit it explicitly.
        let expected = PackedArray::data_byte_len_for(self.word_length, self.n_written);
        let bits_per_entry = self.word_length as u64;
        let bits_written = self.n_written as u64 * bits_per_entry;
        let bytes_already_emitted = bits_written.div_ceil(8) as usize;
        let pad = expected.saturating_sub(bytes_already_emitted);
        if pad > 0 {
            // Stream zeros in a small buffer; `pad` is at most ~7 for
            // any reasonable `word_length`, but we don't rely on that.
            let zeros = [0u8; 64];
            let mut remaining = pad;
            while remaining > 0 {
                let chunk = remaining.min(zeros.len());
                self.writer.write_all(&zeros[..chunk])?;
                remaining -= chunk;
            }
        }
        Ok(self.writer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The streaming writer must produce the same byte sequence as a
    /// fully-materialised [`PackedArray`] for any valid input. Drive
    /// it on randomised inputs at several `word_length` values.
    fn assert_matches_packed_array(word_length: u32, values: &[u64]) {
        // Reference: build a PackedArray, write each value, dump data().
        let mut reference = PackedArray::new(word_length, values.len());
        for (i, &v) in values.iter().enumerate() {
            reference.write(i, v);
        }
        let want = reference.data().to_vec();

        // Streaming: pack values into a Vec<u8> through the writer.
        let mut got: Vec<u8> = Vec::new();
        let mut sw = PackedStreamWriter::new(&mut got, word_length);
        for &v in values {
            sw.write_one(v).unwrap();
        }
        let _w = sw.finish().unwrap();

        assert_eq!(
            got,
            want,
            "stream-write mismatch at word_length={word_length}, len={}",
            values.len()
        );
    }

    #[test]
    fn matches_packed_array_8bit() {
        // Whole-byte alignment — easiest case.
        let values: Vec<u64> = (0u8..=255).map(|v| v as u64).collect();
        assert_matches_packed_array(8, &values);
    }

    #[test]
    fn matches_packed_array_33bit() {
        // STAR's human-genome word width. Awkward (odd) and the most
        // common production setting.
        let values: Vec<u64> = (0..1_000)
            .map(|i| (i * 12345) & ((1u64 << 33) - 1))
            .collect();
        assert_matches_packed_array(33, &values);
    }

    #[test]
    fn matches_packed_array_various_widths() {
        // Sweep widths to exercise both byte-aligned and unaligned
        // boundaries and check the carry logic in `pending`.
        //
        // **Width cap of 57**: `PackedArray::write` shifts the masked
        // value left by `bit_shift` (up to 7) inside a `u64`; for
        // `word_length + 7 > 64` (i.e. `word_length ≥ 58`) the shift
        // truncates the entry's MSBs. STAR's production widths are
        // 32-37 (`gstrand_bit + 1` for the genome size + the strand
        // bit), so this pre-existing `PackedArray` limitation is
        // never hit in practice. The streaming writer uses `u128`
        // and is correct for all widths 1..=64.
        use rand::{RngExt, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(0x5A_C0DE);
        for &wl in &[1u32, 5, 12, 32, 33, 35, 48, 57] {
            for &n in &[0usize, 1, 7, 8, 9, 127, 128, 129, 1000] {
                let mask = if wl == 64 { u64::MAX } else { (1u64 << wl) - 1 };
                let values: Vec<u64> = (0..n).map(|_| rng.random::<u64>() & mask).collect();
                assert_matches_packed_array(wl, &values);
            }
        }
    }

    #[test]
    fn empty_writer_produces_empty_file() {
        let mut got: Vec<u8> = Vec::new();
        let sw = PackedStreamWriter::new(&mut got, 33);
        let _w = sw.finish().unwrap();
        assert!(got.is_empty());
    }
}
