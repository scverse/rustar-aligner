/// Variable-width bit-packed array matching STAR's PackedArray format.
///
/// Stores integers with a specified bit width, packing them at bit-level
/// granularity (LSB-first, little-endian bit packing).
#[derive(Clone)]
pub struct PackedArray {
    /// Number of bits per element
    word_length: u32,

    /// Complement bits (64 - word_length) — kept for STAR compatibility
    #[allow(dead_code)]
    word_comp_length: u32,

    /// Mask for extracting an element (word_length bits set)
    bit_rec_mask: u64,

    /// Number of elements
    length: usize,

    /// Raw byte storage
    data: Vec<u8>,
}

impl PackedArray {
    /// Create a new PackedArray with specified bit width and length.
    ///
    /// # Arguments
    /// * `word_length` - Bits per element (1-64)
    /// * `length` - Number of elements
    pub fn new(word_length: u32, length: usize) -> Self {
        assert!(word_length > 0 && word_length <= 64);

        let word_comp_length = 64 - word_length;
        let bit_rec_mask = if word_length == 64 {
            u64::MAX
        } else {
            (1u64 << word_length) - 1
        };

        // Calculate bytes needed (matching STAR's formula)
        let length_byte = if length == 0 {
            0
        } else {
            ((length - 1) as u64 * word_length as u64) / 8 + 8
        };

        let data = vec![0u8; length_byte as usize];

        Self {
            word_length,
            word_comp_length,
            bit_rec_mask,
            length,
            data,
        }
    }

    /// Write an element at the specified index.
    ///
    /// # Arguments
    /// * `index` - Element index
    /// * `value` - Value to write (will be masked to word_length bits)
    pub fn write(&mut self, index: usize, value: u64) {
        assert!(index < self.length);

        let b = (index as u64 * self.word_length as u64) as usize; // bit offset
        let byte_offset = b / 8;
        let bit_shift = (b % 8) as u32;

        let masked_value = (value & self.bit_rec_mask) << bit_shift;
        let mask = self.bit_rec_mask << bit_shift;

        // Read current 8-byte word, update bits, write back
        let mut word = u64::from_le_bytes([
            self.data.get(byte_offset).copied().unwrap_or(0),
            self.data.get(byte_offset + 1).copied().unwrap_or(0),
            self.data.get(byte_offset + 2).copied().unwrap_or(0),
            self.data.get(byte_offset + 3).copied().unwrap_or(0),
            self.data.get(byte_offset + 4).copied().unwrap_or(0),
            self.data.get(byte_offset + 5).copied().unwrap_or(0),
            self.data.get(byte_offset + 6).copied().unwrap_or(0),
            self.data.get(byte_offset + 7).copied().unwrap_or(0),
        ]);

        word = (word & !mask) | masked_value;

        let bytes = word.to_le_bytes();
        for (i, &byte) in bytes.iter().enumerate() {
            if byte_offset + i < self.data.len() {
                self.data[byte_offset + i] = byte;
            }
        }
    }

    /// Read an element at the specified index.
    ///
    /// # Arguments
    /// * `index` - Element index
    ///
    /// # Returns
    /// The value at the specified index
    pub fn read(&self, index: usize) -> u64 {
        assert!(index < self.length);

        let b = (index as u64 * self.word_length as u64) as usize; // bit offset
        let byte_offset = b / 8;
        let bit_shift = (b % 8) as u32;

        let word = if byte_offset + 8 <= self.data.len() {
            // Fast path: read 8 bytes directly (no per-byte bounds checks)
            // SAFETY: We just verified byte_offset + 8 <= data.len()
            let bytes = &self.data[byte_offset..byte_offset + 8];
            u64::from_le_bytes(bytes.try_into().unwrap())
        } else {
            // Slow path: near end of array, read byte-by-byte with bounds checks
            u64::from_le_bytes([
                self.data.get(byte_offset).copied().unwrap_or(0),
                self.data.get(byte_offset + 1).copied().unwrap_or(0),
                self.data.get(byte_offset + 2).copied().unwrap_or(0),
                self.data.get(byte_offset + 3).copied().unwrap_or(0),
                self.data.get(byte_offset + 4).copied().unwrap_or(0),
                self.data.get(byte_offset + 5).copied().unwrap_or(0),
                self.data.get(byte_offset + 6).copied().unwrap_or(0),
                self.data.get(byte_offset + 7).copied().unwrap_or(0),
            ])
        };

        // Extract and mask the value
        (word >> bit_shift) & self.bit_rec_mask
    }

    /// Get the number of elements.
    pub fn len(&self) -> usize {
        self.length
    }

    /// Storage size in bytes for a [`PackedArray`] of `length` entries
    /// at the given bit width — equivalent to `PackedArray::new(...)
    /// .data().len()` but without allocating. Used by the streaming
    /// writer ([`PackedStreamWriter`]) to know when to stop emitting
    /// padding zeros.
    pub fn data_byte_len_for(word_length: u32, length: usize) -> usize {
        if length == 0 {
            0
        } else {
            // STAR's formula: `((length - 1) * word_length) / 8 + 8`.
            // The +8 reserves the 8-byte read window the `read` path
            // uses on the last entry.
            ((length as u64 - 1) * word_length as u64) as usize / 8 + 8
        }
    }

    /// Check if the array is empty.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Get the number of bits per element.
    pub fn word_length(&self) -> u32 {
        self.word_length
    }

    /// Get a reference to the raw byte data.
    pub fn data(&self) -> &[u8] {
        &self.data
    }

    /// Create a PackedArray from raw byte data.
    ///
    /// # Arguments
    /// * `word_length` - Bits per element
    /// * `length` - Number of elements
    /// * `data` - Raw byte data
    pub fn from_bytes(word_length: u32, length: usize, data: Vec<u8>) -> Self {
        assert!(word_length > 0 && word_length <= 64);

        let word_comp_length = 64 - word_length;
        let bit_rec_mask = if word_length == 64 {
            u64::MAX
        } else {
            (1u64 << word_length) - 1
        };

        Self {
            word_length,
            word_comp_length,
            bit_rec_mask,
            length,
            data,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_single_byte() {
        let mut arr = PackedArray::new(5, 10);

        arr.write(0, 31); // 11111 in 5 bits
        arr.write(1, 0); // 00000
        arr.write(2, 17); // 10001

        assert_eq!(arr.read(0), 31);
        assert_eq!(arr.read(1), 0);
        assert_eq!(arr.read(2), 17);
    }

    #[test]
    fn round_trip_cross_byte_boundary() {
        let mut arr = PackedArray::new(33, 100); // Human genome SA width

        let test_values = [
            0x0001_FFFF_FFFF, // All 33 bits set
            0x0001_0000_0000, // Bit 32 set (strand bit)
            0x0000_FFFF_FFFF, // Bits 0-31 set (max forward position)
            0,
            12_345_678,
        ];

        for (i, &val) in test_values.iter().enumerate() {
            arr.write(i, val);
        }

        for (i, &expected) in test_values.iter().enumerate() {
            assert_eq!(arr.read(i), expected);
        }
    }

    #[test]
    fn masking() {
        let mut arr = PackedArray::new(10, 5);

        // Write value larger than 10 bits — should be masked
        arr.write(0, 0xFFFF); // All bits set
        assert_eq!(arr.read(0), 0x3FF); // Only 10 bits = 1023
    }

    #[test]
    fn bit_width_32() {
        let mut arr = PackedArray::new(32, 10);

        arr.write(0, 0xDEAD_BEEF);
        arr.write(1, 0x1234_5678);
        arr.write(5, 0xCAFE_BABE);

        assert_eq!(arr.read(0), 0xDEAD_BEEF);
        assert_eq!(arr.read(1), 0x1234_5678);
        assert_eq!(arr.read(5), 0xCAFE_BABE);
    }

    #[test]
    fn sequential_writes() {
        let mut arr = PackedArray::new(7, 1000);

        for i in 0..1000 {
            arr.write(i, (i % 128) as u64);
        }

        for i in 0..1000 {
            assert_eq!(arr.read(i), (i % 128) as u64);
        }
    }
}
