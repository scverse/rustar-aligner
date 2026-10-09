//! Format-independent owned read types and transformations.

/// An owned read with bases in ruSTAR's internal encoding.
#[derive(Debug, Clone)]
pub struct EncodedRead {
    /// Read identifier.
    pub name: String,
    /// Base sequence encoded as 0=A, 1=C, 2=G, 3=T, 4=N.
    pub sequence: Vec<u8>,
    /// FASTQ ASCII quality bytes (Phred+33 encoded). Empty when the input format
    /// omitted qualities — nothing in the alignment, clipping or scoring path
    /// reads qualities, so no placeholder is materialized. SAM/BAM builders treat
    /// the empty slice as the format-level missing value (`*` / 0xff bytes).
    pub quality: Vec<u8>,
}

/// An owned paired-end read.
#[derive(Debug, Clone)]
pub struct PairedRead {
    /// Base read name (without a mate suffix).
    pub name: String,
    /// First mate in the pair.
    pub mate1: EncodedRead,
    /// Second mate in the pair.
    pub mate2: EncodedRead,
}

/// Strip a common mate suffix from a read name.
#[allow(clippy::case_sensitive_file_extension_comparisons)] // false positive
pub fn strip_mate_suffix(name: &str) -> String {
    let name = if let Some(pos) = name.find(' ') {
        &name[..pos]
    } else {
        name
    };

    if name.ends_with("/1") || name.ends_with("/2") {
        name[..name.len() - 2].to_string()
    } else if name.ends_with(".R1") || name.ends_with(".R2") {
        name[..name.len() - 3].to_string()
    } else if name.ends_with("_1") || name.ends_with("_2") {
        name[..name.len() - 2].to_string()
    } else {
        name.to_string()
    }
}

/// Convert an ASCII nucleotide to ruSTAR's genome encoding.
pub fn encode_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// Decode a ruSTAR-encoded nucleotide to ASCII.
pub fn decode_base(encoded: u8) -> u8 {
    match encoded {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

/// Complement an encoded base (A=0↔T=3, C=1↔G=2, N=4→N=4).
pub fn complement_base(encoded: u8) -> u8 {
    match encoded {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => encoded,
    }
}

/// Apply fixed 5' and 3' clipping to a sequence and its processing qualities.
///
/// An empty `qual` (an input format that omits qualities) clips to an empty
/// quality, preserving the SAM/BAM missing-quality sentinel through clipping.
pub fn clip_read(seq: &[u8], qual: &[u8], clip5p: usize, clip3p: usize) -> (Vec<u8>, Vec<u8>) {
    let len = seq.len();
    if clip5p + clip3p >= len {
        return (Vec::new(), Vec::new());
    }

    let start = clip5p;
    let end = len - clip3p;
    let clipped_qual = if qual.is_empty() {
        Vec::new()
    } else {
        qual[start..end].to_vec()
    };
    (seq[start..end].to_vec(), clipped_qual)
}
