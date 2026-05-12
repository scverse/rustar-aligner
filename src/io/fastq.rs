/// FASTQ reader with base encoding and decompression support
use crate::error::Error;
use flate2::read::GzDecoder;
use noodles::fastq;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};

/// Writer for unmapped reads in FASTQ format (`--outReadsUnmapped Fastx`).
pub struct UnmappedFastqWriter {
    writer: BufWriter<File>,
}

impl UnmappedFastqWriter {
    pub fn create(path: &Path) -> Result<Self, Error> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| Error::io(e, parent))?;
        }
        let file = File::create(path).map_err(|e| Error::io(e, path))?;
        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    /// Write one FASTQ record. `seq` is in genome encoding (0=A,1=C,2=G,3=T,4=N).
    /// `qual` is raw FASTQ quality bytes.
    pub fn write_record(&mut self, name: &str, seq: &[u8], qual: &[u8]) -> Result<(), Error> {
        self.writer.write_all(b"@").map_err(Error::from)?;
        self.writer
            .write_all(name.as_bytes())
            .map_err(Error::from)?;
        self.writer.write_all(b"\n").map_err(Error::from)?;
        for &b in seq {
            self.writer
                .write_all(&[decode_base(b)])
                .map_err(Error::from)?;
        }
        self.writer.write_all(b"\n+\n").map_err(Error::from)?;
        self.writer.write_all(qual).map_err(Error::from)?;
        self.writer.write_all(b"\n").map_err(Error::from)
    }

    pub fn flush(&mut self) -> Result<(), Error> {
        self.writer.flush().map_err(Error::from)
    }
}

/// A read from a FASTQ file with encoded bases
#[derive(Debug, Clone)]
pub struct EncodedRead {
    /// Read identifier
    pub name: String,
    /// Base sequence encoded as 0=A, 1=C, 2=G, 3=T, 4=N
    pub sequence: Vec<u8>,
    /// Quality scores (raw FASTQ quality values)
    pub quality: Vec<u8>,
}

/// A paired-end read from two FASTQ files
#[derive(Debug, Clone)]
pub struct PairedRead {
    /// Base read name (without /1 or /2 suffix)
    pub name: String,
    /// First mate in pair
    pub mate1: EncodedRead,
    /// Second mate in pair
    pub mate2: EncodedRead,
}

/// FASTQ reader that handles decompression and base encoding
pub struct FastqReader {
    inner: fastq::Reader<Box<dyn BufRead + Send>>,
}

impl FastqReader {
    /// Open a FASTQ file (plain or gzip compressed)
    ///
    /// # Arguments
    /// * `path` - Path to FASTQ file
    /// * `decompress_cmd` - Optional decompression command (e.g., "zcat" for .gz files)
    ///
    /// # Returns
    /// A FastqReader that iterates over encoded reads
    pub fn open(path: &Path, decompress_cmd: Option<&str>) -> Result<Self, Error> {
        let reader: Box<dyn BufRead + Send> = if let Some(cmd) = decompress_cmd {
            // Use external decompression command
            Self::open_with_command(path, cmd)?
        } else {
            // Auto-detect compression by file extension
            let path_str = path.to_string_lossy();
            let is_gzipped = path_str.ends_with(".gz") || path_str.ends_with(".gzip");

            let file = File::open(path).map_err(|e| Error::io(e, path))?;

            if is_gzipped {
                // Gzipped file
                Box::new(BufReader::new(GzDecoder::new(file)))
            } else {
                // Plain text FASTQ
                Box::new(BufReader::new(file))
            }
        };

        let fastq_reader = fastq::Reader::new(reader);

        Ok(Self {
            inner: fastq_reader,
        })
    }

    /// Open FASTQ file using external decompression command
    fn open_with_command(path: &Path, cmd: &str) -> Result<Box<dyn BufRead + Send>, Error> {
        let mut child = Command::new(cmd)
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|e| Error::io(e, path))?;

        let stdout = child.stdout.take().ok_or_else(|| {
            Error::from(std::io::Error::other(
                "failed to capture stdout from decompression command",
            ))
        })?;

        Ok(Box::new(BufReader::new(stdout)))
    }

    /// Get next read with encoded bases
    pub fn next_encoded(&mut self) -> Result<Option<EncodedRead>, Error> {
        match self.inner.records().next() {
            Some(Ok(record)) => {
                let name = std::str::from_utf8(record.name())
                    .map_err(|e| {
                        Error::from(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("invalid UTF-8 in read name: {}", e),
                        ))
                    })?
                    .to_string();

                let sequence = record.sequence().iter().map(|&b| encode_base(b)).collect();

                let quality = record.quality_scores().to_vec();

                Ok(Some(EncodedRead {
                    name,
                    sequence,
                    quality,
                }))
            }
            Some(Err(e)) => Err(Error::from(e)),
            None => Ok(None),
        }
    }

    /// Read a batch of encoded reads for parallel processing
    ///
    /// # Arguments
    /// * `batch_size` - Maximum number of reads to return
    ///
    /// # Returns
    /// Vector of encoded reads (may be shorter than batch_size at end of file)
    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<EncodedRead>, Error> {
        let mut batch = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            match self.next_encoded()? {
                Some(read) => batch.push(read),
                None => break,
            }
        }
        Ok(batch)
    }
}

/// Paired-end FASTQ reader that reads from two files synchronously
pub struct PairedFastqReader {
    reader1: FastqReader,
    reader2: FastqReader,
}

impl PairedFastqReader {
    /// Open two FASTQ files for paired-end reading
    ///
    /// # Arguments
    /// * `path1` - Path to first mate FASTQ file
    /// * `path2` - Path to second mate FASTQ file
    /// * `decompress_cmd` - Optional decompression command
    ///
    /// # Returns
    /// A PairedFastqReader that iterates over paired reads with name validation
    pub fn open(path1: &Path, path2: &Path, decompress_cmd: Option<&str>) -> Result<Self, Error> {
        let reader1 = FastqReader::open(path1, decompress_cmd)?;
        let reader2 = FastqReader::open(path2, decompress_cmd)?;

        Ok(Self { reader1, reader2 })
    }

    /// Get next paired read with name validation
    ///
    /// # Returns
    /// - Ok(Some(PairedRead)) if both mates available and names match
    /// - Ok(None) if both files are exhausted
    /// - Err if only one file exhausted or names don't match
    pub fn next_paired(&mut self) -> Result<Option<PairedRead>, Error> {
        let read1_opt = self.reader1.next_encoded()?;
        let read2_opt = self.reader2.next_encoded()?;

        match (read1_opt, read2_opt) {
            (Some(read1), Some(read2)) => {
                // Strip mate suffixes for comparison
                let name1_base = strip_mate_suffix(&read1.name);
                let name2_base = strip_mate_suffix(&read2.name);

                if name1_base != name2_base {
                    return Err(Error::from(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "Paired FASTQ read names do not match: '{}' vs '{}'",
                            read1.name, read2.name
                        ),
                    )));
                }

                Ok(Some(PairedRead {
                    name: name1_base,
                    mate1: read1,
                    mate2: read2,
                }))
            }
            (None, None) => Ok(None),
            (Some(_), None) => Err(Error::from(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Paired FASTQ files have different lengths: mate1 file has more reads",
            ))),
            (None, Some(_)) => Err(Error::from(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Paired FASTQ files have different lengths: mate2 file has more reads",
            ))),
        }
    }

    /// Read a batch of paired reads for parallel processing
    ///
    /// # Arguments
    /// * `batch_size` - Maximum number of pairs to return
    ///
    /// # Returns
    /// Vector of paired reads (may be shorter than batch_size at end of file)
    pub fn read_paired_batch(&mut self, batch_size: usize) -> Result<Vec<PairedRead>, Error> {
        let mut batch = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            match self.next_paired()? {
                Some(paired) => batch.push(paired),
                None => break,
            }
        }
        Ok(batch)
    }
}

/// Strip mate suffix from read name for pairing
///
/// Removes common paired-end suffixes:
/// - /1 or /2 (Illumina convention)
/// - .R1 or .R2 (alternative convention)
/// - _1 or _2 (another convention)
/// - space and everything after (e.g., "READ_NAME 1:N:0:0" -> "READ_NAME")
///
/// # Arguments
/// * `name` - Original read name from FASTQ
///
/// # Returns
/// Base name with mate suffix removed
pub fn strip_mate_suffix(name: &str) -> String {
    // First, strip space and everything after (Illumina format)
    let name = if let Some(pos) = name.find(' ') {
        &name[..pos]
    } else {
        name
    };

    // Strip common mate suffixes
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

/// Convert FASTQ base character to genome encoding
///
/// # Arguments
/// * `base` - ASCII base character (A, C, G, T, N, or lowercase variants)
///
/// # Returns
/// Encoded base: 0=A, 1=C, 2=G, 3=T, 4=N (or any ambiguous base)
pub fn encode_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // N or any ambiguous base (R, Y, S, W, K, M, etc.)
    }
}

/// Decode genome encoding to ASCII base character
///
/// # Arguments
/// * `encoded` - Encoded base (0-4)
///
/// # Returns
/// ASCII base character (A, C, G, T, or N)
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
        0 => 3,       // A -> T
        1 => 2,       // C -> G
        2 => 1,       // G -> C
        3 => 0,       // T -> A
        _ => encoded, // N -> N
    }
}

/// Apply read clipping from 5' and 3' ends
///
/// # Arguments
/// * `seq` - Original sequence
/// * `qual` - Original quality scores
/// * `clip5p` - Number of bases to clip from 5' end
/// * `clip3p` - Number of bases to clip from 3' end
///
/// # Returns
/// Tuple of (clipped_sequence, clipped_quality)
pub fn clip_read(seq: &[u8], qual: &[u8], clip5p: usize, clip3p: usize) -> (Vec<u8>, Vec<u8>) {
    let len = seq.len();

    // Handle edge cases
    if clip5p + clip3p >= len {
        // Clipping removes entire read
        return (Vec::new(), Vec::new());
    }

    let start = clip5p;
    let end = len - clip3p;

    (seq[start..end].to_vec(), qual[start..end].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'a'), 0);
        assert_eq!(encode_base(b'C'), 1);
        assert_eq!(encode_base(b'c'), 1);
        assert_eq!(encode_base(b'G'), 2);
        assert_eq!(encode_base(b'g'), 2);
        assert_eq!(encode_base(b'T'), 3);
        assert_eq!(encode_base(b't'), 3);
        assert_eq!(encode_base(b'N'), 4);
        assert_eq!(encode_base(b'n'), 4);
        // Ambiguous bases
        assert_eq!(encode_base(b'R'), 4);
        assert_eq!(encode_base(b'Y'), 4);
        assert_eq!(encode_base(b'S'), 4);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), b'A');
        assert_eq!(decode_base(1), b'C');
        assert_eq!(decode_base(2), b'G');
        assert_eq!(decode_base(3), b'T');
        assert_eq!(decode_base(4), b'N');
        assert_eq!(decode_base(5), b'N'); // Invalid -> N
    }

    #[test]
    fn test_clip_read_none() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 0, 0);
        assert_eq!(clipped_seq, seq);
        assert_eq!(clipped_qual, qual);
    }

    #[test]
    fn test_clip_read_5p() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 2, 0);
        assert_eq!(clipped_seq, vec![2, 3, 0]); // GTA
        assert_eq!(clipped_qual, vec![30, 30, 30]);
    }

    #[test]
    fn test_clip_read_3p() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 0, 2);
        assert_eq!(clipped_seq, vec![0, 1, 2]); // ACG
        assert_eq!(clipped_qual, vec![30, 30, 30]);
    }

    #[test]
    fn test_clip_read_both() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 1, 1);
        assert_eq!(clipped_seq, vec![1, 2, 3]); // CGT
        assert_eq!(clipped_qual, vec![30, 30, 30]);
    }

    #[test]
    fn test_clip_read_entire() {
        let seq = vec![0, 1, 2, 3, 0]; // ACGTA
        let qual = vec![30, 30, 30, 30, 30];

        let (clipped_seq, clipped_qual) = clip_read(&seq, &qual, 3, 3);
        assert_eq!(clipped_seq, Vec::<u8>::new());
        assert_eq!(clipped_qual, Vec::<u8>::new());
    }

    #[test]
    fn test_fastq_reader_plain() {
        let mut tmpfile = NamedTempFile::new().unwrap();
        writeln!(tmpfile, "@read1").unwrap();
        writeln!(tmpfile, "ACGTN").unwrap();
        writeln!(tmpfile, "+").unwrap();
        writeln!(tmpfile, "IIIII").unwrap();
        writeln!(tmpfile, "@read2").unwrap();
        writeln!(tmpfile, "TGCA").unwrap();
        writeln!(tmpfile, "+").unwrap();
        writeln!(tmpfile, "HHHH").unwrap();
        tmpfile.flush().unwrap();

        let mut reader = FastqReader::open(tmpfile.path(), None).unwrap();

        let read1 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read1.name, "read1");
        assert_eq!(read1.sequence, vec![0, 1, 2, 3, 4]); // ACGTN
        assert_eq!(read1.quality.len(), 5);

        let read2 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read2.name, "read2");
        assert_eq!(read2.sequence, vec![3, 2, 1, 0]); // TGCA

        let read3 = reader.next_encoded().unwrap();
        assert!(read3.is_none());
    }

    #[test]
    fn test_fastq_reader_gzip() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let tmpfile = tempfile::Builder::new()
            .suffix(".fastq.gz")
            .tempfile()
            .unwrap();
        let mut encoder = GzEncoder::new(tmpfile.as_file(), Compression::default());
        writeln!(encoder, "@read1").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "IIII").unwrap();
        encoder.finish().unwrap();

        let mut reader = FastqReader::open(tmpfile.path(), None).unwrap();

        let read1 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read1.name, "read1");
        assert_eq!(read1.sequence, vec![0, 1, 2, 3]); // ACGT
        assert_eq!(read1.quality.len(), 4);
    }

    #[test]
    fn test_strip_mate_suffix_slash() {
        assert_eq!(strip_mate_suffix("read123/1"), "read123");
        assert_eq!(strip_mate_suffix("read123/2"), "read123");
    }

    #[test]
    fn test_strip_mate_suffix_dot() {
        assert_eq!(strip_mate_suffix("read123.R1"), "read123");
        assert_eq!(strip_mate_suffix("read123.R2"), "read123");
    }

    #[test]
    fn test_strip_mate_suffix_underscore() {
        assert_eq!(strip_mate_suffix("read123_1"), "read123");
        assert_eq!(strip_mate_suffix("read123_2"), "read123");
    }

    #[test]
    fn test_strip_mate_suffix_with_space() {
        assert_eq!(strip_mate_suffix("read123 1:N:0:AGCT"), "read123");
        assert_eq!(strip_mate_suffix("read123/1 1:N:0:AGCT"), "read123");
    }

    #[test]
    fn test_strip_mate_suffix_no_suffix() {
        assert_eq!(strip_mate_suffix("read123"), "read123");
    }

    #[test]
    fn test_paired_reader_matching_names() {
        let mut tmpfile1 = NamedTempFile::new().unwrap();
        writeln!(tmpfile1, "@read1/1").unwrap();
        writeln!(tmpfile1, "ACGT").unwrap();
        writeln!(tmpfile1, "+").unwrap();
        writeln!(tmpfile1, "IIII").unwrap();
        writeln!(tmpfile1, "@read2/1").unwrap();
        writeln!(tmpfile1, "TGCA").unwrap();
        writeln!(tmpfile1, "+").unwrap();
        writeln!(tmpfile1, "HHHH").unwrap();
        tmpfile1.flush().unwrap();

        let mut tmpfile2 = NamedTempFile::new().unwrap();
        writeln!(tmpfile2, "@read1/2").unwrap();
        writeln!(tmpfile2, "GGCC").unwrap();
        writeln!(tmpfile2, "+").unwrap();
        writeln!(tmpfile2, "JJJJ").unwrap();
        writeln!(tmpfile2, "@read2/2").unwrap();
        writeln!(tmpfile2, "AATT").unwrap();
        writeln!(tmpfile2, "+").unwrap();
        writeln!(tmpfile2, "KKKK").unwrap();
        tmpfile2.flush().unwrap();

        let mut reader = PairedFastqReader::open(tmpfile1.path(), tmpfile2.path(), None).unwrap();

        let pair1 = reader.next_paired().unwrap().unwrap();
        assert_eq!(pair1.name, "read1");
        assert_eq!(pair1.mate1.name, "read1/1");
        assert_eq!(pair1.mate1.sequence, vec![0, 1, 2, 3]); // ACGT
        assert_eq!(pair1.mate2.name, "read1/2");
        assert_eq!(pair1.mate2.sequence, vec![2, 2, 1, 1]); // GGCC

        let pair2 = reader.next_paired().unwrap().unwrap();
        assert_eq!(pair2.name, "read2");

        let pair3 = reader.next_paired().unwrap();
        assert!(pair3.is_none());
    }

    #[test]
    fn test_paired_reader_name_mismatch() {
        let mut tmpfile1 = NamedTempFile::new().unwrap();
        writeln!(tmpfile1, "@read1/1").unwrap();
        writeln!(tmpfile1, "ACGT").unwrap();
        writeln!(tmpfile1, "+").unwrap();
        writeln!(tmpfile1, "IIII").unwrap();
        tmpfile1.flush().unwrap();

        let mut tmpfile2 = NamedTempFile::new().unwrap();
        writeln!(tmpfile2, "@read2/2").unwrap();
        writeln!(tmpfile2, "GGCC").unwrap();
        writeln!(tmpfile2, "+").unwrap();
        writeln!(tmpfile2, "JJJJ").unwrap();
        tmpfile2.flush().unwrap();

        let mut reader = PairedFastqReader::open(tmpfile1.path(), tmpfile2.path(), None).unwrap();

        let result = reader.next_paired();
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("read names do not match")
        );
    }

    #[test]
    fn test_paired_reader_length_mismatch_mate1_longer() {
        let mut tmpfile1 = NamedTempFile::new().unwrap();
        writeln!(tmpfile1, "@read1/1").unwrap();
        writeln!(tmpfile1, "ACGT").unwrap();
        writeln!(tmpfile1, "+").unwrap();
        writeln!(tmpfile1, "IIII").unwrap();
        writeln!(tmpfile1, "@read2/1").unwrap();
        writeln!(tmpfile1, "TGCA").unwrap();
        writeln!(tmpfile1, "+").unwrap();
        writeln!(tmpfile1, "HHHH").unwrap();
        tmpfile1.flush().unwrap();

        let mut tmpfile2 = NamedTempFile::new().unwrap();
        writeln!(tmpfile2, "@read1/2").unwrap();
        writeln!(tmpfile2, "GGCC").unwrap();
        writeln!(tmpfile2, "+").unwrap();
        writeln!(tmpfile2, "JJJJ").unwrap();
        tmpfile2.flush().unwrap();

        let mut reader = PairedFastqReader::open(tmpfile1.path(), tmpfile2.path(), None).unwrap();

        // First pair succeeds
        let _ = reader.next_paired().unwrap().unwrap();

        // Second pair fails (mate1 has read but mate2 doesn't)
        let result = reader.next_paired();
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("different lengths")
        );
    }

    #[test]
    fn test_paired_reader_length_mismatch_mate2_longer() {
        let mut tmpfile1 = NamedTempFile::new().unwrap();
        writeln!(tmpfile1, "@read1/1").unwrap();
        writeln!(tmpfile1, "ACGT").unwrap();
        writeln!(tmpfile1, "+").unwrap();
        writeln!(tmpfile1, "IIII").unwrap();
        tmpfile1.flush().unwrap();

        let mut tmpfile2 = NamedTempFile::new().unwrap();
        writeln!(tmpfile2, "@read1/2").unwrap();
        writeln!(tmpfile2, "GGCC").unwrap();
        writeln!(tmpfile2, "+").unwrap();
        writeln!(tmpfile2, "JJJJ").unwrap();
        writeln!(tmpfile2, "@read2/2").unwrap();
        writeln!(tmpfile2, "AATT").unwrap();
        writeln!(tmpfile2, "+").unwrap();
        writeln!(tmpfile2, "KKKK").unwrap();
        tmpfile2.flush().unwrap();

        let mut reader = PairedFastqReader::open(tmpfile1.path(), tmpfile2.path(), None).unwrap();

        // First pair succeeds
        let _ = reader.next_paired().unwrap().unwrap();

        // Second pair fails (mate2 has read but mate1 doesn't)
        let result = reader.next_paired();
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("different lengths")
        );
    }

    #[test]
    fn test_paired_batch_reading() {
        let mut tmpfile1 = NamedTempFile::new().unwrap();
        for i in 1..=5 {
            writeln!(tmpfile1, "@read{}/1", i).unwrap();
            writeln!(tmpfile1, "ACGT").unwrap();
            writeln!(tmpfile1, "+").unwrap();
            writeln!(tmpfile1, "IIII").unwrap();
        }
        tmpfile1.flush().unwrap();

        let mut tmpfile2 = NamedTempFile::new().unwrap();
        for i in 1..=5 {
            writeln!(tmpfile2, "@read{}/2", i).unwrap();
            writeln!(tmpfile2, "GGCC").unwrap();
            writeln!(tmpfile2, "+").unwrap();
            writeln!(tmpfile2, "JJJJ").unwrap();
        }
        tmpfile2.flush().unwrap();

        let mut reader = PairedFastqReader::open(tmpfile1.path(), tmpfile2.path(), None).unwrap();

        // Read batch of 3
        let batch1 = reader.read_paired_batch(3).unwrap();
        assert_eq!(batch1.len(), 3);
        assert_eq!(batch1[0].name, "read1");
        assert_eq!(batch1[2].name, "read3");

        // Read remaining batch (should be 2)
        let batch2 = reader.read_paired_batch(3).unwrap();
        assert_eq!(batch2.len(), 2);
        assert_eq!(batch2[0].name, "read4");
        assert_eq!(batch2[1].name, "read5");

        // EOF batch
        let batch3 = reader.read_paired_batch(3).unwrap();
        assert_eq!(batch3.len(), 0);
    }
}
