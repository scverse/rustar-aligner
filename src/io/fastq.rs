/// FASTQ reader with base encoding and decompression support
use crate::error::Error;
use flate2::read::MultiGzDecoder;
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
    /// FASTQ ASCII quality bytes (Phred+33 encoded) - subtract 33 before
    /// writing to BAM binary QUAL.
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

/// Worker count handed to `rapidgzip-core` for gzip input, or 0 for the
/// single-threaded `flate2` path. Set once per run by
/// [`set_gz_decode_threads`]; read by every `FastqReader::open`.
///
/// A process-wide value rather than a parameter because `FastqReader::open` is
/// called from a dozen places (both mates, solo cDNA/barcode files, SmartSeq
/// manifests) that have no business each deciding a decompression policy.
#[cfg(feature = "rapidgzip")]
static GZ_DECODE_THREADS: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);

/// Decide how many threads gzip decode may use. **Off unless
/// `RUSTAR_GZ_DECODE_THREADS` asks for it**, because at this aligner's current
/// consumption rate parallel decode has nothing to win.
///
/// `rapidgzip-core` decodes a gzip stream in parallel with the marker/window
/// algorithm: blocks are decoded speculatively before their back-references are
/// known, then patched once they are. Its inflate backend is `libz-rs-sys`
/// (zlib-rs), the same one `flate2` uses here, so unlike an FFI decoder it costs
/// nothing at one thread. Decoding a 73 MB level-6 `.fq.gz` (427 MB out) on 16
/// logical cores, median of three:
///
/// | decoder | throughput |
/// |---------|------------|
/// | `flate2` + `zlib-rs`, 1 thread | 1552 MB/s |
/// | `rapidgzip-core`, 1 thread | 1581 MB/s |
/// | `rapidgzip-core`, 4 threads | 1872 MB/s |
/// | `rapidgzip-core`, 8 threads | 3356 MB/s |
/// | `rapidgzip-core`, 12 threads | 4051 MB/s |
///
/// The reason it is still off by default is arithmetic, not doubt about the
/// decoder. Aligning those 2 M reads takes ~3.0 s at `--runThreadN 8`, i.e. we
/// consume 427 MB of decompressed FASTQ at ~140 MB/s. One `flate2` thread
/// supplies ~1550 MB/s, so decode runs at under a tenth of its capacity and is
/// nowhere near the critical path; turning it on measured *slower* (3.15 s vs
/// 3.03 s) and ~360 MB more resident, which is what you would expect from
/// spending threads on a stage that was already idle.
///
/// It becomes the right call when the consumer gets fast enough to drain a
/// single inflate thread, which for this file shape means roughly 7 M reads/s.
/// That is the regime tools like piscem/salmon are in, and it is why they pair
/// a parallel decoder with a broker that splits one thread budget by measured
/// busy time on each side rather than by a fixed ratio.
#[cfg(feature = "rapidgzip")]
pub fn set_gz_decode_threads(_run_thread_n: u32) {
    let threads = std::env::var("RUSTAR_GZ_DECODE_THREADS")
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);
    GZ_DECODE_THREADS.store(threads, std::sync::atomic::Ordering::Relaxed);
}

/// No-op when the feature is off, so callers need no `cfg`.
#[cfg(not(feature = "rapidgzip"))]
pub fn set_gz_decode_threads(_run_thread_n: u32) {}

#[cfg(feature = "rapidgzip")]
fn gz_decode_threads() -> u32 {
    GZ_DECODE_THREADS.load(std::sync::atomic::Ordering::Relaxed)
}

/// FASTQ reader that handles decompression and base encoding
pub struct FastqReader {
    inner: fastq::io::Reader<Box<dyn BufRead + Send>>,
    /// Signed shift applied to every input quality byte so the rest of the
    /// pipeline always sees Phred+33.
    ///
    /// `--readQualityScoreBase 64` contributes `-31`, converting Solexa/Illumina
    /// 1.3 input. `--outQSconversionAdd` contributes its own value, which STAR
    /// documents as an output conversion; applying it here rather than at the
    /// writer means it also reaches anything else that reads qualities. That is
    /// only observable when both it and STARsolo are in use, which STAR itself
    /// does not support either.
    qual_shift: i32,
    /// Characters that terminate a read name (`--readNameSeparator`). Empty
    /// means keep the whole name.
    name_separators: Vec<u8>,
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

            // Larger-than-default (8 KiB) buffers cut read syscalls on the decode
            // hot path. Feed the inflater from a big buffered file, and hand the
            // decoded stream to noodles through a big BufReader.
            const DECODE_BUF: usize = 1 << 19; // 512 KiB
            if is_gzipped {
                if let Some(reader) = Self::open_gz_parallel(path) {
                    reader
                } else {
                    let file = File::open(path).map_err(|e| Error::io(e, path))?;
                    let buffered = BufReader::with_capacity(DECODE_BUF, file);
                    Box::new(BufReader::with_capacity(
                        DECODE_BUF,
                        MultiGzDecoder::new(buffered),
                    ))
                }
            } else {
                // Plain text FASTQ
                let file = File::open(path).map_err(|e| Error::io(e, path))?;
                Box::new(BufReader::with_capacity(DECODE_BUF, file))
            }
        };

        let fastq_reader = fastq::io::Reader::new(reader);

        Ok(Self {
            inner: fastq_reader,
            qual_shift: 0,
            name_separators: vec![b'/'],
        })
    }

    /// Apply the read-input knobs: `--readQualityScoreBase`,
    /// `--outQSconversionAdd` and `--readNameSeparator`.
    #[must_use]
    pub fn with_params(mut self, params: &crate::params::Parameters) -> Self {
        let base = if params.read_quality_score_base == 0 {
            33
        } else {
            params.read_quality_score_base
        };
        self.qual_shift = (33 - base) + params.out_qs_conversion_add;
        self.name_separators = params
            .read_name_separator
            .iter()
            .filter(|s| s.as_str() != "-")
            .filter_map(|s| s.as_bytes().first().copied())
            .collect();
        self
    }

    /// Parallel gzip/BGZF decode via `rapidgzip`, or `None` to use the
    /// single-threaded `flate2` path.
    ///
    /// Returns `None` unless the `rapidgzip` feature is compiled in **and**
    /// [`gz_decode_threads`] resolves to at least 2. Any failure to open the
    /// file through the native decoder also returns `None`, so a run degrades
    /// to the `flate2` path rather than failing: this is an optimisation, never
    /// a requirement.
    #[cfg(feature = "rapidgzip")]
    fn open_gz_parallel(path: &Path) -> Option<Box<dyn BufRead + Send>> {
        let threads = gz_decode_threads();
        if threads < 2 {
            return None;
        }
        // `Decoder::open` owns the file and returns a `Read + Send` stream, so it
        // drops straight into the `Box<dyn BufRead + Send>` the parser wants.
        // Reaching EOF verifies every member footer, which the `flate2` path
        // does too, so error behaviour on a corrupt file is not weakened.
        let decoder = match rapidgzip_core::Decoder::builder()
            .decoder_threads(threads as usize)
            .build()
        {
            Ok(d) => d,
            Err(e) => {
                log::warn!("rapidgzip-core could not be configured: {e}; using flate2");
                return None;
            }
        };
        match decoder.open(path) {
            Ok(reader) => {
                log::info!(
                    "decompressing {} with rapidgzip-core on {threads} threads",
                    path.display()
                );
                Some(Box::new(BufReader::with_capacity(1 << 19, reader)))
            }
            Err(e) => {
                log::warn!(
                    "rapidgzip-core could not open {}: {e}; falling back to flate2",
                    path.display()
                );
                None
            }
        }
    }

    #[cfg(not(feature = "rapidgzip"))]
    fn open_gz_parallel(_path: &Path) -> Option<Box<dyn BufRead + Send>> {
        None
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

        Ok(Box::new(BufReader::with_capacity(1 << 19, stdout)))
    }

    /// Get next read with encoded bases
    pub fn next_encoded(&mut self) -> Result<Option<EncodedRead>, Error> {
        match self.inner.records().next() {
            Some(Ok(record)) => {
                let name = std::str::from_utf8(record.name())
                    .map_err(|e| {
                        Error::from(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("invalid UTF-8 in read name: {e}"),
                        ))
                    })?
                    .to_string();

                let sequence = record.sequence().iter().map(|&b| encode_base(b)).collect();

                let name = match self
                    .name_separators
                    .iter()
                    .filter_map(|&sep| name.as_bytes().iter().position(|&b| b == sep))
                    .min()
                {
                    Some(cut) => name[..cut].to_string(),
                    None => name,
                };

                let quality = if self.qual_shift == 0 {
                    record.quality_scores().to_vec()
                } else {
                    record
                        .quality_scores()
                        .iter()
                        .map(|&b| (b as i32 + self.qual_shift).clamp(33, 126) as u8)
                        .collect()
                };

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

    /// Apply the read-input knobs to both mates. See
    /// [`FastqReader::with_params`].
    #[must_use]
    pub fn with_params(mut self, params: &crate::params::Parameters) -> Self {
        self.reader1 = self.reader1.with_params(params);
        self.reader2 = self.reader2.with_params(params);
        self
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
#[allow(clippy::case_sensitive_file_extension_comparisons)] // false positive
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

    /// A `.gz` written as several concatenated gzip members — what `bcl2fastq`
    /// emits, what `cat a.fq.gz b.fq.gz` produces, and what every BGZF file is.
    /// `flate2::read::GzDecoder` stops after the first member and reports EOF,
    /// so reading such a file used to drop reads with no error at all.
    #[test]
    fn test_fastq_reader_gzip_multi_member() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let mut tmpfile = tempfile::Builder::new()
            .suffix(".fastq.gz")
            .tempfile()
            .unwrap();

        // Member 1: read1. Each `finish()` closes a complete gzip stream, so
        // the next encoder appends a second member rather than continuing.
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        writeln!(encoder, "@read1").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "IIII").unwrap();
        tmpfile.write_all(&encoder.finish().unwrap()).unwrap();

        // Member 2: read2.
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        writeln!(encoder, "@read2").unwrap();
        writeln!(encoder, "TGCA").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "HHHH").unwrap();
        tmpfile.write_all(&encoder.finish().unwrap()).unwrap();
        tmpfile.flush().unwrap();

        let mut reader = FastqReader::open(tmpfile.path(), None).unwrap();

        let read1 = reader.next_encoded().unwrap().unwrap();
        assert_eq!(read1.name, "read1");
        assert_eq!(read1.sequence, vec![0, 1, 2, 3]); // ACGT

        let read2 = reader
            .next_encoded()
            .unwrap()
            .expect("second gzip member must be decoded, not silently truncated");
        assert_eq!(read2.name, "read2");
        assert_eq!(read2.sequence, vec![3, 2, 1, 0]); // TGCA

        assert!(reader.next_encoded().unwrap().is_none());
    }

    /// The parallel decoder must produce exactly what the `flate2` path
    /// produces, including across gzip member boundaries — the case that used
    /// to truncate. Reads the same file twice, once down each path.
    #[cfg(feature = "rapidgzip")]
    #[test]
    fn test_rapidgzip_decodes_identically_to_flate2() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::sync::atomic::Ordering;

        let mut tmpfile = tempfile::Builder::new()
            .suffix(".fastq.gz")
            .tempfile()
            .unwrap();
        // Two members, several reads each, so the parallel decoder has both a
        // member boundary and enough data to split on.
        for member in 0..2 {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            for i in 0..64 {
                writeln!(encoder, "@m{member}_read{i}").unwrap();
                writeln!(encoder, "ACGTACGTNN").unwrap();
                writeln!(encoder, "+").unwrap();
                writeln!(encoder, "IIIIIIIIII").unwrap();
            }
            tmpfile.write_all(&encoder.finish().unwrap()).unwrap();
        }
        tmpfile.flush().unwrap();

        let collect = |threads: u32| {
            GZ_DECODE_THREADS.store(threads, Ordering::Relaxed);
            let mut reader = FastqReader::open(tmpfile.path(), None).unwrap();
            let mut out = Vec::new();
            while let Some(read) = reader.next_encoded().unwrap() {
                out.push((read.name, read.sequence, read.quality));
            }
            out
        };

        let via_flate2 = collect(0);
        let via_rapidgzip = collect(4);
        GZ_DECODE_THREADS.store(0, Ordering::Relaxed);

        assert_eq!(via_flate2.len(), 128, "both members must be decoded");
        assert_eq!(via_rapidgzip, via_flate2);
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
        // Read names are cut at `/`, STAR's default --readNameSeparator, so
        // both mates report the shared name rather than the /1 and /2 forms.
        assert_eq!(pair1.mate1.name, "read1");
        assert_eq!(pair1.mate1.sequence, vec![0, 1, 2, 3]); // ACGT
        assert_eq!(pair1.mate2.name, "read1");
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
            writeln!(tmpfile1, "@read{i}/1").unwrap();
            writeln!(tmpfile1, "ACGT").unwrap();
            writeln!(tmpfile1, "+").unwrap();
            writeln!(tmpfile1, "IIII").unwrap();
        }
        tmpfile1.flush().unwrap();

        let mut tmpfile2 = NamedTempFile::new().unwrap();
        for i in 1..=5 {
            writeln!(tmpfile2, "@read{i}/2").unwrap();
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
