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

/// Record parser used underneath [`FastqReader`].
///
/// Both backends read the same decompressed byte stream (plain file, gzip via
/// `flate2`, or `--readFilesCommand` output) and hand their records to the same
/// post-processing ([`encode_record`]), so the reads they produce are identical.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FastqBackend {
    /// `noodles` FASTQ parser, one record at a time (the default).
    #[default]
    Noodles,
    /// `paraseq` minimal-copy batch parser (issue #95). Needs the `paraseq`
    /// Cargo feature.
    #[cfg(feature = "paraseq")]
    Paraseq,
}

/// Environment variable that selects the FASTQ parser at run time.
pub const FASTQ_BACKEND_ENV: &str = "RUSTAR_FASTQ_BACKEND";

impl FastqBackend {
    /// Backend requested through `RUSTAR_FASTQ_BACKEND` (`noodles` or
    /// `paraseq`). Unset, empty or unknown values fall back to
    /// [`FastqBackend::Noodles`], with a warning for values this build cannot
    /// honour.
    #[must_use]
    pub fn from_env() -> Self {
        let Ok(value) = std::env::var(FASTQ_BACKEND_ENV) else {
            return Self::Noodles;
        };
        match value.trim().to_ascii_lowercase().as_str() {
            "" | "noodles" => Self::Noodles,
            #[cfg(feature = "paraseq")]
            "paraseq" => Self::Paraseq,
            #[cfg(not(feature = "paraseq"))]
            "paraseq" => {
                warn_once(&format!(
                    "{FASTQ_BACKEND_ENV}=paraseq ignored: this build lacks the `paraseq` feature; using noodles"
                ));
                Self::Noodles
            }
            other => {
                warn_once(&format!(
                    "{FASTQ_BACKEND_ENV}={other} not recognised (expected noodles or paraseq); using noodles"
                ));
                Self::Noodles
            }
        }
    }
}

fn warn_once(msg: &str) {
    static WARNED: std::sync::Once = std::sync::Once::new();
    WARNED.call_once(|| log::warn!("{msg}"));
}

enum RecordSource {
    Noodles(fastq::io::Reader<Box<dyn BufRead + Send>>),
    #[cfg(feature = "paraseq")]
    Paraseq(Box<ParaseqSource>),
}

/// Records parsed per `paraseq` fill; converted reads wait in `pending`.
#[cfg(feature = "paraseq")]
const PARASEQ_BATCH: usize = 4096;

#[cfg(feature = "paraseq")]
struct ParaseqSource {
    reader: paraseq::fastq::Reader<Box<dyn BufRead + Send>>,
    record_set: paraseq::fastq::RecordSet,
    pending: std::collections::VecDeque<EncodedRead>,
    done: bool,
}

#[cfg(feature = "paraseq")]
impl ParaseqSource {
    fn new(reader: Box<dyn BufRead + Send>) -> Self {
        Self {
            reader: paraseq::fastq::Reader::new(reader),
            record_set: paraseq::fastq::RecordSet::new(PARASEQ_BATCH),
            pending: std::collections::VecDeque::with_capacity(PARASEQ_BATCH),
            done: false,
        }
    }

    /// Make sure `pending` holds at least one read. Returns `false` once the
    /// input is exhausted.
    fn refill(&mut self, qual_shift: i32, name_separators: &[u8]) -> Result<bool, Error> {
        use paraseq::Record as _;

        if !self.pending.is_empty() {
            return Ok(true);
        }
        if self.done {
            return Ok(false);
        }
        let parse_err = |e: paraseq::Error| {
            Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("FASTQ parse error (paraseq): {e}"),
            ))
        };
        if !self.record_set.fill(&mut self.reader).map_err(parse_err)? {
            self.done = true;
            // paraseq leaves an incomplete trailing record in its overflow
            // buffer and reports a clean end of input; noodles errors on it.
            if !self.reader.exhausted() {
                return Err(Error::from(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "incomplete FASTQ record at end of input",
                )));
            }
            return Ok(false);
        }
        for record in self.record_set.iter() {
            let record = record.map_err(parse_err)?;
            // paraseq's id is the whole header line; noodles' name stops at
            // the first space or tab.
            let id = record.id();
            let name_end = id
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(id.len());
            let read = encode_record(
                qual_shift,
                name_separators,
                &id[..name_end],
                record.seq_raw(),
                record.qual().unwrap_or_default(),
            )?;
            self.pending.push_back(read);
        }
        Ok(true)
    }
}

/// Build an [`EncodedRead`] from one parsed record. Shared by every backend
/// so their output cannot drift apart. `name` is the header up to (not
/// including) the first space or tab.
fn encode_record(
    qual_shift: i32,
    name_separators: &[u8],
    name: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> Result<EncodedRead, Error> {
    let name = std::str::from_utf8(name).map_err(|e| {
        Error::from(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid UTF-8 in read name: {e}"),
        ))
    })?;

    let sequence = seq.iter().map(|&b| encode_base(b)).collect();

    let name = match name_separators
        .iter()
        .filter_map(|&sep| name.as_bytes().iter().position(|&b| b == sep))
        .min()
    {
        Some(cut) => name[..cut].to_string(),
        None => name.to_string(),
    };

    let quality = if qual_shift == 0 {
        qual.to_vec()
    } else {
        qual.iter()
            .map(|&b| (i32::from(b) + qual_shift).clamp(33, 126) as u8)
            .collect()
    };

    Ok(EncodedRead {
        name,
        sequence,
        quality,
    })
}

/// FASTQ reader that handles decompression and base encoding
pub struct FastqReader {
    inner: RecordSource,
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
        Self::open_with_backend(path, decompress_cmd, FastqBackend::from_env())
    }

    /// Like [`FastqReader::open`], with an explicit record parser instead of
    /// the one chosen by `RUSTAR_FASTQ_BACKEND`. Decompression is the same for
    /// every backend.
    pub fn open_with_backend(
        path: &Path,
        decompress_cmd: Option<&str>,
        backend: FastqBackend,
    ) -> Result<Self, Error> {
        let reader: Box<dyn BufRead + Send> = if let Some(cmd) = decompress_cmd {
            // Use external decompression command
            Self::open_with_command(path, cmd)?
        } else {
            // Auto-detect compression by file extension
            let path_str = path.to_string_lossy();
            let is_gzipped = path_str.ends_with(".gz") || path_str.ends_with(".gzip");

            let file = File::open(path).map_err(|e| Error::io(e, path))?;

            // Larger-than-default (8 KiB) buffers cut read syscalls on the decode
            // hot path. Feed the inflater from a big buffered file, and hand the
            // decoded stream to noodles through a big BufReader.
            const DECODE_BUF: usize = 1 << 19; // 512 KiB
            if is_gzipped {
                // Gzipped file
                let buffered = BufReader::with_capacity(DECODE_BUF, file);
                Box::new(BufReader::with_capacity(
                    DECODE_BUF,
                    GzDecoder::new(buffered),
                ))
            } else {
                // Plain text FASTQ
                Box::new(BufReader::with_capacity(DECODE_BUF, file))
            }
        };

        let inner = match backend {
            FastqBackend::Noodles => RecordSource::Noodles(fastq::io::Reader::new(reader)),
            #[cfg(feature = "paraseq")]
            FastqBackend::Paraseq => RecordSource::Paraseq(Box::new(ParaseqSource::new(reader))),
        };

        Ok(Self {
            inner,
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
        match &mut self.inner {
            RecordSource::Noodles(reader) => match reader.records().next() {
                Some(Ok(record)) => encode_record(
                    self.qual_shift,
                    &self.name_separators,
                    record.name(),
                    record.sequence(),
                    record.quality_scores(),
                )
                .map(Some),
                Some(Err(e)) => Err(Error::from(e)),
                None => Ok(None),
            },
            #[cfg(feature = "paraseq")]
            RecordSource::Paraseq(src) => {
                if src.refill(self.qual_shift, &self.name_separators)? {
                    Ok(src.pending.pop_front())
                } else {
                    Ok(None)
                }
            }
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
        #[cfg(feature = "paraseq")]
        if let RecordSource::Paraseq(src) = &mut self.inner {
            while batch.len() < batch_size && src.refill(self.qual_shift, &self.name_separators)? {
                let take = (batch_size - batch.len()).min(src.pending.len());
                batch.extend(src.pending.drain(..take));
            }
            return Ok(batch);
        }
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
        Self::open_with_backend(path1, path2, decompress_cmd, FastqBackend::from_env())
    }

    /// Like [`PairedFastqReader::open`], with an explicit record parser for
    /// both mates.
    pub fn open_with_backend(
        path1: &Path,
        path2: &Path,
        decompress_cmd: Option<&str>,
        backend: FastqBackend,
    ) -> Result<Self, Error> {
        let reader1 = FastqReader::open_with_backend(path1, decompress_cmd, backend)?;
        let reader2 = FastqReader::open_with_backend(path2, decompress_cmd, backend)?;

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

/// Equivalence of the `paraseq` backend against the default `noodles` one:
/// same names, sequences, qualities and PE pairing on identical input.
#[cfg(all(test, feature = "paraseq"))]
mod paraseq_equivalence {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    type Flat = (String, Vec<u8>, Vec<u8>);

    fn flat(r: &EncodedRead) -> Flat {
        (r.name.clone(), r.sequence.clone(), r.quality.clone())
    }

    /// Deterministic FASTQ text with awkward but valid headers: descriptions
    /// after a space or tab, `/1` suffixes, Illumina comments, lowercase and
    /// IUPAC bases, variable lengths, and a `+` line that repeats the name.
    fn fastq_text(n: usize, mate: u8, crlf: bool, trailing_newline: bool) -> String {
        let eol = if crlf { "\r\n" } else { "\n" };
        let mut s = String::new();
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15 ^ u64::from(mate);
        let bases = b"ACGTNacgtnRYK";
        for i in 0..n {
            let len = 20 + (i * 7) % 131;
            let mut seq = String::with_capacity(len);
            let mut qual = String::with_capacity(len);
            for _ in 0..len {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                seq.push(bases[(state % bases.len() as u64) as usize] as char);
                qual.push((b'!' + ((state >> 8) % 42) as u8) as char);
            }
            let header = match i % 4 {
                0 => format!("@r{i}/{mate}"),
                1 => format!("@r{i} {mate}:N:0:ACGT"),
                2 => format!("@r{i}\tdesc with spaces"),
                _ => format!("@r{i}"),
            };
            let plus = if i % 3 == 0 {
                format!("+{}", &header[1..])
            } else {
                "+".to_string()
            };
            s.push_str(&header);
            s.push_str(eol);
            s.push_str(&seq);
            s.push_str(eol);
            s.push_str(&plus);
            s.push_str(eol);
            s.push_str(&qual);
            if trailing_newline || i + 1 < n {
                s.push_str(eol);
            }
        }
        s
    }

    fn write_plain(text: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(text.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    fn write_gz(text: &str) -> NamedTempFile {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let f = tempfile::Builder::new()
            .suffix(".fq.gz")
            .tempfile()
            .unwrap();
        let mut enc = GzEncoder::new(f.reopen().unwrap(), Compression::default());
        enc.write_all(text.as_bytes()).unwrap();
        enc.finish().unwrap();
        f
    }

    fn read_all_se(path: &Path, backend: FastqBackend, batch: usize) -> Vec<Flat> {
        let mut reader = FastqReader::open_with_backend(path, None, backend).unwrap();
        let mut out = Vec::new();
        loop {
            let b = reader.read_batch(batch).unwrap();
            if b.is_empty() {
                break;
            }
            out.extend(b.iter().map(flat));
        }
        out
    }

    fn read_all_pe(p1: &Path, p2: &Path, backend: FastqBackend) -> Vec<(String, Flat, Flat)> {
        let mut reader = PairedFastqReader::open_with_backend(p1, p2, None, backend).unwrap();
        let mut out = Vec::new();
        loop {
            let b = reader.read_paired_batch(1000).unwrap();
            if b.is_empty() {
                break;
            }
            out.extend(
                b.iter()
                    .map(|p| (p.name.clone(), flat(&p.mate1), flat(&p.mate2))),
            );
        }
        out
    }

    #[test]
    fn se_plain_gzip_crlf_and_batch_boundaries() {
        // 10_007 records crosses several paraseq fills (4096) and read_batch
        // sizes that do not divide it.
        for (crlf, trailing) in [(false, true), (true, true), (false, false)] {
            let text = fastq_text(10_007, 1, crlf, trailing);
            for file in [write_plain(&text), write_gz(&text)] {
                let expected = read_all_se(file.path(), FastqBackend::Noodles, 10_000);
                assert_eq!(expected.len(), 10_007);
                for batch in [1, 999, 10_000] {
                    let got = read_all_se(file.path(), FastqBackend::Paraseq, batch);
                    assert_eq!(
                        got, expected,
                        "crlf={crlf} trailing={trailing} batch={batch}"
                    );
                }
            }
        }
    }

    #[test]
    fn se_next_encoded_matches() {
        let text = fastq_text(5000, 1, false, true);
        let file = write_plain(&text);
        let mut a =
            FastqReader::open_with_backend(file.path(), None, FastqBackend::Noodles).unwrap();
        let mut b =
            FastqReader::open_with_backend(file.path(), None, FastqBackend::Paraseq).unwrap();
        loop {
            let (x, y) = (a.next_encoded().unwrap(), b.next_encoded().unwrap());
            assert_eq!(x.as_ref().map(flat), y.as_ref().map(flat));
            if x.is_none() {
                break;
            }
        }
    }

    #[test]
    fn se_with_params_matches() {
        let text = fastq_text(3000, 1, false, true);
        let file = write_plain(&text);
        let mut params =
            crate::params::Parameters::parse_from(["rustar-aligner", "--readFilesIn", "x.fq"]);
        params.read_quality_score_base = 64;
        params.out_qs_conversion_add = 2;
        params.read_name_separator = vec!["_".to_string(), "/".to_string()];
        let collect = |backend| {
            let mut r = FastqReader::open_with_backend(file.path(), None, backend)
                .unwrap()
                .with_params(&params);
            r.read_batch(10_000)
                .unwrap()
                .iter()
                .map(flat)
                .collect::<Vec<_>>()
        };
        assert_eq!(
            collect(FastqBackend::Paraseq),
            collect(FastqBackend::Noodles)
        );
    }

    #[test]
    fn pe_plain_and_gzip_match() {
        let t1 = fastq_text(6001, 1, false, true);
        let t2 = fastq_text(6001, 2, false, true);
        let plain = (write_plain(&t1), write_plain(&t2));
        let gz = (write_gz(&t1), write_gz(&t2));
        for (f1, f2) in [&plain, &gz] {
            let expected = read_all_pe(f1.path(), f2.path(), FastqBackend::Noodles);
            assert_eq!(expected.len(), 6001);
            let got = read_all_pe(f1.path(), f2.path(), FastqBackend::Paraseq);
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn pe_length_mismatch_errors_on_both() {
        let f1 = write_plain(&fastq_text(10, 1, false, true));
        let f2 = write_plain(&fastq_text(9, 2, false, true));
        for backend in [FastqBackend::Noodles, FastqBackend::Paraseq] {
            let mut r =
                PairedFastqReader::open_with_backend(f1.path(), f2.path(), None, backend).unwrap();
            assert!(r.read_paired_batch(100).is_err(), "{backend:?}");
        }
    }

    #[test]
    fn truncated_record_errors_on_both() {
        let mut text = fastq_text(10, 1, false, true);
        text.push_str("@cut\nACGT");
        let file = write_plain(&text);
        for backend in [FastqBackend::Noodles, FastqBackend::Paraseq] {
            let mut r = FastqReader::open_with_backend(file.path(), None, backend).unwrap();
            assert!(r.read_batch(100).is_err(), "{backend:?}");
        }
    }
}
