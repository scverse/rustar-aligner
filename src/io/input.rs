//! Input format detection and layout planning.

use crate::error::Error;
use crate::io::cbq::{
    CbqProducer, CbqReadOpts, PairedEndConverter, SingleEndConverter, decoder_threads,
    paired_end_producer, single_end_producer,
};
use crate::io::fastq::{FastqReader, PairedFastqReader};
use crate::io::reads::{EncodedRead, PairedRead};
use crate::params::{OutReadsUnmapped, Parameters};
use binseq::write::Format;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::mpsc::SyncSender;

const MAGIC_PEEK_LEN: usize = 7;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadLayout {
    SingleEnd,
    PairedEnd,
}

/// A validated ordinary-alignment input configuration.
///
/// Enum variants encode the only supported path/layout combinations, so callers
/// cannot accidentally construct (for example) paired CBQ from two paths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadInputPlan {
    FastqSingle { path: PathBuf },
    FastqPaired { mate1: PathBuf, mate2: PathBuf },
    Cbq { path: PathBuf, layout: ReadLayout },
}

#[derive(Debug)]
enum DetectedInput {
    Fastq,
    Cbq {
        layout: ReadLayout,
        has_qualities: bool,
    },
    Unsupported(Format),
}

pub enum SingleEndProducer {
    Fastq(FastqReader),
    Cbq(CbqProducer<SingleEndConverter>),
}

pub enum PairedEndProducer {
    Fastq(PairedFastqReader),
    Cbq(CbqProducer<PairedEndConverter>),
}

impl ReadInputPlan {
    pub fn resolve(params: &Parameters) -> Result<Self, Error> {
        if params.read_files_in.is_empty() {
            return Err(Error::ReadInput(
                "no read files specified (--readFilesIn)".to_string(),
            ));
        }

        let detected = params
            .read_files_in
            .iter()
            .map(|path| detect_input(path))
            .collect::<Result<Vec<_>, _>>()?;

        for (path, format) in params.read_files_in.iter().zip(&detected) {
            if let DetectedInput::Unsupported(format) = format {
                return Err(Error::ReadInput(format!(
                    "unsupported BINSEQ format {} for {}; only CBQ is supported",
                    format_name(*format),
                    path.display()
                )));
            }
        }

        if (params.solo_enabled() || params.read_files_manifest.is_some())
            && detected
                .iter()
                .any(|input| !matches!(input, DetectedInput::Fastq))
        {
            return Err(Error::ReadInput(
                "CBQ is not supported with STARsolo or --readFilesManifest".to_string(),
            ));
        }

        match (params.read_files_in.as_slice(), detected.as_slice()) {
            ([path], [DetectedInput::Fastq]) => Ok(Self::FastqSingle { path: path.clone() }),
            (
                [path],
                [DetectedInput::Cbq {
                    layout,
                    has_qualities,
                    ..
                }],
            ) => {
                if params.read_files_command.is_some() {
                    return Err(Error::ReadInput(
                        "CBQ is incompatible with --readFilesCommand".to_string(),
                    ));
                }
                if !has_qualities && params.out_reads_unmapped == OutReadsUnmapped::Fastx {
                    return Err(Error::ReadInput(
                        "quality-less CBQ cannot be used with --outReadsUnmapped Fastx"
                            .to_string(),
                    ));
                }
                Ok(Self::Cbq {
                    path: path.clone(),
                    layout: *layout,
                })
            }
            ([mate1, mate2], [DetectedInput::Fastq, DetectedInput::Fastq]) => {
                Ok(Self::FastqPaired {
                    mate1: mate1.clone(),
                    mate2: mate2.clone(),
                })
            }
            ([_, _], _) => Err(Error::ReadInput(
                "paired CBQ is stored in one file; two CBQ paths and mixed FASTQ/CBQ inputs are unsupported"
                    .to_string(),
            )),
            (paths, _) => Err(Error::ReadInput(format!(
                "invalid number of read files: {} (expected 1 or 2)",
                paths.len()
            ))),
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Self::FastqSingle { .. } | Self::FastqPaired { .. } => "Fastq",
            Self::Cbq { .. } => "Cbq",
        }
    }

    pub fn layout(&self) -> ReadLayout {
        match self {
            Self::FastqSingle { .. }
            | Self::Cbq {
                layout: ReadLayout::SingleEnd,
                ..
            } => ReadLayout::SingleEnd,
            Self::FastqPaired { .. }
            | Self::Cbq {
                layout: ReadLayout::PairedEnd,
                ..
            } => ReadLayout::PairedEnd,
        }
    }

    pub fn open_single_end_producer(
        &self,
        params: &Parameters,
    ) -> Result<SingleEndProducer, Error> {
        match self {
            Self::FastqSingle { path } => Ok(SingleEndProducer::Fastq(
                FastqReader::open(path, params.read_files_command.as_deref())?.with_params(params),
            )),
            Self::Cbq {
                path,
                layout: ReadLayout::SingleEnd,
                ..
            } => Ok(SingleEndProducer::Cbq(single_end_producer(
                path.clone(),
                decoder_threads(params),
                CbqReadOpts::from_params(params),
            ))),
            Self::FastqPaired { .. }
            | Self::Cbq {
                layout: ReadLayout::PairedEnd,
                ..
            } => Err(Error::ReadInput(
                "paired input cannot be opened as single-end".to_string(),
            )),
        }
    }

    pub fn open_paired_end_producer(
        &self,
        params: &Parameters,
    ) -> Result<PairedEndProducer, Error> {
        match self {
            Self::FastqPaired { mate1, mate2 } => Ok(PairedEndProducer::Fastq(
                PairedFastqReader::open(mate1, mate2, params.read_files_command.as_deref())?
                    .with_params(params),
            )),
            Self::Cbq {
                path,
                layout: ReadLayout::PairedEnd,
                ..
            } => Ok(PairedEndProducer::Cbq(paired_end_producer(
                path.clone(),
                decoder_threads(params),
                CbqReadOpts::from_params(params),
            ))),
            Self::FastqSingle { .. }
            | Self::Cbq {
                layout: ReadLayout::SingleEnd,
                ..
            } => Err(Error::ReadInput(
                "single-end input cannot be opened as paired-end".to_string(),
            )),
        }
    }
}

impl SingleEndProducer {
    pub fn produce(
        self,
        batch_size: usize,
        max_records: usize,
        sender: &SyncSender<Result<Vec<EncodedRead>, Error>>,
    ) -> Result<(), Error> {
        match self {
            Self::Cbq(producer) => producer.produce(batch_size, max_records, sender),
            Self::Fastq(mut reader) => {
                produce_fastq(|n| reader.read_batch(n), batch_size, max_records, sender)
            }
        }
    }
}

impl PairedEndProducer {
    pub fn produce(
        self,
        batch_size: usize,
        max_records: usize,
        sender: &SyncSender<Result<Vec<PairedRead>, Error>>,
    ) -> Result<(), Error> {
        match self {
            Self::Cbq(producer) => producer.produce(batch_size, max_records, sender),
            Self::Fastq(mut reader) => produce_fastq(
                |n| reader.read_paired_batch(n),
                batch_size,
                max_records,
                sender,
            ),
        }
    }
}

/// Shared FASTQ producer loop: bounded batches, stopping at `max_records`.
fn produce_fastq<T: Send>(
    mut read_batch: impl FnMut(usize) -> Result<Vec<T>, Error>,
    batch_size: usize,
    max_records: usize,
    sender: &SyncSender<Result<Vec<T>, Error>>,
) -> Result<(), Error> {
    let mut produced = 0usize;
    loop {
        let remaining = max_records.saturating_sub(produced);
        if remaining == 0 {
            let _ = sender.send(Ok(Vec::new()));
            return Ok(());
        }
        let batch = read_batch(batch_size.min(remaining))?;
        let finished = batch.is_empty();
        produced += batch.len();
        if sender.send(Ok(batch)).is_err() || finished {
            return Ok(());
        }
    }
}

fn detect_input(path: &Path) -> Result<DetectedInput, Error> {
    let metadata = std::fs::metadata(path).map_err(|source| Error::io(source, path))?;
    if !metadata.is_file() {
        // Preserve the existing FASTQ behavior for named pipes and other streams.
        // CBQ itself is only recognized from regular files because it requires mmap.
        return Ok(DetectedInput::Fastq);
    }
    let mut file = File::open(path).map_err(|source| Error::io(source, path))?;
    let mut magic = [0u8; MAGIC_PEEK_LEN];
    let n = file
        .read(&mut magic)
        .map_err(|source| Error::io(source, path))?;

    match Format::sniff(&magic[..n]) {
        None => Ok(DetectedInput::Fastq),
        Some(Format::Bq) => Ok(DetectedInput::Unsupported(Format::Bq)),
        Some(Format::Vbq) => Ok(DetectedInput::Unsupported(Format::Vbq)),
        Some(Format::Cbq) => {
            let reader = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                binseq::cbq::MmapReader::new(path)
            }))
            .map_err(|_| Error::CbqWorkerPanic {
                path: path.to_path_buf(),
                start: 0,
                end: 0,
            })?
            .map_err(|source| Error::CbqOpen {
                path: path.to_path_buf(),
                source,
            })?;
            let header = reader.header();
            Ok(DetectedInput::Cbq {
                layout: if reader.is_paired() {
                    ReadLayout::PairedEnd
                } else {
                    ReadLayout::SingleEnd
                },
                has_qualities: header.has_qualities(),
            })
        }
    }
}

fn format_name(format: Format) -> &'static str {
    match format {
        Format::Bq => "BQ",
        Format::Vbq => "VBQ",
        Format::Cbq => "CBQ",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use binseq::SequencingRecordBuilder;
    use binseq::write::BinseqWriterBuilder;
    use std::ffi::OsString;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempPath};

    fn params_for(paths: &[&Path]) -> Parameters {
        let mut args = vec![
            OsString::from("rustar-aligner"),
            OsString::from("--readFilesIn"),
        ];
        args.extend(paths.iter().map(|path| path.as_os_str().to_os_string()));
        Parameters::parse_from(args)
    }

    fn write_cbq(paired: bool, qualities: bool, suffix: &str) -> TempPath {
        let temp = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        let path = temp.into_temp_path();
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(paired)
            .quality(qualities)
            .headers(true)
            .block_size(256)
            .build(File::create(&path).unwrap())
            .unwrap();
        let mut record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTN")
            .s_header(b"read/1");
        if qualities {
            record = record.s_qual(b"IIIII");
        }
        if paired {
            record = record.x_seq(b"TGCAN").x_header(b"read/2");
            if qualities {
                record = record.x_qual(b"IIIII");
            }
        }
        writer.push(record.build().unwrap()).unwrap();
        writer.finish().unwrap();
        path
    }

    #[test]
    fn detects_single_and_paired_cbq_by_magic_not_extension() {
        let single = write_cbq(false, true, ".data");
        let paired = write_cbq(true, true, "");

        assert!(matches!(
            ReadInputPlan::resolve(&params_for(&[&single])).unwrap(),
            ReadInputPlan::Cbq {
                layout: ReadLayout::SingleEnd,
                ..
            }
        ));
        assert!(matches!(
            ReadInputPlan::resolve(&params_for(&[&paired])).unwrap(),
            ReadInputPlan::Cbq {
                layout: ReadLayout::PairedEnd,
                ..
            }
        ));
    }

    #[test]
    fn fastq_named_cbq_remains_fastq() {
        let mut fastq = tempfile::Builder::new().suffix(".cbq").tempfile().unwrap();
        fastq.write_all(b"@r\nACGT\n+\nIIII\n").unwrap();
        let plan = ReadInputPlan::resolve(&params_for(&[fastq.path()])).unwrap();
        assert!(matches!(plan, ReadInputPlan::FastqSingle { .. }));
    }

    #[test]
    fn rejects_bq_and_vbq_magic_explicitly() {
        for (magic, expected) in [
            (binseq::bq::FILE_MAGIC.as_slice(), "BQ"),
            (binseq::vbq::FILE_MAGIC.as_slice(), "VBQ"),
        ] {
            let mut input = NamedTempFile::new().unwrap();
            input.write_all(magic).unwrap();
            let error = ReadInputPlan::resolve(&params_for(&[input.path()])).unwrap_err();
            assert!(error.to_string().contains(expected));
        }
    }

    #[test]
    fn recognized_corrupt_cbq_is_not_treated_as_fastq() {
        let mut input = NamedTempFile::new().unwrap();
        input.write_all(binseq::cbq::FILE_MAGIC).unwrap();
        let error = ReadInputPlan::resolve(&params_for(&[input.path()])).unwrap_err();
        assert!(error.to_string().contains("CBQ"));
    }

    #[test]
    fn rejects_mixed_and_two_path_cbq_inputs() {
        let cbq1 = write_cbq(false, true, ".cbq");
        let cbq2 = write_cbq(false, true, ".cbq");
        let mut fastq = NamedTempFile::new().unwrap();
        fastq.write_all(b"@r\nACGT\n+\nIIII\n").unwrap();

        let mixed = ReadInputPlan::resolve(&params_for(&[&cbq1, fastq.path()])).unwrap_err();
        assert!(mixed.to_string().contains("mixed FASTQ/CBQ"));
        let two = ReadInputPlan::resolve(&params_for(&[&cbq1, &cbq2])).unwrap_err();
        assert!(two.to_string().contains("one file"));
    }

    #[test]
    fn rejects_incompatible_cbq_options_upfront() {
        let qualityless = write_cbq(false, false, ".cbq");
        let mut params = params_for(&[&qualityless]);
        params.out_reads_unmapped = OutReadsUnmapped::Fastx;
        assert!(
            ReadInputPlan::resolve(&params)
                .unwrap_err()
                .to_string()
                .contains("quality-less")
        );

        let quality = write_cbq(false, true, ".cbq");
        let mut params = params_for(&[&quality]);
        params.read_files_command = Some("cat".to_string());
        assert!(
            ReadInputPlan::resolve(&params)
                .unwrap_err()
                .to_string()
                .contains("readFilesCommand")
        );

        let mut params = params_for(&[&quality]);
        params.solo_type = crate::params::SoloType::CbUmiSimple;
        assert!(
            ReadInputPlan::resolve(&params)
                .unwrap_err()
                .to_string()
                .contains("STARsolo")
        );
    }
}
