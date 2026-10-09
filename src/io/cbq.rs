//! Native CBQ decoding for the existing ordered alignment pipeline.

use crate::error::Error;
use crate::io::reads::{EncodedRead, PairedRead, encode_base, strip_mate_suffix};
use binseq::{BinseqRecord, ParallelProcessor, ParallelReader};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{SyncSender, sync_channel};

/// Hard ceiling on decoded records held for one CBQ window.
///
/// The window is the reorder buffer's memory limit (a smaller pending-map limit
/// could deadlock), so this constant — not the decoder-thread count and not the
/// input file's block size — is what bounds peak decoded-read memory. At ~330 B
/// per 150 bp single-end read that is ~10 MiB (~21 MiB for paired records).
///
/// A single CBQ block is the unavoidable floor: binseq decodes a whole block per
/// callback, so a file whose blocks each exceed the budget is bounded by one
/// block rather than by this value.
///
/// Windows must also be large enough to amortize a fixed per-window cost: binseq
/// loads the block that starts exactly at a window's exclusive end
/// (`cbq::MmapReader::process_parallel_range` selects blocks with
/// `iv_start <= range.end`), so every window decompresses one block it discards.
///
/// Selected by measurement, not assumption — decoder-only median wall time for
/// 200k x 150 bp SE records in 1 MiB quality-bearing blocks (~4651 records per
/// block), median of 5, AMD Ryzen 7 5700X:
///
/// ```text
/// window     1 thread   2 threads   4 threads   8 threads
///   9 456      73.4 ms     60.5 ms     62.4 ms     61.1 ms
///  16 384      68.3        42.9        44.9        44.1
///  32 768      61.4        43.0        29.9        29.6
///  49 152      58.4        43.5        33.0        28.7
///  65 536      57.6        41.7        33.6        27.3
/// 200 000      53.9        44.1        30.7        28.2
/// ```
///
/// Cost falls steeply to ~32k records and then flattens, so this sits at the knee:
/// best or statistically tied-best for the 1-4 decoder threads the automatic
/// policy can select, at half the peak memory of the next size up. Against the
/// previous `decoder_threads * 2` block policy on the same file: -20.1% at 1
/// thread, -6.3% at 2, -14.1% at 4, +1.9% at 8 (noise; 8 is not selectable
/// automatically).
const MAX_RECORDS_PER_WINDOW: usize = 32_768;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CbqWindow {
    record_start: usize,
    record_end: usize,
    block_count: usize,
}

#[derive(Debug)]
struct IndexedRecord<T> {
    index: usize,
    record: T,
}

#[derive(Debug)]
struct DecodedBlock<T> {
    first_index: usize,
    end_index: usize,
    records: Vec<IndexedRecord<T>>,
}

pub trait RecordConverter: Clone + Send + 'static {
    type Output: Send + 'static;

    fn convert<R: BinseqRecord>(&self, record: R, path: &Path) -> Result<Self::Output, Error>;
}

/// Read-input knobs that also make sense for CBQ input: the same
/// `--readNameSeparator` and `--outQSconversionAdd` handling the FASTQ reader
/// applies in `FastqReader::with_params`. (`--readQualityScoreBase` is a
/// FASTQ-encoding concept: CBQ qualities are stored as Phred+33.)
#[derive(Clone, Debug, Default)]
pub struct CbqReadOpts {
    name_separators: Vec<u8>,
    qual_shift: i32,
}

impl CbqReadOpts {
    pub fn from_params(params: &crate::params::Parameters) -> Self {
        Self {
            name_separators: params
                .read_name_separator
                .iter()
                .filter(|s| s.as_str() != "-")
                .filter_map(|s| s.as_bytes().first().copied())
                .collect(),
            qual_shift: params.out_qs_conversion_add,
        }
    }
}

#[derive(Clone, Default)]
pub struct SingleEndConverter {
    opts: CbqReadOpts,
}

impl RecordConverter for SingleEndConverter {
    type Output = EncodedRead;

    fn convert<R: BinseqRecord>(&self, record: R, path: &Path) -> Result<Self::Output, Error> {
        let index = record_index(&record, path)?;
        if record.is_paired() {
            return Err(record_error(
                path,
                index,
                None,
                "paired record found in single-end CBQ",
            ));
        }
        convert_mate(&record, path, index, 1, false, &self.opts)
    }
}

#[derive(Clone, Default)]
pub struct PairedEndConverter {
    opts: CbqReadOpts,
}

impl RecordConverter for PairedEndConverter {
    type Output = PairedRead;

    fn convert<R: BinseqRecord>(&self, record: R, path: &Path) -> Result<Self::Output, Error> {
        let index = record_index(&record, path)?;
        if !record.is_paired() {
            return Err(record_error(
                path,
                index,
                None,
                "single-end record found in paired CBQ",
            ));
        }

        let mate1 = convert_mate(&record, path, index, 1, false, &self.opts)?;
        let mate2 = convert_mate(&record, path, index, 2, true, &self.opts)?;
        let name1 = strip_mate_suffix(&mate1.name);
        let name2 = strip_mate_suffix(&mate2.name);
        if name1 != name2 {
            return Err(record_error(
                path,
                index,
                None,
                format!(
                    "paired read names do not match: '{}' vs '{}'",
                    mate1.name, mate2.name
                ),
            ));
        }

        Ok(PairedRead {
            name: name1,
            mate1,
            mate2,
        })
    }
}

struct BlockProcessor<C: RecordConverter> {
    converter: C,
    path: PathBuf,
    current: Vec<IndexedRecord<C::Output>>,
    sender: SyncSender<DecodedBlock<C::Output>>,
}

impl<C: RecordConverter> Clone for BlockProcessor<C> {
    fn clone(&self) -> Self {
        Self {
            converter: self.converter.clone(),
            path: self.path.clone(),
            current: Vec::new(),
            sender: self.sender.clone(),
        }
    }
}

impl<C: RecordConverter> BlockProcessor<C> {
    fn flush_block(&mut self) -> binseq::Result<()> {
        if self.current.is_empty() {
            return Ok(());
        }

        let first_index = self.current[0].index;
        let mut expected = first_index;
        for indexed in &self.current {
            if indexed.index != expected {
                return Err(into_binseq_error(Error::CbqOrdering {
                    path: self.path.clone(),
                    expected_index: expected,
                    observed_index: indexed.index,
                }));
            }
            expected = expected.checked_add(1).ok_or_else(|| {
                into_binseq_error(Error::ReadInput("CBQ record index overflow".to_string()))
            })?;
        }

        let records = std::mem::take(&mut self.current);
        self.sender
            .send(DecodedBlock {
                first_index,
                end_index: expected,
                records,
            })
            .map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::BrokenPipe,
                    "CBQ block coordinator disconnected",
                )
                .into()
            })
    }
}

impl<C: RecordConverter> ParallelProcessor for BlockProcessor<C> {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        let index = record_index(&record, &self.path).map_err(into_binseq_error)?;
        let owned = self
            .converter
            .convert(record, &self.path)
            .map_err(into_binseq_error)?;
        self.current.push(IndexedRecord {
            index,
            record: owned,
        });
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.flush_block()
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.flush_block()
    }
}

pub struct CbqProducer<C: RecordConverter> {
    path: PathBuf,
    converter: C,
    decoder_threads: usize,
    /// Records per decode window; always `MAX_RECORDS_PER_WINDOW` in production.
    /// Tests lower it to reach multi-window decoding on small fixtures.
    window_records: usize,
}

pub fn single_end_producer(
    path: PathBuf,
    decoder_threads: usize,
    opts: CbqReadOpts,
) -> CbqProducer<SingleEndConverter> {
    CbqProducer {
        path,
        converter: SingleEndConverter { opts },
        decoder_threads,
        window_records: MAX_RECORDS_PER_WINDOW,
    }
}

pub fn paired_end_producer(
    path: PathBuf,
    decoder_threads: usize,
    opts: CbqReadOpts,
) -> CbqProducer<PairedEndConverter> {
    CbqProducer {
        path,
        converter: PairedEndConverter { opts },
        decoder_threads,
        window_records: MAX_RECORDS_PER_WINDOW,
    }
}

impl<C: RecordConverter> CbqProducer<C> {
    pub fn produce(
        self,
        batch_size: usize,
        max_records: usize,
        sender: &SyncSender<Result<Vec<C::Output>, Error>>,
    ) -> Result<(), Error> {
        let reader = open_reader(&self.path)?;
        let total_to_decode = reader.num_records().min(max_records);
        if total_to_decode == 0 {
            let _ = sender.send(Ok(Vec::new()));
            return Ok(());
        }

        let block_record_counts = read_block_record_counts(&reader, &self.path)?;
        let windows = plan_windows(&block_record_counts, total_to_decode, self.window_records)?;
        log::info!(
            "CBQ input {}: {} records, {} decoder thread(s), {} bounded window(s) of \u{2264}{} records",
            self.path.display(),
            total_to_decode,
            self.decoder_threads,
            windows.len(),
            windows
                .iter()
                .map(|window| window.record_end - window.record_start)
                .max()
                .unwrap_or(0)
        );
        let mut output_batch = Vec::with_capacity(batch_size);

        for window in windows {
            if !decode_window(
                reader.clone(),
                &self.path,
                self.converter.clone(),
                window,
                self.decoder_threads.min(window.block_count).max(1),
                batch_size,
                &mut output_batch,
                sender,
            )? {
                return Ok(());
            }
        }

        if !output_batch.is_empty()
            && sender
                .send(Ok(std::mem::replace(
                    &mut output_batch,
                    Vec::with_capacity(batch_size),
                )))
                .is_err()
        {
            return Ok(());
        }
        let _ = sender.send(Ok(Vec::new()));
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn decode_window<C: RecordConverter>(
    reader: binseq::cbq::MmapReader,
    path: &Path,
    converter: C,
    window: CbqWindow,
    decoder_threads: usize,
    batch_size: usize,
    output_batch: &mut Vec<C::Output>,
    output_sender: &SyncSender<Result<Vec<C::Output>, Error>>,
) -> Result<bool, Error> {
    let (block_sender, block_receiver) = sync_channel(decoder_threads.max(1));
    let processor = BlockProcessor {
        converter,
        path: path.to_path_buf(),
        current: Vec::new(),
        sender: block_sender,
    };
    let range = window.record_start..window.record_end;
    let runner = std::thread::spawn(move || {
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            reader.process_parallel_range(processor, decoder_threads, range)
        }))
    });

    let mut coordinator = OrderedBlocks::new(path.to_path_buf(), window.record_start);
    let mut coordinator_error = None;
    let mut cancelled = false;

    while let Ok(block) = block_receiver.recv() {
        if coordinator_error.is_some() || cancelled {
            continue;
        }
        match coordinator.push(block) {
            Ok(ready) => {
                for record in ready {
                    output_batch.push(record);
                    if output_batch.len() == batch_size {
                        let full = std::mem::replace(output_batch, Vec::with_capacity(batch_size));
                        if output_sender.send(Ok(full)).is_err() {
                            cancelled = true;
                            break;
                        }
                    }
                }
            }
            Err(error) => coordinator_error = Some(error),
        }
    }

    let runner_result = runner.join().map_err(|_| Error::CbqWorkerPanic {
        path: path.to_path_buf(),
        start: window.record_start,
        end: window.record_end,
    })?;

    if cancelled {
        return Ok(false);
    }
    if let Some(error) = coordinator_error {
        return Err(error);
    }
    match runner_result {
        Err(_) => Err(Error::CbqWorkerPanic {
            path: path.to_path_buf(),
            start: window.record_start,
            end: window.record_end,
        }),
        Ok(Err(source)) => Err(Error::CbqDecode {
            path: path.to_path_buf(),
            start: window.record_start,
            end: window.record_end,
            source,
        }),
        Ok(Ok(())) => {
            coordinator.finish(window.record_end)?;
            Ok(true)
        }
    }
}

struct OrderedBlocks<T> {
    path: PathBuf,
    next_record_index: usize,
    pending: BTreeMap<usize, DecodedBlock<T>>,
}

impl<T> OrderedBlocks<T> {
    fn new(path: PathBuf, next_record_index: usize) -> Self {
        Self {
            path,
            next_record_index,
            pending: BTreeMap::new(),
        }
    }

    fn push(&mut self, block: DecodedBlock<T>) -> Result<Vec<T>, Error> {
        if block.records.is_empty() || block.first_index >= block.end_index {
            return Err(Error::CbqOrdering {
                path: self.path.clone(),
                expected_index: self.next_record_index,
                observed_index: block.first_index,
            });
        }
        let expected_end = block
            .first_index
            .checked_add(block.records.len())
            .ok_or_else(|| Error::ReadInput("CBQ block index overflow".to_string()))?;
        if expected_end != block.end_index {
            return Err(Error::CbqOrdering {
                path: self.path.clone(),
                expected_index: expected_end,
                observed_index: block.end_index,
            });
        }
        if block.first_index < self.next_record_index
            || self.pending.contains_key(&block.first_index)
        {
            return Err(Error::CbqOrdering {
                path: self.path.clone(),
                expected_index: self.next_record_index,
                observed_index: block.first_index,
            });
        }

        self.pending.insert(block.first_index, block);
        let mut ready = Vec::new();
        while let Some(block) = self.pending.remove(&self.next_record_index) {
            for indexed in block.records {
                if indexed.index != self.next_record_index {
                    return Err(Error::CbqOrdering {
                        path: self.path.clone(),
                        expected_index: self.next_record_index,
                        observed_index: indexed.index,
                    });
                }
                ready.push(indexed.record);
                self.next_record_index += 1;
            }
        }
        Ok(ready)
    }

    fn finish(&self, expected_end: usize) -> Result<(), Error> {
        if self.next_record_index != expected_end || !self.pending.is_empty() {
            let observed = self
                .pending
                .first_key_value()
                .map_or(self.next_record_index, |(index, _)| *index);
            return Err(Error::CbqOrdering {
                path: self.path.clone(),
                expected_index: expected_end,
                observed_index: observed,
            });
        }
        Ok(())
    }
}

fn open_reader(path: &Path) -> Result<binseq::cbq::MmapReader, Error> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
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
    })
}

fn read_block_record_counts(
    reader: &binseq::cbq::MmapReader,
    path: &Path,
) -> Result<Vec<usize>, Error> {
    let headers = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        reader
            .iter_block_headers()
            .collect::<binseq::Result<Vec<_>>>()
    }))
    .map_err(|_| Error::CbqWorkerPanic {
        path: path.to_path_buf(),
        start: 0,
        end: reader.num_records(),
    })?
    .map_err(|source| Error::CbqOpen {
        path: path.to_path_buf(),
        source,
    })?;

    let counts: Vec<usize> = headers
        .into_iter()
        .map(|header| {
            usize::try_from(header.num_records).map_err(|_| {
                Error::ReadInput(format!(
                    "CBQ block record count does not fit usize in {}",
                    path.display()
                ))
            })
        })
        .collect::<Result<_, _>>()?;
    let indexed_records = counts.iter().try_fold(0usize, |total, count| {
        total
            .checked_add(*count)
            .ok_or_else(|| Error::ReadInput("CBQ record-count overflow".to_string()))
    })?;
    if indexed_records != reader.num_records() {
        return Err(Error::ReadInput(format!(
            "CBQ block headers contain {indexed_records} records but the index advertises {} in {}",
            reader.num_records(),
            path.display()
        )));
    }
    Ok(counts)
}

/// Split the CBQ block index into bounded, block-aligned decode windows.
///
/// `max_records_per_window` is a hard ceiling with one exception: the first block
/// of a window is always admitted, so an oversized block still makes progress
/// instead of stalling the plan.
fn plan_windows(
    block_record_counts: &[usize],
    total_records: usize,
    max_records_per_window: usize,
) -> Result<Vec<CbqWindow>, Error> {
    if total_records == 0 {
        return Ok(Vec::new());
    }
    let budget = max_records_per_window.max(1);
    let mut windows = Vec::new();
    let mut block_index = 0usize;
    let mut record_start = 0usize;

    while record_start < total_records {
        if block_index >= block_record_counts.len() {
            return Err(Error::ReadInput(format!(
                "CBQ index contains fewer records than advertised ({record_start} < {total_records})"
            )));
        }
        let mut block_count = 0usize;
        let mut boundary_end = record_start;

        while block_index < block_record_counts.len() {
            let block_records = block_record_counts[block_index];
            if block_count > 0
                && boundary_end
                    .saturating_sub(record_start)
                    .saturating_add(block_records)
                    > budget
            {
                break;
            }
            boundary_end = boundary_end
                .checked_add(block_records)
                .ok_or_else(|| Error::ReadInput("CBQ record-count overflow".to_string()))?;
            block_index += 1;
            block_count += 1;
            if boundary_end >= total_records {
                break;
            }
        }

        let record_end = boundary_end.min(total_records);
        if record_end <= record_start {
            return Err(Error::ReadInput(
                "CBQ index contains an empty or non-advancing block window".to_string(),
            ));
        }
        windows.push(CbqWindow {
            record_start,
            record_end,
            block_count,
        });
        record_start = record_end;
        if record_end < boundary_end {
            break;
        }
    }

    Ok(windows)
}

fn record_index<R: BinseqRecord>(record: &R, path: &Path) -> Result<usize, Error> {
    usize::try_from(record.index())
        .map_err(|_| record_error(path, usize::MAX, None, "record index does not fit usize"))
}

fn convert_mate<R: BinseqRecord>(
    record: &R,
    path: &Path,
    index: usize,
    mate: u8,
    extended: bool,
    opts: &CbqReadOpts,
) -> Result<EncodedRead, Error> {
    let (header, sequence, quality) = if extended {
        (record.xheader(), record.xseq(), record.xqual())
    } else {
        (record.sheader(), record.sseq(), record.squal())
    };

    let name = std::str::from_utf8(header)
        .map_err(|source| {
            record_error(
                path,
                index,
                Some(mate),
                format!("invalid UTF-8 in read name: {source}"),
            )
        })?
        .split(|c: char| c.is_ascii_whitespace())
        .next()
        .unwrap_or_default()
        .to_string();
    let name = match opts
        .name_separators
        .iter()
        .filter_map(|&sep| name.as_bytes().iter().position(|&b| b == sep))
        .min()
    {
        Some(cut) => name[..cut].to_string(),
        None => name,
    };
    let encoded_sequence = sequence.iter().map(|&base| encode_base(base)).collect();
    let encoded_quality = if !record.has_quality() {
        // No placeholder: qualities are never read by alignment, clipping or
        // scoring, and `clip_read` preserves an empty quality through clipping.
        Vec::new()
    } else {
        if quality.len() != sequence.len() {
            return Err(record_error(
                path,
                index,
                Some(mate),
                format!(
                    "sequence/quality length mismatch: {} bases, {} qualities",
                    sequence.len(),
                    quality.len()
                ),
            ));
        }
        if opts.qual_shift == 0 {
            quality.to_vec()
        } else {
            quality
                .iter()
                .map(|&b| (i32::from(b) + opts.qual_shift).clamp(33, 126) as u8)
                .collect()
        }
    };

    Ok(EncodedRead {
        name,
        sequence: encoded_sequence,
        quality: encoded_quality,
    })
}

fn record_error(
    path: &Path,
    record_index: usize,
    mate: Option<u8>,
    message: impl Into<String>,
) -> Error {
    Error::CbqRecord {
        path: path.to_path_buf(),
        record_index,
        mate_context: mate.map_or_else(String::new, |mate| format!(" mate {mate}")),
        message: message.into(),
    }
}

fn into_binseq_error(error: Error) -> binseq::Error {
    binseq::Error::GenericError(Box::new(error))
}

/// Resolve CBQ decoder concurrency from CLI values.
pub fn decoder_threads(params: &crate::params::Parameters) -> usize {
    let requested = if params.read_files_n_threads > 0 {
        params.read_files_n_threads
    } else {
        (usize::from(params.run_thread_n) / 4).clamp(1, 4)
    };
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    requested.clamp(1, available)
}

#[cfg(test)]
mod tests {
    use super::*;
    use binseq::SequencingRecordBuilder;
    use binseq::write::{BinseqWriterBuilder, Format};
    use std::fs::{File, OpenOptions};
    use std::io::{Read, Seek, SeekFrom, Write};
    use tempfile::{NamedTempFile, TempPath};

    fn write_cbq(paired: bool, qualities: bool, headers: bool, n_records: usize) -> TempPath {
        let temp = NamedTempFile::new().unwrap();
        let path = temp.into_temp_path();
        let file = File::create(&path).unwrap();
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(paired)
            .quality(qualities)
            .headers(headers)
            .block_size(256)
            .build(file)
            .unwrap();

        for i in 0..n_records {
            let h1 = format!("read{i}/1");
            let h2 = format!("read{i}/2");
            let mut builder = SequencingRecordBuilder::default().s_seq(b"ACGTNACGTN");
            if headers {
                builder = builder.s_header(h1.as_bytes());
            }
            if qualities {
                builder = builder.s_qual(b"ABCDEFGHIJ");
            }
            if paired {
                builder = builder.x_seq(b"TGCANTGCAN");
                if headers {
                    builder = builder.x_header(h2.as_bytes());
                }
                if qualities {
                    builder = builder.x_qual(b"JKLMNOPQRS");
                }
            }
            writer.push(builder.build().unwrap()).unwrap();
        }
        writer.finish().unwrap();
        path
    }

    fn collect_single(path: &Path, threads: usize, max_records: usize) -> Vec<EncodedRead> {
        collect_single_windowed(path, threads, max_records, MAX_RECORDS_PER_WINDOW)
    }

    fn collect_single_windowed(
        path: &Path,
        threads: usize,
        max_records: usize,
        window_records: usize,
    ) -> Vec<EncodedRead> {
        let (sender, receiver) = sync_channel(128);
        CbqProducer {
            path: path.to_path_buf(),
            converter: SingleEndConverter::default(),
            decoder_threads: threads,
            window_records,
        }
        .produce(7, max_records, &sender)
        .unwrap();
        drop(sender);
        receiver
            .into_iter()
            .map(Result::unwrap)
            .take_while(|batch| !batch.is_empty())
            .flatten()
            .collect()
    }

    fn collect_paired(path: &Path, threads: usize) -> Vec<PairedRead> {
        collect_paired_windowed(path, threads, MAX_RECORDS_PER_WINDOW)
    }

    fn collect_paired_windowed(
        path: &Path,
        threads: usize,
        window_records: usize,
    ) -> Vec<PairedRead> {
        let (sender, receiver) = sync_channel(128);
        CbqProducer {
            path: path.to_path_buf(),
            converter: PairedEndConverter::default(),
            decoder_threads: threads,
            window_records,
        }
        .produce(9, usize::MAX, &sender)
        .unwrap();
        drop(sender);
        receiver
            .into_iter()
            .map(Result::unwrap)
            .take_while(|batch| !batch.is_empty())
            .flatten()
            .collect()
    }

    fn write_custom_headers(header1: &[u8], header2: Option<&[u8]>) -> TempPath {
        let temp = NamedTempFile::new().unwrap();
        let path = temp.into_temp_path();
        let paired = header2.is_some();
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(paired)
            .quality(true)
            .headers(true)
            .block_size(256)
            .build(File::create(&path).unwrap())
            .unwrap();
        let mut record = SequencingRecordBuilder::default()
            .s_seq(b"ACGT")
            .s_qual(b"IIII")
            .s_header(header1);
        if let Some(header2) = header2 {
            record = record.x_seq(b"TGCA").x_qual(b"IIII").x_header(header2);
        }
        writer.push(record.build().unwrap()).unwrap();
        writer.finish().unwrap();
        path
    }

    /// Every window must tile `0..total` exactly once, in order, at block
    /// boundaries, and honor the record budget (one oversized block excepted).
    fn assert_windows_tile(windows: &[CbqWindow], counts: &[usize], total: usize, budget: usize) {
        let mut expected_start = 0usize;
        let mut consumed_blocks = 0usize;
        for window in windows {
            assert_eq!(
                window.record_start, expected_start,
                "windows must be gapless"
            );
            assert!(window.record_end > window.record_start, "empty window");
            let records = window.record_end - window.record_start;
            let largest = counts[consumed_blocks..consumed_blocks + window.block_count]
                .iter()
                .copied()
                .max()
                .unwrap_or(0);
            assert!(
                records <= budget.max(largest),
                "window of {records} records exceeds max(budget {budget}, largest block {largest})"
            );
            expected_start = window.record_end;
            consumed_blocks += window.block_count;
        }
        assert_eq!(expected_start, total, "windows must cover every record");
    }

    #[test]
    fn window_planner_fills_windows_up_to_the_record_budget() {
        // 20 blocks of 10 records, budget 45 -> 4 whole blocks per window.
        let windows = plan_windows(&[10; 20], 200, 45).unwrap();
        assert_windows_tile(&windows, &[10; 20], 200, 45);
        assert!(windows.iter().all(|window| window.block_count == 4));
        assert_eq!(windows.len(), 5);
    }

    #[test]
    fn window_planner_is_independent_of_decoder_threads() {
        // The plan is a pure function of the block index and the budget: no
        // decoder-thread input, so peak memory cannot scale with --readFilesNthreads.
        // Counts are the measured shape of a 1 MiB-block, quality-bearing,
        // 150 bp single-end CBQ (~4651 records per block).
        let counts = [4_651usize; 43];
        let total: usize = counts.iter().sum();
        let windows = plan_windows(&counts, total, MAX_RECORDS_PER_WINDOW).unwrap();
        assert_windows_tile(&windows, &counts, total, MAX_RECORDS_PER_WINDOW);
        assert_eq!(windows.first().unwrap().block_count, 7);
        assert!(
            windows.len() <= 7,
            "windows must stay large enough to amortize the extra trailing block \
             decompression (one per window): {windows:?}"
        );
    }

    #[test]
    fn window_planner_caps_records_regardless_of_block_size() {
        // A file whose blocks dwarf the budget: each window is exactly one block,
        // so peak memory is one block rather than a multiple of it.
        let counts = [800_000usize; 4];
        let windows = plan_windows(&counts, 3_200_000, MAX_RECORDS_PER_WINDOW).unwrap();
        assert_windows_tile(&windows, &counts, 3_200_000, MAX_RECORDS_PER_WINDOW);
        assert!(windows.iter().all(|window| window.block_count == 1));

        // Mixed sizes: an oversized block never drags neighbours into its window.
        let counts = [10usize, 500, 10, 10];
        let windows = plan_windows(&counts, 530, 100).unwrap();
        assert_windows_tile(&windows, &counts, 530, 100);
        assert_eq!(
            windows,
            vec![
                CbqWindow {
                    record_start: 0,
                    record_end: 10,
                    block_count: 1
                },
                CbqWindow {
                    record_start: 10,
                    record_end: 510,
                    block_count: 1
                },
                CbqWindow {
                    record_start: 510,
                    record_end: 530,
                    block_count: 2
                },
            ]
        );
    }

    #[test]
    fn window_planner_stops_inside_final_block() {
        let windows = plan_windows(&[10, 10, 10], 17, 25).unwrap();
        assert_eq!(
            windows,
            vec![CbqWindow {
                record_start: 0,
                record_end: 17,
                block_count: 2,
            }]
        );
    }

    #[test]
    fn window_planner_rejects_a_short_block_index() {
        // The index advertises more records than its block headers account for.
        let error = plan_windows(&[10, 10], 30, 25).unwrap_err();
        assert!(error.to_string().contains("fewer records than advertised"));
    }

    #[test]
    fn single_end_decoding_is_ordered_and_honors_limit() {
        let path = write_cbq(false, true, true, 75);
        let one_thread = collect_single(&path, 1, 53);
        let four_threads = collect_single(&path, 4, 53);

        assert_eq!(one_thread.len(), 53);
        assert_eq!(four_threads.len(), 53);
        for (i, (one, four)) in one_thread.iter().zip(&four_threads).enumerate() {
            assert_eq!(one.name, format!("read{i}/1"));
            assert_eq!(one.name, four.name);
            assert_eq!(one.sequence, vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4]);
            assert_eq!(one.sequence, four.sequence);
            assert_eq!(one.quality, b"ABCDEFGHIJ");
        }
    }

    /// The record budget makes windows thread-independent, so multi-window
    /// decoding must be exercised explicitly: the alignment-batch buffer is
    /// carried across window boundaries while blocks still arrive out of order
    /// inside each window. Output must be identical for every combination.
    #[test]
    fn decoding_is_identical_across_decoder_threads_and_window_sizes() {
        let path = write_cbq(false, true, true, 300);
        let counts: Vec<usize> = open_reader(&path)
            .unwrap()
            .iter_block_headers()
            .map(|header| header.unwrap().num_records as usize)
            .collect();
        assert!(counts.len() >= 8, "fixture needs many blocks: {counts:?}");

        let reference = collect_single(&path, 1, usize::MAX);
        assert_eq!(reference.len(), 300);
        for window_records in [1usize, 17, 64, 250, MAX_RECORDS_PER_WINDOW] {
            let windows = plan_windows(&counts, 300, window_records).unwrap();
            for threads in [1usize, 2, 3, 4, 8] {
                let decoded = collect_single_windowed(&path, threads, usize::MAX, window_records);
                assert_eq!(
                    decoded.len(),
                    reference.len(),
                    "threads={threads} window_records={window_records}"
                );
                for (expected, actual) in reference.iter().zip(&decoded) {
                    assert_eq!(
                        (&expected.name, &expected.sequence, &expected.quality),
                        (&actual.name, &actual.sequence, &actual.quality),
                        "threads={threads} window_records={window_records}"
                    );
                }
                // A truncating limit must also not leak later records.
                let limited = collect_single_windowed(&path, threads, 137, window_records);
                assert_eq!(limited.len(), 137);
                assert_eq!(limited.last().unwrap().name, reference[136].name);
            }
            if window_records <= 64 {
                assert!(
                    windows.len() >= 3,
                    "expected multi-window decoding for window_records={window_records}"
                );
            }
        }

        // Same for paired records, where each window carries two mates per index.
        let paired_path = write_cbq(true, true, true, 300);
        let paired_reference = collect_paired(&paired_path, 1);
        assert_eq!(paired_reference.len(), 300);
        for window_records in [1usize, 23, 100] {
            for threads in [1usize, 3, 8] {
                let decoded = collect_paired_windowed(&paired_path, threads, window_records);
                let names: Vec<_> = decoded.iter().map(|read| read.name.clone()).collect();
                let expected: Vec<_> = paired_reference
                    .iter()
                    .map(|read| read.name.clone())
                    .collect();
                assert_eq!(
                    names, expected,
                    "threads={threads} window_records={window_records}"
                );
            }
        }
    }

    #[test]
    fn paired_qualityless_decoding_uses_generated_names_without_exporting_placeholders() {
        let path = write_cbq(true, false, false, 41);
        let reads = collect_paired(&path, 4);
        assert_eq!(reads.len(), 41);
        for (i, read) in reads.iter().enumerate() {
            assert_eq!(read.name, i.to_string());
            assert_eq!(read.mate1.name, i.to_string());
            assert_eq!(read.mate2.name, i.to_string());
            // No placeholder buffer is allocated for a quality-less record.
            assert_eq!(read.mate1.quality, Vec::<u8>::new());
            assert_eq!(read.mate2.quality, Vec::<u8>::new());
        }
    }

    #[test]
    fn paired_headers_are_normalized() {
        let path = write_cbq(true, true, true, 17);
        let reads = collect_paired(&path, 3);
        assert_eq!(reads.len(), 17);
        for (i, read) in reads.iter().enumerate() {
            assert_eq!(read.name, format!("read{i}"));
            assert_eq!(read.mate1.name, format!("read{i}/1"));
            assert_eq!(read.mate2.name, format!("read{i}/2"));
        }
    }

    #[test]
    fn read_opts_apply_name_separator_and_quality_shift() {
        let params = crate::params::Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "x.fq",
            "--readNameSeparator",
            "d",
            "--outQSconversionAdd",
            "-2",
        ]);
        let path = write_cbq(false, true, true, 1);
        let (sender, receiver) = sync_channel(8);
        single_end_producer(path.to_path_buf(), 1, CbqReadOpts::from_params(&params))
            .produce(10, usize::MAX, &sender)
            .unwrap();
        drop(sender);
        let read = receiver
            .into_iter()
            .flat_map(Result::unwrap)
            .next()
            .unwrap();
        assert_eq!(read.name, "rea");
        assert_eq!(read.quality, b"?@ABCDEFGH");
    }

    #[test]
    fn invalid_utf8_and_mismatched_pair_names_are_contextual_errors() {
        let invalid = write_custom_headers(&[0xff], None);
        let (sender, _receiver) = sync_channel(8);
        let error = single_end_producer(invalid.to_path_buf(), 1, CbqReadOpts::default())
            .produce(10, usize::MAX, &sender)
            .unwrap_err();
        let message = error.to_string();
        assert!(message.contains("invalid UTF-8"));
        assert!(message.contains("record 0 mate 1"));

        let mismatched = write_custom_headers(b"left/1", Some(b"right/2"));
        let (sender, _receiver) = sync_channel(8);
        let error = paired_end_producer(mismatched.to_path_buf(), 2, CbqReadOpts::default())
            .produce(10, usize::MAX, &sender)
            .unwrap_err();
        assert!(error.to_string().contains("paired read names do not match"));
    }

    #[test]
    fn dropped_alignment_receiver_cancels_without_deadlock() {
        let path = write_cbq(false, true, true, 75);
        let (sender, receiver) = sync_channel(0);
        drop(receiver);
        single_end_producer(path.to_path_buf(), 4, CbqReadOpts::default())
            .produce(1, usize::MAX, &sender)
            .unwrap();
    }

    #[test]
    fn ordered_coordinator_reassembles_blocks_and_rejects_duplicates() {
        fn block(first: usize, values: &[usize]) -> DecodedBlock<usize> {
            DecodedBlock {
                first_index: first,
                end_index: first + values.len(),
                records: values
                    .iter()
                    .enumerate()
                    .map(|(offset, value)| IndexedRecord {
                        index: first + offset,
                        record: *value,
                    })
                    .collect(),
            }
        }

        let mut coordinator = OrderedBlocks::new(PathBuf::from("test.cbq"), 0);
        assert_eq!(coordinator.push(block(2, &[2, 3])).unwrap().len(), 0);
        assert_eq!(
            coordinator.push(block(0, &[0, 1])).unwrap(),
            vec![0, 1, 2, 3]
        );
        assert!(coordinator.push(block(4, &[4])).is_ok());
        assert!(coordinator.push(block(4, &[4])).is_err());
    }

    #[test]
    fn corrupt_compressed_block_returns_error_without_hanging() {
        let path = write_cbq(false, true, true, 30);
        let data_offset = (std::mem::size_of::<binseq::cbq::FileHeader>()
            + std::mem::size_of::<binseq::cbq::BlockHeader>()) as u64;
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        file.seek(SeekFrom::Start(data_offset)).unwrap();
        let mut byte = [0u8; 1];
        file.read_exact(&mut byte).unwrap();
        byte[0] ^= 0xff;
        file.seek(SeekFrom::Start(data_offset)).unwrap();
        file.write_all(&byte).unwrap();
        file.flush().unwrap();

        let (sender, _receiver) = sync_channel(64);
        let result = single_end_producer(path.to_path_buf(), 3, CbqReadOpts::default()).produce(
            10,
            usize::MAX,
            &sender,
        );
        assert!(result.is_err());
    }
}
