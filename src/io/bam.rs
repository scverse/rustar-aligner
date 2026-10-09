/// BAM output writer with noodles (streaming, unsorted)
use crate::error::Error;
use crate::genome::Genome;
use crate::io::bgzf_writer::BgzfWriter;
use crate::io::encode::{RecordEncoder, RecordFormat, bam_record_spans};
use crate::params::Parameters;
use crate::quant::transcriptome::TranscriptomeIndex;
use byteorder::{LittleEndian, WriteBytesExt};
use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::{bam, bgzf, sam};
use std::ffi::CString;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

/// Buffer for BAM records built by parallel threads
#[derive(Default)]
pub struct BufferedBamRecords {
    pub records: Vec<RecordBuf>,
}

impl BufferedBamRecords {
    /// Create new buffer with capacity
    pub fn new() -> Self {
        Self {
            records: Vec::with_capacity(10000),
        }
    }

    /// Add a record to the buffer
    pub fn push(&mut self, record: RecordBuf) {
        self.records.push(record);
    }
}

/// Convert STAR's `--outBAMcompression` integer (-1..9) to a noodles level.
///
/// STAR mapping: -1 or 0 = uncompressed, 1-9 = deflate levels, default 1.
fn bgzf_compression(level: i32) -> bgzf::io::writer::CompressionLevel {
    use bgzf::io::writer::CompressionLevel;
    match level {
        n if n <= 0 => CompressionLevel::NONE,
        n if n >= 9 => CompressionLevel::BEST,
        n => CompressionLevel::try_from(n as u8).unwrap_or_default(),
    }
}

/// Create a BGZF writer with the given STAR compression level.
///
/// `threads > 1` compresses blocks on that many dedicated worker threads
/// (see [`BgzfWriter`]); the output bytes do not depend on `threads`.
fn make_bgzf_writer<W: Write + Send + 'static>(
    inner: W,
    compression: i32,
    threads: usize,
) -> Result<BgzfWriter<W>, Error> {
    Ok(BgzfWriter::new(
        inner,
        bgzf_compression(compression),
        threads,
    )?)
}

/// BGZF compression workers for a run: one per `--runThreadN`, like STAR,
/// which deflates its BAM output inside every alignment thread. The workers
/// only use CPU while there are blocks to compress.
fn bgzf_threads(params: &Parameters) -> usize {
    params.run_thread_n.get()
}

type BgzfFile = BgzfWriter<BufWriter<File>>;
type BgzfStdout = BgzfWriter<BufWriter<std::io::Stdout>>;

/// BAM file writer (streaming, unsorted)
///
/// This writer streams BAM records directly to disk as they're generated,
/// without buffering or sorting. The output is BGZF-compressed but unsorted.
/// Users can sort the output with `samtools sort` if needed.
pub struct BamWriter {
    writer: bam::io::Writer<BgzfFile>,
    header: sam::Header,
}

/// Coordinate-sort RAM budget when `--limitBAMsortRAM` is 0.
///
/// `0` previously meant "unlimited", so the default configuration buffered the
/// entire output (measured: 723 B/record, i.e. ~116 GB for a 160 M-record human
/// sample). A fixed modest budget replaces that.
///
/// Deliberately *not* derived from the genome index size, even though that is
/// STAR's documented rule: scaling the buffer up with the genome enlarges it
/// exactly when RAM is tightest (a 32 GB human index on a 48 GB host would get a
/// multi-GB sort buffer on top). Spilling is measured to be free — 5.83 s with 13
/// spill runs versus 6.15 s unbounded for 1.2 M records — so there is nothing to
/// buy by sorting more in memory. Output is identical either way; only peak memory
/// and temp-file use change.
const DEFAULT_BAM_SORT_RAM: u64 = 512 << 20;

/// Maximum spill runs merged at once, bounding open file descriptors.
///
/// A small `--limitBAMsortRAM` on a large run produces thousands of runs, and
/// opening them all at once hits `EMFILE` (macOS defaults to 256 descriptors;
/// containers are often lower). Above this, runs are merged in balanced passes,
/// costing one extra read/write of the data per `log64(runs)` level.
const MAX_OPEN_RUNS: usize = 64;

fn resolve_bam_sort_ram(params: &Parameters) -> u64 {
    if params.limit_bam_sort_ram > 0 {
        params.limit_bam_sort_ram
    } else {
        DEFAULT_BAM_SORT_RAM
    }
}

/// Directory for coordinate-sort spill files.
///
/// Deliberately the output directory rather than the system temp dir: `/tmp` is
/// tmpfs on many Linux distributions, so spilling there would keep the records in
/// RAM and defeat the budget entirely. STAR likewise keeps its sort scratch
/// (`_STARtmp`) beside the output.
/// Where spill runs go. STAR puts its sort temporaries under `--outTmpDir` when
/// given (`Parameters.cpp:493-498`), so that is honoured. Otherwise they sit
/// beside the output rather than in STAR's `<prefix>_STARtmp/`; see
/// DIVERGENCE.md §4.3.
fn sort_temp_dir(params: &Parameters) -> std::path::PathBuf {
    if params.out_tmp_dir != "-" {
        return std::path::PathBuf::from(&params.out_tmp_dir);
    }
    let prefix = params.output_path("");
    match prefix.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent.to_path_buf(),
        _ => std::path::PathBuf::from("."),
    }
}

/// One buffered record: `(sort key, offset into the byte buffer, length)`.
type SortEntry = ((usize, usize), usize, usize);

/// Coordinate sort key of one BAM record (including its `block_size` prefix):
/// `(reference id, 1-based start)`, with records lacking either sorted last,
/// matching STAR.
fn bam_sort_key(rec: &[u8]) -> (usize, usize) {
    // After block_size: refID (i32), pos (i32, 0-based, -1 = none).
    let ref_id = i32::from_le_bytes(rec[4..8].try_into().unwrap());
    let pos = i32::from_le_bytes(rec[8..12].try_into().unwrap());
    if ref_id >= 0 && pos >= 0 {
        (ref_id as usize, pos as usize + 1)
    } else {
        (usize::MAX, 0)
    }
}

/// External coordinate sorter with bounded memory.
///
/// Records arrive BAM-encoded (serialized on the align workers, #223) and stay
/// encoded: the buffer is one byte vector plus a small index of
/// `(key, offset, length)`, so sorting moves index entries and writing is a
/// plain copy. Once the buffered bytes plus index reach the RAM budget the
/// buffer is sorted and written to a spill run (a headerless BGZF stream of BAM
/// records) beside the output. `write_sorted` merges every run plus the
/// in-memory tail with a k-way merge, so peak memory is the budget rather than
/// the whole output.
///
/// Ordering is identical to a single in-memory sort: runs are created in input
/// order, the sort is stable within a run, and the merge breaks ties on run
/// index, so records sharing a coordinate keep their input order.
///
/// The budget counts the real encoded bytes and the index entries. The realized
/// footprint is somewhat above the budget because the vectors grow by doubling
/// and the spill compressor and merge readers sit outside this accounting.
struct CoordinateSorter {
    encoder: RecordEncoder,
    bytes: Vec<u8>,
    index: Vec<SortEntry>,
    runs: Vec<tempfile::TempPath>,
    header: sam::Header,
    compression: i32,
    threads: usize,
    ram_limit: u64,
    temp_dir: std::path::PathBuf,
    n_records: u64,
}

impl CoordinateSorter {
    fn new(header: sam::Header, params: &Parameters) -> Self {
        Self {
            encoder: RecordEncoder::new(RecordFormat::Bam, &header),
            bytes: Vec::new(),
            index: Vec::new(),
            runs: Vec::new(),
            header,
            compression: params.out_bam_compression,
            threads: bgzf_threads(params),
            ram_limit: resolve_bam_sort_ram(params),
            temp_dir: sort_temp_dir(params),
            n_records: 0,
        }
    }

    /// Bytes the in-memory buffer accounts against the budget.
    fn buffered_bytes(&self) -> u64 {
        (self.bytes.len() + self.index.len() * std::mem::size_of::<SortEntry>()) as u64
    }

    fn push_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        let mut buf = Vec::new();
        self.encoder.encode(batch, &mut buf)?;
        self.push_encoded(&buf)
    }

    /// Buffer records serialized by `self.encoder`, spilling once over budget.
    fn push_encoded(&mut self, bytes: &[u8]) -> Result<(), Error> {
        let base = self.bytes.len();
        for (off, len) in bam_record_spans(bytes)? {
            let key = bam_sort_key(&bytes[off..off + len]);
            self.index.push((key, base + off, len));
            self.n_records += 1;
        }
        self.bytes.extend_from_slice(bytes);
        if self.buffered_bytes() >= self.ram_limit {
            self.spill()?;
        }
        Ok(())
    }

    /// Create an empty spill-run file beside the output.
    fn new_run(&self) -> Result<tempfile::NamedTempFile, Error> {
        // Created on the first spill only, so a run that never spills leaves
        // no directory behind.
        std::fs::create_dir_all(&self.temp_dir)
            .map_err(|source| Error::io(source, &self.temp_dir))?;
        tempfile::Builder::new()
            .prefix("rustar-bamsort-")
            .suffix(".tmp")
            .tempfile_in(&self.temp_dir)
            .map_err(|source| Error::io(source, &self.temp_dir))
    }

    /// Stable-sort the buffered records by coordinate and write them to `out`.
    fn write_buffered<W: Write>(&mut self, out: &mut W) -> Result<(), Error> {
        self.index.sort_by_key(|&(key, _, _)| key);
        for &(_, off, len) in &self.index {
            out.write_all(&self.bytes[off..off + len])?;
        }
        Ok(())
    }

    /// Sort the in-memory buffer and append it as a spill run.
    fn spill(&mut self) -> Result<(), Error> {
        if self.index.is_empty() {
            return Ok(());
        }
        let temp = self.new_run()?;
        // Headerless: runs carry only record blocks and never re-parse
        // reference names.
        let mut writer = make_bgzf_writer(
            BufWriter::new(temp.as_file().try_clone()?),
            self.compression,
            self.threads,
        )?;
        self.write_buffered(&mut writer)?;
        // Finish and flush explicitly rather than on drop, so a failure to write
        // the run surfaces here instead of being swallowed and read back short.
        writer.finish()?.flush()?;
        log::debug!(
            "Coordinate sort: spilled run {} ({} records, ~{} MiB)",
            self.runs.len(),
            self.index.len(),
            self.buffered_bytes() >> 20
        );
        self.runs.push(temp.into_temp_path());
        // Keep the capacity: it is the budget, so holding one stable allocation
        // for the whole run is cheaper than releasing it and growing back.
        self.bytes.clear();
        self.index.clear();
        Ok(())
    }

    /// Merge every spill run and the in-memory tail into `out` as a sorted BAM.
    fn write_sorted<W: Write + Send + 'static>(
        &mut self,
        out: W,
        destination: &str,
    ) -> Result<(), Error> {
        let spilled_runs = self.runs.len();
        self.reduce_runs()?;
        self.index.sort_by_key(|&(key, _, _)| key);

        let mut bgzf = make_bgzf_writer(out, self.compression, self.threads)?;
        write_bam_header_lenient(&mut bgzf, &self.header, Some("coordinate"))?;

        let written = if self.runs.is_empty() {
            // Nothing spilled: identical to the previous in-memory-only path.
            self.write_buffered(&mut bgzf)?;
            self.index.len() as u64
        } else {
            let runs = std::mem::take(&mut self.runs);
            let written = self.merge(&runs, true, &mut bgzf)?;
            drop(runs);
            written
        };
        bgzf.finish()?.flush()?;

        if written != self.n_records {
            return Err(Error::Alignment(format!(
                "BAM sort merge wrote {written} records but {} were buffered",
                self.n_records
            )));
        }
        log::info!(
            "Sorted BAM written to {destination} ({} records, {spilled_runs} spill run(s), \u{2264}{} MiB sort buffer)",
            self.n_records,
            self.ram_limit >> 20
        );
        Ok(())
    }

    /// Merge runs in balanced passes until at most `MAX_OPEN_RUNS` remain, so the
    /// final merge cannot exhaust the process's file descriptors.
    ///
    /// Each pass merges consecutive groups of `MAX_OPEN_RUNS` into one run and
    /// keeps the groups in order, so a lower run index still means "earlier in the
    /// input" and coordinate ties keep resolving to input order.
    fn reduce_runs(&mut self) -> Result<(), Error> {
        while self.runs.len() > MAX_OPEN_RUNS {
            let mut remaining = std::mem::take(&mut self.runs);
            let mut reduced = Vec::with_capacity(remaining.len().div_ceil(MAX_OPEN_RUNS));
            while !remaining.is_empty() {
                let group: Vec<_> = remaining
                    .drain(..MAX_OPEN_RUNS.min(remaining.len()))
                    .collect();
                if group.len() == 1 {
                    reduced.extend(group);
                    continue;
                }
                let temp = self.new_run()?;
                let mut writer = make_bgzf_writer(
                    BufWriter::new(temp.as_file().try_clone()?),
                    self.compression,
                    self.threads,
                )?;
                self.merge(&group, false, &mut writer)?;
                writer.finish()?.flush()?;
                // Dropping `group` here deletes the consumed runs, so scratch use
                // does not grow across passes.
                drop(group);
                reduced.push(temp.into_temp_path());
            }
            log::debug!(
                "Coordinate sort: reduced spill runs to {} (cap {MAX_OPEN_RUNS})",
                reduced.len()
            );
            self.runs = reduced;
        }
        Ok(())
    }

    /// K-way merge `runs` (and, if `with_tail`, the sorted in-memory buffer, which
    /// sorts last on ties) into `out`. Returns the number of records written.
    fn merge<W: Write>(
        &self,
        runs: &[tempfile::TempPath],
        with_tail: bool,
        out: &mut W,
    ) -> Result<u64, Error> {
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;

        let mut readers = Vec::with_capacity(runs.len());
        for path in runs {
            let file = File::open(path).map_err(|source| Error::io(source, path))?;
            readers.push(bgzf::io::Reader::new(BufReader::new(file)));
        }
        // The unspilled tail is the last run, so it keeps the highest run index
        // and therefore loses coordinate ties to everything written before it.
        let tail_run = readers.len();
        let mut tail = if with_tail {
            self.index.as_slice()
        } else {
            &[]
        }
        .iter();

        let mut heads: Vec<Option<Vec<u8>>> = Vec::with_capacity(tail_run + 1);
        let mut heap: BinaryHeap<Reverse<((usize, usize), usize)>> = BinaryHeap::new();
        for (run, (reader, path)) in readers.iter_mut().zip(runs).enumerate() {
            let record = read_run_record(reader, path)?;
            if let Some(record) = &record {
                heap.push(Reverse((bam_sort_key(record), run)));
            }
            heads.push(record);
        }
        let tail_head = tail.next().map(|&(key, off, len)| (key, off, len));
        if let Some((key, _, _)) = tail_head {
            heap.push(Reverse((key, tail_run)));
        }
        heads.push(None);
        let mut tail_head = tail_head;

        let mut written = 0u64;
        while let Some(Reverse((_, run))) = heap.pop() {
            let next = if run == tail_run {
                let (_, off, len) = tail_head
                    .take()
                    .ok_or_else(|| Error::Alignment("BAM sort merge lost a record".to_string()))?;
                out.write_all(&self.bytes[off..off + len])?;
                tail_head = tail.next().copied();
                if let Some((key, _, _)) = tail_head {
                    heap.push(Reverse((key, run)));
                }
                written += 1;
                continue;
            } else {
                let record = heads[run]
                    .take()
                    .ok_or_else(|| Error::Alignment("BAM sort merge lost a record".to_string()))?;
                out.write_all(&record)?;
                read_run_record(&mut readers[run], &runs[run])?
            };
            written += 1;
            if let Some(record) = &next {
                heap.push(Reverse((bam_sort_key(record), run)));
            }
            heads[run] = next;
        }
        Ok(written)
    }
}

/// Read the next BAM record (with its `block_size` prefix) from a spill run, or
/// `None` at end of run.
fn read_run_record<R: std::io::Read>(
    reader: &mut bgzf::io::Reader<R>,
    path: &Path,
) -> Result<Option<Vec<u8>>, Error> {
    use std::io::Read as _;
    let mut size = [0u8; 4];
    let mut got = 0;
    while got < size.len() {
        match reader
            .read(&mut size[got..])
            .map_err(|source| Error::io(source, path))?
        {
            0 if got == 0 => return Ok(None),
            0 => {
                return Err(Error::Alignment(format!(
                    "truncated BAM sort spill run {}",
                    path.display()
                )));
            }
            n => got += n,
        }
    }
    let len = u32::from_le_bytes(size) as usize;
    let mut record = vec![0u8; 4 + len];
    record[..4].copy_from_slice(&size);
    reader
        .read_exact(&mut record[4..])
        .map_err(|source| Error::io(source, path))?;
    Ok(Some(record))
}

/// BAM file writer that sorts records by coordinate with bounded memory,
/// spilling sorted runs beside the output and merging them on `finish()`.
///
/// The header emits `SO:coordinate`. Unmapped records sort to the end.
pub struct SortedBamWriter {
    sorter: CoordinateSorter,
    output_path: std::path::PathBuf,
}

impl BamWriter {
    /// Create a BAM writer at `output_path` with the given prepared header.
    ///
    /// Uses a lenient header writer that bypasses noodles' SAM-spec
    /// reference-name validator. STAR accepts `(` / `)` in reference names
    /// (yeast tRNA transcripts like `tP(UGG)A`), which noodles' strict
    /// validator rejects. Writing the SAM text portion of the BAM header
    /// manually sidesteps that validation while preserving every other byte.
    fn with_header(
        output_path: &Path,
        header: sam::Header,
        compression: i32,
        threads: usize,
    ) -> Result<Self, Error> {
        let buf_writer = BufWriter::new(File::create(output_path)?);
        let mut bgzf = make_bgzf_writer(buf_writer, compression, threads)?;
        write_bam_header_lenient(&mut bgzf, &header, None)?;
        let writer = bam::io::Writer::from(bgzf);
        Ok(Self { writer, header })
    }

    /// Create a new BAM writer with header from genome index.
    pub fn create(output_path: &Path, genome: &Genome, params: &Parameters) -> Result<Self, Error> {
        Self::with_header(
            output_path,
            crate::io::sam::build_sam_header(genome, params)?,
            params.out_bam_compression,
            bgzf_threads(params),
        )
    }

    /// Create a BAM writer whose @SQ header lists the transcripts (not
    /// chromosomes) from `tr_idx`.  Used for
    /// `Aligned.toTranscriptome.out.bam`.
    pub fn create_transcriptome(
        output_path: &Path,
        tr_idx: &TranscriptomeIndex,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let refs = tr_idx
            .tr_ids
            .iter()
            .zip(tr_idx.tr_length.iter())
            .map(|(id, len)| (id.as_str(), *len as usize));
        Self::with_header(
            output_path,
            crate::io::sam::build_sam_header_from_refs(refs, params)?,
            params.out_bam_compression,
            bgzf_threads(params),
        )
    }

    /// Write batch of buffered records (for parallel processing)
    ///
    /// # Arguments
    /// * `batch` - Slice of records to write
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        for record in batch {
            self.writer.write_alignment_record(&self.header, record)?;
        }
        Ok(())
    }

    /// Encoder producing exactly what [`BamWriter::write_batch`] writes.
    pub fn encoder(&self) -> RecordEncoder {
        RecordEncoder::new(RecordFormat::Bam, &self.header)
    }

    /// Append records serialized by [`BamWriter::encoder`].
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<(), Error> {
        self.writer.get_mut().write_all(bytes)?;
        Ok(())
    }

    /// Write `encoded` if the records were pre-serialized, else `records`.
    pub fn write_encoded_or(
        &mut self,
        encoded: Option<&[u8]>,
        records: &[RecordBuf],
    ) -> Result<(), Error> {
        match encoded {
            Some(bytes) => self.write_encoded(bytes),
            None => self.write_batch(records),
        }
    }

    /// Flush and close BAM file
    pub fn finish(&mut self) -> Result<(), Error> {
        self.writer.get_mut().finish()?;
        log::info!("BAM file written successfully");
        Ok(())
    }
}

impl SortedBamWriter {
    /// Create a sorted BAM writer. Records are buffered up to the
    /// `--limitBAMsortRAM` budget and spilled to sorted runs beyond it.
    pub fn create(
        output_path: &Path,
        genome: &crate::genome::Genome,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        Ok(Self {
            sorter: CoordinateSorter::new(header, params),
            output_path: output_path.to_path_buf(),
        })
    }

    /// Buffer records, spilling a sorted run when the RAM budget is reached.
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        self.sorter.push_batch(batch)
    }

    /// Encoder whose output `write_encoded` accepts.
    pub fn encoder(&self) -> RecordEncoder {
        self.sorter.encoder.clone()
    }

    /// Buffer records serialized by `encoder`, spilling once over budget.
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<(), Error> {
        self.sorter.push_encoded(bytes)
    }

    /// Merge every run plus the in-memory tail into a coordinate-sorted BAM.
    pub fn finish(&mut self) -> Result<(), Error> {
        let file = File::create(&self.output_path)
            .map_err(|source| Error::io(source, &self.output_path))?;
        self.sorter.write_sorted(
            BufWriter::new(file),
            &self.output_path.display().to_string(),
        )
    }
}

/// Write a BAM header that tolerates reference sequence names STAR emits
/// (e.g. `tP(UGG)A`) but that noodles' SAM-spec validator rejects.
///
/// Replicates `noodles_bam::io::writer::header::write_header` byte-for-byte
/// for compliant headers; the only divergence is that the SAM text block
/// between `BAM\x01` and the binary reference list is produced via a local
/// formatter instead of `sam::io::Writer::write_header`, so forbidden-char
/// names pass through unchanged.
///
/// Binary reference list (after the text block) uses `CString::new` — the
/// only constraint there is "no interior nul", which is enforced upstream
/// via the usual UTF-8 input.
///
/// `sort_order`: if `Some("coordinate")`, injects `SO:coordinate` into the
/// @HD line. Pass `None` for unsorted output.
fn write_bam_header_lenient<W: Write>(
    writer: &mut W,
    header: &sam::Header,
    sort_order: Option<&str>,
) -> Result<(), Error> {
    const MAGIC: &[u8; 4] = b"BAM\x01";

    writer.write_all(MAGIC)?;

    // Build the SAM text block byte-for-byte identical to
    // `sam::io::Writer::write_header` minus the name validator:
    // `@HD`, `@SQ` (one per reference), `@RG` (if any), `@PG` (if any),
    // `@CO` (if any), each line terminated by `\n`.
    let text = render_sam_text_lenient(header, sort_order);
    let l_text = i32::try_from(text.len()).map_err(|_| {
        Error::Index(format!(
            "BAM SAM-text header exceeds i32::MAX bytes: {} bytes",
            text.len()
        ))
    })?;
    writer.write_i32::<LittleEndian>(l_text)?;
    writer.write_all(&text)?;

    // Binary reference list: n_ref then (l_name, name\0, l_ref) per ref.
    let refs = header.reference_sequences();
    let n_ref = i32::try_from(refs.len())
        .map_err(|_| Error::Index("BAM reference count exceeds i32::MAX".into()))?;
    writer.write_i32::<LittleEndian>(n_ref)?;
    for (name, rs) in refs {
        let c_name = CString::new(name.to_vec())
            .map_err(|e| Error::Index(format!("reference name contains interior NUL byte: {e}")))?;
        let name_bytes = c_name.as_bytes_with_nul();
        let l_name = u32::try_from(name_bytes.len()).map_err(|_| {
            Error::Index(format!(
                "reference name longer than u32::MAX: {} bytes",
                name_bytes.len()
            ))
        })?;
        writer.write_u32::<LittleEndian>(l_name)?;
        writer.write_all(name_bytes)?;
        let l_ref = i32::try_from(usize::from(rs.length()))
            .map_err(|_| Error::Index("reference length exceeds i32::MAX".into()))?;
        writer.write_i32::<LittleEndian>(l_ref)?;
    }

    Ok(())
}

fn render_sam_text_lenient(header: &sam::Header, sort_order: Option<&str>) -> Vec<u8> {
    let mut buf: Vec<u8> = Vec::new();

    // @HD line. noodles' Map<Header> in this version doesn't expose
    // sort_order/group_order via dedicated accessors; we serialize through
    // the generic `other_fields` map (which includes SO/GO when set).
    if let Some(hd) = header.header() {
        buf.extend_from_slice(b"@HD\tVN:");
        buf.extend_from_slice(hd.version().to_string().as_bytes());
        if let Some(so) = sort_order {
            buf.extend_from_slice(b"\tSO:");
            buf.extend_from_slice(so.as_bytes());
        }
        for (tag, value) in hd.other_fields() {
            buf.push(b'\t');
            buf.extend_from_slice(tag.as_ref());
            buf.push(b':');
            buf.extend_from_slice(value);
        }
        buf.push(b'\n');
    }

    // @SQ lines. Use the raw bytes so forbidden characters pass through.
    for (name, rs) in header.reference_sequences() {
        buf.extend_from_slice(b"@SQ\tSN:");
        buf.extend_from_slice(name);
        buf.extend_from_slice(b"\tLN:");
        buf.extend_from_slice(usize::from(rs.length()).to_string().as_bytes());
        // Other optional @SQ fields (AH, AN, AS, DS, M5, SP, TP, UR) —
        // rustar-aligner doesn't set any today, so skip.
        buf.push(b'\n');
    }

    // @RG lines.
    for (id, rg) in header.read_groups() {
        buf.extend_from_slice(b"@RG\tID:");
        buf.extend_from_slice(id);
        for (tag, value) in rg.other_fields() {
            buf.push(b'\t');
            buf.extend_from_slice(tag.as_ref());
            buf.push(b':');
            buf.extend_from_slice(value);
        }
        buf.push(b'\n');
    }

    // @PG lines — noodles' map doesn't guarantee insertion order; for
    // rustar-aligner we emit a single @PG with id "rustar-aligner". If more are added
    // later, pipe them in here.
    for (id, pg) in header.programs().as_ref() {
        buf.extend_from_slice(b"@PG\tID:");
        buf.extend_from_slice(id);
        for (tag, value) in pg.other_fields() {
            buf.push(b'\t');
            buf.extend_from_slice(tag.as_ref());
            buf.push(b':');
            buf.extend_from_slice(value);
        }
        buf.push(b'\n');
    }

    // @CO lines (comments).
    for comment in header.comments() {
        buf.extend_from_slice(b"@CO\t");
        buf.extend_from_slice(comment);
        buf.push(b'\n');
    }

    buf
}

/// Streaming unsorted BAM writer that writes to stdout.
pub struct BamStdoutWriter {
    writer: bam::io::Writer<BgzfStdout>,
    header: sam::Header,
}

impl BamStdoutWriter {
    pub fn create(genome: &crate::genome::Genome, params: &Parameters) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        let mut bgzf = make_bgzf_writer(
            BufWriter::new(std::io::stdout()),
            params.out_bam_compression,
            bgzf_threads(params),
        )?;
        write_bam_header_lenient(&mut bgzf, &header, None)?;
        let writer = bam::io::Writer::from(bgzf);
        Ok(Self { writer, header })
    }

    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        for record in batch {
            self.writer.write_alignment_record(&self.header, record)?;
        }
        Ok(())
    }

    /// Encoder producing exactly what [`BamStdoutWriter::write_batch`] writes.
    pub fn encoder(&self) -> RecordEncoder {
        RecordEncoder::new(RecordFormat::Bam, &self.header)
    }

    /// Append records serialized by [`BamStdoutWriter::encoder`].
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<(), Error> {
        self.writer.get_mut().write_all(bytes)?;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<(), Error> {
        self.writer.get_mut().finish()?;
        Ok(())
    }
}

/// Coordinate-sorted BAM writer that writes to stdout on `finish()`, with the
/// same bounded-memory spill/merge behavior as [`SortedBamWriter`].
pub struct SortedBamStdoutWriter {
    sorter: CoordinateSorter,
}

impl SortedBamStdoutWriter {
    pub fn create(genome: &crate::genome::Genome, params: &Parameters) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        Ok(Self {
            sorter: CoordinateSorter::new(header, params),
        })
    }

    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        self.sorter.push_batch(batch)
    }

    /// Encoder whose output `write_encoded` accepts.
    pub fn encoder(&self) -> RecordEncoder {
        self.sorter.encoder.clone()
    }

    /// Buffer records serialized by `encoder`, spilling once over budget.
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<(), Error> {
        self.sorter.push_encoded(bytes)
    }

    pub fn finish(&mut self) -> Result<(), Error> {
        self.sorter
            .write_sorted(BufWriter::new(std::io::stdout()), "stdout")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{Exon, Transcript};
    use noodles::sam::alignment::record::cigar;
    use tempfile::NamedTempFile;

    fn create_test_genome() -> Genome {
        Genome {
            transform_blocks: None,
            sequence: vec![0, 1, 2, 3, 0, 1, 2, 3].into(), // ACGTACGT
            n_genome: 8,
            n_genome_real: 8,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![8],
            chr_start: vec![0, 8],
        }
    }

    fn default_params() -> Parameters {
        Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"])
    }

    #[test]
    fn test_bam_writer_creation() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let writer = BamWriter::create(temp_file.path(), &genome, &params);
        assert!(writer.is_ok(), "BAM writer creation should succeed");
    }

    #[test]
    fn test_bam_unmapped_write() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        let read_name = "read1";
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        // Build record using SAM builder
        let record = crate::io::sam::SamWriter::build_unmapped_record(
            read_name,
            &read_seq,
            &read_qual,
            &params,
            crate::stats::UnmappedReason::Other,
        )
        .unwrap();

        let result = writer.write_batch(&[record]);
        assert!(result.is_ok(), "Writing unmapped read should succeed");

        // Finish the file
        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_missing_quality_is_encoded_as_absent() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let record = crate::io::sam::SamWriter::build_unmapped_record(
            "read1",
            &[0, 1, 2, 3],
            &[],
            &params,
            crate::stats::UnmappedReason::Other,
        )
        .unwrap();
        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();
        writer.write_batch(&[record]).unwrap();
        writer.finish().unwrap();
        drop(writer);

        let mut reader = bam::io::Reader::new(File::open(temp_file.path()).unwrap());
        reader.read_header().unwrap();
        let record = reader.records().next().unwrap().unwrap();
        assert_eq!(record.sequence().len(), 4);
        assert!(record.quality_scores().is_empty());
    }

    #[test]
    fn test_bam_alignment_write() {
        use cigar::op::{Kind, Op};

        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        // Create a simple transcript
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 104,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 100,
                genome_end: 104,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_name = "read1";
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        // Build records using SAM builder
        let records = crate::io::sam::SamWriter::build_alignment_records(
            read_name,
            &read_seq,
            &read_qual,
            0,
            0,
            &[transcript],
            &genome,
            &params,
            1,
        )
        .unwrap();

        let result = writer.write_batch(&records);
        assert!(result.is_ok(), "Writing alignment should succeed");

        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_transcriptome_writer_creation() {
        use crate::junction::gtf::GtfRecord;
        use crate::quant::transcriptome::TranscriptomeIndex;
        use std::collections::HashMap;

        let genome = create_test_genome();
        // Stretch genome / exon to fit the tiny test chromosome.
        let mut attrs = HashMap::new();
        attrs.insert("gene_id".to_string(), "G1".to_string());
        attrs.insert("transcript_id".to_string(), "T1".to_string());
        let exons = vec![GtfRecord {
            seqname: "chr1".to_string(),
            feature: "exon".to_string(),
            start: 1,
            end: 8,
            strand: '+',
            attributes: attrs,
        }];
        let tr_idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(tr_idx.n_transcripts(), 1);

        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();
        let writer = BamWriter::create_transcriptome(temp_file.path(), &tr_idx, &params);
        assert!(
            writer.is_ok(),
            "transcriptome BAM writer creation should succeed"
        );

        // Header should contain exactly 1 @SQ entry (matching n_transcripts).
        let writer = writer.unwrap();
        assert_eq!(writer.header.reference_sequences().len(), 1);
    }

    #[test]
    fn test_bam_batch_write() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        // Create a batch of unmapped records
        let records = vec![
            crate::io::sam::SamWriter::build_unmapped_record(
                "read1",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                &params,
                crate::stats::UnmappedReason::Other,
            )
            .unwrap(),
            crate::io::sam::SamWriter::build_unmapped_record(
                "read2",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                &params,
                crate::stats::UnmappedReason::Other,
            )
            .unwrap(),
        ];

        let result = writer.write_batch(&records);
        assert!(result.is_ok(), "Writing batch should succeed");

        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_compression_level_zero() {
        let genome = create_test_genome();
        let mut params = default_params();
        params.out_bam_compression = 0;
        let temp_file = NamedTempFile::new().unwrap();
        let writer = BamWriter::create(temp_file.path(), &genome, &params);
        assert!(
            writer.is_ok(),
            "BAM writer with compression=0 should succeed"
        );
    }

    /// A record at `(chr, pos)` named `name`; `None` position means unmapped.
    fn placed_record(name: &str, chr: Option<usize>, pos: Option<usize>) -> RecordBuf {
        let mut record = RecordBuf::default();
        *record.name_mut() = Some(name.into());
        *record.reference_sequence_id_mut() = chr;
        *record.alignment_start_mut() = pos.map(|pos| pos.try_into().unwrap());
        record
    }

    /// `(name, reference_sequence_id, alignment_start)` of a decoded record.
    type DecodedRecord = (String, Option<usize>, Option<usize>);

    /// Decode a sorted BAM into `(name, reference_sequence_id, alignment_start)`.
    fn read_sorted(path: &Path) -> Vec<DecodedRecord> {
        let mut reader = bam::io::Reader::new(File::open(path).unwrap());
        reader.read_header().unwrap();
        let mut out = Vec::new();
        let mut record = RecordBuf::default();
        let header = sam::Header::default();
        while reader.read_record_buf(&header, &mut record).unwrap() != 0 {
            out.push((
                String::from_utf8(record.name().unwrap().to_vec()).unwrap(),
                record.reference_sequence_id(),
                record.alignment_start().map(|p| p.get()),
            ));
        }
        out
    }

    /// A two-chromosome genome, so cross-reference ordering is exercised too.
    ///
    /// References are long enough to contain every position the fixtures use —
    /// records past `LN` would make an out-of-spec BAM that only round-trips
    /// because noodles does not validate position against reference length.
    fn two_chr_genome() -> Genome {
        Genome {
            transform_blocks: None,
            sequence: vec![0u8; 400].into(),
            n_genome: 400,
            n_genome_real: 400,
            n_chr_real: 2,
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
            chr_length: vec![200, 200],
            chr_start: vec![0, 200, 400],
        }
    }

    /// Records in deliberately unsorted input order, spanning two references,
    /// including coordinate ties (so tie-order is observable) and unmapped
    /// records (which must sort last).
    fn shuffled_records() -> Vec<RecordBuf> {
        let mut records = Vec::new();
        for i in 0..200usize {
            // Positions cycle so input order and sorted order disagree, and every
            // (chr, pos) is hit twice to create ties.
            let pos = (i * 37) % 100 + 1;
            records.push(placed_record(&format!("r{i}"), Some(i % 2), Some(pos)));
        }
        for i in 0..10usize {
            records.push(placed_record(&format!("u{i}"), None, None));
        }
        records
    }

    /// Write `records` in batches of `batch` through a `SortedBamWriter` whose
    /// sort budget is `ram_limit`, returning the decoded output.
    fn sorted_output(
        records: &[RecordBuf],
        ram_limit: u64,
        batch: usize,
    ) -> (Vec<DecodedRecord>, usize) {
        let genome = two_chr_genome();
        let dir = tempfile::tempdir().unwrap();
        let mut params = default_params();
        params.limit_bam_sort_ram = ram_limit;
        params.out_file_name_prefix = format!("{}/", dir.path().display());

        let out = dir.path().join("Aligned.sortedByCoord.out.bam");
        let mut writer = SortedBamWriter::create(&out, &genome, &params).unwrap();
        for chunk in records.chunks(batch) {
            writer.write_batch(chunk).unwrap();
        }
        // Captured before `finish()`, which consumes and deletes the runs.
        let spill_runs = writer.sorter.runs.len();
        writer.finish().unwrap();

        let decoded = read_sorted(&out);
        // Spill files must not outlive the merge.
        let leftover = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(Result::ok)
            .filter(|entry| {
                entry
                    .file_name()
                    .to_string_lossy()
                    .starts_with("rustar-bamsort-")
            })
            .count();
        assert_eq!(leftover, 0, "spill files must be cleaned up");
        (decoded, spill_runs)
    }

    #[test]
    fn test_sorted_bam_spill_merge_matches_in_memory_sort() {
        let records = shuffled_records();

        // A budget far above the data never spills — the previous behavior.
        let (in_memory, in_memory_runs) = sorted_output(&records, 1 << 30, 64);
        assert_eq!(in_memory_runs, 0);
        assert_eq!(in_memory.len(), records.len());

        // Tiny budgets force a spill per batch; output must be unchanged.
        for batch in [1usize, 7, 64] {
            let (spilled, runs) = sorted_output(&records, 1, batch);
            assert!(
                runs >= records.len() / batch,
                "expected a spill run per batch (batch={batch}, runs={runs})"
            );
            assert_eq!(
                spilled, in_memory,
                "spill+merge output must equal the in-memory sort (batch={batch})"
            );
        }
    }

    #[test]
    fn test_sorted_bam_coordinate_ties_keep_input_order_across_runs() {
        // Every record shares one coordinate, so the only thing under test is
        // whether the merge preserves input order across spill runs.
        let records: Vec<_> = (0..64)
            .map(|i| placed_record(&format!("r{i:03}"), Some(0), Some(1)))
            .collect();
        let expected: Vec<_> = records
            .iter()
            .map(|r| String::from_utf8(r.name().unwrap().to_vec()).unwrap())
            .collect();

        for batch in [1usize, 3, 16] {
            let (decoded, runs) = sorted_output(&records, 1, batch);
            assert!(runs > 1, "batch={batch} should produce several runs");
            let names: Vec<_> = decoded.into_iter().map(|(name, ..)| name).collect();
            assert_eq!(
                names, expected,
                "tie order must be input order (batch={batch})"
            );
        }
    }

    /// More runs than `MAX_OPEN_RUNS` must merge in passes rather than opening
    /// every run at once — opening them all hits `EMFILE` where the descriptor
    /// limit is low (macOS defaults to 256), and a small `--limitBAMsortRAM` on a
    /// large run produces thousands of runs.
    #[test]
    fn test_sorted_bam_merges_more_runs_than_the_descriptor_cap() {
        // One record per batch with a 1-byte budget => one spill run each.
        let records: Vec<_> = (0..(MAX_OPEN_RUNS * 5 + 3))
            .map(|i| placed_record(&format!("r{i:05}"), Some(i % 2), Some(i % 150 + 1)))
            .collect();
        let expected = {
            let (decoded, runs) = sorted_output(&records, 1 << 30, records.len());
            assert_eq!(runs, 0);
            decoded
        };

        let (decoded, runs) = sorted_output(&records, 1, 1);
        assert!(
            runs > MAX_OPEN_RUNS,
            "fixture must exceed the cap to exercise the reduction pass: {runs}"
        );
        assert_eq!(
            decoded, expected,
            "multi-pass merge must match the in-memory sort exactly"
        );
    }

    #[test]
    fn test_sorted_bam_empty_output_is_a_valid_bam() {
        let (decoded, runs) = sorted_output(&[], 1 << 30, 64);
        assert_eq!(runs, 0);
        assert_eq!(decoded, []);
    }

    #[test]
    fn test_bam_sort_ram_default_is_bounded() {
        // `--limitBAMsortRAM 0` must resolve to a fixed budget, never "unlimited",
        // and must not scale with the genome (which would enlarge the buffer
        // exactly when RAM is tightest).
        let mut params = default_params();
        params.limit_bam_sort_ram = 0;
        assert_eq!(resolve_bam_sort_ram(&params), DEFAULT_BAM_SORT_RAM);
        params.genome_dir = std::path::PathBuf::from("/nonexistent/huge/index");
        assert_eq!(resolve_bam_sort_ram(&params), DEFAULT_BAM_SORT_RAM);

        params.limit_bam_sort_ram = 12345;
        assert_eq!(resolve_bam_sort_ram(&params), 12345);
    }

    #[test]
    fn test_sort_budget_accounts_real_record_size() {
        let params = default_params();
        let mut small_sorter = CoordinateSorter::new(two_ref_header(), &params);
        let mut large_sorter = CoordinateSorter::new(two_ref_header(), &params);
        let small = placed_record("r", Some(0), Some(1));
        let mut large = placed_record("r", Some(0), Some(1));
        *large.sequence_mut() =
            noodles::sam::alignment::record_buf::Sequence::from(vec![b'A'; 150]);
        *large.quality_scores_mut() =
            noodles::sam::alignment::record_buf::QualityScores::from(vec![30u8; 150]);
        small_sorter.push_batch(&[small]).unwrap();
        large_sorter.push_batch(&[large]).unwrap();
        assert!(
            large_sorter.buffered_bytes() >= small_sorter.buffered_bytes() + 150,
            "budget must account for SEQ and QUAL"
        );
    }

    /// Records covering every sort-key case: two references, ties, a record
    /// with a reference but no position, and fully unmapped records.
    fn sort_fixture() -> Vec<RecordBuf> {
        let mk = |name: &str, r: Option<usize>, p: Option<usize>| {
            let mut b = RecordBuf::builder().set_name(name);
            if let Some(r) = r {
                b = b.set_reference_sequence_id(r);
            }
            if let Some(p) = p {
                b = b.set_alignment_start(p.try_into().unwrap());
            }
            b.build()
        };
        vec![
            mk("u1", None, None),
            mk("a", Some(1), Some(5)),
            mk("b", Some(0), Some(7)),
            mk("rnopos", Some(0), None),
            mk("c", Some(0), Some(7)),
            mk("d", Some(0), Some(1)),
            mk("u2", None, None),
            mk("e", Some(1), Some(1)),
        ]
    }

    fn two_ref_header() -> sam::Header {
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use std::num::NonZeroUsize;
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(100).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(100).unwrap()),
            )
            .build()
    }

    #[test]
    fn sorter_matches_recordbuf_sort() {
        let header = two_ref_header();
        let records = sort_fixture();

        // Reference: the RecordBuf-level stable sort the writers used before.
        let mut expected_records = records.clone();
        expected_records.sort_by_key(|r| match (r.reference_sequence_id(), r.alignment_start()) {
            (Some(chr), Some(pos)) => (chr, pos.get()),
            _ => (usize::MAX, 0),
        });
        let mut expected = bam::io::Writer::from(Vec::new());
        for r in &expected_records {
            expected.write_alignment_record(&header, r).unwrap();
        }

        // Fed half as records, half pre-encoded, as the pipelines may.
        let mut sorter = CoordinateSorter::new(header, &default_params());
        sorter.push_batch(&records[..3]).unwrap();
        let mut enc = Vec::new();
        sorter.encoder.encode(&records[3..], &mut enc).unwrap();
        sorter.push_encoded(&enc).unwrap();
        assert_eq!(sorter.index.len(), records.len());
        let mut got = Vec::new();
        sorter.write_buffered(&mut got).unwrap();
        assert_eq!(got, expected.into_inner());
    }

    #[test]
    fn encoder_matches_writer_bytes() {
        let header = two_ref_header();
        let records = sort_fixture();
        let mut direct = bam::io::Writer::from(Vec::new());
        for r in &records {
            direct.write_alignment_record(&header, r).unwrap();
        }
        let mut enc = Vec::new();
        RecordEncoder::new(RecordFormat::Bam, &header)
            .encode(&records, &mut enc)
            .unwrap();
        assert_eq!(enc, direct.into_inner());

        let mut direct = sam::io::Writer::new(Vec::new());
        for r in &records {
            direct.write_alignment_record(&header, r).unwrap();
        }
        let mut enc = Vec::new();
        RecordEncoder::new(RecordFormat::Sam, &header)
            .encode(&records, &mut enc)
            .unwrap();
        assert_eq!(enc, direct.into_inner());
    }
}
