//! Output-record serialization that can run on the align workers (#223).
//!
//! The align pipelines used to hand `RecordBuf`s to the single writer thread,
//! which then encoded every record (SAM text or BAM binary) itself. A
//! [`RecordEncoder`] captures everything that encoding needs (format, header,
//! `--outSAMmode NoQS`), so the rayon workers can serialize each read's records
//! as soon as they are built and the writer thread only appends bytes.
//!
//! The bytes are produced by the same noodles writers the output writers use,
//! record by record, so the concatenation is byte-identical to what the writer
//! thread would have produced.

use std::sync::Arc;

use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf};
use noodles::{bam, sam};

use crate::error::Error;

/// Serialization format of an output stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordFormat {
    /// SAM text lines.
    Sam,
    /// BAM records (`block_size` prefix included), before BGZF compression.
    Bam,
}

/// Serializes alignment records exactly as an output writer would.
#[derive(Clone)]
pub struct RecordEncoder {
    format: RecordFormat,
    header: Arc<sam::Header>,
    strip_quality: bool,
    check_cigar: bool,
}

impl RecordEncoder {
    pub fn new(format: RecordFormat, header: &sam::Header) -> Self {
        Self {
            format,
            header: Arc::new(header.clone()),
            strip_quality: false,
            check_cigar: false,
        }
    }

    pub fn format(&self) -> RecordFormat {
        self.format
    }

    /// Drop quality strings before encoding (`--outSAMmode NoQS`).
    #[must_use]
    pub fn without_quality(mut self) -> Self {
        self.strip_quality = true;
        self
    }

    /// Run the SAM file writer's CIGAR/SEQ length sanity check on every record.
    #[must_use]
    pub fn with_cigar_check(mut self) -> Self {
        self.check_cigar = true;
        self
    }

    /// Append the serialized form of `records` to `out`.
    pub fn encode(&self, records: &[RecordBuf], out: &mut Vec<u8>) -> Result<(), Error> {
        if records.is_empty() {
            return Ok(());
        }
        match self.format {
            RecordFormat::Sam => self.encode_with(&mut sam::io::Writer::new(out), records),
            RecordFormat::Bam => self.encode_with(&mut bam::io::Writer::from(out), records),
        }
    }

    fn encode_with(
        &self,
        writer: &mut impl sam::alignment::io::Write,
        records: &[RecordBuf],
    ) -> Result<(), Error> {
        for record in records {
            if self.check_cigar {
                crate::io::sam::check_cigar_seq_len(record);
            }
            if self.strip_quality {
                let mut stripped = record.clone();
                *stripped.quality_scores_mut() = QualityScores::default();
                writer.write_alignment_record(&self.header, &stripped)?;
            } else {
                writer.write_alignment_record(&self.header, record)?;
            }
        }
        Ok(())
    }
}

/// Split a buffer of BAM-encoded records into `(record_start, record_len)`
/// spans, each including its 4-byte `block_size` prefix.
pub fn bam_record_spans(bytes: &[u8]) -> Result<Vec<(usize, usize)>, Error> {
    let mut spans = Vec::new();
    let mut off = 0;
    while off < bytes.len() {
        let Some(prefix) = bytes.get(off..off + 4) else {
            return Err(Error::Alignment("truncated encoded BAM record".into()));
        };
        let len = 4 + u32::from_le_bytes(prefix.try_into().unwrap()) as usize;
        if off + len > bytes.len() {
            return Err(Error::Alignment("truncated encoded BAM record".into()));
        }
        spans.push((off, len));
        off += len;
    }
    Ok(spans)
}
