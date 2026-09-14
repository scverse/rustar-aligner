//! SAM/BAM output writer with noodles
use crate::align::read_align::PairedAlignment;
use crate::align::transcript::{Transcript, cigar_to_string};
use crate::error::Error;
use crate::genome::Genome;
use crate::io::fastq::{complement_base, decode_base};
use crate::junction::encode_motif;
use crate::mapq::calculate_mapq;
use crate::params::{Parameters, SamAttributes};
use crate::stats::UnmappedReason;
use bstr::BString;
use noodles::sam;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::cigar;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::{
    Map,
    map::{Program, ReadGroup, program::tag as program_tag, tag::Other as HeaderOtherTag},
};
use std::collections::HashSet;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::BufWriter;
use std::num::NonZeroUsize;
use std::path::Path;

/// Buffer for SAM records built by parallel threads
#[derive(Default)]
pub struct BufferedSamRecords {
    pub records: Vec<RecordBuf>,
}

impl BufferedSamRecords {
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

/// Insert STAR-compatible optional tags on unmapped records.
///
/// NH, HI, AS, nM are gated on `outSAMattributes`; `uT:A:` is always emitted
/// (STAR emits it unconditionally so QC tools can parse unmapping categories).
fn insert_unmapped_tags(record: &mut RecordBuf, attrs: SamAttributes, reason: UnmappedReason) {
    let data = record.data_mut();
    if attrs.contains(SamAttributes::NH) {
        data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(0i32));
    }
    if attrs.contains(SamAttributes::HI) {
        data.insert(Tag::HIT_INDEX, Value::from(0i32));
    }
    if attrs.contains(SamAttributes::AS) {
        data.insert(Tag::ALIGNMENT_SCORE, Value::from(0i32));
    }
    if attrs.contains(SamAttributes::NMM) {
        data.insert(Tag::new(b'n', b'M'), Value::from(0i32));
    }
    let ut = match reason {
        UnmappedReason::Other => b'0',
        UnmappedReason::TooShort => b'1',
        UnmappedReason::TooManyMismatches => b'2',
        UnmappedReason::TooManyLoci => b'3',
    };
    data.insert(Tag::new(b'u', b'T'), Value::Character(ut));
}

/// SAM file writer
pub struct SamWriter {
    writer: sam::io::Writer<BufWriter<File>>,
    header: sam::Header,
}

impl SamWriter {
    /// Create a new SAM writer with header from genome index
    ///
    /// # Arguments
    /// * `output_path` - Path to output SAM file
    /// * `genome` - Genome index with chromosome information
    /// * `params` - Parameters (for @PG header)
    pub fn create(output_path: &Path, genome: &Genome, params: &Parameters) -> Result<Self, Error> {
        let file = File::create(output_path)?;
        let buf_writer = BufWriter::new(file);

        let header = build_sam_header(genome, params)?;
        let mut writer = sam::io::Writer::new(buf_writer);

        writer.write_header(&header)?;

        Ok(Self { writer, header })
    }

    /// Write alignment record(s) for a read
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    /// * `transcripts` - Alignment transcripts (1 or more for multi-mappers)
    /// * `genome` - Genome index
    /// * `params` - Parameters
    #[allow(clippy::too_many_arguments)]
    pub fn write_alignment(
        &mut self,
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
        transcripts: &[Transcript],
        genome: &Genome,
        params: &Parameters,
        n_for_mapq: usize,
    ) -> Result<(), Error> {
        if transcripts.is_empty() {
            return Ok(());
        }

        let n_alignments = transcripts.len();
        let max_output = if params.out_sam_mult_nmax < 0 {
            n_alignments
        } else {
            (params.out_sam_mult_nmax as usize).min(n_alignments)
        };
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let attrs = params.out_sam_attributes;
        let rg_id_owned = params.primary_rg_id()?;
        let rg_id = rg_id_owned.as_deref();

        let best_score = transcripts
            .iter()
            .map(|t| t.score)
            .max()
            .unwrap_or(i32::MIN);
        for (hit_index, transcript) in transcripts.iter().take(max_output).enumerate() {
            let mut record = transcript_to_record(
                transcript,
                read_name,
                read_seq,
                read_qual,
                genome,
                mapq,
                max_output,    // NH = number of reported alignments
                hit_index + 1, // 1-based
                params.out_sam_attr_ih_start,
                attrs,
            )?;
            maybe_insert_rg_tag(&mut record, rg_id);
            apply_sam_flag_or_and(&mut record, params);
            apply_primary_flag(&mut record, transcript.score, best_score, params);

            self.writer.write_alignment_record(&self.header, &record)?;
        }

        Ok(())
    }

    /// Write batch of buffered records (for parallel processing)
    ///
    /// # Arguments
    /// * `batch` - Slice of records to write
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        for record in batch {
            // Debug: validate CIGAR vs SEQ length before writing
            let cigar_ops = record.cigar().as_ref();
            let cigar_query_len: usize = cigar_ops
                .iter()
                .filter(|op| op.kind().consumes_read())
                .map(|op| op.len())
                .sum();
            let seq_len = record.sequence().len();
            if cigar_query_len != seq_len && !cigar_ops.is_empty() {
                let name = record
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
                    .unwrap_or_default();
                panic!(
                    "[SAM-MISMATCH] read={} cigar_query_len={} seq_len={} flags={:?} cigar={}",
                    name,
                    cigar_query_len,
                    seq_len,
                    record.flags(),
                    cigar_to_string(cigar_ops)
                );
            }
            self.writer.write_alignment_record(&self.header, record)?;
        }
        Ok(())
    }

    /// Build unmapped record (without writing)
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    /// * `params` - Parameters (used for attribute gating and RG tag)
    /// * `unmapped_reason` - Why the read was not mapped (drives `uT:A:` tag)
    pub fn build_unmapped_record(
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
        params: &Parameters,
        unmapped_reason: UnmappedReason,
    ) -> Result<RecordBuf, Error> {
        let mut record = RecordBuf::default();

        // Name
        record.name_mut().replace(read_name.into());

        // FLAGS: 0x4 (unmapped)
        let flags = sam::alignment::record::Flags::UNMAPPED;
        *record.flags_mut() = flags;

        // Sequence (decode from genome encoding)
        let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);

        *record.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(read_qual));

        let rg_id_owned = params.primary_rg_id()?;
        maybe_insert_rg_tag(&mut record, rg_id_owned.as_deref());
        insert_unmapped_tags(&mut record, params.out_sam_attributes, unmapped_reason);

        Ok(record)
    }

    /// Build alignment records (without writing) for a read
    ///
    /// # Arguments
    /// * `read_name` - Read identifier
    /// * `read_seq` - Read sequence (encoded)
    /// * `read_qual` - Quality scores
    /// * `transcripts` - Alignment transcripts (1 or more for multi-mappers)
    /// * `genome` - Genome index
    /// * `params` - Parameters
    #[allow(clippy::too_many_arguments)]
    pub fn build_alignment_records(
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
        clip5p: usize,
        clip3p: usize,
        transcripts: &[Transcript],
        genome: &Genome,
        params: &Parameters,
        n_for_mapq: usize,
    ) -> Result<Vec<RecordBuf>, Error> {
        if transcripts.is_empty() {
            return Ok(Vec::new());
        }

        // The transcripts were aligned against the clipped read; the core record
        // (CIGAR core, MD, NM) is built from that aligned slice, then apply_read_clips
        // restores the full read + soft-clips. clip5p/clip3p == 0 => aligned == full.
        let aligned_seq = &read_seq[clip5p..read_seq.len() - clip3p];
        let aligned_qual = &read_qual[clip5p..read_qual.len() - clip3p];

        let n_alignments = transcripts.len();
        let max_output = if params.out_sam_mult_nmax < 0 {
            n_alignments
        } else {
            (params.out_sam_mult_nmax as usize).min(n_alignments)
        };
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let attrs = params.out_sam_attributes;
        let rg_id_owned = params.primary_rg_id()?;
        let rg_id = rg_id_owned.as_deref();

        let mut records = Vec::with_capacity(max_output);
        let best_score = transcripts
            .iter()
            .map(|t| t.score)
            .max()
            .unwrap_or(i32::MIN);
        for (hit_index, transcript) in transcripts.iter().take(max_output).enumerate() {
            let mut record = transcript_to_record(
                transcript,
                read_name,
                aligned_seq,
                aligned_qual,
                genome,
                mapq,
                max_output,    // NH = number of reported alignments
                hit_index + 1, // 1-based
                params.out_sam_attr_ih_start,
                attrs,
            )?;
            if clip5p > 0 || clip3p > 0 {
                apply_read_clips(
                    &mut record,
                    &transcript.cigar,
                    read_seq,
                    read_qual,
                    clip5p,
                    clip3p,
                    transcript.is_reverse,
                );
            }
            maybe_insert_rg_tag(&mut record, rg_id);
            apply_sam_flag_or_and(&mut record, params);
            apply_primary_flag(&mut record, transcript.score, best_score, params);
            records.push(record);
        }

        Ok(records)
    }

    /// Build paired-end SAM records (without writing)
    ///
    /// Returns 2 records per pair (one for each mate)
    ///
    /// # Arguments
    /// * `read_name` - Read identifier (base name without /1 or /2)
    /// * `mate1_seq` - First mate sequence (encoded)
    /// * `mate1_qual` - First mate quality scores
    /// * `mate2_seq` - Second mate sequence (encoded)
    /// * `mate2_qual` - Second mate quality scores
    /// * `paired_alignments` - Paired alignments
    /// * `genome` - Genome index
    /// * `params` - Parameters
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_arguments)]
    pub fn build_paired_records(
        read_name: &str,
        mate1_seq: &[u8],
        mate1_qual: &[u8],
        mate2_seq: &[u8],
        mate2_qual: &[u8],
        m1_clip5p: usize,
        m1_clip3p: usize,
        m2_clip5p: usize,
        m2_clip3p: usize,
        paired_alignments: &[PairedAlignment],
        genome: &Genome,
        params: &Parameters,
        n_for_mapq: usize,
    ) -> Result<Vec<RecordBuf>, Error> {
        if paired_alignments.is_empty() {
            // Both mates unmapped
            return Self::build_paired_unmapped_records(
                read_name,
                mate1_seq,
                mate1_qual,
                mate2_seq,
                mate2_qual,
                params,
                UnmappedReason::Other,
            );
        }

        let n_alignments = paired_alignments.len();
        let max_output = if params.out_sam_mult_nmax < 0 {
            n_alignments
        } else {
            (params.out_sam_mult_nmax as usize).min(n_alignments)
        };
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let attrs = params.out_sam_attributes;
        let rg_id_owned = params.primary_rg_id()?;
        let rg_id = rg_id_owned.as_deref();

        let mut records = Vec::with_capacity(max_output * 2);
        let best_score = paired_alignments
            .iter()
            .map(|p| p.combined_wt_score)
            .max()
            .unwrap_or(i32::MIN);

        for (pair_idx, paired_aln) in paired_alignments.iter().take(max_output).enumerate() {
            let hit_index = pair_idx + 1; // 1-based
            // STAR reports the pre-split combined WT score (with length penalty) as AS.
            // This is stored as combined_wt_score, matching STAR's primaryScore.
            let combined_score = paired_aln.combined_wt_score;

            // Create record for mate1 (this=mate1, mate=mate2)
            let mut rec1 = build_paired_mate_record(
                read_name,
                mate1_seq,
                mate1_qual,
                m1_clip5p,
                m1_clip3p,
                &paired_aln.mate1_transcript,
                &paired_aln.mate2_transcript,
                genome,
                mapq,
                true, // is_first_mate
                paired_aln.is_proper_pair,
                paired_tlen(
                    &paired_aln.mate1_transcript,
                    &paired_aln.mate2_transcript,
                    true,
                    paired_aln.insert_size,
                    params,
                ),
                max_output, // NH = number of reported alignments
                hit_index,
                params.out_sam_attr_ih_start,
                combined_score,
                attrs,
            )?;
            maybe_insert_rg_tag(&mut rec1, rg_id);
            apply_sam_flag_or_and(&mut rec1, params);
            apply_primary_flag(&mut rec1, combined_score, best_score, params);
            records.push(rec1);

            // Create record for mate2 (this=mate2, mate=mate1)
            let mut rec2 = build_paired_mate_record(
                read_name,
                mate2_seq,
                mate2_qual,
                m2_clip5p,
                m2_clip3p,
                &paired_aln.mate2_transcript,
                &paired_aln.mate1_transcript,
                genome,
                mapq,
                false, // is_first_mate
                paired_aln.is_proper_pair,
                paired_tlen(
                    &paired_aln.mate1_transcript,
                    &paired_aln.mate2_transcript,
                    false,
                    -paired_aln.insert_size, // Negative for mate2 (mode 1 default)
                    params,
                ),
                max_output, // NH = number of reported alignments
                hit_index,
                params.out_sam_attr_ih_start,
                combined_score,
                attrs,
            )?;
            maybe_insert_rg_tag(&mut rec2, rg_id);
            apply_sam_flag_or_and(&mut rec2, params);
            apply_primary_flag(&mut rec2, combined_score, best_score, params);
            records.push(rec2);
        }

        Ok(records)
    }

    /// Build SAM records for a half-mapped pair (one mate mapped, one unmapped).
    ///
    /// Returns 2 records: mate1 first, mate2 second (regardless of which is mapped).
    ///
    /// **Mapped mate:** Normal alignment with FLAG 0x8 (mate unmapped).
    ///   RNEXT = own chr, PNEXT = own pos (STAR convention for unmapped mate).
    ///
    /// **Unmapped mate:** FLAG 0x4, co-located at mapped mate's position.
    ///   SEQ/QUAL in forward orientation (no RC).
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_arguments)]
    pub fn build_half_mapped_records(
        read_name: &str,
        mate1_seq: &[u8],
        mate1_qual: &[u8],
        mate2_seq: &[u8],
        mate2_qual: &[u8],
        m1_clip5p: usize,
        m1_clip3p: usize,
        m2_clip5p: usize,
        m2_clip3p: usize,
        mapped_transcript: &Transcript,
        mate1_is_mapped: bool,
        genome: &Genome,
        params: &Parameters,
        n_for_mapq: usize,
    ) -> Result<Vec<RecordBuf>, Error> {
        let mut records = Vec::with_capacity(2);

        let n_alignments = 1usize;
        let effective_n = n_alignments.max(n_for_mapq);
        let mapq = calculate_mapq(effective_n, params.out_sam_mapq_unique);
        let attrs = params.out_sam_attributes;
        let rg_id_owned = params.primary_rg_id()?;
        let rg_id = rg_id_owned.as_deref();

        // Compute mapped mate's per-chr position for co-location
        let chr_start = genome.chr_start[mapped_transcript.chr_idx];
        let mapped_pos = (mapped_transcript.genome_start - chr_start + 1) as usize;

        // Determine which sequences (and clip amounts) go where. mate seqs are the
        // full original reads; the mapped mate's core is built from its aligned slice
        // and soft-clipped below, the unmapped mate keeps the full read.
        let (mapped_seq, mapped_qual, unmapped_seq, unmapped_qual, mapped_clip5p, mapped_clip3p) =
            if mate1_is_mapped {
                (
                    mate1_seq, mate1_qual, mate2_seq, mate2_qual, m1_clip5p, m1_clip3p,
                )
            } else {
                (
                    mate2_seq, mate2_qual, mate1_seq, mate1_qual, m2_clip5p, m2_clip3p,
                )
            };
        // Aligned slice of the mapped mate (core SEQ/CIGAR/MD built from this).
        let aligned_mapped_seq = &mapped_seq[mapped_clip5p..mapped_seq.len() - mapped_clip3p];
        let aligned_mapped_qual = &mapped_qual[mapped_clip5p..mapped_qual.len() - mapped_clip3p];

        // --- Build mapped mate record ---
        let mut mapped_rec = RecordBuf::default();
        mapped_rec.name_mut().replace(read_name.into());

        let mut mapped_flags = sam::alignment::record::Flags::SEGMENTED // 0x1
            | sam::alignment::record::Flags::MATE_UNMAPPED; // 0x8
        if mapped_transcript.is_reverse {
            mapped_flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED; // 0x10
        }
        if mate1_is_mapped {
            mapped_flags |= sam::alignment::record::Flags::FIRST_SEGMENT; // 0x40
        } else {
            mapped_flags |= sam::alignment::record::Flags::LAST_SEGMENT; // 0x80
        }
        *mapped_rec.flags_mut() = mapped_flags;

        *mapped_rec.reference_sequence_id_mut() = Some(mapped_transcript.chr_idx);
        *mapped_rec.alignment_start_mut() = Some(
            mapped_pos
                .try_into()
                .map_err(|e| Error::Alignment(format!("invalid position {mapped_pos}: {e}")))?,
        );
        *mapped_rec.mapping_quality_mut() = MappingQuality::new(mapq);
        *mapped_rec.cigar_mut() = mapped_transcript.cigar.iter().copied().collect();

        // RNEXT = own chr, PNEXT = own pos (STAR convention for unmapped mate)
        *mapped_rec.mate_reference_sequence_id_mut() = Some(mapped_transcript.chr_idx);
        *mapped_rec.mate_alignment_start_mut() =
            Some(mapped_pos.try_into().map_err(|e| {
                Error::Alignment(format!("invalid mate position {mapped_pos}: {e}"))
            })?);
        *mapped_rec.template_length_mut() = 0;

        // SEQ/QUAL (core from the aligned slice; apply_read_clips restores full read)
        if mapped_transcript.is_reverse {
            let seq_bytes: Vec<u8> = aligned_mapped_seq
                .iter()
                .rev()
                .map(|&b| decode_base(complement_base(b)))
                .collect();
            *mapped_rec.sequence_mut() = Sequence::from(seq_bytes);
            let mut qual = fastq_qual_to_phred(aligned_mapped_qual);
            qual.reverse();
            *mapped_rec.quality_scores_mut() = QualityScores::from(qual);
        } else {
            let seq_bytes: Vec<u8> = aligned_mapped_seq.iter().map(|&b| decode_base(b)).collect();
            *mapped_rec.sequence_mut() = Sequence::from(seq_bytes);
            *mapped_rec.quality_scores_mut() =
                QualityScores::from(fastq_qual_to_phred(aligned_mapped_qual));
        }

        // Optional tags on mapped mate
        let data = mapped_rec.data_mut();
        if attrs.contains(SamAttributes::NH) {
            data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(n_alignments as i32));
        }
        if attrs.contains(SamAttributes::HI) {
            data.insert(
                Tag::HIT_INDEX,
                Value::from(params.out_sam_attr_ih_start as i32),
            );
        }
        if attrs.contains(SamAttributes::AS) {
            data.insert(Tag::ALIGNMENT_SCORE, Value::from(mapped_transcript.score));
        }
        if attrs.contains(SamAttributes::NMM) {
            data.insert(
                Tag::new(b'n', b'M'),
                Value::from(mapped_transcript.n_mismatch as i32),
            );
        }
        if attrs.contains(SamAttributes::NM) {
            data.insert(
                Tag::EDIT_DISTANCE,
                Value::from(sam_spec_nm(
                    mapped_transcript.n_mismatch,
                    &mapped_transcript.cigar,
                )),
            );
        }
        if attrs.contains(SamAttributes::XS)
            && let Some(xs_strand) = derive_xs_strand(mapped_transcript)
        {
            data.insert(Tag::new(b'X', b'S'), Value::Character(xs_strand as u8));
        }
        if attrs.contains(SamAttributes::JM)
            && let Some(jm) = build_jm_tag(mapped_transcript)
        {
            data.insert(Tag::new(b'j', b'M'), jm);
        }
        if attrs.contains(SamAttributes::JI)
            && let Some(ji) = build_ji_tag(mapped_transcript, chr_start)
        {
            data.insert(Tag::new(b'j', b'I'), ji);
        }
        if attrs.contains(SamAttributes::MD) {
            let md = build_md_tag(
                mapped_transcript,
                aligned_mapped_seq,
                genome,
                mapped_transcript.is_reverse,
            );
            data.insert(Tag::new(b'M', b'D'), Value::String(BString::from(md)));
        }
        maybe_insert_rg_tag(&mut mapped_rec, rg_id);
        apply_sam_flag_or_and(&mut mapped_rec, params);
        // Soft-clip the mapped mate's clipped bases (STAR convention), matching the
        // both-mapped path. No-op when the mapped mate has no clip.
        if mapped_clip5p > 0 || mapped_clip3p > 0 {
            apply_read_clips(
                &mut mapped_rec,
                &mapped_transcript.cigar,
                mapped_seq,
                mapped_qual,
                mapped_clip5p,
                mapped_clip3p,
                mapped_transcript.is_reverse,
            );
        }

        // --- Build unmapped mate record ---
        let mut unmapped_rec = RecordBuf::default();
        unmapped_rec.name_mut().replace(read_name.into());

        let mut unmapped_flags = sam::alignment::record::Flags::SEGMENTED // 0x1
            | sam::alignment::record::Flags::UNMAPPED; // 0x4
        // Mate reverse flag from mapped mate's strand
        if mapped_transcript.is_reverse {
            unmapped_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED; // 0x20
        }
        if mate1_is_mapped {
            // Unmapped is mate2
            unmapped_flags |= sam::alignment::record::Flags::LAST_SEGMENT; // 0x80
        } else {
            // Unmapped is mate1
            unmapped_flags |= sam::alignment::record::Flags::FIRST_SEGMENT; // 0x40
        }
        *unmapped_rec.flags_mut() = unmapped_flags;

        // Co-locate unmapped mate at mapped mate's position
        *unmapped_rec.reference_sequence_id_mut() = Some(mapped_transcript.chr_idx);
        *unmapped_rec.alignment_start_mut() = Some(
            mapped_pos
                .try_into()
                .map_err(|e| Error::Alignment(format!("invalid position {mapped_pos}: {e}")))?,
        );
        *unmapped_rec.mapping_quality_mut() = MappingQuality::new(0);
        // CIGAR = * (default empty cigar)
        // RNEXT = mapped mate's chr
        *unmapped_rec.mate_reference_sequence_id_mut() = Some(mapped_transcript.chr_idx);
        *unmapped_rec.mate_alignment_start_mut() =
            Some(mapped_pos.try_into().map_err(|e| {
                Error::Alignment(format!("invalid mate position {mapped_pos}: {e}"))
            })?);
        *unmapped_rec.template_length_mut() = 0;

        // SEQ/QUAL: forward orientation (no RC for unmapped)
        let unmapped_seq_bytes: Vec<u8> = unmapped_seq.iter().map(|&b| decode_base(b)).collect();
        *unmapped_rec.sequence_mut() = Sequence::from(unmapped_seq_bytes);
        *unmapped_rec.quality_scores_mut() =
            QualityScores::from(fastq_qual_to_phred(unmapped_qual));
        maybe_insert_rg_tag(&mut unmapped_rec, rg_id);
        insert_unmapped_tags(&mut unmapped_rec, attrs, UnmappedReason::Other);

        // Order: mate1 first, mate2 second
        if mate1_is_mapped {
            records.push(mapped_rec);
            records.push(unmapped_rec);
        } else {
            records.push(unmapped_rec);
            records.push(mapped_rec);
        }

        Ok(records)
    }

    /// Build transcriptome-space SAM records for `--quantMode TranscriptomeSAM`.
    ///
    /// Each projected `Transcript` is converted to a record where:
    ///   * `chr_idx` is the transcript index (matches the transcriptome
    ///     header's @SQ order),
    ///   * `genome_start` is the 0-based transcript-space position (→ POS =
    ///     t-space_pos + 1),
    ///   * splice-aware tags (`jM`, `jI`, `XS`) are not emitted (splices
    ///     collapse in t-space and have no meaning there),
    ///   * standard tags (`NH`, `HI`, `AS`, `NM`/`nM`, `MD`) are emitted per
    ///     the `--outSAMattributes` set.
    ///
    /// `primary_hit_idx` (0-based) is the projected alignment selected as
    /// primary (randomly among ties per STAR's `rngUniformReal0to1`).  All
    /// other records get the SECONDARY flag (0x100).
    #[allow(clippy::too_many_arguments)]
    pub fn build_transcriptome_records(
        read_name: &str,
        read_seq: &[u8],
        read_qual: &[u8],
        projected: &[Transcript],
        mapq: u8,
        params: &Parameters,
        primary_hit_idx: usize,
    ) -> Result<Vec<RecordBuf>, Error> {
        if projected.is_empty() {
            return Ok(Vec::new());
        }
        // Splice tags are meaningless in t-space; MD would require the
        // transcript's t-space reference which we do not precompute, and
        // STAR also does not emit MD for transcriptome SAM.
        let attrs = params.out_sam_attributes
            - (SamAttributes::JM | SamAttributes::JI | SamAttributes::XS | SamAttributes::MD);

        let n_alignments = projected.len();
        let mut records = Vec::with_capacity(n_alignments);

        let rg_id_owned = params.primary_rg_id()?;
        let rg_id = rg_id_owned.as_deref();

        for (hit_idx, t) in projected.iter().enumerate() {
            let mut record = RecordBuf::default();
            record.name_mut().replace(read_name.into());

            // FLAGS: SECONDARY if not the primary; REVERSE if is_reverse.
            let mut flags = sam::alignment::record::Flags::empty();
            if t.is_reverse {
                flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
            }
            if hit_idx != primary_hit_idx {
                flags |= sam::alignment::record::Flags::SECONDARY;
            }
            *record.flags_mut() = flags;

            // RNAME = transcript index (maps to transcriptome header).
            *record.reference_sequence_id_mut() = Some(t.chr_idx);

            // POS = t-space position + 1 (1-based).
            let pos = (t.genome_start + 1) as usize;
            *record.alignment_start_mut() =
                Some(pos.try_into().map_err(|e| {
                    Error::Alignment(format!("invalid t-space position {pos}: {e}"))
                })?);

            // MAPQ
            *record.mapping_quality_mut() = MappingQuality::new(mapq);

            // CIGAR (already has N ops stripped by align_to_transcripts)
            *record.cigar_mut() = t.cigar.iter().copied().collect();

            // SEQ / QUAL — STAR writes the original-orientation sequence when
            // FLAG 0x10 is unset (forward alignment in t-space) and RC'd seq
            // when 0x10 is set.  We follow SAM spec: SEQ matches the CIGAR's
            // read orientation, so we mirror `transcript_to_record`.
            if t.is_reverse {
                let seq_bytes: Vec<u8> = read_seq
                    .iter()
                    .rev()
                    .map(|&b| decode_base(complement_base(b)))
                    .collect();
                *record.sequence_mut() = Sequence::from(seq_bytes);
                let mut qual = fastq_qual_to_phred(read_qual);
                qual.reverse();
                *record.quality_scores_mut() = QualityScores::from(qual);
            } else {
                let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
                *record.sequence_mut() = Sequence::from(seq_bytes);
                *record.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(read_qual));
            }

            // Optional tags
            let data = record.data_mut();
            if attrs.contains(SamAttributes::NH) {
                data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(n_alignments as i32));
            }
            if attrs.contains(SamAttributes::HI) {
                // HI defaults to 1-based (primary = 1, secondaries > 1 in emission order);
                // `--outSAMattrIHstart` shifts the whole sequence (0 = CellRanger convention).
                data.insert(
                    Tag::HIT_INDEX,
                    Value::from(hit_idx as i32 + params.out_sam_attr_ih_start as i32),
                );
            }
            if attrs.contains(SamAttributes::AS) {
                data.insert(Tag::ALIGNMENT_SCORE, Value::from(t.score));
            }
            if attrs.contains(SamAttributes::NMM) {
                data.insert(Tag::new(b'n', b'M'), Value::from(t.n_mismatch as i32));
            }
            if attrs.contains(SamAttributes::NM) {
                data.insert(
                    Tag::EDIT_DISTANCE,
                    Value::from(sam_spec_nm(t.n_mismatch, &t.cigar)),
                );
            }

            maybe_insert_rg_tag(&mut record, rg_id);

            records.push(record);
        }

        Ok(records)
    }

    /// Build unmapped paired records (both mates unmapped)
    pub fn build_paired_unmapped_records(
        read_name: &str,
        mate1_seq: &[u8],
        mate1_qual: &[u8],
        mate2_seq: &[u8],
        mate2_qual: &[u8],
        params: &Parameters,
        unmapped_reason: UnmappedReason,
    ) -> Result<Vec<RecordBuf>, Error> {
        let mut records = Vec::with_capacity(2);
        let rg_id_owned = params.primary_rg_id()?;
        let rg_id = rg_id_owned.as_deref();
        let attrs = params.out_sam_attributes;

        // Mate1 record
        let mut rec1 = RecordBuf::default();
        rec1.name_mut().replace(read_name.into());

        // FLAGS: 0x1 (paired) | 0x4 (unmapped) | 0x8 (mate unmapped) | 0x40 (first in pair)
        let flags1 = sam::alignment::record::Flags::SEGMENTED
            | sam::alignment::record::Flags::UNMAPPED
            | sam::alignment::record::Flags::MATE_UNMAPPED
            | sam::alignment::record::Flags::FIRST_SEGMENT;
        *rec1.flags_mut() = flags1;

        let seq1_bytes: Vec<u8> = mate1_seq.iter().map(|&b| decode_base(b)).collect();
        *rec1.sequence_mut() = Sequence::from(seq1_bytes);
        *rec1.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(mate1_qual));
        maybe_insert_rg_tag(&mut rec1, rg_id);
        insert_unmapped_tags(&mut rec1, attrs, unmapped_reason);
        records.push(rec1);

        // Mate2 record
        let mut rec2 = RecordBuf::default();
        rec2.name_mut().replace(read_name.into());

        // FLAGS: 0x1 (paired) | 0x4 (unmapped) | 0x8 (mate unmapped) | 0x80 (second in pair)
        let flags2 = sam::alignment::record::Flags::SEGMENTED
            | sam::alignment::record::Flags::UNMAPPED
            | sam::alignment::record::Flags::MATE_UNMAPPED
            | sam::alignment::record::Flags::LAST_SEGMENT;
        *rec2.flags_mut() = flags2;

        let seq2_bytes: Vec<u8> = mate2_seq.iter().map(|&b| decode_base(b)).collect();
        *rec2.sequence_mut() = Sequence::from(seq2_bytes);
        *rec2.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(mate2_qual));
        maybe_insert_rg_tag(&mut rec2, rg_id);
        insert_unmapped_tags(&mut rec2, attrs, unmapped_reason);
        records.push(rec2);

        Ok(records)
    }
}

/// SAM writer that streams to stdout.
pub struct SamStdoutWriter {
    writer: sam::io::Writer<BufWriter<std::io::Stdout>>,
    header: sam::Header,
}

impl SamStdoutWriter {
    pub fn create(genome: &Genome, params: &Parameters) -> Result<Self, Error> {
        let header = build_sam_header(genome, params)?;
        let mut writer = sam::io::Writer::new(BufWriter::new(std::io::stdout()));
        writer.write_header(&header)?;
        Ok(Self { writer, header })
    }

    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        for record in batch {
            self.writer.write_alignment_record(&self.header, record)?;
        }
        Ok(())
    }
}

/// Build paired SAM header from genome
pub fn build_sam_header(genome: &Genome, params: &Parameters) -> Result<sam::Header, Error> {
    build_sam_header_from_refs(
        (0..genome.n_chr_real)
            .map(|i| (genome.chr_name[i].as_str(), genome.chr_length[i] as usize)),
        params,
    )
}

/// Create a SAM writer for BySJout disk-buffering (temp file). Returns (header, writer).
pub fn create_bysj_writer(
    file: std::fs::File,
    genome: &Genome,
    params: &Parameters,
) -> Result<(sam::Header, sam::io::Writer<BufWriter<std::fs::File>>), Error> {
    let header = build_sam_header(genome, params)?;
    let mut writer = sam::io::Writer::new(BufWriter::new(file));
    writer.write_header(&header)?;
    Ok((header, writer))
}

/// Write a slice of RecordBuf to a SAM writer (for BySJout temp file).
pub fn bysj_write_records<W: std::io::Write>(
    writer: &mut sam::io::Writer<W>,
    header: &sam::Header,
    records: &[RecordBuf],
) -> Result<(), Error> {
    for rec in records {
        writer.write_alignment_record(header, rec)?;
    }
    Ok(())
}

/// Read exactly `n` records from a SAM reader. If `collect` is true, return them in a Vec;
/// otherwise just advance the reader position (discard records).
pub fn bysj_read_n_records<R: std::io::BufRead>(
    reader: &mut sam::io::Reader<R>,
    header: &sam::Header,
    n: u32,
    collect: bool,
) -> Result<Vec<RecordBuf>, Error> {
    let mut out = if collect {
        Vec::with_capacity(n as usize)
    } else {
        Vec::new()
    };
    let mut buf = RecordBuf::default();
    for _ in 0..n {
        reader.read_record_buf(header, &mut buf)?;
        if collect {
            out.push(buf.clone());
        }
    }
    Ok(out)
}

/// Build a SAM header from an iterator of (name, length) reference pairs.
///
/// Used both for the genome header (chromosomes) and the transcriptome header
/// (one @SQ per transcript, length = transcript-space length).
pub fn build_sam_header_from_refs<'a, I>(refs: I, params: &Parameters) -> Result<sam::Header, Error>
where
    I: IntoIterator<Item = (&'a str, usize)>,
{
    let mut builder = sam::Header::builder();

    // @HD line. `--outSAMheaderHD` replaces it wholesale, given as the
    // tab-separated fields STAR expects (`@HD VN:1.4 SO:coordinate`).
    if params.out_sam_header_hd.is_empty() {
        builder = builder.set_header(Map::default());
    } else {
        let mut hd = Map::<sam::header::record::value::map::Header>::default();
        for field in &params.out_sam_header_hd {
            let field = field.trim();
            // STAR takes the leading `@HD` as part of the value list; ignore it.
            if field.is_empty() || field == "@HD" {
                continue;
            }
            let Some((tag, value)) = field.split_once(':') else {
                return Err(Error::Parameter(format!(
                    "--outSAMheaderHD field '{field}' is not TAG:value"
                )));
            };
            if tag.len() != 2 {
                return Err(Error::Parameter(format!(
                    "--outSAMheaderHD tag '{tag}' is not two characters"
                )));
            }
            let tag_bytes: [u8; 2] = tag.as_bytes()[..2].try_into().unwrap();
            let other_tag: HeaderOtherTag<_> = HeaderOtherTag::try_from(tag_bytes)
                .map_err(|e| Error::Parameter(format!("invalid @HD tag '{tag}': {e}")))?;
            hd.other_fields_mut().insert(other_tag, value.into());
        }
        builder = builder.set_header(hd);
    }

    // @SQ lines for each reference
    for (name, length) in refs {
        let length_nz = NonZeroUsize::new(length)
            .ok_or_else(|| Error::Index(format!("reference {name} has zero length")))?;

        builder = builder.add_reference_sequence(
            name,
            Map::<sam::header::record::value::map::ReferenceSequence>::new(length_nz),
        );
    }

    // @RG lines from --outSAMattrRGline. When multiple input files share the
    // same RG ID, only emit one @RG line.
    let rg_lines = params.parsed_rg_lines()?;
    let mut seen_ids: HashSet<String> = HashSet::new();
    for line in &rg_lines {
        let mut fields = line.split('\t');
        let id = fields
            .next()
            .and_then(|f| f.strip_prefix("ID:"))
            .ok_or_else(|| Error::Parameter(format!("malformed RG line '{line}'")))?;
        if !seen_ids.insert(id.to_string()) {
            continue;
        }
        let mut map = Map::<ReadGroup>::default();
        for field in fields {
            if field.len() < 3 || &field[2..3] != ":" {
                return Err(Error::Parameter(format!(
                    "RG field '{field}' is not TAG:value"
                )));
            }
            let tag_bytes: [u8; 2] = field.as_bytes()[..2].try_into().unwrap();
            let other_tag: HeaderOtherTag<_> =
                HeaderOtherTag::try_from(tag_bytes).map_err(|e| {
                    Error::Parameter(format!("invalid RG tag '{}': {}", &field[..2], e))
                })?;
            map.other_fields_mut().insert(other_tag, field[3..].into());
        }
        builder = builder.add_read_group(id, map);
    }

    let mut pg = Map::<Program>::default();
    pg.other_fields_mut()
        .insert(program_tag::NAME, BString::from("rustar-aligner"));
    pg.other_fields_mut().insert(
        program_tag::VERSION,
        BString::from(env!("CARGO_PKG_VERSION")),
    );
    let cl = params
        .command_line
        .clone()
        .unwrap_or_else(|| "rustar-aligner".to_string());
    pg.other_fields_mut()
        .insert(program_tag::COMMAND_LINE, BString::from(cl));
    builder = builder.add_program("rustar-aligner", pg);

    // `--outSAMheaderPG` adds one further @PG line, given the same way.
    if !params.out_sam_header_pg.is_empty() {
        let mut extra = Map::<Program>::default();
        let mut id: Option<String> = None;
        for field in &params.out_sam_header_pg {
            let field = field.trim();
            if field.is_empty() || field == "@PG" {
                continue;
            }
            let Some((tag, value)) = field.split_once(':') else {
                return Err(Error::Parameter(format!(
                    "--outSAMheaderPG field '{field}' is not TAG:value"
                )));
            };
            if tag == "ID" {
                id = Some(value.to_string());
                continue;
            }
            if tag.len() != 2 {
                return Err(Error::Parameter(format!(
                    "--outSAMheaderPG tag '{tag}' is not two characters"
                )));
            }
            let tag_bytes: [u8; 2] = tag.as_bytes()[..2].try_into().unwrap();
            let other_tag: HeaderOtherTag<_> = HeaderOtherTag::try_from(tag_bytes)
                .map_err(|e| Error::Parameter(format!("invalid @PG tag '{tag}': {e}")))?;
            extra.other_fields_mut().insert(other_tag, value.into());
        }
        let id =
            id.ok_or_else(|| Error::Parameter("--outSAMheaderPG needs an ID: field".to_string()))?;
        builder = builder.add_program(id, extra);
    }

    // `--outSAMheaderCommentFile` contributes one @CO line per line of the
    // named file. `-` (the default) means no comments.
    if params.out_sam_header_comment_file != "-" {
        let path = std::path::Path::new(&params.out_sam_header_comment_file);
        let contents = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
        for line in contents.lines() {
            if !line.is_empty() {
                builder = builder.add_comment(line);
            }
        }
    }

    Ok(builder.build())
}

/// Insert `RG:Z:<id>` on the record when an ID is set. `Parameters::try_parse_from`
/// auto-ORs `SamAttributes::RG` into `out_sam_attributes` whenever an RG line is
/// configured, so `rg_id.is_some()` implies the attribute is wanted.
fn maybe_insert_rg_tag(record: &mut RecordBuf, rg_id: Option<&str>) {
    if let Some(id) = rg_id {
        record
            .data_mut()
            .insert(Tag::READ_GROUP, Value::String(BString::from(id)));
    }
}

/// Add STARsolo GX/GN gene tags to each record for the read's `Gene`-feature
/// assignment. `gx`/`gn` are the gene_id / gene_name, or `"-"` when the read is
/// not uniquely assigned to a gene (STARsolo convention). No-ops unless the
/// corresponding attribute bit is set in `--outSAMattributes`.
pub fn add_gene_tags(records: &mut [RecordBuf], gx: &str, gn: &str, attrs: SamAttributes) {
    for rec in records.iter_mut() {
        if attrs.contains(SamAttributes::GX) {
            rec.data_mut()
                .insert(Tag::new(b'G', b'X'), Value::String(BString::from(gx)));
        }
        if attrs.contains(SamAttributes::GN) {
            rec.data_mut()
                .insert(Tag::new(b'G', b'N'), Value::String(BString::from(gn)));
        }
    }
}

/// The per-read STARsolo barcode tag values of one read
/// (`ReadAlign_alignBAM.cpp:389-473`).
///
/// `CB` is filled here only when the barcode is corrected as the read is
/// processed (`--soloType CB_samTagOut`); in a counting run it, and `UB`: are
/// added when the sorted BAM is written, once UMI collapsing has run.
pub struct SoloBarcodeTagValues<'a> {
    /// `CR`: raw cell barcode.
    pub cb_seq: &'a str,
    /// `CY`: raw cell-barcode quality.
    pub cb_qual: &'a str,
    /// `UR`: raw UMI.
    pub umi_seq: &'a str,
    /// `UY`: raw UMI quality.
    pub umi_qual: &'a str,
    /// `CB`: corrected cell barcode, when known at read time.
    pub cb_corrected: Option<&'a str>,
    /// `sM`: STAR's `cbMatch` code.
    pub cb_match: i32,
    /// `sS`: whole barcode read.
    pub barcode_seq: &'a str,
    /// `sQ`: whole barcode-read quality.
    pub barcode_qual: &'a str,
}

/// Add the STARsolo per-read barcode tags to every record of one read (mapped
/// alignments and the unmapped record alike, as in STAR). `attrs` should already
/// be narrowed to the tags this run emits, see `Parameters::solo_sam_tags`.
pub fn add_solo_barcode_tags(
    records: &mut [RecordBuf],
    values: &SoloBarcodeTagValues<'_>,
    attrs: SamAttributes,
) {
    if !attrs.intersects(SamAttributes::SOLO_TAGS) {
        return;
    }
    let str_tags: [(SamAttributes, [u8; 2], &str); 5] = [
        (SamAttributes::CR, *b"CR", values.cb_seq),
        (SamAttributes::CY, *b"CY", values.cb_qual),
        (SamAttributes::UR, *b"UR", values.umi_seq),
        (SamAttributes::UY, *b"UY", values.umi_qual),
        (SamAttributes::SS, *b"sS", values.barcode_seq),
    ];
    for rec in records.iter_mut() {
        // An empty value means the read carries no such sequence at all (a
        // barcode read too short to hold a CB+UMI); STAR has no empty tags, so
        // the tag is left off rather than written blank.
        for (flag, tag, value) in str_tags {
            if attrs.contains(flag) && !value.is_empty() {
                rec.data_mut().insert(
                    Tag::new(tag[0], tag[1]),
                    Value::String(BString::from(value)),
                );
            }
        }
        // sQ is a quality string, written verbatim like sS.
        if attrs.contains(SamAttributes::SQ) && !values.barcode_qual.is_empty() {
            rec.data_mut().insert(
                Tag::new(b's', b'Q'),
                Value::String(BString::from(values.barcode_qual)),
            );
        }
        if attrs.contains(SamAttributes::SM) {
            rec.data_mut()
                .insert(Tag::new(b's', b'M'), Value::Int32(values.cb_match));
        }
        if let (true, Some(cb)) = (attrs.contains(SamAttributes::CB), values.cb_corrected) {
            rec.data_mut()
                .insert(Tag::new(b'C', b'B'), Value::String(BString::from(cb)));
        }
    }
}

/// Add the per-alignment STARsolo gene tags (`gx`, `gn`, `sF`) to a read's
/// records.
///
/// `tags` is one entry per alignment; `records_per_align` is how many records
/// each alignment produced (1 single-end, 2 for a mate pair), so both mates of a
/// pair carry the pair's genes. Records past the end of `tags` (there are none
/// in practice) are left alone.
pub fn add_align_gene_tags(
    records: &mut [RecordBuf],
    tags: &[crate::solo::AlignGeneTag],
    records_per_align: usize,
    attrs: SamAttributes,
) {
    use noodles::sam::alignment::record_buf::data::field::value::Array;

    if !attrs.intersects(SamAttributes::GXM | SamAttributes::GNM | SamAttributes::SF)
        || records_per_align == 0
    {
        return;
    }
    for (i, rec) in records.iter_mut().enumerate() {
        let Some(tag) = tags.get(i / records_per_align) else {
            continue;
        };
        if attrs.contains(SamAttributes::GXM) {
            rec.data_mut().insert(
                Tag::new(b'g', b'x'),
                Value::String(BString::from(tag.gx.as_str())),
            );
        }
        if attrs.contains(SamAttributes::GNM) {
            rec.data_mut().insert(
                Tag::new(b'g', b'n'),
                Value::String(BString::from(tag.gn.as_str())),
            );
        }
        if attrs.contains(SamAttributes::SF) {
            rec.data_mut().insert(
                Tag::new(b's', b'F'),
                Value::Array(Array::Int32(tag.sf.to_vec())),
            );
        }
    }
}

/// Private aux tag carrying the input read index on a buffered record until the
/// sorted BAM is written, where it becomes `CB`/`UB`. STAR does the same thing
/// by encoding `iReadAll` in the record's trailing bytes
/// (`SoloFeature_addBAMtags.cpp:8`); a local `z`-namespace tag survives sorting
/// without a parallel array.
const SOLO_READ_INDEX_TAG: [u8; 2] = *b"zR";

/// Stamp the input read index on each of a read's records, for the CB/UB pass.
pub fn add_solo_read_index(records: &mut [RecordBuf], read_index: u32) {
    for rec in records.iter_mut() {
        rec.data_mut().insert(
            Tag::new(SOLO_READ_INDEX_TAG[0], SOLO_READ_INDEX_TAG[1]),
            Value::UInt32(read_index),
        );
    }
}

/// Replace the private read-index tag with `CB`/`UB`, read out of STAR's
/// readInfo once UMI collapsing has run.
///
/// Both tags are written whenever either was requested, and both fall back to
/// `"-"`, exactly as `SoloFeature::addBAMtags` does: a read that was not counted
/// (no whitelist cell, no valid UMI, no gene) has no cell or molecule to name.
pub fn apply_solo_read_info(
    records: &mut [RecordBuf],
    read_info: &[crate::solo::ReadInfo],
    whitelist: &crate::solo::CbWhitelist,
    umi_len: usize,
) {
    let tag = Tag::new(SOLO_READ_INDEX_TAG[0], SOLO_READ_INDEX_TAG[1]);
    for rec in records.iter_mut() {
        let Some(Value::UInt32(read_index)) = rec.data().get(&tag).cloned() else {
            continue;
        };
        rec.data_mut().remove(&tag);
        let info = read_info
            .get(read_index as usize)
            .copied()
            .unwrap_or_default();
        let cb = (info.cb != u32::MAX)
            .then(|| whitelist.barcode_string(info.cb))
            .flatten()
            .unwrap_or_else(|| "-".to_string());
        let ub = if info.umi == u64::MAX {
            "-".to_string()
        } else {
            crate::solo::whitelist::unpack_barcode(info.umi, umi_len)
        };
        rec.data_mut()
            .insert(Tag::new(b'C', b'B'), Value::String(BString::from(cb)));
        rec.data_mut()
            .insert(Tag::new(b'U', b'B'), Value::String(BString::from(ub)));
    }
}

/// Apply `--outSAMflagOR` / `--outSAMflagAND` to a mapped record's FLAG:
/// `(FLAG & flagAND) | flagOR`. Matches STAR/STAR-rs, which apply this only to
/// mapped-mate records; unmapped and transcriptome-BAM records are untouched.
fn apply_sam_flag_or_and(record: &mut RecordBuf, params: &Parameters) {
    let or_bits = params.out_sam_flag_or as u16;
    let and_bits = params.out_sam_flag_and.min(u16::MAX as u32) as u16;
    let bits = u16::from(record.flags());
    *record.flags_mut() = sam::alignment::record::Flags::from((bits & and_bits) | or_bits);
}

/// `--outSAMprimaryFlag AllBestScore`: clear the SECONDARY bit on every alignment tied for the
/// best score, instead of only the single (already-sorted) best alignment (`OneBestScore`,
/// the default, which the caller's existing hit-index-based SECONDARY assignment already gives).
fn apply_primary_flag(record: &mut RecordBuf, score: i32, best_score: i32, params: &Parameters) {
    if params.out_sam_primary_flag == crate::params::OutSamPrimaryFlag::AllBestScore
        && score == best_score
    {
        let mut flags = record.flags();
        flags.remove(sam::alignment::record::Flags::SECONDARY);
        *record.flags_mut() = flags;
    }
}

/// `--outSAMtlen 2`: per-mate span, signed by whichever mate is genomically leftmost (ties go to
/// mate1). Differs from the default (mode 1, `default_signed`) only for contained/dovetailed
/// pairs, where mode 1's whole-combined-transcript span and mode 2's per-mate max/min span
/// diverge.
fn paired_tlen(
    mate1: &Transcript,
    mate2: &Transcript,
    is_first_mate: bool,
    default_signed: i32,
    params: &Parameters,
) -> i32 {
    if params.out_sam_tlen != 2 {
        return default_signed;
    }
    let span = (mate1.genome_end.max(mate2.genome_end) - mate1.genome_start.min(mate2.genome_start))
        as i64;
    let mate1_is_leftmost = mate1.genome_start <= mate2.genome_start;
    let this_is_leftmost = mate1_is_leftmost == is_first_mate;
    let span = span as i32;
    if this_is_leftmost { span } else { -span }
}

/// Convert FASTQ ASCII quality bytes (Phred+33) to raw Phred values (0-93) for
/// the BAM binary QUAL field, per SAM spec §4.2.3.
///
/// `saturating_sub` clamps malformed bytes < 33 to 0 rather than underflowing.
fn fastq_qual_to_phred(qual: &[u8]) -> Vec<u8> {
    qual.iter().map(|&b| b.saturating_sub(33)).collect()
}

/// Re-attach clipped bases to a mapped record as soft-clips, matching STAR (and
/// STAR-rs): the record is first built from the *clipped* read (so CIGAR core, MD
/// and NM cover only the aligned bases), then this restores the full original read
/// as SEQ/QUAL and wraps the CIGAR with soft-clips for the clipped ends.
///
/// `clip5p`/`clip3p` are in read 5'/3' coordinates. In SAM (genome-forward)
/// orientation the read-5' clip is the leading soft-clip and the read-3' clip the
/// trailing one; on a reverse-strand record the two swap. A clip adjacent to an
/// existing (genomic) soft-clip is merged into a single `S` op (valid SAM). POS is
/// unaffected — soft-clips consume read bases but no reference.
fn apply_read_clips(
    record: &mut RecordBuf,
    core_cigar: &[cigar::Op],
    orig_seq: &[u8],
    orig_qual: &[u8],
    clip5p: usize,
    clip3p: usize,
    is_reverse: bool,
) {
    use cigar::op::Kind;
    // Full original read as SEQ/QUAL (strand-aware), replacing the clipped SEQ.
    if is_reverse {
        let seq: Vec<u8> = orig_seq
            .iter()
            .rev()
            .map(|&b| decode_base(complement_base(b)))
            .collect();
        *record.sequence_mut() = Sequence::from(seq);
        let mut q = fastq_qual_to_phred(orig_qual);
        q.reverse();
        *record.quality_scores_mut() = QualityScores::from(q);
    } else {
        let seq: Vec<u8> = orig_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq);
        *record.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(orig_qual));
    }
    // Wrap the core CIGAR with soft-clips (merging into an adjacent S if present).
    let (lead, trail) = if is_reverse {
        (clip3p, clip5p)
    } else {
        (clip5p, clip3p)
    };
    let mut ops: Vec<cigar::Op> = Vec::with_capacity(core_cigar.len() + 2);
    ops.extend_from_slice(core_cigar);
    if lead > 0 {
        if ops.first().map(|o| o.kind()) == Some(Kind::SoftClip) {
            ops[0] = cigar::Op::new(Kind::SoftClip, ops[0].len() + lead);
        } else {
            ops.insert(0, cigar::Op::new(Kind::SoftClip, lead));
        }
    }
    if trail > 0 {
        let li = ops.len() - 1;
        if ops[li].kind() == Kind::SoftClip {
            ops[li] = cigar::Op::new(Kind::SoftClip, ops[li].len() + trail);
        } else {
            ops.push(cigar::Op::new(Kind::SoftClip, trail));
        }
    }
    *record.cigar_mut() = ops.into_iter().collect();
}

/// Convert Transcript to SAM record
#[allow(clippy::too_many_arguments)]
fn transcript_to_record(
    transcript: &Transcript,
    read_name: &str,
    read_seq: &[u8],
    read_qual: &[u8],
    genome: &Genome,
    mapq: u8,
    n_alignments: usize,
    hit_index: usize,
    ih_start: u32,
    attrs: SamAttributes,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();

    // Name
    record.name_mut().replace(read_name.into());

    // FLAGS
    let mut flags = sam::alignment::record::Flags::empty();
    if transcript.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
    }
    if hit_index > 1 {
        flags |= sam::alignment::record::Flags::SECONDARY;
    }
    *record.flags_mut() = flags;

    // RNAME (reference sequence name)
    if transcript.chr_idx >= genome.n_chr_real {
        return Err(Error::Alignment(format!(
            "invalid chromosome index {} (max {})",
            transcript.chr_idx,
            genome.n_chr_real - 1
        )));
    }
    *record.reference_sequence_id_mut() = Some(transcript.chr_idx);

    // POS (1-based, per-chromosome coordinate)
    // transcript.genome_start is a global genome coordinate, need to convert to per-chr
    let chr_start = genome.chr_start[transcript.chr_idx];
    let pos = (transcript.genome_start - chr_start + 1) as usize;
    *record.alignment_start_mut() = Some(
        pos.try_into()
            .map_err(|e| Error::Alignment(format!("invalid alignment position {pos}: {e}")))?,
    );

    // MAPQ
    *record.mapping_quality_mut() = MappingQuality::new(mapq);

    // CIGAR
    *record.cigar_mut() = transcript.cigar.iter().copied().collect();

    // Sequence and quality scores
    // Per SAM spec: when FLAG & 16 (reverse strand), SEQ is the reverse complement
    // of the original read, and QUAL is reversed.
    if transcript.is_reverse {
        // Reverse complement the sequence
        let seq_bytes: Vec<u8> = read_seq
            .iter()
            .rev()
            .map(|&b| decode_base(complement_base(b)))
            .collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);

        // Reverse the quality scores
        let mut qual = fastq_qual_to_phred(read_qual);
        qual.reverse();
        *record.quality_scores_mut() = QualityScores::from(qual);
    } else {
        let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(read_qual));
    }

    // Optional tags: gated by --outSAMattributes
    let data = record.data_mut();
    if attrs.contains(SamAttributes::NH) {
        data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(n_alignments as i32));
    }
    if attrs.contains(SamAttributes::HI) {
        data.insert(
            Tag::HIT_INDEX,
            Value::from(hit_index as i32 - 1 + ih_start as i32),
        );
    }
    if attrs.contains(SamAttributes::AS) {
        data.insert(Tag::ALIGNMENT_SCORE, Value::from(transcript.score));
    }
    // nM (mismatch count) before NM (edit distance), matching STAR's tag order.
    if attrs.contains(SamAttributes::NMM) {
        data.insert(
            Tag::new(b'n', b'M'),
            Value::from(transcript.n_mismatch as i32),
        );
    }
    if attrs.contains(SamAttributes::NM) {
        data.insert(
            Tag::EDIT_DISTANCE,
            Value::from(sam_spec_nm(transcript.n_mismatch, &transcript.cigar)),
        );
    }
    if attrs.contains(SamAttributes::XS)
        && let Some(xs_strand) = derive_xs_strand(transcript)
    {
        data.insert(Tag::new(b'X', b'S'), Value::Character(xs_strand as u8));
    }
    if attrs.contains(SamAttributes::JM)
        && let Some(jm) = build_jm_tag(transcript)
    {
        data.insert(Tag::new(b'j', b'M'), jm);
    }
    if attrs.contains(SamAttributes::JI)
        && let Some(ji) = build_ji_tag(transcript, chr_start)
    {
        data.insert(Tag::new(b'j', b'I'), ji);
    }
    if attrs.contains(SamAttributes::MD) {
        let md = build_md_tag(transcript, read_seq, genome, transcript.is_reverse);
        data.insert(Tag::new(b'M', b'D'), Value::String(BString::from(md)));
    }

    Ok(record)
}

/// SAM-spec `NM:i:` edit distance: substitutions + inserted bases + deleted bases.
/// Intron `N` skips and clips are excluded per SAM v1 §1.4.
fn sam_spec_nm(n_mismatch: u32, cigar: &[cigar::Op]) -> i32 {
    use cigar::op::Kind;
    let indel_bases: usize = cigar
        .iter()
        .filter_map(|op| match op.kind() {
            Kind::Insertion | Kind::Deletion => Some(op.len()),
            _ => None,
        })
        .sum();
    (n_mismatch + indel_bases as u32) as i32
}

/// Derive XS strand tag from transcript junction motifs.
/// Returns Some('+') or Some('-') if all junctions agree on strand.
/// Returns None if no junctions, all non-canonical, or conflicting strands.
fn derive_xs_strand(transcript: &Transcript) -> Option<char> {
    let mut strand: Option<char> = None;
    for motif in &transcript.junction_motifs {
        if let Some(s) = motif.implied_strand() {
            match strand {
                None => strand = Some(s),
                Some(prev) if prev != s => return None,
                _ => {}
            }
        }
    }
    strand
}

/// Build jM tag: array of junction motif codes (one per intron/RefSkip in CIGAR).
///
/// Encoding: 0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=GT/AT.
/// Add +20 if junction is annotated in GTF.
fn build_jm_tag(transcript: &Transcript) -> Option<Value> {
    if transcript.junction_motifs.is_empty() {
        return None;
    }
    let motifs: Vec<i8> = transcript
        .junction_motifs
        .iter()
        .zip(
            transcript
                .junction_annotated
                .iter()
                .chain(std::iter::repeat(&false)),
        )
        .map(|(motif, &annotated)| {
            let code = encode_motif(*motif) as i8;
            if annotated { code + 20 } else { code }
        })
        .collect();
    Some(Value::Array(Array::Int8(motifs)))
}

/// Build jI tag: array of intron start/end coordinates (1-based, per-chromosome).
///
/// Format: [start1, end1, start2, end2, ...] where start is first intronic base
/// and end is last intronic base (both 1-based, inclusive).
fn build_ji_tag(transcript: &Transcript, chr_start: u64) -> Option<Value> {
    use cigar::op::Kind;
    if transcript.n_junction == 0 {
        return None;
    }
    let mut coords: Vec<i32> = Vec::new();
    let mut genome_pos = transcript.genome_start;
    for op in &transcript.cigar {
        match op.kind() {
            Kind::Skip => {
                let intron_start = (genome_pos - chr_start + 1) as i32; // 1-based
                let intron_end = (genome_pos + op.len() as u64 - chr_start) as i32; // 1-based inclusive
                coords.push(intron_start);
                coords.push(intron_end);
                genome_pos += op.len() as u64;
            }
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                genome_pos += op.len() as u64;
            }
            _ => {} // Ins, SoftClip, HardClip don't consume reference
        }
    }
    Some(Value::Array(Array::Int32(coords)))
}

/// Build MD tag string: matches/mismatches/deletions relative to reference.
///
/// Format: "10A5^AC6" = 10 match, A mismatch, 5 match, 2bp deletion (AC), 6 match.
/// The MD tag describes the reference sequence for positions that differ from the read.
fn build_md_tag(
    transcript: &Transcript,
    read_seq: &[u8],
    genome: &Genome,
    is_reverse: bool,
) -> String {
    use cigar::op::Kind;

    // Build the SAM-order sequence (RC for reverse strand)
    let sam_seq: Vec<u8> = if is_reverse {
        read_seq.iter().rev().map(|&b| complement_base(b)).collect()
    } else {
        read_seq.to_vec()
    };

    let mut md = String::new();
    let mut match_count: u32 = 0;
    let mut genome_pos = transcript.genome_start;
    let mut read_pos: usize = 0;

    for op in &transcript.cigar {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    let ref_base = genome.get_base(genome_pos).unwrap_or(4);
                    let read_base = if read_pos < sam_seq.len() {
                        sam_seq[read_pos]
                    } else {
                        4
                    };
                    if ref_base == read_base {
                        match_count += 1;
                    } else {
                        write!(md, "{match_count}").unwrap();
                        match_count = 0;
                        md.push(decode_base(ref_base) as char);
                    }
                    genome_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Deletion => {
                write!(md, "{match_count}").unwrap();
                match_count = 0;
                md.push('^');
                for _ in 0..op.len() {
                    let ref_base = genome.get_base(genome_pos).unwrap_or(4);
                    md.push(decode_base(ref_base) as char);
                    genome_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += op.len();
            }
            Kind::Skip => {
                genome_pos += op.len() as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    // Emit trailing match count
    write!(md, "{match_count}").unwrap();
    md
}

/// Stamp paired-end mate-bookkeeping (PROPERLY_SEGMENTED, MATE_REVERSE_COMPLEMENTED,
/// RNEXT, PNEXT, TLEN) on two transcriptome-space mate records that share a
/// transcript by construction. The records and transcripts are paired as
/// (mate1_record, mate2_record, mate1_transcript, mate2_transcript).
pub(crate) fn apply_pe_transcriptome_mate_fields(
    rec1: &mut RecordBuf,
    rec2: &mut RecordBuf,
    t1: &Transcript,
    t2: &Transcript,
) -> Result<(), Error> {
    use crate::align::read_align::calculate_insert_size;
    use sam::alignment::record::Flags;

    *rec1.flags_mut() |= Flags::PROPERLY_SEGMENTED;
    *rec2.flags_mut() |= Flags::PROPERLY_SEGMENTED;

    if rec2.flags().is_reverse_complemented() {
        *rec1.flags_mut() |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    if rec1.flags().is_reverse_complemented() {
        *rec2.flags_mut() |= Flags::MATE_REVERSE_COMPLEMENTED;
    }

    *rec1.mate_reference_sequence_id_mut() = Some(t2.chr_idx);
    *rec2.mate_reference_sequence_id_mut() = Some(t1.chr_idx);

    let pos1 = (t1.genome_start + 1) as usize;
    let pos2 = (t2.genome_start + 1) as usize;
    *rec1.mate_alignment_start_mut() = Some(
        pos2.try_into()
            .map_err(|e| Error::Alignment(format!("invalid mate position {pos2}: {e}")))?,
    );
    *rec2.mate_alignment_start_mut() = Some(
        pos1.try_into()
            .map_err(|e| Error::Alignment(format!("invalid mate position {pos1}: {e}")))?,
    );

    let tlen = calculate_insert_size(t1, t2);
    *rec1.template_length_mut() = tlen;
    *rec2.template_length_mut() = -tlen;

    Ok(())
}

/// Build a SAM record for one mate of a paired-end read
#[allow(clippy::too_many_arguments)]
fn build_paired_mate_record(
    read_name: &str,
    mate_seq: &[u8],
    mate_qual: &[u8],
    clip5p: usize,
    clip3p: usize,
    transcript: &Transcript,
    mate_transcript: &Transcript,
    genome: &Genome,
    mapq: u8,
    is_first_mate: bool,
    is_proper_pair: bool,
    insert_size: i32,
    n_alignments: usize,
    hit_index: usize,
    ih_start: u32,
    combined_score: i32,
    attrs: SamAttributes,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();

    // Name
    record.name_mut().replace(read_name.into());

    // FLAGS
    let mut flags = sam::alignment::record::Flags::SEGMENTED; // 0x1 (paired)

    if is_proper_pair {
        flags |= sam::alignment::record::Flags::PROPERLY_SEGMENTED; // 0x2
    }

    if transcript.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED; // 0x10
    }

    // Mate reverse flag from the actual mate's alignment strand
    if mate_transcript.is_reverse {
        flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED; // 0x20
    }

    if is_first_mate {
        flags |= sam::alignment::record::Flags::FIRST_SEGMENT; // 0x40
    } else {
        flags |= sam::alignment::record::Flags::LAST_SEGMENT; // 0x80
    }

    if hit_index > 1 {
        flags |= sam::alignment::record::Flags::SECONDARY; // 0x100
    }

    *record.flags_mut() = flags;

    // RNAME (reference sequence name)
    if transcript.chr_idx >= genome.n_chr_real {
        return Err(Error::Alignment(format!(
            "invalid chromosome index {} (max {})",
            transcript.chr_idx,
            genome.n_chr_real - 1
        )));
    }
    *record.reference_sequence_id_mut() = Some(transcript.chr_idx);

    // POS (1-based, per-chromosome coordinate)
    // transcript.genome_start is a global genome coordinate, need to convert to per-chr
    let chr_start = genome.chr_start[transcript.chr_idx];
    let pos = (transcript.genome_start - chr_start + 1) as usize;
    *record.alignment_start_mut() = Some(
        pos.try_into()
            .map_err(|e| Error::Alignment(format!("invalid alignment position {pos}: {e}")))?,
    );

    // MAPQ
    *record.mapping_quality_mut() = MappingQuality::new(mapq);

    // CIGAR
    *record.cigar_mut() = transcript.cigar.iter().copied().collect();

    // RNEXT (mate reference sequence from the mate's actual alignment)
    *record.mate_reference_sequence_id_mut() = Some(mate_transcript.chr_idx);

    // PNEXT (mate position from the mate's actual alignment, per-chromosome coords)
    let mate_chr_start = genome.chr_start[mate_transcript.chr_idx];
    let mate_pos = (mate_transcript.genome_start - mate_chr_start + 1) as usize;
    *record.mate_alignment_start_mut() = Some(
        mate_pos
            .try_into()
            .map_err(|e| Error::Alignment(format!("invalid mate position {mate_pos}: {e}")))?,
    );

    // TLEN (insert size)
    *record.template_length_mut() = insert_size;

    // Core SEQ/QUAL/MD are built from the aligned (clipped) slice; apply_read_clips
    // below restores the full mate read + soft-clips (clip==0 => aligned == mate_seq).
    let aligned_seq = &mate_seq[clip5p..mate_seq.len() - clip3p];
    let aligned_qual = &mate_qual[clip5p..mate_qual.len() - clip3p];

    // Sequence and quality scores (reverse complement for reverse strand)
    if transcript.is_reverse {
        let seq_bytes: Vec<u8> = aligned_seq
            .iter()
            .rev()
            .map(|&b| decode_base(complement_base(b)))
            .collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);

        let mut qual = fastq_qual_to_phred(aligned_qual);
        qual.reverse();
        *record.quality_scores_mut() = QualityScores::from(qual);
    } else {
        let seq_bytes: Vec<u8> = aligned_seq.iter().map(|&b| decode_base(b)).collect();
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = QualityScores::from(fastq_qual_to_phred(aligned_qual));
    }

    // Optional tags: gated by --outSAMattributes
    let data = record.data_mut();
    if attrs.contains(SamAttributes::NH) {
        data.insert(Tag::ALIGNMENT_HIT_COUNT, Value::from(n_alignments as i32));
    }
    if attrs.contains(SamAttributes::HI) {
        data.insert(
            Tag::HIT_INDEX,
            Value::from(hit_index as i32 - 1 + ih_start as i32),
        );
    }
    if attrs.contains(SamAttributes::AS) {
        // STAR reports combined score (sum of both mates) for PE AS tag
        data.insert(Tag::ALIGNMENT_SCORE, Value::from(combined_score));
    }
    // nM (mismatch count) before NM (edit distance), matching STAR's tag order.
    if attrs.contains(SamAttributes::NMM) {
        data.insert(
            Tag::new(b'n', b'M'),
            Value::from(transcript.n_mismatch as i32),
        );
    }
    if attrs.contains(SamAttributes::NM) {
        data.insert(
            Tag::EDIT_DISTANCE,
            Value::from(sam_spec_nm(transcript.n_mismatch, &transcript.cigar)),
        );
    }
    if attrs.contains(SamAttributes::XS)
        && let Some(xs_strand) = derive_xs_strand(transcript)
    {
        data.insert(Tag::new(b'X', b'S'), Value::Character(xs_strand as u8));
    }
    if attrs.contains(SamAttributes::JM)
        && let Some(jm) = build_jm_tag(transcript)
    {
        data.insert(Tag::new(b'j', b'M'), jm);
    }
    if attrs.contains(SamAttributes::JI)
        && let Some(ji) = build_ji_tag(transcript, chr_start)
    {
        data.insert(Tag::new(b'j', b'I'), ji);
    }
    if attrs.contains(SamAttributes::MD) {
        let md = build_md_tag(transcript, aligned_seq, genome, transcript.is_reverse);
        data.insert(Tag::new(b'M', b'D'), Value::String(BString::from(md)));
    }

    if clip5p > 0 || clip3p > 0 {
        apply_read_clips(
            &mut record,
            &transcript.cigar,
            mate_seq,
            mate_qual,
            clip5p,
            clip3p,
            transcript.is_reverse,
        );
    }

    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::score::SpliceMotif;
    use crate::genome::Genome;
    use noodles::sam::alignment::record::cigar;
    use tempfile::NamedTempFile;

    fn make_test_genome() -> Genome {
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

    /// The solo barcode tags land on every record of the read, and only the
    /// requested ones are written. Values with no sequence behind them (a
    /// barcode read too short to hold a CB+UMI) are left off entirely.
    #[test]
    fn solo_barcode_tags_are_added_per_requested_attribute() {
        let values = SoloBarcodeTagValues {
            cb_seq: "ACGTACGTACGTACGT",
            cb_qual: "IIIIIIIIIIIIIIII",
            umi_seq: "ACGTACGTAC",
            umi_qual: "JJJJJJJJJJ",
            cb_corrected: Some("ACGTACGTACGTACGA"),
            cb_match: 1,
            barcode_seq: "",
            barcode_qual: "",
        };
        let mut records = vec![RecordBuf::default(), RecordBuf::default()];
        let attrs = SamAttributes::CR
            | SamAttributes::CY
            | SamAttributes::UR
            | SamAttributes::SM
            | SamAttributes::SS
            | SamAttributes::CB;
        add_solo_barcode_tags(&mut records, &values, attrs);

        for rec in &records {
            let data = rec.data();
            let get = |t: [u8; 2]| data.get(&Tag::new(t[0], t[1])).cloned();
            assert_eq!(
                get(*b"CR"),
                Some(Value::String(BString::from("ACGTACGTACGTACGT")))
            );
            assert_eq!(
                get(*b"CY"),
                Some(Value::String(BString::from("IIIIIIIIIIIIIIII")))
            );
            assert_eq!(
                get(*b"UR"),
                Some(Value::String(BString::from("ACGTACGTAC")))
            );
            assert_eq!(get(*b"sM"), Some(Value::Int32(1)));
            assert_eq!(
                get(*b"CB"),
                Some(Value::String(BString::from("ACGTACGTACGTACGA")))
            );
            // Not requested, and no value: absent.
            assert_eq!(get(*b"UY"), None);
            assert_eq!(get(*b"sS"), None);
            assert_eq!(get(*b"sQ"), None);
        }

        // No solo attributes requested → untouched records.
        let mut untouched = vec![RecordBuf::default()];
        add_solo_barcode_tags(&mut untouched, &values, SamAttributes::STANDARD);
        assert!(untouched[0].data().is_empty());
    }

    #[test]
    fn test_build_sam_header() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"]);

        let header = build_sam_header(&genome, &params).unwrap();

        // Check that we have reference sequences
        assert_eq!(header.reference_sequences().len(), 1);

        // Check that we have a program line (just check header is valid)
        assert_eq!(header.reference_sequences().len(), 1);
    }

    #[test]
    fn test_build_sam_header_pg_line_populated() {
        let genome = make_test_genome();
        let mut params = Parameters::parse_from(vec!["rustar-aligner", "--readFilesIn", "test.fq"]);
        params.command_line =
            Some("rustar-aligner --readFilesIn test.fq --runThreadN 4".to_string());

        let header = build_sam_header(&genome, &params).unwrap();
        let programs = header.programs().as_ref();
        let pg = programs
            .get(&b"rustar-aligner"[..])
            .expect("@PG line with ID:rustar-aligner must be present");

        let pn: &[u8] = pg
            .other_fields()
            .get(&program_tag::NAME)
            .expect("PN field must be present")
            .as_ref();
        assert_eq!(pn, b"rustar-aligner");

        let vn: &[u8] = pg
            .other_fields()
            .get(&program_tag::VERSION)
            .expect("VN field must be present")
            .as_ref();
        assert_eq!(vn, env!("CARGO_PKG_VERSION").as_bytes());

        let cl: &[u8] = pg
            .other_fields()
            .get(&program_tag::COMMAND_LINE)
            .expect("CL field must be present")
            .as_ref();
        assert!(!cl.is_empty(), "CL field must be non-empty");
        assert_eq!(cl, b"rustar-aligner --readFilesIn test.fq --runThreadN 4");
    }

    #[test]
    fn test_build_sam_header_pg_line_cl_quoted() {
        let genome = make_test_genome();
        let params =
            Parameters::parse_from(vec!["rustar-aligner", "--readFilesIn", "some words.fq"]);

        let header = build_sam_header(&genome, &params).unwrap();
        let programs = header.programs().as_ref();
        let pg = programs.get(&b"rustar-aligner"[..]).unwrap();
        let cl = pg
            .other_fields()
            .get(&program_tag::COMMAND_LINE)
            .expect("CL field must be present");
        assert_eq!(cl, "rustar-aligner --readFilesIn 'some words.fq'");
    }

    #[test]
    fn test_build_sam_header_pg_line_default_cl_when_unset() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec!["rustar-aligner", "--readFilesIn", "test.fq\0"]);
        // when there are null bytes in the args, `command_line` is set to None
        assert!(params.command_line.is_none());

        let header = build_sam_header(&genome, &params).unwrap();
        let programs = header.programs().as_ref();
        let pg = programs.get(&b"rustar-aligner"[..]).unwrap();
        let cl: &[u8] = pg
            .other_fields()
            .get(&program_tag::COMMAND_LINE)
            .expect("CL field must be present even when command_line is None")
            .as_ref();
        assert!(!cl.is_empty());
    }

    #[test]
    fn test_build_sam_header_with_rg() {
        let genome = make_test_genome();
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMattrRGline",
            "ID:rg0",
            "SM:sample0",
            "LB:lib0",
        ]);
        let header = build_sam_header(&genome, &params).unwrap();
        let rgs = header.read_groups();
        assert_eq!(rgs.len(), 1);
        assert!(rgs.contains_key(&b"rg0"[..]));
        let map = rgs.get(&b"rg0"[..]).unwrap();
        // SM and LB should be present as other_fields
        let sm_tag = HeaderOtherTag::<_>::try_from(*b"SM").unwrap();
        let lb_tag = HeaderOtherTag::<_>::try_from(*b"LB").unwrap();
        let sm: &[u8] = map.other_fields().get(&sm_tag).unwrap().as_ref();
        let lb: &[u8] = map.other_fields().get(&lb_tag).unwrap().as_ref();
        assert_eq!(sm, b"sample0");
        assert_eq!(lb, b"lib0");
    }

    #[test]
    fn test_sam_output_includes_rg_header_and_tag() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};
        use std::io::Read;

        let genome = make_test_genome();
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMattrRGline",
            "ID:rg0",
            "SM:sample0",
        ]);

        let tmpfile = NamedTempFile::new().unwrap();
        let mut writer = SamWriter::create(tmpfile.path(), &genome, &params).unwrap();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        writer
            .write_alignment(
                "read1",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                &[transcript],
                &genome,
                &params,
                1,
            )
            .unwrap();
        drop(writer);

        let mut contents = String::new();
        std::fs::File::open(tmpfile.path())
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        assert!(
            contents.contains("@RG\tID:rg0\tSM:sample0"),
            "missing @RG header; got:\n{contents}"
        );
        assert!(
            contents.contains("RG:Z:rg0"),
            "missing RG:Z tag on record; got:\n{contents}"
        );
    }

    #[test]
    fn test_sam_writer_creation() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"]);

        let tmpfile = NamedTempFile::new().unwrap();
        let writer = SamWriter::create(tmpfile.path(), &genome, &params);
        assert!(writer.is_ok());
    }

    /// Regression test for issue #34: BAM binary QUAL must store raw Phred values
    /// (0-93), not FASTQ ASCII bytes (Phred+33). Each FASTQ ASCII byte should be
    /// reduced by 33 before being placed in the BAM QUAL field.
    #[test]
    fn test_qual_phred33_offset_stripped_for_bam() {
        // FASTQ ASCII 'I' = 73 → Phred 40
        let read_qual = b"III";
        let read_seq = vec![0u8, 1, 2]; // ACG
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "t.fq"]);
        let record = SamWriter::build_unmapped_record(
            "read1",
            &read_seq,
            read_qual,
            &params,
            UnmappedReason::Other,
        )
        .unwrap();

        let stored: &[u8] = record.quality_scores().as_ref();
        assert_eq!(stored, &[40u8, 40, 40]);
    }

    /// `saturating_sub` must clamp malformed FASTQ bytes (< 33) to 0 rather
    /// than underflowing.
    #[test]
    fn test_qual_phred33_offset_saturating_sub() {
        // Byte 32 is below the Phred+33 floor; should clamp to 0.
        let read_qual = &[32u8, 33, 34][..];
        let read_seq = vec![0u8, 1, 2];
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "t.fq"]);
        let record = SamWriter::build_unmapped_record(
            "read1",
            &read_seq,
            read_qual,
            &params,
            UnmappedReason::Other,
        )
        .unwrap();

        let stored: &[u8] = record.quality_scores().as_ref();
        assert_eq!(stored, &[0u8, 0, 1]);
    }

    #[test]
    fn test_unmapped_record_tags_emitted() {
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "t.fq"]);
        let read_seq = vec![0u8, 1, 2, 3];
        let read_qual = b"IIII";

        for (reason, expected_ut) in [
            (UnmappedReason::Other, b'0'),
            (UnmappedReason::TooShort, b'1'),
            (UnmappedReason::TooManyMismatches, b'2'),
            (UnmappedReason::TooManyLoci, b'3'),
        ] {
            let record =
                SamWriter::build_unmapped_record("r", &read_seq, read_qual, &params, reason)
                    .unwrap();
            let data = record.data();
            // NH/HI/AS present (Standard attrs enabled by default)
            assert_eq!(
                data.get(&Tag::ALIGNMENT_HIT_COUNT),
                Some(&Value::from(0i32)),
                "NH:i:0 missing for {reason:?}"
            );
            assert_eq!(
                data.get(&Tag::HIT_INDEX),
                Some(&Value::from(0i32)),
                "HI:i:0 missing for {reason:?}"
            );
            assert_eq!(
                data.get(&Tag::ALIGNMENT_SCORE),
                Some(&Value::from(0i32)),
                "AS:i:0 missing for {reason:?}"
            );
            // uT:A: always emitted with correct value
            assert_eq!(
                data.get(&Tag::new(b'u', b'T')),
                Some(&Value::Character(expected_ut)),
                "uT:A: wrong for {reason:?}"
            );
        }
    }

    #[test]
    fn test_transcript_to_record() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 10,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 2,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::STANDARD,
        );
        assert!(record.is_ok());

        let record = record.unwrap();
        assert_eq!(
            record.name().map(ToString::to_string),
            Some("read1".to_string())
        );
        assert_eq!(record.reference_sequence_id(), Some(0));
        assert_eq!(record.alignment_start().map(usize::from), Some(11)); // 1-based
        // hit_index=1, so NOT secondary
        assert!(!record.flags().is_secondary());
    }

    #[test]
    fn test_build_transcriptome_records_stamps_rg_tag() {
        let params = Parameters::parse_from(vec![
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMattrRGline",
            "ID:rg0",
            "SM:sample0",
        ]);

        let projected = vec![
            Transcript {
                chr_idx: 0,
                genome_start: 0,
                genome_end: 4,
                is_reverse: false,
                exons: vec![],
                cigar: vec![cigar::Op::new(cigar::op::Kind::Match, 4)],
                score: 100,
                n_mismatch: 0,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
            Transcript {
                chr_idx: 0,
                genome_start: 2,
                genome_end: 6,
                is_reverse: true,
                exons: vec![],
                cigar: vec![cigar::Op::new(cigar::op::Kind::Match, 4)],
                score: 90,
                n_mismatch: 1,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
        ];

        let records = SamWriter::build_transcriptome_records(
            "read1",
            &[0, 1, 2, 3],
            &[30, 30, 30, 30],
            &projected,
            255,
            &params,
            0,
        )
        .expect("build_transcriptome_records");

        assert_eq!(records.len(), 2, "expected one record per projected hit");
        for (i, rec) in records.iter().enumerate() {
            let rg = rec
                .data()
                .get(&Tag::READ_GROUP)
                .unwrap_or_else(|| panic!("record {i} missing RG tag"));
            match rg {
                Value::String(s) => assert_eq!(
                    s.as_slice(),
                    b"rg0",
                    "record {i}: expected RG:Z:rg0, got {s:?}"
                ),
                other => panic!("record {i}: RG tag is not a string: {other:?}"),
            }
        }
    }

    #[test]
    fn test_build_transcriptome_records_no_rg_when_unset() {
        let params = Parameters::parse_from(vec!["rustar-aligner", "--readFilesIn", "test.fq"]);

        let projected = vec![Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![cigar::Op::new(cigar::op::Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        }];

        let records = SamWriter::build_transcriptome_records(
            "read1",
            &[0, 1, 2, 3],
            &[30, 30, 30, 30],
            &projected,
            255,
            &params,
            0,
        )
        .expect("build_transcriptome_records");

        assert_eq!(records.len(), 1);
        assert!(
            records[0].data().get(&Tag::READ_GROUP).is_none(),
            "RG tag should not be present when --outSAMattrRGline is unset"
        );
    }

    #[test]
    fn test_build_paired_unmapped_records() {
        let mate1_seq = vec![0, 1, 2, 3]; // ACGT
        let mate1_qual = vec![30, 30, 30, 30];
        let mate2_seq = vec![3, 2, 1, 0]; // TGCA
        let mate2_qual = vec![30, 30, 30, 30];
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "t.fq"]);

        let records = SamWriter::build_paired_unmapped_records(
            "read1",
            &mate1_seq,
            &mate1_qual,
            &mate2_seq,
            &mate2_qual,
            &params,
            UnmappedReason::Other,
        )
        .unwrap();

        assert_eq!(records.len(), 2);

        // Check mate1 record
        let rec1 = &records[0];
        assert_eq!(
            rec1.name().map(ToString::to_string),
            Some("read1".to_string())
        );
        assert!(rec1.flags().is_segmented());
        assert!(rec1.flags().is_unmapped());
        assert!(rec1.flags().is_mate_unmapped());
        assert!(rec1.flags().is_first_segment());

        // Check mate2 record
        let rec2 = &records[1];
        assert_eq!(
            rec2.name().map(ToString::to_string),
            Some("read1".to_string())
        );
        assert!(rec2.flags().is_segmented());
        assert!(rec2.flags().is_unmapped());
        assert!(rec2.flags().is_mate_unmapped());
        assert!(rec2.flags().is_last_segment());
    }

    #[test]
    fn test_build_paired_mate_record_flags() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();

        let mate1_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate2_transcript = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 3)],
            score: 90,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate_seq = vec![0, 1, 2, 3];
        let mate_qual = vec![30, 30, 30, 30];

        // Test first mate (forward), mate2 is reverse → 0x20 should be set
        let rec1 = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            0,
            0,
            &mate1_transcript,
            &mate2_transcript,
            &genome,
            255,
            true, // is_first_mate
            true, // is_proper_pair
            300,
            1,   // n_alignments
            1,   // hit_index
            1,   // ih_start
            190, // combined_score (100+90)
            SamAttributes::STANDARD,
        )
        .unwrap();

        assert!(rec1.flags().is_segmented());
        assert!(rec1.flags().is_properly_segmented());
        assert!(rec1.flags().is_first_segment());
        assert!(!rec1.flags().is_last_segment());
        assert!(!rec1.flags().is_reverse_complemented()); // mate1 is forward
        assert!(rec1.flags().is_mate_reverse_complemented()); // mate2 is reverse
        assert!(!rec1.flags().is_secondary());
        assert_eq!(rec1.template_length(), 300);

        // Test second mate (reverse), mate1 is forward → 0x20 should NOT be set
        let rec2 = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            0,
            0,
            &mate2_transcript,
            &mate1_transcript,
            &genome,
            255,
            false, // is_first_mate
            true,  // is_proper_pair
            -300,
            1,   // n_alignments
            1,   // hit_index
            1,   // ih_start
            190, // combined_score (100+90)
            SamAttributes::STANDARD,
        )
        .unwrap();

        assert!(rec2.flags().is_segmented());
        assert!(rec2.flags().is_properly_segmented());
        assert!(!rec2.flags().is_first_segment());
        assert!(rec2.flags().is_last_segment());
        assert!(rec2.flags().is_reverse_complemented()); // mate2 is reverse
        assert!(!rec2.flags().is_mate_reverse_complemented()); // mate1 is forward
        assert!(!rec2.flags().is_secondary());
        assert_eq!(rec2.template_length(), -300);
    }

    #[test]
    fn test_build_paired_mate_record_mate_fields() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();

        // Mate1 at position 0 (chr_start=0, so per-chr pos = 1)
        let this_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 200,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Mate2 at position 4 (chr_start=0, so per-chr pos = 5)
        let mate_transcript = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 3)],
            score: 150,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate_seq = vec![0; 4];
        let mate_qual = vec![30; 4];

        let rec = build_paired_mate_record(
            "read1",
            &mate_seq,
            &mate_qual,
            0,
            0,
            &this_transcript,
            &mate_transcript,
            &genome,
            60,
            true,
            true,
            250,
            1,   // n_alignments
            1,   // hit_index
            1,   // ih_start
            350, // combined_score (200+150)
            SamAttributes::STANDARD,
        )
        .unwrap();

        // RNEXT = mate's chr_idx
        assert_eq!(rec.mate_reference_sequence_id(), Some(0));

        // PNEXT = mate's per-chr position (genome_start=4, chr_start=0 → pos=5)
        assert_eq!(rec.mate_alignment_start().map(usize::from), Some(5));

        // Check TLEN
        assert_eq!(rec.template_length(), 250);

        // AS is the combined score (STAR behavior); nM is per-mate mismatches
        let data = rec.data();
        assert_eq!(
            data.get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(350_i32)),
            "AS should be combined score (200+150=350)"
        );
        assert_eq!(
            data.get(&Tag::new(b'n', b'M')),
            Some(&Value::from(0_i32)),
            "nM should be 0 (no mismatches in this mate)"
        );
    }

    #[test]
    fn test_apply_pe_transcriptome_mate_fields() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};
        use sam::alignment::record::Flags;

        let chr_idx = 7;

        let t1 = Transcript {
            chr_idx,
            genome_start: 100,
            genome_end: 130,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 100,
                genome_end: 130,
                read_start: 0,
                read_end: 30,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 30)],
            score: 60,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let t2 = Transcript {
            chr_idx,
            genome_start: 200,
            genome_end: 230,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 200,
                genome_end: 230,
                read_start: 0,
                read_end: 30,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 30)],
            score: 60,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mut rec1 = RecordBuf::default();
        *rec1.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let mut rec2 = RecordBuf::default();
        *rec2.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED;

        apply_pe_transcriptome_mate_fields(&mut rec1, &mut rec2, &t1, &t2).unwrap();

        assert!(rec1.flags().is_properly_segmented());
        assert!(rec2.flags().is_properly_segmented());

        assert!(rec1.flags().is_mate_reverse_complemented());
        assert!(!rec2.flags().is_mate_reverse_complemented());

        assert_eq!(rec1.mate_reference_sequence_id(), Some(chr_idx));
        assert_eq!(rec2.mate_reference_sequence_id(), Some(chr_idx));

        let pnext1 = rec1.mate_alignment_start().unwrap();
        let pnext2 = rec2.mate_alignment_start().unwrap();
        assert_eq!(usize::from(pnext1), 201);
        assert_eq!(usize::from(pnext2), 101);

        assert_ne!(rec1.template_length(), 0);
        assert_ne!(rec2.template_length(), 0);
        assert!(rec1.template_length() > 0);
        assert!(rec2.template_length() < 0);
        assert_eq!(rec1.template_length(), -rec2.template_length());
    }

    #[test]
    fn test_apply_pe_transcriptome_mate_fields_reverse_cluster() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};
        use sam::alignment::record::Flags;

        let chr_idx = 3;

        let t1 = Transcript {
            chr_idx,
            genome_start: 200,
            genome_end: 230,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 200,
                genome_end: 230,
                read_start: 0,
                read_end: 30,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 30)],
            score: 60,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let t2 = Transcript {
            chr_idx,
            genome_start: 100,
            genome_end: 130,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 100,
                genome_end: 130,
                read_start: 0,
                read_end: 30,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 30)],
            score: 60,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mut rec1 = RecordBuf::default();
        *rec1.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::REVERSE_COMPLEMENTED;
        let mut rec2 = RecordBuf::default();
        *rec2.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT;

        apply_pe_transcriptome_mate_fields(&mut rec1, &mut rec2, &t1, &t2).unwrap();

        assert!(rec1.template_length() < 0);
        assert!(rec2.template_length() > 0);
        assert_eq!(rec1.template_length(), -rec2.template_length());
    }

    #[test]
    fn test_tags_nh_hi_as_nm() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        // Transcript with 2 mismatches and a 3bp deletion → NM = 2 + 3 = 5
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 20),
                Op::new(Kind::Deletion, 3),
                Op::new(Kind::Match, 30),
            ],
            score: 100,
            n_mismatch: 2,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            3, // n_alignments
            2, // hit_index
            1, // ih_start
            // Standard + NM so both nM (mismatch count) and NM (edit distance) emit.
            SamAttributes::STANDARD | SamAttributes::NM,
        )
        .unwrap();

        // hit_index=2 → secondary
        assert!(record.flags().is_secondary());

        let data = record.data();
        assert_eq!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT),
            Some(&Value::from(3_i32)),
            "NH tag should be 3"
        );
        assert_eq!(
            data.get(&Tag::HIT_INDEX),
            Some(&Value::from(2_i32)),
            "HI tag should be 2"
        );
        assert_eq!(
            data.get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(100_i32)),
            "AS tag should be 100"
        );
        assert_eq!(
            data.get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(5_i32)),
            "NM should be edit distance (2 mismatches + 3 deleted bases)"
        );
        assert_eq!(
            data.get(&Tag::new(b'n', b'M')),
            Some(&Value::from(2_i32)),
            "nM should be mismatches only (2)"
        );
    }

    #[test]
    fn test_edit_distance_computation() {
        use cigar::op::{Kind, Op};
        // Pure match: NM = n_mismatch only
        assert_eq!(sam_spec_nm(3, &[Op::new(Kind::Match, 50)]), 3);

        // Match + Ins + Del: NM = n_mismatch + ins_len + del_len
        assert_eq!(
            sam_spec_nm(
                1,
                &[
                    Op::new(Kind::Match, 20),
                    Op::new(Kind::Insertion, 5),
                    Op::new(Kind::Match, 10),
                    Op::new(Kind::Deletion, 7),
                    Op::new(Kind::Match, 20),
                ],
            ),
            13
        );

        // RefSkip (splice junction) should NOT count toward NM
        assert_eq!(
            sam_spec_nm(
                0,
                &[
                    Op::new(Kind::Match, 25),
                    Op::new(Kind::Skip, 1000),
                    Op::new(Kind::Match, 25),
                ],
            ),
            0
        );

        // Soft clips should NOT count toward NM
        assert_eq!(
            sam_spec_nm(
                2,
                &[
                    Op::new(Kind::SoftClip, 10),
                    Op::new(Kind::Match, 40),
                    Op::new(Kind::SoftClip, 10),
                ],
            ),
            2
        );
    }

    #[test]
    fn test_transcript_to_record_has_tags() {
        use cigar::op::{Kind, Op};
        // Verify the existing test_transcript_to_record scenario also has tags
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 10,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 2,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1, // unique mapper
            1,
            1, // ih_start
            SamAttributes::STANDARD,
        )
        .unwrap();

        let data = record.data();
        // Unique mapper: NH=1, HI=1
        assert_eq!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT),
            Some(&Value::from(1_i32))
        );
        assert_eq!(data.get(&Tag::HIT_INDEX), Some(&Value::from(1_i32)));
        assert_eq!(data.get(&Tag::ALIGNMENT_SCORE), Some(&Value::from(100_i32)));
        // Standard = NH HI AS nM: the mismatch count nM, not the edit-distance NM.
        assert_eq!(data.get(&Tag::EDIT_DISTANCE), None);
        assert_eq!(data.get(&Tag::new(b'n', b'M')), Some(&Value::from(2_i32)));
        assert_eq!(data.get(&Tag::new(b'X', b'S')), None);
    }

    #[test]
    fn test_secondary_flag() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"]);

        let transcripts = vec![
            Transcript {
                chr_idx: 0,
                genome_start: 0,
                genome_end: 50,
                is_reverse: false,
                exons: vec![],
                cigar: vec![Op::new(Kind::Match, 50)],
                score: 100,
                n_mismatch: 0,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
            Transcript {
                chr_idx: 0,
                genome_start: 2,
                genome_end: 52,
                is_reverse: false,
                exons: vec![],
                cigar: vec![Op::new(Kind::Match, 50)],
                score: 98,
                n_mismatch: 1,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
            Transcript {
                chr_idx: 0,
                genome_start: 4,
                genome_end: 54,
                is_reverse: true,
                exons: vec![],
                cigar: vec![Op::new(Kind::Match, 50)],
                score: 96,
                n_mismatch: 2,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
        ];

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 3);

        // Record 0 (HI=1): NOT secondary
        assert!(!records[0].flags().is_secondary());
        // Record 1 (HI=2): IS secondary
        assert!(records[1].flags().is_secondary());
        // Record 2 (HI=3): IS secondary + reverse complemented
        assert!(records[2].flags().is_secondary());
        assert!(records[2].flags().is_reverse_complemented());
    }

    #[test]
    fn test_out_sam_flag_or_and_applied_to_mapped_record() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();
        // OR in 0x200 (QC-fail); AND out 0x10 (reverse) via 65535 & !16 = 65519.
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMflagOR",
            "512",
            "--outSAMflagAND",
            "65519",
        ]);

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: true,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            std::slice::from_ref(&transcript),
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        // REVERSE_COMPLEMENTED (0x10) was set by the transcript but AND-ed out;
        // QC_FAIL (0x200) was OR-ed in.
        assert!(!records[0].flags().is_reverse_complemented());
        assert!(records[0].flags().is_qc_fail());
    }

    #[test]
    fn test_out_sam_primary_flag_all_best_score() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMprimaryFlag",
            "AllBestScore",
        ]);

        let mk = |genome_start: u64, score: i32| Transcript {
            chr_idx: 0,
            genome_start,
            genome_end: genome_start + 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };
        // Two alignments tied for the best score (100), one strictly worse (98).
        let transcripts = vec![mk(0, 100), mk(2, 98), mk(4, 100)];

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 3);
        assert!(!records[0].flags().is_secondary(), "best-score #1 primary");
        assert!(
            records[1].flags().is_secondary(),
            "worse-score stays secondary"
        );
        assert!(
            !records[2].flags().is_secondary(),
            "best-score #2 also primary"
        );
    }

    #[test]
    fn test_out_sam_attr_ih_start_shifts_hi_tag() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMattrIHstart",
            "0",
        ]);

        let transcripts = vec![
            Transcript {
                chr_idx: 0,
                genome_start: 0,
                genome_end: 50,
                is_reverse: false,
                exons: vec![],
                cigar: vec![Op::new(Kind::Match, 50)],
                score: 100,
                n_mismatch: 0,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
            Transcript {
                chr_idx: 0,
                genome_start: 2,
                genome_end: 52,
                is_reverse: false,
                exons: vec![],
                cigar: vec![Op::new(Kind::Match, 50)],
                score: 98,
                n_mismatch: 1,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            },
        ];

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 2);
        // --outSAMattrIHstart 0: HI starts at 0 instead of the STAR default 1.
        assert_eq!(
            records[0].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(0_i32))
        );
        assert_eq!(
            records[1].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(1_i32))
        );
        // The secondary FLAG bit is unaffected by IHstart (still rank-based, not tag-value-based).
        assert!(!records[0].flags().is_secondary());
        assert!(records[1].flags().is_secondary());
    }

    #[test]
    fn test_out_sam_tlen_mode2_contained_pair() {
        use cigar::op::{Kind, Op};
        // mate2 fully contained within mate1's span: [100,200) vs [120,180).
        // Mode 1 (default) just returns whatever the caller already computed (opaque here).
        // Mode 2 uses the per-mate max/min span (100) and signs by genomic leftmost (mate1).
        let mate1 = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };
        let mate2 = Transcript {
            genome_start: 120,
            genome_end: 180,
            ..mate1.clone()
        };

        let params_mode1 = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"]);
        assert_eq!(
            paired_tlen(&mate1, &mate2, true, 42, &params_mode1),
            42,
            "mode 1 (default) passes the caller-computed value through unchanged"
        );

        let params_mode2 = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMtlen",
            "2",
        ]);
        assert_eq!(
            paired_tlen(&mate1, &mate2, true, 42, &params_mode2),
            100,
            "mate1 is leftmost -> positive per-mate span"
        );
        assert_eq!(
            paired_tlen(&mate1, &mate2, false, -42, &params_mode2),
            -100,
            "mate2 is not leftmost -> negative per-mate span"
        );
    }

    #[test]
    fn test_xs_tag_spliced() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::ALL,
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            Some(&Value::Character(b'+')),
            "XS should be '+' for GT/AG motif"
        );
    }

    #[test]
    fn test_xs_tag_unspliced() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::ALL,
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should NOT be present for unspliced reads"
        );
    }

    #[test]
    fn test_xs_tag_reverse_strand() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::CtAc],
            junction_annotated: vec![false],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::ALL,
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            Some(&Value::Character(b'-')),
            "XS should be '-' for CT/AC motif"
        );
    }

    #[test]
    fn test_xs_tag_conflicting_motifs() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 300,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc], // +strand and -strand
            junction_annotated: vec![false, false],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::ALL,
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should NOT be present when junction motifs conflict on strand"
        );
    }

    #[test]
    fn test_xs_not_emitted_when_disabled() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1,                       // ih_start
            SamAttributes::STANDARD, // XS not in standard attrs
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should NOT be present when not in attribute set"
        );
    }

    #[test]
    fn test_out_sam_mult_nmax() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "test.fq",
            "--outSAMmultNmax",
            "3",
        ]);

        let transcripts: Vec<Transcript> = (0..5)
            .map(|i| Transcript {
                chr_idx: 0,
                genome_start: i as u64,
                genome_end: (i + 50) as u64,
                is_reverse: false,
                exons: vec![],
                cigar: vec![Op::new(Kind::Match, 50)],
                score: 100 - i,
                n_mismatch: 0,
                n_gap: 0,
                n_junction: 0,
                junction_motifs: vec![],
                junction_annotated: vec![],
            })
            .collect();

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        // Only 3 records output despite 5 transcripts
        assert_eq!(records.len(), 3);

        // NH should be 3 (number of reported alignments)
        for rec in &records {
            let data = rec.data();
            assert_eq!(
                data.get(&Tag::ALIGNMENT_HIT_COUNT),
                Some(&Value::from(3_i32)),
                "NH should be 3 (capped by outSAMmultNmax)"
            );
        }

        // HI should be 1, 2, 3
        assert_eq!(
            records[0].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(1_i32))
        );
        assert_eq!(
            records[1].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(2_i32))
        );
        assert_eq!(
            records[2].data().get(&Tag::HIT_INDEX),
            Some(&Value::from(3_i32))
        );

        // First is primary, rest are secondary
        assert!(!records[0].flags().is_secondary());
        assert!(records[1].flags().is_secondary());
        assert!(records[2].flags().is_secondary());
    }

    #[test]
    fn test_build_jm_tag_basic() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
        };

        let jm = build_jm_tag(&transcript);
        assert!(jm.is_some());
        // GT/AG = motif code 1, not annotated → 1
        assert_eq!(jm.unwrap(), Value::Array(Array::Int8(vec![1])));
    }

    #[test]
    fn test_build_jm_tag_annotated() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![true],
        };

        let jm = build_jm_tag(&transcript);
        assert!(jm.is_some());
        // GT/AG = motif code 1, annotated → 1 + 20 = 21
        assert_eq!(jm.unwrap(), Value::Array(Array::Int8(vec![21])));
    }

    #[test]
    fn test_build_jm_tag_empty() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        assert!(build_jm_tag(&transcript).is_none());
    }

    #[test]
    fn test_build_jm_tag_multiple_junctions() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 400,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc],
            junction_annotated: vec![true, false],
        };

        let jm = build_jm_tag(&transcript);
        assert!(jm.is_some());
        // GT/AG annotated=21, CT/AC not annotated=2
        assert_eq!(jm.unwrap(), Value::Array(Array::Int8(vec![21, 2])));
    }

    #[test]
    fn test_build_ji_tag_basic() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 325,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 200),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
        };

        // chr_start=0, genome_start=100, intron starts at 125, ends at 324
        let ji = build_ji_tag(&transcript, 0);
        assert!(ji.is_some());
        // Intron start: 100+25 - 0 + 1 = 126, Intron end: 100+25+200 - 0 = 325
        assert_eq!(ji.unwrap(), Value::Array(Array::Int32(vec![126, 325])));
    }

    #[test]
    fn test_build_ji_tag_empty() {
        use cigar::op::{Kind, Op};
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 50,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        assert!(build_ji_tag(&transcript, 0).is_none());
    }

    #[test]
    fn test_build_md_tag_perfect_match() {
        use cigar::op::{Kind, Op};
        // Genome: ACGTACGT (A=0,C=1,G=2,T=3)
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 4,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Read exactly matches genome[0..4] = ACGT
        let read_seq = vec![0, 1, 2, 3];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "4");
    }

    #[test]
    fn test_build_md_tag_mismatches() {
        use cigar::op::{Kind, Op};
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 2,
            n_mismatch: 2,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Position 1: read=A, ref=C → mismatch (C in MD)
        // Position 3: read=A, ref=T → mismatch (T in MD)
        let read_seq = vec![0, 0, 2, 0]; // AAGA
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "1C1T0");
    }

    #[test]
    fn test_build_md_tag_deletion() {
        use cigar::op::{Kind, Op};
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 6,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 2),
            ],
            score: 4,
            n_mismatch: 0,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Read: A,C,[del G,T],A,C
        let read_seq = vec![0, 1, 0, 1];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "2^GT2");
    }

    #[test]
    fn test_build_md_tag_insertion() {
        use cigar::op::{Kind, Op};
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::Match, 2),
            ],
            score: 4,
            n_mismatch: 0,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Insertions are invisible in MD — just match counts
        let read_seq = vec![0, 1, 3, 3, 2, 3];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "4");
    }

    #[test]
    fn test_build_md_tag_soft_clip() {
        use cigar::op::{Kind, Op};
        // Genome: ACGTACGT
        let genome = make_test_genome();
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 2, // Starts at G
            genome_end: 6,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::SoftClip, 2),
                Op::new(Kind::Match, 4),
                Op::new(Kind::SoftClip, 2),
            ],
            score: 4,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Soft clips don't appear in MD
        let read_seq = vec![0, 0, 2, 3, 0, 1, 0, 0];
        let md = build_md_tag(&transcript, &read_seq, &genome, false);
        assert_eq!(md, "4");
    }

    #[test]
    fn test_tags_jm_ji_md_in_record() {
        use cigar::op::{Kind, Op};
        // Verify tags appear in a full transcript_to_record call
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 4,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::ALL,
        )
        .unwrap();

        let data = record.data();
        // No junctions → no jM/jI tags
        assert!(data.get(&Tag::new(b'j', b'M')).is_none());
        assert!(data.get(&Tag::new(b'j', b'I')).is_none());
        // MD should be present when in attrs
        assert_eq!(
            data.get(&Tag::new(b'M', b'D')),
            Some(&Value::String(BString::from("4")))
        );
    }

    #[test]
    fn test_build_paired_mate_record_cross_strand() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();

        // Mate1: forward, chr 0, pos 0
        let mate1_trans = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Mate2: reverse, chr 0, pos 4
        let mate2_trans = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: true,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 3)],
            score: 90,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let seq = vec![0, 1, 2, 3];
        let qual = vec![30, 30, 30, 30];

        // Mate1 record: mate is reverse → 0x20 set, PNEXT=5
        let rec1 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            0,
            0,
            &mate1_trans,
            &mate2_trans,
            &genome,
            255,
            true,
            true,
            7,
            1,
            1,
            1,   // ih_start
            190, // combined_score (100+90)
            SamAttributes::STANDARD,
        )
        .unwrap();

        assert!(rec1.flags().is_mate_reverse_complemented());
        assert!(!rec1.flags().is_reverse_complemented());
        assert_eq!(rec1.mate_alignment_start().map(usize::from), Some(5)); // genome_start=4, chr_start=0 → 5

        // Mate2 record: mate is forward → 0x20 NOT set, PNEXT=1
        let rec2 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            0,
            0,
            &mate2_trans,
            &mate1_trans,
            &genome,
            255,
            false,
            true,
            -7,
            1,
            1,
            1,   // ih_start
            190, // combined_score (100+90)
            SamAttributes::STANDARD,
        )
        .unwrap();

        assert!(!rec2.flags().is_mate_reverse_complemented());
        assert!(rec2.flags().is_reverse_complemented());
        assert_eq!(rec2.mate_alignment_start().map(usize::from), Some(1)); // genome_start=0, chr_start=0 → 1
    }

    #[test]
    fn test_build_paired_mate_record_per_mate_tags() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();

        // Mate1: score=100, 0 mismatches, no junctions
        let mate1_trans = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Mate2: score=80, 2 mismatches, 1 deletion
        let mate2_trans = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
                i_frag: 0,
            }],
            cigar: vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, 1),
            ],
            score: 80,
            n_mismatch: 2,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let seq1 = vec![0, 1, 2, 3];
        let qual1 = vec![30, 30, 30, 30];
        let seq2 = vec![0, 1, 2];
        let qual2 = vec![30, 30, 30];

        // Both mates should have combined AS (STAR behavior)
        let rec1 = build_paired_mate_record(
            "read1",
            &seq1,
            &qual1,
            0,
            0,
            &mate1_trans,
            &mate2_trans,
            &genome,
            255,
            true,
            true,
            7,
            1,
            1,
            1,   // ih_start
            180, // combined_score (100+80)
            SamAttributes::STANDARD,
        )
        .unwrap();
        assert_eq!(
            rec1.data().get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(180_i32)),
            "Mate1 AS should be combined score (100+80=180)"
        );
        // NM attribute maps to nM tag (mismatches only)
        assert_eq!(
            rec1.data().get(&Tag::new(b'n', b'M')),
            Some(&Value::from(0_i32)),
            "Mate1 nM should be 0"
        );

        // Mate2 also gets combined AS
        let rec2 = build_paired_mate_record(
            "read1",
            &seq2,
            &qual2,
            0,
            0,
            &mate2_trans,
            &mate1_trans,
            &genome,
            255,
            false,
            true,
            -7,
            1,
            1,
            1,   // ih_start
            180, // combined_score (100+80)
            SamAttributes::STANDARD,
        )
        .unwrap();
        assert_eq!(
            rec2.data().get(&Tag::ALIGNMENT_SCORE),
            Some(&Value::from(180_i32)),
            "Mate2 AS should be combined score (100+80=180)"
        );
        // nM = mismatches only (no indel contribution)
        assert_eq!(
            rec2.data().get(&Tag::new(b'n', b'M')),
            Some(&Value::from(2_i32)),
            "Mate2 nM should be 2 (mismatches only, not edit distance)"
        );
    }

    #[test]
    fn test_build_paired_mate_record_both_forward() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();

        // Both mates forward
        let mate1_trans = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate2_trans = Transcript {
            chr_idx: 0,
            genome_start: 4,
            genome_end: 7,
            is_reverse: false, // Also forward
            exons: vec![Exon {
                genome_start: 4,
                genome_end: 7,
                read_start: 0,
                read_end: 3,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 3)],
            score: 90,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let seq = vec![0, 1, 2, 3];
        let qual = vec![30, 30, 30, 30];

        // Both forward → neither should have 0x20 set
        let rec1 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            0,
            0,
            &mate1_trans,
            &mate2_trans,
            &genome,
            255,
            true,
            true,
            7,
            1,
            1,
            1,   // ih_start
            190, // combined_score (100+90)
            SamAttributes::STANDARD,
        )
        .unwrap();
        assert!(!rec1.flags().is_mate_reverse_complemented());
        assert!(!rec1.flags().is_reverse_complemented());

        let rec2 = build_paired_mate_record(
            "read1",
            &seq,
            &qual,
            0,
            0,
            &mate2_trans,
            &mate1_trans,
            &genome,
            255,
            false,
            true,
            -7,
            1,
            1,
            1,   // ih_start
            190, // combined_score (100+90)
            SamAttributes::STANDARD,
        )
        .unwrap();
        assert!(!rec2.flags().is_mate_reverse_complemented());
        assert!(!rec2.flags().is_reverse_complemented());
    }

    #[test]
    fn test_out_sam_attributes_standard() {
        use cigar::op::{Kind, Op};
        // Default (Standard) → NH, HI, AS, NM present; XS, jM, jI, MD absent
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::STANDARD,
        )
        .unwrap();

        let data = record.data();
        // Standard tags present
        assert!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).is_some(),
            "NH should be present"
        );
        assert!(data.get(&Tag::HIT_INDEX).is_some(), "HI should be present");
        assert!(
            data.get(&Tag::ALIGNMENT_SCORE).is_some(),
            "AS should be present"
        );
        // STAR's Standard preset is NH HI AS nM — the mismatch count, NOT the
        // edit-distance NM (which is opt-in via `NM` / `All`).
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_none(),
            "NM (edit distance) should be absent from Standard"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_some(),
            "nM should be present in Standard"
        );
        // Non-standard tags absent
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'M')).is_none(),
            "jM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'I')).is_none(),
            "jI should be absent"
        );
        assert!(
            data.get(&Tag::new(b'M', b'D')).is_none(),
            "MD should be absent"
        );
    }

    #[test]
    fn test_out_sam_attributes_none() {
        use cigar::op::{Kind, Op};
        // None → no optional tags at all
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 1,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let empty_attrs = SamAttributes::empty();
        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            empty_attrs,
        )
        .unwrap();

        let data = record.data();
        assert!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).is_none(),
            "NH should be absent"
        );
        assert!(data.get(&Tag::HIT_INDEX).is_none(), "HI should be absent");
        assert!(
            data.get(&Tag::ALIGNMENT_SCORE).is_none(),
            "AS should be absent"
        );
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_none(),
            "NM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_none(),
            "nM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'M')).is_none(),
            "jM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'I')).is_none(),
            "jI should be absent"
        );
        assert!(
            data.get(&Tag::new(b'M', b'D')).is_none(),
            "MD should be absent"
        );
    }

    #[test]
    fn test_out_sam_attributes_explicit() {
        use cigar::op::{Kind, Op};
        // Explicit ["NH", "MD"] → only NH and MD present
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let attrs = SamAttributes::NH | SamAttributes::MD;
        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            attrs,
        )
        .unwrap();

        let data = record.data();
        // Only NH and MD present
        assert!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).is_some(),
            "NH should be present"
        );
        assert!(
            data.get(&Tag::new(b'M', b'D')).is_some(),
            "MD should be present"
        );
        // Others absent
        assert!(data.get(&Tag::HIT_INDEX).is_none(), "HI should be absent");
        assert!(
            data.get(&Tag::ALIGNMENT_SCORE).is_none(),
            "AS should be absent"
        );
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_none(),
            "NM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_none(),
            "nM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'M')).is_none(),
            "jM should be absent"
        );
        assert!(
            data.get(&Tag::new(b'j', b'I')).is_none(),
            "jI should be absent"
        );
    }

    #[test]
    fn test_sam_attribute_set_expansion() {
        use crate::params::Parameters;

        // Standard
        let p = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "r.fq"]);
        assert_eq!(p.out_sam_attributes, SamAttributes::STANDARD);

        // All — note: XS is stripped at parse-time because the default
        // --outSAMstrandField is "None" (STAR only emits XS in intronMotif mode).
        let p = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "All",
        ]);
        assert_eq!(p.out_sam_attributes, SamAttributes::ALL - SamAttributes::XS);

        // All + intronMotif keeps XS.
        let p = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "All",
            "--outSAMstrandField",
            "intronMotif",
        ]);
        assert_eq!(p.out_sam_attributes, SamAttributes::ALL);

        // None
        let p = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "None",
        ]);
        assert!(p.out_sam_attributes.is_empty());

        // Explicit subset
        let p = Parameters::parse_from([
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "NH",
            "AS",
            "MD",
        ]);
        assert_eq!(
            p.out_sam_attributes,
            SamAttributes::NH | SamAttributes::AS | SamAttributes::MD,
        );
    }

    #[test]
    fn test_nm_vs_nm_mismatch_difference() {
        use cigar::op::{Kind, Op};
        // Verify NM (edit distance) ≠ nM (mismatches only) when indels present
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 60,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 20),
                Op::new(Kind::Insertion, 5),
                Op::new(Kind::Match, 10),
                Op::new(Kind::Deletion, 3),
                Op::new(Kind::Match, 20),
            ],
            score: 80,
            n_mismatch: 2,
            n_gap: 2,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            SamAttributes::ALL,
        )
        .unwrap();

        let data = record.data();
        assert_eq!(
            data.get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(10_i32)),
            "NM should be edit distance: 2 mismatches + 5 ins + 3 del = 10"
        );
        assert_eq!(
            data.get(&Tag::new(b'n', b'M')),
            Some(&Value::from(2_i32)),
            "nM should be 2 (mismatches only)"
        );
    }

    #[test]
    fn test_nm_tag_emits_both_edit_distance_and_mismatch_count() {
        use cigar::op::{Kind, Op};
        let genome = make_test_genome();

        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 103,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 80),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 21),
            ],
            score: 100,
            n_mismatch: 1,
            n_gap: 1,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        // NM (edit distance) and NMM (nM, mismatch count) are distinct flags now;
        // request both to emit both tags.
        let attrs = SamAttributes::NM | SamAttributes::NMM;
        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            attrs,
        )
        .unwrap();
        let data = record.data();
        assert_eq!(
            data.get(&Tag::EDIT_DISTANCE),
            Some(&Value::from(3_i32)),
            "NM:i: should be 1 mismatch + 2 deleted bases = 3"
        );
        assert_eq!(
            data.get(&Tag::new(b'n', b'M')),
            Some(&Value::from(1_i32)),
            "nM:i: should be 1 (mismatches only)"
        );

        // When NM bit is absent, neither tag is emitted.
        let no_nm = SamAttributes::NH | SamAttributes::HI;
        let record = transcript_to_record(
            &transcript,
            "read1",
            &read_seq,
            &read_qual,
            &genome,
            255,
            1,
            1,
            1, // ih_start
            no_nm,
        )
        .unwrap();
        let data = record.data();
        assert!(
            data.get(&Tag::EDIT_DISTANCE).is_none(),
            "NM:i: absent when not requested"
        );
        assert!(
            data.get(&Tag::new(b'n', b'M')).is_none(),
            "nM:i: absent when not requested"
        );
    }

    #[test]
    fn test_build_half_mapped_flags() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "r1.fq", "r2.fq"]);

        let mapped_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate1_seq = vec![0, 1, 2, 3]; // ACGT
        let mate1_qual = vec![30, 30, 30, 30];
        let mate2_seq = vec![3, 2, 1, 0]; // TGCA
        let mate2_qual = vec![30, 30, 30, 30];

        // mate1 is mapped
        let records = SamWriter::build_half_mapped_records(
            "read1",
            &mate1_seq,
            &mate1_qual,
            &mate2_seq,
            &mate2_qual,
            0,
            0,
            0,
            0,
            &mapped_transcript,
            true, // mate1_is_mapped
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 2);

        // Mate1 record (mapped): 0x1 | 0x8 | 0x40 = paired + mate_unmapped + first
        let rec1 = &records[0];
        assert!(rec1.flags().is_segmented());
        assert!(rec1.flags().is_mate_unmapped());
        assert!(rec1.flags().is_first_segment());
        assert!(!rec1.flags().is_unmapped());
        assert!(!rec1.flags().is_last_segment());

        // Mate2 record (unmapped): 0x1 | 0x4 | 0x80 = paired + unmapped + last
        let rec2 = &records[1];
        assert!(rec2.flags().is_segmented());
        assert!(rec2.flags().is_unmapped());
        assert!(rec2.flags().is_last_segment());
        assert!(!rec2.flags().is_first_segment());
        assert!(!rec2.flags().is_mate_unmapped());
    }

    #[test]
    fn test_build_half_mapped_rnext_pnext() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "r1.fq", "r2.fq"]);

        let mapped_transcript = Transcript {
            chr_idx: 0,
            genome_start: 2,
            genome_end: 6,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 2,
                genome_end: 6,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate1_seq = vec![0, 1, 2, 3];
        let mate1_qual = vec![30, 30, 30, 30];
        let mate2_seq = vec![3, 2, 1, 0];
        let mate2_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_half_mapped_records(
            "read1",
            &mate1_seq,
            &mate1_qual,
            &mate2_seq,
            &mate2_qual,
            0,
            0,
            0,
            0,
            &mapped_transcript,
            true,
            &genome,
            &params,
            1,
        )
        .unwrap();

        // mapped_pos = genome_start(2) - chr_start(0) + 1 = 3
        let expected_pos = 3usize;

        // Mapped mate: RNAME = chr_idx(0), POS = 3
        let mapped = &records[0];
        assert_eq!(mapped.reference_sequence_id(), Some(0));
        assert_eq!(
            mapped.alignment_start().map(usize::from),
            Some(expected_pos)
        );
        // RNEXT and PNEXT should point to own position (STAR convention)
        assert_eq!(mapped.mate_reference_sequence_id(), Some(0));
        assert_eq!(
            mapped.mate_alignment_start().map(usize::from),
            Some(expected_pos)
        );

        // Unmapped mate: co-located at mapped mate's position
        let unmapped = &records[1];
        assert_eq!(unmapped.reference_sequence_id(), Some(0));
        assert_eq!(
            unmapped.alignment_start().map(usize::from),
            Some(expected_pos)
        );
        assert_eq!(unmapped.mate_reference_sequence_id(), Some(0));
        assert_eq!(
            unmapped.mate_alignment_start().map(usize::from),
            Some(expected_pos)
        );
        // MAPQ = 0 for unmapped
        assert_eq!(unmapped.mapping_quality().map(u8::from), Some(0));
    }

    #[test]
    fn test_build_half_mapped_mate_order() {
        use crate::align::transcript::Exon;
        use cigar::op::{Kind, Op};

        let genome = make_test_genome();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "r1.fq", "r2.fq"]);

        let mapped_transcript = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 4,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 0,
                genome_end: 4,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        let mate1_seq = vec![0, 1, 2, 3];
        let mate1_qual = vec![30, 30, 30, 30];
        let mate2_seq = vec![3, 2, 1, 0];
        let mate2_qual = vec![30, 30, 30, 30];

        // When mate1 is mapped: mate1 comes first, mate2 second
        let records_m1 = SamWriter::build_half_mapped_records(
            "read1",
            &mate1_seq,
            &mate1_qual,
            &mate2_seq,
            &mate2_qual,
            0,
            0,
            0,
            0,
            &mapped_transcript,
            true,
            &genome,
            &params,
            1,
        )
        .unwrap();
        assert!(records_m1[0].flags().is_first_segment()); // First record = mate1
        assert!(records_m1[1].flags().is_last_segment()); // Second record = mate2

        // When mate2 is mapped: mate1 still comes first, mate2 second
        let records_m2 = SamWriter::build_half_mapped_records(
            "read1",
            &mate1_seq,
            &mate1_qual,
            &mate2_seq,
            &mate2_qual,
            0,
            0,
            0,
            0,
            &mapped_transcript,
            false,
            &genome,
            &params,
            1,
        )
        .unwrap();
        assert!(records_m2[0].flags().is_first_segment()); // First record = mate1 (unmapped)
        assert!(records_m2[1].flags().is_last_segment()); // Second record = mate2 (mapped)
        assert!(records_m2[0].flags().is_unmapped()); // mate1 is unmapped
        assert!(!records_m2[1].flags().is_unmapped()); // mate2 is mapped
    }

    fn spliced_gtag_transcript() -> Transcript {
        use cigar::op::{Kind, Op};
        Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 200,
            is_reverse: false,
            exons: vec![],
            cigar: vec![
                Op::new(Kind::Match, 25),
                Op::new(Kind::Skip, 100),
                Op::new(Kind::Match, 25),
            ],
            score: 50,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 1,
            junction_motifs: vec![SpliceMotif::GtAg],
            junction_annotated: vec![false],
        }
    }

    #[test]
    fn xs_tag_emitted_when_strand_field_intron_motif_only() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec![
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--outSAMstrandField",
            "intronMotif",
        ]);

        let transcripts = vec![spliced_gtag_transcript()];
        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        let data = records[0].data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            Some(&Value::Character(b'+')),
            "XS:A:+ should be emitted when --outSAMstrandField intronMotif is set"
        );
    }

    #[test]
    fn xs_tag_absent_when_xs_in_attrs_but_no_intron_motif() {
        // XS listed in --outSAMattributes without --outSAMstrandField intronMotif:
        // the bit is stripped at parse time, so no XS tag is emitted.
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec![
            "rustar-aligner",
            "--readFilesIn",
            "r.fq",
            "--outSAMattributes",
            "NH",
            "HI",
            "XS",
        ]);

        let transcripts = vec![spliced_gtag_transcript()];
        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        let data = records[0].data();
        assert!(
            data.get(&Tag::new(b'X', b'S')).is_none(),
            "XS:A: should be absent when --outSAMstrandField is not intronMotif"
        );
    }

    #[test]
    fn xs_tag_absent_when_neither_xs_nor_intron_motif_requested() {
        let genome = make_test_genome();
        let params = Parameters::parse_from(vec!["rustar-aligner", "--readFilesIn", "r.fq"]);

        let transcripts = vec![spliced_gtag_transcript()];
        let read_seq = vec![0, 1, 2, 3];
        let read_qual = vec![30, 30, 30, 30];

        let records = SamWriter::build_alignment_records(
            "read1",
            &read_seq,
            &read_qual,
            0,
            0,
            &transcripts,
            &genome,
            &params,
            1,
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        let data = records[0].data();
        assert_eq!(
            data.get(&Tag::new(b'X', b'S')),
            None,
            "XS should not be emitted under default attributes and strand field"
        );
    }
}
