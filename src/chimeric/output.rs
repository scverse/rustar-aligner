// Chimeric.out.junction file writer and WithinBAM record builder

use crate::align::transcript::cigar_to_string;
use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
use crate::error::Error;
use crate::genome::Genome;
use bstr::BString;
use noodles::sam;
use noodles::sam::alignment::record::{MappingQuality, cigar};
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Writer for Chimeric.out.junction file
pub struct ChimericJunctionWriter {
    writer: BufWriter<File>,
}

impl ChimericJunctionWriter {
    /// Create a new chimeric junction writer
    ///
    /// Creates file: {prefix}Chimeric.out.junction
    pub fn new(prefix: &str) -> Result<Self, Error> {
        let path = PathBuf::from(format!("{prefix}Chimeric.out.junction"));

        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent).map_err(|e| Error::io(e, parent))?;
        }

        let file = File::create(&path).map_err(|e| Error::io(e, &path))?;

        let writer = BufWriter::new(file);
        Ok(Self { writer })
    }

    /// Write a chimeric alignment to the file
    ///
    /// Format: 14 tab-separated columns
    /// 1. Donor chromosome
    /// 2. Donor breakpoint (1-based)
    /// 3. Donor strand (+/-)
    /// 4. Acceptor chromosome
    /// 5. Acceptor breakpoint (1-based)
    /// 6. Acceptor strand (+/-)
    /// 7. Junction type (0-6)
    /// 8. Repeat length donor
    /// 9. Repeat length acceptor
    /// 10. Read name
    /// 11. First segment start (1-based)
    /// 12. First segment CIGAR
    /// 13. Second segment start (1-based)
    /// 14. Second segment CIGAR
    pub fn write_alignment(
        &mut self,
        alignment: &ChimericAlignment,
        chr_names: &[String],
        read_name: &str,
    ) -> Result<(), Error> {
        // Get chromosome names
        let donor_chr = &chr_names[alignment.donor.chr_idx];
        let acceptor_chr = &chr_names[alignment.acceptor.chr_idx];

        // Get breakpoints (1-based)
        let donor_bp = alignment.donor_breakpoint();
        let acceptor_bp = alignment.acceptor_breakpoint();

        // Get strand symbols
        let donor_strand = alignment.donor_strand();
        let acceptor_strand = alignment.acceptor_strand();

        // Get junction type
        let junction_type = alignment.junction_type;

        // Get repeat lengths
        let repeat_donor = alignment.repeat_len_donor;
        let repeat_acceptor = alignment.repeat_len_acceptor;

        // Get segment start positions (1-based)
        let donor_start = alignment.donor.genome_start + 1;
        let acceptor_start = alignment.acceptor.genome_start + 1;

        // Convert CIGAR to string
        let donor_cigar = alignment.donor.cigar_string();
        let acceptor_cigar = alignment.acceptor.cigar_string();

        // Write line
        writeln!(
            self.writer,
            "{donor_chr}\t{donor_bp}\t{donor_strand}\t{acceptor_chr}\t{acceptor_bp}\t{acceptor_strand}\t{junction_type}\t{repeat_donor}\t{repeat_acceptor}\t{read_name}\t{donor_start}\t{donor_cigar}\t{acceptor_start}\t{acceptor_cigar}",
        )
        .map_err(|e| Error::Chimeric(format!("Failed to write chimeric junction: {e}")))?;

        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<(), Error> {
        self.writer
            .flush()
            .map_err(|e| Error::Chimeric(format!("Failed to flush chimeric junction file: {e}")))
    }

    /// `--chimOutJunctionFormat 1` (STAR-Fusion): append a comment header line (version + command
    /// line) then a `# Nreads .. NreadsUnique .. NreadsMulti ..` count summary, after all junction
    /// lines have been written.
    pub fn write_format1_trailer(
        &mut self,
        command_line: &str,
        n_reads: u64,
        n_reads_unique: u64,
        n_reads_multi: u64,
    ) -> Result<(), Error> {
        writeln!(
            self.writer,
            "# {}   {command_line}",
            env!("CARGO_PKG_VERSION")
        )
        .map_err(|e| Error::Chimeric(format!("Failed to write chimeric format-1 header: {e}")))?;
        writeln!(
            self.writer,
            "# Nreads {n_reads}\tNreadsUnique {n_reads_unique}\tNreadsMulti {n_reads_multi}"
        )
        .map_err(|e| Error::Chimeric(format!("Failed to write chimeric format-1 counts: {e}")))
    }
}

/// Build the SAM records for `--chimOutType WithinBAM`.
///
/// Mirrors STAR's `ChimericAlign::chimericBAMoutput` + `ReadAlign::alignBAM`.
/// Returns two records: `[donor, acceptor]` for paired-end, and for
/// single-end the two segments in read order (STAR's `trChim[0]`,
/// `trChim[1]`):
///
/// - Every CIGAR accounts for the whole read (the part outside the segment is
///   clipped), so CIGAR query length always matches SEQ.
/// - Single-end (`chimType==3`): the segment with the higher score is the
///   representative record (donor only if its score is strictly higher); the
///   other one is supplementary (FLAG 0x800). With `hard_clip` (STAR's
///   default `HardClip`) the supplementary record is hard-clipped on its
///   chimeric-junction side (alignType -11/-12) and its SEQ drops those
///   bases; with `SoftClip` (alignType -13) it keeps soft clips and the full
///   SEQ.
/// - Paired-end: donor is representative, acceptor is supplementary with an
///   empty SEQ (unchanged legacy behavior).
/// - Each of the two records carries an SA tag describing the other one
///   (`chr,pos,strand,CIGAR,MAPQ,NM;`), using that record's final CIGAR.
pub fn build_within_bam_records(
    alignment: &ChimericAlignment,
    genome: &Genome,
    mapq: u8,
    single_end: bool,
    hard_clip: bool,
) -> Result<Vec<RecordBuf>, Error> {
    use cigar::op::Kind;

    let read_len = alignment.read_seq.len();
    let mut segs = [&alignment.donor, &alignment.acceptor];
    let mut cigars = segs.map(|seg| full_length_cigar(seg, read_len));

    if single_end {
        // STAR orders trChim by read position (roStart); order by the 5' clip of
        // the full-length CIGAR so this holds whichever detection tier built the pair.
        let ro_start = |seg: &ChimericSegment, ops: &[cigar::Op]| {
            let clip = |op: Option<&cigar::Op>| {
                op.filter(|op| op.kind() == Kind::SoftClip)
                    .map_or(0, |op| op.len())
            };
            if seg.is_reverse {
                clip(ops.last())
            } else {
                clip(ops.first())
            }
        };
        if ro_start(segs[0], &cigars[0]) > ro_start(segs[1], &cigars[1]) {
            segs.swap(0, 1);
            cigars.swap(0, 1);
        }
    }

    // STAR SE: chimRepresent = (trChim[0].maxScore > trChim[1].maxScore) ? 0 : 1
    let represent = usize::from(single_end && segs[0].score <= segs[1].score);
    // Bases hard-clipped from the start/end of SEQ (in CIGAR orientation).
    let mut seq_trim = [(0usize, 0usize); 2];

    if single_end && hard_clip {
        let isuppl = 1 - represent;
        let seg = segs[isuppl];
        // STAR: alignType = (itr%2 == Str) ? -12 (hard clip right) : -11 (hard clip left)
        let hard_right = (isuppl == 1) == seg.is_reverse;
        let ops = &mut cigars[isuppl];
        if hard_right {
            if let Some(last) = ops.last_mut()
                && last.kind() == Kind::SoftClip
            {
                seq_trim[isuppl].1 = last.len();
                *last = cigar::Op::new(Kind::HardClip, last.len());
            }
        } else if let Some(first) = ops.first_mut()
            && first.kind() == Kind::SoftClip
        {
            seq_trim[isuppl].0 = first.len();
            *first = cigar::Op::new(Kind::HardClip, first.len());
        }
    }

    let sa = [0, 1].map(|i| format_sa_entry(segs[i], &cigars[i], genome, mapq));

    let mut records = Vec::with_capacity(2);
    for i in 0..2 {
        let is_supplementary = i != represent;
        let with_seq = single_end || !is_supplementary;
        records.push(build_segment_record(
            &alignment.read_name,
            &alignment.read_seq,
            segs[i],
            &cigars[i],
            if with_seq { Some(seq_trim[i]) } else { None },
            genome,
            mapq,
            is_supplementary,
            &sa[1 - i],
        )?);
    }

    Ok(records)
}

/// Segment CIGAR padded with soft clips so that it covers the whole read.
///
/// Segments built from full-read transcripts already do. Segments re-seeded
/// from a sub-sequence (soft-clip / residual re-mapping) may not; their
/// `read_start` is the left clip in CIGAR orientation.
fn full_length_cigar(seg: &ChimericSegment, read_len: usize) -> Vec<cigar::Op> {
    use cigar::op::Kind;

    let consumes_read = |k: Kind| {
        matches!(
            k,
            Kind::Match
                | Kind::Insertion
                | Kind::SoftClip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        )
    };
    let query_len: usize = seg
        .cigar
        .iter()
        .filter(|op| consumes_read(op.kind()))
        .map(|op| op.len())
        .sum();
    if query_len == read_len {
        return seg.cigar.clone();
    }

    // Strip existing clips, then re-pad to the full read length.
    let core: Vec<cigar::Op> = seg
        .cigar
        .iter()
        .copied()
        .filter(|op| !matches!(op.kind(), Kind::SoftClip | Kind::HardClip))
        .collect();
    let core_len: usize = core
        .iter()
        .filter(|op| consumes_read(op.kind()))
        .map(|op| op.len())
        .sum();
    let left = seg.read_start.min(read_len.saturating_sub(core_len));
    let right = read_len.saturating_sub(left + core_len);

    let mut ops = Vec::with_capacity(core.len() + 2);
    if left > 0 {
        ops.push(cigar::Op::new(Kind::SoftClip, left));
    }
    ops.extend(core);
    if right > 0 {
        ops.push(cigar::Op::new(Kind::SoftClip, right));
    }
    ops
}

/// Format one SA tag entry: `chr,pos,strand,CIGAR,mapQ,NM;`
fn format_sa_entry(seg: &ChimericSegment, ops: &[cigar::Op], genome: &Genome, mapq: u8) -> String {
    let chr = &genome.chr_name[seg.chr_idx];
    let pos = seg.genome_start - genome.chr_start[seg.chr_idx] + 1; // 1-based per-chr
    let strand = if seg.is_reverse { '-' } else { '+' };
    let cigar = cigar_to_string(ops);
    format!(
        "{},{},{},{},{},{};",
        chr, pos, strand, cigar, mapq, seg.n_mismatch
    )
}

/// Build one SAM record for a chimeric segment.
///
/// `seq_trim` is `None` for an empty SEQ (`*`), or the number of bases to drop
/// from the start / end of the CIGAR-oriented read (hard clips).
#[allow(clippy::too_many_arguments)]
fn build_segment_record(
    read_name: &str,
    read_seq: &[u8],
    seg: &ChimericSegment,
    ops: &[cigar::Op],
    seq_trim: Option<(usize, usize)>,
    genome: &Genome,
    mapq: u8,
    is_supplementary: bool,
    sa_tag: &str,
) -> Result<RecordBuf, Error> {
    use crate::io::fastq::{complement_base, decode_base};
    use noodles::sam::alignment::record::data::field::Tag;

    let mut record = RecordBuf::default();
    record.name_mut().replace(read_name.into());

    let mut flags = sam::alignment::record::Flags::empty();
    if seg.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
    }
    if is_supplementary {
        flags |= sam::alignment::record::Flags::SUPPLEMENTARY;
    }
    *record.flags_mut() = flags;

    *record.reference_sequence_id_mut() = Some(seg.chr_idx);

    let chr_start = genome.chr_start[seg.chr_idx];
    let pos = (seg.genome_start - chr_start + 1) as usize;
    *record.alignment_start_mut() = Some(
        pos.try_into()
            .map_err(|e| Error::Chimeric(format!("invalid chimeric position {pos}: {e}")))?,
    );

    *record.mapping_quality_mut() = MappingQuality::new(mapq);

    *record.cigar_mut() = ops.iter().copied().collect();

    if let Some((trim_left, trim_right)) = seq_trim {
        let oriented: Vec<u8> = if seg.is_reverse {
            read_seq
                .iter()
                .rev()
                .map(|&b| decode_base(complement_base(b)))
                .collect()
        } else {
            read_seq.iter().map(|&b| decode_base(b)).collect()
        };
        let end = oriented.len().saturating_sub(trim_right).max(trim_left);
        *record.sequence_mut() = Sequence::from(oriented[trim_left..end].to_vec());
        // Leave QUAL empty (not available for chimeric segments)
        *record.quality_scores_mut() = QualityScores::default();
    }

    let data = record.data_mut();
    data.insert(Tag::new(b'S', b'A'), Value::String(BString::from(sa_tag)));
    data.insert(Tag::new(b'N', b'M'), Value::from(seg.n_mismatch as i32));
    data.insert(Tag::ALIGNMENT_SCORE, Value::from(seg.score));

    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
    use noodles::sam::alignment::record::cigar;
    use std::io::Read;
    use tempfile::tempdir;

    #[test]
    fn test_chimeric_junction_writer_creation() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());

        let writer = ChimericJunctionWriter::new(&prefix);
        assert!(writer.is_ok());

        let path = PathBuf::from(format!("{prefix}Chimeric.out.junction"));
        assert!(path.exists());
    }

    #[test]
    fn test_chimeric_junction_writer_bare_dot_prefix() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/SAMPLE.", dir.path().display());

        let writer = ChimericJunctionWriter::new(&prefix);
        assert!(writer.is_ok());

        let path = PathBuf::from(format!("{prefix}Chimeric.out.junction"));
        assert!(path.exists(), "expected {} to exist", path.display());
        assert!(
            path.file_name().unwrap().to_str().unwrap() == "SAMPLE.Chimeric.out.junction",
            "expected literal concatenation, got {}",
            path.display()
        );
    }

    #[test]
    fn test_chimeric_junction_writer_creates_missing_parent() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/sample/", dir.path().display());
        let prefix_path = PathBuf::from(&prefix);

        assert!(!prefix_path.exists(), "parent dir should not exist yet");

        let writer = ChimericJunctionWriter::new(&prefix);
        assert!(
            writer.is_ok(),
            "writer should create missing parent dir, got: {:?}",
            writer.err()
        );

        let mut path = prefix_path.clone();
        path.push("Chimeric.out.junction");
        assert!(path.exists(), "chim output file should exist at {path:?}");
    }

    #[test]
    fn test_chim_out_junction_format1_trailer() {
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());
        let path = PathBuf::from(format!("{prefix}Chimeric.out.junction"));

        let mut writer = ChimericJunctionWriter::new(&prefix).unwrap();
        writer
            .write_format1_trailer("star --runMode alignReads", 100, 5, 0)
            .unwrap();
        writer.flush().unwrap();

        let mut contents = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("# "));
        assert!(lines[0].contains("star --runMode alignReads"));
        assert_eq!(lines[1], "# Nreads 100\tNreadsUnique 5\tNreadsMulti 0");
    }

    #[test]
    fn test_write_inter_chromosomal() {
        use cigar::op::{Kind, Op};
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());

        let mut writer = ChimericJunctionWriter::new(&prefix).unwrap();

        // Create mock chimeric alignment (chr9 -> chr22, BCR-ABL fusion)
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 133_738_300,
            genome_end: 133_738_363,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            cigar: vec![Op::new(Kind::Match, 63)],
            score: 100,
            n_mismatch: 2,
        };

        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 23_632_600,
            genome_end: 23_632_637,
            is_reverse: false,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::Match, 37)],
            score: 80,
            n_mismatch: 1,
        };

        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            1, // GT/AG
            0,
            0,
            vec![0; 100],
            "READ_001".to_string(),
        );

        let chr_names = vec!["chr9".to_string(), "chr22".to_string()];

        writer
            .write_alignment(&alignment, &chr_names, "READ_001")
            .unwrap();
        writer.flush().unwrap();

        // Read file and verify
        let path = PathBuf::from(format!("{prefix}Chimeric.out.junction"));

        let mut content = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();

        let line = content.trim();
        let fields: Vec<&str> = line.split('\t').collect();

        assert_eq!(fields.len(), 14);
        assert_eq!(fields[0], "chr9"); // donor chr
        assert_eq!(fields[1], "133738363"); // donor breakpoint
        assert_eq!(fields[2], "+"); // donor strand
        assert_eq!(fields[3], "chr22"); // acceptor chr
        assert_eq!(fields[4], "23632601"); // acceptor breakpoint
        assert_eq!(fields[5], "+"); // acceptor strand
        assert_eq!(fields[6], "1"); // junction type
        assert_eq!(fields[7], "0"); // repeat donor
        assert_eq!(fields[8], "0"); // repeat acceptor
        assert_eq!(fields[9], "READ_001"); // read name
        assert_eq!(fields[10], "133738301"); // donor start (1-based)
        assert_eq!(fields[11], "63M"); // donor CIGAR
        assert_eq!(fields[12], "23632601"); // acceptor start (1-based)
        assert_eq!(fields[13], "37M"); // acceptor CIGAR
    }

    #[test]
    fn test_write_strand_break() {
        use cigar::op::{Kind, Op};
        let dir = tempdir().unwrap();
        let prefix = format!("{}/", dir.path().display());

        let mut writer = ChimericJunctionWriter::new(&prefix).unwrap();

        // Create mock chimeric alignment (same chr, opposite strands)
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1050,
            is_reverse: false,
            read_start: 0,
            read_end: 50,
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 1,
        };

        let acceptor = ChimericSegment {
            chr_idx: 0,
            genome_start: 2000,
            genome_end: 2050,
            is_reverse: true,
            read_start: 50,
            read_end: 100,
            cigar: vec![Op::new(Kind::Match, 50)],
            score: 100,
            n_mismatch: 1,
        };

        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            0, // non-canonical
            0,
            0,
            vec![0; 100],
            "READ_002".to_string(),
        );

        let chr_names = vec!["chr1".to_string()];

        writer
            .write_alignment(&alignment, &chr_names, "READ_002")
            .unwrap();
        writer.flush().unwrap();

        // Read file and verify
        let path = PathBuf::from(format!("{prefix}Chimeric.out.junction"));

        let mut content = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();

        let line = content.trim();
        let fields: Vec<&str> = line.split('\t').collect();

        assert_eq!(fields.len(), 14);
        assert_eq!(fields[0], "chr1"); // donor chr
        assert_eq!(fields[2], "+"); // donor strand
        assert_eq!(fields[3], "chr1"); // acceptor chr
        assert_eq!(fields[5], "-"); // acceptor strand (reverse)
        assert_eq!(fields[6], "0"); // junction type (non-canonical)
    }

    // --- build_within_bam_records tests ---

    fn make_genome_2chr() -> crate::genome::Genome {
        use crate::genome::Genome;
        Genome {
            transform_blocks: None,
            sequence: vec![0u8; 2048].into(),
            n_genome: 1024,
            n_genome_real: 1024,
            n_chr_real: 2,
            chr_name: vec!["chr9".to_string(), "chr22".to_string()],
            chr_length: vec![512, 512],
            chr_start: vec![0, 512, 1024],
        }
    }

    #[test]
    fn test_within_bam_returns_two_records() {
        use cigar::op::{Kind, Op};
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 163,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            cigar: vec![Op::new(Kind::Match, 63)],
            score: 63,
            n_mismatch: 0,
        };
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: false,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 1,
        };
        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            0,
            0,
            0,
            vec![0u8; 100],
            "READ_001".to_string(),
        );
        let genome = make_genome_2chr();
        let records = build_within_bam_records(&alignment, &genome, 255, false, true).unwrap();

        assert_eq!(records.len(), 2);
    }

    #[test]
    fn test_within_bam_donor_not_supplementary() {
        use cigar::op::{Kind, Op};
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 163,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            cigar: vec![Op::new(Kind::Match, 63)],
            score: 63,
            n_mismatch: 0,
        };
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: false,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 1,
        };
        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            0,
            0,
            0,
            vec![0u8; 100],
            "READ_001".to_string(),
        );
        let genome = make_genome_2chr();
        let records = build_within_bam_records(&alignment, &genome, 255, false, true).unwrap();

        let donor_flags = records[0].flags();
        let acceptor_flags = records[1].flags();

        assert!(
            !donor_flags.is_supplementary(),
            "donor must not be supplementary"
        );
        assert!(
            acceptor_flags.is_supplementary(),
            "acceptor must be supplementary (0x800)"
        );
    }

    #[test]
    fn test_within_bam_sa_tag_format() {
        use cigar::op::{Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 163,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            cigar: vec![Op::new(Kind::Match, 63)],
            score: 63,
            n_mismatch: 2,
        };
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: true,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 1,
        };
        let alignment = ChimericAlignment::new(
            donor,
            acceptor,
            0,
            0,
            0,
            vec![0u8; 100],
            "READ_001".to_string(),
        );
        let genome = make_genome_2chr();
        let records = build_within_bam_records(&alignment, &genome, 255, false, true).unwrap();

        // Donor record's SA tag should point to acceptor
        let sa_tag = Tag::new(b'S', b'A');
        let donor_sa = records[0].data().get(&sa_tag).unwrap();
        let donor_sa_str = format!("{donor_sa:?}");
        // SA tag: chr22,89,-,37M,255,1; (pos = 600-512+1=89, strand=-, nm=1)
        assert!(
            donor_sa_str.contains("chr22"),
            "SA tag must name acceptor chr"
        );
        assert!(
            donor_sa_str.contains("89"),
            "SA tag must have per-chr position"
        );
        assert!(
            donor_sa_str.contains('-'),
            "SA tag must reflect reverse strand"
        );

        // Acceptor record's SA tag should point to donor
        let acceptor_sa = records[1].data().get(&sa_tag).unwrap();
        let acceptor_sa_str = format!("{acceptor_sa:?}");
        assert!(
            acceptor_sa_str.contains("chr9"),
            "SA tag must name donor chr"
        );
    }

    #[test]
    fn test_within_bam_donor_has_sequence() {
        use cigar::op::{Kind, Op};
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 163,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            cigar: vec![Op::new(Kind::Match, 63)],
            score: 63,
            n_mismatch: 0,
        };
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: false,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 0,
        };
        let read_seq = vec![0u8; 100]; // 100 A bases
        let alignment =
            ChimericAlignment::new(donor, acceptor, 0, 0, 0, read_seq, "READ_001".to_string());
        let genome = make_genome_2chr();
        let records = build_within_bam_records(&alignment, &genome, 255, false, true).unwrap();

        // Donor has sequence, acceptor has empty sequence (*)
        assert!(
            !records[0].sequence().is_empty(),
            "donor record must have SEQ"
        );
        assert!(
            records[1].sequence().is_empty(),
            "supplementary record must have empty SEQ"
        );
    }

    fn cigar_str(rec: &RecordBuf) -> String {
        cigar_to_string(rec.cigar().as_ref())
    }

    fn sa_str(rec: &RecordBuf) -> String {
        match rec.data().get(b"SA") {
            Some(Value::String(s)) => s.to_string(),
            other => panic!("missing SA tag: {other:?}"),
        }
    }

    /// Issue #279: SE segments whose CIGAR only covers the segment (e.g. from
    /// soft-clip re-seeding) must be padded to the read length, and the
    /// supplementary one hard-clipped on its junction side, like STAR.
    fn se_alignment(acceptor_reverse: bool) -> ChimericAlignment {
        use cigar::op::{Kind, Op};
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 163,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            cigar: vec![Op::new(Kind::Match, 63), Op::new(Kind::SoftClip, 37)],
            score: 63,
            n_mismatch: 0,
        };
        // Partial CIGAR (37M): only the segment, as built from a sub-sequence.
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: acceptor_reverse,
            read_start: if acceptor_reverse { 0 } else { 63 },
            read_end: if acceptor_reverse { 37 } else { 100 },
            cigar: vec![Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 1,
        };
        ChimericAlignment::new(
            donor,
            acceptor,
            0,
            0,
            0,
            (0..100u8).map(|i| i % 4).collect(),
            "READ_SE".to_string(),
        )
    }

    #[test]
    fn test_within_bam_se_hard_clip_forward() {
        let genome = make_genome_2chr();
        let records =
            build_within_bam_records(&se_alignment(false), &genome, 255, true, true).unwrap();
        // Representative = higher-scoring donor, full SEQ, soft clips.
        assert!(!records[0].flags().is_supplementary());
        assert_eq!(cigar_str(&records[0]), "63M37S");
        assert_eq!(records[0].sequence().len(), 100);
        // Supplementary acceptor: 5' side (junction) hard-clipped, SEQ trimmed.
        assert!(records[1].flags().is_supplementary());
        assert_eq!(cigar_str(&records[1]), "63H37M");
        assert_eq!(records[1].sequence().len(), 37);
        assert_eq!(sa_str(&records[0]), "chr22,89,+,63H37M,255,1;");
        assert_eq!(sa_str(&records[1]), "chr9,101,+,63M37S,255,0;");
    }

    #[test]
    fn test_within_bam_se_hard_clip_reverse_and_soft_clip() {
        let genome = make_genome_2chr();
        // Reverse acceptor: CIGAR laid out on the reverse-complemented read,
        // so the junction side is on the right.
        let records =
            build_within_bam_records(&se_alignment(true), &genome, 255, true, true).unwrap();
        assert_eq!(cigar_str(&records[1]), "37M63H");
        assert_eq!(records[1].sequence().len(), 37);

        // SoftClip mode (STAR alignType -13): soft clips and full SEQ.
        let records =
            build_within_bam_records(&se_alignment(true), &genome, 255, true, false).unwrap();
        assert!(records[1].flags().is_supplementary());
        assert_eq!(cigar_str(&records[1]), "37M63S");
        assert_eq!(records[1].sequence().len(), 100);
    }

    #[test]
    fn test_within_bam_se_representative_is_higher_score() {
        let genome = make_genome_2chr();
        let mut aln = se_alignment(false);
        aln.acceptor.score = 70;
        let records = build_within_bam_records(&aln, &genome, 255, true, true).unwrap();
        // STAR: chimRepresent = trChim[0].maxScore > trChim[1].maxScore ? 0 : 1
        assert!(records[0].flags().is_supplementary());
        assert!(!records[1].flags().is_supplementary());
        // Donor (5' segment, forward) hard-clipped on its 3' (right) side.
        assert_eq!(cigar_str(&records[0]), "63M37H");
        assert_eq!(records[0].sequence().len(), 63);
        assert_eq!(cigar_str(&records[1]), "63S37M");
    }
}
