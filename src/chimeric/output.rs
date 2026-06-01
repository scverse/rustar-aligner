// Chimeric.out.junction file writer and WithinBAM record builder

use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
use crate::error::Error;
use crate::genome::Genome;
use bstr::BString;
use noodles::sam;
use noodles::sam::alignment::record::MappingQuality;
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
}

/// Build two SAM records for `--chimOutType WithinBAM`.
///
/// Returns `[donor_record, acceptor_record]`:
/// - Donor: normal FLAGS; full read sequence; SA tag pointing to acceptor.
/// - Acceptor: FLAG 0x0800 (supplementary); empty SEQ/QUAL; SA tag pointing to donor.
pub fn build_within_bam_records(
    alignment: &ChimericAlignment,
    genome: &Genome,
    mapq: u8,
) -> Result<Vec<RecordBuf>, Error> {
    let donor = &alignment.donor;
    let acceptor = &alignment.acceptor;

    let donor_sa = format_sa_entry(donor, &genome.chr_name, &genome.chr_start, mapq);
    let acceptor_sa = format_sa_entry(acceptor, &genome.chr_name, &genome.chr_start, mapq);

    let donor_record = build_segment_record(
        &alignment.read_name,
        &alignment.read_seq,
        donor,
        genome,
        mapq,
        false,
        &acceptor_sa,
    )?;
    let acceptor_record = build_segment_record(
        &alignment.read_name,
        &alignment.read_seq,
        acceptor,
        genome,
        mapq,
        true,
        &donor_sa,
    )?;

    Ok(vec![donor_record, acceptor_record])
}

/// Format one SA tag entry: `chr,pos,strand,CIGAR,mapQ,NM;`
fn format_sa_entry(
    seg: &ChimericSegment,
    chr_names: &[String],
    chr_starts: &[u64],
    mapq: u8,
) -> String {
    let chr = &chr_names[seg.chr_idx];
    let chr_start = chr_starts[seg.chr_idx];
    let pos = seg.genome_start - chr_start + 1; // 1-based per-chr
    let strand = if seg.is_reverse { '-' } else { '+' };
    let cigar = seg.cigar_string();
    format!(
        "{},{},{},{},{},{};",
        chr, pos, strand, cigar, mapq, seg.n_mismatch
    )
}

/// Build one SAM record for a chimeric segment.
fn build_segment_record(
    read_name: &str,
    read_seq: &[u8],
    seg: &ChimericSegment,
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

    *record.cigar_mut() = seg.cigar.iter().copied().collect();

    // Primary record carries the full read sequence; supplementary uses * (empty).
    if !is_supplementary {
        if seg.is_reverse {
            let seq_bytes: Vec<u8> = read_seq
                .iter()
                .rev()
                .map(|&b| decode_base(complement_base(b)))
                .collect();
            *record.sequence_mut() = Sequence::from(seq_bytes);
        } else {
            let seq_bytes: Vec<u8> = read_seq.iter().map(|&b| decode_base(b)).collect();
            *record.sequence_mut() = Sequence::from(seq_bytes);
        }
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
            sequence: vec![0u8; 2048],
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
        let records = build_within_bam_records(&alignment, &genome, 255).unwrap();

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
        let records = build_within_bam_records(&alignment, &genome, 255).unwrap();

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
        let records = build_within_bam_records(&alignment, &genome, 255).unwrap();

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
        let records = build_within_bam_records(&alignment, &genome, 255).unwrap();

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
}
