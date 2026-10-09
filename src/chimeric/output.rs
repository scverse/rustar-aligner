// Chimeric.out.junction file writer and WithinBAM record builder

use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
use crate::error::Error;
use crate::genome::Genome;
use bstr::BString;
use noodles::sam;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::data::field::Tag;
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
    /// `multimap` selects STAR's wider file: `--chimMultimapNmax > 0` writes a
    /// header line and six extra columns per record
    /// (`ParametersChimeric_initialize.cpp:48-71`).
    pub fn new_with_multimap(prefix: &str, multimap: bool) -> Result<Self, Error> {
        let mut w = Self::new(prefix)?;
        if multimap {
            // STAR names 21 columns here, `readgrp` included, even though the
            // record itself only carries it when a read group is configured.
            writeln!(
                w.writer,
                "chr_donorA\tbrkpt_donorA\tstrand_donorA\tchr_acceptorB\tbrkpt_acceptorB\t\
                 strand_acceptorB\tjunction_type\trepeat_left_lenA\trepeat_right_lenB\t\
                 read_name\tstart_alnA\tcigar_alnA\tstart_alnB\tcigar_alnB\tnum_chim_aln\t\
                 max_poss_aln_score\tnon_chim_aln_score\tthis_chim_aln_score\t\
                 bestall_chim_aln_score\tPEmerged_bool\treadgrp"
            )
            .map_err(|e| Error::Chimeric(format!("Failed to write chimeric header: {e}")))?;
        }
        Ok(w)
    }

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
        chr_starts: &[u64],
        read_name: &str,
    ) -> Result<(), Error> {
        // The STAR detector records the line as STAR writes it, including
        // two-mate segments whose start and CIGAR no single-mate segment holds.
        // Every coordinate in this file is per-chromosome, as STAR's is, so the
        // chromosome's padded start comes off all four positions.
        let (
            donor_chr_idx,
            donor_bp,
            donor_strand,
            acceptor_chr_idx,
            acceptor_bp,
            acceptor_strand,
            donor_start,
            donor_cigar,
            acceptor_start,
            acceptor_cigar,
        ) = if let Some(j) = &alignment.junction_line {
            let (dc, ac) = (chr_starts[j.donor_chr], chr_starts[j.acceptor_chr]);
            let strand = |rev: bool| if rev { '-' } else { '+' };
            (
                j.donor_chr,
                j.donor_break.wrapping_sub(dc).wrapping_add(1),
                strand(j.donor_reverse),
                j.acceptor_chr,
                j.acceptor_break.wrapping_sub(ac).wrapping_add(1),
                strand(j.acceptor_reverse),
                j.donor_start - dc + 1,
                j.donor_cigar.clone(),
                j.acceptor_start - ac + 1,
                j.acceptor_cigar.clone(),
            )
        } else {
            let (d, a) = (&alignment.donor, &alignment.acceptor);
            let (dc, ac) = (chr_starts[d.chr_idx], chr_starts[a.chr_idx]);
            (
                d.chr_idx,
                alignment.donor_breakpoint() - dc,
                alignment.donor_strand(),
                a.chr_idx,
                alignment.acceptor_breakpoint() - ac,
                alignment.acceptor_strand(),
                d.genome_start - dc + 1,
                d.cigar_string(),
                a.genome_start - ac + 1,
                a.cigar_string(),
            )
        };
        let donor_chr = &chr_names[donor_chr_idx];
        let acceptor_chr = &chr_names[acceptor_chr_idx];
        let junction_type = alignment.junction_type;
        let repeat_donor = alignment.repeat_len_donor;
        let repeat_acceptor = alignment.repeat_len_acceptor;

        // Write line. Under `--chimMultimapNmax` STAR appends six run-level
        // columns (`ChimericAlign_chimericJunctionOutput.cpp:14-19`); without it
        // the file stays at the classic 14.
        let extra = match &alignment.multimap {
            Some(m) => format!(
                "\t{}\t{}\t{}\t{}\t{}\t{}",
                m.chim_n,
                m.max_possible_score,
                m.max_non_chim_score,
                m.chim_score,
                m.best_chim_score,
                u8::from(m.pe_merged)
            ),
            None => String::new(),
        };
        writeln!(
            self.writer,
            "{donor_chr}\t{donor_bp}\t{donor_strand}\t{acceptor_chr}\t{acceptor_bp}\t{acceptor_strand}\t{junction_type}\t{repeat_donor}\t{repeat_acceptor}\t{read_name}\t{donor_start}\t{donor_cigar}\t{acceptor_start}\t{acceptor_cigar}{extra}",
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

/// Build two SAM records for `--chimOutType WithinBAM`.
///
/// Returns `[donor_record, acceptor_record]`:
/// - Donor: normal FLAGS; full read sequence; SA tag pointing to acceptor.
/// - Acceptor: FLAG 0x0800 (supplementary); empty SEQ/QUAL; SA tag pointing to donor.
///
/// `ch_tag` adds `ch:A:1` to both (`--outSAMattributes` has `ch`).
pub fn build_within_bam_records(
    alignment: &ChimericAlignment,
    genome: &Genome,
    mapq: u8,
    ch_tag: bool,
) -> Result<Vec<RecordBuf>, Error> {
    let donor = &alignment.donor;
    let acceptor = &alignment.acceptor;

    // Which end of the supplementary segment's CIGAR covers the OTHER segment's
    // bases. STAR hard-clips exactly that end (`chimericBAMoutput.cpp:55`
    // picks -11 for a left junction, -12 for a right one) and drops those bases
    // from SEQ (`ReadAlign_alignBAM.cpp:503-511`), which is the ordinary SAM
    // rule: hard-clipped bases are absent from SEQ, soft-clipped ones present.
    //
    // In original-read order the other segment is 5' of this one when this one
    // trails. The CIGAR is in reference orientation, so a reverse segment has
    // that end on the opposite side.
    let read_len = alignment.read_seq.len();
    let ro = |seg: &ChimericSegment| {
        if seg.is_reverse {
            read_len.saturating_sub(seg.read_end)
        } else {
            seg.read_start
        }
    };
    let acceptor_trails = ro(acceptor) > ro(donor);
    let hard_leading = acceptor_trails != acceptor.is_reverse;

    let donor_sa = format_sa_entry(donor, &genome.chr_name, &genome.chr_start, mapq);
    let acceptor_sa = format_sa_entry(acceptor, &genome.chr_name, &genome.chr_start, mapq);

    let donor_record = build_segment_record(
        &alignment.read_name,
        &alignment.read_seq,
        donor,
        genome,
        mapq,
        None,
        &acceptor_sa,
    )?;
    let acceptor_record = build_segment_record(
        &alignment.read_name,
        &alignment.read_seq,
        acceptor,
        genome,
        mapq,
        Some(hard_leading),
        &donor_sa,
    )?;

    let mut records = vec![donor_record, acceptor_record];
    if ch_tag {
        // `ch:A:1` marks every chimeric record (ReadAlign_alignBAM.cpp, ATTR_ch,
        // alignType<=-10); it is only valid in BAM output.
        for r in &mut records {
            r.data_mut()
                .insert(Tag::new(b'c', b'h'), Value::Character(b'1'));
        }
    }
    Ok(records)
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
    // `None` for the representative record; `Some(hard_leading)` for the
    // supplementary one, naming which CIGAR end to hard-clip.
    supplementary: Option<bool>,
    sa_tag: &str,
) -> Result<RecordBuf, Error> {
    use crate::io::fastq::{complement_base, decode_base};

    let mut record = RecordBuf::default();
    record.name_mut().replace(read_name.into());

    let mut flags = sam::alignment::record::Flags::empty();
    if seg.is_reverse {
        flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
    }
    if supplementary.is_some() {
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

    // SEQ in reference orientation, matching the CIGAR.
    let mut seq_bytes: Vec<u8> = if seg.is_reverse {
        read_seq
            .iter()
            .rev()
            .map(|&b| decode_base(complement_base(b)))
            .collect()
    } else {
        read_seq.iter().map(|&b| decode_base(b)).collect()
    };

    let mut cigar_ops = seg.cigar.clone();

    // Pad the CIGAR out to the full read with soft clips where the segment does
    // not already carry them. STAR always emits these (`trimL1`/`trimR1`,
    // `ReadAlign_alignBAM.cpp:235,269`), but our segments sometimes arrive as a
    // bare match block, and a record whose CIGAR claims fewer query bases than
    // SEQ holds cannot be written. Read coordinates are in original-read order,
    // so a reverse segment's leading clip is the one past its 3' end.
    {
        use noodles::sam::alignment::record::cigar::{Op, op::Kind};
        let q: usize = cigar_ops
            .iter()
            .filter(|op| op.kind().consumes_read())
            .map(|op| op.len())
            .sum();
        if q < read_seq.len() {
            let before = seg.read_start;
            let after = read_seq.len().saturating_sub(seg.read_end);
            let (lead, trail) = if seg.is_reverse {
                (after, before)
            } else {
                (before, after)
            };
            let lead_present = matches!(cigar_ops.first(), Some(o) if o.kind() == Kind::SoftClip);
            let trail_present = matches!(cigar_ops.last(), Some(o) if o.kind() == Kind::SoftClip);
            if trail > 0 && !trail_present {
                cigar_ops.push(Op::new(Kind::SoftClip, trail));
            }
            if lead > 0 && !lead_present {
                cigar_ops.insert(0, Op::new(Kind::SoftClip, lead));
            }
        }
    }

    // Supplementary: hard-clip the end covering the other segment and drop
    // those bases from SEQ. Leaving them soft-clipped with an empty SEQ, as
    // this did before, is not representable -- soft clips consume the query, so
    // a 55S65M CIGAR asserts a 120-base SEQ and the record fails to write.
    if let Some(hard_leading) = supplementary {
        use noodles::sam::alignment::record::cigar::{Op, op::Kind};
        let idx = if hard_leading { 0 } else { cigar_ops.len() - 1 };
        if let Some(op) = cigar_ops.get(idx).copied()
            && op.kind() == Kind::SoftClip
        {
            let n = op.len();
            cigar_ops[idx] = Op::new(Kind::HardClip, n);
            if hard_leading {
                seq_bytes.drain(..n.min(seq_bytes.len()));
            } else {
                let keep = seq_bytes.len().saturating_sub(n);
                seq_bytes.truncate(keep);
            }
        }
    }

    *record.cigar_mut() = cigar_ops.iter().copied().collect();
    *record.sequence_mut() = Sequence::from(seq_bytes);
    // QUAL is not carried for chimeric segments.
    *record.quality_scores_mut() = QualityScores::default();

    let data = record.data_mut();
    data.insert(Tag::new(b'S', b'A'), Value::String(BString::from(sa_tag)));
    data.insert(Tag::new(b'N', b'M'), Value::from(seg.n_mismatch as i32));
    data.insert(Tag::ALIGNMENT_SCORE, Value::from(seg.score));

    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chimeric::segment::ExonSpan;
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
            first_exon: ExonSpan {
                genome_start: 133_738_300,
                genome_end: 133_738_363,
                read_start: 0,
                read_end: 63,
            },
            last_exon: ExonSpan {
                genome_start: 133_738_300,
                genome_end: 133_738_363,
                read_start: 0,
                read_end: 63,
            },
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
            first_exon: ExonSpan {
                genome_start: 23_632_600,
                genome_end: 23_632_637,
                read_start: 63,
                read_end: 100,
            },
            last_exon: ExonSpan {
                genome_start: 23_632_600,
                genome_end: 23_632_637,
                read_start: 63,
                read_end: 100,
            },
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
        // Distinct, non-zero padded starts. With both at 0 the absolute and
        // per-chromosome coordinates coincide, which is what let the missing
        // chrStart subtraction go unnoticed; distinct values also catch using
        // one chromosome's offset for both.
        let chr_starts = vec![1_000_000u64, 20_000_000u64];

        writer
            .write_alignment(&alignment, &chr_names, &chr_starts, "READ_001")
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
        assert_eq!(fields[1], "132738364"); // donor breakpoint (chimJ0), per-chr
        assert_eq!(fields[2], "+"); // donor strand
        assert_eq!(fields[3], "chr22"); // acceptor chr
        assert_eq!(fields[4], "3632600"); // acceptor breakpoint (chimJ1), per-chr
        assert_eq!(fields[5], "+"); // acceptor strand
        assert_eq!(fields[6], "1"); // junction type
        assert_eq!(fields[7], "0"); // repeat donor
        assert_eq!(fields[8], "0"); // repeat acceptor
        assert_eq!(fields[9], "READ_001"); // read name
        assert_eq!(fields[10], "132738301"); // donor start (1-based, per-chr)
        assert_eq!(fields[11], "63M"); // donor CIGAR
        assert_eq!(fields[12], "3632601"); // acceptor start (1-based, per-chr)
        assert_eq!(fields[13], "37M"); // acceptor CIGAR
    }

    fn chr_names() -> Vec<String> {
        vec!["chr9".to_string(), "chr22".to_string()]
    }

    fn chr_starts() -> Vec<u64> {
        vec![1_000_000, 20_000_000]
    }

    /// A minimal two-segment chimera for output-format tests.
    fn mock_alignment() -> ChimericAlignment {
        let seg = |chr_idx: usize, gs: u64, ge: u64, rs: usize, re: usize| ChimericSegment {
            chr_idx,
            genome_start: gs,
            genome_end: ge,
            is_reverse: false,
            read_start: rs,
            read_end: re,
            cigar: vec![cigar::Op::new(cigar::op::Kind::Match, re - rs)],
            score: (re - rs) as i32,
            n_mismatch: 0,
            first_exon: ExonSpan {
                genome_start: gs,
                genome_end: ge,
                read_start: rs,
                read_end: re,
            },
            last_exon: ExonSpan {
                genome_start: gs,
                genome_end: ge,
                read_start: rs,
                read_end: re,
            },
        };
        ChimericAlignment::new(
            seg(0, 1_100_000, 1_100_060, 0, 60),
            seg(1, 20_100_000, 20_100_060, 60, 120),
            1,
            0,
            0,
            vec![0u8; 120],
            "READ_001".to_string(),
        )
    }

    /// `--chimMultimapNmax` selects STAR's wider file: a 21-name header and six
    /// extra per-record columns. Without it the file stays at the classic 14
    /// with no header. STAR-Fusion parses on that header, so the shape matters
    /// as much as the values.
    #[test]
    fn multimap_mode_writes_stars_header_and_six_extra_columns() {
        use crate::chimeric::MultimapInfo;
        let dir = tempdir().unwrap();

        let build = |multimap: bool| -> String {
            let prefix = format!("{}/mm{}_", dir.path().display(), u8::from(multimap));
            let mut w = ChimericJunctionWriter::new_with_multimap(&prefix, multimap).unwrap();
            let mut aln = mock_alignment();
            if multimap {
                aln = aln.with_multimap(MultimapInfo {
                    chim_n: 3,
                    max_possible_score: 120,
                    max_non_chim_score: 61,
                    chim_score: 116,
                    best_chim_score: 118,
                    pe_merged: false,
                });
            }
            w.write_alignment(&aln, &chr_names(), &chr_starts(), "READ_001")
                .unwrap();
            w.flush().unwrap();
            std::fs::read_to_string(format!("{prefix}Chimeric.out.junction")).unwrap()
        };

        // Off: no header, 14 columns.
        let plain = build(false);
        let plain_lines: Vec<&str> = plain.lines().collect();
        assert_eq!(plain_lines.len(), 1, "no header without the flag");
        assert_eq!(plain_lines[0].split('\t').count(), 14);

        // On: STAR's header verbatim, then 20 columns.
        let multi = build(true);
        let lines: Vec<&str> = multi.lines().collect();
        assert_eq!(lines.len(), 2, "header plus one record");
        assert_eq!(
            lines[0],
            "chr_donorA\tbrkpt_donorA\tstrand_donorA\tchr_acceptorB\tbrkpt_acceptorB\t\
             strand_acceptorB\tjunction_type\trepeat_left_lenA\trepeat_right_lenB\tread_name\t\
             start_alnA\tcigar_alnA\tstart_alnB\tcigar_alnB\tnum_chim_aln\tmax_poss_aln_score\t\
             non_chim_aln_score\tthis_chim_aln_score\tbestall_chim_aln_score\tPEmerged_bool\treadgrp",
            "header must match STAR's byte for byte"
        );
        let f: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(f.len(), 20, "STAR emits 20 columns; readgrp only with a RG");
        assert_eq!(&f[14..], &["3", "120", "61", "116", "118", "0"]);
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
            first_exon: ExonSpan {
                genome_start: 1000,
                genome_end: 1050,
                read_start: 0,
                read_end: 50,
            },
            last_exon: ExonSpan {
                genome_start: 1000,
                genome_end: 1050,
                read_start: 0,
                read_end: 50,
            },
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
            first_exon: ExonSpan {
                genome_start: 2000,
                genome_end: 2050,
                read_start: 50,
                read_end: 100,
            },
            last_exon: ExonSpan {
                genome_start: 2000,
                genome_end: 2050,
                read_start: 50,
                read_end: 100,
            },
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
        let chr_starts = vec![500u64];

        writer
            .write_alignment(&alignment, &chr_names, &chr_starts, "READ_002")
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
        // chrStart 500 comes off every coordinate, on both strands.
        assert_eq!(fields[1], "551"); // donor breakpoint (forward: genome_end + 1)
        assert_eq!(fields[4], "1551"); // acceptor breakpoint (reverse: genome_end + 1)
        assert_eq!(fields[10], "501"); // donor start
        assert_eq!(fields[12], "1501"); // acceptor start
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
            // Real segments span the whole read: the part belonging to the
            // other segment is soft-clipped, not absent.
            cigar: vec![Op::new(Kind::Match, 63)],
            score: 63,
            n_mismatch: 0,
            first_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
            last_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
        };
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: false,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::SoftClip, 63), Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 1,
            first_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
            last_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
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
        let records = build_within_bam_records(&alignment, &genome, 255, false).unwrap();

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
            first_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
            last_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
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
            first_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
            last_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
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
        let records = build_within_bam_records(&alignment, &genome, 255, false).unwrap();

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
            first_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
            last_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
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
            first_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
            last_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
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
        let records = build_within_bam_records(&alignment, &genome, 255, false).unwrap();

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

    /// A segment whose CIGAR is a bare match block must still produce a record
    /// whose CIGAR spans the read.
    ///
    /// Detection sometimes hands us a segment with no clip ops at all (a bare
    /// `59M` on a 120-base read). STAR always emits the surrounding soft clips
    /// (`trimL1`/`trimR1`, `ReadAlign_alignBAM.cpp:235,269`); without them the
    /// record claims fewer query bases than SEQ carries and BAM writing fails
    /// with "read length-sequence length mismatch".
    #[test]
    fn test_within_bam_pads_a_bare_cigar_to_span_the_read() {
        use cigar::op::{Kind, Op};
        let bare = |chr_idx: usize, gs: u64, rs: usize, re: usize, rev: bool| ChimericSegment {
            chr_idx,
            genome_start: gs,
            genome_end: gs + (re - rs) as u64,
            is_reverse: rev,
            read_start: rs,
            read_end: re,
            // No clips, deliberately.
            cigar: vec![Op::new(Kind::Match, re - rs)],
            score: (re - rs) as i32,
            n_mismatch: 0,
            first_exon: ExonSpan {
                genome_start: gs,
                genome_end: gs + (re - rs) as u64,
                read_start: rs,
                read_end: re,
            },
            last_exon: ExonSpan {
                genome_start: gs,
                genome_end: gs + (re - rs) as u64,
                read_start: rs,
                read_end: re,
            },
        };
        // Reverse acceptor, mirroring the fixture read that first tripped this.
        let alignment = ChimericAlignment::new(
            bare(0, 100, 0, 61, false),
            bare(1, 600, 61, 120, true),
            0,
            0,
            0,
            vec![0u8; 120],
            "READ_BARE".to_string(),
        );
        let records =
            build_within_bam_records(&alignment, &make_genome_2chr(), 255, false).unwrap();

        for rec in &records {
            let q: usize = rec
                .cigar()
                .as_ref()
                .iter()
                .filter(|op| op.kind().consumes_read())
                .map(|op| op.len())
                .sum();
            assert_eq!(
                q,
                rec.sequence().len(),
                "CIGAR query length must match SEQ: {}",
                crate::align::transcript::cigar_to_string(rec.cigar().as_ref())
            );
        }
    }

    #[test]
    fn test_within_bam_supplementary_is_hard_clipped_with_trimmed_seq() {
        use cigar::op::{Kind, Op};
        let donor = ChimericSegment {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 163,
            is_reverse: false,
            read_start: 0,
            read_end: 63,
            // Real segments span the whole read: the part belonging to the
            // other segment is soft-clipped, not absent.
            cigar: vec![Op::new(Kind::Match, 63), Op::new(Kind::SoftClip, 37)],
            score: 63,
            n_mismatch: 0,
            first_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
            last_exon: ExonSpan {
                genome_start: 100,
                genome_end: 163,
                read_start: 0,
                read_end: 63,
            },
        };
        let acceptor = ChimericSegment {
            chr_idx: 1,
            genome_start: 600,
            genome_end: 637,
            is_reverse: false,
            read_start: 63,
            read_end: 100,
            cigar: vec![Op::new(Kind::SoftClip, 63), Op::new(Kind::Match, 37)],
            score: 37,
            n_mismatch: 0,
            first_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
            last_exon: ExonSpan {
                genome_start: 600,
                genome_end: 637,
                read_start: 63,
                read_end: 100,
            },
        };
        let read_seq = vec![0u8; 100]; // 100 A bases
        let alignment =
            ChimericAlignment::new(donor, acceptor, 0, 0, 0, read_seq, "READ_001".to_string());
        let genome = make_genome_2chr();
        let records = build_within_bam_records(&alignment, &genome, 255, false).unwrap();

        // The representative record carries the whole read against a CIGAR that
        // spans it.
        assert_eq!(records[0].sequence().len(), 100, "donor keeps the full SEQ");
        assert_eq!(
            crate::align::transcript::cigar_to_string(records[0].cigar().as_ref()),
            "63M37S"
        );

        // The supplementary hard-clips the end covering the donor and drops
        // those bases from SEQ, as STAR does (chimericBAMoutput.cpp:55,
        // alignBAM.cpp:503-511). Leaving them soft-clipped with an empty SEQ is
        // not representable: soft clips consume the query, so the record would
        // claim 100 bases and carry none.
        assert_eq!(
            crate::align::transcript::cigar_to_string(records[1].cigar().as_ref()),
            "63H37M",
            "supplementary hard-clips the other segment's bases"
        );
        assert_eq!(
            records[1].sequence().len(),
            37,
            "SEQ excludes hard-clipped bases"
        );

        // The invariant the SAM writer enforces: CIGAR query length == SEQ length.
        for rec in &records {
            let q: usize = rec
                .cigar()
                .as_ref()
                .iter()
                .filter(|op| op.kind().consumes_read())
                .map(|op| op.len())
                .sum();
            assert_eq!(q, rec.sequence().len(), "CIGAR query length must match SEQ");
        }
    }
}
