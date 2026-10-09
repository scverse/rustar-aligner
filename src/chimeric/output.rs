// Chimeric.out.junction file writer (WithinBAM records are built in `bam.rs`)

use crate::chimeric::segment::ChimericAlignment;
use crate::error::Error;
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
}
