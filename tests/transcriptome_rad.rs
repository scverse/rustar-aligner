//! Integration test for `--quantTranscriptomeRAD Yes`.
//!
//! `Aligned.toTranscriptome.out.rad` must hold the same alignments as
//! `Aligned.toTranscriptome.out.bam`, which is what makes `salmon quant --rad`
//! agree with `salmon quant -a` on the BAM, and asking for it must not change
//! the BAM.

use assert_cmd::cargo::cargo_bin_cmd;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::RecordBuf;
use rustar_aligner::io::rad::{
    FRAG_MAPPED_PAIR, FRAG_SINGLE, RadFile, RadFragment, RadHit, RadTagValue, read_rad,
};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const READ_LEN: usize = 50;
const FRAGMENT_LEN: usize = 180;

fn generate_genome_seq(seed: u32, length: usize) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    let mut state = seed;
    let mut seq = String::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

struct Fixture {
    dir: TempDir,
    genome_dir: PathBuf,
    mate1: PathBuf,
    mate2: PathBuf,
}

const N_FRAGMENTS: usize = 30;
/// One more read than maps, to tell `num_processed` from the fragments written.
const N_READS: usize = N_FRAGMENTS + 1;

/// T1 (+, spliced) and T1b (+, its first exon) share a locus, so reads there have
/// two placements; T2 (-) gives the other strand.
fn fixture() -> Fixture {
    let dir = TempDir::new().unwrap();
    let seq = generate_genome_seq(98765, 4000);

    let fasta = dir.path().join("genome.fa");
    fs::write(&fasta, format!(">chr1\n{seq}\n")).unwrap();

    let gtf = dir.path().join("annotations.gtf");
    let mut f = fs::File::create(&gtf).unwrap();
    for (start, end, strand, gene, tr) in [
        (101, 500, '+', "G1", "T1"),
        (801, 1100, '+', "G1", "T1"),
        (101, 500, '+', "G1", "T1b"),
        (2001, 2500, '-', "G2", "T2"),
    ] {
        writeln!(
            f,
            "chr1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{tr}\";"
        )
        .unwrap();
    }

    let mate1 = dir.path().join("r1.fq");
    let mate2 = dir.path().join("r2.fq");
    let mut w1 = fs::File::create(&mate1).unwrap();
    let mut w2 = fs::File::create(&mate2).unwrap();
    let qual = "I".repeat(READ_LEN);
    for i in 0..N_FRAGMENTS {
        let off = i * 7 % 200;
        // Both read-1 strands on both transcripts: an unstranded library.
        let (start, reverse) = match i % 4 {
            0 => (110 + off, false),
            1 => (2010 + off, true),
            2 => (110 + off, true),
            _ => (2010 + off, false),
        };
        let frag = &seq[start..start + FRAGMENT_LEN];
        let left = frag[..READ_LEN].to_string();
        let right = revcomp(&frag[FRAGMENT_LEN - READ_LEN..]);
        let (r1, r2) = if reverse {
            (right, left)
        } else {
            (left, right)
        };
        writeln!(w1, "@read{i}\n{r1}\n+\n{qual}").unwrap();
        writeln!(w2, "@read{i}\n{r2}\n+\n{qual}").unwrap();
    }
    let noise = "ACGT".repeat(READ_LEN / 4 + 1);
    writeln!(w1, "@noise\n{0}\n+\n{qual}", &noise[..READ_LEN]).unwrap();
    writeln!(w2, "@noise\n{0}\n+\n{qual}", &noise[..READ_LEN]).unwrap();

    let genome_dir = dir.path().join("genome");
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--sjdbOverhang",
            "49",
            "--genomeSAindexNbases",
            "5",
            "--outFileNamePrefix",
            genome_dir.join("run_").to_str().unwrap(),
        ])
        .assert()
        .success();
    Fixture {
        dir,
        genome_dir,
        mate1,
        mate2,
    }
}

fn run(fx: &Fixture, name: &str, paired: bool, opts: &[&str]) -> PathBuf {
    let out = fx.dir.path().join(name);
    fs::create_dir_all(&out).unwrap();
    let prefix = format!("{}/", out.display());
    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.args(["--genomeDir", fx.genome_dir.to_str().unwrap()])
        .arg("--readFilesIn")
        .arg(&fx.mate1);
    if paired {
        cmd.arg(&fx.mate2);
    }
    cmd.args([
        "--runThreadN",
        "1",
        "--outSAMtype",
        "None",
        "--outFileNamePrefix",
        prefix.as_str(),
    ])
    .args(opts)
    .assert()
    .success();
    out
}

fn read_bam(path: &Path) -> Vec<RecordBuf> {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .map(|r| r.unwrap())
        .collect::<Vec<_>>()
}

fn ref_span(r: &RecordBuf) -> u16 {
    r.cigar()
        .as_ref()
        .iter()
        .filter(|op| matches!(op.kind(), Kind::Match | Kind::Deletion))
        .map(|op| op.len() as u16)
        .sum()
}

fn pos0(r: &RecordBuf) -> u32 {
    usize::from(r.alignment_start().unwrap()) as u32 - 1
}

/// The placements each read's BAM records imply, as salmon derives them
/// (`salmon-align/src/bam_rad.rs`), with the score left out: the BAM has none.
fn fragments_from_bam(records: &[RecordBuf], paired: bool) -> Vec<RadFragment> {
    let mut out: Vec<(Vec<u8>, RadFragment)> = Vec::new();
    for group in records.chunks(if paired { 2 } else { 1 }) {
        let name: Vec<u8> = group[0].name().unwrap().to_vec();
        let hit = if paired {
            let r1 = group.iter().find(|r| r.flags().is_first_segment()).unwrap();
            let r2 = group.iter().find(|r| r.flags().is_last_segment()).unwrap();
            assert_eq!(r1.reference_sequence_id(), r2.reference_sequence_id());
            RadHit {
                tid: r1.reference_sequence_id().unwrap() as u32,
                is_fw: !r1.flags().is_reverse_complemented(),
                mate_fw: !r2.flags().is_reverse_complemented(),
                pos: pos0(r1).min(pos0(r2)),
                frag_len: r1.template_length().unsigned_abs() as u16,
                score: 0,
            }
        } else {
            let r = &group[0];
            RadHit {
                tid: r.reference_sequence_id().unwrap() as u32,
                is_fw: !r.flags().is_reverse_complemented(),
                mate_fw: false,
                pos: pos0(r),
                frag_len: ref_span(r),
                score: 0,
            }
        };
        match out.last_mut() {
            Some((n, f)) if *n == name => f.hits.push(hit),
            _ => out.push((
                name,
                RadFragment {
                    frag_map_type: if paired {
                        FRAG_MAPPED_PAIR
                    } else {
                        FRAG_SINGLE
                    },
                    hits: vec![hit],
                },
            )),
        }
    }
    out.into_iter().map(|(_, f)| f).collect()
}

fn without_scores(rad: &RadFile) -> Vec<RadFragment> {
    rad.fragments
        .iter()
        .map(|f| RadFragment {
            frag_map_type: f.frag_map_type,
            hits: f.hits.iter().map(|h| RadHit { score: 0, ..*h }).collect(),
        })
        .collect()
}

fn assert_rad_matches_bam(paired: bool) {
    let fx = fixture();
    let rad_flag = ["--quantTranscriptomeRAD", "Yes"];
    let sam_flag = ["--quantMode", "TranscriptomeSAM"];
    let both = run(
        &fx,
        "both",
        paired,
        &[&sam_flag[..], &rad_flag[..]].concat(),
    );
    let bam_only = run(&fx, "bam", paired, &sam_flag);
    let rad_only = run(&fx, "rad", paired, &rad_flag);
    let bam_name = "Aligned.toTranscriptome.out.bam";
    let rad_name = "Aligned.toTranscriptome.out.rad";

    // The RAD option leaves the BAM alone, and the BAM option does not touch the RAD.
    assert_eq!(
        fs::read(both.join(bam_name)).unwrap(),
        fs::read(bam_only.join(bam_name)).unwrap(),
        "the transcriptome BAM must not depend on --quantTranscriptomeRAD"
    );
    assert!(!rad_only.join(bam_name).exists(), "RAD alone writes no BAM");
    assert!(!bam_only.join(rad_name).exists(), "the RAD is opt-in");
    let rad = read_rad(&both.join(rad_name)).unwrap();
    let alone = read_rad(&rad_only.join(rad_name)).unwrap();
    assert_eq!(rad.fragments, alone.fragments);

    // Same alignments as the BAM, in the same order.
    let records = read_bam(&both.join(bam_name));
    assert!(!records.is_empty());
    let expected = fragments_from_bam(&records, paired);
    assert_eq!(without_scores(&rad), expected);

    // The fixture has to exercise what it is there for.
    let hits = || rad.fragments.iter().flat_map(|f| &f.hits);
    assert!(rad.fragments.iter().any(|f| f.hits.len() > 1));
    assert!(hits().any(|h| h.is_fw) && hits().any(|h| !h.is_fw));
    assert!(hits().all(|h| h.score > 0), "the score of the alignment");

    // Header: the transcript table of the BAM, and every input read counted.
    // STAR's transcript order: by start, then by end.
    assert_eq!(rad.ref_names, ["T1b", "T1", "T2"]);
    assert_eq!(
        rad.file_tag("ref_lengths"),
        Some(&RadTagValue::ArrayU32(vec![400, 700, 500]))
    );
    assert_eq!(rad.is_paired, paired);
    assert_eq!(
        rad.file_tag("num_processed"),
        Some(&RadTagValue::U64(N_READS as u64))
    );
    assert!(
        rad.fragments.len() < N_READS,
        "the unmapped read is left out"
    );
    let Some(RadTagValue::String(program)) = rad.file_tag("source_programs") else {
        panic!("no source_programs tag");
    };
    assert!(program.starts_with("@PG\tID:rustar-aligner"));
    assert!(program.contains("--quantTranscriptomeRAD"));
}

#[test]
fn single_end_rad_matches_the_transcriptome_bam() {
    assert_rad_matches_bam(false);
}

#[test]
fn paired_end_rad_matches_the_transcriptome_bam() {
    assert_rad_matches_bam(true);
}

#[test]
fn pair_score_is_the_whole_pair_and_single_end_the_read() {
    let fx = fixture();
    let flag = ["--quantTranscriptomeRAD", "Yes"];
    let name = "Aligned.toTranscriptome.out.rad";
    let se = read_rad(&run(&fx, "se", false, &flag).join(name)).unwrap();
    let pe = read_rad(&run(&fx, "pe", true, &flag).join(name)).unwrap();
    let max = |r: &RadFile| {
        r.fragments
            .iter()
            .flat_map(|f| &f.hits)
            .map(|h| h.score)
            .max()
            .unwrap()
    };
    // STAR's AS: read length - 1 for a clean read, sum of the mates' lengths - 2
    // for a clean pair (before the genomic length penalty).
    assert!(max(&se) <= READ_LEN as i32 && max(&se) > READ_LEN as i32 - 3);
    assert!(max(&pe) <= 2 * READ_LEN as i32 && max(&pe) > 2 * READ_LEN as i32 - 6);
}

#[test]
fn rad_without_gtf_aware_index_fails_clearly() {
    let dir = TempDir::new().unwrap();
    let seq = generate_genome_seq(1, 2000);
    let fasta = dir.path().join("g.fa");
    fs::write(&fasta, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.path().join("r.fq");
    fs::write(
        &fq,
        format!("@r\n{}\n+\n{}\n", &seq[100..150], "I".repeat(50)),
    )
    .unwrap();
    let genome_dir = dir.path().join("genome");
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--genomeSAindexNbases",
            "5",
            "--outFileNamePrefix",
            genome_dir.join("run_").to_str().unwrap(),
        ])
        .assert()
        .success();
    let prefix = format!("{}/", dir.path().display());
    let output = cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", genome_dir.to_str().unwrap()])
        .args(["--readFilesIn", fq.to_str().unwrap()])
        .args(["--outFileNamePrefix", prefix.as_str()])
        .args(["--quantTranscriptomeRAD", "Yes"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--quantTranscriptomeRAD"));
}
