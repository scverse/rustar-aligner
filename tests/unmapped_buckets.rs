//! `Log.final.out` has to tell "too short" and "other" apart.
//!
//! STAR splits them in `ReadAlign_mappedFilter.cpp`: a read with no good
//! window at all is `unmappedOther`, and only a read that *had* a window whose
//! best transcript failed the score or length thresholds is `unmappedShort`.
//! Folding both into "too short" leaves the "other" bucket permanently zero,
//! which is what issue #48 measured on paired-end data.

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const READ_LEN: usize = 100;

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let bases: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut seq = Vec::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

fn write_fasta(dir: &Path, genome: &[u8]) -> PathBuf {
    let path = dir.join("genome.fa");
    let mut f = fs::File::create(&path).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(genome).unwrap();
    writeln!(f).unwrap();
    path
}

fn build_index(fasta: &Path, genome_dir: &Path) {
    fs::create_dir_all(genome_dir).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--genomeSAindexNbases",
            "7",
        ])
        .assert()
        .success();
}

fn parse_final_log(path: &Path) -> Vec<(String, u64)> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| {
            let (label, value) = l.split_once('|')?;
            let value = value.trim();
            if value.ends_with('%') {
                return None;
            }
            Some((label.trim().to_string(), value.parse::<u64>().ok()?))
        })
        .collect()
}

fn get(rows: &[(String, u64)], label: &str) -> u64 {
    rows.iter()
        .find(|(l, _)| l == label)
        .map(|(_, v)| *v)
        .unwrap_or_else(|| panic!("{label} missing from Log.final.out"))
}

/// Paired reads of three kinds: mapping cleanly, matching the genome over a
/// stretch too short to pass the score filter, and matching nothing at all.
fn write_reads(dir: &Path, genome: &[u8]) -> (PathBuf, PathBuf) {
    let p1 = dir.join("r1.fq");
    let p2 = dir.join("r2.fq");
    let mut f1 = fs::File::create(&p1).unwrap();
    let mut f2 = fs::File::create(&p2).unwrap();

    let mut emit = |name: &str, s1: &[u8], s2: &[u8]| {
        writeln!(f1, "@{name}").unwrap();
        f1.write_all(s1).unwrap();
        writeln!(f1, "\n+\n{}", "I".repeat(s1.len())).unwrap();
        writeln!(f2, "@{name}").unwrap();
        f2.write_all(s2).unwrap();
        writeln!(f2, "\n+\n{}", "I".repeat(s2.len())).unwrap();
    };

    // 20 pairs that map: mate2 is the reverse complement of a downstream slice.
    let rc = |s: &[u8]| -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    };
    for i in 0..20usize {
        let start = 2_000 + i * 300;
        let m1 = &genome[start..start + READ_LEN];
        let m2 = rc(&genome[start + 150..start + 150 + READ_LEN]);
        emit(&format!("ok{i}"), m1, &m2);
    }

    // 20 pairs that seed but cannot clear the score threshold: 30 genomic
    // bases then 70 bases belonging to no chromosome. A window exists, so
    // STAR calls these "too short", not "other".
    for i in 0..20usize {
        let start = 9_000 + i * 137;
        let mut m1 = genome[start..start + 30].to_vec();
        m1.extend_from_slice(&lcg_seq(4_242 + i as u32, READ_LEN - 30));
        let mut m2 = genome[start + 200..start + 230].to_vec();
        m2.extend_from_slice(&lcg_seq(9_191 + i as u32, READ_LEN - 30));
        emit(&format!("short{i}"), &m1, &m2);
    }

    // 20 pairs with no genomic seed at all, which is the "other" bucket.
    for i in 0..20usize {
        let m1 = lcg_seq(700_000 + i as u32, READ_LEN);
        let m2 = lcg_seq(900_000 + i as u32, READ_LEN);
        emit(&format!("other{i}"), &m1, &m2);
    }

    (p1, p2)
}

#[test]
fn paired_end_unmapped_reads_split_between_too_short_and_other() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = lcg_seq(88888, 20_000);
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);
    let (r1, r2) = write_reads(root, &genome);

    let prefix = format!("{}/pe_", root.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let rows = parse_final_log(Path::new(&format!("{prefix}Log.final.out")));
    let input = get(&rows, "Number of input reads");
    let unique = get(&rows, "Uniquely mapped reads number");
    let too_short = get(&rows, "Number of reads unmapped: too short");
    let other = get(&rows, "Number of reads unmapped: other");

    assert_eq!(input, 60, "fixture size");
    assert!(unique > 0, "the mapping pairs must map: {rows:?}");
    assert!(
        other > 0,
        "reads with no genomic seed belong in `other`, which was permanently \
         zero before this split: too_short={too_short}, other={other}"
    );
    assert!(
        too_short > 0,
        "reads that seed but fail the score threshold belong in `too short`: \
         too_short={too_short}, other={other}"
    );
    // Every unmapped read lands in exactly one bucket.
    let multi = get(&rows, "Number of reads mapped to multiple loci");
    let too_many = get(&rows, "Number of reads mapped to too many loci");
    let mismatches = get(&rows, "Number of reads unmapped: too many mismatches");
    assert_eq!(
        unique + multi + too_many + mismatches + too_short + other,
        input,
        "the buckets have to add up to the input: {rows:?}"
    );
}

#[test]
fn single_end_unmapped_reads_split_between_too_short_and_other() {
    // The single-end path already distinguished the two; this keeps it that
    // way, and makes the pair of tests say which path is being checked.
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = lcg_seq(88888, 20_000);
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);
    let (r1, _r2) = write_reads(root, &genome);

    let prefix = format!("{}/se_", root.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            r1.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let rows = parse_final_log(Path::new(&format!("{prefix}Log.final.out")));
    assert!(
        get(&rows, "Number of reads unmapped: other") > 0,
        "single-end `other` bucket: {rows:?}"
    );
    assert!(
        get(&rows, "Number of reads unmapped: too short") > 0,
        "single-end `too short` bucket: {rows:?}"
    );
}
