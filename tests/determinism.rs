//! Output must not depend on how many threads produced it.
//!
//! Junction counts live in a `DashMap`, whose iteration order varies with
//! hashing and with concurrent insertion. Every path that emits an order sorts
//! first, and these tests are what keeps that true: they align the same reads
//! at one thread and at eight and compare the output files byte for byte
//! (issue #210).

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// Two exons with a GT-AG intron between them, so reads spanning the junction
/// produce SJ.out.tab rows rather than an empty file.
const EXON1: (usize, usize) = (2_000, 2_300);
const INTRON: (usize, usize) = (2_300, 2_800);
const EXON2: (usize, usize) = (2_800, 3_100);
/// A second junction close to the first, to give the neighbour-distance filter
/// something to compute and the sort something to order.
const EXON3: (usize, usize) = (6_000, 6_200);
const INTRON2: (usize, usize) = (6_200, 6_600);
const EXON4: (usize, usize) = (6_600, 6_900);

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

fn build_genome() -> Vec<u8> {
    let mut genome = lcg_seq(88888, 20_000);
    for (start, end) in [INTRON, INTRON2] {
        genome[start] = b'G';
        genome[start + 1] = b'T';
        genome[end - 2] = b'A';
        genome[end - 1] = b'G';
    }
    genome
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

/// Reads spanning both junctions plus unspliced filler, enough of them that
/// several threads each hold some.
fn write_fastq(dir: &Path, genome: &[u8]) -> PathBuf {
    let path = dir.join("reads.fq");
    let mut f = fs::File::create(&path).unwrap();
    let mut n = 0usize;

    let spliced = |f: &mut fs::File, e_end: usize, i_end: usize, overhang: usize, n: &mut usize| {
        let left = &genome[e_end - overhang..e_end];
        let right = &genome[i_end..i_end + (100 - overhang)];
        let mut seq = left.to_vec();
        seq.extend_from_slice(right);
        writeln!(f, "@sj{n}").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f, "\n+\n{}", "I".repeat(seq.len())).unwrap();
        *n += 1;
    };

    for overhang in 30..70 {
        spliced(&mut f, INTRON.0, INTRON.1, overhang, &mut n);
        spliced(&mut f, INTRON2.0, INTRON2.1, overhang, &mut n);
    }
    for i in 0..200usize {
        let start = 8_000 + i * 40;
        writeln!(f, "@u{i}").unwrap();
        f.write_all(&genome[start..start + 100]).unwrap();
        writeln!(f, "\n+\n{}", "I".repeat(100)).unwrap();
    }
    let _ = (EXON1, EXON2, EXON3, EXON4);
    path
}

/// Align with `threads` threads and return the bytes of every output file that
/// carries an order.
fn align(
    genome_dir: &Path,
    fastq: &Path,
    prefix: &str,
    threads: &str,
    extra: &[&str],
) -> Vec<(String, Vec<u8>)> {
    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.args([
        "--runMode",
        "alignReads",
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--readFilesIn",
        fastq.to_str().unwrap(),
        "--runThreadN",
        threads,
        "--outFileNamePrefix",
        prefix,
    ]);
    cmd.args(extra);
    cmd.assert().success();

    ["SJ.out.tab", "Aligned.out.sam"]
        .iter()
        .map(|name| {
            let bytes = fs::read(format!("{prefix}{name}")).unwrap_or_default();
            // Drop the @PG header line: it records the command line, which
            // differs by the thread count itself.
            let filtered: Vec<u8> = if *name == "Aligned.out.sam" {
                String::from_utf8_lossy(&bytes)
                    .lines()
                    .filter(|l| !l.starts_with("@PG"))
                    .collect::<Vec<_>>()
                    .join("\n")
                    .into_bytes()
            } else {
                bytes
            };
            ((*name).to_string(), filtered)
        })
        .collect()
}

#[test]
fn output_is_identical_at_one_and_eight_threads() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = build_genome();
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);
    let fq = write_fastq(root, &genome);

    let one = align(
        &genome_dir,
        &fq,
        &format!("{}/t1_", root.display()),
        "1",
        &[],
    );
    let eight = align(
        &genome_dir,
        &fq,
        &format!("{}/t8_", root.display()),
        "8",
        &[],
    );

    // The fixture has to actually produce junctions, or this test would pass
    // on two empty files.
    let sj = &one[0].1;
    let rows = String::from_utf8_lossy(sj).lines().count();
    println!("fixture produced {rows} SJ.out.tab rows");
    assert!(
        rows >= 2,
        "the fixture has to produce at least two junction rows for an order to \
         exist at all, got {rows}"
    );

    for ((name, a), (_, b)) in one.iter().zip(eight.iter()) {
        assert_eq!(
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(b),
            "{name} differs between 1 and 8 threads"
        );
    }
}

#[test]
fn output_is_identical_across_two_runs_at_eight_threads() {
    // Guards against the pair above agreeing only because both runs happened
    // to hit the same map order.
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = build_genome();
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);
    let fq = write_fastq(root, &genome);

    let a = align(
        &genome_dir,
        &fq,
        &format!("{}/a_", root.display()),
        "8",
        &[],
    );
    let b = align(
        &genome_dir,
        &fq,
        &format!("{}/b_", root.display()),
        "8",
        &[],
    );

    for ((name, x), (_, y)) in a.iter().zip(b.iter()) {
        assert_eq!(
            String::from_utf8_lossy(x),
            String::from_utf8_lossy(y),
            "{name} differs between two 8-thread runs"
        );
    }
}

#[test]
fn two_pass_output_is_identical_at_one_and_eight_threads() {
    // Two-pass feeds pass 1's junctions back into the alignment, so any order
    // escaping the junction map would show up here rather than in a single
    // pass.
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = build_genome();
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);
    let fq = write_fastq(root, &genome);

    let extra = ["--twopassMode", "Basic"];
    let one = align(
        &genome_dir,
        &fq,
        &format!("{}/tp1_", root.display()),
        "1",
        &extra,
    );
    let eight = align(
        &genome_dir,
        &fq,
        &format!("{}/tp8_", root.display()),
        "8",
        &extra,
    );

    for ((name, a), (_, b)) in one.iter().zip(eight.iter()) {
        assert_eq!(
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(b),
            "two-pass {name} differs between 1 and 8 threads"
        );
    }

    // The pass-1 junction file is written from the same map and is equally
    // order-sensitive.
    let p1 = fs::read(format!("{}/tp1__STARpass1/SJ.out.tab", root.display())).unwrap();
    let p8 = fs::read(format!("{}/tp8__STARpass1/SJ.out.tab", root.display())).unwrap();
    assert_eq!(
        String::from_utf8_lossy(&p1),
        String::from_utf8_lossy(&p8),
        "pass-1 SJ.out.tab differs between 1 and 8 threads"
    );
}
