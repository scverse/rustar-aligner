//! `--readFilesCommand` takes what STAR takes: a shell command line in one
//! value or several, or `-` for no command. Every form must give the output of
//! the plain file.
#![cfg(unix)]

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Output;
use std::sync::OnceLock;
use std::time::Duration;
use tempfile::TempDir;

const TIMEOUT: Duration = Duration::from_secs(60);
const READ_LEN: usize = 50;

fn genome() -> &'static [u8] {
    static GENOME: OnceLock<Vec<u8>> = OnceLock::new();
    GENOME.get_or_init(|| {
        let mut state = 88888u32;
        (0..20_000)
            .map(|_| {
                state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((state >> 16) & 3) as usize]
            })
            .collect()
    })
}

fn index() -> &'static (TempDir, PathBuf) {
    static INDEX: OnceLock<(TempDir, PathBuf)> = OnceLock::new();
    INDEX.get_or_init(|| {
        let dir = TempDir::new().unwrap();
        let mut fasta = b">chr1\n".to_vec();
        fasta.extend_from_slice(genome());
        fasta.push(b'\n');
        fs::write(dir.path().join("genome.fa"), fasta).unwrap();
        let genome_dir = dir.path().join("idx");
        fs::create_dir(&genome_dir).unwrap();
        cargo_bin_cmd!("rustar-aligner")
            .current_dir(dir.path())
            .args(["--runMode", "genomeGenerate", "--genomeDir", "idx"])
            .args(["--genomeFastaFiles", "genome.fa"])
            .args(["--genomeSAindexNbases", "7"])
            .assert()
            .success();
        (dir, genome_dir)
    })
}

fn fastq() -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..40 {
        let p = 1_000 + 7 * i;
        writeln!(out, "@r{i}").unwrap();
        out.extend_from_slice(&genome()[p..p + READ_LEN]);
        write!(out, "\n+\n{}\n", "I".repeat(READ_LEN)).unwrap();
    }
    out
}

fn gz(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

fn align(work: &Path, args: &[&str]) -> Output {
    let (_, genome_dir) = index();
    cargo_bin_cmd!("rustar-aligner")
        .current_dir(work)
        .args(["--runMode", "alignReads", "--runThreadN", "1"])
        .args(["--genomeDir", genome_dir.to_str().unwrap()])
        .args(args)
        .timeout(TIMEOUT)
        .output()
        .unwrap()
}

fn records(work: &Path) -> Vec<String> {
    fs::read_to_string(work.join("Aligned.out.sam"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_owned)
        .collect()
}

/// The SAM records of a run on `name`, a file holding `data`.
fn run(args: &[&str], name: &str, data: &[u8]) -> Vec<String> {
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join(name), data).unwrap();
    let mut all = args.to_vec();
    all.extend(["--readFilesIn", name]);
    let out = align(dir.path(), &all);
    assert!(
        out.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    records(dir.path())
}

fn plain() -> Vec<String> {
    let records = run(&[], "r.fq", &fastq());
    assert_eq!(records.len(), 40);
    records
}

#[test]
fn multi_word_command_in_one_value() {
    let compressed = gz(&fastq());
    let args = ["--readFilesCommand", "gzip -dc"];
    assert_eq!(run(&args, "r.fq.gz", &compressed), plain());
}

#[test]
fn multi_word_command_in_several_values() {
    let compressed = gz(&fastq());
    let args = ["--readFilesCommand", "gzip", "-d", "-c"];
    assert_eq!(run(&args, "r.fq.gz", &compressed), plain());
    let args = ["--readFilesCommand=gzip", "-dc"];
    assert_eq!(run(&args, "r.fq.gz", &compressed), plain());
}

/// A long option is a new parameter unless it is inside a quoted value.
#[test]
fn long_options_of_the_command_need_one_value() {
    let compressed = gz(&fastq());
    let args = ["--readFilesCommand", "gzip --decompress --stdout"];
    assert_eq!(run(&args, "r.fq.gz", &compressed), plain());

    let dir = TempDir::new().unwrap();
    let args = ["--readFilesCommand", "gzip", "--decompress"];
    let out = align(dir.path(), &args);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--decompress"));
}

#[test]
fn command_line_is_run_by_the_shell() {
    let compressed = gz(&fastq());
    let args = ["--readFilesCommand", "unset CDPATH; gzip -dc"];
    assert_eq!(run(&args, "r.fq.gz", &compressed), plain());
}

#[test]
fn dash_is_no_command() {
    let data = fastq();
    assert_eq!(run(&["--readFilesCommand", "-"], "r.fq", &data), plain());
    // What follows a first `-` is ignored.
    let args = ["--readFilesCommand", "-", "false"];
    assert_eq!(run(&args, "r.fq", &data), plain());
}

#[test]
fn dash_reads_several_files() {
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("a.fq"), fastq()).unwrap();
    fs::write(dir.path().join("b.fq"), fastq()).unwrap();
    let args = ["--readFilesCommand", "-", "--readFilesIn", "a.fq,b.fq"];
    let out = align(dir.path(), &args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(records(dir.path()).len(), 80);
}

/// Through the shell, a file name is still that file, whatever it holds.
#[test]
fn file_names_reach_the_command_unchanged() {
    let compressed = gz(&fastq());
    for name in ["with space.gz", "cost$HOME.gz", "a\"b'c.gz", "*.gz"] {
        let args = ["--readFilesCommand", "gzip -dc"];
        assert_eq!(run(&args, name, &compressed), plain(), "{name}");
    }
}

/// With a shell in between, a command that does not exist is the shell's exit
/// status 127, which must stay an error rather than an empty input.
#[test]
fn command_not_found_is_an_error() {
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r.fq"), fastq()).unwrap();
    let args = [
        "--readFilesCommand",
        "rustar-no-such-command -x",
        "--readFilesIn",
        "r.fq",
    ];
    let out = align(dir.path(), &args);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("exit status: 127"));
}

/// Neither can a failure further in pass for the end of the reads, and a run
/// that stops reading first does not make a failure out of the kill.
#[test]
fn failing_command_is_an_error_and_early_stop_is_not() {
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r.fq"), fastq()).unwrap();
    let args = ["--readFilesCommand", "gzip -dc", "--readFilesIn", "r.fq"];
    let out = align(dir.path(), &args);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("'gzip -dc' failed"));

    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r.fq"), fastq()).unwrap();
    // The file name is `$0` of the inner shell. `true;` keeps the outer shell
    // from replacing itself with it, so the kill misses the writer, which must
    // not hold up the run.
    let args = [
        "--readMapNumber",
        "5",
        "--readFilesCommand",
        "true; sh -c 'while cat \"$0\"; do :; done'",
        "--readFilesIn",
        "r.fq",
    ];
    let out = align(dir.path(), &args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(records(dir.path()).len(), 5);
}
