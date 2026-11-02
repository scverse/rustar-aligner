//! FASTQ input read from streams (FIFOs, pipes, process substitution): the
//! output must be what the same reads give from a file, a mode that would read
//! the input twice must refuse a stream, and a cut stream must not end the run
//! with exit code 0.
#![cfg(unix)]

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::OnceLock;
use std::time::Duration;
use tempfile::TempDir;

/// A command that hangs fails here instead of blocking the test suite.
const TIMEOUT: Duration = Duration::from_secs(60);

const READ_LEN: usize = 50;
const INTRON: (usize, usize) = (2_300, 2_800);
const CB_WHITELIST: [&str; 2] = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT"];

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((state >> 16) & 3) as usize]
        })
        .collect()
}

fn genome() -> &'static [u8] {
    static GENOME: OnceLock<Vec<u8>> = OnceLock::new();
    GENOME.get_or_init(|| {
        let mut g = lcg_seq(88888, 20_000);
        g[INTRON.0..INTRON.0 + 2].copy_from_slice(b"GT");
        g[INTRON.1 - 2..INTRON.1].copy_from_slice(b"AG");
        g
    })
}

/// Index shared by every test: genome, one spliced gene.
fn index() -> &'static (TempDir, PathBuf) {
    static INDEX: OnceLock<(TempDir, PathBuf)> = OnceLock::new();
    INDEX.get_or_init(|| {
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("genome.fa");
        let mut f = fs::File::create(&fasta).unwrap();
        f.write_all(b">chr1\n").unwrap();
        f.write_all(genome()).unwrap();
        f.write_all(b"\n").unwrap();
        let gtf = dir.path().join("genes.gtf");
        fs::write(
            &gtf,
            "chr1\tt\texon\t2001\t2300\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\"; gene_name \"G\";\n\
             chr1\tt\texon\t2801\t3100\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\"; gene_name \"G\";\n",
        )
        .unwrap();
        let genome_dir = dir.path().join("idx");
        fs::create_dir(&genome_dir).unwrap();
        cargo_bin_cmd!("rustar-aligner")
            .current_dir(dir.path())
            .args(["--runMode", "genomeGenerate", "--genomeDir", "idx"])
            .args(["--genomeFastaFiles", "genome.fa", "--sjdbGTFfile", "genes.gtf"])
            .args(["--sjdbOverhang", "49", "--genomeSAindexNbases", "7"])
            .assert()
            .success();
        (dir, genome_dir)
    })
}

/// `n` reads per mate: every other mate 1 crosses the junction.
fn fastq(n: usize) -> (Vec<u8>, Vec<u8>) {
    let g = genome();
    let (mut r1, mut r2) = (Vec::new(), Vec::new());
    for i in 0..n {
        let k = 10 + i % 30;
        let seq1: Vec<u8> = if i % 2 == 0 {
            [
                &g[INTRON.0 - k..INTRON.0],
                &g[INTRON.1..INTRON.1 + READ_LEN - k],
            ]
            .concat()
        } else {
            let p = 4_000 + i;
            g[p..p + READ_LEN].to_vec()
        };
        let p = if i % 2 == 0 {
            INTRON.1 + 100
        } else {
            4_200 + i
        };
        let seq2: Vec<u8> = g[p..p + READ_LEN]
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        for (out, seq) in [(&mut r1, seq1), (&mut r2, seq2)] {
            writeln!(out, "@r{i}").unwrap();
            out.extend_from_slice(&seq);
            write!(out, "\n+\n{}\n", "I".repeat(READ_LEN)).unwrap();
        }
    }
    (r1, r2)
}

fn gz(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

/// Create a FIFO and feed it once from a background thread, as a download would.
/// The thread stays blocked if nothing ever opens the FIFO, and dies with the
/// test process.
fn fifo(path: &Path, data: Vec<u8>) {
    assert!(
        Command::new("mkfifo").arg(path).status().unwrap().success(),
        "mkfifo"
    );
    let path = path.to_path_buf();
    std::thread::spawn(move || {
        if let Ok(mut f) = fs::OpenOptions::new().write(true).open(path) {
            // A reader that stops early (--readMapNumber) closes the pipe.
            let _ = f.write_all(&data);
        }
    });
}

/// Run an alignment in a directory of its own, with `args` as given.
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

/// Same run from a shell, for what only a shell provides: `<(...)`.
fn align_in_bash(work: &Path, args: &str) -> Output {
    let (_, genome_dir) = index();
    let bin = assert_cmd::cargo::cargo_bin!("rustar-aligner");
    // `exec`: on a timeout, the aligner itself is what gets killed.
    let script = format!(
        "exec {} --runMode alignReads --runThreadN 1 --genomeDir {} {args}",
        bin.display(),
        genome_dir.display(),
    );
    assert_cmd::Command::new("bash")
        .current_dir(work)
        .args(["-c", &script])
        .timeout(TIMEOUT)
        .output()
        .unwrap()
}

fn assert_ok(out: &Output) {
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// The run failed by itself, with `needle` in its message, and did not present
/// a result as complete.
fn assert_failed(out: &Output, work: &Path, needle: &str) {
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.code().is_some_and(|c| c != 0),
        "expected an error exit: {:?}\nstderr: {stderr}",
        out.status
    );
    assert!(stderr.contains(needle), "stderr: {stderr}");
    assert!(!work.join("Log.final.out").exists());
}

/// Output files compared between runs. The header's command line differs by
/// construction, and the log's timings by nature.
fn results(work: &Path) -> Vec<String> {
    let sam = fs::read_to_string(work.join("Aligned.out.sam")).unwrap();
    let sam: Vec<&str> = sam
        .lines()
        .filter(|l| !l.starts_with("@PG") && !l.starts_with("@CO"))
        .collect();
    let log = fs::read_to_string(work.join("Log.final.out")).unwrap();
    let log: Vec<&str> = log
        .lines()
        .filter(|l| !l.contains(" on ") && !l.contains("time") && !l.contains("speed"))
        .collect();
    vec![
        sam.join("\n"),
        fs::read_to_string(work.join("SJ.out.tab")).unwrap_or_default(),
        log.join("\n"),
    ]
}

fn input_reads(work: &Path) -> u64 {
    fs::read_to_string(work.join("Log.final.out"))
        .unwrap()
        .lines()
        .find(|l| l.contains("Number of input reads"))
        .and_then(|l| l.rsplit('\t').next()?.trim().parse().ok())
        .unwrap()
}

/// Reference results for `extra` options, from regular files.
fn from_files(extra: &[&str], r1: &[u8], r2: Option<&[u8]>) -> Vec<String> {
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("in1"), r1).unwrap();
    let mut args = extra.to_vec();
    args.extend(["--readFilesIn", "in1"]);
    if let Some(r2) = r2 {
        fs::write(dir.path().join("in2"), r2).unwrap();
        args.push("in2");
    }
    assert_ok(&align(dir.path(), &args));
    results(dir.path())
}

#[test]
fn fifo_gives_the_output_of_the_file() {
    let (r1, r2) = fastq(500);
    let (gz1, gz2) = (gz(&r1), gz(&r2));
    for (label, a, b) in [("plain", &r1, &r2), ("gzip", &gz1, &gz2)] {
        for extra in [
            &[][..],
            &["--outFilterType", "BySJout"][..],
            &["--readFilesCommand", "cat"][..],
        ] {
            // `cat` hands gzip to the parser, so that command is plain only.
            if label == "gzip" && extra.contains(&"cat") {
                continue;
            }
            let dir = TempDir::new().unwrap();
            fifo(&dir.path().join("in1"), a.clone());
            let mut args = extra.to_vec();
            args.extend(["--readFilesIn", "in1"]);
            assert_ok(&align(dir.path(), &args));
            assert_eq!(
                results(dir.path()),
                from_files(extra, a, None),
                "SE {label} {extra:?}"
            );

            let dir = TempDir::new().unwrap();
            fifo(&dir.path().join("in1"), a.clone());
            fifo(&dir.path().join("in2"), b.clone());
            let mut args = extra.to_vec();
            args.extend(["--readFilesIn", "in1", "in2"]);
            assert_ok(&align(dir.path(), &args));
            assert_eq!(
                results(dir.path()),
                from_files(extra, a, Some(b)),
                "PE {label} {extra:?}"
            );
        }
    }
}

#[test]
fn process_substitution_gives_the_output_of_the_file() {
    let (r1, r2) = fastq(500);
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r1.gz"), gz(&r1)).unwrap();
    fs::write(dir.path().join("r2.gz"), gz(&r2)).unwrap();
    assert_ok(&align_in_bash(
        dir.path(),
        "--readFilesIn <(cat r1.gz) <(cat r2.gz)",
    ));
    assert_eq!(results(dir.path()), from_files(&[], &r1, Some(&r2)));
}

#[test]
fn readfilesprefix_with_a_dev_fd_number_reads_the_stream() {
    let (r1, _) = fastq(500);
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r1"), &r1).unwrap();
    assert_ok(&align_in_bash(
        dir.path(),
        "--readFilesPrefix /dev/fd/ --readFilesIn 9 9< <(cat r1)",
    ));
    assert_eq!(results(dir.path()), from_files(&[], &r1, None));
}

#[test]
fn stdin_pipe_gives_the_output_of_the_file() {
    let (r1, _) = fastq(500);
    let dir = TempDir::new().unwrap();
    let (_, genome_dir) = index();
    let out = cargo_bin_cmd!("rustar-aligner")
        .current_dir(dir.path())
        .args(["--runMode", "alignReads", "--runThreadN", "1"])
        .args(["--genomeDir", genome_dir.to_str().unwrap()])
        .args(["--readFilesIn", "/dev/stdin"])
        .write_stdin(gz(&r1))
        .timeout(TIMEOUT)
        .output()
        .unwrap();
    assert_ok(&out);
    assert_eq!(results(dir.path()), from_files(&[], &r1, None));
}

/// Files of a list are opened in turn, so a later FIFO's writer only waits.
#[test]
fn comma_separated_fifos_give_the_output_of_the_files() {
    let (r1, r2) = fastq(500);
    // Every read of the fixture has the same number of lines, four.
    let split = |d: &[u8]| {
        let lines: Vec<&[u8]> = d.split_inclusive(|&b| b == b'\n').collect();
        let (a, b) = lines.split_at(lines.len() / 8 * 4);
        (a.concat(), b.concat())
    };
    let (a1, b1) = split(&r1);
    let (a2, b2) = split(&r2);
    let dir = TempDir::new().unwrap();
    for (name, data) in [("a1", a1), ("b1", b1), ("a2", a2), ("b2", b2)] {
        fifo(&dir.path().join(name), gz(&data));
    }
    assert_ok(&align(dir.path(), &["--readFilesIn", "a1,b1", "a2,b2"]));
    assert_eq!(results(dir.path()), from_files(&[], &r1, Some(&r2)));
}

/// Stopping early leaves the writer with a closed pipe, and never opens the
/// next FIFO of a list; neither stops the run.
#[test]
fn readmapnumber_stops_on_a_fifo() {
    let (r1, _) = fastq(500);
    let limit = ["--readMapNumber", "10"];
    let dir = TempDir::new().unwrap();
    fifo(&dir.path().join("a"), gz(&r1));
    fifo(&dir.path().join("never_opened"), gz(&r1));
    let mut args = limit.to_vec();
    args.extend(["--readFilesIn", "a,never_opened"]);
    assert_ok(&align(dir.path(), &args));
    assert_eq!(input_reads(dir.path()), 10);
    assert_eq!(results(dir.path()), from_files(&limit, &r1, None));
}

fn barcodes(n: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..n {
        let umi = String::from_utf8(lcg_seq(i as u32 + 1, 12)).unwrap();
        write!(
            out,
            "@r{i}\n{}{umi}\n+\n{}\n",
            CB_WHITELIST[i % 2],
            "I".repeat(28)
        )
        .unwrap();
    }
    out
}

fn solo_args(whitelist: &str) -> Vec<&str> {
    vec![
        "--soloType",
        "CB_UMI_Simple",
        "--soloCBwhitelist",
        whitelist,
        "--soloCBlen",
        "16",
        "--soloUMIstart",
        "17",
        "--soloUMIlen",
        "12",
        "--soloBarcodeReadLength",
        "0",
        "--readFilesIn",
        "cdna",
        "barcode",
    ]
}

fn write_whitelist(dir: &Path) {
    fs::write(dir.join("wl"), CB_WHITELIST.join("\n") + "\n").unwrap();
}

fn solo_matrix(dir: &Path) -> String {
    fs::read_to_string(dir.join("Solo.out/Gene/raw/matrix.mtx")).unwrap()
}

#[test]
fn solo_fifos_give_the_counts_of_the_files() {
    let (r1, _) = fastq(500);
    let bc = barcodes(500);

    let files = TempDir::new().unwrap();
    write_whitelist(files.path());
    fs::write(files.path().join("cdna"), &r1).unwrap();
    fs::write(files.path().join("barcode"), &bc).unwrap();
    assert_ok(&align(files.path(), &solo_args("wl")));

    let streams = TempDir::new().unwrap();
    write_whitelist(streams.path());
    fifo(&streams.path().join("cdna"), gz(&r1));
    fifo(&streams.path().join("barcode"), gz(&bc));
    assert_ok(&align(streams.path(), &solo_args("wl")));

    assert_eq!(solo_matrix(streams.path()), solo_matrix(files.path()));
    assert_eq!(results(streams.path()), results(files.path()));
}

#[test]
fn smartseq_manifest_of_fifos_gives_the_counts_of_the_files() {
    let (r1, _) = fastq(200);
    let run = |make_input: &dyn Fn(&Path, Vec<u8>)| {
        let dir = TempDir::new().unwrap();
        make_input(&dir.path().join("cell1"), r1.clone());
        make_input(&dir.path().join("cell2"), r1.clone());
        let manifest = format!("{0}/cell1\t-\tc1\n{0}/cell2\t-\tc2\n", dir.path().display());
        fs::write(dir.path().join("manifest"), manifest).unwrap();
        assert_ok(&align(
            dir.path(),
            &[
                "--soloType",
                "SmartSeq",
                "--soloUMIdedup",
                "Exact",
                "--outSAMtype",
                "None",
                "--readFilesManifest",
                "manifest",
            ],
        ));
        solo_matrix(dir.path())
    };
    let from_fifos = run(&|path, data| fifo(path, data));
    let from_files = run(&|path, data| fs::write(path, data).unwrap());
    assert_eq!(from_fifos, from_files);
}

const TWOPASS_REFUSED: &str = "reads the input twice";

/// Without the check this blocks forever on the second open of the FIFO.
#[test]
fn twopass_refuses_a_fifo() {
    let (r1, r2) = fastq(50);
    let dir = TempDir::new().unwrap();
    fifo(&dir.path().join("in1"), r1.clone());
    let out = align(
        dir.path(),
        &["--twopassMode", "Basic", "--readFilesIn", "in1"],
    );
    assert_failed(&out, dir.path(), TWOPASS_REFUSED);

    // Mate 2 alone, and a stream that is not the first of its list.
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("in1"), &r1).unwrap();
    fifo(&dir.path().join("in2"), r2);
    let out = align(
        dir.path(),
        &["--twopassMode", "Basic", "--readFilesIn", "in1", "in2"],
    );
    assert_failed(&out, dir.path(), TWOPASS_REFUSED);

    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("a"), &r1).unwrap();
    fifo(&dir.path().join("b"), r1);
    let out = align(
        dir.path(),
        &["--twopassMode", "Basic", "--readFilesIn", "a,b"],
    );
    assert_failed(&out, dir.path(), TWOPASS_REFUSED);
}

/// Without the check this exits 0 with a second pass of 0 reads.
#[test]
fn twopass_refuses_process_substitution() {
    let (r1, _) = fastq(50);
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r1"), &r1).unwrap();
    let out = align_in_bash(dir.path(), "--twopassMode Basic --readFilesIn <(cat r1)");
    assert_failed(&out, dir.path(), TWOPASS_REFUSED);
    let out = align_in_bash(
        dir.path(),
        "--twopassMode Basic --readFilesCommand cat --readFilesIn <(cat r1)",
    );
    assert_failed(&out, dir.path(), TWOPASS_REFUSED);
}

/// `/dev/stdin` redirected from a file stats as a regular file, yet on macOS
/// opening it again does not rewind it.
#[test]
fn twopass_refuses_dev_stdin_redirected_from_a_file() {
    let (r1, _) = fastq(50);
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r1"), &r1).unwrap();
    let out = align_in_bash(
        dir.path(),
        "--twopassMode Basic --readFilesIn /dev/stdin < r1",
    );
    assert_failed(&out, dir.path(), TWOPASS_REFUSED);
}

#[test]
fn twopass_still_reads_regular_files() {
    let (r1, _) = fastq(50);
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("r1"), &r1).unwrap();
    assert_ok(&align(
        dir.path(),
        &["--twopassMode", "Basic", "--readFilesIn", "r1"],
    ));
    assert_eq!(input_reads(dir.path()), 50);
}

/// Solo aligns only the cDNA in pass 1, so a barcode stream is read once.
#[test]
fn twopass_solo_accepts_a_barcode_stream_but_not_a_cdna_stream() {
    let (r1, _) = fastq(200);
    let bc = barcodes(200);
    let mut args = vec!["--twopassMode", "Basic"];
    args.extend(solo_args("wl"));

    let dir = TempDir::new().unwrap();
    write_whitelist(dir.path());
    fs::write(dir.path().join("cdna"), &r1).unwrap();
    fifo(&dir.path().join("barcode"), bc.clone());
    assert_ok(&align(dir.path(), &args));
    assert_eq!(input_reads(dir.path()), 200);

    let dir = TempDir::new().unwrap();
    write_whitelist(dir.path());
    fifo(&dir.path().join("cdna"), r1);
    fs::write(dir.path().join("barcode"), &bc).unwrap();
    assert_failed(&align(dir.path(), &args), dir.path(), TWOPASS_REFUSED);
}

fn run_on_fifo(data: Vec<u8>) -> (TempDir, Output) {
    let dir = TempDir::new().unwrap();
    fifo(&dir.path().join("in"), data);
    let out = align(dir.path(), &["--readFilesIn", "in"]);
    (dir, out)
}

#[test]
fn gzip_stream_cut_is_an_error() {
    let (r1, _) = fastq(200);
    let whole = gz(&r1);
    for (label, cut) in [
        ("in the header", 5),
        ("in the deflate data", whole.len() / 2),
        ("in the trailer", whole.len() - 4),
        ("before the trailer", whole.len() - 8),
    ] {
        let (dir, out) = run_on_fifo(whole[..cut].to_vec());
        assert_failed(&out, dir.path(), "truncated gzip input");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains("(in)"), "{label}: names the file: {stderr}");
    }
}

/// Paired-end: the mate with the cut stream is what the error names.
#[test]
fn gzip_cut_on_one_mate_is_an_error() {
    let (r1, r2) = fastq(200);
    let dir = TempDir::new().unwrap();
    let cut = gz(&r2);
    fifo(&dir.path().join("in1"), gz(&r1));
    fifo(&dir.path().join("in2"), cut[..cut.len() / 2].to_vec());
    let out = align(dir.path(), &["--readFilesIn", "in1", "in2"]);
    assert_failed(&out, dir.path(), "truncated gzip input");
    assert!(String::from_utf8_lossy(&out.stderr).contains("(in2)"));
}

/// Every cut inside the last read fails. Only the cuts that leave complete
/// reads pass: at the boundary and just before the final newline.
#[test]
fn plain_stream_cut_mid_read_is_an_error_and_on_a_boundary_is_not() {
    let (r1, _) = fastq(3);
    let last_start = r1.len() - r1.rsplit(|&b| b == b'@').next().unwrap().len() - 1;
    for cut in (last_start + 1)..(r1.len() - 1) {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("in"), &r1[..cut]).unwrap();
        let out = align(dir.path(), &["--readFilesIn", "in"]);
        assert!(
            out.status.code().is_some_and(|c| c != 0),
            "cut {} bytes into the read exited {:?}",
            cut - last_start,
            out.status
        );
        assert!(!dir.path().join("Log.final.out").exists());
    }

    // The limit of what the aligner can see, recorded in the docs: a cut
    // between reads looks like a shorter input and only the read count shows it.
    let (dir, out) = run_on_fifo(r1[..last_start].to_vec());
    assert_ok(&out);
    assert_eq!(input_reads(dir.path()), 2);
    let (dir, out) = run_on_fifo(r1.clone());
    assert_ok(&out);
    assert_eq!(input_reads(dir.path()), 3);
}

/// An executable script standing for a decompressor.
fn command_script(dir: &Path, body: &str) -> String {
    use std::os::unix::fs::PermissionsExt;
    let path = dir.join("decompress.sh");
    fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
    fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
    path.to_str().unwrap().to_owned()
}

/// `zcat` on a cut file writes what it could, then exits non-zero. Complete
/// reads followed by a failure status must not pass for the whole input.
#[test]
fn readfilescommand_failure_is_an_error() {
    let (r1, _) = fastq(50);
    let dir = TempDir::new().unwrap();
    let cmd = command_script(dir.path(), "cat \"$1\"\nexit 3");
    fs::write(dir.path().join("in"), &r1).unwrap();
    let out = align(
        dir.path(),
        &["--readFilesCommand", &cmd, "--readFilesIn", "in"],
    );
    assert_failed(&out, dir.path(), "failed");
    assert!(String::from_utf8_lossy(&out.stderr).contains("exit status: 3"));

    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("in"), &r1).unwrap();
    let out = align(
        dir.path(),
        &["--readFilesCommand", "false", "--readFilesIn", "in"],
    );
    assert_failed(&out, dir.path(), "'false' failed");
}

#[test]
fn readfilescommand_success_and_early_stop_are_not_errors() {
    let (r1, _) = fastq(50);
    let dir = TempDir::new().unwrap();
    let cmd = command_script(dir.path(), "cat \"$1\"");
    fs::write(dir.path().join("in"), &r1).unwrap();
    assert_ok(&align(
        dir.path(),
        &["--readFilesCommand", &cmd, "--readFilesIn", "in"],
    ));
    assert_eq!(input_reads(dir.path()), 50);

    // Stopped by --readMapNumber: the command, still writing and bound to fail
    // in the end, is killed, and that is not a failure.
    let dir = TempDir::new().unwrap();
    let cmd = command_script(
        dir.path(),
        "i=0; while [ $i -lt 20000 ]; do cat \"$1\"; i=$((i+1)); done; exit 3",
    );
    fs::write(dir.path().join("in"), &r1).unwrap();
    assert_ok(&align(
        dir.path(),
        &[
            "--readMapNumber",
            "10",
            "--readFilesCommand",
            &cmd,
            "--readFilesIn",
            "in",
        ],
    ));
    assert_eq!(input_reads(dir.path()), 10);
}
