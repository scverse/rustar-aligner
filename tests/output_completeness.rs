//! A run that exits 0 must have written everything it reports: `SJ.out.tab` even
//! with no junction (STAR writes it empty), and SAM output's last buffer.

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

fn lcg_seq(seed: u32, length: usize) -> String {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((state >> 16) & 3) as usize] as char
        })
        .collect()
}

/// One 5 kb contig indexed in `dir/idx`, and unspliced reads from it in `dir/r.fq`.
fn fixture(dir: &Path) -> (String, String) {
    let genome = lcg_seq(11, 5000);
    fs::write(dir.join("g.fa"), format!(">chr1\n{genome}\n")).unwrap();
    let reads: String = (0..20)
        .map(|i| {
            let s = &genome[100 + i * 50..200 + i * 50];
            format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(100))
        })
        .collect();
    fs::write(dir.join("r.fq"), reads).unwrap();
    let idx = dir.join("idx");
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeSAindexNbases", "5"])
        .arg("--genomeDir")
        .arg(&idx)
        .arg("--genomeFastaFiles")
        .arg(dir.join("g.fa"))
        .arg("--outFileNamePrefix")
        .arg(dir.join("gg_"))
        .assert()
        .success();
    (
        idx.to_str().unwrap().to_string(),
        dir.join("r.fq").to_str().unwrap().to_string(),
    )
}

#[test]
fn sj_out_tab_is_written_empty_without_junctions() {
    let dir = TempDir::new().unwrap();
    let (idx, reads) = fixture(dir.path());
    let prefix = dir.path().join("a_");
    cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", &idx, "--readFilesIn", &reads])
        .args(["--outSAMtype", "None"])
        .arg("--outFileNamePrefix")
        .arg(&prefix)
        .assert()
        .success();
    assert_eq!(fs::read(dir.path().join("a_SJ.out.tab")).unwrap(), b"");
}

#[test]
fn out_sj_type_none_still_writes_no_file() {
    let dir = TempDir::new().unwrap();
    let (idx, reads) = fixture(dir.path());
    cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", &idx, "--readFilesIn", &reads])
        .args(["--outSAMtype", "None", "--outSJtype", "None"])
        .arg("--outFileNamePrefix")
        .arg(dir.path().join("a_"))
        .assert()
        .success();
    assert!(!dir.path().join("a_SJ.out.tab").exists());
}

/// `/dev/full` fails every write, so the only write this small run makes, the
/// flush of the SAM buffer, has to report it.
#[cfg(target_os = "linux")]
#[test]
fn failed_sam_flush_fails_the_run() {
    let dir = TempDir::new().unwrap();
    let (idx, reads) = fixture(dir.path());
    std::os::unix::fs::symlink("/dev/full", dir.path().join("a_Aligned.out.sam")).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", &idx, "--readFilesIn", &reads])
        .arg("--outFileNamePrefix")
        .arg(dir.path().join("a_"))
        .assert()
        .failure();
}
