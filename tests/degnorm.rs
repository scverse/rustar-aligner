//! End-to-end DegNorm: align two synthetic samples with
//! `--quantMode GeneCoverage`, then merge them with `--runMode degNorm`.
//!
//! Sample `A` covers both genes evenly. Sample `B` is identical except that
//! gene `G1` is only sequenced over its 5' 40%, which is what 3' transcript
//! degradation looks like in coverage space. The DI score for `G1` in `B` must
//! therefore be clearly higher than in `A`, while `G2` stays undegraded in both.

use assert_cmd::cargo::cargo_bin_cmd;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// LCG pseudo-random sequence generator (same LCG as `tests/alignment_features.rs`).
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

/// Single-exon genes, 0-based half-open: G1 = [2000, 4000), G2 = [6000, 8000).
const G1: (usize, usize) = (2000, 4000);
const G2: (usize, usize) = (6000, 8000);

fn write_fasta(dir: &Path, genome: &[u8]) -> PathBuf {
    let path = dir.join("genome.fa");
    let mut f = fs::File::create(&path).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(genome).unwrap();
    writeln!(f).unwrap();
    path
}

fn write_gtf(dir: &Path) -> PathBuf {
    let path = dir.join("genes.gtf");
    let mut f = fs::File::create(&path).unwrap();
    for (gene, (start, end)) in [("G1", G1), ("G2", G2)] {
        writeln!(
            f,
            "chr1\ttest\texon\t{}\t{}\t.\t+\t.\tgene_id \"{gene}\"; transcript_id \"{gene}_T1\";",
            start + 1,
            end
        )
        .unwrap();
    }
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

/// Tile `[from, to)` with 100 bp reads every `stride` bases.
fn tile_reads(f: &mut fs::File, genome: &[u8], label: &str, from: usize, to: usize, stride: usize) {
    let read_len = 100;
    let mut pos = from;
    let mut i = 0;
    while pos + read_len <= to {
        writeln!(f, "@{label}_{i}").unwrap();
        f.write_all(&genome[pos..pos + read_len]).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(read_len)).unwrap();
        pos += stride;
        i += 1;
    }
}

/// FASTQ for one sample. When `degrade_g1` is set, G1's reads only come from
/// its 5' 40%.
fn write_fastq(dir: &Path, genome: &[u8], name: &str, degrade_g1: bool) -> PathBuf {
    let path = dir.join(format!("{name}.fq"));
    let mut f = fs::File::create(&path).unwrap();
    let g1_end = if degrade_g1 {
        G1.0 + (G1.1 - G1.0) * 2 / 5
    } else {
        G1.1
    };
    tile_reads(&mut f, genome, &format!("{name}g1"), G1.0, g1_end, 25);
    tile_reads(&mut f, genome, &format!("{name}g2"), G2.0, G2.1, 25);
    path
}

fn align(genome_dir: &Path, fastq: &Path, prefix: &str, quant_mode: &[&str]) {
    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.args([
        "--runMode",
        "alignReads",
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--readFilesIn",
        fastq.to_str().unwrap(),
        "--outFileNamePrefix",
        prefix,
        "--quantMode",
    ]);
    cmd.args(quant_mode);
    cmd.arg("--sjdbGTFfile");
    cmd.arg(genome_dir.parent().unwrap().join("genes.gtf"));
    cmd.assert().success();
}

/// Parse a `gene<TAB>sample...` matrix into `gene -> (sample -> value)`.
fn read_matrix(path: &Path) -> (Vec<String>, HashMap<String, Vec<f64>>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .skip(1)
        .map(String::from)
        .collect();
    let mut rows = HashMap::new();
    for line in lines {
        let mut it = line.split('\t');
        let gene = it.next().unwrap().to_string();
        rows.insert(gene, it.map(|v| v.parse::<f64>().unwrap()).collect());
    }
    (header, rows)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn degnorm_end_to_end_flags_the_degraded_sample() {
    let tmpdir = TempDir::new().unwrap();
    let root = tmpdir.path();
    let genome = lcg_seq(88888, 20000);
    let fasta = write_fasta(root, &genome);
    write_gtf(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &genome_dir);

    let fq_a = write_fastq(root, &genome, "A", false);
    let fq_b = write_fastq(root, &genome, "B", true);

    let prefix_a = format!("{}/A_", root.display());
    let prefix_b = format!("{}/B_", root.display());
    align(
        &genome_dir,
        &fq_a,
        &prefix_a,
        &["GeneCounts", "GeneCoverage"],
    );
    align(
        &genome_dir,
        &fq_b,
        &prefix_b,
        &["GeneCounts", "GeneCoverage"],
    );

    let cov_a = root.join("A_GeneCoverage.out.bin");
    let cov_b = root.join("B_GeneCoverage.out.bin");
    assert!(cov_a.exists(), "A_GeneCoverage.out.bin not written");
    assert!(cov_b.exists(), "B_GeneCoverage.out.bin not written");

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "degNorm",
            "--degNormCoverageFiles",
            cov_a.to_str().unwrap(),
            cov_b.to_str().unwrap(),
            "--degNormIter",
            "2",
            "--degNormNmfIter",
            "50",
            "--degNormMinHighCoverage",
            "20",
            "--outFileNamePrefix",
            &format!("{}/", root.display()),
        ])
        .assert()
        .success();

    let out_dir = root.join("DegNorm.out");
    for name in [
        "DegradationIndex.tab",
        "AdjustedCounts.tab",
        "RawCounts.tab",
        "ScaleFactors.tab",
        "Summary.txt",
    ] {
        assert!(out_dir.join(name).exists(), "{name} missing");
    }

    let (samples, di) = read_matrix(&out_dir.join("DegradationIndex.tab"));
    assert_eq!(samples, vec!["A".to_string(), "B".to_string()]);

    let g1 = di.get("G1").expect("G1 missing from DI table");
    let g2 = di.get("G2").expect("G2 missing from DI table");
    assert!(
        g1[1] > g1[0] + 0.2,
        "G1 DI should be much higher in the degraded sample: A={}, B={}",
        g1[0],
        g1[1]
    );
    assert!(
        g2[0] < 0.1 && g2[1] < 0.1,
        "G2 is undegraded in both samples but got DI A={}, B={}",
        g2[0],
        g2[1]
    );

    // Degradation-adjusted counts scale the degraded sample's G1 count up.
    let (_, adjusted) = read_matrix(&out_dir.join("AdjustedCounts.tab"));
    let (_, raw) = read_matrix(&out_dir.join("RawCounts.tab"));
    assert!(raw["G1"][1] > 0.0, "G1 has no reads in sample B");
    assert!(
        adjusted["G1"][1] > raw["G1"][1],
        "adjusted count {} should exceed raw count {}",
        adjusted["G1"][1],
        raw["G1"][1]
    );
}

#[test]
fn gene_coverage_does_not_change_alignment_output() {
    let tmpdir = TempDir::new().unwrap();
    let root = tmpdir.path();
    let genome = lcg_seq(88888, 20000);
    let fasta = write_fasta(root, &genome);
    write_gtf(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &genome_dir);
    let fq = write_fastq(root, &genome, "A", false);

    let with_prefix = format!("{}/with_", root.display());
    let without_prefix = format!("{}/without_", root.display());
    align(&genome_dir, &fq, &without_prefix, &["GeneCounts"]);
    align(
        &genome_dir,
        &fq,
        &with_prefix,
        &["GeneCounts", "GeneCoverage"],
    );

    // Compare alignment records only: the @PG header records the command line,
    // which necessarily differs by the extra --quantMode word.
    let records = |p: PathBuf| -> Vec<String> {
        fs::read_to_string(p)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(String::from)
            .collect()
    };
    let a = records(root.join("without_Aligned.out.sam"));
    let b = records(root.join("with_Aligned.out.sam"));
    assert!(!a.is_empty(), "no alignments produced");
    assert_eq!(a, b, "GeneCoverage changed the SAM output");

    // ReadsPerGene.out.tab must be identical too.
    let counts_a = fs::read(root.join("without_ReadsPerGene.out.tab")).unwrap();
    let counts_b = fs::read(root.join("with_ReadsPerGene.out.tab")).unwrap();
    assert_eq!(counts_a, counts_b, "GeneCoverage changed ReadsPerGene");
}
