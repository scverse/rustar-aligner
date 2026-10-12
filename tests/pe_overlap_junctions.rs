//! A pair whose first mate is unspliced across a short intron that the second
//! mate splices over is a valid STAR alignment (STAR only compares junctions
//! that both overlapping mates have). Expected alignments are STAR 2.7.11b's on
//! the same genome and reads (`--genomeSAindexNbases 6`, default parameters).

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let bases: [u8; 4] = *b"ACGT";
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            bases[((state >> 16) & 3) as usize]
        })
        .collect()
}

fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b,
        })
        .collect()
}

fn write_fastq(path: &std::path::Path, reads: &[(String, Vec<u8>)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in reads {
        writeln!(f, "@{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f, "\n+\n{}", "I".repeat(seq.len())).unwrap();
    }
}

#[test]
fn unspliced_mate_may_overlap_the_junction_of_the_other_mate() {
    let tmp = TempDir::new().unwrap();
    // 8 kb background with two planted 21 bp GT..AG introns.
    let mut genome = lcg_seq(424242, 8000);
    let introns = [3000usize, 5000];
    for &p in &introns {
        genome[p] = b'G';
        genome[p + 1] = b'T';
        genome[p + 19] = b'A';
        genome[p + 20] = b'G';
    }
    let fasta = tmp.path().join("genome.fa");
    let mut f = fs::File::create(&fasta).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(&genome).unwrap();
    writeln!(f).unwrap();

    let genome_dir = tmp.path().join("genome");
    fs::create_dir_all(&genome_dir).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeDir"])
        .arg(&genome_dir)
        .arg("--genomeFastaFiles")
        .arg(&fasta)
        .args(["--genomeSAindexNbases", "6", "--outFileNamePrefix"])
        .arg(genome_dir.join("run_"))
        .assert()
        .success();

    // Mate 1: 150 contiguous bases across the intron. Mate 2: 11 bases of the
    // upstream exon, then 139 bases past the intron.
    let mut mate1 = Vec::new();
    let mut mate2 = Vec::new();
    for (i, &p) in introns.iter().enumerate() {
        mate1.push((format!("p{i}/1"), genome[p - 100..p + 50].to_vec()));
        let mut right = genome[p - 11..p].to_vec();
        right.extend_from_slice(&genome[p + 21..p + 160]);
        mate2.push((format!("p{i}/2"), rc(&right)));
    }
    let (r1, r2) = (tmp.path().join("r1.fq"), tmp.path().join("r2.fq"));
    write_fastq(&r1, &mate1);
    write_fastq(&r2, &mate2);

    let out = tmp.path().join("out");
    fs::create_dir_all(&out).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "alignReads", "--genomeDir"])
        .arg(&genome_dir)
        .arg("--readFilesIn")
        .arg(&r1)
        .arg(&r2)
        .arg("--outFileNamePrefix")
        .arg(format!("{}/", out.display()))
        .assert()
        .success();

    let sam = fs::read_to_string(out.join("Aligned.out.sam")).unwrap();
    let got: Vec<String> = sam
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let score = c.iter().find(|t| t.starts_with("AS:i:")).unwrap();
            format!("{} {} {} {} {}", c[0], c[1], c[3], c[5], score)
        })
        .collect();
    assert_eq!(
        got,
        vec![
            "p0 99 2901 150M AS:i:298",
            "p0 147 2990 11M21N139M AS:i:298",
            "p1 99 4901 150M AS:i:298",
            "p1 147 4990 11M21N139M AS:i:298",
        ]
    );
}
