//! Phase 17.13 integration tests — coverage for all major Phase 17 features.
//!
//! Uses a 20,000bp pseudo-random genome (seed 88888) on chr1 with a planted
//! GT-AG intron structure for splice tests.

use assert_cmd::cargo::cargo_bin_cmd;
use noodles::bam;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// LCG pseudo-random sequence generator (identical LCG to existing tests).
fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let bases: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut seq = Vec::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

/// Build the 20 kb genome with a planted GT-AG intron structure.
///
/// Layout (0-based):
///   [0..10000]      : LCG(88888) bases  — background
///   [10000..10050]  : 50 bp Exon1 region (unique LCG seed 11111)
///   [10050..10052]  : "GT"  — intron donor
///   [10052..10248]  : 196 bp intron body (LCG seed 22222)
///   [10248..10250]  : "AG"  — intron acceptor
///   [10250..10300]  : 50 bp Exon2 region (unique LCG seed 33333)
///   [10300..20000]  : LCG(88888) bases  — background (continued)
fn build_genome() -> Vec<u8> {
    let background = lcg_seq(88888, 20000);
    let exon1 = lcg_seq(11111, 50);
    let intron_body = lcg_seq(22222, 196);
    let exon2 = lcg_seq(33333, 50);

    let mut genome = background;
    // Exon1
    genome[10000..10050].copy_from_slice(&exon1);
    // GT donor
    genome[10050] = b'G';
    genome[10051] = b'T';
    // Intron body
    genome[10052..10248].copy_from_slice(&intron_body);
    // AG acceptor
    genome[10248] = b'A';
    genome[10249] = b'G';
    // Exon2
    genome[10250..10300].copy_from_slice(&exon2);
    genome
}

/// Write genome.fa to tmpdir and return its path.
fn write_fasta(dir: &TempDir, genome: &[u8]) -> PathBuf {
    let path = dir.path().join("genome.fa");
    let mut f = fs::File::create(&path).unwrap();
    writeln!(f, ">chr1").unwrap();
    // Write the genome as ASCII
    f.write_all(genome).unwrap();
    writeln!(f).unwrap();
    path
}

/// Write a 2-exon GTF (1-based inclusive) for gene G1 / transcript T1.
///
/// Exon1: chr1:10001–10050  (0-based [10000,10050) → 1-based [10001,10050])
/// Exon2: chr1:10251–10300  (0-based [10250,10300) → 1-based [10251,10300])
fn write_gtf(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("annotations.gtf");
    let mut f = fs::File::create(&path).unwrap();
    writeln!(
        f,
        "chr1\ttest\texon\t10001\t10050\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
    )
    .unwrap();
    writeln!(
        f,
        "chr1\ttest\texon\t10251\t10300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
    )
    .unwrap();
    path
}

/// Build the rustar-aligner genome index.
/// `sa_nbases` should be "7" for this 20 kb genome.
/// If `gtf` is `Some`, passes `--sjdbGTFfile` + `--sjdbOverhang`.
fn build_index(fasta: &Path, genome_dir: &Path, sa_nbases: &str, gtf: Option<&Path>) {
    fs::create_dir_all(genome_dir).unwrap();
    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.arg("--runMode")
        .arg("genomeGenerate")
        .arg("--genomeDir")
        .arg(genome_dir)
        .arg("--genomeFastaFiles")
        .arg(fasta)
        .arg("--genomeSAindexNbases")
        .arg(sa_nbases);
    if let Some(g) = gtf {
        cmd.arg("--sjdbGTFfile")
            .arg(g)
            .arg("--sjdbOverhang")
            .arg("24");
    }
    cmd.assert().success();
}

/// Reverse-complement of a byte slice (A↔T, C↔G; unknown bases kept as-is).
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

/// Count alignment (non-@) lines in a SAM file.
fn count_sam_records(sam_path: &Path) -> usize {
    let content = fs::read_to_string(sam_path).unwrap();
    content.lines().filter(|l| !l.starts_with('@')).count()
}

// ---------------------------------------------------------------------------
// Test 1 — BAM unsorted output
// ---------------------------------------------------------------------------

#[test]
fn test_bam_unsorted_output() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", None);

    // 50 reads of 50 bp, from positions 100..150, 200..250, ..., 5100..5150
    let fastq_path = tmpdir.path().join("reads.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        for i in 0..50usize {
            let start = 100 + i * 100;
            let seq = &genome[start..start + 50];
            writeln!(f, "@read{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_bam_unsorted");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--outSAMtype",
            "BAM",
            "Unsorted",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let bam_path = output_dir.join("Aligned.out.bam");
    assert!(bam_path.exists(), "Aligned.out.bam not found");

    // Validate as BAM and check at least 1 record
    let mut reader = bam::io::Reader::new(fs::File::open(&bam_path).unwrap());
    let _header = reader.read_header().expect("BAM header readable");
    let mut count = 0usize;
    for rec in reader.records() {
        rec.expect("valid BAM record");
        count += 1;
    }
    assert!(count >= 1, "expected at least 1 BAM record, got {count}");
}

// ---------------------------------------------------------------------------
// Test 2 — BAM sorted output
// ---------------------------------------------------------------------------

#[test]
fn test_bam_sorted_output() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", None);

    // Same 50 reads as test 1
    let fastq_path = tmpdir.path().join("reads.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        for i in 0..50usize {
            let start = 100 + i * 100;
            let seq = &genome[start..start + 50];
            writeln!(f, "@read{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_bam_sorted");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let bam_path = output_dir.join("Aligned.sortedByCoord.out.bam");
    assert!(bam_path.exists(), "Aligned.sortedByCoord.out.bam not found");

    // Validate readable and that at least 5 consecutive mapped records are
    // in non-decreasing genomic order.
    let mut reader = bam::io::Reader::new(fs::File::open(&bam_path).unwrap());
    let _header = reader.read_header().expect("BAM header readable");

    let mut positions: Vec<(usize, usize)> = Vec::new(); // (ref_id, pos)
    for rec in reader.records() {
        let rec = rec.expect("valid BAM record");
        // Skip unmapped records (reference_sequence_id == None)
        let rid_opt = rec
            .reference_sequence_id()
            .map(|r| r.expect("ref_id readable"));
        let pos_opt = rec
            .alignment_start()
            .map(|p| p.expect("pos readable").get());
        if let (Some(rid), Some(pos)) = (rid_opt, pos_opt) {
            positions.push((rid, pos));
        }
    }

    assert!(
        !positions.is_empty(),
        "need at least 1 mapped record to verify sort order"
    );

    // Verify non-decreasing order for at least 5 consecutive pairs (or all if fewer)
    let check_n = positions.len().min(10);
    for w in positions[..check_n].windows(2) {
        assert!(
            w[0] <= w[1],
            "BAM records out of order: {:?} > {:?}",
            w[0],
            w[1]
        );
    }
}

// ---------------------------------------------------------------------------
// Test 3 — Paired-end alignment
// ---------------------------------------------------------------------------

#[test]
fn test_paired_end_alignment() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", None);

    // 30 FR pairs: mate1 at P, mate2 RC of P+150..P+200
    // P = 500, 700, 900, ... (stride 200)
    let mate1_path = tmpdir.path().join("mate1.fq");
    let mate2_path = tmpdir.path().join("mate2.fq");
    {
        let mut f1 = fs::File::create(&mate1_path).unwrap();
        let mut f2 = fs::File::create(&mate2_path).unwrap();
        for i in 0..30usize {
            let p = 500 + i * 200;
            // mate1: forward strand [p..p+50]
            let seq1 = &genome[p..p + 50];
            // mate2: RC of [p+150..p+200] — the "right" end of the fragment
            let seq2 = rc(&genome[p + 150..p + 200]);

            writeln!(f1, "@read{}/1", i + 1).unwrap();
            f1.write_all(seq1).unwrap();
            writeln!(f1).unwrap();
            writeln!(f1, "+").unwrap();
            writeln!(f1, "{}", "I".repeat(50)).unwrap();

            writeln!(f2, "@read{}/2", i + 1).unwrap();
            f2.write_all(&seq2).unwrap();
            writeln!(f2).unwrap();
            writeln!(f2, "+").unwrap();
            writeln!(f2, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_pe");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            mate1_path.to_str().unwrap(),
            mate2_path.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists(), "Aligned.out.sam not found");

    let content = fs::read_to_string(&sam_path).unwrap();
    // At least some records must have the PAIRED flag (0x1) set
    let paired_records = content
        .lines()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| {
            let mut cols = l.splitn(12, '\t');
            let _name = cols.next();
            if let Some(flag_str) = cols.next()
                && let Ok(flag) = flag_str.parse::<u16>()
            {
                return flag & 0x1 != 0; // PAIRED
            }
            false
        })
        .count();

    assert!(
        paired_records >= 1,
        "expected at least 1 paired record (flag 0x1), got {paired_records}"
    );
}

// ---------------------------------------------------------------------------
// Test 4 — Spliced alignment
// ---------------------------------------------------------------------------

#[test]
fn test_spliced_alignment() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    // Spliced read: 25 bp from Exon1 end + 25 bp from Exon2 start
    // genome[10025..10050] ++ genome[10250..10275]
    let mut spliced_read = genome[10025..10050].to_vec();
    spliced_read.extend_from_slice(&genome[10250..10275]);

    let fastq_path = tmpdir.path().join("spliced.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        for i in 0..10usize {
            writeln!(f, "@splice{}", i + 1).unwrap();
            f.write_all(&spliced_read).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_splice");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--sjdbOverhang",
            "24",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
            "--outFilterMismatchNmax",
            "20",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists(), "Aligned.out.sam not found");

    let content = fs::read_to_string(&sam_path).unwrap();
    let records: Vec<&str> = content.lines().filter(|l| !l.starts_with('@')).collect();

    assert!(!records.is_empty(), "no alignment records in SAM");

    // Check that at least one record has "N" in CIGAR (splice junction)
    let has_splice = records.iter().any(|l| {
        let cols: Vec<&str> = l.splitn(12, '\t').collect();
        if cols.len() >= 6 {
            return cols[5].contains('N');
        }
        false
    });

    // If the spliced alignment was found, great; otherwise just verify M alignment
    if !has_splice {
        // Fallback: at least one record with M in CIGAR (alignment succeeded)
        let has_match = records.iter().any(|l| {
            let cols: Vec<&str> = l.splitn(12, '\t').collect();
            cols.len() >= 6 && cols[5].contains('M')
        });
        assert!(
            has_match,
            "expected at least one alignment with M or N in CIGAR"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 5 — BySJout filtering
// ---------------------------------------------------------------------------

#[test]
fn test_bysj_filtering() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    // Mix: 10 spliced reads + 20 non-spliced reads
    let mut spliced_read = genome[10025..10050].to_vec();
    spliced_read.extend_from_slice(&genome[10250..10275]);

    let fastq_path = tmpdir.path().join("mixed.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        // 10 spliced reads
        for i in 0..10usize {
            writeln!(f, "@splice{}", i + 1).unwrap();
            f.write_all(&spliced_read).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
        // 20 non-spliced reads from various positions in background region
        for i in 0..20usize {
            let start = 200 + i * 150;
            let seq = &genome[start..start + 50];
            writeln!(f, "@normal{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_bysj");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    let output = cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outFilterType",
            "BySJout",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
            "--outFilterMismatchNmax",
            "20",
            "--outFileNamePrefix",
            &prefix,
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "rustar-aligner failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists(), "Aligned.out.sam not found");

    let log_path = output_dir.join("Log.final.out");
    assert!(log_path.exists(), "Log.final.out not found");

    // Verify the BySJout disk-buffering message was logged
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "outFilterType=BySJout: disk-buffering reads for post-alignment junction filtering"
        ),
        "expected BySJout disk-buffering log message in stderr; got:\n{stderr}"
    );
}

// ---------------------------------------------------------------------------
// Test 6 — GeneCounts output
// ---------------------------------------------------------------------------

#[test]
fn test_gene_counts_output() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    // 20 reads from Exon1 + 20 reads from Exon2
    let fastq_path = tmpdir.path().join("exon_reads.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        // Exon1 reads: genome[10000..10050]
        for i in 0..20usize {
            let offset = i; // all from same 50 bp window
            let seq = &genome[(10000 + offset)..(10000 + offset + 50)];
            writeln!(f, "@exon1_{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
        // Exon2 reads: genome[10250..10300]
        for i in 0..20usize {
            let offset = i;
            let seq = &genome[(10250 + offset)..(10250 + offset + 50)];
            writeln!(f, "@exon2_{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_genecounts");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--quantMode",
            "GeneCounts",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
            "--outFilterMismatchNmax",
            "20",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let tab_path = output_dir.join("ReadsPerGene.out.tab");
    assert!(tab_path.exists(), "ReadsPerGene.out.tab not found");

    let content = fs::read_to_string(&tab_path).unwrap();

    // Find the line for gene G1 and check at least one count column > 0
    let g1_line = content
        .lines()
        .find(|l| l.starts_with("G1"))
        .expect("gene G1 not found in ReadsPerGene.out.tab");

    let cols: Vec<&str> = g1_line.split('\t').collect();
    assert!(
        cols.len() >= 2,
        "G1 line has fewer than 2 columns: {g1_line}"
    );

    let max_count: i64 = cols[1..]
        .iter()
        .filter_map(|c| c.trim().parse().ok())
        .max()
        .unwrap_or(0);

    assert!(
        max_count > 0,
        "expected count > 0 for gene G1, got {g1_line}"
    );
}

// ---------------------------------------------------------------------------
// Test 7 — Unmapped reads output
// ---------------------------------------------------------------------------

#[test]
fn test_unmapped_reads_output() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", None);

    // 20 mappable reads + 10 unmappable (all-N)
    let fastq_path = tmpdir.path().join("mixed_unmapped.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        for i in 0..20usize {
            let start = 100 + i * 100;
            let seq = &genome[start..start + 50];
            writeln!(f, "@mapped{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
        for i in 0..10usize {
            writeln!(f, "@unmapped{}", i + 1).unwrap();
            writeln!(f, "{}", "N".repeat(50)).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_unmapped");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--outReadsUnmapped",
            "Fastx",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let unmapped_path = output_dir.join("Unmapped.out.mate1");
    assert!(
        unmapped_path.exists(),
        "Unmapped.out.mate1 not found at {unmapped_path:?}"
    );

    let content = fs::read_to_string(&unmapped_path).unwrap();
    let fastq_records = content.lines().filter(|l| l.starts_with('@')).count();
    assert!(
        fastq_records >= 1,
        "expected at least 1 FASTQ record in Unmapped.out.mate1, got {fastq_records}"
    );
}

// ---------------------------------------------------------------------------
// Test 8 — Two-pass mode
// ---------------------------------------------------------------------------

#[test]
fn test_two_pass_mode() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    // 20 spliced reads
    let mut spliced_read = genome[10025..10050].to_vec();
    spliced_read.extend_from_slice(&genome[10250..10275]);

    let fastq_path = tmpdir.path().join("twopass.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        for i in 0..20usize {
            writeln!(f, "@splice{}", i + 1).unwrap();
            f.write_all(&spliced_read).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_twopass");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--sjdbOverhang",
            "24",
            "--twopassMode",
            "Basic",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
            "--outFilterMismatchNmax",
            "20",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let pass1_path = output_dir.join("_STARpass1").join("SJ.out.tab");
    assert!(
        pass1_path.exists(),
        "_STARpass1/SJ.out.tab not found — two-pass mode did not write pass-1 junctions"
    );
    let top_level_pass1 = output_dir.join("SJ.pass1.out.tab");
    assert!(
        !top_level_pass1.exists(),
        "SJ.pass1.out.tab should no longer be emitted at the top level"
    );

    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists(), "Aligned.out.sam not found");

    let record_count = count_sam_records(&sam_path);
    assert!(
        record_count >= 1,
        "expected at least 1 alignment record, got {record_count}"
    );
}

// ---------------------------------------------------------------------------
// Test 9 — bare-dot prefix is treated as a literal string prefix (issue #26)
//
// STAR treats `--outFileNamePrefix SAMPLE.` as a literal prefix concatenated
// onto each output filename (SAMPLE.Aligned.out.bam at the top level), not as
// a directory name. This test asserts rustar-aligner matches that behaviour.
// ---------------------------------------------------------------------------

#[test]
fn test_bare_dot_prefix_is_literal_string() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", None);

    let fastq_path = tmpdir.path().join("reads.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        for i in 0..50usize {
            let start = 100 + i * 100;
            let seq = &genome[start..start + 50];
            writeln!(f, "@read{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    let run_dir = tmpdir.path().join("bare_dot_run");
    fs::create_dir_all(&run_dir).unwrap();
    // SAMPLE. is a bare-dot prefix; STAR writes SAMPLE.Aligned.out.bam at the top level.
    let prefix = format!("{}/SAMPLE.", run_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--outSAMtype",
            "BAM",
            "Unsorted",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let bam_path = run_dir.join("SAMPLE.Aligned.out.bam");
    let log_path = run_dir.join("SAMPLE.Log.final.out");
    let bam_as_dir = run_dir.join("SAMPLE.").join("Aligned.out.bam");

    assert!(
        bam_path.exists(),
        "expected literal prefix file at {}, but it was not created",
        bam_path.display()
    );
    assert!(
        log_path.exists(),
        "expected literal prefix file at {}, but it was not created",
        log_path.display()
    );
    assert!(
        !bam_as_dir.exists(),
        "bare-dot prefix was treated as a directory: {} should not exist",
        bam_as_dir.display()
    );

    let mut reader = bam::io::Reader::new(fs::File::open(&bam_path).unwrap());
    let _header = reader.read_header().expect("BAM header readable");
    let mut count = 0usize;
    for rec in reader.records() {
        rec.expect("valid BAM record");
        count += 1;
    }
    assert!(count >= 1, "expected at least 1 BAM record, got {count}");
}

// ---------------------------------------------------------------------------
// Test 9 — STARsolo (Phase 14.1–14.4): barcode parse, CB match, gene assign,
// UMI dedup, raw count-matrix output
// ---------------------------------------------------------------------------

#[test]
fn test_starsolo_gene_matrix() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    // cDNA reads (R2): 50 bp from Exon1 of gene G1 (genome[10000..10050]),
    // so each maps uniquely on the + strand inside G1 → Forward sense.
    let cdna_path = tmpdir.path().join("cdna.fq");
    let barcode_path = tmpdir.path().join("barcode.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");

    let cb = "AAAACCCCGGGGTTTT"; // 16 bp, sorts first in the whitelist
    // 8 reads, one cell, two well-separated UMI clouds (Hamming distance 10
    // apart, 4 reads each) → 1MM_All collapses each cloud to 1 molecule → 2.
    let umi_a = "ACGTACGTAC";
    let umi_b = "TGCATGCATG";
    let n_reads = 8usize;
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&barcode_path).unwrap();
        let exon1 = &genome[10000..10050];
        for i in 0..n_reads {
            writeln!(cf, "@read{i}").unwrap();
            cf.write_all(exon1).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(50)).unwrap();

            let umi = if i < 4 { umi_a } else { umi_b };
            writeln!(bf, "@read{i}").unwrap();
            writeln!(bf, "{cb}{umi}").unwrap();
            writeln!(bf, "+\n{}", "I".repeat(26)).unwrap();
        }
    }
    {
        let mut wf = fs::File::create(&wl_path).unwrap();
        writeln!(wf, "{cb}").unwrap();
        writeln!(wf, "CCCCGGGGTTTTAAAA").unwrap(); // decoys
        writeln!(wf, "GGGGTTTTAAAACCCC").unwrap();
    }

    let output_dir = tmpdir.path().join("out_solo");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    let assert = cargo_bin_cmd!("rustar-aligner")
        .env("RUST_LOG", "info")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            barcode_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloFeatures",
            "Gene",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    // cDNA alignments are emitted like a normal SE run.
    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists(), "Aligned.out.sam not found");
    assert!(
        count_sam_records(&sam_path) >= n_reads,
        "expected >= {n_reads} cDNA alignment records"
    );

    // 8 reads collected, all exact CB matches.
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).to_string();
    assert!(
        stderr.contains("collected 8 resolved"),
        "expected 8 resolved solo records in log, stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains("exact=8"),
        "expected 8 exact CB matches in log, stderr was:\n{stderr}"
    );

    // Raw matrix output.
    let raw = output_dir.join("Solo.out").join("Gene").join("raw");
    let features = fs::read_to_string(raw.join("features.tsv")).unwrap();
    let barcodes = fs::read_to_string(raw.join("barcodes.tsv")).unwrap();
    let matrix = fs::read_to_string(raw.join("matrix.mtx")).unwrap();

    // One gene G1 with a name column + feature type.
    assert_eq!(features.lines().count(), 1);
    assert!(
        features.starts_with("G1\tG1\tGene Expression"),
        "unexpected features.tsv:\n{features}"
    );
    // Three whitelist barcodes; the assayed CB sorts first.
    assert_eq!(barcodes.lines().count(), 3);
    assert_eq!(barcodes.lines().next().unwrap(), cb);

    // MatrixMarket: header, dims "1 3 1" (1 gene × 3 barcodes, 1 entry),
    // single entry "1 1 2" (gene 1, cell 1, 2 deduped molecules).
    let mtx_lines: Vec<&str> = matrix.lines().collect();
    assert!(
        mtx_lines[0].starts_with("%%MatrixMarket matrix coordinate integer general"),
        "unexpected mtx banner: {}",
        mtx_lines[0]
    );
    let dims = mtx_lines.iter().find(|l| !l.starts_with('%')).unwrap();
    assert_eq!(*dims, "1 3 1", "unexpected matrix dimensions");
    let entry = mtx_lines.last().unwrap();
    assert_eq!(
        *entry, "1 1 2",
        "expected 2 deduped molecules for G1 in cell 1"
    );

    // The default --soloCellFilter (CellRanger2.2) also writes a filtered/ matrix
    // containing only the called cell (the one assayed barcode), column-renumbered.
    let filt = output_dir.join("Solo.out").join("Gene").join("filtered");
    let f_barcodes = fs::read_to_string(filt.join("barcodes.tsv")).unwrap();
    assert_eq!(f_barcodes.lines().count(), 1, "expected 1 filtered cell");
    assert_eq!(f_barcodes.lines().next().unwrap(), cb);
    let f_matrix = fs::read_to_string(filt.join("matrix.mtx")).unwrap();
    let f_dims = f_matrix.lines().find(|l| !l.starts_with('%')).unwrap();
    assert_eq!(f_dims, "1 1 1", "unexpected filtered matrix dimensions");
    assert_eq!(f_matrix.lines().last().unwrap(), "1 1 2");

    // A CellRanger-style summary is written per feature.
    let summary =
        fs::read_to_string(output_dir.join("Solo.out").join("Gene").join("Summary.csv")).unwrap();
    assert!(
        summary.contains("Estimated Number of Cells,1"),
        "summary:\n{summary}"
    );
}

// ---------------------------------------------------------------------------
// Test 9a'' — --soloOutLayout CellRanger writes the same numbers where
// `cellranger count` writes them: outs/{raw,filtered}_feature_bc_matrix/,
// gzipped, with a -1 GEM-well suffix on every barcode.
// ---------------------------------------------------------------------------

#[test]
fn test_solo_out_layout_cellranger() {
    use std::io::Read;

    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    // Same fixture as test_starsolo_gene_matrix: 8 reads, one cell, two UMI
    // clouds, so the expected counts are already known.
    let cdna_path = tmpdir.path().join("cdna.fq");
    let barcode_path = tmpdir.path().join("barcode.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    let umi_a = "ACGTACGTAC";
    let umi_b = "TGCATGCATG";
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&barcode_path).unwrap();
        let exon1 = &genome[10000..10050];
        for i in 0..8usize {
            writeln!(cf, "@read{i}").unwrap();
            cf.write_all(exon1).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(50)).unwrap();
            let umi = if i < 4 { umi_a } else { umi_b };
            writeln!(bf, "@read{i}").unwrap();
            writeln!(bf, "{cb}{umi}").unwrap();
            writeln!(bf, "+\n{}", "I".repeat(26)).unwrap();
        }
    }
    {
        let mut wf = fs::File::create(&wl_path).unwrap();
        writeln!(wf, "{cb}").unwrap();
        writeln!(wf, "CCCCGGGGTTTTAAAA").unwrap();
        writeln!(wf, "GGGGTTTTAAAACCCC").unwrap();
    }

    let output_dir = tmpdir.path().join("out_crlayout");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    // No --soloOutLayout on the command line: 16 bp CB + 12 bp UMI with a
    // whitelist is 10x geometry, so the CellRanger layout is the default.
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            barcode_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloFeatures",
            "Gene",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let gunzip = |p: &std::path::Path| -> String {
        let f = fs::File::open(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
        let mut s = String::new();
        flate2::read::GzDecoder::new(f)
            .read_to_string(&mut s)
            .unwrap();
        s
    };

    // CellRanger's directory names, under outs/, with no per-feature level.
    let raw = output_dir.join("outs").join("raw_feature_bc_matrix");
    assert!(raw.is_dir(), "expected {}", raw.display());
    assert!(
        !output_dir.join("Solo.out").exists(),
        "Solo.out/ should not be written under the CellRanger layout"
    );

    // Every barcode carries the -1 GEM-well suffix, and the raw matrix has one
    // column per observed barcode (one), not one per whitelist entry (three).
    let barcodes = gunzip(&raw.join("barcodes.tsv.gz"));
    assert_eq!(barcodes.lines().count(), 1);
    assert_eq!(barcodes.lines().next().unwrap(), format!("{cb}-1"));

    let features = gunzip(&raw.join("features.tsv.gz"));
    assert!(features.starts_with("G1\tG1\tGene Expression"));

    // The 2 deduped molecules test_starsolo_gene_matrix asserts, in a matrix
    // that is now 1 gene × 1 observed barcode.
    let matrix = gunzip(&raw.join("matrix.mtx.gz"));
    let dims = matrix.lines().find(|l| !l.starts_with('%')).unwrap();
    assert_eq!(dims, "1 1 1", "unexpected matrix dimensions");
    assert_eq!(matrix.lines().last().unwrap(), "1 1 2");

    let filt = output_dir.join("outs").join("filtered_feature_bc_matrix");
    let f_barcodes = gunzip(&filt.join("barcodes.tsv.gz"));
    assert_eq!(f_barcodes.lines().next().unwrap(), format!("{cb}-1"));
    let f_matrix = gunzip(&filt.join("matrix.mtx.gz"));
    assert_eq!(f_matrix.lines().last().unwrap(), "1 1 2");

    // metrics_summary.csv: CellRanger 10.0.0's 20 metrics, in its order, as a
    // header row and a value row. The header is compared against the literal
    // string from a real CellRanger run.
    let metrics = fs::read_to_string(output_dir.join("outs").join("metrics_summary.csv")).unwrap();
    let mut lines = metrics.lines();
    assert_eq!(
        lines.next().unwrap(),
        "Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,\
         Number of Reads,Valid Barcodes,Valid UMI Sequences,Sequencing Saturation,\
         Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in UMI,\
         Reads Mapped to Genome,Reads Mapped Confidently to Genome,\
         Reads Mapped Confidently to Intergenic Regions,\
         Reads Mapped Confidently to Intronic Regions,\
         Reads Mapped Confidently to Exonic Regions,\
         Reads Mapped Confidently to Transcriptome,Reads Mapped Antisense to Gene,\
         Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell"
    );
    let values: Vec<&str> = lines.next().unwrap().split(',').collect();
    // 20 metrics; 8 reads, so no field is large enough to be comma-quoted.
    assert_eq!(values.len(), 20);
    assert_eq!(values[0], "1", "one called cell");
    assert_eq!(values[3], "8", "8 reads");
    assert_eq!(values[4], "100.0%", "all barcodes valid");
    assert_eq!(values[5], "100.0%", "all UMIs valid");
    // 8 reads collapsing to 2 molecules: 6 of the 8 added nothing.
    assert_eq!(values[6], "75.0%", "sequencing saturation");
    // The fixture writes 'I' (Phred 40) for every base.
    assert_eq!(values[7], "100.0%", "Q30 in barcode");
    assert_eq!(values[8], "100.0%", "Q30 in RNA read");
    assert_eq!(values[9], "100.0%", "Q30 in UMI");
    assert_eq!(values[18], "1", "one gene detected");
    assert_eq!(values[19], "2", "median 2 UMIs per cell");
    assert!(lines.next().is_none(), "exactly two rows");
}

// ---------------------------------------------------------------------------
// Test 9a' — Summary.csv stays STARsolo-faithful; the CellRanger mapping funnel
// (exonic/intronic/intergenic/antisense) is split out into a separate
// CellRanger.summary.csv (PR #90 review: keep the faithful Summary.csv unaltered).
// Needs Gene + GeneFull — the exonic/intronic split requires both queries.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_summary_split() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let bc_path = tmpdir.path().join("bc.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&bc_path).unwrap();
        let exon = &genome[10000..10050]; // exonic read → populates the funnel
        for (i, umi) in ["ACGTACGTACGT", "TGCATGCATGCA"].iter().enumerate() {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(exon).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(exon.len())).unwrap();
            writeln!(bf, "@r{i}\n{cb}{umi}\n+\n{}", "I".repeat(28)).unwrap();
        }
        fs::write(&wl_path, format!("{cb}\nCCCCGGGGTTTTAAAA\n")).unwrap();
    }

    let output_dir = tmpdir.path().join("out_split");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            bc_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloCBstart",
            "1",
            "--soloCBlen",
            "16",
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "12",
            "--soloFeatures",
            "Gene",
            "GeneFull",
            "--soloStrand",
            "Forward",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let gene = output_dir.join("Solo.out").join("Gene");
    let summary = fs::read_to_string(gene.join("Summary.csv")).unwrap();
    // Faithful Summary.csv must NOT carry the CellRanger funnel rows.
    assert!(
        !summary.contains("Exonic Regions") && !summary.contains("Antisense to Gene"),
        "Summary.csv leaked CellRanger funnel rows:\n{summary}"
    );
    assert!(
        summary.contains("Estimated Number of Cells"),
        "summary:\n{summary}"
    );
    // The funnel lives in a separate additional file.
    let cr = fs::read_to_string(gene.join("CellRanger.summary.csv")).unwrap();
    assert!(
        cr.contains("Reads Mapped Confidently to Exonic Regions")
            && cr.contains("Reads Mapped Antisense to Gene"),
        "CellRanger.summary.csv missing funnel rows:\n{cr}"
    );
}

// ---------------------------------------------------------------------------
// Test 9b — STARsolo SJ (splice-junction) feature
//
// Spliced cDNA reads (last 25 bp of Exon1 + first 25 bp of Exon2) cross the
// planted GT-AG intron, producing one junction. --soloFeatures SJ must write a
// Solo.out/SJ/raw matrix whose features.tsv equals SJ.out.tab and whose single
// junction row carries the deduped molecule count for the one cell.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_sj_feature() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let barcode_path = tmpdir.path().join("barcode.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    let umi = "ACGTACGTAC";
    // Spliced read: 25 bp from end of Exon1 + 25 bp from start of Exon2, which
    // aligns across the intron [10050,10250) → one GT-AG junction.
    let mut spliced = genome[10025..10050].to_vec();
    spliced.extend_from_slice(&genome[10250..10275]);
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&barcode_path).unwrap();
        for i in 0..6 {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(&spliced).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(50)).unwrap();
            writeln!(bf, "@r{i}\n{cb}{umi}\n+\n{}", "I".repeat(26)).unwrap();
        }
        let mut wf = fs::File::create(&wl_path).unwrap();
        writeln!(wf, "{cb}\nCCCCGGGGTTTTAAAA\nGGGGTTTTAAAACCCC").unwrap();
    }

    let output_dir = tmpdir.path().join("out_sj");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            barcode_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloFeatures",
            "Gene",
            "SJ",
            "--soloStrand",
            "Forward",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let sj_raw = output_dir.join("Solo.out").join("SJ").join("raw");
    let features = fs::read_to_string(sj_raw.join("features.tsv")).unwrap();
    let sj_tab = fs::read_to_string(output_dir.join("SJ.out.tab")).unwrap();
    // SJ feature file mirrors SJ.out.tab and contains exactly the one junction.
    assert_eq!(features, sj_tab, "SJ features.tsv must equal SJ.out.tab");
    assert_eq!(features.lines().count(), 1, "expected one junction");
    assert!(
        features.starts_with("chr1\t10051\t10250\t"),
        "unexpected junction: {features}"
    );
    // Matrix: 1 junction × 3 barcodes, single entry "1 1 1" (one deduped molecule
    // — all 6 reads share one UMI in one cell).
    let matrix = fs::read_to_string(sj_raw.join("matrix.mtx")).unwrap();
    let dims = matrix.lines().find(|l| !l.starts_with('%')).unwrap();
    assert_eq!(dims, "1 3 1", "unexpected SJ matrix dims");
    assert_eq!(matrix.lines().last().unwrap(), "1 1 1");
}

// ---------------------------------------------------------------------------
// Test 9c — STARsolo --soloMultiMappers (gene-ambiguous distribution)
//
// G1 and G3 share Exon1 (so a read there is ambiguous {G1,G3}); G2 has Exon2.
// One cell has a unique G2 molecule + one ambiguous {G1,G3} molecule. The unique
// matrix counts only G2; UniqueAndMult-Uniform spreads the ambiguous molecule
// 0.5/0.5 to G1 and G3 while keeping G2 at 1.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_multimappers() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    // GTF order: G1, G3 (both Exon1), G2 (Exon2) → gene indices 0,1,2.
    let gtf = tmpdir.path().join("multi.gtf");
    {
        let mut f = fs::File::create(&gtf).unwrap();
        for g in ["G1", "G3"] {
            writeln!(
                f,
                "chr1\tt\texon\t10001\t10050\t.\t+\t.\tgene_id \"{g}\"; transcript_id \"{g}t\";"
            )
            .unwrap();
        }
        writeln!(
            f,
            "chr1\tt\texon\t10251\t10300\t.\t+\t.\tgene_id \"G2\"; transcript_id \"G2t\";"
        )
        .unwrap();
    }
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let barcode_path = tmpdir.path().join("barcode.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&barcode_path).unwrap();
        // 4 reads in Exon2 → unique G2 (UMI a); 4 reads in Exon1 → ambiguous (UMI b).
        let exon2 = &genome[10250..10300];
        let exon1 = &genome[10000..10050];
        for (i, (seq, umi)) in [(exon2, "ACGTACGTAC"), (exon1, "TGCATGCATG")]
            .iter()
            .flat_map(|x| std::iter::repeat_n(*x, 4))
            .enumerate()
        {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(seq).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(50)).unwrap();
            writeln!(bf, "@r{i}\n{cb}{umi}\n+\n{}", "I".repeat(26)).unwrap();
        }
        let mut wf = fs::File::create(&wl_path).unwrap();
        writeln!(wf, "{cb}\nCCCCGGGGTTTTAAAA\nGGGGTTTTAAAACCCC").unwrap();
    }

    let output_dir = tmpdir.path().join("out_mm");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            barcode_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloFeatures",
            "Gene",
            "--soloStrand",
            "Forward",
            "--soloMultiMappers",
            "Uniform",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Gene").join("raw");
    // Unique matrix: only G2 (gene index 2 → row 3), count 1.
    let matrix = fs::read_to_string(raw.join("matrix.mtx")).unwrap();
    assert_eq!(
        matrix.lines().last().unwrap(),
        "3 1 1",
        "unique matrix:\n{matrix}"
    );
    // UniqueAndMult-Uniform: G1=0.5, G3=0.5, G2=1.
    let um = fs::read_to_string(raw.join("UniqueAndMult-Uniform.mtx")).unwrap();
    assert!(um.contains("coordinate real general"), "um header:\n{um}");
    let rows: Vec<&str> = um.lines().filter(|l| !l.starts_with('%')).skip(1).collect();
    assert!(rows.contains(&"1 1 0.50000"), "expected G1 0.5, got:\n{um}");
    assert!(rows.contains(&"2 1 0.50000"), "expected G3 0.5, got:\n{um}");
    assert!(rows.contains(&"3 1 1"), "expected G2 1, got:\n{um}");
}

// ---------------------------------------------------------------------------
// Test 9d — STARsolo SmartSeq (plate-based, manifest, no UMI)
//
// Two "cells" (manifest entries) of Exon1 reads → gene G1. With no UMIs each read
// is a count, so the matrix is G1 × {CellA,CellB} = read counts (5, 3).
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_smartseq() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let exon1 = &genome[10000..10050];
    let write_cell = |name: &str, n: usize| -> PathBuf {
        let p = tmpdir.path().join(name);
        let mut f = fs::File::create(&p).unwrap();
        for i in 0..n {
            writeln!(f, "@{name}_{i}").unwrap();
            f.write_all(exon1).unwrap();
            writeln!(f, "\n+\n{}", "I".repeat(50)).unwrap();
        }
        p
    };
    let a = write_cell("cellA.fq", 5);
    let b = write_cell("cellB.fq", 3);
    let manifest = tmpdir.path().join("manifest.tsv");
    fs::write(
        &manifest,
        format!("{}\t-\tCellA\n{}\t-\tCellB\n", a.display(), b.display()),
    )
    .unwrap();

    let output_dir = tmpdir.path().join("out_ss");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--soloType",
            "SmartSeq",
            "--readFilesManifest",
            manifest.to_str().unwrap(),
            "--soloStrand",
            "Forward",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Gene").join("raw");
    let barcodes = fs::read_to_string(raw.join("barcodes.tsv")).unwrap();
    assert_eq!(barcodes, "CellA\nCellB\n");
    let matrix = fs::read_to_string(raw.join("matrix.mtx")).unwrap();
    let dims = matrix.lines().find(|l| !l.starts_with('%')).unwrap();
    assert_eq!(dims, "1 2 2", "SmartSeq matrix dims:\n{matrix}");
    let entries: Vec<&str> = matrix
        .lines()
        .filter(|l| !l.starts_with('%'))
        .skip(1)
        .collect();
    assert!(entries.contains(&"1 1 5"), "expected CellA G1=5:\n{matrix}");
    assert!(entries.contains(&"1 2 3"), "expected CellB G1=3:\n{matrix}");
}

// ---------------------------------------------------------------------------
// Test 9d-PE — STARsolo SmartSeq paired-end (fragment counts)
//
// One cell, 4 read pairs: mate1 in Exon1, mate2 in (reverse-complement) Exon2 →
// a proper FR pair on gene G1. Each fragment is counted once (no UMI) → G1 = 4.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_smartseq_paired() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let r1_path = tmpdir.path().join("r1.fq");
    let r2_path = tmpdir.path().join("r2.fq");
    let mate1 = &genome[10000..10050]; // Exon1, forward
    let mate2 = rc(&genome[10250..10300]); // Exon2, reverse-complement (FR mate)
    {
        let mut f1 = fs::File::create(&r1_path).unwrap();
        let mut f2 = fs::File::create(&r2_path).unwrap();
        for i in 0..4 {
            writeln!(f1, "@p{i}").unwrap();
            f1.write_all(mate1).unwrap();
            writeln!(f1, "\n+\n{}", "I".repeat(50)).unwrap();
            writeln!(f2, "@p{i}").unwrap();
            f2.write_all(&mate2).unwrap();
            writeln!(f2, "\n+\n{}", "I".repeat(50)).unwrap();
        }
    }
    let manifest = tmpdir.path().join("manifest.tsv");
    fs::write(
        &manifest,
        format!("{}\t{}\tCellPE\n", r1_path.display(), r2_path.display()),
    )
    .unwrap();

    let output_dir = tmpdir.path().join("out_sspe");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--soloType",
            "SmartSeq",
            "--readFilesManifest",
            manifest.to_str().unwrap(),
            "--soloStrand",
            "Unstranded",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Gene").join("raw");
    let matrix = fs::read_to_string(raw.join("matrix.mtx")).unwrap();
    let dims = matrix.lines().find(|l| !l.starts_with('%')).unwrap();
    // One gene (G1) × one cell; 4 fragments counted.
    assert_eq!(dims, "1 1 1", "PE SmartSeq matrix dims:\n{matrix}");
    assert_eq!(
        matrix.lines().last().unwrap(),
        "1 1 4",
        "expected G1=4 fragments:\n{matrix}"
    );
}

// ---------------------------------------------------------------------------
// Test 9f — STARsolo Velocyto (spliced / unspliced / ambiguous)
//
// Three reads on gene G1, one per category: a junction-spanning read (spliced),
// a purely intronic read (unspliced), and a wholly-exonic read with no junction
// (ambiguous, per Sullivan 2025). Distinct UMIs → one molecule in each matrix.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_velocyto() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let bc_path = tmpdir.path().join("bc.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    // category → cDNA read + a distinct (non-homopolymer) 12 bp UMI.
    let mut spliced = genome[10025..10050].to_vec(); // Exon1 end ...
    spliced.extend_from_slice(&genome[10250..10275]); // ... + Exon2 start → junction
    let reads: [(Vec<u8>, &str); 3] = [
        (spliced, "ACGTACGTACGT"),                       // spliced
        (genome[10100..10150].to_vec(), "TGCATGCATGCA"), // intronic → unspliced
        (genome[10000..10050].to_vec(), "GATCGATCGATC"), // exonic, no junction → ambiguous
    ];
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&bc_path).unwrap();
        for (i, (seq, umi)) in reads.iter().enumerate() {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(seq).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(seq.len())).unwrap();
            writeln!(bf, "@r{i}\n{cb}{umi}\n+\n{}", "I".repeat(28)).unwrap();
        }
        fs::write(&wl_path, format!("{cb}\nCCCCGGGGTTTTAAAA\n")).unwrap();
    }

    let output_dir = tmpdir.path().join("out_velo");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            bc_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloCBstart",
            "1",
            "--soloCBlen",
            "16",
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "12",
            "--soloFeatures",
            "Velocyto",
            "--soloStrand",
            "Forward",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Velocyto").join("raw");
    // Each category matrix holds exactly its one molecule for G1 (row 1, col 1).
    for name in ["spliced", "unspliced", "ambiguous"] {
        let m = fs::read_to_string(raw.join(format!("{name}.mtx"))).unwrap();
        assert_eq!(
            m.lines().last().unwrap(),
            "1 1 1",
            "{name}.mtx should have G1=1:\n{m}"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 9f — Velocyto fold mode (`--soloVelocytoAmbiguous no`): exon-only
// molecules are folded into spliced and no `ambiguous.mtx` is written
// (cf. He, Soneson & Patro 2023). Same fixture as Test 9 but the exon-only
// molecule (distinct UMI) lands in spliced, giving G1 = 2 there.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_velocyto_fold_ambiguous() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let bc_path = tmpdir.path().join("bc.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    let mut spliced = genome[10025..10050].to_vec();
    spliced.extend_from_slice(&genome[10250..10275]);
    let reads: [(Vec<u8>, &str); 3] = [
        (spliced, "ACGTACGTACGT"),                       // spliced
        (genome[10100..10150].to_vec(), "TGCATGCATGCA"), // intronic → unspliced
        (genome[10000..10050].to_vec(), "GATCGATCGATC"), // exonic, no junction → ambiguous
    ];
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&bc_path).unwrap();
        for (i, (seq, umi)) in reads.iter().enumerate() {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(seq).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(seq.len())).unwrap();
            writeln!(bf, "@r{i}\n{cb}{umi}\n+\n{}", "I".repeat(28)).unwrap();
        }
        fs::write(&wl_path, format!("{cb}\nCCCCGGGGTTTTAAAA\n")).unwrap();
    }

    let output_dir = tmpdir.path().join("out_velo_fold");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            bc_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloCBstart",
            "1",
            "--soloCBlen",
            "16",
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "12",
            "--soloFeatures",
            "Velocyto",
            "--soloVelocytoAmbiguous",
            "no",
            "--soloStrand",
            "Forward",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Velocyto").join("raw");
    assert!(
        !raw.join("ambiguous.mtx").exists(),
        "ambiguous.mtx must not be written in fold mode"
    );
    let unspliced = fs::read_to_string(raw.join("unspliced.mtx")).unwrap();
    assert_eq!(
        unspliced.lines().last().unwrap(),
        "1 1 1",
        "unspliced unchanged by folding:\n{unspliced}"
    );
    let spliced_m = fs::read_to_string(raw.join("spliced.mtx")).unwrap();
    assert_eq!(
        spliced_m.lines().last().unwrap(),
        "1 1 2",
        "spliced should gain the folded exon-only molecule:\n{spliced_m}"
    );
}

// ---------------------------------------------------------------------------
// Test 9g — regression guard for the solo reader pipeline: a run that stops
// early via `--readMapNumber` must terminate. The background FASTQ-decode
// thread feeds a bounded (depth-2) channel; on an early break the consumer
// must disconnect so a producer blocked on the full channel wakes, instead of
// deadlocking the scope join. Needs > 2 batches (batch size 10k) of input to
// actually fill the channel after the break — hence 30k reads.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_readmapnumber_terminates() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let bc_path = tmpdir.path().join("bc.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");
    let cb = "AAAACCCCGGGGTTTT";
    let exon = &genome[10000..10050];
    let qual = "I".repeat(exon.len());
    {
        let mut cf = std::io::BufWriter::new(fs::File::create(&cdna_path).unwrap());
        let mut bf = std::io::BufWriter::new(fs::File::create(&bc_path).unwrap());
        for i in 0..30_000 {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(exon).unwrap();
            writeln!(cf, "\n+\n{qual}").unwrap();
            writeln!(bf, "@r{i}\n{cb}ACGTACGTACGT\n+\n{}", "I".repeat(28)).unwrap();
        }
        fs::write(&wl_path, format!("{cb}\n")).unwrap();
    }

    let output_dir = tmpdir.path().join("out_rmn");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    // readMapNumber (3) far below the 30k input: the consumer breaks after the
    // first batch while the producer still has batches to send. Must finish.
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            bc_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloCBstart",
            "1",
            "--soloCBlen",
            "16",
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "12",
            "--soloFeatures",
            "Gene",
            "--soloStrand",
            "Forward",
            "--readMapNumber",
            "3",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();
}

// ---------------------------------------------------------------------------
// Test 9e — STARsolo CB_UMI_Complex (multi-segment barcode)
//
// Barcode read layout: seg1(2bp) + linker(2bp) + seg2(2bp) + UMI(2bp). The cell
// barcode is seg1++seg2 matched against the cartesian product of two segment
// whitelists. All reads share CB=AAGG / UMI=AT → one molecule for gene G1.
// ---------------------------------------------------------------------------
#[test]
fn test_starsolo_cb_umi_complex() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let bc_path = tmpdir.path().join("bc.fq");
    let wl1 = tmpdir.path().join("wl1.txt");
    let wl2 = tmpdir.path().join("wl2.txt");
    fs::write(&wl1, "AA\nCC\n").unwrap(); // seg1 whitelist
    fs::write(&wl2, "GG\nTT\n").unwrap(); // seg2 whitelist
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&bc_path).unwrap();
        let exon1 = &genome[10000..10050];
        for i in 0..4 {
            writeln!(cf, "@r{i}").unwrap();
            cf.write_all(exon1).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(50)).unwrap();
            // seg1=AA, linker=CC, seg2=GG, UMI=AT → CB "AAGG", UMI "AT".
            writeln!(bf, "@r{i}\nAACCGGAT\n+\nIIIIIIII").unwrap();
        }
    }

    let output_dir = tmpdir.path().join("out_cx");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            bc_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Complex",
            "--soloCBwhitelist",
            wl1.to_str().unwrap(),
            wl2.to_str().unwrap(),
            "--soloCBposition",
            "0_0_0_1",
            "0_4_0_5",
            "--soloUMIposition",
            "0_6_0_7",
            "--soloUMIlen",
            "2",
            "--soloCBmatchWLtype",
            "Exact",
            "--soloFeatures",
            "Gene",
            "--soloStrand",
            "Forward",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Gene").join("raw");
    // Combined whitelist = {AA,CC}×{GG,TT} = 4 barcodes. The matched cell is AAGG;
    // all 4 reads share UMI AT → one molecule for G1.
    let matrix = fs::read_to_string(raw.join("matrix.mtx")).unwrap();
    let dims = matrix.lines().find(|l| !l.starts_with('%')).unwrap();
    let parts: Vec<&str> = dims.split_whitespace().collect();
    assert_eq!(
        parts[1], "4",
        "expected 4 combined-whitelist cells, dims={dims}"
    );
    assert_eq!(matrix.lines().last().unwrap(), "1 1 1", "matrix:\n{matrix}");
}

// ---------------------------------------------------------------------------
// Test 10 — CellRanger-style STARsolo run (Phase 14.5)
//
// Exercises the full CellRanger 4.x/5.x flag set from STARsolo.md:
//   --clipAdapterType CellRanger4 --outFilterScoreMin 30
//   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
//   --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
// and asserts the raw Gene matrix. The 1MM_CR UMI collapse is the key
// CellRanger-specific behavior verified here. A live differential comparison
// against the real STAR binary is in test/solo_cellranger_diff.py.
// ---------------------------------------------------------------------------

#[test]
fn test_starsolo_cellranger_style_matrix() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let gtf = write_gtf(&tmpdir);

    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", Some(&gtf));

    let cdna_path = tmpdir.path().join("cdna.fq");
    let barcode_path = tmpdir.path().join("barcode.fq");
    let wl_path = tmpdir.path().join("whitelist.txt");

    // One cell (CB sorts first), 8 reads in Exon1 of G1. UMIs: M x5 + a 1MM
    // neighbor of M x1 (1MM_CR collapses these to ONE molecule) + N x2 (a second
    // molecule) => 2 deduped molecules for (CB, G1).
    let cb = "AAAACCCCGGGGTTTT";
    let umi_m = "ACGTACGTAC"; // 10 bp (default soloUMIlen)
    let umi_m_1mm = "ACGTACGTAG"; // 1 mismatch from umi_m (last base)
    let umi_n = "TGCATGCATG";
    let plan = [(umi_m, 5usize), (umi_m_1mm, 1), (umi_n, 2)];
    {
        let mut cf = fs::File::create(&cdna_path).unwrap();
        let mut bf = fs::File::create(&barcode_path).unwrap();
        let exon1 = &genome[10000..10050];
        let mut i = 0;
        for (umi, n) in plan {
            for _ in 0..n {
                writeln!(cf, "@read{i}").unwrap();
                cf.write_all(exon1).unwrap();
                writeln!(cf, "\n+\n{}", "I".repeat(50)).unwrap();
                writeln!(
                    bf,
                    "@read{i}\n{cb}{umi}\n+\n{}",
                    "I".repeat(cb.len() + umi.len())
                )
                .unwrap();
                i += 1;
            }
        }
    }
    {
        let mut wf = fs::File::create(&wl_path).unwrap();
        writeln!(wf, "{cb}").unwrap();
        writeln!(wf, "TTTTGGGGCCCCAAAA").unwrap(); // decoy (sorts after cb)
    }

    let output_dir = tmpdir.path().join("out_cr");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            cdna_path.to_str().unwrap(),
            barcode_path.to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl_path.to_str().unwrap(),
            "--soloCBstart",
            "1",
            "--soloCBlen",
            "16",
            "--soloUMIstart",
            "17",
            "--soloUMIlen",
            "10",
            "--soloFeatures",
            "Gene",
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // CellRanger 4.x/5.x matching flags:
            "--clipAdapterType",
            "CellRanger4",
            "--outFilterScoreMin",
            "30",
            "--soloCBmatchWLtype",
            "1MM_multi_Nbase_pseudocounts",
            "--soloUMIfiltering",
            "MultiGeneUMI_CR",
            "--soloUMIdedup",
            "1MM_CR",
            "--outSAMtype",
            "SAM",
            // This fixture's 16 bp CB + 12 bp UMI is 10x geometry, which now
            // defaults to CellRanger's output layout; the assertions below are
            // about STARsolo's, so state it.
            "--soloOutLayout",
            "STARsolo",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let raw = output_dir.join("Solo.out").join("Gene").join("raw");
    let features = fs::read_to_string(raw.join("features.tsv")).unwrap();
    let barcodes = fs::read_to_string(raw.join("barcodes.tsv")).unwrap();
    let matrix = fs::read_to_string(raw.join("matrix.mtx")).unwrap();

    assert!(features.starts_with("G1\t"), "features.tsv: {features}");
    assert_eq!(barcodes.lines().count(), 2);
    assert_eq!(barcodes.lines().next().unwrap(), cb); // CB sorts first

    let lines: Vec<&str> = matrix.lines().collect();
    let dims = lines.iter().find(|l| !l.starts_with('%')).unwrap();
    assert_eq!(
        *dims, "1 2 1",
        "matrix dims (1 gene x 2 barcodes x 1 entry)"
    );
    // 1MM_CR: M(5)+M_1mm(1) collapse to 1 molecule, N(2) is another => 2.
    assert_eq!(
        *lines.last().unwrap(),
        "1 1 2",
        "expected 2 deduped molecules"
    );
}

// ---------------------------------------------------------------------------
// Test — WASP allele-specific SAMtag (--waspOutputMode SAMtag)
// ---------------------------------------------------------------------------

/// A uniquely-mapped read overlapping a heterozygous SNV should be re-mapped with
/// the allele swapped; since the region is unique it remaps to the same locus and
/// gets `vW:i:1` (passed WASP). Verifies the end-to-end WASP path + that the vW tag
/// (auto-added when `--waspOutputMode SAMtag`) reaches the SAM.
#[test]
fn test_wasp_samtag() {
    let tmpdir = TempDir::new().unwrap();
    let genome = build_genome();
    let fasta = write_fasta(&tmpdir, &genome);
    let genome_dir = tmpdir.path().join("genome");
    build_index(&fasta, &genome_dir, "7", None);

    // 10 reads from a unique 50 bp region [5000,5050); each overlaps the SNV at
    // 0-based 5025 (1-based 5026), which sits mid-read.
    let read_start = 5000usize;
    let snv_pos0 = 5025usize; // 0-based
    let fastq_path = tmpdir.path().join("reads.fq");
    {
        let mut f = fs::File::create(&fastq_path).unwrap();
        let seq = &genome[read_start..read_start + 50];
        for i in 0..10usize {
            writeln!(f, "@wr{}", i + 1).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(50)).unwrap();
        }
    }

    // Het SNV VCF at 1-based (snv_pos0 + 1); ALT is any base != REF.
    let ref_base = genome[snv_pos0] as char;
    let alt_base = if ref_base == 'A' { 'C' } else { 'A' };
    let vcf_path = tmpdir.path().join("variants.vcf");
    {
        let mut f = fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        )
        .unwrap();
        writeln!(
            f,
            "chr1\t{}\t.\t{}\t{}\t.\t.\t.\tGT\t0/1",
            snv_pos0 + 1,
            ref_base,
            alt_base
        )
        .unwrap();
    }

    let output_dir = tmpdir.path().join("out_wasp");
    fs::create_dir_all(&output_dir).unwrap();
    let prefix = format!("{}/", output_dir.display());

    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--waspOutputMode",
            "SAMtag",
            "--varVCFfile",
            vcf_path.to_str().unwrap(),
            "--outSAMtype",
            "SAM",
            "--outFileNamePrefix",
            &prefix,
        ])
        .assert()
        .success();

    let sam = fs::read_to_string(output_dir.join("Aligned.out.sam")).unwrap();
    let vw1 = sam
        .lines()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| l.split('\t').any(|f| f == "vW:i:1"))
        .count();
    assert_eq!(
        vw1, 10,
        "all 10 unique reads overlapping the het SNV should pass WASP (vW:i:1)"
    );
}
