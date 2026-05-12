/// Phase 9: Threading integration tests
use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use tempfile::TempDir;

/// Generate a pseudo-random genome sequence using an LCG PRNG.
/// Avoids repetitive patterns so seeds have manageable SA ranges
/// (compatible with bin-based windowing's anchor multimap filtering).
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

/// Helper to create a simple test genome (larger for proper indexing)
fn create_test_genome(dir: &TempDir) -> (PathBuf, String, String) {
    let fasta_path = dir.path().join("genome.fa");
    let mut file = fs::File::create(&fasta_path).unwrap();

    // Generate pseudo-random 1000bp chromosomes (unique 20-mers → SA ranges of 1-2)
    let chr1_seq = generate_genome_seq(12345, 1000);
    let chr2_seq = generate_genome_seq(67890, 1000);

    writeln!(file, ">chr1").unwrap();
    writeln!(file, "{}", chr1_seq).unwrap();
    writeln!(file, ">chr2").unwrap();
    writeln!(file, "{}", chr2_seq).unwrap();

    (fasta_path, chr1_seq, chr2_seq)
}

/// Helper to create test FASTQ with known reads extracted from genome sequences
fn create_test_fastq(dir: &TempDir, n_reads: usize, chr1_seq: &str, chr2_seq: &str) -> PathBuf {
    let fastq_path = dir.path().join("reads.fq");
    let mut file = fs::File::create(&fastq_path).unwrap();

    for i in 0..n_reads {
        writeln!(file, "@read{}", i + 1).unwrap();
        // Extract 20bp reads from known positions in the genome
        let read_start = (i * 7) % 980; // Stagger across the genome
        if i % 2 == 0 {
            writeln!(file, "{}", &chr1_seq[read_start..read_start + 20]).unwrap();
        } else {
            writeln!(file, "{}", &chr2_seq[read_start..read_start + 20]).unwrap();
        }
        writeln!(file, "+").unwrap();
        writeln!(file, "IIIIIIIIIIIIIIIIIIII").unwrap();
    }

    fastq_path
}

#[test]
fn test_single_thread_alignment() {
    let tmpdir = TempDir::new().unwrap();

    // Create genome
    let (fasta_path, chr1_seq, chr2_seq) = create_test_genome(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");

    // Generate genome index
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("genomeGenerate")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--genomeFastaFiles")
        .arg(&fasta_path)
        .arg("--genomeSAindexNbases")
        .arg("5")
        .assert()
        .success();

    // Create test reads
    let fastq_path = create_test_fastq(&tmpdir, 100, &chr1_seq, &chr2_seq);
    let output_dir = tmpdir.path().join("output_1t");

    // Align with 1 thread
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("alignReads")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--readFilesIn")
        .arg(&fastq_path)
        .arg("--runThreadN")
        .arg("1")
        .arg("--outFileNamePrefix")
        .arg(format!("{}/", output_dir.display()))
        .assert()
        .success()
        .stderr(predicate::str::contains("Alignment complete!"))
        .stderr(predicate::str::contains("Number of input reads: 100"));

    // Verify SAM output exists
    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists());

    // Verify SAM has correct number of lines (header + 100 records)
    let sam_content = fs::read_to_string(&sam_path).unwrap();
    let lines: Vec<&str> = sam_content.lines().collect();

    // Count header lines (@HD, @SQ, @PG)
    let header_lines = lines.iter().filter(|l| l.starts_with('@')).count();
    assert!(header_lines >= 3); // At least @HD, @SQ for chr1, @SQ for chr2, @PG

    // Count alignment lines (should be 100)
    let alignment_lines = lines.iter().filter(|l| !l.starts_with('@')).count();
    assert_eq!(alignment_lines, 100);
}

#[test]
fn test_multi_thread_alignment() {
    let tmpdir = TempDir::new().unwrap();

    // Create genome
    let (fasta_path, chr1_seq, chr2_seq) = create_test_genome(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");

    // Generate genome index
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("genomeGenerate")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--genomeFastaFiles")
        .arg(&fasta_path)
        .arg("--genomeSAindexNbases")
        .arg("5")
        .assert()
        .success();

    // Create test reads
    let fastq_path = create_test_fastq(&tmpdir, 100, &chr1_seq, &chr2_seq);
    let output_dir = tmpdir.path().join("output_4t");

    // Align with 4 threads
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("alignReads")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--readFilesIn")
        .arg(&fastq_path)
        .arg("--runThreadN")
        .arg("4")
        .arg("--outFileNamePrefix")
        .arg(format!("{}/", output_dir.display()))
        .assert()
        .success()
        .stderr(predicate::str::contains("Alignment complete!"))
        .stderr(predicate::str::contains("Number of input reads: 100"))
        .stderr(predicate::str::contains("Using 4 threads for alignment"));

    // Verify SAM output exists
    let sam_path = output_dir.join("Aligned.out.sam");
    assert!(sam_path.exists());

    // Verify SAM has correct number of lines
    let sam_content = fs::read_to_string(&sam_path).unwrap();
    let lines: Vec<&str> = sam_content.lines().collect();

    // Count alignment lines (should be 100)
    let alignment_lines = lines.iter().filter(|l| !l.starts_with('@')).count();
    assert_eq!(alignment_lines, 100);
}

#[test]
fn test_thread_count_consistency() {
    let tmpdir = TempDir::new().unwrap();

    // Create genome
    let (fasta_path, chr1_seq, chr2_seq) = create_test_genome(&tmpdir);
    let genome_dir = tmpdir.path().join("genome");

    // Generate genome index
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("genomeGenerate")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--genomeFastaFiles")
        .arg(&fasta_path)
        .arg("--genomeSAindexNbases")
        .arg("5")
        .assert()
        .success();

    // Create test reads
    let fastq_path = create_test_fastq(&tmpdir, 200, &chr1_seq, &chr2_seq);

    // Run with 1 thread
    let output_1t = tmpdir.path().join("output_1t");
    let result_1t = Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("alignReads")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--readFilesIn")
        .arg(&fastq_path)
        .arg("--runThreadN")
        .arg("1")
        .arg("--outFileNamePrefix")
        .arg(format!("{}/", output_1t.display()))
        .output()
        .unwrap();

    assert!(result_1t.status.success());

    // Run with 4 threads
    let output_4t = tmpdir.path().join("output_4t");
    let result_4t = Command::cargo_bin("rustar-aligner")
        .unwrap()
        .arg("--runMode")
        .arg("alignReads")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--readFilesIn")
        .arg(&fastq_path)
        .arg("--runThreadN")
        .arg("4")
        .arg("--outFileNamePrefix")
        .arg(format!("{}/", output_4t.display()))
        .output()
        .unwrap();

    assert!(result_4t.status.success());

    // Parse stats from both runs (logs go to stderr)
    let stderr_1t = String::from_utf8_lossy(&result_1t.stderr);
    let stderr_4t = String::from_utf8_lossy(&result_4t.stderr);

    // Both should process 200 reads
    assert!(stderr_1t.contains("Number of input reads: 200"));
    assert!(stderr_4t.contains("Number of input reads: 200"));

    // Extract mapped percentages (should be identical)
    let extract_percentage = |s: &str, pattern: &str| -> f64 {
        s.lines()
            .find(|l| l.contains(pattern))
            .and_then(|l| l.split('(').nth(1)?.split('%').next()?.trim().parse().ok())
            .unwrap_or(0.0)
    };

    let unique_1t = extract_percentage(&stderr_1t, "Uniquely mapped");
    let unique_4t = extract_percentage(&stderr_4t, "Uniquely mapped");

    // Percentages should match (allowing small floating point differences)
    assert!((unique_1t - unique_4t).abs() < 0.1);

    // SAM files should have same number of alignments
    let sam_1t = fs::read_to_string(output_1t.join("Aligned.out.sam")).unwrap();
    let sam_4t = fs::read_to_string(output_4t.join("Aligned.out.sam")).unwrap();

    let count_alignments = |s: &str| s.lines().filter(|l| !l.starts_with('@')).count();

    assert_eq!(count_alignments(&sam_1t), count_alignments(&sam_4t));
}

#[test]
fn test_batched_fastq_reading() {
    // Unit test for batched reading
    use rustar_aligner::io::fastq::FastqReader;
    use std::io::Write;

    let tmpdir = TempDir::new().unwrap();
    let fastq_path = tmpdir.path().join("test.fq");
    let mut file = fs::File::create(&fastq_path).unwrap();

    // Write 25 reads
    for i in 0..25 {
        writeln!(file, "@read{}", i + 1).unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "IIII").unwrap();
    }
    drop(file);

    // Open reader and read in batches of 10
    let mut reader = FastqReader::open(&fastq_path, None).unwrap();

    let batch1 = reader.read_batch(10).unwrap();
    assert_eq!(batch1.len(), 10);
    assert_eq!(batch1[0].name, "read1");
    assert_eq!(batch1[9].name, "read10");

    let batch2 = reader.read_batch(10).unwrap();
    assert_eq!(batch2.len(), 10);
    assert_eq!(batch2[0].name, "read11");
    assert_eq!(batch2[9].name, "read20");

    let batch3 = reader.read_batch(10).unwrap();
    assert_eq!(batch3.len(), 5); // Only 5 reads left
    assert_eq!(batch3[0].name, "read21");
    assert_eq!(batch3[4].name, "read25");

    let batch4 = reader.read_batch(10).unwrap();
    assert_eq!(batch4.len(), 0); // No more reads
}
