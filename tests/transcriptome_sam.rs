//! Integration test for `--quantMode TranscriptomeSAM`.
//!
//! Builds a tiny genome + 2-transcript GTF + a small FASTQ, runs rustar-aligner
//! with `--quantMode TranscriptomeSAM`, and asserts that
//! `Aligned.toTranscriptome.out.bam` is produced, is a valid BAM file,
//! and contains at least one record.  Acts as a smoke test for the
//! end-to-end pipeline.

use assert_cmd::Command;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use tempfile::TempDir;

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

fn create_fasta(dir: &TempDir) -> (PathBuf, String) {
    let fasta_path = dir.path().join("genome.fa");
    let mut file = fs::File::create(&fasta_path).unwrap();
    // Single chromosome, 2000 bp pseudo-random content.
    let chr1_seq = generate_genome_seq(98765, 2000);
    writeln!(file, ">chr1").unwrap();
    writeln!(file, "{}", chr1_seq).unwrap();
    (fasta_path, chr1_seq)
}

fn create_gtf(dir: &TempDir) -> PathBuf {
    let gtf_path = dir.path().join("annotations.gtf");
    let mut file = fs::File::create(&gtf_path).unwrap();

    // Transcript T1: one exon [101, 400] forward (1-based inclusive)
    writeln!(
        file,
        "chr1\ttest\texon\t101\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
    )
    .unwrap();
    // Transcript T2: one exon [601, 900] forward
    writeln!(
        file,
        "chr1\ttest\texon\t601\t900\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";"
    )
    .unwrap();

    gtf_path
}

fn create_fastq(dir: &TempDir, n_reads: usize, chr1_seq: &str) -> PathBuf {
    let fastq_path = dir.path().join("reads.fq");
    let mut file = fs::File::create(&fastq_path).unwrap();

    // Alternate reads from T1 region and T2 region.
    for i in 0..n_reads {
        writeln!(file, "@read{}", i + 1).unwrap();
        let start = if i % 2 == 0 { 120 } else { 620 };
        let off = (i * 3) % 180;
        let s = start + off;
        writeln!(file, "{}", &chr1_seq[s..s + 30]).unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
    }

    fastq_path
}

/// Confirm the 5 STAR-compatible transcriptome index files are present
/// and structurally correct. Uses the test GTF (T1 [101,400) + /
/// T2 [601,900), both forward, G1/G2 single-exon transcripts) to derive
/// expected byte contents.
fn assert_star_transcriptome_files(genome_dir: &std::path::Path) {
    // transcriptInfo.tab: header = 2, then two lines in sorted order.
    // T1: trStart=100, trEnd=399 (inclusive), trEmax=399 (first), strand=1 (+),
    //     trExN=1, trExI=0, trGene=0.
    // T2: trStart=600, trEnd=899, trEmax=399 (running max excludes current),
    //     strand=1, trExN=1, trExI=1, trGene=1.
    let tr_info = fs::read_to_string(genome_dir.join("transcriptInfo.tab")).unwrap();
    assert_eq!(
        tr_info,
        "2\n\
         T1\t100\t399\t399\t1\t1\t0\t0\n\
         T2\t600\t899\t399\t1\t1\t1\t1\n",
        "transcriptInfo.tab byte format divergent from STAR"
    );

    // exonInfo.tab: header = 2, two single-exon transcripts with relative
    // [0, len-1] coords and exLenCum=0 for each first exon.
    let ex_info = fs::read_to_string(genome_dir.join("exonInfo.tab")).unwrap();
    assert_eq!(
        ex_info,
        "2\n\
         0\t299\t0\n\
         0\t299\t0\n",
        "exonInfo.tab byte format divergent from STAR"
    );

    // geneInfo.tab: header = 2, gene IDs in first-seen order. STAR falls
    // back to gene_id for missing gene_name, and the literal string
    // "MissingGeneType" for missing gene_biotype.
    let ge_info = fs::read_to_string(genome_dir.join("geneInfo.tab")).unwrap();
    assert_eq!(
        ge_info,
        "2\n\
         G1\tG1\tMissingGeneType\n\
         G2\tG2\tMissingGeneType\n",
        "geneInfo.tab byte format divergent from STAR"
    );

    // exonGeTrInfo.tab: header = 2, exons sorted by (start, end, strand, ...).
    // Both are forward (+strand == 1).
    let ge_tr_info = fs::read_to_string(genome_dir.join("exonGeTrInfo.tab")).unwrap();
    assert_eq!(
        ge_tr_info,
        "2\n\
         100\t399\t1\t0\t0\n\
         600\t899\t1\t1\t1\n",
        "exonGeTrInfo.tab byte format divergent from STAR"
    );

    // sjdbList.fromGTF.out.tab: both transcripts are single-exon, so no
    // splice junctions are produced. File must exist but be empty.
    let sj_path = genome_dir.join("sjdbList.fromGTF.out.tab");
    assert!(sj_path.exists(), "sjdbList.fromGTF.out.tab was not written");
    let sj = fs::read_to_string(&sj_path).unwrap();
    assert!(
        sj.is_empty(),
        "single-exon-only transcripts should yield empty sjdbList.fromGTF.out.tab, got: {:?}",
        sj
    );
}

#[test]
fn transcriptome_sam_end_to_end_smoke_test() {
    let tmpdir = TempDir::new().unwrap();

    let (fasta_path, chr1_seq) = create_fasta(&tmpdir);
    let gtf_path = create_gtf(&tmpdir);
    let fastq_path = create_fastq(&tmpdir, 20, &chr1_seq);

    let genome_dir = tmpdir.path().join("genome");
    let output_dir = tmpdir.path().join("output");
    fs::create_dir_all(&output_dir).unwrap();
    // STAR's prefix semantics: if prefix ends with `/`, it is treated as a
    // directory and files are placed inside it.  rustar-aligner's `PathBuf::join`
    // also handles the trailing slash convention.
    let output_prefix = format!("{}/", output_dir.display());

    // Build genome index, passing the GTF so transcriptInfo.tab + friends
    // are persisted (matches STAR's workflow).
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta_path.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf_path.to_str().unwrap(),
            "--genomeSAindexNbases",
            "5",
        ])
        .assert()
        .success();

    // Run alignment with --quantMode TranscriptomeSAM.
    Command::cargo_bin("rustar-aligner")
        .unwrap()
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--runThreadN",
            "1",
            "--outFileNamePrefix",
            output_prefix.as_str(),
            "--sjdbGTFfile",
            gtf_path.to_str().unwrap(),
            "--quantMode",
            "TranscriptomeSAM",
            // Permissive mismatch filter so our tiny read set gets through.
            "--outFilterMismatchNmax",
            "20",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
        ])
        .assert()
        .success();

    // All 5 STAR-compatible transcriptome index files must be written at
    // genomeGenerate time. Verify their byte formats match STAR exactly
    // (see src/quant/transcriptome.rs::write_*).
    assert_star_transcriptome_files(&genome_dir);

    // The transcriptome BAM must exist.
    let tr_bam = output_dir.join("Aligned.toTranscriptome.out.bam");
    assert!(
        tr_bam.exists(),
        "Aligned.toTranscriptome.out.bam was not created at {:?}",
        tr_bam
    );

    // File size sanity: BAM header + at least BGZF EOF marker > 100 bytes.
    let meta = fs::metadata(&tr_bam).unwrap();
    assert!(
        meta.len() > 100,
        "transcriptome BAM is suspiciously small: {} bytes",
        meta.len()
    );

    // Validate it's readable as a BAM and check the header has our
    // two @SQ lines (T1 + T2).
    use noodles::bam;
    use std::fs::File;
    let mut reader = bam::io::Reader::new(File::open(&tr_bam).unwrap());
    let header = reader.read_header().expect("valid BAM header");
    let refs = header.reference_sequences();
    assert_eq!(
        refs.len(),
        2,
        "expected 2 @SQ entries (T1, T2), got {}",
        refs.len()
    );
    let keys: Vec<String> = refs
        .keys()
        .map(|k| String::from_utf8_lossy(k).to_string())
        .collect();
    assert!(keys.iter().any(|k| k == "T1"), "T1 missing from @SQ");
    assert!(keys.iter().any(|k| k == "T2"), "T2 missing from @SQ");
}
