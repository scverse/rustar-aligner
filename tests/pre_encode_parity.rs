//! Output parity for #223: records serialized on the align workers
//! (`RUSTAR_PRE_ENCODE=1`) must produce byte-identical files to records
//! serialized on the writer thread (`RUSTAR_PRE_ENCODE=0`).

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::Path;
use tempfile::TempDir;

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((state >> 16) & 3) as usize]
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
            _ => b'C',
        })
        .collect()
}

/// Two chromosomes; chr1 carries a GT-AG intron at [10050, 10250).
fn build_genome() -> (Vec<u8>, Vec<u8>) {
    let mut chr1 = lcg_seq(88888, 20000);
    chr1[10050..10052].copy_from_slice(b"GT");
    chr1[10248..10250].copy_from_slice(b"AG");
    (chr1, lcg_seq(4242, 12000))
}

fn write_fastq(path: &Path, reads: &[Vec<u8>]) {
    let mut f = fs::File::create(path).unwrap();
    for (i, seq) in reads.iter().enumerate() {
        // Varied qualities so --outSAMmode NoQS has something to strip.
        let qual: String = (0..seq.len())
            .map(|j| char::from(b"#5?FIJ"[(i + j) % 6]))
            .collect();
        writeln!(f, "@r{i}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f, "\n+\n{qual}").unwrap();
    }
}

struct Fixture {
    _tmp: TempDir,
    root: std::path::PathBuf,
}

fn fixture() -> Fixture {
    let tmp = TempDir::new().unwrap();
    let root = tmp.path().to_path_buf();
    let (chr1, chr2) = build_genome();

    let fasta = root.join("genome.fa");
    let mut f = fs::File::create(&fasta).unwrap();
    for (name, seq) in [("chr1", &chr1), ("chr2", &chr2)] {
        writeln!(f, ">{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f).unwrap();
    }
    let gtf = root.join("a.gtf");
    fs::write(
        &gtf,
        "chr1\tt\texon\t9901\t10050\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
         chr1\tt\texon\t10251\t10400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();

    let genome_dir = root.join("idx");
    fs::create_dir_all(&genome_dir).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeSAindexNbases", "7"])
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--genomeFastaFiles")
        .arg(&fasta)
        .arg("--sjdbGTFfile")
        .arg(&gtf)
        .args(["--sjdbOverhang", "49"])
        .arg("--outFileNamePrefix")
        .arg(genome_dir.join("run_"))
        .assert()
        .success();

    // SE: genomic (both strands), spliced across the intron, chimeric
    // (chr1 + chr2 halves) and random (unmapped) reads.
    let mut se = Vec::new();
    for i in 0..120usize {
        let s = 200 + i * 150;
        let r = chr1[s..s + 60].to_vec();
        se.push(if i % 2 == 0 { r } else { rc(&r) });
    }
    for k in 10..40usize {
        let mut r = chr1[10050 - k..10050].to_vec();
        r.extend_from_slice(&chr1[10250..10250 + 60 - k]);
        se.push(r);
    }
    for i in 0..20usize {
        let mut r = chr1[1000 + i * 300..1030 + i * 300].to_vec();
        r.extend_from_slice(&chr2[500 + i * 400..530 + i * 400]);
        se.push(r);
    }
    for i in 0..20u32 {
        se.push(lcg_seq(7 + i, 60));
    }
    write_fastq(&root.join("se.fq"), &se);

    // PE: 250 bp fragments, mate 2 reverse-complemented.
    let (mut m1, mut m2) = (Vec::new(), Vec::new());
    for i in 0..80usize {
        let (src, s) = if i % 3 == 0 {
            (&chr2, 100 + i * 120)
        } else {
            (&chr1, 100 + i * 200)
        };
        let frag = &src[s..s + 250];
        m1.push(frag[..60].to_vec());
        m2.push(rc(&frag[190..]));
    }
    // Chimeric pairs for --chimOutType WithinBAM: mate 1 on chr1, mate 2 on
    // chr2 (inter-mate), and mate 1 itself split across chr1 / chr2.
    for i in 0..30usize {
        let a = &chr1[12000 + i * 200..12000 + i * 200 + 250];
        let b = &chr2[6000 + i * 150..6000 + i * 150 + 250];
        if i % 2 == 0 {
            m1.push(a[..100].to_vec());
            m2.push(rc(&b[150..]));
        } else {
            let mut r = a[..50].to_vec();
            r.extend_from_slice(&b[..50]);
            m1.push(r);
            m2.push(rc(&b[150..]));
        }
    }
    write_fastq(&root.join("m1.fq"), &m1);
    write_fastq(&root.join("m2.fq"), &m2);

    Fixture { _tmp: tmp, root }
}

/// Run one configuration with pre-encoding forced on and off, from two
/// sibling directories with prefix `./` (the `@PG CL:` line records argv, so
/// both runs need identical arguments), and compare every output file.
fn assert_parity(fx: &Fixture, label: &str, args: &[&str]) {
    let mut outputs = Vec::new();
    for mode in ["0", "1"] {
        let dir = fx.root.join(format!("{label}_{mode}"));
        fs::create_dir_all(&dir).unwrap();
        let out = cargo_bin_cmd!("rustar-aligner")
            .current_dir(&dir)
            .env("RUSTAR_PRE_ENCODE", mode)
            .args(["--genomeDir", "../idx", "--outFileNamePrefix", "./"])
            .args(["--runThreadN", "3"])
            .args(args)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{label} (pre-encode {mode}) failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        fs::write(dir.join("stdout.bin"), &out.stdout).unwrap();
        outputs.push(dir);
    }
    let mut compared = 0;
    for entry in fs::read_dir(&outputs[0]).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap().to_string();
        if name.starts_with("Log.") || !path.is_file() {
            continue;
        }
        let a = fs::read(&path).unwrap();
        let b = fs::read(outputs[1].join(&name))
            .unwrap_or_else(|_| panic!("{label}: {name} missing with pre-encoding"));
        assert!(a == b, "{label}: {name} differs with pre-encoding");
        compared += 1;
    }
    assert!(compared >= 2, "{label}: nothing compared");
}

#[test]
fn pre_encoded_output_is_byte_identical() {
    let fx = fixture();
    let se = ["--readFilesIn", "../se.fq"];
    let pe = ["--readFilesIn", "../m1.fq", "../m2.fq"];
    let uw = ["--outSAMunmapped", "Within"];
    let cat = |parts: &[&[&'static str]]| -> Vec<&'static str> { parts.concat() };

    let configs: Vec<(&str, Vec<&str>)> = vec![
        ("se_sam", cat(&[&se, &uw])),
        (
            "se_bam",
            cat(&[&se, &uw, &["--outSAMtype", "BAM", "Unsorted"]]),
        ),
        (
            "se_sorted",
            cat(&[&se, &uw, &["--outSAMtype", "BAM", "SortedByCoordinate"]]),
        ),
        ("se_noqs", cat(&[&se, &uw, &["--outSAMmode", "NoQS"]])),
        ("se_stdout_sam", cat(&[&se, &["--outStd", "SAM"]])),
        (
            "se_stdout_bam",
            cat(&[
                &se,
                &uw,
                &[
                    "--outSAMtype",
                    "BAM",
                    "Unsorted",
                    "--outStd",
                    "BAM_Unsorted",
                ],
            ]),
        ),
        (
            "se_stdout_sorted",
            cat(&[
                &se,
                &[
                    "--outSAMtype",
                    "BAM",
                    "SortedByCoordinate",
                    "--outStd",
                    "BAM_SortedByCoordinate",
                ],
            ]),
        ),
        (
            "se_transcriptome",
            cat(&[
                &se,
                &[
                    "--quantMode",
                    "TranscriptomeSAM",
                    "--outSAMtype",
                    "BAM",
                    "Unsorted",
                ],
            ]),
        ),
        (
            "se_chimeric",
            cat(&[
                &se,
                &["--chimSegmentMin", "20", "--outSAMtype", "BAM", "Unsorted"],
            ]),
        ),
        (
            "se_unmapped_fastx",
            cat(&[&se, &["--outReadsUnmapped", "Fastx"]]),
        ),
        ("pe_sam", cat(&[&pe, &uw])),
        (
            "pe_sorted",
            cat(&[&pe, &uw, &["--outSAMtype", "BAM", "SortedByCoordinate"]]),
        ),
        (
            "pe_within_bam",
            cat(&[
                &pe,
                &[
                    "--chimSegmentMin",
                    "20",
                    "--chimOutType",
                    "Junctions",
                    "WithinBAM",
                    "--outSAMtype",
                    "BAM",
                    "Unsorted",
                ],
            ]),
        ),
        (
            "pe_transcriptome_noqs",
            cat(&[
                &pe,
                &[
                    "--quantMode",
                    "TranscriptomeSAM",
                    "--outSAMmode",
                    "NoQS",
                    "--outSAMtype",
                    "BAM",
                    "Unsorted",
                ],
            ]),
        ),
    ];
    for (label, args) in &configs {
        assert_parity(&fx, label, args);
    }
}
