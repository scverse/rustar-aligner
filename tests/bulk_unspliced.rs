//! Integration tests for the bulk total-RNA options (rustar-aligner
//! extensions, not in STAR):
//!
//! - `--quantMode GeneVelocyto` (per-gene spliced / unspliced / ambiguous);
//! - `--quantTranscriptomePreMRNA BanRetainedIntron` (TranscriptomeSAM).
//!
//! A synthetic genome carries three genes, each with a fully spliced isoform
//! `A` and a retained-intron isoform `R` (exon 1 + intron 1 + exon 2 as one
//! exon). Reads are simulated from a known mixture of mature `A`, mature `R`
//! and unspliced pre-mRNA, as a reverse-stranded (dUTP) single-end library.

use assert_cmd::cargo::cargo_bin_cmd;
use noodles::bam;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const READ_LEN: usize = 50;
const CHR_LEN: usize = 40_000;

/// Gene layout relative to the gene start: exons `[0,300) [800,1000)
/// [1600,1900)`, introns `[300,800) [1000,1600)`.
const EXONS: [(usize, usize); 3] = [(0, 300), (800, 1000), (1600, 1900)];
const GENE_LEN: usize = 1900;

struct Gene {
    name: &'static str,
    start: usize,
    reverse: bool,
    /// Reads simulated from mature `A`, mature `R`, and pre-mRNA.
    n_a: usize,
    n_r: usize,
    n_pre: usize,
}

const GENES: [Gene; 3] = [
    Gene {
        name: "G1",
        start: 5_000,
        reverse: false,
        n_a: 300,
        n_r: 0,
        n_pre: 400,
    },
    Gene {
        name: "G2",
        start: 15_000,
        reverse: false,
        n_a: 300,
        n_r: 150,
        n_pre: 400,
    },
    Gene {
        name: "G3",
        start: 25_000,
        reverse: true,
        n_a: 300,
        n_r: 0,
        n_pre: 400,
    },
];

struct Lcg(u32);
impl Lcg {
    fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        self.0 >> 16
    }
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

fn build_genome() -> Vec<u8> {
    let mut rng = Lcg(424_242);
    let mut g: Vec<u8> = (0..CHR_LEN)
        .map(|_| b"ACGT"[(rng.next() & 3) as usize])
        .collect();
    // Plant canonical motifs at every intron: GT..AG on the transcribed
    // strand, i.e. CT..AC in forward genome coordinates for a - gene.
    for gene in &GENES {
        for w in EXONS.windows(2) {
            let (d, a) = (gene.start + w[0].1, gene.start + w[1].0);
            let (left, right) = if gene.reverse {
                (b"CT", b"AC")
            } else {
                (b"GT", b"AG")
            };
            g[d..d + 2].copy_from_slice(left);
            g[a - 2..a].copy_from_slice(right);
        }
    }
    g
}

/// Returns the FASTA, GTF and FASTQ paths, plus the gene-relative start of
/// every simulated pre-mRNA read.
fn write_inputs(dir: &Path, genome: &[u8]) -> (PathBuf, PathBuf, PathBuf, Vec<usize>) {
    let fasta = dir.join("genome.fa");
    let mut f = fs::File::create(&fasta).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(genome).unwrap();
    writeln!(f).unwrap();

    let gtf = dir.join("ann.gtf");
    let mut f = fs::File::create(&gtf).unwrap();
    for gene in &GENES {
        let strand = if gene.reverse { '-' } else { '+' };
        let a_exons = EXONS.to_vec();
        let r_exons = vec![(EXONS[0].0, EXONS[1].1), EXONS[2]];
        for (tr, exons) in [("A", a_exons), ("R", r_exons)] {
            for (s, e) in exons {
                writeln!(
                    f,
                    "chr1\tsim\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"{g}\"; transcript_id \"{g}_{tr}\";",
                    gene.start + s + 1,
                    gene.start + e,
                    g = gene.name
                )
                .unwrap();
            }
        }
    }

    // Reads: fragments of the molecule in genome-forward orientation; a
    // reverse-stranded read 1 is the reverse complement of the RNA, i.e. the
    // reverse complement of the fragment for a + gene and the fragment itself
    // for a - gene.
    let fq = dir.join("reads.fq");
    let mut f = fs::File::create(&fq).unwrap();
    let mut rng = Lcg(7);
    let mut pre_starts = Vec::new();
    for gene in &GENES {
        let body = [(0, GENE_LEN)];
        let r_exons = [(EXONS[0].0, EXONS[1].1), EXONS[2]];
        for (kind, blocks, n) in [
            ("A", &EXONS[..], gene.n_a),
            ("R", &r_exons[..], gene.n_r),
            ("P", &body[..], gene.n_pre),
        ] {
            let seq: Vec<u8> = blocks
                .iter()
                .flat_map(|&(s, e)| genome[gene.start + s..gene.start + e].iter().copied())
                .collect();
            for i in 0..n {
                let pos = rng.next() as usize % (seq.len() - READ_LEN + 1);
                if kind == "P" {
                    pre_starts.push(pos);
                }
                let frag = &seq[pos..pos + READ_LEN];
                let read = if gene.reverse {
                    frag.to_vec()
                } else {
                    rc(frag)
                };
                writeln!(f, "@{}_{kind}_{i}", gene.name).unwrap();
                f.write_all(&read).unwrap();
                writeln!(f, "\n+\n{}", "I".repeat(READ_LEN)).unwrap();
            }
        }
    }
    (fasta, gtf, fq, pre_starts)
}

/// Build the index once per test and return (tmpdir, genome dir, gtf,
/// fastq, pre-mRNA read starts).
fn setup() -> (TempDir, PathBuf, PathBuf, PathBuf, Vec<usize>) {
    let tmp = TempDir::new().unwrap();
    let genome = build_genome();
    let (fasta, gtf, fq, pre_starts) = write_inputs(tmp.path(), &genome);
    let gdir = tmp.path().join("idx");
    fs::create_dir_all(&gdir).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeDir"])
        .arg(&gdir)
        .arg("--genomeFastaFiles")
        .arg(&fasta)
        .args(["--genomeSAindexNbases", "8", "--sjdbOverhang", "49"])
        .arg("--sjdbGTFfile")
        .arg(&gtf)
        .arg("--outFileNamePrefix")
        .arg(gdir.join("gen_"))
        .assert()
        .success();
    (tmp, gdir, gtf, fq, pre_starts)
}

fn align(tmp: &TempDir, gdir: &Path, gtf: &Path, fq: &Path, name: &str, extra: &[&str]) -> PathBuf {
    let out = tmp.path().join(name);
    fs::create_dir_all(&out).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .arg("--genomeDir")
        .arg(gdir)
        .arg("--readFilesIn")
        .arg(fq)
        .arg("--sjdbGTFfile")
        .arg(gtf)
        .arg("--outFileNamePrefix")
        .arg(format!("{}/", out.display()))
        .args(extra)
        .assert()
        .success();
    out
}

/// SAM body plus header lines that do not embed the command line.
fn sam_without_command_line(path: &Path) -> String {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG") && !l.starts_with("@CO"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Every BAM record, in file order (the header carries the command line).
fn bam_records(path: &Path) -> Vec<String> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    reader.read_header().unwrap();
    reader
        .records()
        .map(|r| format!("{:?}", r.unwrap()))
        .collect()
}

/// `Log.final.out` minus the wall-clock lines.
fn log_final_without_times(path: &Path) -> String {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| {
            !(l.contains("Started") || l.contains("Finished") || l.contains("Mapping speed"))
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap()
}

#[test]
fn new_options_leave_existing_outputs_unchanged() {
    let (tmp, gdir, gtf, fq, _) = setup();
    let base = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "base",
        &["--quantMode", "GeneCounts", "TranscriptomeSAM"],
    );
    let velo = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "velo",
        &[
            "--quantMode",
            "GeneCounts",
            "TranscriptomeSAM",
            "GeneVelocyto",
            "--quantTranscriptomePreMRNA",
            "Keep",
        ],
    );
    let ban = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "ban",
        &[
            "--quantMode",
            "GeneCounts",
            "TranscriptomeSAM",
            "--quantTranscriptomePreMRNA",
            "BanRetainedIntron",
        ],
    );

    let base_sam = sam_without_command_line(&base.join("Aligned.out.sam"));
    assert!(base_sam.lines().count() > 2000, "too few alignments");
    for other in [&velo, &ban] {
        assert_eq!(
            base_sam,
            sam_without_command_line(&other.join("Aligned.out.sam"))
        );
        for f in ["ReadsPerGene.out.tab", "SJ.out.tab"] {
            assert_eq!(read(&base.join(f)), read(&other.join(f)), "{f} differs");
        }
        assert_eq!(
            log_final_without_times(&base.join("Log.final.out")),
            log_final_without_times(&other.join("Log.final.out"))
        );
    }
    // GeneVelocyto and PreMRNA Keep do not touch the transcriptome BAM ...
    let base_tr = bam_records(&base.join("Aligned.toTranscriptome.out.bam"));
    assert_eq!(
        base_tr,
        bam_records(&velo.join("Aligned.toTranscriptome.out.bam"))
    );
    // ... BanRetainedIntron removes records and adds none.
    let ban_tr = bam_records(&ban.join("Aligned.toTranscriptome.out.bam"));
    assert!(ban_tr.len() < base_tr.len());
    // The new output files exist only when requested.
    assert!(velo.join("ReadsPerGeneVelocyto.out.tab").exists());
    assert!(!base.join("ReadsPerGeneVelocyto.out.tab").exists());
    assert!(!ban.join("ReadsPerGeneVelocyto.summary.tsv").exists());
}
