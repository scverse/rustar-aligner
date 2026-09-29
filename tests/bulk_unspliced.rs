//! Integration tests for the bulk total-RNA options (rustar-aligner
//! extensions, not in STAR):
//!
//! - `--quantMode GeneSplicing` (per-gene spliced / unspliced / ambiguous);
//! - `--quantTranscriptomeUnspliced Intron|PreMRNA` (`<gene_id>-I` targets
//!   in TranscriptomeSAM).
//!
//! A synthetic genome carries three genes, each with a fully spliced isoform
//! `A` and a retained-intron isoform `R` (exon 1 + intron 1 + exon 2 as one
//! exon). Reads are simulated from a known mixture of mature `A`, mature `R`
//! and unspliced pre-mRNA, as a reverse-stranded (dUTP) single-end library.

use assert_cmd::cargo::cargo_bin_cmd;
use noodles::bam;
use std::collections::HashMap;
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

/// Target names (in `@SQ` order) and, per read, the targets it was written
/// to in the transcriptome BAM.
/// `(target name, length)` in `@SQ` order, and read name -> target indices.
type TargetHits = (Vec<(String, usize)>, HashMap<String, Vec<usize>>);

fn read_targets(path: &Path) -> TargetHits {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let refs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(k, r)| {
            (
                String::from_utf8_lossy(k).to_string(),
                usize::from(r.length()),
            )
        })
        .collect();
    let mut per_read: HashMap<String, Vec<usize>> = HashMap::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        let q = String::from_utf8_lossy(rec.name().unwrap()).to_string();
        let tid = rec.reference_sequence_id().unwrap().unwrap();
        per_read.entry(q).or_default().push(tid);
    }
    (refs, per_read)
}

/// Reads per target estimated from the transcriptome BAM by a minimal
/// Salmon-like EM: a read compatible with targets `S` is shared in
/// proportion to `alpha_t / efflen_t`, `efflen_t = len_t - READ_LEN + 1`.
fn target_counts(path: &Path) -> HashMap<String, f64> {
    let (refs, per_read) = read_targets(path);
    let efflen: Vec<f64> = refs
        .iter()
        .map(|(_, l)| (l.saturating_sub(READ_LEN) + 1) as f64)
        .collect();
    let classes: Vec<Vec<usize>> = per_read.into_values().collect();
    let n = refs.len();
    let mut alpha = vec![1.0 / n as f64; n];
    let mut counts = vec![0.0; n];
    for _ in 0..2000 {
        counts.iter_mut().for_each(|c| *c = 0.0);
        for tids in &classes {
            let z: f64 = tids.iter().map(|&t| alpha[t] / efflen[t]).sum();
            if z > 0.0 {
                for &t in tids {
                    counts[t] += alpha[t] / efflen[t] / z;
                }
            }
        }
        let total: f64 = counts.iter().sum();
        for (a, c) in alpha.iter_mut().zip(&counts) {
            *a = c / total;
        }
    }
    refs.into_iter().map(|(n, _)| n).zip(counts).collect()
}

fn summary(path: &Path) -> HashMap<String, Vec<String>> {
    read(path)
        .lines()
        .skip(1)
        .map(|l| {
            let mut it = l.split('\t').map(str::to_string);
            (it.next().unwrap(), it.collect())
        })
        .collect()
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
    let splicing = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "splicing",
        &[
            "--quantMode",
            "GeneCounts",
            "TranscriptomeSAM",
            "GeneSplicing",
            "--quantTranscriptomeUnspliced",
            "None",
        ],
    );
    let unspliced = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "unspliced",
        &[
            "--quantMode",
            "GeneCounts",
            "TranscriptomeSAM",
            "--quantTranscriptomeUnspliced",
            "Intron",
        ],
    );

    let base_sam = sam_without_command_line(&base.join("Aligned.out.sam"));
    assert!(base_sam.lines().count() > 2000, "too few alignments");
    for other in [&splicing, &unspliced] {
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
    // GeneSplicing and --quantTranscriptomeUnspliced None leave the
    // transcriptome BAM untouched.
    let base_bam = base.join("Aligned.toTranscriptome.out.bam");
    assert_eq!(
        bam_records(&base_bam),
        bam_records(&splicing.join("Aligned.toTranscriptome.out.bam"))
    );
    // Intron adds the <gene_id>-I targets after the annotated transcripts.
    let (base_refs, _) = read_targets(&base_bam);
    let (refs, _) = read_targets(&unspliced.join("Aligned.toTranscriptome.out.bam"));
    let names: Vec<&str> = refs.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(refs[..base_refs.len()], base_refs[..]);
    assert_eq!(names[base_refs.len()..], ["G1-I", "G2-I", "G3-I"]);
    // New output files exist only when requested.
    assert!(splicing.join("ReadsPerGeneSplicing.out.tab").exists());
    assert!(!base.join("ReadsPerGeneSplicing.out.tab").exists());
    assert!(
        unspliced
            .join("Aligned.toTranscriptome.targets.tsv")
            .exists()
    );
    assert!(!base.join("Aligned.toTranscriptome.targets.tsv").exists());
    assert!(
        !unspliced
            .join("Aligned.toTranscriptome.unspliced.fa")
            .exists()
    );
}

#[test]
fn simulated_total_rna_mixture() {
    let (tmp, gdir, gtf, fq, pre_starts) = setup();
    let star = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "star",
        &["--quantMode", "TranscriptomeSAM", "GeneSplicing"],
    );
    let intron = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "intron",
        &[
            "--quantMode",
            "TranscriptomeSAM",
            "--quantTranscriptomeUnspliced",
            "Intron",
        ],
    );
    let premrna = align(
        &tmp,
        &gdir,
        &gtf,
        &fq,
        "premrna",
        &[
            "--quantMode",
            "TranscriptomeSAM",
            "--quantTranscriptomeUnspliced",
            "PreMRNA",
        ],
    );

    // --- Every simulated pre-mRNA read lying inside intron 1 (the intron
    // that isoform R retains) is written to the gene's unspliced target, and
    // no pre-mRNA read is written to a spliced target unless an isoform
    // contains it.
    let (refs, per_read) = read_targets(&intron.join("Aligned.toTranscriptome.out.bam"));
    let (i1s, i1e) = (EXONS[0].1, EXONS[1].0);
    let mut idx = 0;
    for gene in &GENES {
        for i in 0..gene.n_pre {
            let p = pre_starts[idx];
            idx += 1;
            if p >= i1s && p + READ_LEN <= i1e {
                let targets = &per_read[&format!("{}_P_{i}", gene.name)];
                let names: Vec<&str> = targets.iter().map(|&t| refs[t].0.as_str()).collect();
                assert!(
                    names.contains(&format!("{}-I", gene.name).as_str()),
                    "{names:?}"
                );
                assert!(
                    names.contains(&format!("{}_R", gene.name).as_str()),
                    "{names:?}"
                );
            }
        }
    }

    // --- EM estimates: share of the mature reads given to the retained-
    // intron isoform, and reads given to the spliced isoform A.
    let counts: Vec<HashMap<String, f64>> = [&star, &intron, &premrna]
        .iter()
        .map(|d| target_counts(&d.join("Aligned.toTranscriptome.out.bam")))
        .collect();
    let get = |c: &HashMap<String, f64>, t: String| c.get(&t).copied().unwrap_or(0.0);
    let mut share_err = [0.0f64; 3];
    let mut a_err = [0.0f64; 3];
    for gene in &GENES {
        let truth = gene.n_r as f64 / (gene.n_a + gene.n_r) as f64;
        let mut line = format!("{}: RI share truth {truth:.3}", gene.name);
        for (m, c) in counts.iter().enumerate() {
            let a = get(c, format!("{}_A", gene.name));
            let r = get(c, format!("{}_R", gene.name));
            let u = get(c, format!("{}-I", gene.name));
            share_err[m] += (r / (a + r) - truth).abs();
            a_err[m] += (a - gene.n_a as f64).abs() / gene.n_a as f64;
            line.push_str(&format!(
                " | {} A {a:.0} R {r:.0} -I {u:.0} share {:.3}",
                ["STAR", "Intron", "PreMRNA"][m],
                r / (a + r)
            ));
        }
        eprintln!(
            "{line} (simulated A {}, R {}, pre-mRNA {})",
            gene.n_a, gene.n_r, gene.n_pre
        );
    }
    eprintln!(
        "summed |RI share error|: STAR {:.3}, Intron {:.3}, PreMRNA {:.3}",
        share_err[0], share_err[1], share_err[2]
    );
    eprintln!(
        "summed relative error on A: STAR {:.3}, Intron {:.3}, PreMRNA {:.3}",
        a_err[0], a_err[1], a_err[2]
    );
    // Both modes stop the retained-intron isoform from absorbing pre-mRNA.
    for m in [1, 2] {
        assert!(
            share_err[m] < share_err[0] / 2.0,
            "RI share error {share_err:?}"
        );
    }
    // PreMRNA also improves the spliced isoform and recovers the pre-mRNA
    // reads on <gene_id>-I. Intron (splici-style) cannot: the exonic part of
    // pre-mRNA has no unspliced target to go to, so the spliced isoform
    // over-counts it, as it does with splici.
    assert!(a_err[2] < a_err[0], "A error {a_err:?}");
    for gene in &GENES {
        let u = get(&counts[2], format!("{}-I", gene.name));
        assert!(
            (u - gene.n_pre as f64).abs() < 0.15 * gene.n_pre as f64,
            "{}: -I {u}",
            gene.name
        );
    }

    // --- GeneSplicing: the library is reverse-stranded, so the forward
    // columns see nothing and the reverse columns see every gene read.
    let s = summary(&star.join("ReadsPerGeneSplicing.summary.tsv"));
    let n = |k: &str, col: usize| s[k][col].parse::<u64>().unwrap();
    let total_reads: usize = GENES.iter().map(|g| g.n_a + g.n_r + g.n_pre).sum();
    assert_eq!(
        n("N_spliced", 1) + n("N_unspliced", 1) + n("N_ambiguous", 1),
        0
    );
    let assigned_rev = n("N_spliced", 2) + n("N_unspliced", 2) + n("N_ambiguous", 2);
    assert!(assigned_rev as f64 > 0.95 * total_reads as f64);
    // A read is unspliced when some model needs pre-mRNA and none is only
    // exonic: here, exactly the pre-mRNA reads reaching more than 6 bases
    // (STAR's minOverlapMinusOne) into intron 2, which no isoform retains.
    let (i2s, i2e) = (EXONS[1].1, EXONS[2].0);
    let expected_unspliced = pre_starts
        .iter()
        .filter(|&&p| (p + READ_LEN).min(i2e).saturating_sub(p.max(i2s)) > 6)
        .count() as u64;
    let got = n("N_unspliced", 2);
    eprintln!("unspliced: expected {expected_unspliced}, got {got}");
    assert!(got.abs_diff(expected_unspliced) * 50 <= expected_unspliced);
    assert!(n("N_ambiguous", 2) > 0);
    let mature: usize = GENES.iter().map(|g| g.n_a).sum();
    assert!(n("N_spliced", 2) as usize >= mature * 95 / 100);
}
