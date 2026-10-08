//! Each CellRanger-matching solo flag has to change the count matrix.
//!
//! Issue #172 measured two of them doing nothing at all while still being
//! accepted, which is the worst shape a flag can have: the run succeeds and
//! the numbers are someone else's. The existing CellRanger test cannot see
//! that, because its fixture has three non-zero matrix entries and none of the
//! populations these flags act on.
//!
//! These fixtures are built so that a flag *must* move a count: a UMI shared
//! between two genes with tied read support (which `MultiGeneUMI_CR` gives to
//! nobody), and a barcode carrying an `N` (which the `Nbase` match types
//! decide the fate of). If a flag goes inert again, a count changes here.

use assert_cmd::cargo::cargo_bin_cmd;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const CB: &str = "AAAACCCCGGGGTTTT";
const DECOY_CB: &str = "TTTTGGGGCCCCAAAA";
/// Two genes, far enough apart that no read spans both.
const G1: (usize, usize) = (2_000, 2_600);
const G2: (usize, usize) = (8_000, 8_600);
const READ_LEN: usize = 60;

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

fn write_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let genome = lcg_seq(88888, 20_000);
    let fasta = dir.join("genome.fa");
    let mut f = fs::File::create(&fasta).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(&genome).unwrap();
    writeln!(f).unwrap();

    let gtf = dir.join("genes.gtf");
    let mut g = fs::File::create(&gtf).unwrap();
    for (name, (s, e)) in [("G1", G1), ("G2", G2)] {
        writeln!(
            g,
            "chr1\tsyn\texon\t{}\t{}\t.\t+\t.\tgene_id \"{name}\"; transcript_id \"{name}_T1\";",
            s + 1,
            e
        )
        .unwrap();
    }
    (fasta, gtf)
}

fn build_index(fasta: &Path, gtf: &Path, genome_dir: &Path) {
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
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--sjdbOverhang",
            "59",
        ])
        .assert()
        .success();
}

/// Reads for one cell: `plan` gives (gene start, barcode, UMI, read count).
fn write_reads(dir: &Path, plan: &[(usize, String, &str, usize)]) {
    let genome = lcg_seq(88888, 20_000);
    let mut cf = fs::File::create(dir.join("cdna.fq")).unwrap();
    let mut bf = fs::File::create(dir.join("barcode.fq")).unwrap();
    let mut i = 0usize;
    for (gene_start, cb, umi, n) in plan {
        for k in 0..*n {
            // Stagger the start so the reads are not identical duplicates.
            let start = gene_start + 20 + (k % 5) * 7;
            writeln!(cf, "@read{i}").unwrap();
            cf.write_all(&genome[start..start + READ_LEN]).unwrap();
            writeln!(cf, "\n+\n{}", "I".repeat(READ_LEN)).unwrap();
            writeln!(
                bf,
                "@read{i}\n{cb}{umi}\n+\n{}",
                "I".repeat(cb.len() + umi.len())
            )
            .unwrap();
            i += 1;
        }
    }
    let mut wf = fs::File::create(dir.join("whitelist.txt")).unwrap();
    writeln!(wf, "{CB}").unwrap();
    writeln!(wf, "{DECOY_CB}").unwrap();
}

/// Run solo and return `gene_name -> count` from the raw matrix.
fn run_solo(
    dir: &Path,
    genome_dir: &Path,
    gtf: &Path,
    tag: &str,
    extra: &[&str],
) -> HashMap<String, u64> {
    let out = dir.join(format!("out_{tag}"));
    fs::create_dir_all(&out).unwrap();
    let prefix = format!("{}/", out.display());

    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.args([
        "--runMode",
        "alignReads",
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--readFilesIn",
        dir.join("cdna.fq").to_str().unwrap(),
        dir.join("barcode.fq").to_str().unwrap(),
        "--soloType",
        "CB_UMI_Simple",
        "--soloCBwhitelist",
        dir.join("whitelist.txt").to_str().unwrap(),
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
        "--outFileNamePrefix",
        &prefix,
    ]);
    cmd.args(extra);
    cmd.assert().success();

    let raw = out.join("Solo.out/Gene/raw");
    let features: Vec<String> = fs::read_to_string(raw.join("features.tsv"))
        .unwrap()
        .lines()
        .map(|l| l.split('\t').next().unwrap_or("").to_string())
        .collect();

    let mut counts: HashMap<String, u64> = HashMap::new();
    for line in fs::read_to_string(raw.join("matrix.mtx"))
        .unwrap()
        .lines()
        .skip(3)
    {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 {
            continue;
        }
        let gene = f[0].parse::<usize>().unwrap() - 1;
        let count = f[2].parse::<u64>().unwrap();
        *counts.entry(features[gene].clone()).or_insert(0) += count;
    }
    counts
}

#[test]
fn multi_gene_umi_cr_removes_a_tied_umi_from_the_matrix() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let (fasta, gtf) = write_fixture(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &gtf, &genome_dir);

    // One UMI shared between both genes with equal read support (STAR gives a
    // tie to nobody), plus one clean UMI per gene so the matrix is not empty
    // once the tied one is dropped.
    let shared = "ACGTACGTAC";
    let only1 = "TGCATGCATG";
    let only2 = "GGTTCCAAGG";
    write_reads(
        root,
        &[
            (G1.0, CB.to_string(), shared, 3),
            (G2.0, CB.to_string(), shared, 3),
            (G1.0, CB.to_string(), only1, 2),
            (G2.0, CB.to_string(), only2, 2),
        ],
    );

    let plain = run_solo(root, &genome_dir, &gtf, "plain", &[]);
    let cr = run_solo(
        root,
        &genome_dir,
        &gtf,
        "cr",
        &[
            "--soloUMIfiltering",
            "MultiGeneUMI_CR",
            "--soloUMIdedup",
            "1MM_CR",
        ],
    );

    let plain_total: u64 = plain.values().sum();
    let cr_total: u64 = cr.values().sum();
    assert!(
        plain_total > 0,
        "the fixture has to produce counts at all: {plain:?}"
    );
    assert!(
        cr_total < plain_total,
        "--soloUMIfiltering MultiGeneUMI_CR must drop the tied multi-gene UMI: \
         {plain_total} without it, {cr_total} with it ({plain:?} vs {cr:?})"
    );
    // Precisely: the tied UMI is the only multi-gene one, so it is the only
    // molecule that disappears, from each of the two genes.
    assert_eq!(
        plain_total - cr_total,
        2,
        "exactly the tied UMI's two molecules should go: {plain:?} vs {cr:?}"
    );
}

#[test]
fn an_n_containing_barcode_is_decided_by_the_cb_match_type() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let (fasta, gtf) = write_fixture(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &gtf, &genome_dir);

    // Half the reads carry a barcode with an `N` one substitution from the
    // whitelist entry; `Exact` cannot place them, the 1MM types can. The two
    // halves carry *different* UMIs on purpose: with a shared UMI a rescued
    // read would fold into the molecule the clean read already contributes,
    // and the totals would agree whether or not the flag did anything.
    let n_cb = format!("{}N{}", &CB[..8], &CB[9..]);
    write_reads(
        root,
        &[
            (G1.0, CB.to_string(), "ACGTACGTAC", 4),
            (G1.0, n_cb, "TGCATGCATG", 4),
        ],
    );

    let exact = run_solo(
        root,
        &genome_dir,
        &gtf,
        "exact",
        &["--soloCBmatchWLtype", "Exact"],
    );
    let onemm = run_solo(
        root,
        &genome_dir,
        &gtf,
        "1mm",
        &["--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts"],
    );

    let exact_total: u64 = exact.values().sum();
    let onemm_total: u64 = onemm.values().sum();
    println!("Exact total {exact_total}, 1MM_multi_Nbase_pseudocounts total {onemm_total}");
    assert!(
        onemm_total > exact_total,
        "an N-containing barcode cannot be placed by Exact but can be by the \
         1MM types: Exact {exact_total}, 1MM_multi_Nbase_pseudocounts {onemm_total}"
    );
    assert!(
        exact_total > 0,
        "the exact-barcode half of the fixture must still count: {exact:?}"
    );
}
