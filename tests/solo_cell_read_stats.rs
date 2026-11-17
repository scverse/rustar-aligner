//! `--soloCellReadStats Standard` writes STARsolo's `CellReads.stats`, one table
//! per feature.
//!
//! The fixture is built so that the features disagree about a read: an
//! intronic read is a gene for `GeneFull` and nothing for `Gene`, which a single
//! table shared by both features cannot show. It also carries the reads STAR
//! sums into `CBnotInPasslist` (a barcode off the whitelist, a rejected UMI, a
//! correction `oneExact` drops) and a read on `chrM`, which STAR's default
//! `--genomeChrSetMitochondrial` (`chrM M MT`) counts as mitochondrial.

use assert_cmd::cargo::cargo_bin_cmd;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const CB: &str = "AAAACCCCGGGGTTTT";
const DECOY_CB: &str = "TTTTGGGGCCCCAAAA";
/// Two exons of G1 with an intron between them, and a one-exon G2.
const G1_EXON1: (usize, usize) = (2_000, 2_600);
const G1_EXON2: (usize, usize) = (3_400, 4_000);
const G1_INTRON: usize = 2_900;
const G2: (usize, usize) = (8_000, 8_600);
const INTERGENIC: usize = 12_000;
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

fn chr1() -> Vec<u8> {
    lcg_seq(88888, 20_000)
}

fn chr_m() -> Vec<u8> {
    lcg_seq(4242, 2_000)
}

fn write_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let fasta = dir.join("genome.fa");
    let mut f = fs::File::create(&fasta).unwrap();
    for (name, seq) in [("chr1", chr1()), ("chrM", chr_m())] {
        writeln!(f, ">{name}").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
    }

    let gtf = dir.join("genes.gtf");
    let mut g = fs::File::create(&gtf).unwrap();
    for (name, (s, e)) in [("G1", G1_EXON1), ("G1", G1_EXON2), ("G2", G2)] {
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

/// Reads: `(contig, start, barcode, UMI, count)`.
fn write_reads(dir: &Path, plan: &[(&[u8], usize, &str, &str, usize)]) {
    let mut cf = fs::File::create(dir.join("cdna.fq")).unwrap();
    let mut bf = fs::File::create(dir.join("barcode.fq")).unwrap();
    let mut i = 0usize;
    for (contig, start, cb, umi, n) in plan {
        for k in 0..*n {
            // Stagger the start so the reads are not identical duplicates.
            let s = start + (k % 5) * 7;
            writeln!(cf, "@read{i}").unwrap();
            cf.write_all(&contig[s..s + READ_LEN]).unwrap();
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
    for bc in [CB, DECOY_CB] {
        writeln!(wf, "{bc}").unwrap();
    }
}

/// One feature's file: barcode -> column -> value, plus the header.
type Table = HashMap<String, HashMap<String, String>>;

fn run(dir: &Path, genome_dir: &Path, tag: &str, extra: &[&str]) -> (Vec<String>, Table, Table) {
    let out = dir.join(format!("out_{tag}"));
    fs::create_dir_all(&out).unwrap();
    let prefix = format!("{}/", out.display());
    cargo_bin_cmd!("rustar-aligner")
        .args([
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
            "GeneFull",
            "--outFileNamePrefix",
            &prefix,
        ])
        .args(extra)
        .assert()
        .success();
    let load = |feature: &str| -> (Vec<String>, Table) {
        let text =
            fs::read_to_string(out.join(format!("Solo.out/{feature}/CellReads.stats"))).unwrap();
        let mut lines = text.lines();
        let header: Vec<String> = lines
            .next()
            .unwrap()
            .split('\t')
            .map(String::from)
            .collect();
        let rows = lines
            .map(|l| {
                let v: Vec<&str> = l.split('\t').collect();
                assert_eq!(v.len(), header.len(), "short row {l}");
                let row = header[1..]
                    .iter()
                    .cloned()
                    .zip(v[1..].iter().map(|s| (*s).to_string()))
                    .collect();
                (v[0].to_string(), row)
            })
            .collect();
        (header, rows)
    };
    let (header, gene) = load("Gene");
    let (header_full, full) = load("GeneFull");
    assert_eq!(header, header_full);
    (header, gene, full)
}

fn fixture() -> (TempDir, PathBuf) {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let (fasta, gtf) = write_fixture(root);
    let genome_dir = root.join("genome");
    build_index(&fasta, &gtf, &genome_dir);
    let (c1, cm) = (chr1(), chr_m());
    // One substitution from DECOY_CB, a barcode no read matches exactly.
    let near_decoy = format!("A{}", &DECOY_CB[1..]);
    write_reads(
        root,
        &[
            (&c1, G1_EXON1.0 + 20, CB, "ACGTACGTAC", 3),
            (&c1, G1_INTRON, CB, "TGCATGCATG", 2),
            (&c1, INTERGENIC, CB, "CCCCAAAAGT", 1),
            (&cm, 500, CB, "GATTACAGAT", 1),
            (&c1, G2.0 + 20, &near_decoy, "AACCGGTTAA", 2),
            (&c1, G2.0 + 20, "CCCCAAAAGGGGTTTT", "TTGGCCAATT", 1),
            (&c1, G2.0 + 20, CB, "ACGTNCGTAC", 1),
        ],
    );
    (dir, genome_dir)
}

fn get<'a>(t: &'a Table, cb: &str, col: &str) -> &'a str {
    &t[cb][col]
}

#[test]
fn each_feature_counts_its_own_flags() {
    let (dir, genome_dir) = fixture();
    let (header, gene, full) = run(
        dir.path(),
        &genome_dir,
        "std",
        &["--soloCellReadStats", "Standard"],
    );
    assert_eq!(header.len(), 20);
    assert_eq!(header[0], "CB");
    assert_eq!(header[1], "cbMatch");

    // STAR's rows: the passlist misses plus the one cell that resolved.
    for t in [&gene, &full] {
        let mut cbs: Vec<&str> = t.keys().map(String::as_str).collect();
        cbs.sort_unstable();
        assert_eq!(cbs, [CB, "CBnotInPasslist"]);
        // Every read lands in one row or the other.
        let total: u64 = t
            .values()
            .map(|r| r["cbMatch"].parse::<u64>().unwrap())
            .sum();
        assert_eq!(total, 11);
        assert_eq!(get(t, CB, "cbMatch"), "7");
        assert_eq!(get(t, CB, "cbPerfect"), "7");
        assert_eq!(get(t, CB, "genomeU"), "7");
        assert_eq!(get(t, CB, "mito"), "1", "chrM is mitochondrial by default");
        assert_eq!(
            get(t, CB, "intronic"),
            "0",
            "no feature sets an intronic ovType"
        );
    }

    // Gene: only the exonic reads are a feature, and `exonic` follows them.
    assert_eq!(get(&gene, CB, "featureU"), "3");
    assert_eq!(get(&gene, CB, "exonic"), "3");
    assert_eq!(get(&gene, CB, "countedU"), "3");
    assert_eq!(get(&gene, CB, "nUMIunique"), "1");
    // GeneFull: the intronic reads are G1's too, and GeneFull has no ovType.
    assert_eq!(get(&full, CB, "featureU"), "5");
    assert_eq!(get(&full, CB, "exonic"), "0");
    assert_eq!(get(&full, CB, "countedU"), "5");
    assert_eq!(get(&full, CB, "nUMIunique"), "2");

    // CBnotInPasslist: the two corrections oneExact drops, the barcode off the
    // whitelist, and the rejected UMI, all on G2, none counted.
    for t in [&gene, &full] {
        let miss = "CBnotInPasslist";
        assert_eq!(get(t, miss, "cbMatch"), "4");
        assert_eq!(get(t, miss, "cbPerfect"), "0");
        assert_eq!(get(t, miss, "cbMMunique"), "0");
        assert_eq!(get(t, miss, "genomeU"), "4");
        assert_eq!(get(t, miss, "featureU"), "4");
        assert_eq!(get(t, miss, "countedU"), "0");
    }
}

/// STAR's value is `Standard`; the `CB` rustar-aligner took before is kept as
/// an alias of it, and anything else is refused.
#[test]
fn standard_is_the_value_and_cb_an_alias() {
    let (dir, genome_dir) = fixture();
    let std = run(
        dir.path(),
        &genome_dir,
        "std",
        &["--soloCellReadStats", "Standard"],
    );
    let cb = run(
        dir.path(),
        &genome_dir,
        "cb",
        &["--soloCellReadStats", "CB"],
    );
    assert_eq!(std, cb);

    let out = dir.path().join("out_bad");
    fs::create_dir_all(&out).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            dir.path().join("cdna.fq").to_str().unwrap(),
            dir.path().join("barcode.fq").to_str().unwrap(),
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            dir.path().join("whitelist.txt").to_str().unwrap(),
            "--soloUMIlen",
            "10",
            "--soloCellReadStats",
            "Yes",
            "--outFileNamePrefix",
            &format!("{}/", out.display()),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("expected Standard or None"));
}

/// `--genomeChrSetMitochondrial -` names no mitochondrial chromosome.
#[test]
fn mito_follows_the_chromosome_set() {
    let (dir, genome_dir) = fixture();
    let (_, gene, _) = run(
        dir.path(),
        &genome_dir,
        "nomito",
        &[
            "--soloCellReadStats",
            "Standard",
            "--genomeChrSetMitochondrial",
            "-",
        ],
    );
    assert_eq!(get(&gene, CB, "mito"), "0");
}
