//! A damaged index is refused at load, and junctions near a genome edge or
//! dropped by an on-the-fly insertion leave a working index.

use assert_cmd::cargo::cargo_bin_cmd;
use predicates::str::contains;
use std::fs;
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

fn write_genome(dir: &Path, seq: &[u8]) {
    let mut fa = b">c1\n".to_vec();
    fa.extend_from_slice(seq);
    fa.push(b'\n');
    fs::write(dir.join("g.fa"), fa).unwrap();
}

fn generate(dir: &Path, extra: &[&str]) -> assert_cmd::assert::Assert {
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeSAindexNbases", "8"])
        .arg("--genomeDir")
        .arg(dir.join("idx"))
        .arg("--genomeFastaFiles")
        .arg(dir.join("g.fa"))
        .arg("--outFileNamePrefix")
        .arg(dir.join("gg_"))
        .args(extra)
        .assert()
}

fn write_reads(dir: &Path, seq: &[u8]) {
    let mut out = String::new();
    for i in 0..50 {
        let s = std::str::from_utf8(&seq[100 + i * 40..200 + i * 40]).unwrap();
        out += &format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(100));
    }
    fs::write(dir.join("r.fq"), out).unwrap();
}

fn align(dir: &Path, extra: &[&str]) -> assert_cmd::assert::Assert {
    cargo_bin_cmd!("rustar-aligner")
        .arg("--genomeDir")
        .arg(dir.join("idx"))
        .arg("--readFilesIn")
        .arg(dir.join("r.fq"))
        .args(["--outSAMtype", "None"])
        .arg("--outFileNamePrefix")
        .arg(dir.join("a_"))
        .args(extra)
        .assert()
}

fn built_index() -> TempDir {
    let dir = TempDir::new().unwrap();
    let seq = lcg_seq(3, 5000);
    write_genome(dir.path(), &seq);
    write_reads(dir.path(), &seq);
    generate(dir.path(), &[]).success();
    align(dir.path(), &[]).success();
    dir
}

fn truncate(path: &Path, drop_bytes: u64) {
    let f = fs::OpenOptions::new().write(true).open(path).unwrap();
    f.set_len(f.metadata().unwrap().len() - drop_bytes).unwrap();
}

#[test]
fn short_sa_is_refused() {
    let dir = built_index();
    truncate(&dir.path().join("idx/SA"), 40);
    align(dir.path(), &[]).failure().stderr(contains("SA is"));
}

#[test]
fn short_saindex_is_refused() {
    let dir = built_index();
    truncate(&dir.path().join("idx/SAindex"), 40);
    align(dir.path(), &[])
        .failure()
        .stderr(contains("SAindex is"));
}

#[test]
fn garbage_chr_length_is_refused_without_a_panic() {
    let dir = built_index();
    fs::write(dir.path().join("idx/chrLength.txt"), "not-a-number\n").unwrap();
    align(dir.path(), &[])
        .failure()
        .stderr(contains("not a number"));
}

#[test]
fn short_chr_start_is_refused_without_a_panic() {
    let dir = built_index();
    fs::write(dir.path().join("idx/chrStart.txt"), "0\n").unwrap();
    align(dir.path(), &[]).failure().stderr(contains("starts"));
}

/// An exon at the very start of the genome, overhang 100: the donor flank runs
/// off the start. STAR reads its padding there and so builds the index.
#[test]
fn junction_at_the_genome_start_builds_and_maps() {
    let dir = TempDir::new().unwrap();
    let seq = lcg_seq(5, 5000);
    write_genome(dir.path(), &seq);
    write_reads(dir.path(), &seq);
    fs::write(dir.path().join("sj.tab"), "c1\t51\t499\t+\n").unwrap();
    let sj = dir.path().join("sj.tab");
    generate(
        dir.path(),
        &[
            "--sjdbFileChrStartEnd",
            sj.to_str().unwrap(),
            "--sjdbOverhang",
            "100",
        ],
    )
    .success();
    align(dir.path(), &[]).success();
}

/// Two listed junctions that are one repeat collapse when the index is reloaded,
/// so an on-the-fly insertion drops one old block from the SA. The SAindex has to
/// follow: it must equal the one a fresh build of the surviving junctions writes.
#[test]
fn insertion_that_drops_an_old_junction_keeps_a_valid_saindex() {
    let dir = TempDir::new().unwrap();
    let mut seq = lcg_seq(9, 40_000);
    // A GT..AG intron at 5001-6000, whose flanking base repeats, so that
    // 5000-5999 is the same junction shifted by one.
    seq[5000..5002].copy_from_slice(b"GT");
    seq[5998..6000].copy_from_slice(b"AG");
    seq[4999] = b'G';
    write_genome(dir.path(), &seq);
    write_reads(dir.path(), &seq);

    let old = dir.path().join("old.tab");
    fs::write(&old, "c1\t5001\t6000\t.\nc1\t5000\t5999\t+\n").unwrap();
    generate(
        dir.path(),
        &["--sjdbFileChrStartEnd", old.to_str().unwrap()],
    )
    .success();

    let new = dir.path().join("new.tab");
    fs::write(&new, "c1\t20001\t21000\t+\n").unwrap();
    align(
        dir.path(),
        &[
            "--sjdbFileChrStartEnd",
            new.to_str().unwrap(),
            "--sjdbInsertSave",
            "All",
        ],
    )
    .success();

    let fresh = TempDir::new().unwrap();
    write_genome(fresh.path(), &seq);
    let merged = fresh.path().join("merged.tab");
    fs::write(&merged, "c1\t5001\t6000\t+\nc1\t20001\t21000\t+\n").unwrap();
    generate(
        fresh.path(),
        &["--sjdbFileChrStartEnd", merged.to_str().unwrap()],
    )
    .success();

    let inserted = dir.path().join("a__STARgenome");
    for name in ["SA", "SAindex"] {
        assert!(
            fs::read(inserted.join(name)).unwrap()
                == fs::read(fresh.path().join("idx").join(name)).unwrap(),
            "{name} differs from a fresh build"
        );
    }
}

/// An index built with a GTF of one three-exon transcript: two junctions and
/// the transcript tables.
fn annotated_index() -> TempDir {
    let dir = TempDir::new().unwrap();
    let seq = lcg_seq(11, 5000);
    write_genome(dir.path(), &seq);
    write_reads(dir.path(), &seq);
    let gtf: String = [(1001, 1200), (2001, 2200), (3001, 3200)]
        .iter()
        .map(|(s, e)| {
            format!("c1\tt\texon\t{s}\t{e}\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n")
        })
        .collect();
    fs::write(dir.path().join("a.gtf"), gtf).unwrap();
    let gtf = dir.path().join("a.gtf");
    generate(
        dir.path(),
        &[
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--sjdbOverhang",
            "50",
        ],
    )
    .success();
    dir
}

/// A table cut mid-line still parses, with a wrong last value; STAR reads it
/// silently, except `transcriptInfo.tab`, so every one is refused.
#[test]
fn table_with_a_cut_last_line_is_refused() {
    let dir = annotated_index();
    let new = dir.path().join("new.tab");
    fs::write(&new, "c1\t4001\t4100\t+\n").unwrap();
    for name in [
        "sjdbInfo.txt",
        "sjdbList.out.tab",
        "transcriptInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "exonGeTrInfo.tab",
    ] {
        let path = dir.path().join("idx").join(name);
        let intact = fs::read(&path).unwrap();
        truncate(&path, 2);
        // sjdbList.out.tab is read only when junctions are inserted.
        align(
            dir.path(),
            &["--sjdbFileChrStartEnd", new.to_str().unwrap()],
        )
        .failure()
        .stderr(contains(format!("{name} is empty or its last line is cut")));
        fs::write(&path, intact).unwrap();
    }
    align(dir.path(), &[]).success();
}

/// STAR reloads the index's junctions from sjdbList.out.tab to insert new ones
/// and hangs when one is missing; without insertion it never reads the file.
#[test]
fn short_sjdb_list_is_refused_when_junctions_are_inserted() {
    let dir = annotated_index();
    let path = dir.path().join("idx/sjdbList.out.tab");
    let list = fs::read_to_string(&path).unwrap();
    fs::write(&path, list.lines().next().unwrap().to_string() + "\n").unwrap();
    align(dir.path(), &[]).success();
    align(dir.path(), &["--twopassMode", "Basic"])
        .failure()
        .stderr(contains("has 1 junctions but sjdbInfo.txt has 2"));
}

/// The junction count must fill the junction part of the Genome; STAR maps
/// reads wrongly when it does not.
#[test]
fn sjdb_info_that_disagrees_with_the_genome_is_refused() {
    let dir = annotated_index();
    let path = dir.path().join("idx/sjdbInfo.txt");
    let info = fs::read_to_string(&path).unwrap();
    fs::write(&path, info.replacen("\t50\n", "\t49\n", 1)).unwrap();
    align(dir.path(), &[])
        .failure()
        .stderr(contains("do not fill"));
}
