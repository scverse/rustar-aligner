//! A pair with one aligned mate goes through the same selection and filters as a
//! pair with two. Expected values are STAR 2.7.11b's observed output on this
//! fixture (same genome and reads, `--genomeSAindexNbases 7`).

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

const COPIES: [usize; 3] = [5000, 10000, 15000];
const COPY_LEN: usize = 150;

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((state >> 16) & 3) as usize]
        })
        .collect()
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn genome() -> Vec<u8> {
    let mut g = lcg_seq(11, 20_900);
    let copy = lcg_seq(23, COPY_LEN);
    for start in COPIES {
        g[start..start + COPY_LEN].copy_from_slice(&copy);
    }
    g
}

/// A run of one pair, its `Aligned.out.sam` split into fields.
struct Run {
    records: Vec<Vec<String>>,
    log: String,
    genes: Option<String>,
}

fn run(mate1: &[u8], mate2: &[u8], options: &[&str]) -> Run {
    let dir = TempDir::new().unwrap();
    let write_read = |name: &str, seq: &[u8]| {
        let path = dir.path().join(name);
        let qual = "I".repeat(seq.len());
        fs::write(
            &path,
            format!("@r\n{}\n+\n{qual}\n", String::from_utf8_lossy(seq)),
        )
        .unwrap();
        path
    };
    let r1 = write_read("r1.fq", mate1);
    let r2 = write_read("r2.fq", mate2);
    let fasta = dir.path().join("g.fa");
    fs::write(
        &fasta,
        format!(">c1\n{}\n", String::from_utf8_lossy(&genome())),
    )
    .unwrap();
    let gtf = dir.path().join("g.gtf");
    fs::write(
        &gtf,
        "c1\tt\texon\t4001\t16000\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
    )
    .unwrap();
    let index = dir.path().join("i");
    fs::create_dir_all(&index).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeSAindexNbases", "7"])
        .args(["--genomeDir", index.to_str().unwrap()])
        .args(["--genomeFastaFiles", fasta.to_str().unwrap()])
        .args(["--sjdbGTFfile", gtf.to_str().unwrap()])
        .assert()
        .success();
    let prefix = format!("{}/", dir.path().display());
    cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", index.to_str().unwrap()])
        .args(["--readFilesIn", r1.to_str().unwrap(), r2.to_str().unwrap()])
        .args(["--outFileNamePrefix", &prefix, "--outSAMtype", "SAM"])
        .args(options)
        .assert()
        .success();
    let read = |name: &str| fs::read_to_string(Path::new(&prefix).join(name));
    Run {
        records: read("Aligned.out.sam")
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').map(String::from).collect())
            .collect(),
        log: read("Log.final.out").unwrap(),
        genes: read("ReadsPerGene.out.tab").ok(),
    }
}

impl Run {
    /// FLAG, POS, MAPQ, PNEXT of every record, in file order.
    fn core(&self) -> Vec<(u32, u32, u32, u32)> {
        self.records
            .iter()
            .map(|r| {
                let n = |i: usize| r[i].parse().unwrap();
                (n(1), n(3), n(4), n(7))
            })
            .collect()
    }

    fn tag(&self, record: usize, tag: &str) -> &str {
        let found = self.records[record]
            .iter()
            .find(|f| f.starts_with(tag))
            .unwrap();
        &found[tag.len() + 3..]
    }

    fn count(&self, label: &str) -> u64 {
        self.log
            .lines()
            .find_map(|l| {
                let (name, value) = l.split_once('|')?;
                (name.trim() == label).then(|| value.trim().parse().unwrap())
            })
            .unwrap()
    }
}

fn copy() -> Vec<u8> {
    lcg_seq(23, COPY_LEN)
}

fn no_mate() -> Vec<u8> {
    vec![b'N'; 15]
}

/// STAR on this fixture: the copies in the order 5001, 15001, 10001, MAPQ 1.
const ON_ALL_COPIES: [(u32, u32, u32, u32); 3] =
    [(73, 5001, 1, 0), (329, 15001, 1, 0), (329, 10001, 1, 0)];

#[test]
fn mate_in_a_repeat_is_a_multimapper_on_every_copy() {
    let run = run(&copy(), &no_mate(), &[]);
    assert_eq!(run.core(), ON_ALL_COPIES);
    for (i, hi) in ["1", "2", "3"].iter().enumerate() {
        assert_eq!(run.tag(i, "NH"), "3");
        assert_eq!(run.tag(i, "HI"), *hi);
        assert_eq!(run.tag(i, "AS"), "148");
    }
    assert_eq!(run.count("Uniquely mapped reads number"), 0);
    assert_eq!(run.count("Number of reads mapped to multiple loci"), 1);
}

#[test]
fn unmapped_mate_record_only_on_request() {
    let within = run(&copy(), &no_mate(), &["--outSAMunmapped", "Within"]);
    let mut expected = ON_ALL_COPIES.to_vec();
    expected.push((133, 0, 0, 5001));
    assert_eq!(within.core(), expected);
    assert_eq!(within.tag(3, "uT"), "4");

    let keep = run(
        &copy(),
        &no_mate(),
        &["--outSAMunmapped", "Within", "KeepPairs"],
    );
    assert_eq!(
        keep.core(),
        [
            (73, 5001, 1, 0),
            (133, 0, 0, 5001),
            (329, 15001, 1, 0),
            (389, 0, 0, 15001),
            (329, 10001, 1, 0),
            (389, 0, 0, 10001),
        ]
    );
}

#[test]
fn single_mates_are_held_to_the_loci_and_mismatch_limits() {
    let loci = run(&copy(), &no_mate(), &["--outFilterMultimapNmax", "2"]);
    assert!(loci.records.is_empty());
    assert_eq!(loci.count("Number of reads mapped to too many loci"), 1);

    let mut mutated = copy();
    for p in (40..120).step_by(30) {
        mutated[p] = if mutated[p] == b'A' { b'C' } else { b'A' };
    }
    let mismatches = run(&mutated, &no_mate(), &["--outFilterMismatchNmax", "1"]);
    assert!(mismatches.records.is_empty());
    // The mismatch limit also stops the extensions, so STAR calls this too short.
    assert_eq!(mismatches.count("Number of reads unmapped: too short"), 1);
}

#[test]
fn pair_and_single_mates_share_the_score_range() {
    let g = genome();
    let near = reverse_complement(&g[5250..5270]);
    let run = run(&copy(), &near, &["--outFilterMultimapScoreRange", "30"]);
    assert_eq!(
        run.core(),
        [
            (99, 5001, 1, 5251),
            (147, 5251, 1, 5001),
            (329, 15001, 1, 0),
            (329, 10001, 1, 0),
        ]
    );
    assert_eq!(run.tag(0, "AS"), "168");
    assert_eq!(run.tag(2, "AS"), "148");
    assert_eq!(run.tag(3, "HI"), "3");
}

#[test]
fn nh_counts_alignments_past_the_output_limit() {
    let g = genome();
    let near = reverse_complement(&g[5250..5270]);
    let run = run(
        &copy(),
        &near,
        &[
            "--outFilterMultimapScoreRange",
            "30",
            "--outSAMmultNmax",
            "2",
        ],
    );
    assert_eq!(
        run.core(),
        [
            (99, 5001, 1, 5251),
            (147, 5251, 1, 5001),
            (329, 15001, 1, 0),
        ]
    );
    for i in 0..3 {
        assert_eq!(run.tag(i, "NH"), "3");
    }
}

#[test]
fn mates_that_cannot_pair_are_two_alignments() {
    let g = genome();
    let run = run(
        &g[1000..1100],
        &g[2000..2100],
        &[
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
        ],
    );
    assert_eq!(run.core(), [(73, 1001, 3, 0), (393, 2001, 3, 0)]);
    assert_eq!(run.tag(1, "NH"), "2");
}

/// Not compared with STAR: its `--quantMode` fails on this machine. By the rule
/// that an alignment count above one is a multimapper.
#[test]
fn repeat_mate_counts_as_multimapping_in_gene_counts() {
    let run = run(&copy(), &no_mate(), &["--quantMode", "GeneCounts"]);
    let genes = run.genes.unwrap();
    assert!(genes.contains("N_multimapping\t1\t1\t1"), "{genes}");
    assert!(genes.contains("N_unmapped\t0\t0\t0"), "{genes}");
}
