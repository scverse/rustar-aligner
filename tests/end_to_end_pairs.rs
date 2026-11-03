//! `--alignEndsType` on the two ends of a mate that meet in a pair. Expected values
//! are STAR 2.7.11b's observed output on this fixture (same genome and reads,
//! `--genomeSAindexNbases 7`).

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

fn genome() -> Vec<u8> {
    let mut state = 11u32;
    (0..20_000)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((state >> 16) & 3) as usize]
        })
        .collect()
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    }
}

/// `seq` with a mismatch at each of `positions`.
fn mutated(seq: &[u8], positions: &[usize]) -> Vec<u8> {
    let mut out = seq.to_vec();
    for &p in positions {
        out[p] = complement(out[p]);
    }
    out
}

const LAST_FOUR: [usize; 4] = [96, 97, 98, 99];
const FIRST_FOUR: [usize; 4] = [0, 1, 2, 3];

/// A pair 400 bases long at 2001, with 4 mismatches at each of the mates' `first`
/// and `second` positions on the forward strand.
fn pair(first: &[usize], second: &[usize]) -> (Vec<u8>, Vec<u8>) {
    let g = genome();
    (
        mutated(&g[2000..2100], first),
        reverse_complement(&mutated(&g[2300..2400], second)),
    )
}

/// FLAG, POS, CIGAR, AS, nM.
type Record = (u32, u32, String, i32, u32);

/// The alignment records of a run, and its `Log.final.out`.
fn run(reads: &[&[u8]], options: &[&str]) -> (Vec<Record>, String) {
    let fasta = format!(">c1\n{}\n", String::from_utf8_lossy(&genome()));
    run_on(&fasta, reads, options)
}

/// [`run`] on the genome in `fasta_text`.
fn run_on(fasta_text: &str, reads: &[&[u8]], options: &[&str]) -> (Vec<Record>, String) {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("g.fa");
    fs::write(&fasta, fasta_text).unwrap();
    let index = dir.path().join("i");
    fs::create_dir_all(&index).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "genomeGenerate", "--genomeSAindexNbases", "7"])
        .args(["--genomeDir", index.to_str().unwrap()])
        .args(["--genomeFastaFiles", fasta.to_str().unwrap()])
        .assert()
        .success();
    let files: Vec<String> = reads
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            let path = dir.path().join(format!("r{i}.fq"));
            let qual = "I".repeat(seq.len());
            fs::write(
                &path,
                format!("@r\n{}\n+\n{qual}\n", String::from_utf8_lossy(seq)),
            )
            .unwrap();
            path.to_str().unwrap().to_string()
        })
        .collect();
    let prefix = format!("{}/", dir.path().display());
    cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", index.to_str().unwrap()])
        .arg("--readFilesIn")
        .args(&files)
        .args(["--outFileNamePrefix", &prefix, "--outSAMtype", "SAM"])
        .args(options)
        .assert()
        .success();
    let sam = fs::read_to_string(Path::new(&prefix).join("Aligned.out.sam")).unwrap();
    let tag = |fields: &[&str], tag: &str| -> i64 {
        let found = fields.iter().find(|f| f.starts_with(tag)).unwrap();
        found[tag.len() + 3..].parse().unwrap()
    };
    let records = sam
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[1].parse::<u32>().unwrap() & 4 == 0)
        .map(|f| {
            (
                f[1].parse().unwrap(),
                f[3].parse().unwrap(),
                f[5].to_string(),
                tag(&f, "AS") as i32,
                tag(&f, "nM") as u32,
            )
        })
        .collect();
    let log = fs::read_to_string(Path::new(&prefix).join("Log.final.out")).unwrap();
    (records, log)
}

fn record(flag: u32, pos: u32, cigar: &str, score: i32, nm: u32) -> Record {
    (flag, pos, cigar.to_string(), score, nm)
}

fn end_to_end_pair(mate1: &[u8], mate2: &[u8]) -> Vec<Record> {
    run(&[mate1, mate2], &["--alignEndsType", "EndToEnd"]).0
}

#[test]
fn inner_ends_are_extended_and_charged() {
    let (m1, m2) = pair(&LAST_FOUR, &FIRST_FOUR);
    assert_eq!(
        end_to_end_pair(&m1, &m2),
        [
            record(99, 2001, "100M", 182, 8),
            record(147, 2301, "100M", 182, 8),
        ]
    );
    // The same pair with the mates swapped, on the other strand.
    assert_eq!(
        end_to_end_pair(&m2, &m1),
        [
            record(163, 2001, "100M", 182, 8),
            record(83, 2301, "100M", 182, 8),
        ]
    );
}

#[test]
fn outer_ends_are_charged_to_the_pair() {
    let (m1, m2) = pair(&FIRST_FOUR, &LAST_FOUR);
    assert_eq!(
        end_to_end_pair(&m1, &m2),
        [
            record(99, 2001, "100M", 182, 8),
            record(147, 2301, "100M", 182, 8),
        ]
    );
}

#[test]
fn a_mate_past_the_chromosome_start_cannot_be_clipped_into_a_pair() {
    let g = genome();
    let mut over_the_edge = b"ACGTT".to_vec();
    over_the_edge.extend_from_slice(&g[0..95]);
    let mate2 = reverse_complement(&g[300..400]);
    assert!(end_to_end_pair(&over_the_edge, &mate2).is_empty());

    // With filters low enough for one mate, only the mate that reaches its ends maps.
    let (records, _) = run(
        &[&over_the_edge, &mate2],
        &[
            "--alignEndsType",
            "EndToEnd",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
        ],
    );
    assert_eq!(records, [record(153, 301, "100M", 98, 0)]);
}

/// STAR checks the chromosome edge before the mate spacer (`extendAlign.cpp`),
/// so a mate whose 3' end sits on the edge cannot be extended into a pair.
#[test]
fn a_mate_ending_on_the_chromosome_end_cannot_be_extended_into_a_pair() {
    let g = genome();
    let at_end = &g[g.len() - 100..];
    let low = [
        "--outFilterScoreMinOverLread",
        "0.3",
        "--outFilterMatchNminOverLread",
        "0.3",
    ];
    for mode in ["EndToEnd", "Extend3pOfRead1"] {
        let options = [&["--alignEndsType", mode][..], &low[..]].concat();
        let (records, _) = run(&[at_end, &reverse_complement(at_end)], &options);
        assert_eq!(records, [record(153, 19_901, "100M", 98, 0)], "{mode}");
    }
}

#[test]
fn one_sided_modes_extend_only_their_inner_end() {
    let (m1, m2) = pair(&LAST_FOUR, &FIRST_FOUR);
    // Read 1's 3' end is forced, mate 2's 3' end stays clipped.
    let (records, _) = run(&[&m1, &m2], &["--alignEndsType", "Extend3pOfRead1"]);
    assert_eq!(
        records,
        [
            record(99, 2001, "100M", 186, 4),
            record(147, 2305, "4S96M", 186, 4),
        ]
    );
}

#[test]
fn single_end_read_past_the_chromosome_start_has_no_alignment() {
    let g = genome();
    let mut over_the_edge = b"ACGTT".to_vec();
    over_the_edge.extend_from_slice(&g[0..95]);
    let (records, log) = run(&[&over_the_edge], &["--alignEndsType", "EndToEnd"]);
    assert!(records.is_empty());
    let count = |label: &str| {
        log.lines()
            .find(|l| l.contains(label))
            .and_then(|l| l.split('|').nth(1))
            .unwrap()
            .trim()
            .to_string()
    };
    assert_eq!(count("Number of reads unmapped: other"), "1");
    assert_eq!(count("Number of reads unmapped: too short"), "0");
}

/// `--peOverlapNbasesMin` maps overlapping mates as one merged read. When that read
/// hangs over a chromosome edge, STAR keeps its end-to-end transcript at score 0
/// instead of dropping it (`stitchWindowAligns.cpp:44-55, 224`), and once another
/// window scores above 0 it splits every window transcript back into mates and
/// rescores them (`ReadAlign_peOverlapMergeMap.cpp:276-292`). The kept one then wins,
/// soft-clipped at the edge. Here the other window is a diverged copy on c2 of the
/// merged read's edge end.
#[test]
fn overlapping_mates_over_a_chromosome_edge_are_clipped_like_star() {
    let g = genome();
    let diverged = |s: &[u8]| -> Vec<u8> {
        let positions: Vec<usize> = (5..s.len()).step_by(10).collect();
        mutated(s, &positions)
    };
    let options = [
        "--alignEndsType",
        "EndToEnd",
        "--peOverlapNbasesMin",
        "10",
        "--outFilterScoreMinOverLread",
        "0.3",
        "--outFilterMatchNminOverLread",
        "0.3",
    ];
    for over in [&b"T"[..], b"GA"] {
        let k = over.len() as u32;
        let decoy_start = [over, &g[..40], &diverged(&g[40..150])].concat();
        let decoy_end = [
            &diverged(&g[g.len() - 150..g.len() - 40])[..],
            &g[g.len() - 40..],
            over,
        ]
        .concat();
        let c2 = [
            &g[5000..6000],
            &decoy_start,
            &g[6000..7000],
            &decoy_end,
            &g[7000..8000],
        ]
        .concat();
        let fasta = format!(
            ">c1\n{}\n>c2\n{}\n",
            String::from_utf8_lossy(&g),
            String::from_utf8_lossy(&c2)
        );
        for fragment in [130u32, 169] {
            let f = fragment as usize;
            // Over the start: mate 2 forward from the overhang, mate 1 reverse.
            let frag = [over, &g[..f - over.len()]].concat();
            let (records, _) = run_on(
                &fasta,
                &[&reverse_complement(&frag[f - 100..]), &frag[..100]],
                &options,
            );
            assert_eq!(
                records,
                [
                    record(163, 1, &format!("{k}S{}M", 100 - k), 198 - k as i32, 0),
                    record(83, fragment - 99 - k, "100M", 198 - k as i32, 0),
                ],
                "start, overhang {k}, fragment {fragment}"
            );
            // Over the end: mate 1 forward, mate 2 reverse into the overhang.
            let frag = [&g[g.len() - f + over.len()..], over].concat();
            let (records, _) = run_on(
                &fasta,
                &[&frag[..100], &reverse_complement(&frag[f - 100..])],
                &options,
            );
            assert_eq!(
                records,
                [
                    record(99, 20_001 - fragment + k, "100M", 198 - k as i32, 0),
                    record(
                        147,
                        19_901 + k,
                        &format!("{}M{k}S", 100 - k),
                        198 - k as i32,
                        0
                    ),
                ],
                "end, overhang {k}, fragment {fragment}"
            );
        }
    }
}
