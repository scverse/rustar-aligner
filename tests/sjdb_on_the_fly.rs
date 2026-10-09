//! `--sjdbGTFfile` at mapping time on an unannotated index (STAR's
//! sjdbInsertJunctions) must build exactly what `genomeGenerate` builds.

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

fn seq(seed: u32, length: usize) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            bases[((state >> 16) & 3) as usize]
        })
        .collect()
}

fn generate(dir: &std::path::Path, fasta: &std::path::Path, gtf: Option<&std::path::Path>) {
    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.args([
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        dir.to_str().unwrap(),
        "--genomeFastaFiles",
        fasta.to_str().unwrap(),
        "--genomeSAindexNbases",
        "5",
        "--sjdbOverhang",
        "30",
        "--outFileNamePrefix",
        dir.join("run_").to_str().unwrap(),
    ]);
    if let Some(g) = gtf {
        cmd.args(["--sjdbGTFfile", g.to_str().unwrap()]);
    }
    cmd.assert().success();
}

#[test]
fn on_the_fly_insertion_matches_genome_generate() {
    let tmp = TempDir::new().unwrap();
    let chr = seq(4242, 3000);
    let fasta = tmp.path().join("g.fa");
    writeln!(fs::File::create(&fasta).unwrap(), ">chr1\n{chr}").unwrap();
    let gtf = tmp.path().join("a.gtf");
    // Two spliced transcripts on opposite strands; random sequence makes the
    // introns non-canonical, so sjdbPrepare shifts them.
    let mut f = fs::File::create(&gtf).unwrap();
    for (s, e, st, t) in [
        (101, 400, "+", "T1"),
        (801, 1100, "+", "T1"),
        (1501, 1800, "-", "T2"),
        (2201, 2500, "-", "T2"),
    ] {
        writeln!(
            f,
            "chr1\tt\texon\t{s}\t{e}\t.\t{st}\t.\tgene_id \"G{t}\"; transcript_id \"{t}\";"
        )
        .unwrap();
    }
    drop(f);

    let annotated = tmp.path().join("idx_gtf");
    let plain = tmp.path().join("idx_plain");
    generate(&annotated, &fasta, Some(&gtf));
    generate(&plain, &fasta, None);

    let fq = tmp.path().join("r.fq");
    let read = format!("{}{}", &chr[370..400], &chr[800..830]);
    writeln!(
        fs::File::create(&fq).unwrap(),
        "@r1\n{read}\n+\n{}",
        "I".repeat(60)
    )
    .unwrap();
    let prefix = tmp.path().join("out_");
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--genomeDir",
            plain.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--sjdbOverhang",
            "30",
            "--sjdbInsertSave",
            "All",
            "--readFilesIn",
            fq.to_str().unwrap(),
            "--outFileNamePrefix",
            prefix.to_str().unwrap(),
        ])
        .assert()
        .success();

    let otf = tmp.path().join("out__STARgenome");
    for name in [
        "Genome",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.out.tab",
        "exonInfo.tab",
        "transcriptInfo.tab",
    ] {
        assert_eq!(
            fs::read(otf.join(name)).unwrap(),
            fs::read(annotated.join(name)).unwrap(),
            "{name} differs between on-the-fly insertion and genomeGenerate"
        );
    }
}
