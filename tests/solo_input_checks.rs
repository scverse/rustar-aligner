//! Solo options that used to be accepted and then give an empty or wrong result
//! with exit 0, and the per-feature and gene-name outputs that depend on them.

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const CB: &str = "AAAACCCCGGGGTTTT";
const READ_LEN: usize = 60;
/// One gene with an intron: exons at 2001-2300 and 2901-3200.
const EXONIC: usize = 2_020;
const INTRONIC: usize = 2_500;

fn genome() -> Vec<u8> {
    let bases = *b"ACGT";
    let mut state = 12_345u32;
    (0..12_000)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            bases[((state >> 16) & 3) as usize]
        })
        .collect()
}

struct Fixture {
    dir: TempDir,
    genome_dir: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("g.fa");
        let mut f = fs::File::create(&fasta).unwrap();
        writeln!(f, ">chr1").unwrap();
        f.write_all(&genome()).unwrap();
        writeln!(f).unwrap();
        let genome_dir = dir.path().join("idx");
        fs::create_dir_all(&genome_dir).unwrap();
        cargo_bin_cmd!("rustar-aligner")
            .args(["--runMode", "genomeGenerate", "--genomeSAindexNbases", "5"])
            .arg("--genomeDir")
            .arg(&genome_dir)
            .arg("--genomeFastaFiles")
            .arg(&fasta)
            .arg("--outFileNamePrefix")
            .arg(dir.path().join("gen_"))
            .assert()
            .success();
        let gtf = dir.path().join("a.gtf");
        let attrs = r#"gene_id "G1"; transcript_id "T1"; gene_name "NAME"; gene_symbol "SYMBOL";"#;
        fs::write(
            &gtf,
            format!(
                "chr1\tt\texon\t2001\t2300\t.\t+\t.\t{attrs}\nchr1\tt\texon\t2901\t3200\t.\t+\t.\t{attrs}\n"
            ),
        )
        .unwrap();
        Self { dir, genome_dir }
    }

    /// Write `barcodes` paired with reads at `at` (0-based), and return the paths.
    fn reads(&self, at: &[usize], barcodes: &[String]) -> (PathBuf, PathBuf) {
        let g = genome();
        let (cdna, bc) = (self.dir.path().join("c.fq"), self.dir.path().join("b.fq"));
        let (mut fc, mut fb) = (
            fs::File::create(&cdna).unwrap(),
            fs::File::create(&bc).unwrap(),
        );
        for (i, (&p, b)) in at.iter().zip(barcodes).enumerate() {
            let seq = std::str::from_utf8(&g[p..p + READ_LEN]).unwrap();
            writeln!(fc, "@r{i}\n{seq}\n+\n{}", "I".repeat(READ_LEN)).unwrap();
            writeln!(fb, "@r{i}\n{b}\n+\n{}", "I".repeat(b.len())).unwrap();
        }
        (cdna, bc)
    }

    fn whitelist(&self, name: &str, barcode: &str) -> PathBuf {
        let p = self.dir.path().join(name);
        fs::write(&p, format!("{barcode}\n")).unwrap();
        p
    }

    fn run(&self, cdna: &Path, bc: &Path, args: &[&str]) -> assert_cmd::assert::Assert {
        cargo_bin_cmd!("rustar-aligner")
            .arg("--genomeDir")
            .arg(&self.genome_dir)
            .arg("--readFilesIn")
            .arg(cdna)
            .arg(bc)
            .arg("--sjdbGTFfile")
            .arg(self.dir.path().join("a.gtf"))
            .arg("--outFileNamePrefix")
            .arg(self.out())
            .args(args)
            .assert()
    }

    fn out(&self) -> PathBuf {
        self.dir.path().join("out/")
    }
}

fn stderr(a: &assert_cmd::assert::Assert) -> String {
    String::from_utf8_lossy(&a.get_output().stderr).into_owned()
}

/// Refused while parsing, so nothing needs to exist on disk.
fn refused(args: &[&str]) -> String {
    let a = cargo_bin_cmd!("rustar-aligner")
        .args(["--genomeDir", "none", "--readFilesIn", "a", "b"])
        .args(args)
        .assert()
        .failure();
    stderr(&a)
}

#[test]
fn solo_refuses_side_outputs_its_read_loops_do_not_write() {
    for (args, named) in [
        (&["--quantMode", "GeneCounts"][..], "GeneCounts"),
        (&["--quantMode", "TranscriptomeSAM"][..], "TranscriptomeSAM"),
        (&["--outReadsUnmapped", "Fastx"][..], "outReadsUnmapped"),
        (&["--outFilterType", "BySJout"][..], "BySJout"),
        (&["--chimSegmentMin", "12"][..], "chimSegmentMin"),
    ] {
        let mut all = vec!["--soloType", "CB_UMI_Simple"];
        all.extend_from_slice(args);
        let err = refused(&all);
        assert!(err.contains(named), "{args:?}: {err}");
    }
}

#[test]
fn umi_longer_than_a_u64_is_refused() {
    let err = refused(&["--soloType", "CB_UMI_Simple", "--soloUMIlen", "33"]);
    assert!(err.contains("--soloUMIlen 33"), "{err}");
}

#[test]
fn several_umi_dedup_methods_are_refused() {
    let err = refused(&[
        "--soloType",
        "CB_UMI_Simple",
        "--soloUMIdedup",
        "1MM_All",
        "Exact",
    ]);
    assert!(err.contains("--soloUMIdedup"), "{err}");
}

#[test]
fn complex_needs_a_valid_umi_position() {
    let base = [
        "--soloType",
        "CB_UMI_Complex",
        "--soloCBwhitelist",
        "wl",
        "--soloCBposition",
        "0_0_0_15",
    ];
    assert!(refused(&base).contains("--soloUMIposition"));
    let mut bad = base.to_vec();
    bad.extend_from_slice(&["--soloUMIposition", "0_16_nonsense"]);
    assert!(refused(&bad).contains("0_16_nonsense"));
    let mut wide = base.to_vec();
    wide.extend_from_slice(&["--soloUMIposition", "0_16_0_48"]);
    assert!(refused(&wide).contains("longer than 32"));
    let mut bad_cb = base[..4].to_vec();
    bad_cb.extend_from_slice(&["--soloCBposition", "oops", "--soloUMIposition", "0_16_0_25"]);
    assert!(refused(&bad_cb).contains("oops"));
}

#[test]
fn whitelist_of_another_length_is_fatal() {
    let fx = Fixture::new();
    let (c, b) = fx.reads(&[EXONIC], &[format!("{CB}AAAAACCCCC")]);
    let wl = fx.whitelist("wl14.txt", &CB[..14]);
    let a = fx
        .run(
            &c,
            &b,
            &[
                "--soloType",
                "CB_UMI_Simple",
                "--soloCBwhitelist",
                wl.to_str().unwrap(),
            ],
        )
        .failure();
    assert!(stderr(&a).contains("length 14"), "{}", stderr(&a));
}

#[test]
fn complex_segment_whitelist_of_another_length_is_fatal() {
    let fx = Fixture::new();
    let (c, b) = fx.reads(&[EXONIC], &[format!("{CB}AAAAACCCCC")]);
    let (w1, w2) = (
        fx.whitelist("a.txt", "AAAAAA"),
        fx.whitelist("b.txt", "CCCCCCCCCC"),
    );
    // 6 + 10 is the right total; each file is the wrong way round.
    let a = fx
        .run(
            &c,
            &b,
            &[
                "--soloType",
                "CB_UMI_Complex",
                "--soloCBwhitelist",
                w1.to_str().unwrap(),
                w2.to_str().unwrap(),
                "--soloCBposition",
                "0_0_0_9",
                "0_10_0_15",
                "--soloUMIposition",
                "0_16_0_25",
            ],
        )
        .failure();
    assert!(stderr(&a).contains("length 6"), "{}", stderr(&a));
}

#[test]
fn complex_umi_width_comes_from_its_position() {
    let fx = Fixture::new();
    let (c, b) = fx.reads(&[EXONIC], &[format!("{CB}ACGTACGT")]);
    let wl = fx.whitelist("wl.txt", CB);
    fx.run(
        &c,
        &b,
        &[
            "--soloType",
            "CB_UMI_Complex",
            "--soloCBwhitelist",
            wl.to_str().unwrap(),
            "--soloCBposition",
            "0_0_0_15",
            "--soloUMIposition",
            "0_16_0_23",
            "--soloBarcodeReadLength",
            "0",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outSAMattributes",
            "NH",
            "CB",
            "UB",
        ],
    )
    .success();
    let bam = fs::File::open(fx.out().join("Aligned.sortedByCoord.out.bam")).unwrap();
    let mut decoded = Vec::new();
    flate2::read::MultiGzDecoder::new(bam)
        .read_to_end(&mut decoded)
        .unwrap();
    let tag = b"UBZACGTACGT\0";
    assert!(
        decoded.windows(tag.len()).any(|w| w == tag),
        "the UB tag must be the eight bases of --soloUMIposition"
    );
}

/// An intronic read counts to GeneFull and not to Gene, so the two features'
/// files must disagree about it.
#[test]
fn cell_reads_stats_are_per_feature() {
    let fx = Fixture::new();
    let (c, b) = fx.reads(
        &[EXONIC, INTRONIC],
        &[format!("{CB}AAAAACCCCC"), format!("{CB}GGGGGTTTTT")],
    );
    let wl = fx.whitelist("wl.txt", CB);
    fx.run(
        &c,
        &b,
        &[
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl.to_str().unwrap(),
            "--soloFeatures",
            "Gene",
            "GeneFull",
            "--soloCellReadStats",
            "CB",
        ],
    )
    .success();
    let counted_u = |feature: &str| -> String {
        let text = fs::read_to_string(fx.out().join(format!("Solo.out/{feature}/CellReads.stats")))
            .unwrap();
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let col = header.iter().position(|h| *h == "countedU").unwrap();
        let row = lines.find(|l| l.starts_with(CB)).unwrap();
        row.split('\t').nth(col).unwrap().to_string()
    };
    assert_eq!(counted_u("Gene"), "1");
    assert_eq!(counted_u("GeneFull"), "2");
}

#[test]
fn gene_name_tag_applies_to_a_gtf_given_at_mapping() {
    let fx = Fixture::new();
    let (c, b) = fx.reads(&[EXONIC], &[format!("{CB}AAAAACCCCC")]);
    let wl = fx.whitelist("wl.txt", CB);
    fx.run(
        &c,
        &b,
        &[
            "--soloType",
            "CB_UMI_Simple",
            "--soloCBwhitelist",
            wl.to_str().unwrap(),
            "--sjdbGTFtagExonParentGeneName",
            "gene_symbol",
        ],
    )
    .success();
    let features = fs::read_to_string(fx.out().join("Solo.out/Gene/raw/features.tsv")).unwrap();
    assert!(features.contains("G1\tSYMBOL"), "{features}");
}

/// `soloCellFiltering` writes a plain matrix, so a gzipped feature list has to be
/// decompressed on the way through, not copied into a `.tsv`.
#[test]
fn cell_filtering_writes_plain_features_from_gzipped_ones() {
    let dir = TempDir::new().unwrap();
    let raw = dir.path().join("raw");
    fs::create_dir_all(&raw).unwrap();
    fs::write(
        raw.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n%\n2 2 2\n1 1 500\n2 2 1\n",
    )
    .unwrap();
    fs::write(raw.join("barcodes.tsv"), "AAAA\nCCCC\n").unwrap();
    let mut gz = flate2::write::GzEncoder::new(
        fs::File::create(raw.join("features.tsv.gz")).unwrap(),
        flate2::Compression::default(),
    );
    gz.write_all(b"g1\tG1\tGene Expression\ng2\tG2\tGene Expression\n")
        .unwrap();
    gz.finish().unwrap();
    let out = dir.path().join("filtered/");
    cargo_bin_cmd!("rustar-aligner")
        .args(["--runMode", "soloCellFiltering"])
        .arg(&raw)
        .arg(&out)
        .args(["--soloCellFilter", "TopCells", "1"])
        .assert()
        .success();
    let written = fs::read(out.join("features.tsv")).unwrap();
    assert!(written.starts_with(b"g1\tG1"), "{written:?}");
}
