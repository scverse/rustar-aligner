//! Compare `solo::cr_annot` against CellRanger BAM annotation tags.
//!
//! Usage: samtools view bam | cargo run --release --example crannot_check -- <index_dir>
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use rustar_aligner::align::transcript::Transcript;
use rustar_aligner::genome::{Genome, GenomeSeq};
use rustar_aligner::quant::transcriptome::TranscriptomeIndex;
use rustar_aligner::solo::cr_annot::CrModel;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

fn main() {
    let dir = std::env::args().nth(1).expect("index dir");
    let dir = Path::new(&dir);
    let names: Vec<(String, u64)> = std::fs::read_to_string(dir.join("chrNameLength.txt"))
        .unwrap()
        .lines()
        .map(|l| {
            let mut f = l.split('\t');
            (
                f.next().unwrap().to_string(),
                f.next().unwrap().parse().unwrap(),
            )
        })
        .collect();
    let chr_start: Vec<u64> = std::fs::read_to_string(dir.join("chrStart.txt"))
        .unwrap()
        .lines()
        .map(|l| l.trim().parse().unwrap())
        .collect();
    let n = names.len();
    let genome = Genome {
        sequence: GenomeSeq::Owned(Vec::new()),
        n_genome: *chr_start.last().unwrap(),
        n_genome_real: *chr_start.last().unwrap(),
        n_chr_real: n,
        chr_name: names.iter().map(|x| x.0.clone()).collect(),
        chr_length: names.iter().map(|x| x.1).collect(),
        chr_start: chr_start.clone(),
        transform_blocks: None,
    };
    let tr = TranscriptomeIndex::from_index_dir(dir, &genome).unwrap();
    let model = CrModel::from_transcriptome(&tr, n);
    let chr_idx: HashMap<&str, usize> = genome
        .chr_name
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let (mut tot, mut skipped) = (0u64, 0u64);
    let mut bad: HashMap<&str, u64> = HashMap::new();
    let mut shown: HashMap<&str, u32> = HashMap::new();
    let (mut multi, mut multi_bad) = (0u64, 0u64);
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        let flag: u32 = f[1].parse().unwrap();
        if f[2] == "*" || flag & 4 != 0 || !chr_idx.contains_key(f[2]) {
            skipped += 1;
            continue;
        }
        let ci = chr_idx[f[2]];
        let pos: u64 = f[3].parse::<u64>().unwrap() - 1 + chr_start[ci];
        let mut cig = Vec::new();
        let mut num = 0usize;
        for c in f[5].chars() {
            if let Some(d) = c.to_digit(10) {
                num = num * 10 + d as usize;
            } else {
                let k = match c {
                    'M' => Kind::Match,
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    'N' => Kind::Skip,
                    'S' => Kind::SoftClip,
                    'H' => Kind::HardClip,
                    'P' => Kind::Pad,
                    '=' => Kind::SequenceMatch,
                    _ => Kind::SequenceMismatch,
                };
                cig.push(Op::new(k, num));
                num = 0;
            }
        }
        let mut tags: HashMap<&str, &str> = HashMap::new();
        for t in &f[11..] {
            tags.insert(&t[..2], &t[5..]);
        }
        let t = Transcript {
            chr_idx: ci,
            genome_start: pos,
            genome_end: pos,
            is_reverse: flag & 16 != 0,
            exons: Vec::new(),
            cigar: cig.clone(),
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            star_order: 0,
        };
        let r = model.annotate_read(&[t]).unwrap();
        let a = &r.alns[0];
        let nh: u32 = tags.get("NH").map_or(1, |v| v.parse().unwrap());
        let (gx, gn) = model.gx_gn(&r.genes);
        let opt = |s: String| if s.is_empty() { None } else { Some(s) };
        let mut mism: Vec<&'static str> = Vec::new();
        let chk = |name: &'static str, ours: Option<String>, mism: &mut Vec<&'static str>| {
            if ours.as_deref() != tags.get(name).copied() {
                mism.push(name);
            }
        };
        chk("RE", Some(a.region.tag_char().to_string()), &mut mism);
        chk("TX", opt(model.tx_tag(&cig, a)), &mut mism);
        chk("AN", opt(model.an_tag(&cig, a)), &mut mism);
        chk("GX", opt(gx), &mut mism);
        chk("GN", opt(gn), &mut mism);
        let xf: u32 = tags.get("xf").map_or(0, |v| v.parse().unwrap());
        if nh == 1 {
            tot += 1;
            if (xf & 17 == 17) != r.gene.is_some() {
                mism.push("xf");
            }
            let mq: u32 = f[4].parse().unwrap();
            if mq != 255 {
                mism.push("MAPQ");
            }
            if tags.contains_key("mm") {
                mism.push("mm");
            }
        } else {
            multi += 1;
            // Cannot rebuild the other loci: only check consistency.
            let mm = tags.contains_key("mm");
            let mq: u32 = f[4].parse().unwrap();
            if mm != (mq == 255) {
                mism.push("mm/MAPQ");
            }
            if mm && xf & 1 == 0 && r.gene.is_some() {
                mism.push("xf(multi)");
            }
            if !mism.is_empty() {
                multi_bad += 1;
            }
        }
        for m in &mism {
            *bad.entry(m).or_default() += 1;
            let c = shown.entry(m).or_default();
            if *c < 4 {
                *c += 1;
                eprintln!("MISMATCH {m} NH={nh}: {}", &line[..line.len().min(400)]);
            }
        }
    }
    println!(
        "NH==1 records {tot}, skipped {skipped}, multimappers {multi} (with a mismatch: {multi_bad})"
    );
    for (k, v) in &bad {
        println!("  {k}: {v}");
    }
}
