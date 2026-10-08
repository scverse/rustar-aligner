//! Micro-benchmarks for the paths the open "measure first" questions land on.
//!
//! Several dependency decisions (#162, #202, #205, #208) all begin with a
//! measurement, and each was about to invent its own. These are the shared
//! ones, deliberately small: pure functions with no genome index to build, so
//! a run takes seconds and a regression is attributable to one function.
//!
//! ```text
//! cargo bench                       # everything
//! cargo bench -- seed_scan          # one group
//! cargo bench -- --sample-count 200 # more samples
//! ```
//!
//! Divan reports allocation counts alongside wall time, which matters here:
//! two of the open questions (`sufr`/`libsais` against `caps-sa`) are about
//! peak RSS as much as speed.

use std::collections::HashMap;

use rustar_aligner::align::simd_scan::find_stop;
use rustar_aligner::genome::Genome;
use rustar_aligner::junction::gtf::GtfRecord;
use rustar_aligner::quant::GeneAnnotation;

fn main() {
    divan::main();
}

/// Deterministic pseudo-random bases (0..=3), the same generator the tests use.
fn lcg_bases(seed: u32, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            ((state >> 16) & 3) as u8
        })
        .collect()
}

// ── Seed extension ──────────────────────────────────────────────────────────
//
// `find_stop` is the innermost loop of seed extension: it walks read against
// genome and returns the first mismatch or padding byte. It runs once per
// candidate extension, which is millions of times per million reads, and it is
// the function a portable-SIMD crate (#205) would replace.

#[divan::bench(args = [50, 100, 250])]
fn seed_scan_full_match(bencher: divan::Bencher, len: usize) {
    let read = lcg_bases(1, len);
    let genome = read.clone();
    bencher.bench(|| find_stop(divan::black_box(&read), divan::black_box(&genome)));
}

#[divan::bench(args = [50, 100, 250])]
fn seed_scan_stop_at_half(bencher: divan::Bencher, len: usize) {
    // The common case in practice: a long run of matches, then a mismatch.
    let read = lcg_bases(1, len);
    let mut genome = read.clone();
    genome[len / 2] = (genome[len / 2] + 1) % 4;
    bencher.bench(|| find_stop(divan::black_box(&read), divan::black_box(&genome)));
}

// ── Gene overlap ────────────────────────────────────────────────────────────
//
// The segment-tree overlap query runs per read (twice for `GeneFull`) and was
// the top solo hotspot before it replaced a linear scan. It is what the
// interval-crate question (#208) would replace.

fn synthetic_annotation(n_genes: usize) -> (GeneAnnotation, Genome) {
    let genome_len = (n_genes as u64 + 2) * 1_000;
    let genome = Genome {
        transform_blocks: None,
        sequence: vec![0u8; genome_len as usize].into(),
        n_genome: genome_len,
        n_genome_real: genome_len,
        n_chr_real: 1,
        chr_start: vec![0, genome_len],
        chr_length: vec![genome_len],
        chr_name: vec!["chr1".to_string()],
    };

    let exons: Vec<GtfRecord> = (0..n_genes)
        .flat_map(|g| {
            // Two exons per gene, genes 1 kb apart with a little overlap
            // between neighbours so the query cannot stop at the first hit.
            let base = g as u64 * 1_000 + 100;
            [(base, base + 400), (base + 600, base + 1_100)]
                .into_iter()
                .map(move |(s, e)| {
                    let mut attributes = HashMap::new();
                    attributes.insert("gene_id".to_string(), format!("G{g}"));
                    attributes.insert("transcript_id".to_string(), format!("G{g}_T1"));
                    GtfRecord {
                        seqname: "chr1".to_string(),
                        feature: "exon".to_string(),
                        start: s + 1,
                        end: e,
                        strand: '+',
                        attributes,
                    }
                })
        })
        .collect();

    (GeneAnnotation::from_gtf_exons(&exons, &genome), genome)
}

#[divan::bench(args = [100, 2_000, 20_000])]
fn gene_overlap_query(bencher: divan::Bencher, n_genes: usize) {
    use rustar_aligner::align::transcript::{Exon, Transcript};

    let (ann, _genome) = synthetic_annotation(n_genes);
    // A read landing in the middle of the annotation, spanning two exons.
    let mid = (n_genes as u64 / 2) * 1_000 + 150;
    let transcript = Transcript {
        chr_idx: 0,
        genome_start: mid,
        genome_end: mid + 800,
        is_reverse: false,
        exons: vec![
            Exon {
                genome_start: mid,
                genome_end: mid + 200,
                read_start: 0,
                read_end: 200,
                i_frag: 0,
            },
            Exon {
                genome_start: mid + 600,
                genome_end: mid + 800,
                read_start: 200,
                read_end: 400,
                i_frag: 0,
            },
        ],
        cigar: Vec::new(),
        score: 0,
        n_mismatch: 0,
        n_gap: 0,
        n_junction: 1,
        junction_motifs: Vec::new(),
        junction_annotated: Vec::new(),
    };

    let mut out = Vec::new();
    bencher.bench_local(|| {
        ann.overlapping_genes_into(divan::black_box(&transcript), &mut out);
        out.len()
    });
}

// ── Annotation build ────────────────────────────────────────────────────────
//
// Building the annotation is once per run, but it is O(exons log exons) and
// shows up on the solo startup path; it is also the allocation-heavy half of
// the interval question.

#[divan::bench(args = [2_000, 20_000])]
fn gene_annotation_build(bencher: divan::Bencher, n_genes: usize) {
    bencher
        .with_inputs(|| n_genes)
        .bench_values(|n| synthetic_annotation(n).0.n_genes());
}
