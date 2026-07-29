//! `--genomeType SuperTranscriptome`: condense a genome to the union of its
//! annotated exons (STAR `GTF_superTranscript.cpp` +
//! `SuperTranscriptome::sjCollapse`).
//!
//! Instead of indexing the whole genome, the index covers only exonic
//! sequence: overlapping exons are merged, the merged intervals are
//! concatenated, and the transcripts that overlap are grouped into
//! superTranscripts, one per condensed chromosome (`st0`, `st1`, ...). Reads
//! then align against a reference a fraction of the size, and introns cannot
//! be crossed because they are not there.
//!
//! The output is a genome and an annotation in condensed coordinates, which
//! the ordinary build path consumes unchanged, plus four files STAR writes for
//! downstream use: the superTranscript and transcript sequences, the collapsed
//! junction table, and the condensed-to-full coordinate map.
//!
//! # D20: minus-strand exons
//!
//! STAR condenses the superTranscript sequence *before* it fills the
//! reverse-complement half of its genome buffer, about twenty-five lines later
//! in `Genome_genomeGenerate.cpp`. Every minus-strand exon therefore reads the
//! allocation's spacer bytes, and minus-strand superTranscripts come out as
//! uninitialised sequence.
//!
//! This does not reproduce that. Here the reverse-complement half is already
//! filled when condensing runs, so a minus-strand exon yields its actual
//! reverse complement. The divergence is deliberate: STAR's output on that path
//! is not a rule to follow, it is a read of uninitialised memory.
//! `minus_strand_supertranscript_is_reverse_complement_not_spacer` locks the
//! correct answer, derived by hand.

use std::collections::HashMap;

use crate::error::Error;
use crate::genome::Genome;
use crate::genome::fasta::Chromosome;
use crate::junction::gtf::GtfRecord;

/// The condensed genome and everything derived from it.
pub struct SuperTranscriptome {
    /// Condensed chromosomes `st0..`, base codes, ready for `Genome::from_chromosomes`.
    pub chromosomes: Vec<Chromosome>,
    /// The annotation in condensed coordinates: `seqname` is the
    /// superTranscript, positions are 1-based within it, strand is always `+`
    /// because a superTranscript has no orientation left to speak of.
    pub exons: Vec<GtfRecord>,
    /// `superTranscriptSequences.fasta`.
    pub supertranscript_fasta: Vec<u8>,
    /// `transcriptSequences.fasta`.
    pub transcript_fasta: Vec<u8>,
    /// `superTranscriptSJcollapsed.tsv`.
    pub sj_collapsed_tsv: Vec<u8>,
    /// `fullGenome/conversionToFullGenome.tsv`.
    pub conversion_tsv: Vec<u8>,
}

/// One exon in global doubled-genome coordinates, tagged with its transcript.
#[derive(Clone, Copy)]
struct Exon {
    transcript: usize,
    start: u64,
    end: u64,
    /// Index of the GTF record this came from, so the remapped annotation can
    /// keep its attributes — gene id, gene name, gene type — rather than
    /// dropping everything the transcriptome tables need.
    source: usize,
}

/// Condense `genome` to the union of the exons in `records`.
///
/// `transcript_tag` is `--sjdbGTFtagExonParentTranscript`. Records naming a
/// chromosome the genome does not have, or carrying no transcript tag, are
/// skipped, matching what the junction extractor does with them.
pub fn condense(
    genome: &Genome,
    records: &[GtfRecord],
    transcript_tag: &str,
    chr_bin_nbits: u32,
) -> Result<SuperTranscriptome, Error> {
    let chr_index: HashMap<&str, usize> = genome
        .chr_name
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    // Transcripts in first-appearance order, so a rebuild of the same GTF
    // produces the same superTranscript numbering.
    let mut transcript_ids: Vec<String> = Vec::new();
    let mut transcript_of: HashMap<String, usize> = HashMap::new();
    let mut transcript_strand: Vec<char> = Vec::new();
    let mut exons: Vec<Exon> = Vec::new();

    for (src, rec) in records.iter().enumerate() {
        let Some(&ci) = chr_index.get(rec.seqname.as_str()) else {
            continue;
        };
        let Some(id) = rec.attributes.get(transcript_tag) else {
            continue;
        };
        let t = *transcript_of.entry(id.clone()).or_insert_with(|| {
            transcript_ids.push(id.clone());
            transcript_strand.push(rec.strand);
            transcript_ids.len() - 1
        });
        // GTF is 1-based inclusive; the genome is 0-based.
        exons.push(Exon {
            transcript: t,
            start: genome.chr_start[ci] + rec.start - 1,
            end: genome.chr_start[ci] + rec.end - 1,
            source: src,
        });
    }

    if exons.is_empty() {
        return Err(Error::Gtf(
            "--genomeType SuperTranscriptome needs at least one exon on a known chromosome"
                .to_string(),
        ));
    }

    // (1) Reflect minus-strand exons into the reverse-complement half. See D20
    // in the module docs for why this reads real sequence here and does not in
    // STAR.
    let n = genome.n_genome;
    for e in &mut exons {
        if transcript_strand[e.transcript] == '-' {
            let (s, t) = (e.start, e.end);
            e.start = 2 * n - 1 - t;
            e.end = 2 * n - 1 - s;
        }
    }

    // (2) Sort by start, (3) merge overlapping or touching exons and rewrite
    // every exon into condensed coordinates as we go.
    exons.sort_by_key(|e| e.start);
    let mut merged: Vec<[u64; 2]> = Vec::new();
    let mut gap = exons[0].start;
    let mut curr = [exons[0].start, exons[0].end];
    exons[0].start = 0;
    exons[0].end -= gap;
    for e in exons.iter_mut().skip(1) {
        if e.start <= curr[1] + 1 {
            curr[1] = curr[1].max(e.end);
        } else {
            gap += e.start - curr[1] - 1;
            merged.push(curr);
            curr = [e.start, e.end];
        }
        e.start -= gap;
        e.end -= gap;
    }
    merged.push(curr);

    // The condensed sequence: the merged intervals' bases, back to back.
    let mut concat: Vec<u8> = Vec::new();
    for m in &merged {
        for j in m[0]..=m[1] {
            concat.push(genome.get_base(j).unwrap_or(4));
        }
    }

    // (4) Each transcript's span in condensed coordinates.
    let n_transcripts = transcript_ids.len();
    let mut tr_span = vec![[u64::MAX, 0u64]; n_transcripts];
    for e in &exons {
        tr_span[e.transcript][0] = tr_span[e.transcript][0].min(e.start);
        tr_span[e.transcript][1] = tr_span[e.transcript][1].max(e.end);
    }

    // (5) Group overlapping transcript spans into superTranscripts, aligned to
    // merged-interval boundaries so a superTranscript never starts mid-interval.
    let mut merged_super = vec![0u64; merged.len()];
    let mut spans = tr_span.clone();
    spans.sort_by_key(|s| s[0]);
    let mut super_span: Vec<[u64; 2]> = Vec::new();
    let mut mi = 0usize;
    let mut mi_len = 0u64;
    let mut curr = spans[0];
    for s in &spans {
        while mi_len < curr[1] && mi < merged.len() {
            mi_len += merged[mi][1] - merged[mi][0] + 1;
            merged_super[mi] = super_span.len() as u64;
            mi += 1;
        }
        curr[1] = curr[1].max(mi_len.saturating_sub(1));
        if s[0] <= curr[1] {
            curr[1] = curr[1].max(s[1]);
        } else {
            super_span.push(curr);
            curr = *s;
        }
    }
    super_span.push(curr);
    // Intervals past the last transcript span still belong somewhere.
    for m in merged_super.iter_mut().skip(mi) {
        *m = super_span.len() as u64 - 1;
    }

    // (6) The superTranscript sequences.
    let super_seq: Vec<Vec<u8>> = super_span
        .iter()
        .map(|s| concat[s[0] as usize..=s[1] as usize].to_vec())
        .collect();

    // (7) Which superTranscript each transcript landed in.
    let mut tr_super = vec![0u64; n_transcripts];
    {
        let mut ist = 0usize;
        for e in &exons {
            while ist + 1 < super_span.len() && e.start > super_span[ist][1] {
                ist += 1;
            }
            tr_super[e.transcript] = ist as u64;
        }
    }

    // (8) Transcript sequences: exons in transcript order, concatenated.
    exons.sort_by_key(|e| (e.transcript, e.start));
    let mut transcript_seq: Vec<Vec<u8>> = vec![Vec::new(); n_transcripts];
    for e in &exons {
        transcript_seq[e.transcript].extend_from_slice(&concat[e.start as usize..=e.end as usize]);
    }

    let mut transcript_fasta = Vec::new();
    for (i, seq) in transcript_seq.iter().enumerate() {
        transcript_fasta.push(b'>');
        transcript_fasta.extend_from_slice(transcript_ids[i].as_bytes());
        transcript_fasta.push(b'\n');
        transcript_fasta.extend(seq.iter().map(|&c| code_to_char(c)));
        transcript_fasta.push(b'\n');
    }

    let mut supertranscript_fasta = Vec::new();
    for (i, seq) in super_seq.iter().enumerate() {
        supertranscript_fasta.extend_from_slice(format!(">st{i}\n").as_bytes());
        supertranscript_fasta.extend(seq.iter().map(|&c| code_to_char(c)));
        supertranscript_fasta.push(b'\n');
    }

    // (9) Junctions, in superTranscript-relative coordinates.
    let mut sj: Vec<(u64, u64, u64)> = Vec::new(); // (superTranscript, start, end)
    for i in 1..exons.len() {
        if exons[i].transcript == exons[i - 1].transcript && exons[i].start > exons[i - 1].end + 1 {
            let st = tr_super[exons[i].transcript];
            let base = super_span[st as usize][0];
            sj.push((st, exons[i - 1].end - base, exons[i].start - base));
        }
    }
    let sj_collapsed_tsv = collapse_sj(&mut sj);

    // (10) The condensed chromosomes, and the annotation rewritten onto them.
    let chromosomes: Vec<Chromosome> = super_seq
        .iter()
        .enumerate()
        .map(|(i, s)| Chromosome {
            name: format!("st{i}"),
            sequence: s.clone(),
        })
        .collect();

    exons.sort_by_key(|e| e.start);
    let mut out_records: Vec<GtfRecord> = Vec::with_capacity(exons.len());
    {
        let mut ist = 0usize;
        for e in &exons {
            while ist + 1 < super_span.len() && e.start > super_span[ist][1] {
                ist += 1;
            }
            let base = super_span[ist][0];
            let source = &records[e.source];
            out_records.push(GtfRecord {
                seqname: format!("st{ist}"),
                feature: source.feature.clone(),
                attributes: source.attributes.clone(),
                start: e.start - base + 1, // back to 1-based, superTranscript-local
                end: e.end - base + 1,
                strand: '+',
            });
        }
    }

    // (11) conversionToFullGenome.tsv: header, then one line per merged block.
    let chr_start = super::compute_chr_starts(
        &chromosomes
            .iter()
            .map(|c| c.sequence.len() as u64)
            .collect::<Vec<_>>(),
        chr_bin_nbits,
    );
    let mut conversion_tsv = Vec::new();
    conversion_tsv.extend_from_slice(format!("{}\t{}\n", merged.len(), n).as_bytes());
    let mut cond = 0u64;
    for (ib, b) in merged.iter().enumerate() {
        let len = b[1] - b[0] + 1;
        let st = merged_super[ib] as usize;
        let start_in_super = chr_start[st] + cond - super_span[st][0];
        conversion_tsv.extend_from_slice(format!("{start_in_super}\t{len}\t{}\n", b[0]).as_bytes());
        cond += len;
    }

    Ok(SuperTranscriptome {
        chromosomes,
        exons: out_records,
        supertranscript_fasta,
        transcript_fasta,
        sj_collapsed_tsv,
        conversion_tsv,
    })
}

/// Genome code to FASTA letter. Codes above 4 cannot occur here: every base
/// comes from a merged interval inside the real genome.
fn code_to_char(code: u8) -> u8 {
    match code {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

/// STAR's `SuperTranscriptome::sjCollapse`: sort by (superTranscript, start,
/// end), drop exact duplicates, emit one line each.
fn collapse_sj(sj: &mut [(u64, u64, u64)]) -> Vec<u8> {
    sj.sort_unstable();
    let mut out = Vec::new();
    let mut last: Option<(u64, u64, u64)> = None;
    for &j in sj.iter() {
        if last != Some(j) {
            out.extend_from_slice(format!("{}\t{}\t{}\n", j.0, j.1, j.2).as_bytes());
            last = Some(j);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome_of(seq: &str) -> Genome {
        let codes: Vec<u8> = seq
            .bytes()
            .map(|b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();
        Genome::from_chromosomes(
            &[Chromosome {
                name: "chr1".to_string(),
                sequence: codes,
            }],
            18,
            None,
        )
        .unwrap()
    }

    fn exon(start: u64, end: u64, strand: char, transcript: &str) -> GtfRecord {
        let mut attributes = HashMap::new();
        attributes.insert("transcript_id".to_string(), transcript.to_string());
        GtfRecord {
            seqname: "chr1".to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes,
        }
    }

    fn fasta_text(bytes: &[u8]) -> String {
        String::from_utf8(bytes.to_vec()).unwrap()
    }

    /// Two exons with an intron between them collapse to one superTranscript
    /// holding only exonic sequence, and the intron is gone rather than spliced.
    #[test]
    fn exons_are_concatenated_and_the_intron_disappears() {
        //            1234567890123456789012345678901234567890
        let genome = genome_of("AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTT");
        // Exons 1-8 and 25-32 (1-based), one transcript.
        let records = vec![exon(1, 8, '+', "t1"), exon(25, 32, '+', "t1")];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();

        assert_eq!(st.chromosomes.len(), 1);
        assert_eq!(
            fasta_text(&st.supertranscript_fasta),
            ">st0\nAAAACCCCAAAACCCC\n"
        );
        assert_eq!(fasta_text(&st.transcript_fasta), ">t1\nAAAACCCCAAAACCCC\n");
        // No junction: with nothing else annotated between them, the two
        // exons end up contiguous in the condensed genome, so there is no gap
        // left for a read to cross.
        assert_eq!(fasta_text(&st.sj_collapsed_tsv), "");
    }

    /// A junction survives condensing only when another transcript's exon sits
    /// between the two, keeping them apart in the condensed genome.
    #[test]
    fn a_junction_is_recorded_when_another_exon_separates_the_two() {
        let genome = genome_of("AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTT");
        let records = vec![
            exon(1, 8, '+', "t1"),   // t1 exon 1
            exon(13, 20, '+', "t2"), // t2, between them
            exon(25, 32, '+', "t1"), // t1 exon 2
        ];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();
        // Condensed: t1[0..7], t2[8..15], t1[16..23]. t1 jumps 7 -> 16.
        assert_eq!(fasta_text(&st.sj_collapsed_tsv), "0\t7\t16\n");
    }

    /// The condensed annotation addresses the superTranscript, one-based, so
    /// the ordinary GTF consumers need no special case.
    #[test]
    fn the_annotation_comes_back_in_condensed_coordinates() {
        let genome = genome_of("AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTT");
        let records = vec![exon(1, 8, '+', "t1"), exon(25, 32, '+', "t1")];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();

        assert_eq!(st.exons.len(), 2);
        assert_eq!(st.exons[0].seqname, "st0");
        assert_eq!((st.exons[0].start, st.exons[0].end), (1, 8));
        assert_eq!((st.exons[1].start, st.exons[1].end), (9, 16));
        assert!(st.exons.iter().all(|e| e.strand == '+'));
    }

    /// Overlapping exons from different transcripts are merged once, not
    /// duplicated, which is the whole point of condensing.
    #[test]
    fn overlapping_exons_are_merged_once() {
        let genome = genome_of("AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTT");
        let records = vec![exon(1, 12, '+', "t1"), exon(5, 16, '+', "t2")];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();
        assert_eq!(
            fasta_text(&st.supertranscript_fasta),
            ">st0\nAAAACCCCGGGGTTTT\n",
            "the union of the two exons, each base once"
        );
    }

    /// D20. STAR reads its unfilled reverse strand here and emits spacer bytes;
    /// this emits the reverse complement, derived by hand below.
    #[test]
    fn minus_strand_supertranscript_is_reverse_complement_not_spacer() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let genome = genome_of(seq);
        // One minus-strand exon at 1-based 11..20 — forward offsets 10..=19.
        let records = vec![exon(11, 20, '-', "tM")];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();

        let complement = |b: u8| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        };
        let expected: String = seq.as_bytes()[10..20]
            .iter()
            .rev()
            .map(|&b| complement(b) as char)
            .collect();

        assert_eq!(
            fasta_text(&st.supertranscript_fasta),
            format!(">st0\n{expected}\n"),
            "a minus-strand exon must yield real sequence, not the spacer STAR reads"
        );
        assert!(
            !st.supertranscript_fasta.contains(&0x00),
            "no uninitialised bytes may reach the output"
        );
    }

    /// Two transcripts that share no sequence become two superTranscripts, so
    /// the aligner cannot stitch across a boundary that does not exist.
    #[test]
    fn disjoint_transcripts_become_separate_supertranscripts() {
        let genome = genome_of("AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTT");
        let records = vec![exon(1, 8, '+', "t1"), exon(33, 40, '+', "t2")];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();
        assert_eq!(st.chromosomes.len(), 2);
        assert_eq!(st.chromosomes[0].name, "st0");
        assert_eq!(st.chromosomes[1].name, "st1");
    }

    /// The conversion map is what makes condensed coordinates translatable
    /// back, so its blocks must cover every merged interval.
    #[test]
    fn the_conversion_map_has_one_block_per_merged_interval() {
        let genome = genome_of("AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTT");
        let records = vec![exon(1, 8, '+', "t1"), exon(25, 32, '+', "t1")];
        let st = condense(&genome, &records, "transcript_id", 18).unwrap();
        let text = fasta_text(&st.conversion_tsv);
        let mut lines = text.lines();
        assert_eq!(lines.next().unwrap().split('\t').next().unwrap(), "2");
        let blocks: Vec<&str> = lines.collect();
        assert_eq!(blocks.len(), 2);
        // condensed start, length, full-genome start
        assert_eq!(blocks[0], "0\t8\t0");
        assert_eq!(blocks[1], "8\t8\t24");
    }

    #[test]
    fn a_gtf_with_no_usable_exon_is_an_error() {
        let genome = genome_of("AAAACCCC");
        let mut rec = exon(1, 4, '+', "t1");
        rec.seqname = "chrUnknown".to_string();
        assert!(condense(&genome, &[rec], "transcript_id", 18).is_err());
    }
}
