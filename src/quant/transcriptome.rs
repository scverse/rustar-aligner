//! Transcript-level quantification (`--quantMode TranscriptomeSAM`).
//!
//! Builds a per-transcript exon map from GTF records and projects genome-space
//! alignments onto the transcriptome coordinate system so downstream tools
//! (Salmon, RSEM) can read `Aligned.toTranscriptome.out.bam`.
//!
//! This module mirrors STAR's `Transcriptome` / `alignToTranscript` /
//! `quantAlign` logic (see `source/Transcriptome.cpp` and
//! `source/Transcriptome_quantAlign.cpp` in the upstream repo), with the
//! only substantive divergence that rustar-aligner builds the transcript tables on
//! the fly from the input GTF instead of loading persisted
//! `transcriptInfo.tab`/`exonInfo.tab` files.
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;

use crate::align::transcript::{CigarOp, Exon, Transcript};
use crate::error::Error;
use crate::genome::Genome;
use crate::junction::gtf::GtfRecord;
use crate::params::Parameters;

/// GTF attribute names the transcriptome index reads.
const GTF_ATTR_TRANSCRIPT_ID: &str = "transcript_id";
const GTF_ATTR_GENE_ID: &str = "gene_id";
const GTF_ATTR_GENE_NAME: &str = "gene_name";
const GTF_ATTR_GENE_BIOTYPE: &str = "gene_biotype";

/// STAR's fallback for GTF records missing `gene_biotype`
/// (`source/GTF.cpp`).
const MISSING_GENE_TYPE: &str = "MissingGeneType";

// ---------------------------------------------------------------------------
// Filter-mode enum + softclip extension (subtask 3)
// ---------------------------------------------------------------------------

/// STAR's `--quantTranscriptomeSAMoutput` value.
///
/// Controls how alignments are filtered / adjusted before being projected
/// onto the transcriptome:
///   * `BanSingleEndBanIndelsExtendSoftclip` (STAR default, RSEM-compatible) —
///     reject alignments containing indels; reject PE alignments where both
///     sides of the read come from a single mate; extend soft-clipped bases
///     back into matched bases (re-counting mismatches).
///   * `BanSingleEnd` — only reject single-mate-only PE alignments; keep
///     indels and soft-clips as-is.
///   * `BanSingleEndExtendSoftclip` — reject single-mate-only PE alignments;
///     keep indels; extend soft-clips.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QuantTranscriptomeSAMoutput {
    /// Default: ban indels, ban single-mate PE hits, extend soft-clips.
    #[default]
    BanSingleEndBanIndelsExtendSoftclip,
    /// Ban single-mate PE hits only; keep indels and soft-clips.
    BanSingleEnd,
    /// Ban single-mate PE hits; keep indels; extend soft-clips.
    BanSingleEndExtendSoftclip,
}

impl std::str::FromStr for QuantTranscriptomeSAMoutput {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "BanSingleEnd_BanIndels_ExtendSoftclip" => {
                Ok(Self::BanSingleEndBanIndelsExtendSoftclip)
            }
            "BanSingleEnd" => Ok(Self::BanSingleEnd),
            "BanSingleEnd_ExtendSoftclip" => Ok(Self::BanSingleEndExtendSoftclip),
            _ => Err(format!(
                "unknown --quantTranscriptomeSAMoutput '{s}'; expected \
                 'BanSingleEnd_BanIndels_ExtendSoftclip', 'BanSingleEnd', or \
                 'BanSingleEnd_ExtendSoftclip'"
            )),
        }
    }
}

impl QuantTranscriptomeSAMoutput {
    /// Whether indels are permitted in the output.
    pub fn allow_indels(self) -> bool {
        !matches!(self, Self::BanSingleEndBanIndelsExtendSoftclip)
    }
    /// Whether soft-clipped bases should be kept (true) or extended back into
    /// matched bases (false).
    pub fn allow_softclip(self) -> bool {
        matches!(self, Self::BanSingleEnd)
    }
}

/// Per-transcript exon in absolute genome coordinates (0-based half-open),
/// paired with the cumulative transcript-space length of all preceding exons
/// (STAR's `exLenCum`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrExon {
    /// Absolute genome start (0-based, inclusive).
    pub genome_start: u64,
    /// Absolute genome end (0-based, exclusive).
    pub genome_end: u64,
    /// Cumulative transcript-space length of all prior exons (0 for the first).
    pub ex_len_cum: u32,
}

/// Transcriptome metadata built from GTF exon records.
///
/// Each transcript is a contiguous set of exons on a single chromosome and
/// strand, ordered by ascending genome start. `tr_order` / `tr_starts_sorted`
/// / `tr_end_max_sorted` support STAR's `binarySearch1a` + running-max
/// early-exit in `quantAlign`.
#[derive(Debug, Clone)]
pub struct TranscriptomeIndex {
    /// Transcript IDs in insertion order (index = transcript_idx).
    pub tr_ids: Vec<String>,
    /// Genome chromosome index per transcript.
    pub tr_chr_idx: Vec<usize>,
    /// Strand per transcript: 1 = forward/`+`, 2 = reverse/`-`.
    pub tr_strand: Vec<u8>,
    /// Index into `gene_ids` per transcript — matches STAR's `trGene` column.
    /// Look up the string form via `gene_ids[tr_gene_idx[tr] as usize]`.
    pub tr_gene_idx: Vec<u32>,
    /// Unique gene IDs in first-seen order (STAR's `geneInfo.tab` column 1).
    pub gene_ids: Vec<String>,
    /// gene_name per gene. Falls back to the `gene_id` string when the GTF
    /// has no `gene_name` attribute — matches STAR's `GTF.cpp` behavior.
    pub gene_names: Vec<String>,
    /// gene_biotype per gene. Falls back to the literal string
    /// `"MissingGeneType"` when the GTF has no `gene_biotype` attribute —
    /// matches STAR's `GTF.cpp` behavior for `geneAttr[ig][1]`.
    pub gene_biotypes: Vec<String>,
    /// Absolute genome start of the first exon (0-based).
    pub tr_start: Vec<u64>,
    /// Absolute genome end of the last exon (0-based half-open).
    pub tr_end: Vec<u64>,
    /// Exons per transcript (sorted by genome_start, with `ex_len_cum`).
    pub tr_exons: Vec<Vec<TrExon>>,
    /// Total transcript-space length (sum of exon spans).
    pub tr_length: Vec<u32>,
    /// First-exon offset into a flat global exon array, per transcript,
    /// computed in STAR's SORTED-BY-(trStart, trEnd) order — matches the
    /// `trExI` column in `transcriptInfo.tab`. Indexed by insertion position.
    pub tr_exi: Vec<u32>,
    /// Permutation of `0..n_transcripts` sorted by `(tr_start, tr_end)`.
    pub tr_order: Vec<usize>,
    /// `tr_start[tr_order[i]]` — for `binary_search`.
    pub tr_starts_sorted: Vec<u64>,
    /// Running max of `tr_end` along `tr_order` (STAR's `trEmax`).
    pub tr_end_max_sorted: Vec<u64>,
}

impl TranscriptomeIndex {
    /// Build from already-parsed GTF exon records with configurable GTF attribute names.
    ///
    /// `transcript_tag`: STAR `sjdbGTFtagExonParentTranscript` (default `"transcript_id"`).
    /// `gene_tag`: STAR `sjdbGTFtagExonParentGene` (default `"gene_id"`).
    pub fn from_gtf_exons_configured(
        exons: &[GtfRecord],
        genome: &Genome,
        transcript_tag: &str,
        gene_tag: &str,
    ) -> Result<Self, Error> {
        // Group exons by transcript_tag, preserving first-seen insertion order.
        let mut order: Vec<String> = Vec::new();
        let mut groups: HashMap<String, Vec<&GtfRecord>> = HashMap::new();
        for rec in exons {
            let tid = match rec.attributes.get(transcript_tag) {
                Some(id) => id.clone(),
                None => continue,
            };
            if !groups.contains_key(&tid) {
                order.push(tid.clone());
            }
            groups.entry(tid).or_default().push(rec);
        }

        let mut tr_ids: Vec<String> = Vec::new();
        let mut tr_chr_idx: Vec<usize> = Vec::new();
        let mut tr_strand: Vec<u8> = Vec::new();
        let mut tr_gene_idx: Vec<u32> = Vec::new();
        let mut tr_start: Vec<u64> = Vec::new();
        let mut tr_end: Vec<u64> = Vec::new();
        let mut tr_exons_vec: Vec<Vec<TrExon>> = Vec::new();
        let mut tr_length: Vec<u32> = Vec::new();

        // Gene-interning pass: each unique gene_id gets a position-based
        // integer index, matching STAR's geneInfo.tab / trGene column.
        let mut gene_ids: Vec<String> = Vec::new();
        let mut gene_names: Vec<String> = Vec::new();
        let mut gene_biotypes: Vec<String> = Vec::new();
        let mut gene_id_to_idx: HashMap<String, u32> = HashMap::new();

        for tid in &order {
            let mut recs = groups.remove(tid).unwrap();
            if recs.is_empty() {
                continue;
            }

            // Validate consistent chromosome + strand.
            let first = recs[0];
            let chr_name = &first.seqname;
            let strand_char = first.strand;
            let mut inconsistent = false;
            for r in &recs {
                if &r.seqname != chr_name || r.strand != strand_char {
                    inconsistent = true;
                    break;
                }
            }
            if inconsistent {
                log::warn!(
                    "quantMode TranscriptomeSAM: transcript {} has inconsistent chromosome/strand across exons — skipping",
                    tid
                );
                continue;
            }

            // Map chromosome name to genome index.
            let chr_idx = match genome.chr_name.iter().position(|n| n == chr_name) {
                Some(i) if i < genome.n_chr_real => i,
                _ => {
                    log::warn!(
                        "quantMode TranscriptomeSAM: transcript {} on unknown chromosome {} — skipping",
                        tid,
                        chr_name
                    );
                    continue;
                }
            };
            let chr_offset = genome.chr_start[chr_idx];

            // Sort exons by GTF start ASC.
            recs.sort_by_key(|r| r.start);

            // Build TrExon list in absolute genome coords (0-based half-open).
            let mut tr_exons: Vec<TrExon> = Vec::with_capacity(recs.len());
            let mut ex_len_cum: u32 = 0;
            for r in &recs {
                // GTF is 1-based inclusive; rustar-aligner uses 0-based half-open.
                let abs_start = chr_offset + r.start.saturating_sub(1);
                let abs_end = chr_offset + r.end;
                if abs_end <= abs_start {
                    log::warn!(
                        "quantMode TranscriptomeSAM: transcript {} has invalid exon {}-{} — skipping",
                        tid,
                        r.start,
                        r.end
                    );
                    tr_exons.clear();
                    break;
                }
                let len = (abs_end - abs_start) as u32;
                tr_exons.push(TrExon {
                    genome_start: abs_start,
                    genome_end: abs_end,
                    ex_len_cum,
                });
                ex_len_cum = ex_len_cum.saturating_add(len);
            }
            if tr_exons.is_empty() {
                continue;
            }

            let first_ex = tr_exons.first().unwrap();
            let last_ex = tr_exons.last().unwrap();
            let start_abs = first_ex.genome_start;
            let end_abs = last_ex.genome_end;
            let total_len = ex_len_cum;

            let strand_u8 = match strand_char {
                '+' => 1u8,
                '-' => 2u8,
                _ => 1u8, // unknown → treat as forward (STAR default)
            };

            let gene_id = first.attributes.get(gene_tag).cloned().unwrap_or_default();
            // STAR-faithful fallbacks: when the GTF record omits gene_name,
            // STAR's GTF.cpp writes geneAttr[ig][0] = gene_id (not empty).
            // When gene_biotype is omitted, it writes `MISSING_GENE_TYPE`.
            let gene_name = first
                .attributes
                .get(GTF_ATTR_GENE_NAME)
                .cloned()
                .unwrap_or_else(|| gene_id.clone());
            let gene_biotype = first
                .attributes
                .get(GTF_ATTR_GENE_BIOTYPE)
                .cloned()
                .unwrap_or_else(|| MISSING_GENE_TYPE.to_string());

            // Intern the gene — first transcript for each gene_id wins the
            // name/biotype slot. Subsequent transcripts with a richer name or
            // biotype do NOT overwrite (STAR's Transcriptome writer is
            // first-seen-wins too).
            let gene_idx = match gene_id_to_idx.get(&gene_id) {
                Some(&i) => i,
                None => {
                    let i = gene_ids.len() as u32;
                    gene_id_to_idx.insert(gene_id.clone(), i);
                    gene_ids.push(gene_id.clone());
                    gene_names.push(gene_name);
                    gene_biotypes.push(gene_biotype);
                    i
                }
            };

            tr_ids.push(tid.clone());
            tr_chr_idx.push(chr_idx);
            tr_strand.push(strand_u8);
            tr_gene_idx.push(gene_idx);
            tr_start.push(start_abs);
            tr_end.push(end_abs);
            tr_exons_vec.push(tr_exons);
            tr_length.push(total_len);
        }

        let n_tr = tr_ids.len();

        // Sorted view by (tr_start, tr_end) for binary-search +
        // running-max early-exit.
        let mut tr_order: Vec<usize> = (0..n_tr).collect();
        tr_order.sort_by(|&a, &b| {
            tr_start[a]
                .cmp(&tr_start[b])
                .then_with(|| tr_end[a].cmp(&tr_end[b]))
        });

        // Cumulative exon offset in SORTED order — matches STAR's `trExI`,
        // where exons in `exonInfo.tab` are grouped by transcript in sorted
        // transcript order. `tr_exi[i]` (insertion position i) = sum of exon
        // counts of transcripts that come BEFORE `i` in sorted order.
        let mut tr_exi: Vec<u32> = vec![0; n_tr];
        let mut cum: u32 = 0;
        for &sorted_idx in &tr_order {
            tr_exi[sorted_idx] = cum;
            cum = cum.saturating_add(tr_exons_vec[sorted_idx].len() as u32);
        }

        let tr_starts_sorted: Vec<u64> = tr_order.iter().map(|&i| tr_start[i]).collect();
        let mut tr_end_max_sorted: Vec<u64> = Vec::with_capacity(n_tr);
        let mut running_max: u64 = 0;
        for &i in &tr_order {
            running_max = running_max.max(tr_end[i]);
            tr_end_max_sorted.push(running_max);
        }

        Ok(TranscriptomeIndex {
            tr_ids,
            tr_chr_idx,
            tr_strand,
            tr_gene_idx,
            gene_ids,
            gene_names,
            gene_biotypes,
            tr_start,
            tr_end,
            tr_exons: tr_exons_vec,
            tr_length,
            tr_exi,
            tr_order,
            tr_starts_sorted,
            tr_end_max_sorted,
        })
    }

    /// Build from already-parsed GTF exon records using default STAR attribute names.
    pub fn from_gtf_exons(exons: &[GtfRecord], genome: &Genome) -> Result<Self, Error> {
        Self::from_gtf_exons_configured(exons, genome, GTF_ATTR_TRANSCRIPT_ID, GTF_ATTR_GENE_ID)
    }

    /// Number of transcripts indexed.
    pub fn n_transcripts(&self) -> usize {
        self.tr_ids.len()
    }

    /// Load from STAR-compatible index files in `dir`.
    ///
    /// Requires `transcriptInfo.tab`, `exonInfo.tab`, `geneInfo.tab` to be
    /// present — all written by either STAR's `genomeGenerate` or rustar-aligner's
    /// `write_*` methods. Returns a fully-populated index equivalent to
    /// what `from_gtf_exons` would produce, except transcripts come back
    /// in STAR's SORTED `(tr_start, tr_end)` order (not GTF insertion
    /// order), because that's how `transcriptInfo.tab` is written.
    pub fn from_index_dir(dir: &Path, genome: &Genome) -> Result<Self, Error> {
        let (tr_ids, tr_start, tr_end, tr_strand, tr_exn, tr_exi, tr_gene_idx) =
            read_transcript_info(&dir.join("transcriptInfo.tab"))?;
        let (ex_start_rel, ex_end_rel_incl, ex_len_cum_flat) =
            read_exon_info(&dir.join("exonInfo.tab"))?;
        let (gene_ids, gene_names, gene_biotypes) = read_gene_info(&dir.join("geneInfo.tab"))?;

        let n_tr = tr_ids.len();
        let n_exons_total = ex_start_rel.len();

        // Sanity: exon count in header must equal sum of trExN.
        let sum_exn: u64 = tr_exn.iter().map(|&n| n as u64).sum();
        if sum_exn != n_exons_total as u64 {
            return Err(Error::Index(format!(
                "transcriptome index inconsistent: sum(trExN)={} but exonInfo has {} rows",
                sum_exn, n_exons_total
            )));
        }

        // Reconstruct per-transcript TrExon lists + tr_length from flat exon
        // arrays. `ex_start_rel[k]` is transcript-relative; add tr_start[i] to
        // get absolute 0-based coords. `ex_end_rel_incl` is inclusive; add 1
        // for rustar-aligner's exclusive convention.
        let mut tr_exons: Vec<Vec<TrExon>> = Vec::with_capacity(n_tr);
        let mut tr_length: Vec<u32> = Vec::with_capacity(n_tr);
        for i in 0..n_tr {
            let start = tr_exi[i] as usize;
            let end = start + tr_exn[i] as usize;
            let tr_s = tr_start[i];
            let mut exons: Vec<TrExon> = Vec::with_capacity(end - start);
            let mut total: u32 = 0;
            for k in start..end {
                let gs = tr_s + ex_start_rel[k];
                let ge = tr_s + ex_end_rel_incl[k] + 1;
                exons.push(TrExon {
                    genome_start: gs,
                    genome_end: ge,
                    ex_len_cum: ex_len_cum_flat[k],
                });
                total = total.saturating_add((ge - gs) as u32);
            }
            tr_exons.push(exons);
            tr_length.push(total);
        }

        // Derive tr_chr_idx from each transcript's first exon genome_start
        // using the existing `Genome::position_to_chr` binary-search helper.
        let mut tr_chr_idx: Vec<usize> = Vec::with_capacity(n_tr);
        for exs in &tr_exons {
            let pos = exs
                .first()
                .map(|e| e.genome_start)
                .ok_or_else(|| Error::Index("transcript with zero exons".into()))?;
            let (chr_idx, _) = genome.position_to_chr(pos).ok_or_else(|| {
                Error::Index(format!(
                    "transcript first-exon position {pos} does not fall inside any chromosome"
                ))
            })?;
            tr_chr_idx.push(chr_idx);
        }

        // Sorted view: transcripts are already in sorted order on disk, so
        // tr_order = identity.
        let tr_order: Vec<usize> = (0..n_tr).collect();
        let tr_starts_sorted: Vec<u64> = tr_start.clone();
        let mut tr_end_max_sorted: Vec<u64> = Vec::with_capacity(n_tr);
        let mut m: u64 = 0;
        for &i in &tr_order {
            m = m.max(tr_end[i]);
            tr_end_max_sorted.push(m);
        }

        Ok(TranscriptomeIndex {
            tr_ids,
            tr_chr_idx,
            tr_strand,
            tr_gene_idx,
            gene_ids,
            gene_names,
            gene_biotypes,
            tr_start,
            tr_end,
            tr_exons,
            tr_length,
            tr_exi,
            tr_order,
            tr_starts_sorted,
            tr_end_max_sorted,
        })
    }

    /// Write `transcriptInfo.tab` into `dir`, byte-for-byte matching STAR's
    /// `GTF_transcriptGeneSJ.cpp:86-112` format:
    ///
    /// - Header line: integer transcript count.
    /// - Per-transcript line (in sorted order by `(tr_start, tr_end)`):
    ///   `trID\ttrStart\ttrEnd\ttrEmax\ttrStrand\ttrExN\ttrExI\ttrGene\n`
    ///
    /// `trStart` / `trEnd` are 0-based absolute genome coordinates
    /// (STAR's `exonLoci` convention: chrStart + 1-based GTF pos − 1),
    /// with `trEnd` INCLUSIVE (last base of transcript). `trEmax` is the
    /// running max of prior transcripts' `trEnd` in sorted order — with
    /// the initial value set to the first sorted transcript's `trEnd`
    /// (so sorted position 0 has `trEmax == trEnd`, matching STAR's
    /// `trend=extrLoci[GTF_extrTrEnd(0)]` initializer). `trStrand` uses
    /// STAR's GTF.cpp encoding: `+`→1, `-`→2, other→0.
    pub fn write_transcript_info(&self, dir: &Path) -> Result<(), Error> {
        let path = dir.join("transcriptInfo.tab");
        let file = std::fs::File::create(&path).map_err(|e| Error::io(e, &path))?;
        let mut out = std::io::BufWriter::new(file);

        let n_tr = self.tr_ids.len();
        writeln!(out, "{n_tr}").map_err(|e| Error::io(e, &path))?;

        let first_end_inclusive = self.tr_order.first().map(|&i| self.tr_end_inclusive(i));
        let mut tr_emax = first_end_inclusive.unwrap_or(0);

        for &i in &self.tr_order {
            let tr_end_incl = self.tr_end_inclusive(i);
            let tr_exn = self.tr_exons[i].len();
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.tr_ids[i],
                self.tr_start[i],
                tr_end_incl,
                tr_emax,
                self.tr_strand[i],
                tr_exn,
                self.tr_exi[i],
                self.tr_gene_idx[i],
            )
            .map_err(|e| Error::io(e, &path))?;

            // Running-max update matches STAR's `trend=max(trend, ...)`
            // that runs AFTER the write.
            tr_emax = tr_emax.max(tr_end_incl);
        }

        Ok(())
    }

    /// 0-based inclusive end of transcript `i` — STAR's `exE` convention.
    fn tr_end_inclusive(&self, i: usize) -> u64 {
        self.tr_end[i].saturating_sub(1)
    }

    /// Write `sjdbList.fromGTF.out.tab` into `dir`, byte-for-byte matching
    /// STAR's `GTF_transcriptGeneSJ.cpp:115-171` format:
    ///
    /// - No header line.
    /// - Per-junction line: `chr\tstart(1-based)\tend(1-based)\tstrand\tgenes\n`
    ///
    /// Junctions come from consecutive exon pairs within each transcript
    /// where an intron exists. `start` and `end` are 1-based chromosome-local
    /// positions of the first and last base of the intron. `strand` is `.`,
    /// `+`, or `-` (mapping STAR's `transcriptStrand` 0/1/2). Duplicate
    /// junctions (same chr/start/end/strand) across multiple transcripts
    /// share a single record with comma-separated 1-based gene indices.
    pub fn write_sjdb_list_from_gtf(&self, dir: &Path, genome: &Genome) -> Result<(), Error> {
        let path = dir.join("sjdbList.fromGTF.out.tab");
        let file = std::fs::File::create(&path).map_err(|e| Error::io(e, &path))?;
        let mut out = std::io::BufWriter::new(file);

        // Collect (sjStart_abs_0based_excl, sjEnd_abs_0based_incl, chr_idx,
        // strand, gene_idx_1based) for every intron within every transcript.
        let mut junctions: Vec<(u64, u64, usize, u8, u32)> = Vec::new();
        for (tr_idx, exs) in self.tr_exons.iter().enumerate() {
            let chr_idx = self.tr_chr_idx[tr_idx];
            let strand = self.tr_strand[tr_idx];
            let gene1 = self.tr_gene_idx[tr_idx] + 1;
            for pair in exs.windows(2) {
                let e0 = &pair[0];
                let e1 = &pair[1];
                // STAR: if e1.start <= e0.end+1, exons touch (no intron). In
                // rustar-aligner's 0-based half-open: e1.genome_start <= e0.genome_end.
                if e1.genome_start <= e0.genome_end {
                    continue;
                }
                junctions.push((
                    e0.genome_end,       // sjStart: 0-based first intron base
                    e1.genome_start - 1, // sjEnd: 0-based last intron base
                    chr_idx,
                    strand,
                    gene1,
                ));
            }
        }
        // STAR sorts by sjStart only (funCompareUint2 on the first uint64).
        // Keep it stable so gene list order across duplicates matches
        // transcript-insertion order.
        junctions.sort_by_key(|&(s, _, _, _, _)| s);

        let strand_char = |s: u8| match s {
            1 => '+',
            2 => '-',
            _ => '.',
        };

        // Dedup pass: merge genes across identical (chr, start, end, strand).
        let mut i = 0;
        while i < junctions.len() {
            let (sj_start, sj_end, chr_idx, strand, gene1) = junctions[i];
            let chr_offset = genome.chr_start[chr_idx];
            let start_1based = sj_start + 1 - chr_offset;
            let end_1based = (sj_end + 1) - chr_offset;
            write!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                genome.chr_name[chr_idx],
                start_1based,
                end_1based,
                strand_char(strand),
                gene1
            )
            .map_err(|e| Error::io(e, &path))?;

            // Append genes from subsequent entries with the same key.
            let mut j = i + 1;
            let mut seen: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
            seen.insert(gene1);
            while j < junctions.len() {
                let (s2, e2, c2, st2, g2) = junctions[j];
                if s2 == sj_start && e2 == sj_end && c2 == chr_idx && st2 == strand {
                    if seen.insert(g2) {
                        write!(out, ",{g2}").map_err(|e| Error::io(e, &path))?;
                    }
                    j += 1;
                } else {
                    break;
                }
            }
            writeln!(out).map_err(|e| Error::io(e, &path))?;
            i = j;
        }

        Ok(())
    }

    /// Write `exonGeTrInfo.tab` into `dir`, byte-for-byte matching STAR's
    /// `GTF_transcriptGeneSJ.cpp:33-53` format:
    ///
    /// - Header: total exon count.
    /// - Per-exon line: `exStart\texEnd\texStrand\tgeID\ttrID\n`
    ///
    /// Coordinates are 0-based absolute with `exEnd` INCLUSIVE (STAR's
    /// `exE`). Exons are sorted by `(exStart, exEnd, exStrand, geID, trID)`
    /// — STAR's `funCompareArrays<uint64,5>` order.
    ///
    /// Note: STAR writes the insertion-order transcript index here, not the
    /// sorted one. Its own comment at line 51 calls this "wrong" because
    /// transcripts are re-sorted in `transcriptInfo.tab`. For byte parity
    /// we replicate the behavior.
    pub fn write_exon_ge_tr_info(&self, dir: &Path) -> Result<(), Error> {
        let path = dir.join("exonGeTrInfo.tab");
        let file = std::fs::File::create(&path).map_err(|e| Error::io(e, &path))?;
        let mut out = std::io::BufWriter::new(file);

        // Flatten (exStart, exEnd_inclusive, strand, gene_idx, tr_insertion_idx).
        let mut records: Vec<(u64, u64, u8, u32, u32)> =
            Vec::with_capacity(self.tr_exons.iter().map(|e| e.len()).sum());
        for (tr_idx, exs) in self.tr_exons.iter().enumerate() {
            let strand = self.tr_strand[tr_idx];
            let gene_idx = self.tr_gene_idx[tr_idx];
            for ex in exs {
                records.push((
                    ex.genome_start,
                    ex.genome_end.saturating_sub(1),
                    strand,
                    gene_idx,
                    tr_idx as u32,
                ));
            }
        }
        records.sort_unstable();

        writeln!(out, "{}", records.len()).map_err(|e| Error::io(e, &path))?;
        for (ex_start, ex_end, strand, gene_idx, tr_idx) in records {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                ex_start, ex_end, strand, gene_idx, tr_idx
            )
            .map_err(|e| Error::io(e, &path))?;
        }

        Ok(())
    }

    /// Write `geneInfo.tab` into `dir`, byte-for-byte matching STAR's
    /// `GTF_transcriptGeneSJ.cpp:55-60` format:
    ///
    /// - Header: integer gene count.
    /// - Per-gene line: `geneID\tgeneName\tgeneBiotype\n`
    ///
    /// Order is first-seen (same as STAR, which accumulates genes during
    /// GTF parsing). Empty string if the GTF record had no `gene_name` /
    /// `gene_biotype` attribute.
    pub fn write_gene_info(&self, dir: &Path) -> Result<(), Error> {
        let path = dir.join("geneInfo.tab");
        let file = std::fs::File::create(&path).map_err(|e| Error::io(e, &path))?;
        let mut out = std::io::BufWriter::new(file);

        writeln!(out, "{}", self.gene_ids.len()).map_err(|e| Error::io(e, &path))?;
        for ((id, name), biotype) in self
            .gene_ids
            .iter()
            .zip(&self.gene_names)
            .zip(&self.gene_biotypes)
        {
            writeln!(out, "{}\t{}\t{}", id, name, biotype).map_err(|e| Error::io(e, &path))?;
        }

        Ok(())
    }

    /// Write `exonInfo.tab` into `dir`, byte-for-byte matching STAR's
    /// `GTF_transcriptGeneSJ.cpp:86-112` format:
    ///
    /// - Header: total exon count across all transcripts.
    /// - Per-exon line: `exStart_relative\texEnd_relative\texLenCum\n`
    ///
    /// Coordinates are transcript-relative (exon start/end minus transcript
    /// start), 0-based, with `exEnd` INCLUSIVE (matches STAR's `exE`).
    /// Exons are emitted in STAR's sort order: transcripts in `(tr_start,
    /// tr_end)` order, exons within each transcript in genome order.
    pub fn write_exon_info(&self, dir: &Path) -> Result<(), Error> {
        let path = dir.join("exonInfo.tab");
        let file = std::fs::File::create(&path).map_err(|e| Error::io(e, &path))?;
        let mut out = std::io::BufWriter::new(file);

        let total_exons: usize = self.tr_exons.iter().map(|e| e.len()).sum();
        writeln!(out, "{total_exons}").map_err(|e| Error::io(e, &path))?;

        for &i in &self.tr_order {
            let tr_s = self.tr_start[i];
            for ex in &self.tr_exons[i] {
                // STAR's exEnd is 0-based INCLUSIVE; rustar-aligner stores exclusive.
                let ex_start_rel = ex.genome_start - tr_s;
                let ex_end_rel = ex.genome_end.saturating_sub(1) - tr_s;
                writeln!(out, "{}\t{}\t{}", ex_start_rel, ex_end_rel, ex.ex_len_cum)
                    .map_err(|e| Error::io(e, &path))?;
            }
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Loaders for STAR-compatible transcriptome index files.
// ---------------------------------------------------------------------------

/// Flat exon columns: `(starts_transcript_relative, ends_transcript_relative_inclusive, ex_len_cum)`.
type ExonInfoColumns = (Vec<u64>, Vec<u64>, Vec<u32>);

/// Gene table columns: `(gene_ids, gene_names, gene_biotypes)`.
type GeneInfoColumns = (Vec<String>, Vec<String>, Vec<String>);

type TrInfoColumns = (
    Vec<String>, // tr_ids
    Vec<u64>,    // tr_start (0-based absolute)
    Vec<u64>,    // tr_end   (0-based exclusive, converted from STAR's inclusive)
    Vec<u8>,     // tr_strand
    Vec<u32>,    // tr_exn
    Vec<u32>,    // tr_exi
    Vec<u32>,    // tr_gene_idx
);

/// Read `transcriptInfo.tab` back into column vectors. Converts STAR's
/// 0-based-inclusive `trEnd` to rustar-aligner's 0-based-exclusive convention.
/// Ignores the `trEmax` column (reconstructed from `tr_end` + sort order).
fn read_transcript_info(path: &Path) -> Result<TrInfoColumns, Error> {
    let body = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
    let mut lines = body.lines();

    let header = lines
        .next()
        .ok_or_else(|| Error::Index(format!("{}: missing header line", path.display())))?;
    let n_tr: usize = header
        .trim()
        .parse()
        .map_err(|_| Error::Index(format!("{}: header is not an integer", path.display())))?;

    let mut tr_ids = Vec::with_capacity(n_tr);
    let mut tr_start = Vec::with_capacity(n_tr);
    let mut tr_end = Vec::with_capacity(n_tr);
    let mut tr_strand = Vec::with_capacity(n_tr);
    let mut tr_exn = Vec::with_capacity(n_tr);
    let mut tr_exi = Vec::with_capacity(n_tr);
    let mut tr_gene_idx = Vec::with_capacity(n_tr);

    for (row, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 8 {
            return Err(Error::Index(format!(
                "{}: row {} has {} fields, expected 8",
                path.display(),
                row + 1,
                fields.len()
            )));
        }
        let parse_u64 = |s: &str, col: &str| -> Result<u64, Error> {
            s.parse().map_err(|_| {
                Error::Index(format!(
                    "{}: row {} column {} not an integer: '{}'",
                    path.display(),
                    row + 1,
                    col,
                    s
                ))
            })
        };
        tr_ids.push(fields[0].to_string());
        tr_start.push(parse_u64(fields[1], "trStart")?);
        // STAR's trEnd is 0-based INCLUSIVE; rustar-aligner stores exclusive.
        tr_end.push(parse_u64(fields[2], "trEnd")? + 1);
        tr_strand.push(parse_u64(fields[4], "trStrand")? as u8);
        tr_exn.push(parse_u64(fields[5], "trExN")? as u32);
        tr_exi.push(parse_u64(fields[6], "trExI")? as u32);
        tr_gene_idx.push(parse_u64(fields[7], "trGene")? as u32);
    }

    if tr_ids.len() != n_tr {
        return Err(Error::Index(format!(
            "{}: header says {} transcripts but found {}",
            path.display(),
            n_tr,
            tr_ids.len()
        )));
    }

    Ok((
        tr_ids,
        tr_start,
        tr_end,
        tr_strand,
        tr_exn,
        tr_exi,
        tr_gene_idx,
    ))
}

/// Read `exonInfo.tab` into flat arrays. Start/end are transcript-relative,
/// 0-based; exEnd is INCLUSIVE (STAR's convention).
fn read_exon_info(path: &Path) -> Result<ExonInfoColumns, Error> {
    let body = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
    let mut lines = body.lines();

    let header = lines
        .next()
        .ok_or_else(|| Error::Index(format!("{}: missing header line", path.display())))?;
    let n_ex: usize = header
        .trim()
        .parse()
        .map_err(|_| Error::Index(format!("{}: header is not an integer", path.display())))?;

    let mut starts = Vec::with_capacity(n_ex);
    let mut ends = Vec::with_capacity(n_ex);
    let mut len_cums = Vec::with_capacity(n_ex);

    for (row, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 3 {
            return Err(Error::Index(format!(
                "{}: row {} has {} fields, expected 3",
                path.display(),
                row + 1,
                fields.len()
            )));
        }
        let parse = |s: &str, col: &str| -> Result<u64, Error> {
            s.parse().map_err(|_| {
                Error::Index(format!(
                    "{}: row {} column {} not an integer: '{}'",
                    path.display(),
                    row + 1,
                    col,
                    s
                ))
            })
        };
        starts.push(parse(fields[0], "exStart")?);
        ends.push(parse(fields[1], "exEnd")?);
        len_cums.push(parse(fields[2], "exLenCum")? as u32);
    }

    if starts.len() != n_ex {
        return Err(Error::Index(format!(
            "{}: header says {} exons but found {}",
            path.display(),
            n_ex,
            starts.len()
        )));
    }

    Ok((starts, ends, len_cums))
}

/// Read `geneInfo.tab`. Each record: geneID \t geneName \t geneBiotype.
fn read_gene_info(path: &Path) -> Result<GeneInfoColumns, Error> {
    let body = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
    let mut lines = body.lines();

    let header = lines
        .next()
        .ok_or_else(|| Error::Index(format!("{}: missing header line", path.display())))?;
    let n_ge: usize = header
        .trim()
        .parse()
        .map_err(|_| Error::Index(format!("{}: header is not an integer", path.display())))?;

    let mut ids = Vec::with_capacity(n_ge);
    let mut names = Vec::with_capacity(n_ge);
    let mut biotypes = Vec::with_capacity(n_ge);

    for (row, line) in lines.enumerate() {
        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        if fields.len() != 3 {
            return Err(Error::Index(format!(
                "{}: row {} has {} fields, expected 3",
                path.display(),
                row + 1,
                fields.len()
            )));
        }
        ids.push(fields[0].to_string());
        names.push(fields[1].to_string());
        biotypes.push(fields[2].to_string());
    }

    if ids.len() != n_ge {
        return Err(Error::Index(format!(
            "{}: header says {} genes but found {}",
            path.display(),
            n_ge,
            ids.len()
        )));
    }

    Ok((ids, names, biotypes))
}

// ---------------------------------------------------------------------------
// Projection — STAR's `Transcriptome::quantAlign` + `alignToTranscript`.
// ---------------------------------------------------------------------------

/// Project a genome-space `Transcript` onto every transcript whose coordinates
/// fully contain the alignment.  Returns transcript-space alignments ready for
/// BAM emission.
///
/// Mirrors STAR `Transcriptome_quantAlign.cpp`:
///   * binary-search `tr_starts_sorted` for the greatest `tr_start <=
///     align.genome_start`,
///   * walk back while `tr_end_max_sorted[i] >= align.genome_end`,
///   * for each candidate whose `[tr_start, tr_end]` fully contains the align,
///     call `align_to_one_transcript`.
///
/// `lread` is the total read length and is needed to flip read-space
/// coordinates on reverse-strand transcripts (`read_pos' = Lread - (read_pos +
/// len)`).  For paired-end inputs pass the sum of per-mate lengths.
pub fn align_to_transcripts(
    align: &Transcript,
    idx: &TranscriptomeIndex,
    lread: u32,
) -> Vec<Transcript> {
    let mut out: Vec<Transcript> = Vec::new();
    if idx.n_transcripts() == 0 || align.exons.is_empty() {
        return out;
    }

    let a_start = align.genome_start;
    let a_end = align.genome_end;

    // Binary-search `tr_starts_sorted` for greatest tr_start <= a_start.
    // partition_point returns the first index with tr_start > a_start; subtract 1.
    let upper = idx.tr_starts_sorted.partition_point(|&s| s <= a_start);
    if upper == 0 {
        return out; // a_start is to the left of all transcripts
    }
    let mut sorted_i = upper - 1;

    // Walk backwards while the running-max end still covers this alignment.
    loop {
        if idx.tr_end_max_sorted[sorted_i] < a_end {
            break;
        }
        let tr_idx = idx.tr_order[sorted_i];
        if idx.tr_chr_idx[tr_idx] == align.chr_idx
            && idx.tr_start[tr_idx] <= a_start
            && idx.tr_end[tr_idx] >= a_end
            && let Some(projected) = align_to_one_transcript(align, tr_idx, idx, lread)
        {
            out.push(projected);
        }
        if sorted_i == 0 {
            break;
        }
        sorted_i -= 1;
    }
    out
}

/// Project a single alignment onto a single transcript.  Returns `None` if the
/// alignment is inconsistent with the transcript's exon structure (block
/// extends past an exon boundary, splice boundary doesn't match a transcript
/// junction, etc.).
///
/// Direct port of STAR's `alignToTranscript` (see
/// `source/Transcriptome_quantAlign.cpp:5-89`).
fn align_to_one_transcript(
    align: &Transcript,
    tr_idx: usize,
    idx: &TranscriptomeIndex,
    lread: u32,
) -> Option<Transcript> {
    let tr_exons = &idx.tr_exons[tr_idx];
    let tr_strand = idx.tr_strand[tr_idx];
    let tr_length = idx.tr_length[tr_idx];

    // Find exon containing first block's start.
    let first_block = &align.exons[0];
    let g1 = first_block.genome_start;

    let mut ex = find_containing_exon(tr_exons, g1)?;

    // Build projected exons (in t-space) as we walk the alignment blocks.
    let mut proj_exons: Vec<Exon> = Vec::new();

    for iab in 0..align.exons.len() {
        let block = &align.exons[iab];
        let block_start = block.genome_start;
        let block_end = block.genome_end; // half-open exclusive

        // Mate boundary (STAR canonSJ == -3): the previous block and this
        // block belong to different mates. The next block may sit in any
        // transcript exon — re-locate `ex` via the binary-search helper.
        let crossed_mate_boundary = iab > 0 && align.exons[iab - 1].i_frag != block.i_frag;
        if crossed_mate_boundary {
            ex = find_containing_exon(tr_exons, block_start)?;
        }

        // Block must not extend past the current transcript exon's end.
        if block_end > tr_exons[ex].genome_end {
            return None;
        }
        if block_start < tr_exons[ex].genome_start {
            return None;
        }

        // STAR starts a new projected exon on the first block, after a
        // preceding canonSJ junction (>= 0), or at a mate boundary;
        // insertions coalesce into the previous block.
        let start_new = iab == 0 || crossed_mate_boundary || is_splice_boundary_before(align, iab);

        if start_new {
            // t-space position = ex_len_cum + (block_start - exon_start)
            let tr_offset =
                tr_exons[ex].ex_len_cum as u64 + (block_start - tr_exons[ex].genome_start);
            let len = (block_end - block_start) as usize;
            proj_exons.push(Exon {
                genome_start: tr_offset,
                genome_end: tr_offset + len as u64,
                read_start: block.read_start,
                read_end: block.read_start + len,
                i_frag: block.i_frag,
            });
        } else {
            // Coalesce: extend the last projected exon by this block's length
            // (STAR adds block EX_L directly — does NOT update start position).
            if let Some(last) = proj_exons.last_mut() {
                let len = block_end - block_start;
                last.genome_end += len;
                last.read_end = block.read_end;
            } else {
                return None;
            }
        }

        // Advance `ex` across any splice boundary BEFORE the next block,
        // but only when both blocks belong to the same mate (mate boundaries
        // are handled above at the start of the next iteration).
        if iab + 1 < align.exons.len()
            && align.exons[iab + 1].i_frag == block.i_frag
            && is_splice_boundary_before(align, iab + 1)
        {
            // Require the junction to match a transcript junction.
            let next_block = &align.exons[iab + 1];
            if ex + 1 >= tr_exons.len() {
                return None;
            }
            if block_end != tr_exons[ex].genome_end {
                return None;
            }
            if next_block.genome_start != tr_exons[ex + 1].genome_start {
                return None;
            }
            ex += 1;
        }
    }

    // Apply reverse-strand coordinate flip (STAR canonSJ == -999 branch).
    let projected_is_reverse = if tr_strand == 2 {
        !align.is_reverse
    } else {
        align.is_reverse
    };

    if tr_strand == 2 {
        let tr_len = tr_length as u64;
        let lread_u = lread as u64;
        for e in proj_exons.iter_mut() {
            let len = e.genome_end - e.genome_start;
            let new_g = tr_len - (e.genome_start + len);
            e.genome_start = new_g;
            e.genome_end = new_g + len;

            let read_len = e.read_end - e.read_start;
            let new_r = (lread_u as usize).saturating_sub(e.read_start + read_len);
            e.read_start = new_r;
            e.read_end = new_r + read_len;
        }
        proj_exons.reverse();
    }

    // Build projected CIGAR: drop N operations (splices collapse in t-space);
    // for reverse-strand transcripts, reverse the resulting op sequence.
    let mut proj_cigar: Vec<CigarOp> = align
        .cigar
        .iter()
        .filter(|op| !matches!(op, CigarOp::RefSkip(_)))
        .copied()
        .collect();
    if tr_strand == 2 {
        proj_cigar.reverse();
    }

    // Projected genome bounds = outermost t-space exon positions.
    let proj_start = proj_exons.first().map(|e| e.genome_start).unwrap_or(0);
    let proj_end = proj_exons.last().map(|e| e.genome_end).unwrap_or(0);

    Some(Transcript {
        chr_idx: tr_idx,
        genome_start: proj_start,
        genome_end: proj_end,
        is_reverse: projected_is_reverse,
        exons: proj_exons,
        cigar: proj_cigar,
        score: align.score,
        n_mismatch: align.n_mismatch,
        n_gap: align.n_gap,
        // Splices collapse in t-space; projected records don't carry junction
        // metadata or a read_seq (the SAM writer uses the caller-supplied read).
        n_junction: 0,
        junction_motifs: Vec::new(),
        junction_annotated: Vec::new(),
        read_seq: Vec::new(),
    })
}

/// Filter a genome-space alignment per `--quantTranscriptomeSAMoutput`
/// semantics and return the resulting transcriptome-space projections.
///
/// Mirrors STAR `ReadAlign::quantTranscriptome` (see
/// `source/ReadAlign_quantTranscriptome.cpp:9-66`): reject indels if banned,
/// reject single-mate-only PE hits (caller enforces — rustar-aligner lacks per-block
/// mate tags), extend leading/trailing soft-clips if requested, then project.
///
/// `read_bases_align_orientation` is the read in genome base encoding
/// (A=0,C=1,G=2,T=3,N=4) already reversed/complemented when the alignment is
/// on the reverse strand — STAR's `Read1[roStr==0 ? 0 : 2]`.
pub fn filter_and_project(
    align: &Transcript,
    read_bases_align_orientation: &[u8],
    genome: &Genome,
    idx: &TranscriptomeIndex,
    lread: u32,
    mode: QuantTranscriptomeSAMoutput,
    params: &Parameters,
) -> Vec<Transcript> {
    if !mode.allow_indels() && align.n_gap > 0 {
        return Vec::new();
    }

    let align_for_projection = if mode.allow_softclip() || !has_soft_clip(align) {
        align.clone()
    } else {
        match extend_softclips(align, read_bases_align_orientation, genome, lread, params) {
            Some(extended) => extended,
            None => return Vec::new(),
        }
    };

    align_to_transcripts(&align_for_projection, idx, lread)
}

fn has_soft_clip(align: &Transcript) -> bool {
    align
        .cigar
        .iter()
        .any(|op| matches!(op, CigarOp::SoftClip(n) if *n > 0))
}

/// Extend the 5'/3' soft-clips of `align` back into matched bases, counting
/// mismatches.  Returns `None` if the extension exceeds the mismatch budget.
fn extend_softclips(
    align: &Transcript,
    read_bases_align_orientation: &[u8],
    genome: &Genome,
    lread: u32,
    params: &Parameters,
) -> Option<Transcript> {
    // Determine left / right clip sizes from the CIGAR.
    let (left_clip, right_clip) = align.count_soft_clips();

    let mut n_mm_extra: u32 = 0;

    // Walk the left-clip: bases `[first_exon.read_start - left_clip,
    // first_exon.read_start)` in read space, `[first_exon.genome_start -
    // left_clip, first_exon.genome_start)` in genome space.
    if left_clip > 0
        && let Some(first) = align.exons.first()
    {
        for b in 1..=left_clip as usize {
            if b > first.read_start {
                break;
            }
            if (first.genome_start as usize) < b {
                break;
            }
            let r_idx = first.read_start - b;
            let g_idx = (first.genome_start as usize) - b;
            if r_idx >= read_bases_align_orientation.len() || g_idx >= genome.sequence.len() {
                break;
            }
            let r1 = read_bases_align_orientation[r_idx];
            let g1 = genome.sequence[g_idx];
            if r1 != g1 && r1 < 4 && g1 < 4 {
                n_mm_extra += 1;
            }
        }
    }

    // Walk the right-clip: bases `[last_exon.read_end, last_exon.read_end +
    // right_clip)` in read space, `[last_exon.genome_end, last_exon.genome_end
    // + right_clip)` in genome space.
    if right_clip > 0
        && let Some(last) = align.exons.last()
    {
        for b in 0..right_clip as usize {
            let r_idx = last.read_end + b;
            let g_idx = (last.genome_end as usize) + b;
            if r_idx >= read_bases_align_orientation.len() || g_idx >= genome.sequence.len() {
                break;
            }
            let r1 = read_bases_align_orientation[r_idx];
            let g1 = genome.sequence[g_idx];
            if r1 != g1 && r1 < 4 && g1 < 4 {
                n_mm_extra += 1;
            }
        }
    }

    // Apply STAR's mismatch budget.
    let mismatch_nmax_abs = params.out_filter_mismatch_nmax;
    let mismatch_nmax_rel =
        ((params.out_filter_mismatch_nover_lmax * (lread.saturating_sub(1) as f64)).floor()) as u32;
    let budget = mismatch_nmax_abs.min(mismatch_nmax_rel);
    if align.n_mismatch.saturating_add(n_mm_extra) > budget {
        return None;
    }

    // Construct the extended alignment: remove the soft-clip CIGAR ops and
    // extend the leading/trailing match blocks.
    let mut ext = align.clone();
    ext.n_mismatch = ext.n_mismatch.saturating_add(n_mm_extra);
    if left_clip > 0
        && let Some(first) = ext.exons.first_mut()
    {
        let shift = (left_clip as u64)
            .min(first.read_start as u64)
            .min(first.genome_start);
        first.read_start -= shift as usize;
        first.genome_start -= shift;
    }
    if right_clip > 0
        && let Some(last) = ext.exons.last_mut()
    {
        last.read_end += right_clip as usize;
        last.genome_end += right_clip as u64;
    }

    ext.cigar = rebuild_cigar_without_softclips(&align.cigar, left_clip, right_clip);
    if let Some(first) = ext.exons.first() {
        ext.genome_start = first.genome_start;
    }
    if let Some(last) = ext.exons.last() {
        ext.genome_end = last.genome_end;
    }
    Some(ext)
}

/// Strip the leading/trailing `SoftClip` ops and fold their lengths into the
/// adjacent `Match` ops.  Interior soft-clips (rare — only appear in chimeric
/// contexts) are left alone.
fn rebuild_cigar_without_softclips(
    cigar: &[CigarOp],
    left_clip: u32,
    right_clip: u32,
) -> Vec<CigarOp> {
    let mut out: Vec<CigarOp> = Vec::with_capacity(cigar.len());
    let mut start_idx = 0;
    let mut end_idx = cigar.len();
    if left_clip > 0 && matches!(cigar.first(), Some(CigarOp::SoftClip(_))) {
        start_idx = 1;
    }
    if right_clip > 0 && end_idx > start_idx && matches!(cigar[end_idx - 1], CigarOp::SoftClip(_)) {
        end_idx -= 1;
    }

    let body = &cigar[start_idx..end_idx];
    for (i, op) in body.iter().enumerate() {
        if i == 0 && left_clip > 0 {
            // Fold left_clip into the first op if it's match-like.
            match op {
                CigarOp::Match(n) => out.push(CigarOp::Match(n + left_clip)),
                CigarOp::Equal(n) => out.push(CigarOp::Equal(n + left_clip)),
                _ => {
                    // Extension landed on a non-match op (shouldn't normally
                    // happen).  Emit as Match.
                    out.push(CigarOp::Match(left_clip));
                    out.push(*op);
                }
            }
        } else if i + 1 == body.len() && right_clip > 0 {
            match op {
                CigarOp::Match(n) => out.push(CigarOp::Match(n + right_clip)),
                CigarOp::Equal(n) => out.push(CigarOp::Equal(n + right_clip)),
                _ => {
                    out.push(*op);
                    out.push(CigarOp::Match(right_clip));
                }
            }
        } else {
            out.push(*op);
        }
    }
    out
}

/// Find the transcript exon (by index in `tr_exons`) that contains position
/// `pos` (0-based genome coord).  Returns `None` if `pos` is in an intron or
/// outside the transcript.
fn find_containing_exon(tr_exons: &[TrExon], pos: u64) -> Option<usize> {
    let mut lo = 0usize;
    let mut hi = tr_exons.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if tr_exons[mid].genome_end <= pos {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo < tr_exons.len() && pos >= tr_exons[lo].genome_start && pos < tr_exons[lo].genome_end {
        Some(lo)
    } else {
        None
    }
}

/// Return true if the boundary between `align.exons[iab-1]` and
/// `align.exons[iab]` is a splice (`RefSkip`) rather than an indel.
///
/// rustar-aligner's `Transcript.exons` is a list of read-contiguous match blocks;
/// splices / insertions / deletions all create block boundaries.  We
/// discriminate based on the read-side gap:
///   * `read_end_prev == read_start_curr` AND `genome gap` → deletion OR splice.
///     We call it a splice — the caller verifies the boundary matches a
///     transcript junction and rejects the alignment if it does not.
///   * read gap (insertion) → not a splice (coalesce).
///   * Pure deletion (read contiguous, small genome gap that does NOT match a
///     transcript junction) → handled by the caller rejecting the alignment
///     when the junction check fails.  In practice rustar-aligner produces `Del` ops
///     inside a single exon (no block split for pure deletions because the
///     stitch merge coalesces across Del), so this branch rarely fires.
fn is_splice_boundary_before(align: &Transcript, iab: usize) -> bool {
    if iab == 0 {
        return false;
    }
    let prev = &align.exons[iab - 1];
    let cur = &align.exons[iab];
    // Insertion: read gap between blocks with no genome gap → coalesce.
    if prev.read_end < cur.read_start && prev.genome_end == cur.genome_start {
        return false;
    }
    // Splice / large gap on the genome side, with read-contiguous: treat as
    // potential splice junction (caller validates).
    prev.genome_end < cur.genome_start
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_genome() -> Genome {
        Genome {
            sequence: vec![0u8; 3000],
            n_genome: 3000,
            n_chr_real: 2,
            chr_start: vec![0, 1000, 3000],
            chr_length: vec![1000, 2000],
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        }
    }

    fn make_exon(
        seqname: &str,
        start: u64,
        end: u64,
        strand: char,
        gene_id: &str,
        transcript_id: &str,
    ) -> GtfRecord {
        let mut attrs = HashMap::new();
        attrs.insert("gene_id".to_string(), gene_id.to_string());
        attrs.insert("transcript_id".to_string(), transcript_id.to_string());
        GtfRecord {
            seqname: seqname.to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes: attrs,
        }
    }

    #[test]
    fn single_exon_transcript_metadata() {
        let genome = make_genome();
        let exons = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_ids[0], "T1");
        assert_eq!(idx.gene_ids[idx.tr_gene_idx[0] as usize], "G1");
        assert_eq!(idx.tr_strand[0], 1);
        assert_eq!(idx.tr_chr_idx[0], 0);
        // GTF 1-based 101 → absolute 0-based 100; end 200 → absolute 200 (exclusive)
        assert_eq!(idx.tr_start[0], 100);
        assert_eq!(idx.tr_end[0], 200);
        assert_eq!(idx.tr_length[0], 100);
        assert_eq!(idx.tr_exons[0].len(), 1);
        assert_eq!(idx.tr_exons[0][0].ex_len_cum, 0);
    }

    #[test]
    fn multi_exon_transcript_ex_len_cum() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"), // len 100
            make_exon("chr1", 301, 400, '+', "G1", "T1"), // len 100
            make_exon("chr1", 501, 650, '+', "G1", "T1"), // len 150
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_length[0], 350);

        let ex = &idx.tr_exons[0];
        assert_eq!(ex[0].ex_len_cum, 0);
        assert_eq!(ex[1].ex_len_cum, 100);
        assert_eq!(ex[2].ex_len_cum, 200);

        // Ex_len_cum must be monotonically non-decreasing.
        for w in ex.windows(2) {
            assert!(w[0].ex_len_cum <= w[1].ex_len_cum);
        }

        // tr_start / tr_end from first/last absolute exon.
        assert_eq!(idx.tr_start[0], 100);
        assert_eq!(idx.tr_end[0], 650);
    }

    #[test]
    fn reverse_strand_transcript() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '-', "G2", "T2"),
            make_exon("chr1", 301, 400, '-', "G2", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_strand[0], 2);
        assert_eq!(idx.tr_length[0], 200);
    }

    #[test]
    fn two_transcripts_same_gene() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
            make_exon("chr1", 101, 200, '+', "G1", "T2"),
            make_exon("chr1", 301, 500, '+', "G1", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 2);
        assert_eq!(idx.gene_ids[idx.tr_gene_idx[0] as usize], "G1");
        assert_eq!(idx.gene_ids[idx.tr_gene_idx[1] as usize], "G1");
        // Different lengths: T1 = 100+100 = 200, T2 = 100+200 = 300
        assert_eq!(idx.tr_length[0], 200);
        assert_eq!(idx.tr_length[1], 300);
    }

    fn make_exon_with_attrs(
        seqname: &str,
        start: u64,
        end: u64,
        strand: char,
        transcript_id: &str,
        gene_attrs: &[(&str, &str)],
    ) -> GtfRecord {
        let mut attrs = HashMap::new();
        attrs.insert("transcript_id".to_string(), transcript_id.to_string());
        for (k, v) in gene_attrs {
            attrs.insert(k.to_string(), v.to_string());
        }
        GtfRecord {
            seqname: seqname.to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes: attrs,
        }
    }

    #[test]
    fn gene_interning_first_seen_wins() {
        let genome = make_genome();
        // Two transcripts on gene G1, one on gene G2. G1 name set via T1;
        // T2 has no name — interning should NOT overwrite G1's name.
        let exons = vec![
            make_exon_with_attrs(
                "chr1",
                101,
                200,
                '+',
                "T1",
                &[
                    ("gene_id", "G1"),
                    ("gene_name", "GENE1"),
                    ("gene_biotype", "protein_coding"),
                ],
            ),
            make_exon_with_attrs("chr1", 301, 400, '+', "T2", &[("gene_id", "G1")]),
            make_exon_with_attrs(
                "chr1",
                501,
                600,
                '+',
                "T3",
                &[("gene_id", "G2"), ("gene_name", "GENE2")],
            ),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.gene_ids, vec!["G1".to_string(), "G2".to_string()]);
        assert_eq!(
            idx.gene_names,
            vec!["GENE1".to_string(), "GENE2".to_string()]
        );
        // G2 has no gene_biotype attribute — STAR-faithful fallback is
        // the literal "MissingGeneType".
        assert_eq!(
            idx.gene_biotypes,
            vec!["protein_coding".to_string(), "MissingGeneType".to_string()]
        );
        assert_eq!(idx.tr_gene_idx, vec![0, 0, 1]);
    }

    #[test]
    fn transcript_info_tab_byte_format() {
        let genome = make_genome();
        // Two non-overlapping transcripts, T1 before T2. Same gene G1.
        // T1: chr1 [100..200) forward, 1 exon of length 100. End inclusive = 199.
        // T2: chr1 [300..500) forward, 1 exon of length 200. End inclusive = 499.
        let exons = vec![
            make_exon_with_attrs("chr1", 101, 200, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 301, 500, '+', "T2", &[("gene_id", "G1")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_transcript_info(dir.path()).unwrap();

        let body = std::fs::read_to_string(dir.path().join("transcriptInfo.tab")).unwrap();
        // Header + two records.
        // Format: trID \t trStart \t trEnd(incl) \t trEmax \t trStrand \t trExN \t trExI \t trGene
        // Sorted position 0 (T1): trEnd = 199, trEmax = 199 (initial = first end).
        // Sorted position 1 (T2): trEnd = 499, trEmax = 199 (running max excludes current).
        assert_eq!(
            body,
            "2\n\
             T1\t100\t199\t199\t1\t1\t0\t0\n\
             T2\t300\t499\t199\t1\t1\t1\t0\n"
        );
    }

    #[test]
    fn transcript_info_tab_reverse_strand() {
        let genome = make_genome();
        let exons = vec![make_exon_with_attrs(
            "chr1",
            101,
            200,
            '-',
            "T1",
            &[("gene_id", "G1")],
        )];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_transcript_info(dir.path()).unwrap();
        let body = std::fs::read_to_string(dir.path().join("transcriptInfo.tab")).unwrap();
        // Reverse strand in STAR's GTF.cpp encoding = 2.
        assert_eq!(body, "1\nT1\t100\t199\t199\t2\t1\t0\t0\n");
    }

    #[test]
    fn transcript_info_tab_emax_running_max_excludes_current() {
        let genome = make_genome();
        // Three transcripts ordered by sort: T_small (end 199), T_big (end 699), T_mid (end 499).
        // Sort by (trStart, trEnd) on different start positions so order is
        // deterministic by start alone.
        let exons = vec![
            make_exon_with_attrs("chr1", 101, 200, '+', "T_small", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 301, 700, '+', "T_big", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 801, 900, '+', "T_mid", &[("gene_id", "G1")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_transcript_info(dir.path()).unwrap();
        let body = std::fs::read_to_string(dir.path().join("transcriptInfo.tab")).unwrap();
        // Sorted-order emax values:
        //   pos 0 (T_small, end 199): emax = 199 (initial = first end)
        //   pos 1 (T_big, end 699):   emax = 199 (running max of {T_small})
        //   pos 2 (T_mid, end 899):   emax = 699 (running max of {T_small, T_big})
        assert_eq!(
            body,
            "3\n\
             T_small\t100\t199\t199\t1\t1\t0\t0\n\
             T_big\t300\t699\t199\t1\t1\t1\t0\n\
             T_mid\t800\t899\t699\t1\t1\t2\t0\n"
        );
    }

    #[test]
    fn roundtrip_write_then_load_matches_in_memory_index() {
        let genome = make_genome();
        let exons = vec![
            make_exon_with_attrs(
                "chr1",
                101,
                200,
                '+',
                "T1",
                &[
                    ("gene_id", "G1"),
                    ("gene_name", "GENE1"),
                    ("gene_biotype", "protein_coding"),
                ],
            ),
            make_exon_with_attrs("chr1", 301, 400, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs(
                "chr1",
                601,
                700,
                '-',
                "T2",
                &[("gene_id", "G2"), ("gene_name", "GENE2")],
            ),
        ];
        let built = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();

        let dir = tempfile::tempdir().unwrap();
        built.write_transcript_info(dir.path()).unwrap();
        built.write_exon_info(dir.path()).unwrap();
        built.write_gene_info(dir.path()).unwrap();

        let loaded = TranscriptomeIndex::from_index_dir(dir.path(), &genome).unwrap();

        // Same transcripts, same exon structure. After load the order is
        // sorted (which equals insertion order here — starts are distinct).
        assert_eq!(loaded.tr_ids, built.tr_ids);
        assert_eq!(loaded.tr_start, built.tr_start);
        assert_eq!(loaded.tr_end, built.tr_end);
        assert_eq!(loaded.tr_strand, built.tr_strand);
        assert_eq!(loaded.tr_gene_idx, built.tr_gene_idx);
        assert_eq!(loaded.tr_chr_idx, built.tr_chr_idx);
        assert_eq!(loaded.tr_length, built.tr_length);
        assert_eq!(loaded.tr_exons, built.tr_exons);
        assert_eq!(loaded.gene_ids, built.gene_ids);
        assert_eq!(loaded.gene_names, built.gene_names);
        assert_eq!(loaded.gene_biotypes, built.gene_biotypes);
    }

    #[test]
    fn load_handles_empty_gene_name_and_biotype() {
        let genome = make_genome();
        let exons = vec![make_exon_with_attrs(
            "chr1",
            101,
            200,
            '+',
            "T1",
            &[("gene_id", "G1")],
        )];
        let built = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        built.write_transcript_info(dir.path()).unwrap();
        built.write_exon_info(dir.path()).unwrap();
        built.write_gene_info(dir.path()).unwrap();

        let loaded = TranscriptomeIndex::from_index_dir(dir.path(), &genome).unwrap();
        // STAR fallbacks: gene_name → gene_id, gene_biotype → "MissingGeneType".
        assert_eq!(loaded.gene_names, vec!["G1".to_string()]);
        assert_eq!(loaded.gene_biotypes, vec!["MissingGeneType".to_string()]);
    }

    #[test]
    fn sjdb_list_from_gtf_tab_byte_format() {
        let genome = make_genome();
        // T1: exons [100,200) and [300,400) → intron [200,300), 1-based [201,299].
        // T2: exons [600,700) and [800,900) → intron [700,800), 1-based [701,799].
        let exons = vec![
            make_exon_with_attrs("chr1", 101, 200, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 301, 400, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 601, 700, '-', "T2", &[("gene_id", "G2")]),
            make_exon_with_attrs("chr1", 801, 900, '-', "T2", &[("gene_id", "G2")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_sjdb_list_from_gtf(dir.path(), &genome).unwrap();
        let body = std::fs::read_to_string(dir.path().join("sjdbList.fromGTF.out.tab")).unwrap();
        // Gene indices are 1-based: G1 → 1, G2 → 2.
        assert_eq!(
            body,
            "chr1\t201\t300\t+\t1\n\
             chr1\t701\t800\t-\t2\n"
        );
    }

    #[test]
    fn sjdb_list_dedups_and_merges_genes() {
        let genome = make_genome();
        // Two transcripts on DIFFERENT genes share the same junction.
        // T1 gene G1 (idx 0, 1-based 1): exons [100,200) + [300,400)
        // T2 gene G2 (idx 1, 1-based 2): exons [100,200) + [300,400)
        let exons = vec![
            make_exon_with_attrs("chr1", 101, 200, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 301, 400, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 101, 200, '+', "T2", &[("gene_id", "G2")]),
            make_exon_with_attrs("chr1", 301, 400, '+', "T2", &[("gene_id", "G2")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_sjdb_list_from_gtf(dir.path(), &genome).unwrap();
        let body = std::fs::read_to_string(dir.path().join("sjdbList.fromGTF.out.tab")).unwrap();
        // Single line, genes comma-separated.
        assert_eq!(body, "chr1\t201\t300\t+\t1,2\n");
    }

    #[test]
    fn sjdb_list_single_exon_transcript_has_no_junctions() {
        let genome = make_genome();
        let exons = vec![make_exon_with_attrs(
            "chr1",
            101,
            200,
            '+',
            "T1",
            &[("gene_id", "G1")],
        )];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_sjdb_list_from_gtf(dir.path(), &genome).unwrap();
        let body = std::fs::read_to_string(dir.path().join("sjdbList.fromGTF.out.tab")).unwrap();
        assert_eq!(body, "");
    }

    #[test]
    fn exon_ge_tr_info_tab_byte_format() {
        let genome = make_genome();
        // T1 (insertion 0): forward strand, 2 exons at [100,200) and [300,400).
        // T2 (insertion 1): reverse strand, 1 exon at [150,250).
        // T3 (insertion 2): forward, 1 exon at [100,200) — overlaps T1's first exon.
        // Sort by (start, end, strand, gene, trIdx):
        //   (100, 199, 1, 0, 0) — T1's exon 1
        //   (100, 199, 1, 2, 2) — T3 (same start/end/strand, different gene)
        //   (150, 249, 2, 1, 1) — T2
        //   (300, 399, 1, 0, 0) — T1's exon 2
        let exons = vec![
            make_exon_with_attrs("chr1", 101, 200, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 301, 400, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 151, 250, '-', "T2", &[("gene_id", "G2")]),
            make_exon_with_attrs("chr1", 101, 200, '+', "T3", &[("gene_id", "G3")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_exon_ge_tr_info(dir.path()).unwrap();
        let body = std::fs::read_to_string(dir.path().join("exonGeTrInfo.tab")).unwrap();
        assert_eq!(
            body,
            "4\n\
             100\t199\t1\t0\t0\n\
             100\t199\t1\t2\t2\n\
             150\t249\t2\t1\t1\n\
             300\t399\t1\t0\t0\n"
        );
    }

    #[test]
    fn gene_info_tab_byte_format() {
        let genome = make_genome();
        let exons = vec![
            make_exon_with_attrs(
                "chr1",
                101,
                200,
                '+',
                "T1",
                &[
                    ("gene_id", "G1"),
                    ("gene_name", "GENE1"),
                    ("gene_biotype", "protein_coding"),
                ],
            ),
            make_exon_with_attrs(
                "chr1",
                301,
                400,
                '+',
                "T2",
                &[("gene_id", "G2"), ("gene_name", "GENE2")],
            ),
            make_exon_with_attrs("chr1", 501, 600, '+', "T3", &[("gene_id", "G3")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_gene_info(dir.path()).unwrap();
        let body = std::fs::read_to_string(dir.path().join("geneInfo.tab")).unwrap();
        // STAR-faithful fallbacks: missing gene_name → gene_id string,
        // missing gene_biotype → literal "MissingGeneType".
        assert_eq!(
            body,
            "3\n\
             G1\tGENE1\tprotein_coding\n\
             G2\tGENE2\tMissingGeneType\n\
             G3\tG3\tMissingGeneType\n"
        );
    }

    #[test]
    fn exon_info_tab_byte_format() {
        let genome = make_genome();
        // T1: 2 exons at [100,200) and [300,400), total 200 bases.
        // T2: 1 exon at [500,650), 150 bases.
        let exons = vec![
            make_exon_with_attrs("chr1", 101, 200, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 301, 400, '+', "T1", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 501, 650, '+', "T2", &[("gene_id", "G1")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_exon_info(dir.path()).unwrap();
        let body = std::fs::read_to_string(dir.path().join("exonInfo.tab")).unwrap();
        // Order: T1 (sorted pos 0) then T2. Each exon: start_rel \t end_rel(incl) \t exLenCum
        //   T1 tr_start = 100. Exon 1: 100-100=0, 200-1-100=99, lenCum=0.
        //                      Exon 2: 300-100=200, 400-1-100=299, lenCum=100.
        //   T2 tr_start = 500. Exon 1: 500-500=0, 650-1-500=149, lenCum=0.
        assert_eq!(
            body,
            "3\n\
             0\t99\t0\n\
             200\t299\t100\n\
             0\t149\t0\n"
        );
    }

    #[test]
    fn exon_info_respects_sort_order() {
        let genome = make_genome();
        // Insert T_late first, T_early second. In sorted order T_early comes first.
        let exons = vec![
            make_exon_with_attrs("chr1", 501, 600, '+', "T_late", &[("gene_id", "G1")]),
            make_exon_with_attrs("chr1", 101, 200, '+', "T_early", &[("gene_id", "G1")]),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        let dir = tempfile::tempdir().unwrap();
        idx.write_exon_info(dir.path()).unwrap();
        let body = std::fs::read_to_string(dir.path().join("exonInfo.tab")).unwrap();
        // T_early first (sorted pos 0): 0 \t 99 \t 0
        // T_late second:                0 \t 99 \t 0
        assert_eq!(
            body,
            "2\n\
             0\t99\t0\n\
             0\t99\t0\n"
        );
    }

    #[test]
    fn tr_exi_sorted_cumulative() {
        let genome = make_genome();
        // Three transcripts with non-overlapping starts: T1 [100..600),
        // T2 [700..900), T3 [1100..1300) — sorted order = insertion order.
        // Counts: T1=3, T2=2, T3=1. tr_exi in sorted order = [0, 3, 5].
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
            make_exon("chr1", 501, 600, '+', "G1", "T1"),
            make_exon("chr1", 701, 800, '+', "G1", "T2"),
            make_exon("chr1", 801, 900, '+', "G1", "T2"),
            make_exon("chr2", 101, 200, '+', "G2", "T3"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.tr_exi, vec![0, 3, 5]);
    }

    #[test]
    fn tr_exi_respects_sort_order_not_insertion_order() {
        let genome = make_genome();
        // T_late inserted first but sorts LAST (starts 501); T_early inserted
        // second and sorts FIRST (starts 101). tr_exi must reflect sorted
        // cumulative count: T_early gets 0, T_late gets 1 (T_early's exon count).
        let exons = vec![
            make_exon("chr1", 501, 600, '+', "G1", "T_late"),
            make_exon("chr1", 101, 200, '+', "G1", "T_early"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        // Insertion order: [T_late, T_early]
        assert_eq!(
            idx.tr_ids,
            vec!["T_late".to_string(), "T_early".to_string()]
        );
        // Sorted order: T_early (pos 0), T_late (pos 1)
        // tr_order maps sorted → insertion: tr_order[0] = 1 (T_early), tr_order[1] = 0 (T_late)
        assert_eq!(idx.tr_order, vec![1, 0]);
        // tr_exi by insertion position: T_late (insertion 0, sort pos 1) = 1,
        // T_early (insertion 1, sort pos 0) = 0.
        assert_eq!(idx.tr_exi, vec![1, 0]);
    }

    #[test]
    fn unknown_chromosome_skipped() {
        let genome = make_genome();
        let exons = vec![make_exon("chrX", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 0);
    }

    #[test]
    fn inconsistent_strand_skipped() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '-', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 0);
    }

    fn make_align(
        chr_idx: usize,
        is_reverse: bool,
        exons: Vec<(u64, u64, usize, usize)>,
        cigar: Vec<CigarOp>,
    ) -> Transcript {
        let proj_exons: Vec<Exon> = exons
            .into_iter()
            .map(|(gs, ge, rs, re)| Exon {
                genome_start: gs,
                genome_end: ge,
                read_start: rs,
                read_end: re,
                i_frag: 0,
            })
            .collect();
        let gs = proj_exons.first().map(|e| e.genome_start).unwrap_or(0);
        let ge = proj_exons.last().map(|e| e.genome_end).unwrap_or(0);
        Transcript {
            chr_idx,
            genome_start: gs,
            genome_end: ge,
            is_reverse,
            exons: proj_exons,
            cigar,
            score: 100,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        }
    }

    #[test]
    fn project_single_exon_align_into_single_exon_transcript() {
        let genome = make_genome();
        // Transcript: chr1 [100, 200) forward.
        let gtf = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align fully inside exon: genome [110, 150), read [0, 40).
        let align = make_align(0, false, vec![(110, 150, 0, 40)], vec![CigarOp::Match(40)]);
        let results = align_to_transcripts(&align, &idx, 40);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.chr_idx, 0); // transcript index
        assert_eq!(r.is_reverse, false);
        assert_eq!(r.exons.len(), 1);
        // t-space offset = 0 (ex_len_cum) + (110 - 100) = 10
        assert_eq!(r.exons[0].genome_start, 10);
        assert_eq!(r.exons[0].genome_end, 50);
        assert_eq!(r.exons[0].read_start, 0);
        assert_eq!(r.exons[0].read_end, 40);
        // CIGAR must have no N
        assert!(r.cigar.iter().all(|op| !matches!(op, CigarOp::RefSkip(_))));
        assert_eq!(r.cigar.len(), 1);
    }

    #[test]
    fn project_two_exon_align_matching_junction() {
        let genome = make_genome();
        // Transcript T1: chr1 [100, 200) + [300, 400) forward. tr_length = 200.
        let gtf = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align: two blocks (50 M each), junction matches transcript junction.
        // Genome: [150, 200) + [300, 350); read [0, 50) + [50, 100).
        let align = make_align(
            0,
            false,
            vec![(150, 200, 0, 50), (300, 350, 50, 100)],
            vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(100),
                CigarOp::Match(50),
            ],
        );
        let results = align_to_transcripts(&align, &idx, 100);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        // Two t-space exons, no N in CIGAR.
        assert_eq!(r.exons.len(), 2);
        // First exon: t-space [50, 100)
        assert_eq!(r.exons[0].genome_start, 50);
        assert_eq!(r.exons[0].genome_end, 100);
        // Second exon: t-space [100, 150) — starts right after prev (splice collapsed)
        assert_eq!(r.exons[1].genome_start, 100);
        assert_eq!(r.exons[1].genome_end, 150);
        // CIGAR: no N
        assert!(r.cigar.iter().all(|op| !matches!(op, CigarOp::RefSkip(_))));
        assert_eq!(r.cigar.len(), 2);
    }

    #[test]
    fn project_mismatched_junction_fails() {
        let genome = make_genome();
        let gtf = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align with junction NOT matching transcript — splice ends at 195 instead of 200.
        let align = make_align(
            0,
            false,
            vec![(150, 195, 0, 45), (305, 350, 45, 90)],
            vec![
                CigarOp::Match(45),
                CigarOp::RefSkip(110),
                CigarOp::Match(45),
            ],
        );
        let results = align_to_transcripts(&align, &idx, 90);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn project_onto_reverse_strand_transcript() {
        let genome = make_genome();
        // Reverse transcript: chr1 [100, 200) - strand. tr_length = 100.
        let gtf = vec![make_exon("chr1", 101, 200, '-', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align forward on genome [120, 160), read [0, 40).
        let align = make_align(0, false, vec![(120, 160, 0, 40)], vec![CigarOp::Match(40)]);
        let results = align_to_transcripts(&align, &idx, 40);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        // is_reverse flipped (transcript strand == 2, align was false → result true)
        assert_eq!(r.is_reverse, true);
        // t-space position after flip:
        //   pre-flip: genome_start=20, length=40 → new_g = 100 - (20 + 40) = 40
        //   read_start=0, read_len=40 → new_r = 40 - (0 + 40) = 0
        assert_eq!(r.exons.len(), 1);
        assert_eq!(r.exons[0].genome_start, 40);
        assert_eq!(r.exons[0].genome_end, 80);
        assert_eq!(r.exons[0].read_start, 0);
        assert_eq!(r.exons[0].read_end, 40);
    }

    #[test]
    fn project_multi_exon_align_onto_longer_transcript() {
        let genome = make_genome();
        // 3-exon transcript: [100,200) + [300,400) + [500,600). tr_length = 300.
        let gtf = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
            make_exon("chr1", 501, 600, '+', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align spans just exon 2 and 3 (skipping first exon): genome [350, 400) + [500, 550).
        let align = make_align(
            0,
            false,
            vec![(350, 400, 0, 50), (500, 550, 50, 100)],
            vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(100),
                CigarOp::Match(50),
            ],
        );
        let results = align_to_transcripts(&align, &idx, 100);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.exons.len(), 2);
        // First t-space exon: ex_len_cum[1]=100, offset within exon = 350-300=50 → t-space start 150
        assert_eq!(r.exons[0].genome_start, 150);
        assert_eq!(r.exons[0].genome_end, 200);
        // Second t-space exon: ex_len_cum[2]=200, offset = 500-500=0 → t-space start 200
        assert_eq!(r.exons[1].genome_start, 200);
        assert_eq!(r.exons[1].genome_end, 250);
    }

    #[test]
    fn project_past_transcript_end_fails() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align extends past transcript end (genome [150, 250), transcript ends at 200).
        let align = make_align(
            0,
            false,
            vec![(150, 250, 0, 100)],
            vec![CigarOp::Match(100)],
        );
        let results = align_to_transcripts(&align, &idx, 100);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn project_before_all_transcripts_returns_empty() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 501, 600, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        let align = make_align(0, false, vec![(100, 150, 0, 50)], vec![CigarOp::Match(50)]);
        let results = align_to_transcripts(&align, &idx, 50);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn project_onto_multiple_overlapping_transcripts() {
        let genome = make_genome();
        // Two transcripts both containing the align:
        //  T1: [100, 400) single exon
        //  T2: [100, 400) single exon, different ID
        let gtf = vec![
            make_exon("chr1", 101, 400, '+', "G1", "T1"),
            make_exon("chr1", 101, 400, '+', "G2", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        let align = make_align(0, false, vec![(200, 250, 0, 50)], vec![CigarOp::Match(50)]);
        let results = align_to_transcripts(&align, &idx, 50);
        assert_eq!(results.len(), 2);
    }

    // ---- Subtask 3: filter-mode tests ----

    fn default_params() -> Parameters {
        use clap::Parser;
        Parameters::parse_from(vec!["rustar-aligner", "--readFilesIn", "r.fq"])
    }

    #[test]
    fn mode_from_str_all_three() {
        use std::str::FromStr;
        assert_eq!(
            QuantTranscriptomeSAMoutput::from_str("BanSingleEnd_BanIndels_ExtendSoftclip").unwrap(),
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip,
        );
        assert_eq!(
            QuantTranscriptomeSAMoutput::from_str("BanSingleEnd").unwrap(),
            QuantTranscriptomeSAMoutput::BanSingleEnd,
        );
        assert_eq!(
            QuantTranscriptomeSAMoutput::from_str("BanSingleEnd_ExtendSoftclip").unwrap(),
            QuantTranscriptomeSAMoutput::BanSingleEndExtendSoftclip,
        );
        assert!(QuantTranscriptomeSAMoutput::from_str("garbage").is_err());
    }

    #[test]
    fn mode_flags() {
        assert!(!QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip.allow_indels());
        assert!(!QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip.allow_softclip());
        assert!(QuantTranscriptomeSAMoutput::BanSingleEnd.allow_indels());
        assert!(QuantTranscriptomeSAMoutput::BanSingleEnd.allow_softclip());
        assert!(QuantTranscriptomeSAMoutput::BanSingleEndExtendSoftclip.allow_indels());
        assert!(!QuantTranscriptomeSAMoutput::BanSingleEndExtendSoftclip.allow_softclip());
    }

    #[test]
    fn filter_default_rejects_indels() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 101, 300, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();
        let mut align = make_align(0, false, vec![(110, 150, 0, 40)], vec![CigarOp::Match(40)]);
        align.n_gap = 1; // simulate insertion/deletion
        let params = default_params();
        let read = vec![0u8; 40];
        let results = filter_and_project(
            &align,
            &read,
            &genome,
            &idx,
            40,
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip,
            &params,
        );
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn filter_keeps_indels_when_allowed() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 101, 300, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();
        let mut align = make_align(0, false, vec![(110, 150, 0, 40)], vec![CigarOp::Match(40)]);
        align.n_gap = 1;
        let params = default_params();
        let read = vec![0u8; 40];
        let results = filter_and_project(
            &align,
            &read,
            &genome,
            &idx,
            40,
            QuantTranscriptomeSAMoutput::BanSingleEnd,
            &params,
        );
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn filter_extends_left_softclip_with_zero_mismatches() {
        // Build a custom genome with known content so we can construct a
        // read whose soft-clipped bases match the adjacent genome.
        let mut seq = vec![4u8; 1000];
        // Place pattern "AAAA" at genome [100, 104) — clip region
        for i in 100..104 {
            seq[i] = 0; // A
        }
        // Aligned region [104, 144) — fill with zeros (A) so read bases match
        for i in 104..144 {
            seq[i] = 0;
        }
        let genome = Genome {
            sequence: seq,
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        let gtf = vec![make_exon("chr1", 1, 1000, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        // Align: 4S + 40M starting at genome 104, read [4..44).
        let align = make_align(
            0,
            false,
            vec![(104, 144, 4, 44)],
            vec![CigarOp::SoftClip(4), CigarOp::Match(40)],
        );
        // Read is 44 bases of A (0s) — matches all of genome [100, 144).
        let read = vec![0u8; 44];

        let params = default_params();
        let results = filter_and_project(
            &align,
            &read,
            &genome,
            &idx,
            44,
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip,
            &params,
        );
        assert_eq!(results.len(), 1);
        let r = &results[0];
        // CIGAR should now be a single 44M (4S folded into leading M).
        assert_eq!(r.cigar.len(), 1);
        match r.cigar[0] {
            CigarOp::Match(n) => assert_eq!(n, 44),
            _ => panic!("expected Match(44)"),
        }
    }

    #[test]
    fn filter_extends_softclip_too_many_mismatches_rejects() {
        // Left clip is 4 bases, all mismatches.  n_mismatch = 0 to start.
        // With very tight out_filter_mismatch_nmax, the alignment is rejected.
        let mut seq = vec![4u8; 1000];
        // Clip region [100, 104): all zeros (A)
        for i in 100..104 {
            seq[i] = 0;
        }
        // Aligned region [104, 144): all zeros
        for i in 104..144 {
            seq[i] = 0;
        }
        let genome = Genome {
            sequence: seq,
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };
        let gtf = vec![make_exon("chr1", 1, 1000, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        let align = make_align(
            0,
            false,
            vec![(104, 144, 4, 44)],
            vec![CigarOp::SoftClip(4), CigarOp::Match(40)],
        );
        // Read clip bases are all T (3) — mismatch against genome A (0)
        let mut read = vec![0u8; 44];
        for b in read.iter_mut().take(4) {
            *b = 3;
        }

        let mut params = default_params();
        params.out_filter_mismatch_nmax = 2; // budget = 2
        let results = filter_and_project(
            &align,
            &read,
            &genome,
            &idx,
            44,
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip,
            &params,
        );
        // 4 extension mismatches > 2 → rejected
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn filter_mode_single_end_keeps_softclip_as_is() {
        let genome = make_genome();
        let gtf = vec![make_exon("chr1", 101, 300, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&gtf, &genome).unwrap();

        let align = make_align(
            0,
            false,
            vec![(110, 150, 4, 44)],
            vec![CigarOp::SoftClip(4), CigarOp::Match(40)],
        );
        let read = vec![0u8; 44];
        let params = default_params();

        // Mode BanSingleEnd → keep soft-clips as-is (no extension).
        let results = filter_and_project(
            &align,
            &read,
            &genome,
            &idx,
            44,
            QuantTranscriptomeSAMoutput::BanSingleEnd,
            &params,
        );
        assert_eq!(results.len(), 1);
        // CIGAR preserves the soft-clip.
        assert!(
            results[0]
                .cigar
                .iter()
                .any(|op| matches!(op, CigarOp::SoftClip(_)))
        );
    }

    #[test]
    fn tr_end_max_sorted_is_running_max() {
        let genome = make_genome();
        let exons = vec![
            // T1 chr1 starts at 100, ends at 200
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            // T2 chr1 starts at 150, ends at 500 (encloses T1)
            make_exon("chr1", 151, 500, '+', "G1", "T2"),
            // T3 chr1 starts at 200, ends at 300 (nested)
            make_exon("chr1", 201, 300, '+', "G1", "T3"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 3);

        // tr_order must be sorted by (start, end)
        let sorted_starts: Vec<u64> = idx.tr_starts_sorted.clone();
        let mut check = sorted_starts.clone();
        check.sort();
        assert_eq!(sorted_starts, check);

        // tr_end_max_sorted must be monotonically non-decreasing
        for w in idx.tr_end_max_sorted.windows(2) {
            assert!(w[0] <= w[1]);
        }
        // Last entry equals overall maximum tr_end
        let overall_max = *idx.tr_end.iter().max().unwrap();
        assert_eq!(*idx.tr_end_max_sorted.last().unwrap(), overall_max);
    }
}
