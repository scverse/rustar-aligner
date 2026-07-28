/// Splice junction annotation and tracking
///
/// This module handles:
/// - GTF file parsing for gene/transcript/exon annotations
/// - Building a junction database from annotated exons
/// - Junction lookup during alignment (annotated vs novel)
/// - Junction statistics collection for SJ.out.tab output
pub(crate) mod chr_start_end;
pub(crate) mod gtf;
mod sj_output;
pub mod sjdb_insert;

pub use sj_output::SpliceJunctionStats;
pub(crate) use sj_output::{SjKey, encode_motif};
pub use sjdb_insert::PreparedJunction;

use crate::params::Parameters;

use crate::error::Error;
use crate::genome::Genome;
use std::collections::HashMap;
use std::path::Path;

/// Key for junction lookup: (chr_idx, intron_start, intron_end, strand).
///
/// `intron_start` / `intron_end` are genome-absolute 0-based positions of
/// the first and last intronic bases, matching the convention used by
/// `PreparedJunction`, `SpliceJunctionStats`, and the alignment-time
/// `genome_pos` variables.
#[derive(Hash, Eq, PartialEq, Clone, Debug)]
struct JunctionKey {
    chr_idx: usize,
    intron_start: u64,
    intron_end: u64,
    strand: u8, // 0=unknown, 1=+, 2=-
}

/// Information about a splice junction
#[derive(Debug, Clone)]
pub struct JunctionInfo {
    pub annotated: bool,
    // Future: gene_id, transcript_ids for provenance tracking
}

/// Key for novel junction insertion (public for two-pass mode)
#[derive(Hash, Eq, PartialEq, Clone, Debug)]
pub struct NovelJunctionKey {
    pub chr_idx: usize,
    pub intron_start: u64,
    pub intron_end: u64,
    pub strand: u8, // 0=unknown, 1=+, 2=-
}

/// Splice junction database built from GTF annotations
#[derive(Clone)]
pub struct SpliceJunctionDb {
    /// Map: (chr_idx, intron_start, intron_end, strand) → annotated
    junctions: HashMap<JunctionKey, JunctionInfo>,
    /// STAR's `mapGen.sjdb*` arrays: the annotated junctions in the exact
    /// order they occupy the Gsj buffer, sorted by `(stored_start,
    /// stored_end)`. Index `i` here is STAR's `sjA` / `sjdbInd`, i.e. the
    /// value [`sjdb_insert::decode_gsj_hit`] derives from a Gsj SA hit.
    ///
    /// Empty when the database was built without motif/shift metadata; in
    /// that case [`find`](Self::find) always returns `None` and the
    /// annotated-junction fast paths simply do not engage.
    table: Vec<PreparedJunction>,
}

impl SpliceJunctionDb {
    /// Create empty database (for no-GTF mode)
    pub fn empty() -> Self {
        Self {
            junctions: HashMap::new(),
            table: Vec::new(),
        }
    }

    /// Build from an already sorted-and-deduplicated [`PreparedJunction`]
    /// list — the form `sjdbInfo.txt` is read back into and the form
    /// `genomeGenerate` produces before building the Gsj buffer.
    ///
    /// The HashMap is keyed on the *stored* (post-`sjdbPrepare`) donor and
    /// acceptor coordinates, which is what the stitch-time scan produces.
    pub fn from_prepared(prepared: Vec<PreparedJunction>) -> Self {
        let mut junctions = HashMap::with_capacity(prepared.len());
        for j in &prepared {
            junctions.insert(
                JunctionKey {
                    chr_idx: j.chr_idx,
                    intron_start: j.stored_start(),
                    intron_end: j.stored_end(),
                    strand: j.strand,
                },
                JunctionInfo { annotated: true },
            );
        }
        Self {
            junctions,
            table: prepared,
        }
    }

    /// Replace the `sjA`-addressable table without touching the annotated
    /// lookup map.
    ///
    /// Needed because `sj_a` tags come from [`sjdb_insert::decode_gsj_hit`],
    /// which indexes the junction list stored in the index (`sjdbInfo.txt`).
    /// When a GTF is *also* supplied at align time the map is rebuilt from
    /// that GTF, but the table must keep addressing the index's array or the
    /// tags would point at the wrong junctions.
    pub fn set_table(&mut self, prepared: Vec<PreparedJunction>) {
        self.table = prepared;
    }

    /// STAR's `binarySearch2` (`ReadAlign_stitchAlignToTranscript.cpp` →
    /// `binarySearch2.h`): the index of the junction whose stored donor and
    /// acceptor are exactly `x` and `y`, or `None`.
    ///
    /// Coordinates are genome-absolute, so no chromosome index is needed:
    /// STAR's `sjdbStart` / `sjdbEnd` are already unique across contigs.
    pub fn find(&self, x: u64, y: u64) -> Option<usize> {
        let n = self.table.len();
        if n == 0 || x > self.table[n - 1].stored_start() || x < self.table[0].stored_start() {
            return None;
        }
        let (mut i1, mut i2) = (0usize, n - 1);
        while i2 > i1 + 1 {
            let i3 = usize::midpoint(i1, i2);
            if self.table[i3].stored_start() > x {
                i2 = i3;
            } else {
                i1 = i3;
            }
        }
        let i3 = if x == self.table[i1].stored_start() {
            i1
        } else if x == self.table[i2].stored_start() {
            i2
        } else {
            return None;
        };
        // Scan the run of equal `stored_start` values (backward then forward)
        // for a matching `stored_end`.
        for jj in (0..=i3).rev() {
            if x != self.table[jj].stored_start() {
                break;
            } else if y == self.table[jj].stored_end() {
                return Some(jj);
            }
        }
        for jj in i3..n {
            if x != self.table[jj].stored_start() {
                return None;
            } else if y == self.table[jj].stored_end() {
                return Some(jj);
            }
        }
        None
    }

    /// The junction at `sjA` index `i`, or `None` when the table is absent
    /// or the index is out of range.
    pub fn entry(&self, i: usize) -> Option<&PreparedJunction> {
        self.table.get(i)
    }

    /// Number of entries in the `sjA`-addressable table.
    pub fn table_len(&self) -> usize {
        self.table.len()
    }

    /// Build junction database from GTF file with configurable GTF attribute names.
    pub fn from_gtf_configured(
        gtf_path: &Path,
        genome: &Genome,
        feature_exon: &str,
        chr_prefix: &str,
        transcript_tag: &str,
    ) -> Result<Self, Error> {
        log::info!("Loading GTF annotations from: {}", gtf_path.display());

        let exons = gtf::parse_gtf_configured(gtf_path, feature_exon, chr_prefix)?;
        log::debug!("Parsed {} exon features from GTF", exons.len());

        let raw = gtf::extract_junctions_configured(exons, genome, transcript_tag)?;
        log::info!("Extracted {} annotated junctions from GTF", raw.len());

        // Keep the historical (raw-coordinate) annotated lookup map, but also
        // derive the motif/shift/strand table so the annotated-junction fast
        // paths have something to address. Building the table here makes the
        // align-time GTF path carry the same metadata the index-loaded path
        // already had.
        let mut db = Self::from_raw_junctions(&raw);
        db.set_table(Self::prepare_table(&raw, genome));
        Ok(db)
    }

    /// Run every raw `(chr_idx, intron_start, intron_end, strand)` through
    /// `sjdbPrepare`'s motif detection and micro-repeat shift computation,
    /// then apply STAR's post-dedup sort. Produces exactly the array
    /// `genomeGenerate` writes to `sjdbInfo.txt`.
    fn prepare_table(raw: &[(usize, u64, u64, u8)], genome: &Genome) -> Vec<PreparedJunction> {
        let n_genome_real = genome.n_genome_real;
        let prepared: Vec<PreparedJunction> = raw
            .iter()
            .map(|&(chr_idx, intron_start, intron_end, strand)| {
                sjdb_insert::prepare_junction(
                    chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                    genome,
                    n_genome_real,
                )
            })
            .collect();
        sjdb_insert::sort_and_dedup(prepared)
    }

    /// Build junction database from GTF file (default STAR attribute names).
    pub fn from_gtf(gtf_path: &Path, genome: &Genome) -> Result<Self, Error> {
        Self::from_gtf_configured(gtf_path, genome, "exon", "", "transcript_id")
    }

    /// Build junction database from a pre-extracted list of annotated
    /// junctions `(chr_idx, intron_start, intron_end, strand)`. Used by
    /// the `genomeGenerate` path so it can share the parsed GTF with
    /// `TranscriptomeIndex` and the `sjdb_insert` pipeline without
    /// re-parsing the file.
    ///
    /// The resulting database has **no** `sjA` table; callers that need
    /// [`find`](Self::find) must follow up with [`set_table`](Self::set_table)
    /// or use [`from_prepared`](Self::from_prepared) instead.
    pub fn from_raw_junctions(raw: &[(usize, u64, u64, u8)]) -> Self {
        let mut junctions = HashMap::with_capacity(raw.len());
        for &(chr_idx, intron_start, intron_end, strand) in raw {
            let key = JunctionKey {
                chr_idx,
                intron_start,
                intron_end,
                strand,
            };
            junctions.insert(key, JunctionInfo { annotated: true });
        }
        Self {
            junctions,
            table: Vec::new(),
        }
    }

    /// Check if a junction is annotated in the GTF.
    ///
    /// # Arguments
    /// * `chr_idx` - Chromosome index
    /// * `start` - Genome-absolute 0-based position of the first intronic base
    /// * `end` - Genome-absolute 0-based position of the last intronic base
    /// * `strand` - Strand (0=unknown, 1=+, 2=-)
    ///
    /// # Returns
    /// `true` if junction is annotated, `false` otherwise
    pub fn is_annotated(&self, chr_idx: usize, start: u64, end: u64, strand: u8) -> bool {
        let key = JunctionKey {
            chr_idx,
            intron_start: start,
            intron_end: end,
            strand,
        };
        self.junctions.get(&key).is_some_and(|info| info.annotated)
    }

    /// Get the number of annotated junctions in the database
    pub fn len(&self) -> usize {
        self.junctions.len()
    }

    /// Check if the database is empty
    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }

    /// Insert novel junctions discovered during two-pass mode
    ///
    /// # Arguments
    /// * `novel_junctions` - Vector of (key, info) pairs for novel junctions
    pub fn insert_novel(&mut self, novel_junctions: Vec<(NovelJunctionKey, JunctionInfo)>) {
        for (key, info) in novel_junctions {
            let junction_key = JunctionKey {
                chr_idx: key.chr_idx,
                intron_start: key.intron_start,
                intron_end: key.intron_end,
                strand: key.strand,
            };
            self.junctions.insert(junction_key, info);
        }
    }
}

/// Filter novel junctions by coverage and overhang thresholds (for two-pass mode)
///
/// # Arguments
/// * `sj_stats` - Junction statistics from pass 1
/// * `params` - Parameters (for thresholds)
///
/// # Returns
/// Vector of novel junctions that meet filtering criteria
pub fn filter_novel_junctions(
    sj_stats: &SpliceJunctionStats,
    params: &Parameters,
) -> Vec<(NovelJunctionKey, JunctionInfo)> {
    use crate::align::score::SpliceMotif;
    use std::sync::atomic::Ordering;

    let max_intron = if params.align_intron_max == 0 {
        params.win_bin_window_dist()
    } else {
        params.align_intron_max as u64
    };

    sj_stats
        .iter()
        .filter_map(|entry| {
            let key = entry.key();
            let counts = entry.value();

            // Skip if already annotated (from GTF)
            if counts.annotated {
                return None;
            }

            let unique = counts.unique_count.load(Ordering::Relaxed);
            let multi = counts.multi_count.load(Ordering::Relaxed);
            let max_overhang = counts.max_overhang.load(Ordering::Relaxed);

            // Use motif-specific thresholds from outSJfilter* params
            let cat = SpliceMotif::filter_category_from_encoded(key.motif);

            // Overhang threshold (motif-specific)
            let min_overhang = params.out_sj_filter_overhang_min[cat] as u32;
            let has_overhang = max_overhang >= min_overhang;

            // Coverage threshold (motif-specific). STAR keeps a junction if EITHER
            // the unique-read count OR the total (unique+multi) count meets its
            // threshold — an OR, not an AND (STAR manual: "Junctions are output if
            // one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
            // conditions are satisfied"; confirmed against STAR's source and the
            // byte-faithful STAR-rs `star_sj.rs`). Using AND here dropped every
            // junction supported only by multi-mapping reads (unique==0), which is
            // why rustar-aligner reported far fewer novel junctions than STAR.
            let min_unique = params.out_sj_filter_count_unique_min[cat] as u32;
            let min_total = params.out_sj_filter_count_total_min[cat] as u32;
            let total = unique + multi;
            let has_coverage = unique >= min_unique || total >= min_total;

            // Intron length threshold
            let intron_len = key.intron_end.saturating_sub(key.intron_start) + 1;
            let within_intron_limit = intron_len <= max_intron;

            if has_coverage && has_overhang && within_intron_limit {
                let novel_key = NovelJunctionKey {
                    chr_idx: key.chr_idx,
                    intron_start: key.intron_start,
                    intron_end: key.intron_end,
                    strand: key.strand,
                };
                let info = JunctionInfo {
                    annotated: false, // Novel junctions are not annotated
                };
                Some((novel_key, info))
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_junction_db_empty() {
        let db = SpliceJunctionDb::empty();
        assert_eq!(db.len(), 0);
        assert!(db.is_empty());
        assert!(!db.is_annotated(0, 100, 200, 1));
    }

    #[test]
    fn test_junction_key_equality() {
        let key1 = JunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
        };
        let key2 = JunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
        };
        let key3 = JunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 2,
        };

        assert_eq!(key1, key2);
        assert_ne!(key1, key3); // Different strand
    }

    #[test]
    fn test_junction_lookup() {
        let mut db = SpliceJunctionDb::empty();

        // Manually insert a junction
        db.junctions.insert(
            JunctionKey {
                chr_idx: 0,
                intron_start: 100,
                intron_end: 200,
                strand: 1,
            },
            JunctionInfo { annotated: true },
        );

        // Should find annotated junction
        assert!(db.is_annotated(0, 100, 200, 1));

        // Should not find with different strand
        assert!(!db.is_annotated(0, 100, 200, 2));

        // Should not find with different coordinates
        assert!(!db.is_annotated(0, 101, 200, 1));
        assert!(!db.is_annotated(0, 100, 201, 1));
    }

    #[test]
    fn test_junction_strand_specific() {
        let mut db = SpliceJunctionDb::empty();

        // Add same junction coordinates but different strands
        db.junctions.insert(
            JunctionKey {
                chr_idx: 0,
                intron_start: 100,
                intron_end: 200,
                strand: 1,
            },
            JunctionInfo { annotated: true },
        );
        db.junctions.insert(
            JunctionKey {
                chr_idx: 0,
                intron_start: 100,
                intron_end: 200,
                strand: 2,
            },
            JunctionInfo { annotated: true },
        );

        assert_eq!(db.len(), 2);
        assert!(db.is_annotated(0, 100, 200, 1));
        assert!(db.is_annotated(0, 100, 200, 2));
        assert!(!db.is_annotated(0, 100, 200, 0)); // Unknown strand
    }

    #[test]
    fn test_insert_novel_junctions() {
        let mut db = SpliceJunctionDb::empty();

        // Insert a novel junction
        let key = NovelJunctionKey {
            chr_idx: 0,
            intron_start: 100,
            intron_end: 200,
            strand: 1,
        };
        let info = JunctionInfo { annotated: false };
        db.insert_novel(vec![(key, info)]);

        assert_eq!(db.len(), 1);
        assert!(!db.is_annotated(0, 100, 200, 1)); // Novel, not annotated

        // Insert another novel junction
        let key2 = NovelJunctionKey {
            chr_idx: 0,
            intron_start: 300,
            intron_end: 400,
            strand: 2,
        };
        let info2 = JunctionInfo { annotated: false };
        db.insert_novel(vec![(key2, info2)]);

        assert_eq!(db.len(), 2);
    }

    #[test]
    fn test_filter_novel_junctions() {
        use crate::align::score::SpliceMotif;

        let sj_stats = SpliceJunctionStats::new();

        // Add a high-quality novel canonical junction (should pass filter)
        // Needs overhang >= 12 (default outSJfilterOverhangMin for GT/AG)
        // Needs unique >= 1 (default outSJfilterCountUniqueMin for GT/AG)
        sj_stats.record_junction(0, 100, 200, 1, SpliceMotif::GtAg, true, 20, false);

        // Add a low-overhang novel junction (should fail filter: overhang 2 < 12)
        sj_stats.record_junction(0, 300, 400, 1, SpliceMotif::GtAg, true, 2, false);

        // Add an annotated junction (should be excluded from novel list)
        sj_stats.record_junction(0, 500, 600, 1, SpliceMotif::GtAg, true, 20, true);

        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "reads.fq"]);
        let novel_junctions = filter_novel_junctions(&sj_stats, &params);

        // Should only get the high-quality novel junction
        assert_eq!(novel_junctions.len(), 1);
        assert_eq!(novel_junctions[0].0.intron_start, 100);
        assert_eq!(novel_junctions[0].0.intron_end, 200);
        assert!(!novel_junctions[0].1.annotated);
    }

    #[test]
    fn test_filter_novel_junctions_noncanonical_strict() {
        use crate::align::score::SpliceMotif;

        let sj_stats = SpliceJunctionStats::new();

        // Non-canonical junction with moderate overhang (20 < 30 default for non-canonical)
        // Record 5 unique reads (>= 3 count threshold)
        for _ in 0..5 {
            sj_stats.record_junction(0, 100, 200, 1, SpliceMotif::NonCanonical, true, 20, false);
        }

        // Non-canonical junction with enough overhang (35 >= 30)
        for _ in 0..5 {
            sj_stats.record_junction(0, 300, 400, 1, SpliceMotif::NonCanonical, true, 35, false);
        }

        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "reads.fq"]);
        let novel_junctions = filter_novel_junctions(&sj_stats, &params);

        // Only the 35-overhang junction should pass (30bp minimum for non-canonical)
        assert_eq!(novel_junctions.len(), 1);
        assert_eq!(novel_junctions[0].0.intron_start, 300);
    }

    /// Build a canonical (motif != 0) prepared junction whose stored
    /// coordinates are exactly `(start, end)`.
    fn pj(start: u64, end: u64, strand: u8) -> PreparedJunction {
        PreparedJunction {
            chr_idx: 0,
            start_pos: start,
            end_pos: end,
            motif: 1,
            shift_left: 0,
            shift_right: 0,
            strand,
        }
    }

    #[test]
    fn find_round_trips_every_prepared_junction() {
        // Deliberately unsorted input: `sort_and_dedup` establishes the
        // (stored_start, stored_end) order `find`'s binary search needs.
        let prepared = sjdb_insert::sort_and_dedup(vec![
            pj(900, 1000, 1),
            pj(100, 200, 1),
            pj(500, 600, 2),
            pj(300, 400, 1),
        ]);
        let db = SpliceJunctionDb::from_prepared(prepared.clone());
        assert_eq!(db.table_len(), 4);

        // The invariant the sjAB fast path and the annotated snap both rely
        // on: index i in the table is addressable by its own stored coords.
        for (i, j) in prepared.iter().enumerate() {
            assert_eq!(
                db.find(j.stored_start(), j.stored_end()),
                Some(i),
                "junction {i} did not round-trip"
            );
            assert_eq!(db.entry(i).unwrap().stored_start(), j.stored_start());
        }

        // Misses in every direction: below the first, above the last, a
        // start that exists with a wrong end, and an end that exists with a
        // wrong start.
        assert_eq!(db.find(50, 200), None);
        assert_eq!(db.find(2000, 3000), None);
        assert_eq!(db.find(100, 201), None);
        assert_eq!(db.find(101, 200), None);
    }

    #[test]
    fn find_disambiguates_a_run_of_equal_starts() {
        // Several junctions sharing a donor: the backward/forward scan around
        // the binary-search landing point must pick the right acceptor.
        let prepared = sjdb_insert::sort_and_dedup(vec![
            pj(100, 700, 1),
            pj(100, 200, 1),
            pj(100, 500, 1),
            pj(100, 300, 1),
            pj(900, 1000, 1),
        ]);
        let db = SpliceJunctionDb::from_prepared(prepared.clone());
        for (i, j) in prepared.iter().enumerate() {
            assert_eq!(db.find(j.stored_start(), j.stored_end()), Some(i));
        }
        assert_eq!(db.find(100, 400), None);
    }

    #[test]
    fn find_is_inert_without_a_table() {
        // `from_raw_junctions` carries no motif/shift metadata, so the
        // annotated fast paths must simply not engage rather than misfire.
        let db = SpliceJunctionDb::from_raw_junctions(&[(0, 100, 200, 1)]);
        assert!(db.is_annotated(0, 100, 200, 1));
        assert_eq!(db.table_len(), 0);
        assert_eq!(db.find(100, 200), None);
        assert!(db.entry(0).is_none());
    }

    #[test]
    fn from_prepared_keys_the_map_on_stored_coordinates() {
        // Non-canonical junction: stored coords are the shifted ones, which is
        // what the stitch-time scan produces.
        let noncan = PreparedJunction {
            chr_idx: 0,
            start_pos: 100,
            end_pos: 200,
            motif: 0,
            shift_left: 3,
            shift_right: 0,
            strand: 0,
        };
        assert_eq!(noncan.stored_start(), 100);
        assert_eq!(noncan.original_start(), 103);

        let db = SpliceJunctionDb::from_prepared(vec![noncan]);
        assert!(db.is_annotated(0, 100, 200, 0));
        assert!(!db.is_annotated(0, 103, 203, 0));
        assert_eq!(db.find(100, 200), Some(0));
    }

    #[test]
    fn test_db_keyed_in_genome_absolute_zero_based_multi_chr() {
        use crate::junction::gtf::{GtfRecord, extract_junctions_configured};
        use std::collections::HashMap;

        // Two-chromosome toy genome so chr_start[1] != 0.
        let genome = Genome {
            transform_blocks: None,
            sequence: vec![0; 4000].into(),
            n_genome: 2000,
            n_genome_real: 2000,
            n_chr_real: 2,
            chr_start: vec![0, 1000, 2000],
            chr_length: vec![1000, 1000],
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        };

        let make_exon = |seqname: &str, start: u64, end: u64, transcript: &str| -> GtfRecord {
            let mut attrs: HashMap<String, String> = HashMap::new();
            attrs.insert("gene_id".to_string(), "G".to_string());
            attrs.insert("transcript_id".to_string(), transcript.to_string());
            GtfRecord {
                seqname: seqname.to_string(),
                feature: "exon".to_string(),
                start,
                end,
                strand: '+',
                attributes: attrs,
            }
        };

        let exons = vec![
            make_exon("chr1", 100, 200, "T1"),
            make_exon("chr1", 300, 400, "T1"),
            make_exon("chr2", 100, 200, "T2"),
            make_exon("chr2", 300, 400, "T2"),
        ];

        let raw = extract_junctions_configured(exons, &genome, "transcript_id").unwrap();
        let db = SpliceJunctionDb::from_raw_junctions(&raw);
        assert_eq!(db.len(), 2);

        // Junction on chr1: intron local 1-based 201..299
        // → genome-absolute 0-based: chr_start[0] + 200 .. chr_start[0] + 298
        assert!(db.is_annotated(0, 200, 298, 1));
        // Off-by-one in either direction must miss.
        assert!(!db.is_annotated(0, 201, 299, 1));
        assert!(!db.is_annotated(0, 199, 297, 1));

        // Junction on chr2: same chr-local coords, but chr_start[1] = 1000
        // → genome-absolute 0-based: 1200 .. 1298
        assert!(db.is_annotated(1, 1200, 1298, 1));
        // The pre-fix chr-local 1-based key (201, 299) must not match on chr2.
        assert!(!db.is_annotated(1, 201, 299, 1));
        // The pre-fix stitch-time off-by-one (1199, 1297) must not match either.
        assert!(!db.is_annotated(1, 1199, 1297, 1));
    }
}
