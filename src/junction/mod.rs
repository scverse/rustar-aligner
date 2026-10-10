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

/// Splice junction database built from GTF annotations
#[derive(Clone)]
pub struct SpliceJunctionDb {
    /// Map: (chr_idx, intron_start, intron_end, strand) → annotated
    junctions: HashMap<JunctionKey, JunctionInfo>,
}

impl SpliceJunctionDb {
    /// Create empty database (for no-GTF mode)
    pub fn empty() -> Self {
        Self {
            junctions: HashMap::new(),
        }
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

        Ok(Self::from_raw_junctions(&raw))
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
        Self { junctions }
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
