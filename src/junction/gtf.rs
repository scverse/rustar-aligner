/// GTF file parsing for gene annotations
///
/// Supports standard GTF format (tab-separated, 9 columns):
/// 1. seqname (chromosome)
/// 2. source (ignored)
/// 3. feature (gene, transcript, exon, etc.)
/// 4. start (1-based inclusive)
/// 5. end (1-based inclusive)
/// 6. score (ignored)
/// 7. strand (+, -, .)
/// 8. frame (ignored)
/// 9. attributes (semicolon-separated key-value pairs)
use crate::error::Error;
use crate::genome::Genome;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// GTF record (single line)
#[derive(Debug, Clone)]
pub struct GtfRecord {
    pub seqname: String,
    pub feature: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub attributes: HashMap<String, String>,
}

/// Parse GTF file, returning only records matching `feature_exon`.
///
/// `chr_prefix` is prepended to every seqname (STAR: `sjdbGTFchrPrefix`).
/// `feature_exon` is the feature column value to keep (STAR: `sjdbGTFfeatureExon`, default `"exon"`).
pub fn parse_gtf_configured(
    path: &Path,
    feature_exon: &str,
    chr_prefix: &str,
) -> Result<Vec<GtfRecord>, Error> {
    let file =
        File::open(path).map_err(|e| Error::Gtf(format!("Failed to open GTF file: {}", e)))?;
    let reader = BufReader::new(file);

    let mut exons = Vec::new();
    let mut line_num = 0;

    for line in reader.lines() {
        line_num += 1;
        let line =
            line.map_err(|e| Error::Gtf(format!("Failed to read line {}: {}", line_num, e)))?;

        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        match parse_gtf_line(line) {
            Ok(mut record) => {
                if record.feature.eq_ignore_ascii_case(feature_exon) {
                    if !chr_prefix.is_empty() {
                        record.seqname = format!("{}{}", chr_prefix, record.seqname);
                    }
                    exons.push(record);
                }
            }
            Err(e) => {
                log::warn!("Skipping malformed GTF line {}: {}", line_num, e);
            }
        }
    }

    Ok(exons)
}

/// Parse GTF file and extract exon features (default: feature `"exon"`, no chr prefix).
#[allow(dead_code)]
pub fn parse_gtf(path: &Path) -> Result<Vec<GtfRecord>, Error> {
    parse_gtf_configured(path, "exon", "")
}

/// Parse a single GTF line
fn parse_gtf_line(line: &str) -> Result<GtfRecord, Error> {
    let fields: Vec<&str> = line.split('\t').collect();

    if fields.len() < 9 {
        return Err(Error::Gtf(format!(
            "GTF line has {} fields, expected 9",
            fields.len()
        )));
    }

    let seqname = fields[0].to_string();
    let feature = fields[2].to_string();
    let start = fields[3]
        .parse::<u64>()
        .map_err(|e| Error::Gtf(format!("Invalid start position: {}", e)))?;
    let end = fields[4]
        .parse::<u64>()
        .map_err(|e| Error::Gtf(format!("Invalid end position: {}", e)))?;
    let strand = fields[6]
        .chars()
        .next()
        .ok_or_else(|| Error::Gtf("Empty strand field".to_string()))?;

    // Parse attributes (semicolon-separated key-value pairs)
    let attributes = parse_attributes(fields[8])?;

    Ok(GtfRecord {
        seqname,
        feature,
        start,
        end,
        strand,
        attributes,
    })
}

/// Parse GTF attributes field
///
/// Format: key1 "value1"; key2 "value2";
fn parse_attributes(attr_str: &str) -> Result<HashMap<String, String>, Error> {
    let mut attributes = HashMap::new();

    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }

        // Split on first space to separate key and value
        let parts: Vec<&str> = pair.splitn(2, ' ').collect();
        if parts.len() != 2 {
            continue; // Skip malformed attributes
        }

        let key = parts[0].trim().to_string();
        let value = parts[1]
            .trim()
            .trim_matches('"') // Remove quotes
            .to_string();

        attributes.insert(key, value);
    }

    Ok(attributes)
}

/// Extract junctions from exon records, grouping by `transcript_tag`.
///
/// `transcript_tag` is the GTF attribute key for the parent transcript
/// (STAR: `sjdbGTFtagExonParentTranscript`, default `"transcript_id"`).
///
/// Returns: Vec<(chr_idx, intron_start, intron_end, strand)>
pub fn extract_junctions_configured(
    exons: Vec<GtfRecord>,
    genome: &Genome,
    transcript_tag: &str,
) -> Result<Vec<(usize, u64, u64, u8)>, Error> {
    let mut transcripts: HashMap<String, Vec<GtfRecord>> = HashMap::new();

    for exon in exons {
        let transcript_id = exon
            .attributes
            .get(transcript_tag)
            .ok_or_else(|| Error::Gtf(format!("Exon missing {} attribute", transcript_tag)))?
            .clone();

        transcripts.entry(transcript_id).or_default().push(exon);
    }

    // Extract junctions from each transcript
    let mut junctions = Vec::new();

    for (_transcript_id, mut exons) in transcripts {
        if exons.len() < 2 {
            // Single-exon transcript, no junctions
            continue;
        }

        // Get chromosome index
        let chr_name = &exons[0].seqname;
        let chr_idx = genome.chr_name.iter().position(|name| name == chr_name);

        let chr_idx = match chr_idx {
            Some(idx) => idx,
            None => {
                log::warn!("Skipping transcript on unknown chromosome: {}", chr_name);
                continue;
            }
        };

        // Convert strand
        let strand = match exons[0].strand {
            '+' => 1u8,
            '-' => 2u8,
            _ => 0u8, // Unknown strand
        };

        // Sort exons by position
        exons.sort_by_key(|e| e.start);

        // Calculate junction coordinates from consecutive exons
        for i in 0..exons.len() - 1 {
            let exon1 = &exons[i];
            let exon2 = &exons[i + 1];

            // Intron coordinates (1-based, STAR convention)
            let intron_start = exon1.end + 1;
            let intron_end = exon2.start - 1;

            // Validate junction
            if intron_end <= intron_start {
                log::warn!(
                    "Invalid junction coordinates: {}-{} (possibly overlapping exons)",
                    intron_start,
                    intron_end
                );
                continue;
            }

            junctions.push((chr_idx, intron_start, intron_end, strand));
        }
    }

    // Deduplicate junctions (same junction can appear in multiple transcripts)
    junctions.sort_unstable();
    junctions.dedup();

    Ok(junctions)
}

/// Extract junctions grouping by `"transcript_id"` (default, backward-compatible wrapper).
#[allow(dead_code)]
pub fn extract_junctions_from_exons(
    exons: Vec<GtfRecord>,
    genome: &Genome,
) -> Result<Vec<(usize, u64, u64, u8)>, Error> {
    extract_junctions_configured(exons, genome, "transcript_id")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_attributes() {
        let attr = r#"gene_id "ENSG001"; transcript_id "ENST001"; gene_name "MYC";"#;
        let attrs = parse_attributes(attr).unwrap();

        assert_eq!(attrs.get("gene_id"), Some(&"ENSG001".to_string()));
        assert_eq!(attrs.get("transcript_id"), Some(&"ENST001".to_string()));
        assert_eq!(attrs.get("gene_name"), Some(&"MYC".to_string()));
    }

    #[test]
    fn test_parse_attributes_no_trailing_semicolon() {
        let attr = r#"gene_id "ENSG001"; transcript_id "ENST001""#;
        let attrs = parse_attributes(attr).unwrap();

        assert_eq!(attrs.get("gene_id"), Some(&"ENSG001".to_string()));
        assert_eq!(attrs.get("transcript_id"), Some(&"ENST001".to_string()));
    }

    #[test]
    fn test_parse_gtf_line_valid() {
        let line = "chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";";
        let record = parse_gtf_line(line).unwrap();

        assert_eq!(record.seqname, "chr1");
        assert_eq!(record.feature, "exon");
        assert_eq!(record.start, 100);
        assert_eq!(record.end, 200);
        assert_eq!(record.strand, '+');
        assert_eq!(record.attributes.get("gene_id"), Some(&"G1".to_string()));
    }

    #[test]
    fn test_parse_gtf_line_invalid_columns() {
        let line = "chr1\ttest\texon"; // Only 3 columns
        assert!(parse_gtf_line(line).is_err());
    }

    #[test]
    fn test_parse_gtf_with_comments() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "# Comment line\nchr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
        )
        .unwrap();

        let exons = parse_gtf(file.path()).unwrap();
        assert_eq!(exons.len(), 1);
    }

    #[test]
    fn test_parse_gtf_filters_non_exons() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "chr1\ttest\tgene\t50\t300\t.\t+\t.\tgene_id \"G1\";").unwrap();
        writeln!(
            file,
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
        )
        .unwrap();

        let exons = parse_gtf(file.path()).unwrap();
        assert_eq!(exons.len(), 1); // Only exon, not gene
        assert_eq!(exons[0].feature, "exon");
    }

    #[test]
    fn test_extract_junctions_single_transcript() {
        // Create a simple genome
        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        // Create two exons for one transcript
        let exons = vec![
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 100,
                end: 200,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 300,
                end: 400,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        ];

        let junctions = extract_junctions_from_exons(exons, &genome).unwrap();

        assert_eq!(junctions.len(), 1);
        let (chr_idx, start, end, strand) = junctions[0];
        assert_eq!(chr_idx, 0);
        assert_eq!(start, 201); // exon1.end + 1
        assert_eq!(end, 299); // exon2.start - 1
        assert_eq!(strand, 1); // + strand
    }

    #[test]
    fn test_extract_junctions_multiple_transcripts() {
        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        // Two transcripts with the same junction
        let exons = vec![
            // Transcript 1
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 100,
                end: 200,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 300,
                end: 400,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            // Transcript 2 (same junction)
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 100,
                end: 200,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T2".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 300,
                end: 500,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T2".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        ];

        let junctions = extract_junctions_from_exons(exons, &genome).unwrap();

        // Should have 1 unique junction (both transcripts share junction 201-299)
        // Note: T1 has 100-200, 300-400 and T2 has 100-200, 300-500
        // They share the junction from exon ending at 200 to exon starting at 300
        assert_eq!(junctions.len(), 1);
        assert_eq!(junctions[0], (0, 201, 299, 1)); // chr0, junction 201-299, strand +
    }

    #[test]
    fn test_extract_junctions_single_exon_transcript() {
        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        // Single exon transcript (no junctions)
        let exons = vec![GtfRecord {
            seqname: "chr1".to_string(),
            feature: "exon".to_string(),
            start: 100,
            end: 200,
            strand: '+',
            attributes: vec![
                ("gene_id".to_string(), "G1".to_string()),
                ("transcript_id".to_string(), "T1".to_string()),
            ]
            .into_iter()
            .collect(),
        }];

        let junctions = extract_junctions_from_exons(exons, &genome).unwrap();

        assert_eq!(junctions.len(), 0); // No junctions
    }

    #[test]
    fn test_extract_junctions_unknown_chromosome() {
        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        // Exons on unknown chromosome
        let exons = vec![
            GtfRecord {
                seqname: "chr99".to_string(), // Unknown
                feature: "exon".to_string(),
                start: 100,
                end: 200,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            GtfRecord {
                seqname: "chr99".to_string(),
                feature: "exon".to_string(),
                start: 300,
                end: 400,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        ];

        let junctions = extract_junctions_from_exons(exons, &genome).unwrap();

        assert_eq!(junctions.len(), 0); // Skipped due to unknown chr
    }

    #[test]
    fn test_junction_coordinate_calculation() {
        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        // Exon 1: 100-200, Exon 2: 300-400
        // Expected junction: 201-299 (intron coordinates)
        let exons = vec![
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 100,
                end: 200,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 300,
                end: 400,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("transcript_id".to_string(), "T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        ];

        let junctions = extract_junctions_from_exons(exons, &genome).unwrap();

        assert_eq!(junctions.len(), 1);
        let (_chr_idx, start, end, _strand) = junctions[0];
        assert_eq!(start, 201);
        assert_eq!(end, 299);
    }

    #[test]
    fn test_parse_gtf_configured_chr_prefix() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
        )
        .unwrap();

        let exons = parse_gtf_configured(file.path(), "exon", "chr").unwrap();
        assert_eq!(exons.len(), 1);
        assert_eq!(exons[0].seqname, "chr1");
    }

    #[test]
    fn test_parse_gtf_configured_custom_feature() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
        )
        .unwrap();
        writeln!(
            file,
            "chr1\ttest\tCDS\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
        )
        .unwrap();

        let exons = parse_gtf_configured(file.path(), "CDS", "").unwrap();
        assert_eq!(exons.len(), 1);
        assert_eq!(exons[0].feature, "CDS");
    }

    #[test]
    fn test_extract_junctions_configured_custom_transcript_tag() {
        let genome = Genome {
            sequence: vec![0; 1000],
            n_genome: 1000,
            n_chr_real: 1,
            chr_start: vec![0, 1000],
            chr_length: vec![1000],
            chr_name: vec!["chr1".to_string()],
        };

        let exons = vec![
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 100,
                end: 200,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("Parent".to_string(), "transcript:T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
            GtfRecord {
                seqname: "chr1".to_string(),
                feature: "exon".to_string(),
                start: 300,
                end: 400,
                strand: '+',
                attributes: vec![
                    ("gene_id".to_string(), "G1".to_string()),
                    ("Parent".to_string(), "transcript:T1".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        ];

        let junctions = extract_junctions_configured(exons, &genome, "Parent").unwrap();
        assert_eq!(junctions.len(), 1);
        assert_eq!(junctions[0], (0, 201, 299, 1));
    }
}
