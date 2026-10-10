//! `--sjdbFileChrStartEnd`: explicit splice junctions given as TSV files, unioned with any
//! `--sjdbGTFfile`-derived junctions before sjdb insertion (STAR's `sjdbFileChrStartEnd`).
//!
//! Each line is `chr  start  end  [strand]`, whitespace-delimited, with `start`/`end` the 1-based
//! coordinates of the intron's first and last base. The strand column is optional; `.` or absent
//! means unknown (the motif is derived from the genome at sjdb-insertion time). Lines with fewer
//! than 3 fields, or a non-integer `start`/`end`, are skipped (STAR's `istringstream >>` silently
//! stops parsing a malformed line). An unresolvable chromosome name is a hard error, matching
//! STAR's fatal exit on the same condition.

use std::path::Path;

use crate::error::Error;
use crate::genome::Genome;

/// Parse `--sjdbFileChrStartEnd` TSV files into raw junction tuples
/// `(chr_idx, intron_start, intron_end, strand)`, in the same genome-absolute 0-based /
/// `0=unknown,1=+,2=-` convention as [`crate::junction::gtf::extract_junctions_configured`].
pub fn parse_sjdb_chr_start_end(
    paths: &[std::path::PathBuf],
    genome: &Genome,
) -> Result<Vec<(usize, u64, u64, u8)>, Error> {
    let mut out = Vec::new();
    for path in paths {
        let text = std::fs::read_to_string(path).map_err(|e| Error::io(e, path.clone()))?;
        for line in text.lines() {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 3 {
                continue;
            }
            let (Ok(start), Ok(end)) = (fields[1].parse::<u64>(), fields[2].parse::<u64>()) else {
                continue;
            };
            let strand = match fields.get(3).and_then(|s| s.chars().next()) {
                // STAR's sjdbLoadFromStream also takes SJ.out.tab's numeric codes.
                Some('+' | '1') => 1u8,
                Some('-' | '2') => 2u8,
                _ => 0u8,
            };
            let chr_idx = chr_index(genome, fields[0], path)?;
            let base = genome.chr_start[chr_idx];
            out.push((chr_idx, base + start - 1, base + end - 1, strand));
        }
    }
    Ok(out)
}

fn chr_index(genome: &Genome, name: &str, path: &Path) -> Result<usize, Error> {
    genome
        .chr_name
        .iter()
        .position(|n| n == name)
        .ok_or_else(|| {
            Error::Parameter(format!(
                "{}: chromosome '{name}' not found in the genome (--sjdbFileChrStartEnd)",
                path.display()
            ))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_genome() -> Genome {
        Genome {
            transform_blocks: None,
            sequence: vec![0u8; 2000].into(),
            n_genome: 2000,
            n_genome_real: 2000,
            n_chr_real: 2,
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
            chr_length: vec![1000, 1000],
            chr_start: vec![0, 1000, 2000],
        }
    }

    #[test]
    fn parses_chr_start_end_strand() {
        let genome = tiny_genome();
        let dir =
            std::env::temp_dir().join(format!("sjdb-chr-start-end-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("junctions.tsv");
        std::fs::write(&path, "chr1\t101\t200\t+\nchr2\t51\t150\t-\n").unwrap();

        let raw = parse_sjdb_chr_start_end(std::slice::from_ref(&path), &genome).unwrap();
        assert_eq!(raw, vec![(0, 100, 199, 1), (1, 1050, 1149, 2)]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn parses_sj_out_tab_numeric_strands() {
        // A 1st-pass SJ.out.tab, as STAR's sjdbLoadFromStream reads it.
        let genome = tiny_genome();
        let dir = std::env::temp_dir().join(format!("sjdb-sjout-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("SJ.out.tab");
        std::fs::write(&path, "chr1\t101\t200\t1\t1\t0\t3\t0\t30\nchr2\t51\t150\t2\t2\t0\t1\t0\t12\nchr2\t300\t400\t0\t0\t0\t1\t0\t9\n").unwrap();
        let raw = parse_sjdb_chr_start_end(std::slice::from_ref(&path), &genome).unwrap();
        assert_eq!(
            raw,
            vec![(0, 100, 199, 1), (1, 1050, 1149, 2), (1, 1299, 1399, 0)]
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn missing_strand_defaults_to_unknown() {
        let genome = tiny_genome();
        let dir =
            std::env::temp_dir().join(format!("sjdb-chr-start-end-test2-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("junctions.tsv");
        std::fs::write(&path, "chr1\t101\t200\n").unwrap();

        let raw = parse_sjdb_chr_start_end(std::slice::from_ref(&path), &genome).unwrap();
        assert_eq!(raw, vec![(0, 100, 199, 0)]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn short_lines_are_skipped() {
        let genome = tiny_genome();
        let dir =
            std::env::temp_dir().join(format!("sjdb-chr-start-end-test3-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("junctions.tsv");
        std::fs::write(&path, "chr1\t101\n# comment too short\nchr1\t101\t200\t+\n").unwrap();

        let raw = parse_sjdb_chr_start_end(std::slice::from_ref(&path), &genome).unwrap();
        assert_eq!(raw, vec![(0, 100, 199, 1)]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn unknown_chromosome_errors() {
        let genome = tiny_genome();
        let dir =
            std::env::temp_dir().join(format!("sjdb-chr-start-end-test4-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("junctions.tsv");
        std::fs::write(&path, "chrZ\t101\t200\t+\n").unwrap();

        let err = parse_sjdb_chr_start_end(std::slice::from_ref(&path), &genome).unwrap_err();
        assert!(err.to_string().contains("chrZ"));
        std::fs::remove_dir_all(&dir).ok();
    }
}
