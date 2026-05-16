/// MAPQ (mapping quality) calculation
///
/// Uses STAR's exact lookup table (ReadAlign_outputTranscriptSAM.cpp):
/// - Unique mappers (n=1): use outSAMmapqUnique (default 255)
/// - n=2: MAPQ 3
/// - n=3 or 4: MAPQ 1
/// - n>=5: MAPQ 0
/// - Unmapped (n=0): 0
///
/// # Arguments
/// * `n_alignments` - Number of valid alignments for this read
/// * `mapq_unique` - MAPQ value for unique mappers (typically 255)
///
/// # Returns
/// MAPQ score (0-255)
pub fn calculate_mapq(n_alignments: usize, mapq_unique: u8) -> u8 {
    match n_alignments {
        1 => mapq_unique,
        2 => 3,
        3 | 4 => 1,
        _ => 0, // 0 or >= 5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapq_unmapped() {
        assert_eq!(calculate_mapq(0, 255), 0);
    }

    #[test]
    fn test_mapq_unique() {
        assert_eq!(calculate_mapq(1, 255), 255);
        assert_eq!(calculate_mapq(1, 60), 60);
    }

    #[test]
    fn test_mapq_multi() {
        // STAR lookup table values
        assert_eq!(calculate_mapq(2, 255), 3);
        assert_eq!(calculate_mapq(10, 255), 0);
        assert_eq!(calculate_mapq(100, 255), 0);
    }

    #[test]
    fn test_mapq_capped() {
        // Even with very high alignment quality, cap at 255
        assert_eq!(calculate_mapq(1, 255), 255);
    }

    #[test]
    fn test_mapq_star_lookup() {
        // Verify STAR's exact lookup table values
        assert_eq!(calculate_mapq(2, 255), 3);
        assert_eq!(calculate_mapq(3, 255), 1);
        assert_eq!(calculate_mapq(4, 255), 1);
        assert_eq!(calculate_mapq(5, 255), 0);
        assert_eq!(calculate_mapq(6, 255), 0);
    }
}
