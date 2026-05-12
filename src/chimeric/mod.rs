// Chimeric alignment detection module
//
// Detects split reads that span distant genomic locations:
// - Inter-chromosomal fusions (e.g., BCR-ABL)
// - Intra-chromosomal strand breaks
// - Circular RNAs (back-splices)
//
// Detection strategy:
// - Tier 1: chimericDetectionOld (post-stitching, transcript-pair search — detect_chimeric_old)
// - Tier 2: Multi-cluster chimeric stitching (during clustering)
// - Tier 1b: detect_from_soft_clips (re-seed primary soft-clips when Tier 1 finds nothing)
// - Tier 3: detect_from_chimeric_residuals (re-seed outer uncovered regions of Tier 1/2 pairs)

mod detect;
mod output;
mod score;
mod segment;

pub use detect::{ChimericDetector, detect_chimeric_old, detect_inter_mate_chimeric};
pub use output::{ChimericJunctionWriter, build_within_bam_records};
pub use segment::{ChimericAlignment, ChimericSegment};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_module_exports() {
        // Ensure all public types are accessible
        let _ = std::mem::size_of::<ChimericAlignment>();
        let _ = std::mem::size_of::<ChimericSegment>();
    }
}
