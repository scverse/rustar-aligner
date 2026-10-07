// Chimeric alignment detection module
//
// Detects split reads that span distant genomic locations:
// - Inter-chromosomal fusions (e.g., BCR-ABL)
// - Intra-chromosomal strand breaks
// - Circular RNAs (back-splices)
//
// Detection (`detect.rs`) is STAR's: `chimericDetectionOld` by default, or
// `chimericDetectionMult` under `--chimMultimapNmax`, both over the read's
// window transcripts — combined two-mate transcripts for a pair.

mod detect;
mod output;
mod segment;

pub use detect::{ChimRead, WinTr, chimeric_detection};
pub use output::{ChimericJunctionWriter, build_within_bam_records};
pub use segment::{ChimericAlignment, ChimericSegment, ExonSpan, JunctionLine, MultimapInfo};

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
