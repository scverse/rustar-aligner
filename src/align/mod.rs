pub mod cr_primary;
pub mod pe_overlap;
pub mod read_align;
pub mod score;
pub mod seed;
mod simd_scan;
pub mod stitch;
pub mod transcript;

// Re-export commonly used types
pub use read_align::{
    AlignReadResult, PairedAlignment, PairedAlignmentResult, align_paired_read, align_read,
    last_unmapped_best,
};
pub use seed::Seed;
pub use stitch::{SeedCluster, WindowAlignment, stitch_seeds};
pub use transcript::{Exon, Transcript};
