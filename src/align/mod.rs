pub mod pe_overlap;
pub mod read_align;
pub mod score;
pub mod seed;
// Private to the crate, except under the `bench` feature: `benches/hot_paths.rs`
// measures `find_stop` directly, because it is the function a portable-SIMD
// crate would replace (#205) and the one whose cost is worth watching.
#[cfg(not(feature = "bench"))]
mod simd_scan;
#[cfg(feature = "bench")]
#[doc(hidden)]
pub mod simd_scan;
pub mod stitch;
pub mod transcript;

// Re-export commonly used types
pub use read_align::{
    AlignReadResult, PairedAlignment, PairedAlignmentResult, align_paired_read, align_read,
};
pub use seed::Seed;
pub use stitch::{SeedCluster, WindowAlignment, stitch_seeds};
pub use transcript::{Exon, Transcript};
