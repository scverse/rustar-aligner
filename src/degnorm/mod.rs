//! DegNorm-style transcript degradation normalization.
//!
//! **Not a STAR feature** (see `DIVERGENCE.md`). Two phases:
//!
//! 1. `--quantMode GeneCoverage` during alignment writes per-gene, per-exonic
//!    base coverage plus raw unique counts to `GeneCoverage.out.bin`
//!    ([`crate::quant::coverage`]).
//! 2. `--runMode degNorm` merges several samples' coverage files and fits
//!    DegNorm's rank-one NMF over-approximation, producing Degradation Index
//!    (DI) scores and degradation-adjusted read counts.
//!
//! The DI is defined *across* samples: the model separates a shared coverage
//! envelope from per-sample abundances, so at least two samples are required
//! and phase 2 cannot run inside a single alignment.
//!
//! Reference: Xiong et al., "Normalization of generalized transcript
//! degradation improves accuracy in RNA-seq analysis", *Genome Biology* 2019;
//! implementation ported from <https://github.com/NUStatBioinfo/DegNorm>.

pub mod baseline;
pub mod nmf;
pub mod run;

use crate::params::Parameters;
use crate::quant::coverage::CoverageFile;

/// `--runMode degNorm` entry point.
pub fn run_mode(params: &Parameters) -> anyhow::Result<()> {
    let mut files = Vec::with_capacity(params.deg_norm_coverage_files.len());
    for path in &params.deg_norm_coverage_files {
        log::info!("degNorm: loading {}", path.display());
        files.push(CoverageFile::read(path)?);
    }

    let cfg = run::DegNormConfig {
        iter: params.deg_norm_iter,
        nmf_iter: params.deg_norm_nmf_iter,
        downsample_rate: params.deg_norm_downsample_rate.max(1),
        minimax_coverage: params.deg_norm_minimax_coverage,
        skip_baseline: params.deg_norm_skip_baseline_selection,
        bins: params.deg_norm_bins,
        min_high_coverage: params.deg_norm_min_high_coverage,
        seed: params.run_rng_seed,
    };

    let out = run::run_degnorm(&files, &cfg)?;
    let dir = params.output_path("DegNorm.out");
    run::write_outputs(&out, &dir, &cfg)?;
    log::info!(
        "degNorm: {} genes x {} samples written to {}",
        out.gene_ids.len(),
        out.sample_ids.len(),
        dir.display()
    );
    Ok(())
}
