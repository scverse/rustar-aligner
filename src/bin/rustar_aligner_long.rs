use rustar_aligner::cpu;
use rustar_aligner::params::Parameters;

/// Global allocator override — see `main.rs` / `Cargo.toml`'s `mimalloc` comment.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// STAR's `STARlong` binary: identical CLI to `rustar-aligner`, but forces the
/// chaining-DP long-read stitcher (STAR's `-DCOMPILE_FOR_LONG_READS`), matching
/// how native STAR selects it by which binary you run, not a CLI flag.
fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    let mut params = Parameters::parse();
    params.long_read = true;
    rustar_aligner::run(&params)
}
