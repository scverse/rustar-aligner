use rustar_aligner::cpu;
use rustar_aligner::params::Parameters;

/// Global allocator override — see the comment in `Cargo.toml`
/// next to the `mimalloc` dependency for the rationale. The
/// `#[global_allocator]` attribute applies to the whole binary; no
/// further work needed at allocation sites.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    rustar_aligner::run(&Parameters::parse())
}
