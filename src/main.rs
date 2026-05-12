use clap::Parser;

use rustar_aligner::cpu;
use rustar_aligner::params::Parameters;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    let params = Parameters::parse();
    rustar_aligner::run(&params)
}
