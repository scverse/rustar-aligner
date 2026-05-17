use rustar_aligner::cpu;
use rustar_aligner::params::Parameters;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let mut params = Parameters::parse();
    params.command_line = Some(command_line);
    rustar_aligner::run(&params)
}
