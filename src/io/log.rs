//! Writers for STAR-compatible `Log.out` and `Log.progress.out`.

use crate::genome::Genome;
use crate::params::Parameters;
use crate::stats::AlignmentStats;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

/// Parse `--param value…` pairs out of the raw command-line string.
///
/// Stops accumulating values for a parameter when the next `--` token is hit.
fn cli_params(cmd: &str) -> Vec<(String, String)> {
    let mut result = Vec::new();
    let tokens: Vec<&str> = cmd.split_whitespace().collect();
    let mut i = 0;
    // skip the binary name
    if i < tokens.len() && !tokens[i].starts_with("--") {
        i += 1;
    }
    while i < tokens.len() {
        if let Some(name) = tokens[i].strip_prefix("--") {
            i += 1;
            let mut vals = Vec::new();
            while i < tokens.len() && !tokens[i].starts_with("--") {
                vals.push(tokens[i]);
                i += 1;
            }
            result.push((name.to_string(), vals.join(" ")));
        } else {
            i += 1;
        }
    }
    result
}

/// Write STAR-compatible `Log.out`.
///
/// Content mirrors STAR's Log.out structure: version header, command line,
/// parameter sections, chromosome table, phase timestamps, and "ALL DONE!".
pub fn write_log_out(
    path: &Path,
    params: &Parameters,
    genome: &Genome,
    time_start: chrono::DateTime<chrono::Local>,
    time_genome_loaded: chrono::DateTime<chrono::Local>,
    time_finish: chrono::DateTime<chrono::Local>,
) -> std::io::Result<()> {
    let f = std::fs::File::create(path)?;
    let mut w = BufWriter::new(f);

    // STAR uses two timestamp formats in Log.out
    let long_fmt = "%a %b %e %H:%M:%S %Y"; // "Tue Feb 10 17:11:24 2026"
    let short_fmt = "%b %e %H:%M:%S"; //      "Feb 10 17:11:26"

    // ── Header ─────────────────────────────────────────────────────────────
    writeln!(w, "STAR version={}", env!("CARGO_PKG_VERSION"))?;
    writeln!(
        w,
        "STAR compilation time,server,dir={} :",
        time_start.format("%Y-%m-%dT%H:%M:%S%:z")
    )?;
    writeln!(w, "STAR git: ")?;

    // ── Command line and parameter sections ────────────────────────────────
    let cmd = params.command_line.as_deref().unwrap_or("");
    let pairs = cli_params(cmd);

    writeln!(w, "##### Command Line:")?;
    writeln!(w, "{cmd}")?;

    writeln!(w, "##### Initial USER parameters from Command Line:")?;
    if let Some((_, v)) = pairs.iter().find(|(k, _)| k == "outFileNamePrefix") {
        writeln!(w, "{:<33}{v}", "outFileNamePrefix")?;
    }

    writeln!(w, "###### All USER parameters from Command Line:")?;
    for (k, v) in &pairs {
        writeln!(w, "{k:<30}{v}     ~RE-DEFINED")?;
    }
    writeln!(w, "##### Finished reading parameters from all sources")?;
    writeln!(w)?;

    writeln!(
        w,
        "##### Final user re-defined parameters-----------------:"
    )?;
    for (k, v) in &pairs {
        writeln!(w, "{k:<34}{v}")?;
    }
    writeln!(w)?;
    writeln!(w, "-------------------------------")?;
    writeln!(w, "##### Final effective command line:")?;
    writeln!(w, "{cmd}")?;
    writeln!(w, "----------------------------------------")?;
    writeln!(w)?;

    // ── Chromosome table ───────────────────────────────────────────────────
    writeln!(
        w,
        "Number of real (reference) chromosomes= {}",
        genome.n_chr_real
    )?;
    for i in 0..genome.n_chr_real {
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            i + 1,
            genome.chr_name[i],
            genome.chr_length[i],
            genome.chr_start[i]
        )?;
    }

    // ── Phase timestamps ───────────────────────────────────────────────────
    writeln!(
        w,
        "Started loading the genome: {}",
        time_genome_loaded.format(long_fmt)
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "Finished loading the genome: {}",
        time_genome_loaded.format(long_fmt)
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "{} ..... finished mapping",
        time_finish.format(short_fmt)
    )?;
    writeln!(w, "ALL DONE!")?;

    Ok(())
}

/// Write STAR-compatible `Log.progress.out`.
///
/// STAR updates this file periodically during alignment; for short runs (and
/// in our current single-pass implementation) it contains only the header,
/// one summary line with final counts, and "ALL DONE!".
pub fn write_log_progress_out(
    path: &Path,
    stats: &Arc<AlignmentStats>,
    time_start: chrono::DateTime<chrono::Local>,
    time_finish: chrono::DateTime<chrono::Local>,
) -> std::io::Result<()> {
    use std::sync::atomic::Ordering;

    let f = std::fs::File::create(path)?;
    let mut w = BufWriter::new(f);

    // STAR's exact two-line header
    writeln!(
        w,
        "           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped"
    )?;
    writeln!(
        w,
        "                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other"
    )?;

    // Final summary line — mirrors what STAR would write in its last progress tick
    let total = stats.total_reads.load(Ordering::Relaxed);
    if total > 0 {
        let elapsed_secs = (time_finish - time_start).num_seconds().max(1) as f64;
        let speed_m_hr = (total as f64 / elapsed_secs) * 3600.0 / 1_000_000.0;

        let read_len = stats.read_bases.load(Ordering::Relaxed) / total;
        let unique = stats.uniquely_mapped.load(Ordering::Relaxed);
        let mapped_len = if unique > 0 {
            stats.mapped_bases.load(Ordering::Relaxed) / unique
        } else {
            0
        };
        let multi = stats.multi_mapped.load(Ordering::Relaxed);
        let mm_rate = if unique + multi > 0 {
            (multi as f64 / (unique + multi) as f64) * 100.0
        } else {
            0.0
        };
        let too_short = stats.unmapped_short.load(Ordering::Relaxed);
        let too_many_mm = stats.unmapped_mismatches.load(Ordering::Relaxed);
        let other = stats.unmapped_other.load(Ordering::Relaxed)
            + stats.too_many_loci.load(Ordering::Relaxed);

        let elapsed = time_finish - time_start;
        let h = elapsed.num_hours();
        let m = elapsed.num_minutes() % 60;
        let s = elapsed.num_seconds() % 60;

        writeln!(
            w,
            "{h:>15}:{m:02}:{s:02}  {speed_m_hr:>11.2}  {total:>11}  {read_len:>7}  {unique:>7}  {mapped_len:>7}  {mm_rate:>7.2}%  {multi:>7}        0  {too_many_mm:>7}  {too_short:>7}  {other:>7}"
        )?;
    }

    writeln!(w, "ALL DONE!")?;

    Ok(())
}
