//! Locks the CLI parameter surface against STAR 2.7.11b.
//!
//! `tests/data/star_2.7.11b_params.txt` lists every column-1 name of STAR
//! 2.7.11b's `parametersDefault`, with none removed. This test asserts that clap recognises each one, so
//! the coverage figure is a machine-checked number rather than a claim, and so
//! surface drift shows up the moment it happens rather than in a bug report.
//!
//! A recognised name here means "the CLI accepts it and does not error", not
//! "it is fully implemented". Parameters that are deliberately accepted but
//! output-neutral are listed in `ACCEPTED_BUT_INERT` below, with the reason.
//! Parameters not yet accepted at all live in `NOT_YET_ACCEPTED`, which is the
//! honest remaining gap.

use std::collections::BTreeSet;

/// Names accepted by the CLI but with no effect on output, each for a stated
/// reason. Keeping the list explicit stops "accepted" from quietly meaning
/// "implemented".
///
/// Only half of that claim is machine-checked: `inert_parameters_are_actually
/// _accepted` verifies each name *is* accepted, but nothing verifies it is
/// inert. Adding a name here is an assertion made by hand, so check that the
/// parameter is genuinely unread before adding one — `limitBAMsortRAM` and
/// `runRNGseed` were both listed here while being implemented.
const ACCEPTED_BUT_INERT: &[(&str, &str)] = &[
    (
        "genomeFileSizes",
        "an input hint for STAR's loader; no output bytes depend on it",
    ),
    (
        "limitGenomeGenerateRAM",
        "genome generation manages its own memory",
    ),
    (
        "limitIObufferSize",
        "buffer sizing is an implementation detail with no output effect",
    ),
    (
        "limitNreadsSoft",
        "a soft warning threshold in STAR; no output bytes depend on it",
    ),
    (
        "limitOutSAMoneReadBytes",
        "an overflow guard on STAR's fixed output buffer",
    ),
    (
        "limitOutSJcollapsed",
        "an allocation cap on STAR's collapsed-junction array",
    ),
    (
        "limitOutSJoneRead",
        "an allocation cap on STAR's per-read junction array",
    ),
    (
        "limitSjdbInsertNsj",
        "an allocation cap on the inserted-junction array",
    ),
    (
        "outBAMsortingBinsN",
        "sorting is not binned yet; output is unaffected",
    ),
    (
        "outBAMsortingThreadN",
        "BGZF writing is single-threaded; output is unaffected",
    ),
    (
        "outTmpDir",
        "no intermediate files are written that a user could observe",
    ),
    ("outTmpKeep", "as outTmpDir: nothing to keep"),
    (
        "readMatesLengthsIn",
        "a read-length hint; lengths are taken from the FASTQ",
    ),
    (
        "runDirPerm",
        "no directories are created whose mode a user could observe",
    ),
];

/// Names STAR accepts that this CLI does not yet accept at all. This is the
/// real remaining surface gap; shrinking it is the point of the port.
///
/// Adding a name here must always be a deliberate act. Removing one is what
/// progress looks like.
const NOT_YET_ACCEPTED: &[&str] = &[
    // STAR meta-parameters, none of them accepted here. clap rejects them, so
    // a user who passes one is told rather than quietly ignored, which is the
    // behaviour these three need most: silently dropping `--parametersFiles`
    // would discard every parameter in that file.
    "parametersFiles",
    "sysShell",
    "versionGenome",
    // Aligner core (annotated-junction stitching, alignEndsType, in-recursion
    // length penalty).
    "alignEndsProtrude",
    "alignInsertionFlush",
    "alignSoftClipAtReferenceEnds",
    "alignTranscriptsPerReadNmax",
    "outFilterMismatchNoverReadLmax",
    "seedNoneLociPerWindow",
    "seedSplitMin",
    // Chimeric multimapping.
    "chimFilter",
    "chimMultimapNmax",
    "chimMultimapScoreRange",
    "chimNonchimScoreDropMin",
    // CellRanger4 adapter clipping.
    "clip5pAdapterMMp",
    "clip5pAdapterSeq",
    // Genome index types and transforms.
    "genomeSuffixLengthMax",
    "genomeTransformOutput",
    "genomeType",
    "sjdbInsertSave",
    // STARsolo barcode chemistry.
    "soloAdapterMismatchesNmax",
    "soloAdapterSequence",
    "soloCBtype",
    "soloOutFormatFeaturesGeneField3",
];

fn star_parameter_names() -> Vec<String> {
    let raw = include_str!("data/star_2.7.11b_params.txt");
    raw.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// Every long flag clap knows about, walked recursively so flattened argument
/// groups are included.
fn recognised_flags() -> BTreeSet<String> {
    use clap::CommandFactory;
    fn walk(cmd: &clap::Command, out: &mut BTreeSet<String>) {
        for arg in cmd.get_arguments() {
            if let Some(long) = arg.get_long() {
                out.insert(long.to_string());
            }
            for alias in arg.get_all_aliases().unwrap_or_default() {
                out.insert(alias.to_string());
            }
        }
        for sub in cmd.get_subcommands() {
            walk(sub, out);
        }
    }
    let cmd = rustar_aligner::params::Parameters::command();
    let mut out = BTreeSet::new();
    walk(&cmd, &mut out);
    out
}

#[test]
fn star_parameter_surface_is_fully_accounted_for() {
    let star = star_parameter_names();
    let ours = recognised_flags();
    let not_yet: BTreeSet<&str> = NOT_YET_ACCEPTED.iter().copied().collect();

    let mut unexpected_missing = Vec::new();
    let mut unexpectedly_present = Vec::new();

    for name in &star {
        let known = ours.contains(name.as_str());
        let declared_missing = not_yet.contains(name.as_str());
        if !known && !declared_missing {
            unexpected_missing.push(name.clone());
        }
        if known && declared_missing {
            unexpectedly_present.push(name.clone());
        }
    }

    assert!(
        unexpected_missing.is_empty(),
        "{} STAR parameter(s) are not recognised and are not declared in \
         NOT_YET_ACCEPTED. Either implement them or add them to that list \
         deliberately:\n  {}",
        unexpected_missing.len(),
        unexpected_missing.join("\n  ")
    );

    assert!(
        unexpectedly_present.is_empty(),
        "{} parameter(s) are now recognised but still listed in \
         NOT_YET_ACCEPTED; remove them from that list:\n  {}",
        unexpectedly_present.len(),
        unexpectedly_present.join("\n  ")
    );
}

#[test]
fn inert_parameters_are_actually_accepted() {
    // A name in ACCEPTED_BUT_INERT that the CLI does not accept would be a
    // documentation lie, so check the claim rather than trusting it.
    let ours = recognised_flags();
    let missing: Vec<&str> = ACCEPTED_BUT_INERT
        .iter()
        .map(|(name, _)| *name)
        .filter(|name| !ours.contains(*name))
        .collect();
    assert!(
        missing.is_empty(),
        "listed as accepted-but-inert, but not accepted: {missing:?}"
    );
}

#[test]
fn report_coverage() {
    // Not an assertion: prints the machine-checked coverage figure so it can be
    // quoted without anyone counting by hand.
    let star = star_parameter_names();
    let ours = recognised_flags();
    let covered = star.iter().filter(|n| ours.contains(n.as_str())).count();
    println!(
        "STAR 2.7.11b parameter surface: {covered}/{} accepted ({} inert), \
         {} not yet accepted",
        star.len(),
        ACCEPTED_BUT_INERT.len(),
        star.len() - covered,
    );
    assert!(covered > 0);
}
