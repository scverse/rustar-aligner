//! End-to-end checks for the SAM/SJ/read-input knobs.
//!
//! These drive `Parameters` directly rather than spawning the binary: the
//! knobs under test are all decided at parse time or at a single output site,
//! so a full alignment run would add minutes without adding coverage.

use rustar_aligner::params::Parameters;

fn parse(args: &[&str]) -> Result<Parameters, clap::Error> {
    let mut full = vec!["rustar-aligner"];
    full.extend_from_slice(args);
    Parameters::try_parse_from(&full)
}

fn with_reads(args: &[&str]) -> Result<Parameters, clap::Error> {
    let mut full = vec!["--readFilesIn", "reads.fq"];
    full.extend_from_slice(args);
    parse(&full)
}

#[test]
fn out_sam_mode_accepts_stars_three_values_and_rejects_others() {
    for value in ["Full", "NoQS", "None"] {
        assert!(
            with_reads(&["--outSAMmode", value]).is_ok(),
            "--outSAMmode {value} should parse"
        );
    }
    assert!(with_reads(&["--outSAMmode", "Quiet"]).is_err());
}

#[test]
fn out_sj_knobs_accept_stars_values_and_reject_others() {
    assert_eq!(with_reads(&[]).unwrap().out_sj_type, "Standard");
    assert!(with_reads(&["--outSJtype", "None"]).is_ok());
    assert!(with_reads(&["--outSJtype", "Compact"]).is_err());

    assert_eq!(with_reads(&[]).unwrap().out_sj_filter_reads, "All");
    assert!(with_reads(&["--outSJfilterReads", "Unique"]).is_ok());
    assert!(with_reads(&["--outSJfilterReads", "Best"]).is_err());
}

#[test]
fn read_files_prefix_is_applied_to_every_input_path() {
    let p = with_reads(&["--readFilesPrefix", "/data/run7/"]).unwrap();
    assert_eq!(p.read_files_in[0].to_str().unwrap(), "/data/run7/reads.fq");

    // Paired input: both mates get the prefix.
    let p = parse(&[
        "--readFilesIn",
        "r1.fq",
        "r2.fq",
        "--readFilesPrefix",
        "/data/run7/",
    ])
    .unwrap();
    assert_eq!(p.read_files_in[0].to_str().unwrap(), "/data/run7/r1.fq");
    assert_eq!(p.read_files_in[1].to_str().unwrap(), "/data/run7/r2.fq");

    // Absent by default.
    let p = with_reads(&[]).unwrap();
    assert_eq!(p.read_files_in[0].to_str().unwrap(), "reads.fq");
}

#[test]
fn read_name_separator_defaults_to_stars_slash() {
    let p = with_reads(&[]).unwrap();
    assert_eq!(p.read_name_separator, vec!["/".to_string()]);

    let p = with_reads(&["--readNameSeparator", "-"]).unwrap();
    assert_eq!(p.read_name_separator, vec!["-".to_string()]);
}

#[test]
fn read_quality_score_base_is_restricted_to_the_two_real_encodings() {
    assert_eq!(with_reads(&[]).unwrap().read_quality_score_base, 33);
    assert!(with_reads(&["--readQualityScoreBase", "64"]).is_ok());
    assert!(with_reads(&["--readQualityScoreBase", "0"]).is_ok());
    assert!(with_reads(&["--readQualityScoreBase", "40"]).is_err());
}

#[test]
fn out_qs_conversion_add_takes_negative_values() {
    // -31 is the Phred+64 to Phred+33 conversion, and needs to survive clap's
    // hyphen handling.
    let p = with_reads(&["--outQSconversionAdd", "-31"]).unwrap();
    assert_eq!(p.out_qs_conversion_add, -31);
}

#[test]
fn unsupported_modes_are_refused_rather_than_ignored() {
    // Both need machinery this aligner does not have. Accepting and ignoring
    // them would silently produce output the user did not ask for.
    assert!(with_reads(&["--outSAMfilter", "KeepOnlyAddedReferences"]).is_err());
    assert!(with_reads(&["--outSAMfilter", "KeepAllAddedReferences"]).is_err());
    assert!(with_reads(&["--readFilesType", "SAM", "SE"]).is_err());

    // The defaults, and the supported input type, still parse.
    assert!(with_reads(&["--outSAMfilter", "None"]).is_ok());
    assert!(with_reads(&["--readFilesType", "Fastx"]).is_ok());
}

#[test]
fn inert_limit_knobs_parse_with_stars_defaults() {
    let p = with_reads(&[]).unwrap();
    assert_eq!(p.limit_out_sam_one_read_bytes, 100_000);
    assert_eq!(p.limit_out_sj_collapsed, 1_000_000);
    assert_eq!(p.limit_out_sj_one_read, 1_000);
    assert_eq!(p.limit_sjdb_insert_nsj, 1_000_000);
    assert_eq!(p.limit_nreads_soft, -1);
    assert_eq!(p.limit_io_buffer_size, vec![30_000_000, 50_000_000]);
    assert_eq!(p.run_dir_perm, "User_RWX");
    assert_eq!(p.out_tmp_dir, "-");
    assert_eq!(p.out_tmp_keep, "None");
    assert_eq!(p.out_bam_sorting_bins_n, 50);
    assert_eq!(p.out_bam_sorting_thread_n, 0);

    // And they accept non-default values without complaint, since STAR does.
    assert!(with_reads(&["--limitOutSJcollapsed", "2000000"]).is_ok());
    assert!(with_reads(&["--limitNreadsSoft", "-1"]).is_ok());
    assert!(with_reads(&["--outBAMsortingThreadN", "8"]).is_ok());
}
