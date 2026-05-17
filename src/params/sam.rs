// ---------------------------------------------------------------------------
// SAM output type enums
// ---------------------------------------------------------------------------

/// STAR's `--outSAMtype` format component.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum OutSamFormat {
    #[default]
    Sam,
    Bam,
    None,
}

/// STAR's `--outSAMtype` sort order component (only applies to BAM).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutSamSortOrder {
    Unsorted,
    SortedByCoordinate,
}

/// Combined `--outSAMtype` value.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct OutSamType {
    pub format: OutSamFormat,
    pub sort_order: Option<OutSamSortOrder>,
}

impl std::fmt::Display for OutSamType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (&self.format, &self.sort_order) {
            (OutSamFormat::Sam, _) => write!(f, "SAM"),
            (OutSamFormat::None, _) => write!(f, "None"),
            (OutSamFormat::Bam, Some(OutSamSortOrder::SortedByCoordinate)) => {
                write!(f, "BAM SortedByCoordinate")
            }
            (OutSamFormat::Bam, _) => write!(f, "BAM Unsorted"),
        }
    }
}

impl clap::FromArgMatches for OutSamType {
    fn from_arg_matches(matches: &clap::ArgMatches) -> Result<Self, clap::Error> {
        let mut s = Self::default();
        s.update_from_arg_matches(matches)?;
        Ok(s)
    }

    fn update_from_arg_matches(&mut self, matches: &clap::ArgMatches) -> Result<(), clap::Error> {
        let Some(values) = matches.get_many::<String>("outSAMtype") else {
            return Ok(());
        };
        let tokens: Vec<&str> = values.map(String::as_str).collect();
        *self = match tokens.as_slice() {
            ["SAM"] => Self {
                format: OutSamFormat::Sam,
                sort_order: None,
            },
            ["None"] => Self {
                format: OutSamFormat::None,
                sort_order: None,
            },
            ["BAM", "Unsorted"] => Self {
                format: OutSamFormat::Bam,
                sort_order: Some(OutSamSortOrder::Unsorted),
            },
            ["BAM", "SortedByCoordinate"] => Self {
                format: OutSamFormat::Bam,
                sort_order: Some(OutSamSortOrder::SortedByCoordinate),
            },
            other => {
                return Err(invalid_multi_arg(
                    other,
                    &["SAM", "None", "BAM Unsorted", "BAM SortedByCoordinate"],
                ));
            }
        };
        Ok(())
    }
}

impl clap::Args for OutSamType {
    fn augment_args(cmd: clap::Command) -> clap::Command {
        cmd.arg(
            clap::Arg::new("outSAMtype")
                .long("outSAMtype")
                .num_args(1..=2)
                .default_values(["SAM"])
                .help(
                    "Output type: SAM, BAM Unsorted, BAM SortedByCoordinate, None. \
                     Provide as space-separated tokens, e.g. `--outSAMtype BAM SortedByCoordinate`.",
                ),
        )
    }

    fn augment_args_for_update(cmd: clap::Command) -> clap::Command {
        Self::augment_args(cmd)
    }
}

// ---------------------------------------------------------------------------
// SAM unmapped output
// ---------------------------------------------------------------------------

/// STAR’s `--outSAMunmapped` value
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum OutSamUnmapped {
    #[default]
    None,
    Within,
    WithinKeepPairs,
}

impl clap::FromArgMatches for OutSamUnmapped {
    fn from_arg_matches(matches: &clap::ArgMatches) -> Result<Self, clap::Error> {
        let mut s = Self::default();
        s.update_from_arg_matches(matches)?;
        Ok(s)
    }

    fn update_from_arg_matches(&mut self, matches: &clap::ArgMatches) -> Result<(), clap::Error> {
        let Some(values) = matches.get_many::<String>("outSAMunmapped") else {
            return Ok(());
        };
        let tokens: Vec<&str> = values.map(String::as_str).collect();
        *self = match tokens.as_slice() {
            ["None"] => Self::None,
            ["Within"] => Self::Within,
            ["Within", "KeepPairs"] => Self::WithinKeepPairs,
            other => {
                return Err(invalid_multi_arg(
                    other,
                    &["None", "Within", "Within KeepPairs"],
                ));
            }
        };
        Ok(())
    }
}

impl clap::Args for OutSamUnmapped {
    fn augment_args(cmd: clap::Command) -> clap::Command {
        cmd.arg(
            clap::Arg::new("outSAMunmapped")
                .long("outSAMunmapped")
                .num_args(1..=2)
                .default_values(["None"])
                .help(
                    "Unmapped reads in SAM output: None, Within, or Within KeepPairs. \
                     Provide as space-separated tokens, e.g. `--outSAMunmapped Within KeepPairs`.",
                ),
        )
    }

    fn augment_args_for_update(cmd: clap::Command) -> clap::Command {
        Self::augment_args(cmd)
    }
}

// Helpers

fn invalid_multi_arg(other: &[&str], valid: &[&str]) -> clap::error::Error {
    use clap::error::{ContextKind, ContextValue, ErrorKind};

    let mut err = clap::Error::new(ErrorKind::InvalidValue);
    err.insert(
        ContextKind::InvalidArg,
        ContextValue::String("--outSAMtype".into()),
    );
    err.insert(
        ContextKind::InvalidValue,
        ContextValue::String(other.join(" ")),
    );
    err.insert(
        ContextKind::ValidValue,
        // replace spaces with an invisible non-whitespace character to prevent clap from adding quotes
        ContextValue::Strings(valid.iter().map(|s| s.replace(' ', "\u{2800}")).collect()),
    );
    err
}
