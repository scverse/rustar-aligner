// ---------------------------------------------------------------------------
// SAM optional-tag attribute set (`--outSAMattributes`)
// ---------------------------------------------------------------------------

use std::str::FromStr;

bitflags::bitflags! {
    /// Optional SAM tags requested via `--outSAMattributes`.
    ///
    /// Each bit corresponds to one tag the writer may emit. `STANDARD` and
    /// `ALL` are convenience aliases matching STAR's preset names.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    pub struct SamAttributes: u32 {
        const NH = 1 << 0;
        const HI = 1 << 1;
        const AS = 1 << 2;
        /// `NM:i` — SAM-standard edit distance (mismatches + inserted + deleted
        /// bases). NOT in STAR's `Standard` preset; opt-in via `NM` / `All`.
        const NM = 1 << 3;
        const MD = 1 << 4;
        const JM = 1 << 5;
        const JI = 1 << 6;
        const XS = 1 << 7;
        const RG = 1 << 8;
        /// STARsolo gene id of the read's Gene-feature assignment (GX:Z).
        const GX = 1 << 9;
        /// STARsolo gene name (symbol) of the Gene-feature assignment (GN:Z).
        const GN = 1 << 10;
        /// `nM:i` — STAR's mismatch count (mismatches only, excluding indels). This
        /// is the tag in STAR's `Standard` preset (distinct from `NM`).
        const NMM = 1 << 11;
        /// WASP allele-specific-mapping tags (require --waspOutputMode SAMtag).
        const VW = 1 << 12;
        const VA = 1 << 13;
        const VG = 1 << 14;

        // ---- STARsolo barcode tags (BAM output only, like STAR) ----
        /// `CR:Z`: raw (uncorrected) cell barcode sequence.
        const CR = 1 << 15;
        /// `CY:Z`: quality string of the raw cell barcode.
        const CY = 1 << 16;
        /// `UR:Z`: raw (uncorrected) UMI sequence.
        const UR = 1 << 17;
        /// `UY:Z`: quality string of the raw UMI.
        const UY = 1 << 18;
        /// `CB:Z`: whitelist-corrected cell barcode. Filled at sorting time from
        /// the solo read info (except `--soloType CB_samTagOut`, which corrects
        /// the barcode as the read is processed).
        const CB = 1 << 19;
        /// `UB:Z`: collapsed (corrected) UMI. Only known after UMI collapsing,
        /// so it is added when the sorted BAM is written.
        const UB = 1 << 20;
        /// `sM:i`: STAR's `cbMatch` code (its barcode/UMI assessment).
        const SM = 1 << 21;
        /// `sS:Z`: full barcode-read sequence (CB + UMI + any adapter).
        const SS = 1 << 22;
        /// `sQ:Z`: full barcode-read quality string.
        const SQ = 1 << 23;
        /// `gx:Z`: gene ids of THIS alignment, `;`-joined (multi-gene allowed,
        /// unlike the read-level unique-gene `GX`).
        const GXM = 1 << 24;
        /// `gn:Z`: gene names of this alignment, `;`-joined.
        const GNM = 1 << 25;
        /// `sF:B:i`: `(overlap type, number of genes)` for the read.
        const SF = 1 << 26;

        /// Every STARsolo barcode/gene tag. STAR emits these in BAM output only.
        const SOLO_TAGS =
            Self::CR.bits() | Self::CY.bits() | Self::UR.bits() | Self::UY.bits()
            | Self::CB.bits() | Self::UB.bits()
            | Self::SM.bits() | Self::SS.bits() | Self::SQ.bits()
            | Self::GX.bits() | Self::GN.bits()
            | Self::GXM.bits() | Self::GNM.bits() | Self::SF.bits();

        /// The tags derived from the gene model rather than the barcode read.
        const SOLO_GENE_TAGS =
            Self::GX.bits() | Self::GN.bits()
            | Self::GXM.bits() | Self::GNM.bits() | Self::SF.bits();

        // STAR `Standard` = NH HI AS nM  (the mismatch count nM, NOT edit-distance NM).
        const STANDARD =
            Self::NH.bits() | Self::HI.bits() | Self::AS.bits()
            | Self::NMM.bits();
        // STAR `All` additionally includes the edit-distance NM, MD, jM, jI (+ XS here).
        const ALL =
            Self::STANDARD.bits()
            | Self::NM.bits()
            | Self::MD.bits() | Self::JM.bits() | Self::JI.bits() | Self::XS.bits();
    }
}

impl FromStr for SamAttributes {
    type Err = String;
    /// Parse a single CLI token into a flag. `Standard`/`All`/`None` expand to
    /// their preset sets; individual tag names map to a single bit.
    fn from_str(s: &str) -> Result<Self, String> {
        Ok(match s {
            "Standard" => Self::STANDARD,
            "All" => Self::ALL,
            "None" => Self::empty(),
            "NH" => Self::NH,
            "HI" => Self::HI,
            "AS" => Self::AS,
            // `NM` = SAM-standard edit distance; `nM` = STAR's mismatch count. These
            // are distinct tags (only `nM` is in the `Standard` preset).
            "NM" => Self::NM,
            "nM" => Self::NMM,
            "MD" => Self::MD,
            "jM" => Self::JM,
            "jI" => Self::JI,
            "XS" => Self::XS,
            "RG" => Self::RG,
            "GX" => Self::GX,
            "GN" => Self::GN,
            "CR" => Self::CR,
            "CY" => Self::CY,
            "UR" => Self::UR,
            "UY" => Self::UY,
            "CB" => Self::CB,
            "UB" => Self::UB,
            "sM" => Self::SM,
            "sS" => Self::SS,
            "sQ" => Self::SQ,
            "gx" => Self::GXM,
            "gn" => Self::GNM,
            "sF" => Self::SF,
            "vW" => Self::VW,
            "vA" => Self::VA,
            "vG" => Self::VG,
            other => return Err(format!("unknown --outSAMattributes token '{other}'")),
        })
    }
}

impl clap::FromArgMatches for SamAttributes {
    fn from_arg_matches(matches: &clap::ArgMatches) -> Result<Self, clap::Error> {
        let mut s = Self::STANDARD;
        s.update_from_arg_matches(matches)?;
        Ok(s)
    }

    fn update_from_arg_matches(&mut self, matches: &clap::ArgMatches) -> Result<(), clap::Error> {
        let Some(values) = matches.get_many::<String>("outSAMattributes") else {
            return Ok(());
        };
        let mut acc = Self::empty();
        for tok in values {
            let flag = tok.parse().map_err(|e| {
                use clap::error::{ContextKind, ContextValue, ErrorKind};
                let mut err = clap::Error::new(ErrorKind::InvalidValue);
                err.insert(
                    ContextKind::InvalidArg,
                    ContextValue::String("--outSAMattributes".into()),
                );
                err.insert(ContextKind::InvalidValue, ContextValue::String(tok.clone()));
                err.insert(ContextKind::Custom, ContextValue::String(e));
                err
            })?;
            acc |= flag;
        }
        *self = acc;
        Ok(())
    }
}

impl clap::Args for SamAttributes {
    fn augment_args(cmd: clap::Command) -> clap::Command {
        cmd.arg(
            clap::Arg::new("outSAMattributes")
                .long("outSAMattributes")
                .num_args(1..)
                .default_values(["Standard"])
                .help(
                    "SAM optional tags: Standard, All, None, or any combination of \
                     NH HI AS NM nM MD jM jI XS RG vW vA vG, plus the STARsolo tags \
                     CR CY UR UY CB UB GX GN gx gn sM sS sQ sF (BAM output only).",
                ),
        )
    }

    fn augment_args_for_update(cmd: clap::Command) -> clap::Command {
        Self::augment_args(cmd)
    }
}

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
