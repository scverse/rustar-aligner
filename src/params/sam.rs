// ---------------------------------------------------------------------------
// SAM optional-tag attribute set (`--outSAMattributes`)
// ---------------------------------------------------------------------------

use std::str::FromStr;

/// One optional SAM tag requested via `--outSAMattributes`.
///
/// The discriminant is the bit position in [`SamAttributes`]'s presence mask.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum SamAttr {
    NH = 0,
    HI = 1,
    AS = 2,
    /// `NM:i`: SAM-standard edit distance (mismatches + inserted + deleted
    /// bases). NOT in STAR's `Standard` preset; opt-in via `NM` / `All`.
    NM = 3,
    MD = 4,
    JM = 5,
    JI = 6,
    XS = 7,
    RG = 8,
    /// STARsolo gene id of the read's Gene-feature assignment (GX:Z).
    GX = 9,
    /// STARsolo gene name (symbol) of the Gene-feature assignment (GN:Z).
    GN = 10,
    /// `nM:i`: STAR's mismatch count (mismatches only, excluding indels). This
    /// is the tag in STAR's `Standard` preset (distinct from `NM`).
    NMM = 11,
    /// WASP allele-specific-mapping tags (require --waspOutputMode SAMtag).
    VW = 12,
    VA = 13,
    VG = 14,
    // ---- STARsolo barcode tags (BAM output only, like STAR) ----
    /// `CR:Z`: raw (uncorrected) cell barcode sequence.
    CR = 15,
    /// `CY:Z`: quality string of the raw cell barcode.
    CY = 16,
    /// `UR:Z`: raw (uncorrected) UMI sequence.
    UR = 17,
    /// `UY:Z`: quality string of the raw UMI.
    UY = 18,
    /// `CB:Z`: whitelist-corrected cell barcode. Filled at sorting time from
    /// the solo read info (except `--soloType CB_samTagOut`, which corrects
    /// the barcode as the read is processed).
    CB = 19,
    /// `UB:Z`: collapsed (corrected) UMI. Only known after UMI collapsing,
    /// so it is added when the sorted BAM is written.
    UB = 20,
    /// `sM:i`: STAR's `cbMatch` code (its barcode/UMI assessment).
    SM = 21,
    /// `sS:Z`: full barcode-read sequence (CB + UMI + any adapter).
    SS = 22,
    /// `sQ:Z`: full barcode-read quality string.
    SQ = 23,
    /// `gx:Z`: gene ids of THIS alignment, `;`-joined (multi-gene allowed,
    /// unlike the read-level unique-gene `GX`).
    GXM = 24,
    /// `gn:Z`: gene names of this alignment, `;`-joined.
    GNM = 25,
    /// `sF:B:i`: `(overlap type, number of genes)` for the read.
    SF = 26,
}

impl SamAttr {
    /// The two-letter SAM tag.
    pub const fn tag(self) -> [u8; 2] {
        match self {
            Self::NH => *b"NH",
            Self::HI => *b"HI",
            Self::AS => *b"AS",
            Self::NM => *b"NM",
            Self::MD => *b"MD",
            Self::JM => *b"jM",
            Self::JI => *b"jI",
            Self::XS => *b"XS",
            Self::RG => *b"RG",
            Self::GX => *b"GX",
            Self::GN => *b"GN",
            Self::NMM => *b"nM",
            Self::VW => *b"vW",
            Self::VA => *b"vA",
            Self::VG => *b"vG",
            Self::CR => *b"CR",
            Self::CY => *b"CY",
            Self::UR => *b"UR",
            Self::UY => *b"UY",
            Self::CB => *b"CB",
            Self::UB => *b"UB",
            Self::SM => *b"sM",
            Self::SS => *b"sS",
            Self::SQ => *b"sQ",
            Self::GXM => *b"gx",
            Self::GNM => *b"gn",
            Self::SF => *b"sF",
        }
    }
}

pub const MAX_ATTRS: usize = 27;

/// Ordered set of optional SAM tags (`--outSAMattributes`).
///
/// STAR writes the optional tags in the order given on the command line
/// (`outSAMattrOrder`, Parameters_samAttributes.cpp), with the derived tags
/// (RG for `--outSAMattrRGline`, XS for `--outSAMstrandField intronMotif`, vW
/// for WASP) appended after it. This type keeps both the presence mask (for
/// `contains` checks) and that order (for the writers). `|` appends the
/// right-hand attributes that are not yet present, in their order, which is
/// how STAR's appends behave. `Copy`, so it is passed around like the former
/// bitflags value; `PartialEq` compares the order too.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamAttributes {
    mask: u32,
    len: u8,
    order: [SamAttr; MAX_ATTRS],
}

impl Default for SamAttributes {
    fn default() -> Self {
        Self::empty()
    }
}

impl SamAttributes {
    pub const NH: Self = Self::of(&[SamAttr::NH]);
    pub const HI: Self = Self::of(&[SamAttr::HI]);
    pub const AS: Self = Self::of(&[SamAttr::AS]);
    pub const NM: Self = Self::of(&[SamAttr::NM]);
    pub const MD: Self = Self::of(&[SamAttr::MD]);
    pub const JM: Self = Self::of(&[SamAttr::JM]);
    pub const JI: Self = Self::of(&[SamAttr::JI]);
    pub const XS: Self = Self::of(&[SamAttr::XS]);
    pub const RG: Self = Self::of(&[SamAttr::RG]);
    pub const GX: Self = Self::of(&[SamAttr::GX]);
    pub const GN: Self = Self::of(&[SamAttr::GN]);
    pub const NMM: Self = Self::of(&[SamAttr::NMM]);
    pub const VW: Self = Self::of(&[SamAttr::VW]);
    pub const VA: Self = Self::of(&[SamAttr::VA]);
    pub const VG: Self = Self::of(&[SamAttr::VG]);
    pub const CR: Self = Self::of(&[SamAttr::CR]);
    pub const CY: Self = Self::of(&[SamAttr::CY]);
    pub const UR: Self = Self::of(&[SamAttr::UR]);
    pub const UY: Self = Self::of(&[SamAttr::UY]);
    pub const CB: Self = Self::of(&[SamAttr::CB]);
    pub const UB: Self = Self::of(&[SamAttr::UB]);
    pub const SM: Self = Self::of(&[SamAttr::SM]);
    pub const SS: Self = Self::of(&[SamAttr::SS]);
    pub const SQ: Self = Self::of(&[SamAttr::SQ]);
    pub const GXM: Self = Self::of(&[SamAttr::GXM]);
    pub const GNM: Self = Self::of(&[SamAttr::GNM]);
    pub const SF: Self = Self::of(&[SamAttr::SF]);

    /// Every STARsolo barcode/gene tag. STAR emits these in BAM output only.
    pub const SOLO_TAGS: Self = Self::of(&[
        SamAttr::CR,
        SamAttr::CY,
        SamAttr::UR,
        SamAttr::UY,
        SamAttr::CB,
        SamAttr::UB,
        SamAttr::SM,
        SamAttr::SS,
        SamAttr::SQ,
        SamAttr::GX,
        SamAttr::GN,
        SamAttr::GXM,
        SamAttr::GNM,
        SamAttr::SF,
    ]);

    /// The tags derived from the gene model rather than the barcode read.
    pub const SOLO_GENE_TAGS: Self = Self::of(&[
        SamAttr::GX,
        SamAttr::GN,
        SamAttr::GXM,
        SamAttr::GNM,
        SamAttr::SF,
    ]);

    /// STAR `Standard` = NH HI AS nM  (the mismatch count nM, NOT edit-distance NM).
    pub const STANDARD: Self = Self::of(&[SamAttr::NH, SamAttr::HI, SamAttr::AS, SamAttr::NMM]);
    /// STAR `All` = NH HI AS nM NM MD jM jI (MC and ch are not implemented),
    /// followed here by XS, which the parameter fold keeps only under
    /// `--outSAMstrandField intronMotif`.
    pub const ALL: Self = Self::of(&[
        SamAttr::NH,
        SamAttr::HI,
        SamAttr::AS,
        SamAttr::NMM,
        SamAttr::NM,
        SamAttr::MD,
        SamAttr::JM,
        SamAttr::JI,
        SamAttr::XS,
    ]);

    /// No attributes.
    pub const fn empty() -> Self {
        Self {
            mask: 0,
            len: 0,
            order: [SamAttr::NH; MAX_ATTRS],
        }
    }

    /// Build from an ordered list (duplicates dropped, first occurrence wins).
    #[must_use]
    pub const fn of(list: &[SamAttr]) -> Self {
        let mut s = Self::empty();
        let mut i = 0;
        while i < list.len() {
            s = s.with(list[i]);
            i += 1;
        }
        s
    }

    /// Append one attribute if absent.
    #[must_use]
    pub const fn with(mut self, a: SamAttr) -> Self {
        let bit = 1u32 << (a as u8);
        if self.mask & bit == 0 {
            self.mask |= bit;
            self.order[self.len as usize] = a;
            self.len += 1;
        }
        self
    }

    /// True when no attribute is present.
    pub const fn is_empty(&self) -> bool {
        self.mask == 0
    }

    /// True when every attribute of `other` is present (order ignored).
    pub const fn contains(&self, other: Self) -> bool {
        self.mask & other.mask == other.mask
    }

    /// True when at least one attribute of `other` is present.
    pub const fn intersects(&self, other: Self) -> bool {
        self.mask & other.mask != 0
    }

    /// The attributes in output order.
    pub fn iter(&self) -> impl Iterator<Item = SamAttr> + '_ {
        self.order[..self.len as usize].iter().copied()
    }

    /// Remove the attributes of `other`, keeping the order of the rest.
    pub fn remove(&mut self, other: Self) {
        let mut out = Self::empty();
        for a in self.iter() {
            if other.mask & (1u32 << (a as u8)) == 0 {
                out = out.with(a);
            }
        }
        *self = out;
    }
}

impl std::ops::BitOr for SamAttributes {
    type Output = Self;
    fn bitor(mut self, rhs: Self) -> Self {
        for a in rhs.iter() {
            self = self.with(a);
        }
        self
    }
}

impl std::ops::BitAnd for SamAttributes {
    type Output = Self;
    /// The attributes of `self` that are also in `rhs`, in `self`'s order.
    fn bitand(self, rhs: Self) -> Self {
        let mut out = Self::empty();
        for a in self.iter() {
            if rhs.mask & (1u32 << (a as u8)) != 0 {
                out = out.with(a);
            }
        }
        out
    }
}

impl std::ops::BitOrAssign for SamAttributes {
    fn bitor_assign(&mut self, rhs: Self) {
        *self = *self | rhs;
    }
}

impl std::ops::Sub for SamAttributes {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self {
        self.remove(rhs);
        self
    }
}

impl FromStr for SamAttributes {
    type Err = String;
    /// Parse a single CLI token into a flag. `Standard`/`All`/`None` expand to
    /// their preset sets; individual tag names map to a single attribute.
    fn from_str(s: &str) -> Result<Self, String> {
        Ok(match s {
            // STAR's `All` has no XS (it is appended later by intronMotif).
            "All" => Self::ALL - Self::XS,
            "Standard" => Self::STANDARD,
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
