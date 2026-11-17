//! `--soloCellReadStats Standard`: the per-cell-barcode read summary STARsolo
//! writes as `Solo.out/<feature>/CellReads.stats`.
//!
//! One row per cell barcode, fifteen counters describing what happened to the
//! reads carrying it (how the barcode matched, whether the read mapped
//! uniquely, whether it landed on a feature, and whether it reached the
//! matrix) plus the per-cell UMI and gene totals. Reads whose barcode never
//! resolved are summed into a single `CBnotInPasslist` row, so `cbMatch`
//! accounts for every input read.
//!
//! # What STAR counts
//!
//! Each feature keeps its own table (`SoloReadFeature::readFlag`), because the
//! flags depend on the feature's annotation of the read:
//!
//! - At the read stage (`SoloReadFeature_record.cpp:25-62`) every read gets
//!   `genomeU`/`genomeM` from its number of loci (0 for an unmapped or
//!   too-many-loci read), `mito` if any alignment is on a
//!   `--genomeChrSetMitochondrial` chromosome, a region bit from the feature's
//!   `ovType`, and `featureU`/`featureM` from the size of the feature's gene
//!   set. Only `Gene` sets an `ovType` (`exonic`, whenever it found a gene,
//!   `Transcriptome_classifyAlign.cpp:262-263`); `GeneFull` sets none, so its
//!   region columns stay zero, and no feature rustar quantifies ever sets
//!   `intronic` or the antisense columns. A read whose barcode did not match
//!   the whitelist, or whose UMI was rejected, goes to `CBnotInPasslist` here.
//! - At the counting stage (`SoloReadFeature_inputRecords.cpp:100-150`) the
//!   rest are resolved: an exact barcode is `cbPerfect`; a single one-mismatch
//!   correction is `cbMMunique`, unless `oneExact` drops it because no read
//!   matched that barcode exactly; several one-mismatch candidates are
//!   `cbMMmultiple` when the posterior picks one. The dropped ones also go to
//!   `CBnotInPasslist`. A resolved read that reached the unique matrix is
//!   `countedU`; one that reached the multi-gene distribution
//!   (`--soloMultiMappers`) is `countedM`.
//!
//! Corrections that wait on the exact counts are held in [`CellReadStats`]
//! until [`CellReadStats::resolve`] is given them.
//!
//! # D24: row order
//!
//! STAR emits these rows by iterating a `std::unordered_map` keyed by
//! whitelist index, so the order is a hash-table walk, not a sort. This writes
//! them in reverse order of first appearance in the reads, which agrees with
//! STAR only on the first few rows of CellRanger's tiny 10x data; what sets
//! STAR's order is not settled (DIVERGENCE.md 3.2). The **values** match, only
//! which line they appear on differs. Read this file by barcode, not by
//! position; sorting by barcode is the only stable thing to do with it.
//!
//! Reads are folded in by the in-order consumer of the alignment loop, so
//! "first appearance" is input order regardless of thread count. A barcode
//! first reached through a deferred correction appears where that read was.

use std::collections::HashMap;
use std::fmt::Write as _;

use crate::solo::whitelist::CbCandidate;

/// STAR's `SoloReadFlagClass` bits, in its column order.
pub mod bit {
    pub const CB_MATCH: u16 = 1 << 0;
    pub const CB_PERFECT: u16 = 1 << 1;
    pub const CB_MM_UNIQUE: u16 = 1 << 2;
    pub const CB_MM_MULTIPLE: u16 = 1 << 3;
    pub const GENOME_U: u16 = 1 << 4;
    pub const GENOME_M: u16 = 1 << 5;
    pub const FEATURE_U: u16 = 1 << 6;
    pub const FEATURE_M: u16 = 1 << 7;
    pub const EXONIC: u16 = 1 << 8;
    pub const INTRONIC: u16 = 1 << 9;
    pub const EXONIC_AS: u16 = 1 << 10;
    pub const INTRONIC_AS: u16 = 1 << 11;
    pub const MITO: u16 = 1 << 12;
    pub const COUNTED_U: u16 = 1 << 13;
    pub const COUNTED_M: u16 = 1 << 14;
    /// Number of bits, i.e. counter columns.
    pub const N: usize = 15;
}

/// The column names, in bit order (STAR's `statNames`).
const STAT_NAMES: [&str; bit::N] = [
    "cbMatch",
    "cbPerfect",
    "cbMMunique",
    "cbMMmultiple",
    "genomeU",
    "genomeM",
    "featureU",
    "featureM",
    "exonic",
    "intronic",
    "exonicAS",
    "intronicAS",
    "mito",
    "countedU",
    "countedM",
];

/// How one read's cell barcode matched, as far as the read stage knows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadCb {
    /// No whitelist match, a rejected UMI, or no barcode at all (`cbMatch < 0`).
    None,
    /// Exact whitelist hit (`cbMatch` 0).
    Exact(u32),
    /// A single one-mismatch correction (`cbMatch` 1).
    Corrected(u32),
    /// Several one-mismatch candidates (`cbMatch` > 1), for the posterior.
    Multi(Vec<CbCandidate>),
}

/// One read's contribution: its barcode and, per feature, the read-stage flag
/// bits (`genomeU`/`genomeM`, `featureU`/`featureM`, region, `mito`).
#[derive(Debug, Clone)]
pub struct CellReadEntry {
    pub cb: ReadCb,
    pub flags: Vec<u16>,
}

/// One barcode's counters and the input position of its first read (D24).
#[derive(Debug, Clone)]
struct CellRow {
    first_seq: u64,
    counts: [u64; bit::N],
}

#[derive(Debug, Default, Clone)]
struct FeatureTable {
    cells: HashMap<u32, CellRow>,
    no_cb: [u64; bit::N],
}

/// A read whose barcode resolution waits on the whole run's exact counts.
#[derive(Debug, Clone)]
struct Pending {
    seq: u64,
    cb: ReadCb,
    flags: Vec<u16>,
}

/// The accumulator behind `CellReads.stats`: one table per quantified feature.
#[derive(Debug, Clone)]
pub struct CellReadStats {
    features: Vec<FeatureTable>,
    next_seq: u64,
    /// `oneExact`: single corrections are kept only for barcodes some read
    /// matched exactly, which is known only at the end.
    defer_corrected: bool,
    /// `--soloMultiMappers` other than `Unique`: gene-ambiguous reads are
    /// counted (`countedM`).
    count_multi_gene: bool,
    /// With a whitelist, unmatched reads are summed into `CBnotInPasslist`;
    /// without one STAR records nothing for them.
    whitelist: bool,
    pending: Vec<Pending>,
}

impl CellReadStats {
    pub fn new(
        n_features: usize,
        defer_corrected: bool,
        count_multi_gene: bool,
        whitelist: bool,
    ) -> Self {
        Self {
            features: vec![FeatureTable::default(); n_features],
            next_seq: 0,
            defer_corrected,
            count_multi_gene,
            whitelist,
            pending: Vec::new(),
        }
    }

    fn add(v: &mut [u64; bit::N], flags: u16) {
        for (i, c) in v.iter_mut().enumerate() {
            *c += u64::from((flags >> i) & 1);
        }
    }

    fn add_no_cb(&mut self, flags: &[u16]) {
        for (t, &f) in self.features.iter_mut().zip(flags) {
            Self::add(&mut t.no_cb, f | bit::CB_MATCH);
        }
    }

    /// `countsAdd(cb)`: the read resolved to whitelist cell `cb`, matched as
    /// `cb_bit`. A resolved read with a feature was counted.
    fn add_cell(&mut self, cb: u32, seq: u64, cb_bit: u16, flags: &[u16]) {
        for (t, &f) in self.features.iter_mut().zip(flags) {
            let mut f = f | bit::CB_MATCH | cb_bit;
            if f & bit::FEATURE_U != 0 {
                f |= bit::COUNTED_U;
            }
            if self.count_multi_gene && f & bit::FEATURE_M != 0 {
                f |= bit::COUNTED_M;
            }
            let row = t.cells.entry(cb).or_insert(CellRow {
                first_seq: seq,
                counts: [0; bit::N],
            });
            row.first_seq = row.first_seq.min(seq);
            Self::add(&mut row.counts, f);
        }
    }

    /// Fold in one read, in input order.
    pub fn fold(&mut self, entry: CellReadEntry) {
        debug_assert_eq!(entry.flags.len(), self.features.len());
        let seq = self.next_seq;
        self.next_seq += 1;
        match entry.cb {
            ReadCb::None => {
                if self.whitelist {
                    self.add_no_cb(&entry.flags);
                }
            }
            ReadCb::Exact(cb) => self.add_cell(cb, seq, bit::CB_PERFECT, &entry.flags),
            ReadCb::Corrected(cb) if !self.defer_corrected => {
                self.add_cell(cb, seq, bit::CB_MM_UNIQUE, &entry.flags);
            }
            cb @ (ReadCb::Corrected(_) | ReadCb::Multi(_)) => self.pending.push(Pending {
                seq,
                cb,
                flags: entry.flags,
            }),
        }
    }

    /// Resolve the held corrections once the exact counts are known:
    /// `has_exact(cb)` is `oneExact`'s test, `posterior(candidates)` the
    /// 1MM_multi choice. Reads neither accepts go to `CBnotInPasslist`.
    pub fn resolve(
        &mut self,
        has_exact: impl Fn(u32) -> bool,
        posterior: impl Fn(&[CbCandidate]) -> Option<u32>,
    ) {
        for p in std::mem::take(&mut self.pending) {
            let resolved = match &p.cb {
                ReadCb::Corrected(cb) => has_exact(*cb).then_some((*cb, bit::CB_MM_UNIQUE)),
                ReadCb::Multi(cands) => posterior(cands).map(|cb| (cb, bit::CB_MM_MULTIPLE)),
                ReadCb::None | ReadCb::Exact(_) => unreachable!("never held"),
            };
            match resolved {
                Some((cb, cb_bit)) => self.add_cell(cb, p.seq, cb_bit, &p.flags),
                None => self.add_no_cb(&p.flags),
            }
        }
    }

    /// Render feature `fi`'s file. `unique` gives a cell's
    /// `(nUMIunique, nGenesUnique)` and `multi` its `(nUMImulti, nGenesMulti)`;
    /// a cell absent from them prints zeros. `barcode_of` renders a whitelist
    /// index as its barcode.
    pub fn render<S: std::hash::BuildHasher>(
        &self,
        fi: usize,
        barcode_of: impl Fn(u32) -> String,
        unique: &HashMap<u32, (u64, u32), S>,
        multi: &HashMap<u32, (f64, u32), S>,
    ) -> String {
        let table = &self.features[fi];
        let mut s = String::from("CB");
        for name in STAT_NAMES {
            s.push('\t');
            s.push_str(name);
        }
        s.push_str("\tnUMIunique\tnGenesUnique\tnUMImulti\tnGenesMulti\n");
        s.push_str("CBnotInPasslist");
        for v in &table.no_cb {
            let _ = write!(s, "\t{v}");
        }
        s.push_str("\t0\t0\t0\t0\n");
        // Reverse first-appearance order (D24).
        let mut rows: Vec<(&u32, &CellRow)> = table.cells.iter().collect();
        rows.sort_unstable_by_key(|(_, r)| std::cmp::Reverse(r.first_seq));
        for (&cb, row) in rows {
            s.push_str(&barcode_of(cb));
            for v in &row.counts {
                let _ = write!(s, "\t{v}");
            }
            let (n_umi, n_gene) = unique.get(&cb).copied().unwrap_or((0, 0));
            let (m_umi, m_gene) = multi.get(&cb).copied().unwrap_or((0.0, 0));
            // `nUMIperCBmulti` is a double (`ostream << double`).
            let _ = writeln!(
                s,
                "\t{n_umi}\t{n_gene}\t{}\t{m_gene}",
                crate::solo::count::fmt_double(m_umi)
            );
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ONE: [u16; 1] = [bit::GENOME_U | bit::FEATURE_U | bit::EXONIC];

    fn entry(cb: ReadCb, flags: &[u16]) -> CellReadEntry {
        CellReadEntry {
            cb,
            flags: flags.to_vec(),
        }
    }

    fn stats(n_features: usize) -> CellReadStats {
        CellReadStats::new(n_features, false, false, true)
    }

    fn render(st: &CellReadStats, fi: usize) -> String {
        st.render(fi, |i| format!("CB{i}"), &HashMap::new(), &HashMap::new())
    }

    fn rows(text: &str) -> Vec<Vec<&str>> {
        text.lines()
            .skip(1)
            .map(|l| l.split('\t').collect())
            .collect()
    }

    /// The column of `name` in a split row.
    fn col(name: &str) -> usize {
        1 + STAT_NAMES.iter().position(|&n| n == name).unwrap()
    }

    /// Every column has a name, and every row has a value under each of them.
    #[test]
    fn each_row_fills_the_header() {
        let mut st = stats(1);
        st.fold(entry(ReadCb::Exact(0), &ONE));
        st.fold(entry(ReadCb::None, &[0]));
        let out = render(&st, 0);
        let header: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
        assert_eq!(header.len(), 20);
        assert_eq!(header[1], "cbMatch");
        assert_eq!(header[15], "countedM");
        let rows = rows(&out);
        assert_eq!(rows.len(), 2, "the passlist-miss row plus one cell");
        for row in rows {
            assert_eq!(row.len(), header.len());
        }
    }

    /// A read is counted once under `cbMatch` and once under each flag it sets.
    #[test]
    fn flags_accumulate_per_read() {
        let mut st = stats(1);
        for _ in 0..3 {
            st.fold(entry(ReadCb::Exact(7), &ONE));
        }
        st.fold(entry(ReadCb::Corrected(7), &[bit::GENOME_M]));
        let out = render(&st, 0);
        let row = &rows(&out)[1];
        assert_eq!(row[0], "CB7");
        assert_eq!(row[col("cbMatch")], "4");
        assert_eq!(row[col("cbPerfect")], "3");
        assert_eq!(row[col("cbMMunique")], "1");
        assert_eq!(row[col("genomeU")], "3");
        assert_eq!(row[col("genomeM")], "1");
        assert_eq!(row[col("countedU")], "3", "a resolved featureU read counts");
    }

    /// Each feature keeps its own table: the same read can be featureU for
    /// one and featureM for the other, as STAR's per-feature `readFlag` does.
    #[test]
    fn features_are_counted_separately() {
        let mut st = stats(2);
        st.fold(entry(
            ReadCb::Exact(1),
            &[
                bit::GENOME_U | bit::FEATURE_U | bit::EXONIC,
                bit::GENOME_U | bit::FEATURE_M,
            ],
        ));
        let gene = render(&st, 0);
        let full = render(&st, 1);
        let (g, f) = (&rows(&gene)[1], &rows(&full)[1]);
        assert_eq!(
            (g[col("featureU")], g[col("featureM")], g[col("exonic")]),
            ("1", "0", "1")
        );
        assert_eq!(
            (f[col("featureU")], f[col("featureM")], f[col("exonic")]),
            ("0", "1", "0")
        );
        assert_eq!(g[col("countedU")], "1");
        assert_eq!(f[col("countedU")], "0");
        assert_eq!(
            f[col("countedM")],
            "0",
            "multi-gene reads count only with --soloMultiMappers"
        );
    }

    /// With `--soloMultiMappers`, a resolved gene-ambiguous read is countedM.
    #[test]
    fn multi_gene_reads_count_with_multimappers() {
        let mut st = CellReadStats::new(1, false, true, true);
        st.fold(entry(ReadCb::Exact(1), &[bit::GENOME_M | bit::FEATURE_M]));
        st.fold(entry(ReadCb::None, &[bit::GENOME_M | bit::FEATURE_M]));
        let out = render(&st, 0);
        let rows = rows(&out);
        assert_eq!(rows[1][col("countedM")], "1");
        assert_eq!(rows[0][col("featureM")], "1");
        assert_eq!(
            rows[0][col("countedM")],
            "0",
            "an unresolved read is never counted"
        );
    }

    /// Reads whose barcode never resolved are summed rather than dropped, and
    /// keep their mapping and feature flags.
    #[test]
    fn unresolved_reads_land_in_the_passlist_miss_row() {
        let mut st = stats(1);
        st.fold(entry(ReadCb::Exact(0), &ONE));
        for _ in 0..5 {
            st.fold(entry(ReadCb::None, &ONE));
        }
        let out = render(&st, 0);
        let row = &rows(&out)[0];
        assert_eq!(row[0], "CBnotInPasslist");
        assert_eq!(row[col("cbMatch")], "5");
        assert_eq!(row[col("cbPerfect")], "0");
        assert_eq!(row[col("featureU")], "5");
        assert_eq!(row[col("countedU")], "0");
        assert_eq!(&row[16..], ["0", "0", "0", "0"], "no UMI columns for it");
    }

    /// Without a whitelist STAR records nothing for unmatched reads.
    #[test]
    fn no_whitelist_records_no_passlist_misses() {
        let mut st = CellReadStats::new(1, false, false, false);
        st.fold(entry(ReadCb::None, &ONE));
        assert_eq!(rows(&render(&st, 0))[0][col("cbMatch")], "0");
    }

    /// `oneExact`: a correction is held until the exact counts are known, then
    /// kept only for a barcode some read matched exactly.
    #[test]
    fn one_exact_corrections_wait_for_the_exact_counts() {
        let mut st = CellReadStats::new(1, true, false, true);
        st.fold(entry(ReadCb::Corrected(3), &ONE));
        st.fold(entry(ReadCb::Corrected(4), &ONE));
        assert_eq!(rows(&render(&st, 0)).len(), 1, "nothing resolved yet");
        st.resolve(|cb| cb == 3, |_| None);
        let out = render(&st, 0);
        let rows = rows(&out);
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows[0][col("cbMatch")],
            "1",
            "the drop goes to CBnotInPasslist"
        );
        assert_eq!(rows[0][col("featureU")], "1");
        assert_eq!(rows[1][0], "CB3");
        assert_eq!(rows[1][col("cbMMunique")], "1");
        assert_eq!(rows[1][col("countedU")], "1");
    }

    /// A 1MM_multi barcode is `cbMMmultiple` when the posterior picks a cell,
    /// `CBnotInPasslist` when it does not.
    #[test]
    fn multi_match_barcodes_follow_the_posterior() {
        let cand = |wl_index| CbCandidate {
            wl_index,
            mismatch_pos: 0,
            mismatch_qual: b'I',
        };
        let mut st = stats(1);
        st.fold(entry(ReadCb::Multi(vec![cand(5), cand(6)]), &ONE));
        st.fold(entry(ReadCb::Multi(vec![cand(8), cand(9)]), &ONE));
        st.resolve(|_| true, |c| (c[0].wl_index == 5).then_some(6));
        let out = render(&st, 0);
        let rows = rows(&out);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0][col("cbMatch")], "1");
        assert_eq!(rows[1][0], "CB6");
        assert_eq!(rows[1][col("cbMMmultiple")], "1");
        assert_eq!(rows[1][col("countedU")], "1");
    }

    /// D24: rows come out in reverse first-appearance order, a held read
    /// placed where it was read rather than where it was resolved.
    #[test]
    fn rows_are_emitted_in_reverse_first_appearance_order() {
        let mut st = CellReadStats::new(1, true, false, true);
        st.fold(entry(ReadCb::Exact(4), &ONE));
        st.fold(entry(ReadCb::Corrected(2), &ONE));
        st.fold(entry(ReadCb::Exact(1), &ONE));
        st.fold(entry(ReadCb::Exact(9), &ONE));
        st.fold(entry(ReadCb::Exact(1), &ONE)); // seen again: order unchanged
        st.resolve(|_| true, |_| None);
        let out = render(&st, 0);
        let cbs: Vec<&str> = rows(&out).iter().skip(1).map(|r| r[0]).collect();
        assert_eq!(cbs, ["CB9", "CB1", "CB2", "CB4"]);
    }

    /// The UMI and gene totals come from the matrices, not from the read
    /// counters, so a cell missing from them prints zeros; `nUMImulti` prints
    /// as the double it is.
    #[test]
    fn umi_and_gene_totals_come_from_the_matrices() {
        let mut st = stats(1);
        st.fold(entry(ReadCb::Exact(2), &ONE));
        st.fold(entry(ReadCb::Exact(3), &ONE));
        let unique = HashMap::from([(2u32, (17u64, 5u32))]);
        let multi = HashMap::from([(2u32, (1.5f64, 1u32))]);
        let out = st.render(0, |i| format!("CB{i}"), &unique, &multi);
        let rows = rows(&out);
        assert_eq!(rows[1][0], "CB3");
        assert_eq!(&rows[1][16..], ["0", "0", "0", "0"]);
        assert_eq!(rows[2][0], "CB2");
        assert_eq!(&rows[2][16..], ["17", "5", "1.5", "1"]);
    }
}
