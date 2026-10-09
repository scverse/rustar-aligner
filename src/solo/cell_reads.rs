//! `--soloCellReadStats CB`: the per-cell-barcode read summary STARsolo writes
//! as `Solo.out/<feature>/CellReads.stats`.
//!
//! One row per cell barcode, fifteen counters describing what happened to the
//! reads carrying it — how the barcode matched, whether the read mapped
//! uniquely, whether it landed on a feature, where in the gene, and whether it
//! reached the matrix — plus the per-cell UMI and gene totals. The reads whose
//! barcode never resolved are not dropped; they are summed into a single
//! `CBnotInPasslist` row, so the columns account for every read rather than
//! only the ones that succeeded.
//!
//! # D24: row order
//!
//! STAR iterates a libc++ `std::unordered_map` to emit these rows, so the order
//! is a hash-table walk, not a sort. For the small maps this produces, libc++
//! chains new entries at the head of their bucket and walks buckets in order,
//! which comes out as the reverse of first appearance in read order. That is
//! what this reproduces.
//!
//! It is not reproducible in general: past the load factor libc++ rehashes, and
//! the order after a rehash depends on the bucket count, which depends on how
//! many distinct barcodes were seen. At that size the order diverges. The
//! **values never do** — only which line they appear on. A consumer that reads
//! this file by barcode rather than by position is unaffected either way, and
//! sorting by barcode is the only stable thing to do with it.
//!
//! Reads are folded in under a mutex held by `SoloContext`, so the order is the
//! order reads were processed in regardless of thread count.

use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;

/// What happened to one read, as the fourteen optional flags STAR tracks.
///
/// `cbMatch` is not here: every read that reaches the accumulator matched a
/// barcode well enough to be attributed somewhere, so it is always set.
///
/// Fourteen bools rather than a bitfield because they are written once per read
/// and read once per fold, and the names are the column names of the file.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone, Copy, Default)]
pub struct CellReadFlag {
    /// The barcode was an exact whitelist hit.
    pub cb_perfect: bool,
    /// Corrected via a single one-mismatch neighbour.
    pub cb_mm_unique: bool,
    /// Corrected via several one-mismatch neighbours, resolved by the posterior.
    pub cb_mm_multiple: bool,
    /// The read mapped to exactly one genomic locus.
    pub genome_u: bool,
    /// The read mapped to several genomic loci.
    pub genome_m: bool,
    /// It fell on exactly one feature (gene).
    pub feature_u: bool,
    /// It fell on several features.
    pub feature_m: bool,
    /// Exonic, on the annotated strand.
    pub exonic: bool,
    /// Intronic, on the annotated strand.
    pub intronic: bool,
    /// Exonic, antisense to the annotation.
    pub exonic_as: bool,
    /// Intronic, antisense to the annotation.
    pub intronic_as: bool,
    /// On a chromosome named by `--genomeChrSetMitochondrial`.
    pub mito: bool,
    /// Counted into the unique-gene matrix.
    pub counted_u: bool,
    /// Counted through the multi-gene distribution.
    pub counted_m: bool,
}

/// The accumulator behind `CellReads.stats`.
#[derive(Debug, Default, Clone)]
pub struct CellReadStats {
    /// Whitelist index to its fifteen counters.
    cells: BTreeMap<u32, [u64; 15]>,
    /// The single bucket for reads whose barcode did not resolve.
    no_cb: [u64; 15],
    /// Whitelist indices in first-appearance order; emitted reversed (D24).
    order: Vec<u32>,
    seen: HashSet<u32>,
}

impl CellReadStats {
    pub fn new() -> Self {
        Self::default()
    }

    fn fold(v: &mut [u64; 15], f: &CellReadFlag) {
        v[0] += 1; // cbMatch: every read that got this far
        for (i, set) in [
            f.cb_perfect,
            f.cb_mm_unique,
            f.cb_mm_multiple,
            f.genome_u,
            f.genome_m,
            f.feature_u,
            f.feature_m,
            f.exonic,
            f.intronic,
            f.exonic_as,
            f.intronic_as,
            f.mito,
            f.counted_u,
            f.counted_m,
        ]
        .into_iter()
        .enumerate()
        {
            if set {
                v[i + 1] += 1;
            }
        }
    }

    /// Record a read that resolved to whitelist cell `cb`.
    pub fn add_cell(&mut self, cb: u32, flag: &CellReadFlag) {
        if self.seen.insert(cb) {
            self.order.push(cb);
        }
        Self::fold(self.cells.entry(cb).or_insert([0; 15]), flag);
    }

    /// Record a read whose barcode did not resolve, or whose UMI was rejected.
    pub fn add_no_cb(&mut self, flag: &CellReadFlag) {
        Self::fold(&mut self.no_cb, flag);
    }

    /// Render the file. `umi_gene` gives each cell its final
    /// `(nUMIunique, nGenesUnique)`; a cell absent from it prints zeros.
    /// `barcode_of` renders a whitelist index as its barcode string.
    pub fn render(
        &self,
        barcode_of: impl Fn(u32) -> String,
        umi_gene: &BTreeMap<u32, (u32, u32)>,
    ) -> String {
        let mut s = String::from(
            "CB\tcbMatch\tcbPerfect\tcbMMunique\tcbMMmultiple\tgenomeU\tgenomeM\tfeatureU\t\
             featureM\texonic\tintronic\texonicAS\tintronicAS\tmito\tcountedU\tcountedM\t\
             nUMIunique\tnGenesUnique\tnUMImulti\tnGenesMulti\n",
        );
        s.push_str("CBnotInPasslist");
        for v in &self.no_cb {
            s.push('\t');
            s.push_str(&v.to_string());
        }
        s.push_str("\t0\t0\t0\t0\n");
        // Reverse first-appearance order — see D24 in the module docs.
        for &cb in self.order.iter().rev() {
            s.push_str(&barcode_of(cb));
            for v in &self.cells[&cb] {
                s.push('\t');
                s.push_str(&v.to_string());
            }
            let (n_umi, n_gene) = umi_gene.get(&cb).copied().unwrap_or((0, 0));
            // nUMImulti and nGenesMulti stay zero: multi-gene UMIs are not
            // collapsed into per-cell totals here.
            let _ = writeln!(s, "\t{n_umi}\t{n_gene}\t0\t0");
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flag_perfect_counted() -> CellReadFlag {
        CellReadFlag {
            cb_perfect: true,
            genome_u: true,
            feature_u: true,
            exonic: true,
            counted_u: true,
            ..Default::default()
        }
    }

    fn header_and_rows(text: &str) -> (Vec<&str>, Vec<&str>) {
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        (header, lines.collect())
    }

    /// Every column has a name, and every row has a value under each of them.
    #[test]
    fn each_row_fills_the_header() {
        let mut st = CellReadStats::new();
        st.add_cell(0, &flag_perfect_counted());
        st.add_no_cb(&CellReadFlag::default());
        let out = st.render(|i| format!("CB{i}"), &BTreeMap::new());
        let (header, rows) = header_and_rows(&out);
        assert_eq!(header.len(), 20);
        assert_eq!(rows.len(), 2, "the passlist-miss row plus one cell");
        for row in rows {
            assert_eq!(row.split('\t').count(), header.len());
        }
    }

    /// A read is counted once under `cbMatch` and once under each flag it sets,
    /// so the flag columns can exceed neither `cbMatch` nor each other's logic.
    #[test]
    fn flags_accumulate_per_read() {
        let mut st = CellReadStats::new();
        for _ in 0..3 {
            st.add_cell(7, &flag_perfect_counted());
        }
        st.add_cell(
            7,
            &CellReadFlag {
                cb_mm_unique: true,
                genome_m: true,
                ..Default::default()
            },
        );
        let out = st.render(|i| format!("CB{i}"), &BTreeMap::new());
        let row: Vec<&str> = out.lines().nth(2).unwrap().split('\t').collect();
        assert_eq!(row[0], "CB7");
        assert_eq!(row[1], "4", "cbMatch counts every read");
        assert_eq!(row[2], "3", "cbPerfect");
        assert_eq!(row[3], "1", "cbMMunique");
        assert_eq!(row[5], "3", "genomeU");
        assert_eq!(row[6], "1", "genomeM");
        assert_eq!(row[14], "3", "countedU");
    }

    /// Reads whose barcode never resolved are summed rather than dropped, so
    /// the file accounts for the whole input.
    #[test]
    fn unresolved_reads_land_in_the_passlist_miss_row() {
        let mut st = CellReadStats::new();
        st.add_cell(0, &flag_perfect_counted());
        for _ in 0..5 {
            st.add_no_cb(&CellReadFlag {
                genome_u: true,
                ..Default::default()
            });
        }
        let out = st.render(|i| format!("CB{i}"), &BTreeMap::new());
        let row: Vec<&str> = out.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[0], "CBnotInPasslist");
        assert_eq!(row[1], "5");
        assert_eq!(row[5], "5", "genomeU");
        assert_eq!(&row[16..], ["0", "0", "0", "0"], "no UMI columns for it");
    }

    /// D24: rows come out in reverse first-appearance order, which is what
    /// STAR's libc++ hash-map walk produces at these sizes.
    #[test]
    fn rows_are_emitted_in_reverse_first_appearance_order() {
        let mut st = CellReadStats::new();
        for cb in [4u32, 1, 9] {
            st.add_cell(cb, &flag_perfect_counted());
        }
        st.add_cell(1, &flag_perfect_counted()); // seen again: order unchanged
        let out = st.render(|i| format!("CB{i}"), &BTreeMap::new());
        let cbs: Vec<&str> = out
            .lines()
            .skip(2)
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(cbs, ["CB9", "CB1", "CB4"]);
    }

    /// The UMI and gene totals come from the final matrix, not from the read
    /// counters, so a cell missing from that map prints zeros rather than
    /// inheriting a read count.
    #[test]
    fn umi_and_gene_totals_come_from_the_matrix() {
        let mut st = CellReadStats::new();
        st.add_cell(2, &flag_perfect_counted());
        st.add_cell(3, &flag_perfect_counted());
        let mut umi_gene = BTreeMap::new();
        umi_gene.insert(2u32, (17u32, 5u32));
        let out = st.render(|i| format!("CB{i}"), &umi_gene);
        let rows: Vec<Vec<&str>> = out
            .lines()
            .skip(2)
            .map(|l| l.split('\t').collect())
            .collect();
        // Reverse order: CB3 first, then CB2.
        assert_eq!(rows[0][0], "CB3");
        assert_eq!(&rows[0][16..18], ["0", "0"]);
        assert_eq!(rows[1][0], "CB2");
        assert_eq!(&rows[1][16..18], ["17", "5"]);
    }
}
