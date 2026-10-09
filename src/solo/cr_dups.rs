//! CellRanger's per-barcode duplicate marking, ported from
//! `lib/rust/tx_annotation/src/mark_dups.rs` (cellranger-10.0.0).
//!
//! For one barcode's confidently mapped reads, each a `(UMI, gene)` pair:
//!
//! 1. every `(UMI, gene)` is moved to its 1-mismatch neighbour with the higher
//!    read count (a tie goes to the larger UMI); the move is not transitive;
//! 2. one read of each corrected pair is counted *first* and the rest *after*
//!    the low-support test, so a corrected UMI weighs as many reads as raw
//!    UMIs it absorbed, not as their reads;
//! 3. per UMI, a gene whose count is below the maximum, or any gene at a tied
//!    maximum, is low support: its reads are never counted;
//! 4. each surviving corrected `(UMI, gene)` is a molecule, represented by the
//!    read with the smallest `(UMI type, read name)` key among the reads of
//!    every raw UMI that was corrected into it; the other reads are duplicates.

use std::collections::{BTreeMap, HashMap, HashSet};

/// One read of a barcode.
#[derive(Debug, Clone, Copy)]
pub struct DupRead {
    pub read_index: u32,
    pub umi: u64,
    pub gene: u32,
    /// `(UMI type << 63) | read-name order`: transcriptomic reads sort first.
    pub key: u64,
}

/// What marking decided for one read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DupRes {
    pub read_index: u32,
    pub processed_umi: u64,
    pub low_support: bool,
    /// The molecule's representative: `UMI_COUNT`, and what the matrix counts.
    pub umi_count: bool,
}

/// A counted molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DupMolecule {
    pub gene: u32,
    pub umi: u64,
    pub reads: u32,
    /// CellRanger's `UmiType`: 1 transcriptomic, 0 not.
    pub utype: u8,
}

fn correct_umis(counts: &HashMap<(u64, u32), u64>, umi_len: usize) -> HashMap<(u64, u32), u64> {
    let mut out = HashMap::new();
    for (&(umi, gene), &orig) in counts {
        let (mut best_count, mut best_umi) = (orig, umi);
        for pos in 0..umi_len {
            let shift = 2 * (umi_len - 1 - pos);
            let cur = (umi >> shift) & 3;
            for b in 0..4u64 {
                if b == cur {
                    continue;
                }
                let test = (umi & !(3 << shift)) | (b << shift);
                let c = counts.get(&(test, gene)).copied().unwrap_or(0);
                if c > best_count || (c == best_count && test > best_umi) {
                    best_count = c;
                    best_umi = test;
                }
            }
        }
        if best_umi != umi {
            out.insert((umi, gene), best_umi);
        }
    }
    out
}

fn low_support(counts: &HashMap<(u64, u32), u64>) -> HashSet<(u64, u32)> {
    let mut by_umi: BTreeMap<u64, Vec<(u32, u64)>> = BTreeMap::new();
    for (&(umi, gene), &c) in counts {
        by_umi.entry(umi).or_default().push((gene, c));
    }
    let mut out = HashSet::new();
    for (umi, genes) in by_umi {
        let max = genes.iter().map(|g| g.1).max().unwrap_or(0);
        let tied = genes.iter().filter(|g| g.1 == max).count() >= 2;
        for (gene, c) in genes {
            if tied || c < max {
                out.insert((umi, gene));
            }
        }
    }
    out
}

/// Mark one barcode's reads. Returns a result per read (same order) and the
/// counted molecules.
pub fn mark_dups(reads: &[DupRead], umi_len: usize) -> (Vec<DupRes>, Vec<DupMolecule>) {
    let mut counts: HashMap<(u64, u32), u64> = HashMap::new();
    let mut min_key: HashMap<(u64, u32), u64> = HashMap::new();
    for r in reads {
        let k = (r.umi, r.gene);
        *counts.entry(k).or_insert(0) += 1;
        let e = min_key.entry(k).or_insert(r.key);
        *e = (*e).min(r.key);
    }
    let corrections = correct_umis(&counts, umi_len);

    // One read of each corrected pair moves before the low-support test.
    type Key = (u64, u32);
    let moves: Vec<(Key, Key, u64)> = corrections
        .iter()
        .map(|(&raw, &cu)| (raw, (cu, raw.1), counts[&raw]))
        .collect();
    for (raw, to, _) in &moves {
        *counts.get_mut(raw).unwrap() -= 1;
        *counts.get_mut(to).unwrap() += 1;
    }
    let low = low_support(&counts);
    // ... and the rest after it.
    for (raw, to, n) in &moves {
        *counts.get_mut(raw).unwrap() -= n - 1;
        *counts.get_mut(to).unwrap() += n - 1;
    }

    // The smallest key among the raw UMIs that end up at a corrected pair. A
    // raw UMI that is itself a correction target keeps its own key too.
    let mut min_raw: HashMap<(u64, u32), u64> = HashMap::new();
    for (&(raw, gene), &cu) in &corrections {
        if raw < cu || corrections.contains_key(&(cu, gene)) {
            let e = min_raw.entry((cu, gene)).or_insert(raw);
            *e = (*e).min(raw);
        }
    }
    let mut key_after = min_key.clone();
    for (&(cu, gene), &raw) in &min_raw {
        key_after.insert((cu, gene), min_key[&(raw, gene)]);
    }

    let mut res = Vec::with_capacity(reads.len());
    let mut mols: HashMap<(u64, u32), DupMolecule> = HashMap::new();
    for r in reads {
        let cu = corrections.get(&(r.umi, r.gene)).copied().unwrap_or(r.umi);
        let ck = (cu, r.gene);
        let low_support = low.contains(&ck);
        let is_min = key_after.get(&ck) == Some(&r.key);
        let umi_count = !low_support && is_min;
        if umi_count {
            mols.entry(ck).or_insert(DupMolecule {
                gene: r.gene,
                umi: cu,
                reads: counts[&ck] as u32,
                utype: u8::from(r.key >> 63 == 0),
            });
        }
        res.push(DupRes {
            read_index: r.read_index,
            processed_umi: cu,
            low_support,
            umi_count,
        });
    }
    let mut mols: Vec<DupMolecule> = mols.into_values().collect();
    mols.sort_unstable_by_key(|m| (m.gene, m.umi));
    (res, mols)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(i: u32, umi: u64, gene: u32, key: u64) -> DupRead {
        DupRead {
            read_index: i,
            umi,
            gene,
            key,
        }
    }

    #[test]
    fn a_lower_count_umi_moves_to_its_neighbour_and_the_first_name_represents_it() {
        // 4-base UMIs 0 and 1 differ in one base.
        let reads = [read(0, 1, 7, 30), read(1, 0, 7, 20), read(2, 0, 7, 10)];
        let (res, mols) = mark_dups(&reads, 4);
        assert_eq!(mols.len(), 1);
        assert_eq!(mols[0].umi, 0);
        assert_eq!(mols[0].reads, 3);
        assert_eq!(res[0].processed_umi, 0);
        // The smallest key, 10, is the representative.
        assert!(res[2].umi_count && !res[0].umi_count && !res[1].umi_count);
    }

    #[test]
    fn equal_counts_move_to_the_larger_umi() {
        let reads = [read(0, 0, 1, 1), read(1, 1, 1, 2)];
        let (res, mols) = mark_dups(&reads, 4);
        assert_eq!(mols.len(), 1);
        assert_eq!(res[0].processed_umi, 1);
    }

    #[test]
    fn a_umi_shared_by_two_genes_keeps_only_the_better_supported_one() {
        let reads = [read(0, 5, 1, 1), read(1, 5, 1, 2), read(2, 5, 2, 3)];
        let (res, mols) = mark_dups(&reads, 4);
        assert_eq!(mols.len(), 1);
        assert_eq!(mols[0].gene, 1);
        assert!(res[2].low_support && !res[0].low_support);
    }

    #[test]
    fn a_tie_between_genes_drops_the_umi() {
        let reads = [read(0, 5, 1, 1), read(1, 5, 2, 2)];
        let (res, mols) = mark_dups(&reads, 4);
        assert_eq!(mols.len(), 0);
        assert!(res.iter().all(|r| r.low_support));
    }

    #[test]
    fn a_non_transcriptomic_read_is_never_preferred() {
        let reads = [read(0, 5, 1, 1 << 63), read(1, 5, 1, 900)];
        let (res, mols) = mark_dups(&reads, 4);
        assert_eq!(mols[0].utype, 1);
        assert!(res[1].umi_count && !res[0].umi_count);
    }
}
