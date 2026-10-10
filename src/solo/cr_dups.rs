//! Per-barcode UMI correction and duplicate marking.
//!
//! Independent implementation from a behavioural specification.
//!
//! All confidently, uniquely gene-assigned reads of one barcode are passed to
//! [`mark_dups`]. UMIs are corrected to a same-gene Hamming-1 neighbour with
//! better support, low-support (gene-ambiguous) UMIs are flagged, and one
//! representative read per molecule is chosen; the other reads are duplicates.

use std::collections::HashMap;

/// One read entering the stage.
#[derive(Debug, Clone, Copy)]
pub struct DupRead {
    pub read_index: u32,
    /// Raw UMI, 2 bits per base, first base most significant.
    pub umi: u64,
    pub gene: u32,
    /// Order key: bit 63 is 1 for non-transcriptomic reads, low bits are the
    /// read-name order.
    pub key: u64,
}

/// Per-read result, in input order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DupRes {
    pub read_index: u32,
    pub processed_umi: u64,
    pub low_support: bool,
    /// True when the read is the representative of its molecule.
    pub umi_count: bool,
}

/// One counted molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DupMolecule {
    pub gene: u32,
    pub umi: u64,
    pub reads: u32,
    /// 1 = transcriptomic representative, 0 = otherwise.
    pub utype: u8,
}

struct Pair {
    umi: u64,
    gene: u32,
    n: u32,
    kmin: u64,
    target: Option<usize>,
}

/// Correct UMIs, flag low-support reads and pick molecule representatives.
pub fn mark_dups(reads: &[DupRead], umi_len: usize) -> (Vec<DupRes>, Vec<DupMolecule>) {
    // Distinct raw pairs.
    let mut index: HashMap<(u64, u32), usize> = HashMap::new();
    let mut pairs: Vec<Pair> = Vec::new();
    let mut read_pair: Vec<usize> = Vec::with_capacity(reads.len());
    for r in reads {
        let id = *index.entry((r.umi, r.gene)).or_insert_with(|| {
            pairs.push(Pair {
                umi: r.umi,
                gene: r.gene,
                n: 0,
                kmin: u64::MAX,
                target: None,
            });
            pairs.len() - 1
        });
        pairs[id].n += 1;
        pairs[id].kmin = pairs[id].kmin.min(r.key);
        read_pair.push(id);
    }

    // R1: correction.
    for id in 0..pairs.len() {
        let (u, g, n) = (pairs[id].umi, pairs[id].gene, pairs[id].n);
        let (mut best_n, mut best_u, mut best_id) = (n, u, id);
        for pos in 0..umi_len {
            let shift = 2 * pos;
            let cur = (u >> shift) & 3;
            for alt in 0..4u64 {
                if alt == cur {
                    continue;
                }
                let v = (u & !(3u64 << shift)) | (alt << shift);
                if let Some(&vid) = index.get(&(v, g)) {
                    let vn = pairs[vid].n;
                    if vn > best_n || (vn == best_n && v > best_u) {
                        best_n = vn;
                        best_u = v;
                        best_id = vid;
                    }
                }
            }
        }
        if best_id != id {
            pairs[id].target = Some(best_id);
        }
    }

    // R2 step 1: one read per corrected pair moves.
    let mut m: Vec<i64> = pairs.iter().map(|p| i64::from(p.n)).collect();
    for (id, p) in pairs.iter().enumerate() {
        if let Some(t) = p.target {
            m[id] -= 1;
            m[t] += 1;
        }
    }

    // R3: low support on the intermediate counts, grouped by UMI.
    let mut by_umi: HashMap<u64, Vec<usize>> = HashMap::new();
    for (id, p) in pairs.iter().enumerate() {
        by_umi.entry(p.umi).or_default().push(id);
    }
    let mut low = vec![false; pairs.len()];
    for ids in by_umi.values() {
        let max = ids.iter().map(|&i| m[i]).max().unwrap_or(0);
        let n_max = ids.iter().filter(|&&i| m[i] == max).count();
        for &i in ids {
            low[i] = n_max >= 2 || m[i] < max;
        }
    }

    // R2 step 3: move the remaining reads.
    for (id, p) in pairs.iter().enumerate() {
        if let Some(t) = p.target {
            let rest = i64::from(p.n) - 1;
            m[id] -= rest;
            m[t] += rest;
        }
    }

    // R4: molecule keys. For each target, the smallest qualifying source UMI.
    let mut min_src: HashMap<usize, usize> = HashMap::new();
    for (id, p) in pairs.iter().enumerate() {
        if let Some(t) = p.target {
            let c = pairs[t].umi;
            if p.umi < c || pairs[t].target.is_some() {
                min_src
                    .entry(t)
                    .and_modify(|cur| {
                        if p.umi < pairs[*cur].umi {
                            *cur = id;
                        }
                    })
                    .or_insert(id);
            }
        }
    }
    let final_key = |id: usize| -> u64 {
        match min_src.get(&id) {
            Some(&s) => pairs[s].kmin,
            None => pairs[id].kmin,
        }
    };

    let mut res = Vec::with_capacity(reads.len());
    let mut mols: HashMap<usize, DupMolecule> = HashMap::new();
    for (r, &pid) in reads.iter().zip(&read_pair) {
        let cid = pairs[pid].target.unwrap_or(pid);
        let is_low = low[cid];
        let rep = !is_low && r.key == final_key(cid);
        if rep {
            mols.entry(cid).or_insert(DupMolecule {
                gene: pairs[cid].gene,
                umi: pairs[cid].umi,
                reads: m[cid].max(0) as u32,
                utype: u8::from(r.key >> 63 == 0),
            });
        }
        res.push(DupRes {
            read_index: r.read_index,
            processed_umi: pairs[cid].umi,
            low_support: is_low,
            umi_count: rep,
        });
    }
    let mut mols: Vec<DupMolecule> = mols.into_values().collect();
    mols.sort_by_key(|x| (x.gene, x.umi));
    (res, mols)
}

#[cfg(test)]
mod tests {
    use super::*;

    const L: usize = 4;

    fn rd(i: u32, umi: u64, gene: u32, key: u64) -> DupRead {
        DupRead {
            read_index: i,
            umi,
            gene,
            key,
        }
    }

    // 4-base UMIs: AAAC=1, AAAG=2, AACA=4, CAAA=64.
    #[test]
    fn single_read() {
        let (r, m) = mark_dups(&[rd(7, 5, 1, 3)], L);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].processed_umi, 5);
        assert!(!r[0].low_support && r[0].umi_count);
        assert_eq!(
            m,
            vec![DupMolecule {
                gene: 1,
                umi: 5,
                reads: 1,
                utype: 1
            }]
        );
    }

    #[test]
    fn singletons_smaller_moves_to_larger() {
        let (r, m) = mark_dups(&[rd(0, 1, 0, 10), rd(1, 2, 0, 20)], L);
        assert_eq!(r[0].processed_umi, 2);
        assert_eq!(r[1].processed_umi, 2);
        // S = {1} (1 < 2): key is that of UMI 1's read (10).
        assert!(r[0].umi_count && !r[1].umi_count);
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].reads, 2);
    }

    #[test]
    fn better_supported_neighbour_wins() {
        let reads = [rd(0, 1, 0, 5), rd(1, 2, 0, 9), rd(2, 2, 0, 8)];
        let (r, m) = mark_dups(&reads, L);
        assert_eq!(r[0].processed_umi, 2);
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].umi, 2);
        assert_eq!(m[0].reads, 3);
        // S={1}: key 5 replaces target's own (8); read 0 is representative.
        assert!(r[0].umi_count && !r[1].umi_count && !r[2].umi_count);
    }

    #[test]
    fn source_with_larger_umi_does_not_influence_key() {
        // target 1 (2 reads), source 2 (1 read, 2 > 1, target not corrected).
        let reads = [rd(0, 1, 0, 50), rd(1, 1, 0, 60), rd(2, 2, 0, 1)];
        let (r, m) = mark_dups(&reads, L);
        assert_eq!(r[2].processed_umi, 1);
        assert!(!r[2].umi_count);
        assert!(r[0].umi_count && !r[1].umi_count);
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].reads, 3);
    }

    #[test]
    fn not_transitive() {
        // A=1 (1 read) -> B=2 (2 reads) -> C=3? 2 and 3 differ at pos0 (AAAG/AAAT).
        // C=3 has 3 reads. 1 and 3 differ (AAAC vs AAAT) so A sees both B and C.
        // Use A=1(1), B=2(2), C=3(3): A picks best (3 reads) = C directly.
        // Make A see only B: A=4 (AACA), B=0? Use chain 4 -> 5 -> 7.
        // 4=AACA,5=AACC,7=AACT: all mutually 1 apart. Use positions instead:
        // A=0 (AAAA,1 read), B=1 (AAAC, 2 reads), C=5 (AACC, 3 reads).
        // 0~1 one apart, 1~5 one apart, 0~5 two apart.
        let mut reads = vec![rd(0, 0, 0, 1)];
        reads.extend([rd(1, 1, 0, 2), rd(2, 1, 0, 3)]);
        reads.extend([rd(3, 5, 0, 4), rd(4, 5, 0, 5), rd(5, 5, 0, 6)]);
        let (r, m) = mark_dups(&reads, L);
        assert_eq!(r[0].processed_umi, 1);
        assert_eq!(r[1].processed_umi, 5);
        assert_eq!(r[3].processed_umi, 5);
        // Intermediate: A:0, B:2-1+1=2, C:3+1=4 -> A (0) low support? Only one
        // gene per UMI, so nothing is low support.
        assert!(r.iter().all(|x| !x.low_support));
        // C: S={1} (1<5), key = kmin(B)=2. A's reads end at B, which is
        // itself corrected: B's S = {0} (B is a corrected pair) -> key 1.
        let c = r.iter().filter(|x| x.umi_count).count();
        assert_eq!(c, m.len());
        // C molecule has all of B's and C's reads except A's.
        assert!(m.iter().any(|x| x.umi == 5 && x.reads == 5));
        // A's read ends at B, whose final count is 1 (A's read only after
        // B's reads moved on to C).
        assert!(m.iter().any(|x| x.umi == 1 && x.reads == 1));
    }

    #[test]
    fn different_genes_never_correct() {
        let (r, _) = mark_dups(&[rd(0, 1, 0, 1), rd(1, 2, 1, 2)], L);
        assert_eq!(r[0].processed_umi, 1);
        assert_eq!(r[1].processed_umi, 2);
    }

    #[test]
    fn shared_umi_tie_is_low_support() {
        let (r, m) = mark_dups(&[rd(0, 9, 0, 1), rd(1, 9, 1, 2)], L);
        assert!(r[0].low_support && r[1].low_support);
        assert!(!r[0].umi_count && !r[1].umi_count);
        assert_eq!(m.len(), 0);
    }

    #[test]
    fn shared_umi_unequal_support() {
        let reads = [rd(0, 9, 0, 1), rd(1, 9, 0, 2), rd(2, 9, 1, 3)];
        let (r, m) = mark_dups(&reads, L);
        assert!(!r[0].low_support && !r[1].low_support && r[2].low_support);
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].gene, 0);
        assert_eq!(m[0].reads, 2);
        assert!(r[0].umi_count && !r[1].umi_count);
    }

    #[test]
    fn step_one_weighting_for_low_support() {
        // Gene0: UMI 1 x5 absorbing UMI 2 x5? Equal counts: smaller moves.
        // UMI 2 (5 reads, gene 0) -> corrected to nothing; UMI 1 (5) -> 2.
        // Intermediate: UMI1: 4, UMI2: 6. Shared UMI 1 with gene 1 x5 reads:
        // UMI 1 gene0 = 4 < gene1 = 5 -> gene0 low support (read weight 1).
        let mut reads = Vec::new();
        for i in 0..5 {
            reads.push(rd(i, 1, 0, u64::from(i)));
            reads.push(rd(10 + i, 2, 0, 100 + u64::from(i)));
            reads.push(rd(20 + i, 1, 1, 200 + u64::from(i)));
        }
        let (r, _) = mark_dups(&reads, L);
        // UMI 1 gene 0 corrected to UMI 2; (2,g0) intermediate 6, no clash.
        let g0_u1 = r.iter().find(|x| x.read_index == 0).unwrap();
        assert_eq!(g0_u1.processed_umi, 2);
        assert!(!g0_u1.low_support);
        // (1,g1) is alone with UMI 1 at intermediate: (1,g0)=4, (1,g1)=5.
        // g1 is max (unique): not low support; g0 pair low support, but its
        // reads are judged by the corrected pair (2,g0), not low.
        let g1 = r.iter().find(|x| x.read_index == 20).unwrap();
        assert!(!g1.low_support);
    }

    #[test]
    fn low_support_of_corrected_pair_applies() {
        // (3,g0) 1 read, (3,g1) 1 read: tied -> low support both.
        // (2,g0) 1 read: neighbour of (3,g0)=AAAT? 2=AAAG,3=AAAT.
        let reads = [rd(0, 3, 0, 1), rd(1, 3, 1, 2), rd(2, 2, 0, 3)];
        let (r, m) = mark_dups(&reads, L);
        // (2,g0) -> (3,g0) (equal count, larger umi). Intermediate: (2,g0)=0,
        // (3,g0)=2 ; (3,g1)=1. UMI3 max=2 unique -> g1 low support only.
        assert_eq!(r[2].processed_umi, 3);
        assert!(!r[0].low_support && r[1].low_support && !r[2].low_support);
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].reads, 2);
    }

    #[test]
    fn nontranscriptomic_loses_representative() {
        let top = 1u64 << 63;
        let reads = [rd(0, 5, 0, top | 1), rd(1, 5, 0, 7)];
        let (r, m) = mark_dups(&reads, L);
        assert!(!r[0].umi_count && r[1].umi_count);
        assert_eq!(m[0].utype, 1);
        let (r, m) = mark_dups(&[rd(0, 5, 0, top | 1)], L);
        assert!(r[0].umi_count);
        assert_eq!(m[0].utype, 0);
    }

    #[test]
    fn equal_keys_all_representative_and_molecule_unique() {
        let (r, m) = mark_dups(&[rd(0, 5, 0, 4), rd(1, 5, 0, 4)], L);
        assert!(r[0].umi_count && r[1].umi_count);
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].reads, 2);
    }

    #[test]
    fn molecule_without_representative() {
        // B=1 (2 reads) absorbs A=0 (1 read, key 9); B's own key 1 is replaced
        // by A's key (9) since 0<1. Then C=5 (3 reads) absorbs B (B corrected):
        // C's S = {1} -> key kmin(B)=1... build instead: S replaced key
        // belongs to reads that moved elsewhere.
        // A=0 (1 read key 9) -> B=1 (2 reads keys 1,2) -> C=5 (3 reads).
        let reads = [
            rd(0, 0, 0, 9),
            rd(1, 1, 0, 1),
            rd(2, 1, 0, 2),
            rd(3, 5, 0, 3),
            rd(4, 5, 0, 4),
            rd(5, 5, 0, 5),
        ];
        let (r, m) = mark_dups(&reads, L);
        // B's key becomes 9 (from A) but B's reads (keys 1,2) end at C whose
        // key is kmin(B)=1. A's read ends at B, with key 9 == key(B): rep of
        // B, but B is not a counted pair for A... A's corrected pair is B.
        assert!(r[0].umi_count);
        assert!(r[1].umi_count && !r[2].umi_count);
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn output_order_and_sorting() {
        let reads = [rd(3, 9, 2, 1), rd(1, 8, 1, 1), rd(2, 6, 1, 1)];
        let (r, m) = mark_dups(&reads, L);
        assert_eq!(
            r.iter().map(|x| x.read_index).collect::<Vec<_>>(),
            vec![3, 1, 2]
        );
        let keys: Vec<_> = m.iter().map(|x| (x.gene, x.umi)).collect();
        assert_eq!(keys, vec![(1, 6), (1, 8), (2, 9)]);
    }

    /// Oracle check against a TSV extracted from a real BAM (name, flag, cb,
    /// UR, gene, UB, xf, transcriptomic 0/1). Set `DUPS_TSV`; optional `DUPS_OUT` receives the
    /// per-(barcode, gene) molecule counts.
    #[test]
    #[ignore = "needs DUPS_TSV extracted from a real BAM"]
    fn oracle_real_bam() {
        use std::io::Write;
        let path = std::env::var("DUPS_TSV").expect("DUPS_TSV");
        let text = std::fs::read_to_string(path).unwrap();
        let mut by_bc: HashMap<&str, Vec<Vec<&str>>> = HashMap::new();
        for l in text.lines() {
            let f: Vec<&str> = l.split(' ').collect();
            // Only confidently mapped reads (xf bit 1) enter the stage.
            if f[6].parse::<u32>().unwrap_or(0) & 1 == 0 {
                continue;
            }
            by_bc.entry(f[2]).or_default().push(f);
        }
        let enc = |s: &str| {
            s.bytes().fold(0u64, |a, b| {
                a * 4
                    + match b {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    }
            })
        };
        let dec = |mut v: u64, l: usize| {
            let mut o = vec![b'A'; l];
            for i in (0..l).rev() {
                o[i] = b"ACGT"[(v & 3) as usize];
                v >>= 2;
            }
            String::from_utf8(o).unwrap()
        };
        let mut out = std::env::var("DUPS_OUT")
            .ok()
            .map(|p| std::fs::File::create(p).unwrap());
        let (mut nreads, mut bad_ub, mut bad_ls, mut bad_c, mut bad_dup) = (0, 0, 0, 0, 0);
        let mut bad_bc = std::collections::BTreeSet::new();
        for (bc, rows) in &by_bc {
            let l = rows[0][3].len();
            let mut names: Vec<&str> = rows.iter().map(|f| f[0]).collect();
            names.sort_unstable();
            names.dedup();
            let mut genes: Vec<&str> = rows.iter().map(|f| f[4]).collect();
            genes.sort_unstable();
            genes.dedup();
            let reads: Vec<DupRead> = rows
                .iter()
                .enumerate()
                .map(|(i, f)| DupRead {
                    read_index: i as u32,
                    umi: enc(f[3]),
                    gene: genes.binary_search(&f[4]).unwrap() as u32,
                    // Type bit: 1 unless the read has a transcript-level (ENST) TX tag.
                    key: names.binary_search(&f[0]).unwrap() as u64
                        | (u64::from(f.get(7) != Some(&"1")) << 63),
                })
                .collect();
            let (res, mols) = mark_dups(&reads, l);
            for (r, f) in res.iter().zip(rows) {
                nreads += 1;
                let xf: u32 = f[6].parse().unwrap_or(0);
                let flag: u32 = f[1].parse().unwrap();
                let mut bad = false;
                if dec(r.processed_umi, l) != f[5] {
                    bad_ub += 1;
                    bad = true;
                }
                if r.low_support != (xf & 2 != 0) {
                    bad_ls += 1;
                    bad = true;
                }
                if r.umi_count != (xf & 8 != 0) {
                    bad_c += 1;
                    bad = true;
                }
                if (!r.low_support && !r.umi_count) != (flag & 1024 != 0) {
                    bad_dup += 1;
                    bad = true;
                }
                if bad {
                    bad_bc.insert(*bc);
                }
            }
            if let Some(o) = out.as_mut() {
                let mut c: HashMap<u32, u32> = HashMap::new();
                for m in &mols {
                    *c.entry(m.gene).or_default() += 1;
                }
                for (g, n) in c {
                    writeln!(o, "{bc}\t{}\t{n}", genes[g as usize]).unwrap();
                }
            }
        }
        eprintln!(
            "barcodes {} reads {nreads} bad_ub {bad_ub} bad_lowsupport {bad_ls} bad_count {bad_c} bad_dupflag {bad_dup} badbc {}",
            by_bc.len(),
            bad_bc.len()
        );
        eprintln!(
            "first bad barcodes: {:?}",
            bad_bc.iter().take(5).collect::<Vec<_>>()
        );
        assert_eq!(bad_bc.len(), 0);
    }
}
