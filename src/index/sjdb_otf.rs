//! On-the-fly junction insertion into a loaded index: a port of STAR's
//! `sjdbInsertJunctions.cpp` + `sjdbBuildIndex.cpp`.
//!
//! STAR inserts junctions into an already-built genome in two situations: an
//! annotation (`--sjdbGTFfile` / `--sjdbFileChrStartEnd`) given at the mapping
//! step, and the junctions of the first pass in `--twopassMode Basic`. Both
//! rebuild the Gsj flanking-sequence buffer for the union of the old and new
//! junctions, keep every existing suffix-array entry in place (shifting the
//! positions that moved), binary-search the suffixes of the *new* junctions
//! into the suffix array, and regenerate the SA index. After that the new
//! junctions behave exactly like junctions built into the genome: seeds can
//! cross them in the Gsj buffer and the stitcher treats them as sjdb.

use std::cmp::Ordering;
use std::collections::HashMap;

use rayon::prelude::*;

use crate::error::Error;
use crate::genome::{Genome, GenomeSeq};
use crate::index::GenomeIndex;
use crate::index::packed_array::PackedArray;
use crate::index::sa_index::SaIndex;
use crate::index::suffix_array::SuffixArray;
use crate::junction::SpliceJunctionDb;
use crate::junction::sjdb_insert::{self, PreparedJunction};

/// STAR's `GENOME_spacingChar`.
const SPACER: u8 = 5;

/// STAR's `-2llu` from suffixArraySearch1: inserted after every existing entry.
const BEYOND_END: usize = usize::MAX - 1;

/// STAR's in-place SAindex update after inserting suffixes into the SA
/// (`sjdbBuildIndex.cpp:197-268`), a literal port. `inserts` are the inserted
/// suffixes as (insertion point in the old SA, offset in `gsj2`), in final
/// order; every recorded SA index moves up by the number of suffixes inserted
/// before it, and prefixes that only the inserted suffixes carry become present.
fn update_sa_index(sai: &mut SaIndex, inserts: &[(usize, usize)], gsj2: &[u8]) {
    let gb = sai.gstrand_bit;
    let n_mark = 1u64 << (gb + 1);
    let absent = 1u64 << (gb + 2);
    let value_mask = !(n_mark | absent);
    let start = sai.genome_sa_index_start.clone();
    let n_ind = inserts.len();
    // funCalcSAi: prefix code of the first iL+1 bases, minus its partial value
    // when a non-ACGT base comes first.
    let calc = |off: usize, il: usize| -> i64 {
        let mut ind1: i64 = 0;
        for k in 0..=il {
            let g = gsj2.get(off + k).copied().unwrap_or(5);
            if g > 3 {
                return -ind1;
            }
            ind1 = (ind1 << 2) + i64::from(g);
        }
        ind1
    };
    // indArray[2*iSJ]-1 and the -999 sentinel after the last entry.
    let ins_m1 = |i: usize| -> u64 {
        if i < n_ind {
            (inserts[i].0 as u64).wrapping_sub(1)
        } else {
            (-1000i64) as u64
        }
    };
    let data = &mut sai.data;
    for il in 0..sai.nbases as usize {
        let mut isj = 0usize;
        let mut ind0 = start[il] as i64 - 1;
        for ii in start[il]..start[il + 1] {
            let rel = (ii - start[il]) as i64;
            let isa1 = data.read(ii as usize);
            let isa2 = isa1 & value_mask;
            if isj < n_ind && isa1 & absent != 0 {
                let isj1 = isj;
                let mut ind1 = calc(inserts[isj].1, il);
                while ind1 < rel && ins_m1(isj) < isa2 {
                    isj += 1;
                    ind1 = if isj < n_ind {
                        calc(inserts[isj].1, il)
                    } else {
                        i64::MIN
                    };
                }
                if ind1 == rel {
                    let v = ins_m1(isj).wrapping_add(isj as u64 + 1);
                    data.write(ii as usize, v);
                    for ii0 in (ind0 + 1) as u64..ii {
                        data.write(ii0 as usize, v | absent);
                    }
                    isj += 1;
                    ind0 = ii as i64;
                } else {
                    isj = isj1;
                }
            } else {
                while isj < n_ind && ins_m1(isj).wrapping_add(1) < isa2 {
                    isj += 1;
                }
                while isj < n_ind && ins_m1(isj).wrapping_add(1) == isa2 {
                    if calc(inserts[isj].1, il) >= rel {
                        break;
                    }
                    isj += 1;
                }
                data.write(ii as usize, isa1 + isj as u64);
                for ii0 in (ind0 + 1) as u64..ii {
                    data.write(ii0 as usize, (isa2 + isj as u64) | absent);
                }
                ind0 = ii as i64;
            }
        }
    }
    // Inserted suffixes with a non-ACGT base inside the prefix mark the last
    // present entry before them (`:270-290`).
    let nbases = sai.nbases as usize;
    for &(_, off) in inserts {
        let mut ind1: i64 = 0;
        for il in 0..nbases {
            let g = gsj2.get(off + il).copied().unwrap_or(5);
            ind1 <<= 2;
            if g > 3 {
                for &level_start in &start[il..nbases] {
                    ind1 += 3;
                    let mut ind2 = level_start as i64 + ind1;
                    while ind2 >= 0 && data.read(ind2 as usize) & absent != 0 {
                        ind2 -= 1;
                    }
                    let i = ind2.max(0) as usize;
                    data.write(i, data.read(i) | n_mark);
                    ind1 <<= 2;
                }
                break;
            }
            ind1 += i64::from(g);
        }
    }
}

/// Insert `added` junctions into `index` (STAR's `sjdbInsertJunctions`).
///
/// `added` must already be prepared (`sjdbPrepare`d); they are merged with the
/// index's own junctions and deduplicated as STAR does. `sjdb_overhang` is used
/// only when the index has no junctions yet; otherwise the index's overhang
/// is kept, since its Gsj blocks are laid out with it. The junction database is
/// rebuilt from the merged list, so every inserted junction is annotated.
///
/// Returns the number of junctions that were not already in the index.
pub fn insert_junctions(
    index: &mut GenomeIndex,
    added: Vec<PreparedJunction>,
    sjdb_overhang: u32,
) -> Result<usize, Error> {
    let overhang = if index.sjdb_overhang > 0 {
        index.sjdb_overhang
    } else {
        sjdb_overhang
    };
    let old_junctions = std::mem::take(&mut index.prepared_junctions);
    let mut union = old_junctions.clone();
    union.extend(added);
    let all = sjdb_insert::sort_and_dedup(union);

    // The Gsj block of a junction is determined by its original coordinates,
    // so an old block keeps its content and only moves to its new index.
    let old_by_coords: HashMap<(u64, u64), usize> = old_junctions
        .iter()
        .enumerate()
        .map(|(i, j)| ((j.original_start(), j.original_end()), i))
        .collect();
    let mut old_to_new: Vec<Option<usize>> = vec![None; old_junctions.len()];
    let mut is_new = vec![true; all.len()];
    for (inew, j) in all.iter().enumerate() {
        if let Some(&iold) = old_by_coords.get(&(j.original_start(), j.original_end())) {
            old_to_new[iold] = Some(inew);
            is_new[inew] = false;
        }
    }
    let n_new = is_new.iter().filter(|&&n| n).count();
    // A dropped old block cannot be shifted anywhere; its suffixes are removed.
    let dropped_old = old_to_new.iter().any(Option::is_none);

    if n_new == 0 && !dropped_old {
        index.junction_db = SpliceJunctionDb::from_prepared(&all);
        index.prepared_junctions = all;
        return Ok(0);
    }

    let genome = &index.genome;
    let n_real = genome.n_genome_real;
    let old_n = genome.n_genome;
    let sjdb_length = (2 * overhang + 1) as u64;
    let n_gsj = all.len() as u64 * sjdb_length;
    let new_n = n_real + n_gsj;

    let gstrand_bit = index.suffix_array.gstrand_bit;
    if SuffixArray::calculate_gstrand_bit(new_n) > gstrand_bit {
        return Err(Error::Index(
            "EXITING because of FATAL ERROR: cannot insert junctions on the fly because of \
             strand GstrandBit problem"
                .to_string(),
        ));
    }
    let n2bit = 1u64 << gstrand_bit;
    let strand_mask = n2bit - 1;

    // New genome text: [real | Gsj] + reverse complement.
    let old_view = genome.sequence.view();
    let real: Vec<u8> = (0..n_real as usize).map(|i| old_view.base(i)).collect();
    let real_genome = Genome {
        transform_blocks: genome.transform_blocks.clone(),
        sequence: GenomeSeq::Owned(real),
        n_genome: n_real,
        n_genome_real: n_real,
        n_chr_real: genome.n_chr_real,
        chr_name: genome.chr_name.clone(),
        chr_length: genome.chr_length.clone(),
        chr_start: genome.chr_start.clone(),
    };
    let gsj = sjdb_insert::build_gsj(&all, &real_genome, n_real, overhang)?;
    // STAR's working `Gsj` array: forward blocks, their reverse complement,
    // then a terminating spacer (`sjdbBuildIndex.cpp:31-40`).
    let mut gsj2 = Vec::with_capacity(2 * gsj.len() + 1);
    gsj2.extend_from_slice(&gsj);
    gsj2.extend(gsj.iter().rev().map(|&b| if b < 4 { 3 - b } else { b }));
    gsj2.push(SPACER);

    // Suffixes of the new junctions, both strands, skipping those that start
    // with a non-ACGT base (`sjdbBuildIndex.cpp:62-80`).
    let n_jn = all.len();
    let l = sjdb_length as usize;
    let starts: Vec<usize> = (0..2 * n_jn)
        .filter(|&b| is_new[if b < n_jn { b } else { 2 * n_jn - 1 - b }])
        .flat_map(|b| (b * l..(b + 1) * l).filter(|&off| gsj2[off] < 4))
        .collect();

    // Insertion point of each suffix in the old SA: the first entry that is
    // greater (`suffixArraySearch1` with gInsert = -1).
    let sa = &index.suffix_array;
    let n_sa = sa.len();
    let old_text = |abs: usize| -> u8 {
        if abs < 2 * old_n as usize {
            old_view.base(abs)
        } else {
            SPACER
        }
    };
    let query_vs_sa = |off: usize, sa_value: u64| -> Ordering {
        let rev = sa_value & n2bit != 0;
        let pos = (sa_value & strand_mask) as usize;
        let abs = if rev { old_n as usize + pos } else { pos };
        let mut ii = 0;
        loop {
            let s = gsj2[off + ii];
            let g = old_text(abs + ii);
            if s != g {
                return s.cmp(&g);
            }
            if s == SPACER {
                // compareRefEnds(strR=true, gInsert=-1): the junction suffix
                // sorts after an equal forward-genome suffix and before an
                // equal reverse one.
                return if rev {
                    Ordering::Less
                } else {
                    Ordering::Greater
                };
            }
            ii += 1;
        }
    };
    let mut inserts: Vec<(usize, usize)> = starts
        .par_iter()
        .map(|&off| {
            let ins = partition_point(n_sa, |i| query_vs_sa(off, sa.get(i)) == Ordering::Greater);
            // suffixArraySearch1 returns -2 for a suffix past the last entry.
            (if ins == n_sa { BEYOND_END } else { ins }, off)
        })
        .collect();
    // funCompareUintAndSuffixes: by insertion point, then by suffix up to the
    // spacer, then by position in the Gsj array.
    inserts.par_sort_unstable_by(|a, b| {
        a.0.cmp(&b.0).then_with(|| {
            let mut ig = 0;
            loop {
                let (x, y) = (gsj2[a.1 + ig], gsj2[b.1 + ig]);
                if x != y {
                    return x.cmp(&y);
                }
                if x == SPACER {
                    return a.1.cmp(&b.1);
                }
                ig += 1;
            }
        })
    });

    // Merge: old entries in order, shifted; new ones before the entry they
    // were found to precede.
    let n_gsj_u = n_gsj as usize;
    let encode_new = |off: usize| -> u64 {
        if off < n_gsj_u {
            n_real + off as u64
        } else {
            (off - n_gsj_u) as u64 | n2bit
        }
    };
    let shift_block = |fwd_pos: u64| -> Option<u64> {
        let rel = fwd_pos - n_real;
        let sj_old = (rel / sjdb_length) as usize;
        old_to_new[sj_old].map(|sj_new| n_real + sj_new as u64 * sjdb_length + rel % sjdb_length)
    };
    let shift_entry = |v: u64| -> Option<u64> {
        let pos = v & strand_mask;
        if v & n2bit != 0 {
            // Reverse strand: offset into the RC half, whose start moves by
            // the Gsj growth; RC Gsj suffixes are mirrored back to their block.
            let fwd_base = old_n - 1 - pos;
            if fwd_base >= n_real {
                shift_block(fwd_base).map(|f| (new_n - 1 - f) | n2bit)
            } else {
                Some((pos + new_n - old_n) | n2bit)
            }
        } else if pos >= n_real {
            shift_block(pos)
        } else {
            Some(v)
        }
    };
    // Entries of dropped old blocks disappear; count them so the new SA can be
    // allocated at its exact size (a direct write, no intermediate vector).
    let n_dropped = if dropped_old {
        (0..n_sa)
            .into_par_iter()
            .filter(|&i| shift_entry(sa.get(i)).is_none())
            .count()
    } else {
        0
    };
    let n_out = n_sa - n_dropped + inserts.len();
    let word = gstrand_bit + 1;
    let data = if dropped_old {
        let mut data = PackedArray::new(word, n_out);
        let mut w = 0usize;
        let mut next = inserts.iter().peekable();
        for isa in 0..n_sa {
            while let Some(&&(ins, off)) = next.peek() {
                if ins != isa {
                    break;
                }
                data.write(w, encode_new(off));
                w += 1;
                next.next();
            }
            if let Some(value) = shift_entry(sa.get(isa)) {
                data.write(w, value);
                w += 1;
            }
        }
        for &(_, off) in next {
            data.write(w, encode_new(off));
            w += 1;
        }
        data
    } else {
        // Every old entry survives, so output position o is reachable by a
        // binary search, and blocks of whole bytes (a multiple of 8 entries)
        // are filled in parallel with a bit writer.
        let ins_at = |j: usize| inserts[j].0.min(n_sa);
        let before = |isa: usize| isa + inserts.partition_point(|&(ins, _)| ins.min(n_sa) < isa);
        let start_state = |o: usize| -> (usize, usize) {
            // Largest isa whose preceding outputs number at most o.
            let (mut lo, mut hi) = (0usize, n_sa);
            while lo < hi {
                let mid = lo + (hi - lo).div_ceil(2);
                if before(mid) <= o {
                    lo = mid;
                } else {
                    hi = mid - 1;
                }
            }
            let j = inserts.partition_point(|&(ins, _)| ins.min(n_sa) < lo) + (o - before(lo));
            (lo, j)
        };
        const BLOCK: usize = 1 << 22; // entries per parallel block, multiple of 8
        let mut bytes = vec![0u8; PackedArray::data_byte_len_for(word, n_out)];
        let block_bytes = BLOCK * word as usize / 8;
        bytes
            .par_chunks_mut(block_bytes)
            .enumerate()
            .for_each(|(bi, chunk)| {
                let o0 = bi * BLOCK;
                if o0 >= n_out {
                    return;
                }
                let o1 = (o0 + BLOCK).min(n_out);
                let (mut isa, mut j) = start_state(o0);
                let (mut acc, mut nbits, mut pos) = (0u128, 0u32, 0usize);
                for _ in o0..o1 {
                    let v = if j < inserts.len() && ins_at(j) == isa {
                        j += 1;
                        encode_new(inserts[j - 1].1)
                    } else {
                        isa += 1;
                        shift_entry(sa.get(isa - 1)).expect("no dropped blocks")
                    };
                    acc |= u128::from(v) << nbits;
                    nbits += word;
                    while nbits >= 8 {
                        chunk[pos] = acc as u8;
                        acc >>= 8;
                        nbits -= 8;
                        pos += 1;
                    }
                }
                if nbits > 0 {
                    chunk[pos] = acc as u8;
                }
            });
        PackedArray::from_bytes(word, n_out, bytes)
    };
    let new_sa = SuffixArray {
        data,
        gstrand_bit,
        gstrand_mask: strand_mask,
    };

    let mut new_genome = real_genome;
    new_genome.append_sjdb(&gsj);
    let mut new_sai = std::mem::replace(
        &mut index.sa_index,
        SaIndex {
            nbases: 0,
            genome_sa_index_start: Vec::new(),
            data: PackedArray::new(1, 0),
            word_length: 1,
            gstrand_bit,
        },
    );
    // A loaded SAindex is memory-mapped: take an owned copy to update.
    new_sai.data = PackedArray::from_bytes(
        new_sai.data.word_length(),
        new_sai.data.len(),
        new_sai.data.data().to_vec(),
    );
    update_sa_index(&mut new_sai, &inserts, &gsj2);

    log::info!(
        "Inserted {} new junctions on the fly ({} total, {} SA entries)",
        n_new,
        all.len(),
        new_sa.len()
    );
    index.genome = new_genome;
    index.suffix_array = new_sa;
    index.sa_index = new_sai;
    index.junction_db = SpliceJunctionDb::from_prepared(&all);
    index.prepared_junctions = all;
    index.sjdb_overhang = overhang;
    Ok(n_new)
}

impl GenomeIndex {
    /// `sjdbPrepare` raw `(chr_idx, start, end, strand)` junctions against the
    /// real genome (motif, repeat shifts, shifted coordinates), from a source
    /// of the given STAR priority (see [`PreparedJunction::priority`]).
    pub fn prepare_raw_junctions(
        &self,
        raw: &[(usize, u64, u64, u8)],
        priority: u8,
    ) -> Vec<PreparedJunction> {
        let n_real = self.genome.n_genome_real;
        raw.iter()
            .map(|&(chr_idx, start, end, strand)| PreparedJunction {
                priority,
                ..sjdb_insert::prepare_junction(chr_idx, start, end, strand, &self.genome, n_real)
            })
            .collect()
    }
}

/// First index in `0..n` for which `pred` is false (`pred` must be monotone).
fn partition_point(n: usize, pred: impl Fn(usize) -> bool) -> usize {
    let (mut lo, mut hi) = (0usize, n);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if pred(mid) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;

    /// Differential check against STAR: inserting a GTF's junctions on the fly
    /// into an unannotated index must give the genome and suffix array that
    /// STAR's `genomeGenerate --sjdbGTFfile` wrote (STAR builds that index by
    /// the same insertion). Needs two STAR indices of one genome; run with
    /// `RUSTAR_OTF_BASE=<plain index> RUSTAR_OTF_REF=<annotated index>
    /// RUSTAR_OTF_GTF=<gtf> cargo test --release -- --ignored otf_matches_star`.
    #[test]
    #[ignore = "needs local STAR indices"]
    fn otf_matches_star_genome_generate() {
        let base = std::env::var("RUSTAR_OTF_BASE").unwrap();
        let reference = std::env::var("RUSTAR_OTF_REF").unwrap();
        let gtf = std::env::var("RUSTAR_OTF_GTF").unwrap();
        let params = Parameters::parse_from([
            "rustar-aligner",
            "--genomeDir",
            &base,
            "--readFilesIn",
            "x.fq",
        ]);
        let mut index = GenomeIndex::load(std::path::Path::new(&base), &params).unwrap();
        let ref_params = Parameters::parse_from([
            "rustar-aligner",
            "--genomeDir",
            &reference,
            "--readFilesIn",
            "x.fq",
        ]);
        let star = GenomeIndex::load(std::path::Path::new(&reference), &ref_params).unwrap();

        let exons =
            crate::junction::gtf::parse_gtf_configured(std::path::Path::new(&gtf), "exon", "")
                .unwrap();
        let raw = crate::junction::gtf::extract_junctions_configured(
            exons,
            &index.genome,
            "transcript_id",
        )
        .unwrap();
        let n_real = index.genome.n_genome_real;
        let prepared = raw
            .iter()
            .map(|&(c, s, e, st)| PreparedJunction {
                priority: 20,
                ..sjdb_insert::prepare_junction(c, s, e, st, &index.genome, n_real)
            })
            .collect();
        insert_junctions(&mut index, prepared, star.sjdb_overhang).unwrap();

        assert_eq!(index.prepared_junctions, star.prepared_junctions);
        assert_eq!(index.genome.n_genome, star.genome.n_genome);
        let (a, b) = (index.genome.sequence.view(), star.genome.sequence.view());
        for i in 0..2 * index.genome.n_genome as usize {
            assert_eq!(a.base(i), b.base(i), "genome byte {i}");
        }
        assert_eq!(index.suffix_array.len(), star.suffix_array.len());
        let mut n_diff = 0usize;
        for i in 0..star.suffix_array.len() {
            if index.suffix_array.get(i) != star.suffix_array.get(i) {
                if n_diff < 10 {
                    eprintln!(
                        "SA[{i}]: ours {:#x} star {:#x}",
                        index.suffix_array.get(i),
                        star.suffix_array.get(i)
                    );
                }
                n_diff += 1;
            }
        }
        assert_eq!(n_diff, 0, "{n_diff} SA entries differ");

        // The SA index is a function of genome + SA; rebuild it at the
        // reference's k-mer length (the plain index may use another one).
        // With the same k-mer length the SA index was updated in place, as
        // STAR does; otherwise rebuild it at the reference length.
        let sai = if index.sa_index.nbases == star.sa_index.nbases {
            index.sa_index.clone()
        } else {
            SaIndex::build(&index.genome, &index.suffix_array, star.sa_index.nbases).unwrap()
        };
        assert_eq!(sai.data.len(), star.sa_index.data.len());
        for i in 0..sai.data.len() {
            assert_eq!(sai.data.read(i), star.sa_index.data.read(i), "SAi[{i}]");
        }
    }
}
