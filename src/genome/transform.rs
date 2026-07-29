//! `--genomeTransformType Haploid` (STAR `Genome_transformGenome.cpp`): substitute VCF alleles into
//! the genome before the suffix array is built, so reads carrying the alternate allele align with
//! fewer mismatches.
//!
//! Ported from STAR-rs `crates/star-index/src/transform.rs`. Both **Haploid** (one allele per site)
//! and **Diploid** (genotype-aware, duplicating the genome into `_h1`/`_h2` haplotypes) are
//! implemented, including indels that shift coordinates and split the `transformGenomeBlocks.tsv`
//! block map. With the only supported `--genomeTransformOutput` (`None`, the default) this is a pure
//! `genomeGenerate` transform: the aligner reports transformed-genome coordinates directly, so no
//! align-time back-transform is implemented (that's a separate follow-up, gated in
//! `Parameters::validate`) — for Diploid that also means no `ha:i:` haplotype tag yet.
//!
//! Unlike STAR-rs, which transforms an already-laid-out, bin-padded genome buffer, this operates on
//! rustar-aligner's `Vec<Chromosome>` (name + unpadded base-code sequence) directly, before
//! [`Genome::from_fasta`](super::Genome::from_fasta)'s own padding pass runs. Blocks are computed in
//! two stages: substitution happens in chromosome-local coordinates (no padding involved), then
//! globalized using [`compute_chr_starts`](super::compute_chr_starts) run once over the original
//! chromosome lengths and once over the transformed lengths — the same deterministic function
//! `Genome::from_fasta` itself uses, so the geometry is guaranteed to match the final on-disk index.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::error::Error;
use crate::io::fastq::encode_base;

use super::compute_chr_starts;
use super::fasta::Chromosome;

/// One VCF variant applied to a chromosome: 0-based chr-local `pos`, the reference length, and the
/// alternate allele bytes (raw ASCII letters; [`encode_base`] is applied at substitution time).
#[derive(Debug, Clone)]
pub struct Variant {
    pub pos: u64,
    pub ref_len: usize,
    pub alt: Vec<u8>,
}

impl Variant {
    /// STAR's `len = alt.size() - ref.size()` (0 for a SNV, `>0` insertion, `<0` deletion).
    fn len_delta(&self) -> i64 {
        self.alt.len() as i64 - self.ref_len as i64
    }
}

/// Parse a VCF for the Haploid transform (`Genome_transformGenome.cpp:41-107`): every record on a
/// known chromosome contributes its **first** alternate allele (no genotype filtering); `#` lines and
/// records on unknown chromosomes are skipped. Returns per-chromosome-index variants, each list sorted
/// by position with STAR's overlap filter applied (keep a variant only if it starts at or after the
/// end of the last kept one).
pub fn parse_vcf_haploid(text: &str, chr_name: &[String]) -> BTreeMap<usize, Vec<Variant>> {
    let index: BTreeMap<&str, usize> = chr_name
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let mut per_chr: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 5 {
            continue;
        }
        let Some(&ci) = index.get(f[0]) else {
            continue;
        };
        let Ok(pos1): Result<u64, _> = f[1].parse() else {
            continue;
        };
        if pos1 == 0 {
            continue;
        }
        let ref_allele = f[3].as_bytes();
        // STAR takes the FIRST alternate allele for Haploid (`altV[0]`).
        let alt = f[4].split(',').next().unwrap_or(f[4]).as_bytes().to_vec();
        per_chr.entry(ci).or_default().push(Variant {
            pos: pos1 - 1,
            ref_len: ref_allele.len(),
            alt,
        });
    }

    for variants in per_chr.values_mut() {
        *variants = filter_sort_variants(std::mem::take(variants));
    }
    per_chr
}

/// Parse a VCF for the Diploid transform (`Genome_transformGenome.cpp`'s diploid path): genotype-aware,
/// producing one variant list per haplotype (`[hap0, hap1]`). For each record and each haplotype `ih`,
/// the genotype digit is read from the sample column (field index 9) at character offset `ih*2` (e.g.
/// `0|1` → hap0 digit at offset 0, hap1 digit at offset 2); a genotype of `0` means no variant on that
/// haplotype, otherwise it 1-indexes into the comma-separated ALT list (`altV[gt-1]`). Records without a
/// sample column, or with a non-digit at the expected offset, contribute no variant for that haplotype.
pub fn parse_vcf_diploid(text: &str, chr_name: &[String]) -> [BTreeMap<usize, Vec<Variant>>; 2] {
    let index: BTreeMap<&str, usize> = chr_name
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let mut per_chr: [BTreeMap<usize, Vec<Variant>>; 2] = [BTreeMap::new(), BTreeMap::new()];
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 10 {
            continue;
        }
        let Some(&ci) = index.get(f[0]) else {
            continue;
        };
        let Ok(pos1): Result<u64, _> = f[1].parse() else {
            continue;
        };
        if pos1 == 0 {
            continue;
        }
        let ref_allele = f[3].as_bytes();
        let alt_alleles: Vec<&[u8]> = f[4].split(',').map(str::as_bytes).collect();
        let gt_field = f[9].split(':').next().unwrap_or(f[9]).as_bytes();

        for (ih, per_hap) in per_chr.iter_mut().enumerate() {
            let offset = ih * 2;
            let Some(&c) = gt_field.get(offset) else {
                continue;
            };
            if !c.is_ascii_digit() {
                continue;
            }
            let gt = (c - b'0') as usize;
            if gt == 0 {
                continue;
            }
            let Some(&alt) = alt_alleles.get(gt - 1) else {
                continue;
            };
            per_hap.entry(ci).or_default().push(Variant {
                pos: pos1 - 1,
                ref_len: ref_allele.len(),
                alt: alt.to_vec(),
            });
        }
    }

    for per_hap in &mut per_chr {
        for variants in per_hap.values_mut() {
            *variants = filter_sort_variants(std::mem::take(variants));
        }
    }
    per_chr
}

/// Sort variants by position and apply STAR's overlap filter: keep a variant only if it starts at or
/// after the end of the previously kept variant's reference span.
fn filter_sort_variants(mut variants: Vec<Variant>) -> Vec<Variant> {
    variants.sort_by_key(|v| v.pos);
    let mut kept: Vec<Variant> = Vec::with_capacity(variants.len());
    let mut g0: u64 = 0;
    for v in variants {
        if v.pos >= g0 {
            g0 = v.pos + v.ref_len as u64;
            kept.push(v);
        }
    }
    kept
}

/// One chromosome's Haploid substitution, in chromosome-LOCAL coordinates (no padding): the
/// transformed sequence, and its blocks `[orig_local_start, len, new_local_start]` (globalized by the
/// caller). A chromosome with no variants yields one identity block spanning the whole sequence.
fn transform_one_chromosome(seq: &[u8], variants: &[Variant]) -> (Vec<u8>, Vec<[u64; 3]>) {
    if variants.is_empty() {
        return (seq.to_vec(), vec![[0, seq.len() as u64, 0]]);
    }

    let cl0 = seq.len() as u64;
    let mut gnew: Vec<u8> = Vec::with_capacity(seq.len());
    let mut blocks: Vec<[u64; 3]> = Vec::new();
    let mut iv = 0usize;
    let mut g0: u64 = 0;
    let mut g1: u64 = 0;
    blocks.push([g0, 0, g1]); // first block

    while g0 < cl0 {
        if g0 == variants[iv].pos {
            let v = &variants[iv];
            for &b in &v.alt {
                gnew.push(encode_base(b));
            }
            g0 += v.ref_len as u64;
            g1 += v.alt.len() as u64;
            if v.len_delta() != 0 {
                // Close the previous block; STAR's length formula, then open a new one.
                let last = blocks.last_mut().unwrap();
                last[1] = g0 - v.ref_len as u64 + v.ref_len.min(v.alt.len()) as u64 - last[0];
                blocks.push([g0, 0, g1]);
            }
            if iv < variants.len() - 1 {
                iv += 1;
            }
        } else {
            gnew.push(seq[g0 as usize]);
            g0 += 1;
            g1 += 1;
        }
    }
    if blocks.last().unwrap()[1] == 0 {
        let last = blocks.last_mut().unwrap();
        last[1] = g0 - last[0];
    }
    (gnew, blocks)
}

/// The result of transforming a genome's chromosome list: the substituted chromosomes, and the global
/// `[orig_start, length, new_start]` block map (STAR's `array<uint64,3>`), in chromosome order.
pub struct TransformedGenome {
    pub chromosomes: Vec<Chromosome>,
    pub blocks: Vec<[u64; 3]>,
}

/// Apply the Haploid VCF transform to a genome's chromosome list. `variants` maps a chromosome index
/// (matching `chromosomes`' order) to its filtered, sorted variants (from [`parse_vcf_haploid`]);
/// `chr_bin_nbits` is `--genomeChrBinNbits`, used only to globalize the block coordinates the same way
/// [`Genome::from_fasta`](super::Genome::from_fasta) will pad the transformed chromosomes.
pub fn transform_chromosomes(
    chromosomes: &[Chromosome],
    variants: &BTreeMap<usize, Vec<Variant>>,
    chr_bin_nbits: u32,
) -> TransformedGenome {
    let orig_lengths: Vec<u64> = chromosomes
        .iter()
        .map(|c| c.sequence.len() as u64)
        .collect();

    let empty: Vec<Variant> = Vec::new();
    let mut new_chromosomes = Vec::with_capacity(chromosomes.len());
    let mut local_blocks_per_chr: Vec<Vec<[u64; 3]>> = Vec::with_capacity(chromosomes.len());
    for (ci, chrom) in chromosomes.iter().enumerate() {
        let vs = variants.get(&ci).unwrap_or(&empty);
        let (seq, blocks) = transform_one_chromosome(&chrom.sequence, vs);
        new_chromosomes.push(Chromosome {
            name: chrom.name.clone(),
            sequence: seq,
        });
        local_blocks_per_chr.push(blocks);
    }

    let new_lengths: Vec<u64> = new_chromosomes
        .iter()
        .map(|c| c.sequence.len() as u64)
        .collect();
    let orig_chr_start = compute_chr_starts(&orig_lengths, chr_bin_nbits);
    let new_chr_start = compute_chr_starts(&new_lengths, chr_bin_nbits);

    let mut blocks = Vec::new();
    for (ci, local) in local_blocks_per_chr.into_iter().enumerate() {
        for [lo, len, ln] in local {
            blocks.push([orig_chr_start[ci] + lo, len, new_chr_start[ci] + ln]);
        }
    }

    TransformedGenome {
        chromosomes: new_chromosomes,
        blocks,
    }
}

/// Apply the Diploid VCF transform: run the Haploid single-genome substitution once per haplotype
/// (`variants[0]`/`variants[1]`), rename chromosomes `<name>_h1`/`<name>_h2`, and concatenate hap0's
/// transformed chromosomes then hap1's into one combined chromosome list. Block `orig_start` stays in
/// the single shared original-genome coordinate space for both halves (both haplotypes substitute the
/// same reference); block `new_start` is globalized per-half by [`transform_chromosomes`] and then
/// hap1's is additionally shifted by `offset` — hap0's own padded genome length, computed the same
/// deterministic way [`Genome::from_fasta`](super::Genome::from_fasta) will lay out the final combined
/// genome, so hap1's blocks land exactly where hap1's chromosomes end up after concatenation.
pub fn transform_genome_diploid(
    chromosomes: &[Chromosome],
    variants: &[BTreeMap<usize, Vec<Variant>>; 2],
    chr_bin_nbits: u32,
) -> TransformedGenome {
    let hap0 = transform_chromosomes(chromosomes, &variants[0], chr_bin_nbits);
    let hap1 = transform_chromosomes(chromosomes, &variants[1], chr_bin_nbits);

    let hap0_lengths: Vec<u64> = hap0
        .chromosomes
        .iter()
        .map(|c| c.sequence.len() as u64)
        .collect();
    let offset = *compute_chr_starts(&hap0_lengths, chr_bin_nbits)
        .last()
        .unwrap();

    let mut new_chromosomes = Vec::with_capacity(hap0.chromosomes.len() + hap1.chromosomes.len());
    for c in hap0.chromosomes {
        new_chromosomes.push(Chromosome {
            name: format!("{}_h1", c.name),
            sequence: c.sequence,
        });
    }
    for c in hap1.chromosomes {
        new_chromosomes.push(Chromosome {
            name: format!("{}_h2", c.name),
            sequence: c.sequence,
        });
    }

    let mut blocks = hap0.blocks;
    blocks.extend(hap1.blocks.into_iter().map(|[o, l, n]| [o, l, n + offset]));

    TransformedGenome {
        chromosomes: new_chromosomes,
        blocks,
    }
}

/// Render `transformGenomeBlocks.tsv` (STAR's `transformBlocksWrite`): header `<nBlocks>\t-1`, then
/// one `new_start\tlength\torig_start` line per block (reverting the stored
/// `[orig_start, length, new_start]` order, for reverse conversion).
pub fn blocks_to_tsv(blocks: &[[u64; 3]]) -> String {
    let mut s = format!("{}\t-1\n", blocks.len());
    for b in blocks {
        let _ = writeln!(s, "{}\t{}\t{}", b[2], b[1], b[0]);
    }
    s
}

/// The stop entry appended to a parsed block map.
///
/// The back-transform walks forward over the blocks an exon reaches into and
/// stops on the first one that starts past its end. A sentinel that no
/// coordinate can reach makes that a plain comparison instead of a bounds test
/// on every iteration, which is how STAR writes it (`genomeOutLoad`).
pub const BLOCK_SENTINEL: [u64; 3] = [u64::MAX, 0, 0];

/// Parse `transformGenomeBlocks.tsv` back into the conversion map used by the
/// align-time back-transform.
///
/// The file is written for reverse conversion, so its columns are already
/// `[transformed_start, length, original_start]` — the order
/// [`crate::genome::transform_align::transform_transcript`] wants, and the
/// reverse of the `[u64; 3]` order the build side keeps in memory. The result
/// is sorted by transformed start and terminated by [`BLOCK_SENTINEL`].
pub fn blocks_from_tsv(text: &str) -> Result<Vec<[u64; 3]>, Error> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| Error::Index("transformGenomeBlocks.tsv is empty".to_string()))?;
    let declared: usize = header
        .split_whitespace()
        .next()
        .and_then(|n| n.parse().ok())
        .ok_or_else(|| {
            Error::Index(format!(
                "transformGenomeBlocks.tsv: bad header line {header:?}"
            ))
        })?;

    let mut blocks = Vec::with_capacity(declared + 1);
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let mut f = line.split_whitespace();
        let mut next = |what: &str| -> Result<u64, Error> {
            f.next().and_then(|v| v.parse().ok()).ok_or_else(|| {
                Error::Index(format!("transformGenomeBlocks.tsv: bad {what} in {line:?}"))
            })
        };
        let new_start = next("new_start")?;
        let length = next("length")?;
        let orig_start = next("orig_start")?;
        blocks.push([new_start, length, orig_start]);
    }

    if blocks.len() != declared {
        return Err(Error::Index(format!(
            "transformGenomeBlocks.tsv declares {declared} blocks but has {}",
            blocks.len()
        )));
    }
    // The map is written in chromosome order, which is already ascending by
    // transformed start; sorting is a cheap guarantee rather than a fix.
    blocks.sort_unstable();
    blocks.push(BLOCK_SENTINEL);
    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_first_alt_and_filters_overlaps() {
        let chr = vec!["chr1".to_string(), "chr2".to_string()];
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10\t.\tG\tA,T\t.\t.\t.
chrX\t5\t.\tG\tA\t.\t.\t.
chr2\t3\t.\tC\tT\t.\t.\t.
";
        let v = parse_vcf_haploid(vcf, &chr);
        assert_eq!(v[&0].len(), 1);
        assert_eq!(v[&0][0].pos, 9); // 0-based
        assert_eq!(v[&0][0].alt, b"A"); // first alt only
        assert_eq!(v[&1][0].pos, 2);
        assert_eq!(v[&1][0].alt, b"T");
        assert!(!v.contains_key(&2)); // chrX not in index; only 2 contigs
    }

    #[test]
    fn snv_substitution_and_identity_blocks() {
        let chromosomes = vec![
            Chromosome {
                name: "chr1".to_string(),
                sequence: vec![2, 2, 2, 2, 2], // GGGGG
            },
            Chromosome {
                name: "chr2".to_string(),
                sequence: vec![2, 2, 2, 2], // GGGG
            },
        ];
        let mut variants = BTreeMap::new();
        variants.insert(
            0usize,
            vec![Variant {
                pos: 2,
                ref_len: 1,
                alt: b"A".to_vec(),
            }],
        );
        let t = transform_chromosomes(&chromosomes, &variants, 18);
        assert_eq!(t.chromosomes[0].sequence[2], 0); // A code
        assert_eq!(t.chromosomes[0].sequence.len(), 5);
        assert_eq!(t.chromosomes[1].sequence.len(), 4);
        // One identity block per chromosome; lengths unchanged.
        let cs = compute_chr_starts(&[5, 4], 18);
        assert_eq!(t.blocks, vec![[cs[0], 5, cs[0]], [cs[1], 4, cs[1]]]);
    }

    #[test]
    fn indel_shifts_coordinates_and_splits_blocks() {
        // chr1 = 6 bases; a 1-base insertion at pos 2 (ref C -> alt CTT, +2) then a deletion at pos 4
        // (ref GA -> alt G, -1). Net length 6 + 2 - 1 = 7.
        let chromosomes = vec![Chromosome {
            name: "chr1".to_string(),
            sequence: vec![0, 0, 1, 2, 2, 0], // AACGGA
        }];
        let mut variants = BTreeMap::new();
        variants.insert(
            0usize,
            vec![
                Variant {
                    pos: 2,
                    ref_len: 1,
                    alt: b"CTT".to_vec(),
                },
                Variant {
                    pos: 4,
                    ref_len: 2,
                    alt: b"G".to_vec(),
                },
            ],
        );
        let t = transform_chromosomes(&chromosomes, &variants, 18);
        assert_eq!(t.chromosomes[0].sequence.len(), 7);
        let seq: Vec<u8> = t.chromosomes[0]
            .sequence
            .iter()
            .map(|&c| b"ACGT"[c as usize])
            .collect();
        assert_eq!(&seq, b"AACTTGG");
        // Blocks split at each indel: [orig_start, length, new_start] (global, chr1 starts at 0).
        assert_eq!(t.blocks[0], [0, 3, 0]);
    }

    #[test]
    fn blocks_to_tsv_reverts_column_order() {
        let tsv = blocks_to_tsv(&[[0, 5, 0], [5, 2, 7]]);
        assert_eq!(tsv, "2\t-1\n0\t5\t0\n7\t2\t5\n");
    }

    #[test]
    fn diploid_genotype_parsing() {
        let chr = vec!["chr1".to_string()];
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t3\t.\tC\tA,T\t.\t.\t.\tGT\t0|1
chr1\t6\t.\tG\tA\t.\t.\t.\tGT\t1/1
";
        let [h0, h1] = parse_vcf_diploid(vcf, &chr);
        // Line 1 (0-based pos 2, het "0|1"): hap0 gt=0 -> no variant; hap1 gt=1 -> alt[0]="A".
        // Line 2 (0-based pos 5, hom "1/1"): both haplotypes get alt[0]="A".
        assert_eq!(h0[&0].len(), 1);
        assert_eq!(h0[&0][0].pos, 5);
        assert_eq!(h1[&0].len(), 2);
        assert_eq!(h1[&0][0].pos, 2);
        assert_eq!(h1[&0][0].alt, b"A");
        assert_eq!(h1[&0][1].pos, 5);
    }

    #[test]
    fn diploid_transform_concatenates_haplotypes_with_offset() {
        let chromosomes = vec![Chromosome {
            name: "chr1".to_string(),
            sequence: vec![0, 0, 0, 0], // AAAA
        }];
        let mut h0 = BTreeMap::new();
        h0.insert(
            0usize,
            vec![Variant {
                pos: 1,
                ref_len: 1,
                alt: b"C".to_vec(),
            }],
        );
        let h1: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        let variants = [h0, h1];

        let t = transform_genome_diploid(&chromosomes, &variants, 18);
        assert_eq!(t.chromosomes.len(), 2);
        assert_eq!(t.chromosomes[0].name, "chr1_h1");
        assert_eq!(t.chromosomes[1].name, "chr1_h2");
        assert_eq!(t.chromosomes[0].sequence[1], 1); // C substituted on hap0
        assert_eq!(t.chromosomes[1].sequence, vec![0, 0, 0, 0]); // hap1 identity

        let offset = *compute_chr_starts(&[4], 18).last().unwrap();
        assert_eq!(t.blocks.len(), 2); // one identity-length block per haplotype (SNV doesn't split)
        assert_eq!(t.blocks[0][0], 0); // orig_start
        assert_eq!(t.blocks[1][0], 0); // same shared orig_start, not offset
        assert_eq!(t.blocks[1][2], offset); // new_start shifted by hap0's padded genome length
    }

    #[test]
    fn block_map_round_trips_through_the_tsv() {
        // In-memory order is [orig_start, length, new_start]; the file stores
        // the reverse, which is also what the back-transform reads.
        let written = vec![[0u64, 50, 0], [55, 100, 50]];
        let parsed = blocks_from_tsv(&blocks_to_tsv(&written)).unwrap();
        assert_eq!(
            parsed,
            vec![[0, 50, 0], [50, 100, 55], BLOCK_SENTINEL],
            "the parsed map is keyed on transformed coordinates, plus the stop entry"
        );
    }

    #[test]
    fn a_block_count_that_disagrees_with_the_body_is_an_error() {
        let text = "3\t-1\n0\t50\t0\n50\t100\t55\n";
        let err = blocks_from_tsv(text).unwrap_err().to_string();
        assert!(err.contains("declares 3 blocks but has 2"), "{err}");
    }
}
