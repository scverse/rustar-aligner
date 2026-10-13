//! `--runMode liftOver`: arithmetic lift-over of a GTF through a UCSC chain file (STAR `Chain.cpp`).
//!
//! A chain file describes how blocks of an old assembly (`tName`) map to a new one (`qName`). This
//! lifts each GTF feature's start/end coordinate through the blocks and re-emits the line under the
//! new chromosome name; features that cannot be lifted (an endpoint falls in an unaligned gap, or the
//! lifted end precedes the lifted start) go to a `.unlifted` sidecar. No genome index is involved.
//!
//! Faithful to STAR's quirks: only a single chain per chromosome is supported, coordinates are
//! unsigned with `-1` as the "impossible" sentinel, and the remainder of each GTF line (everything
//! after the end coordinate, including its leading tab) is preserved verbatim.
//!
//! Near-verbatim port of the sister project STAR-rs's `crates/star-index/src/chain.rs`.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use crate::params::Parameters;

/// One chain's aligned blocks for a single old chromosome: block starts in the old (`b_start1`) and
/// new (`b_start2`) assemblies and the block lengths (`b_len`), plus the new chromosome name.
#[derive(Debug, Default, Clone)]
pub struct OneChain {
    /// New (`qName`) chromosome name the coordinates lift to.
    pub chr2: String,
    /// Block starts in the old assembly.
    pub b_start1: Vec<u64>,
    /// Block starts in the new assembly.
    pub b_start2: Vec<u64>,
    /// Block lengths (shared by both assemblies).
    pub b_len: Vec<u64>,
}

/// STAR `binarySearch1a`: the index of the last element of the sorted `xs` that is `<= x`, or `-1`
/// when `x < xs[0]`. Ties resolve to the last equal element. `xs` must be non-empty.
fn binary_search_1a(x: u64, xs: &[u64]) -> i64 {
    let n = xs.len();
    if x > xs[n - 1] {
        return n as i64 - 1;
    } else if x < xs[0] {
        return -1;
    }
    let (mut i1, mut i2) = (0i64, n as i64 - 1);
    while i2 > i1 + 1 {
        let i3 = i1 + (i2 - i1) / 2;
        if xs[i3 as usize] > x {
            i2 = i3;
        } else {
            i1 = i3;
        }
    }
    while (i1 as usize) < n - 1 && x == xs[(i1 + 1) as usize] {
        i1 += 1;
    }
    i1
}

/// Parse a chain file (STAR `Chain::chainLoad`). Each header line (`chain score tName tSize tStrand
/// tStart tEnd qName qSize qStrand qStart qEnd id`) opens a chromosome; block lines (`size dt dq`)
/// accumulate the aligned blocks; a lone `size` closes the chain. Fields are whitespace-delimited.
pub fn chain_load(text: &str) -> BTreeMap<String, OneChain> {
    let mut chains: BTreeMap<String, OneChain> = BTreeMap::new();
    let mut chr1 = String::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let f = |i: usize| fields.get(i).copied().unwrap_or("");
        if f(0).is_empty() {
            // empty line
        } else if f(1).is_empty() {
            // end of chain: the last block's length
            if let Some(ch) = chains.get_mut(&chr1) {
                ch.b_len.push(parse_u64(f(0)));
            }
        } else if f(3).is_empty() {
            // normal block "size dt dq": push the block length, then the next block's starts
            if let Some(ch) = chains.get_mut(&chr1) {
                let blen = parse_u64(f(0));
                ch.b_len.push(blen);
                let s1 = ch.b_start1.last().copied().unwrap_or(0) + blen + parse_u64(f(1));
                ch.b_start1.push(s1);
                let s2 = ch.b_start2.last().copied().unwrap_or(0) + blen + parse_u64(f(2));
                ch.b_start2.push(s2);
            }
        } else {
            // chain header: open (or extend) the chromosome's chain
            chr1 = f(2).to_string();
            let ch = chains.entry(chr1.clone()).or_default();
            ch.chr2 = f(7).to_string();
            ch.b_start1.push(parse_u64(f(5)));
            ch.b_start2.push(parse_u64(f(10)));
        }
    }
    chains
}

/// STAR `std::stoi`-style leading-integer parse (ignores any trailing non-digits); `0` on no digits.
fn parse_u64(s: &str) -> u64 {
    let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    s[..end].parse().unwrap_or(0)
}

/// Read the first `n` whitespace-delimited tokens of `line` (STAR's `istringstream >>`), returning
/// them plus the unread remainder (from just after the `n`-th token, including its leading
/// whitespace, i.e. `istringstream::rdbuf`). Returns `None` if the line has fewer than `n` tokens.
fn tokens_and_rest(line: &str, n: usize) -> Option<(Vec<&str>, &str)> {
    let bytes = line.as_bytes();
    let mut tokens = Vec::with_capacity(n);
    let mut i = 0;
    for _ in 0..n {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            return None;
        }
        let start = i;
        while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        tokens.push(&line[start..i]);
    }
    Some((tokens, &line[i..]))
}

/// Lift one coordinate (STAR `liftOverGTF`): `is_end` selects the right-end edge rule (map an
/// out-of-block start to the next block's start, an out-of-block end to the current block's end).
/// Returns `None` for STAR's `-1` (impossible to lift).
fn lift_coord(ch: &OneChain, c1: u64, is_end: bool) -> Option<u64> {
    let i1 = binary_search_1a(c1, &ch.b_start1);
    if i1 >= 0 && c1 < ch.b_start1[i1 as usize] + ch.b_len[i1 as usize] {
        // inside the block: a straight shift into the new assembly
        Some(ch.b_start2[i1 as usize] + c1 - ch.b_start1[i1 as usize])
    } else if !is_end && i1 < ch.b_start1.len() as i64 - 1 {
        // left end outside a block: the start of the next block (i1 == -1 gives block 0)
        Some(ch.b_start2[(i1 + 1) as usize])
    } else if is_end && i1 >= 0 {
        // right end outside a block: the end of the current block
        Some(ch.b_start2[i1 as usize] + ch.b_len[i1 as usize] - 1)
    } else {
        None
    }
}

/// Lift a whole GTF (STAR `Chain::liftOverGTF`). Returns `(lifted, unlifted)`: the lifted GTF (each
/// feature re-emitted under the new chromosome with lifted start/end and its remaining fields kept
/// verbatim) and the sidecar of lines that could not be lifted. `Err` names the offending chromosome
/// when the GTF references one the chain file lacks (STAR exits with that error).
pub fn lift_over_gtf(
    chains: &BTreeMap<String, OneChain>,
    gtf: &str,
    chain_file: &str,
) -> Result<(String, String), String> {
    let mut lifted = String::new();
    let mut unlifted = String::new();
    for line in gtf.lines() {
        let first = line.split_whitespace().next().unwrap_or("");
        if first.is_empty() || first.starts_with('#') {
            continue; // empty or comment line
        }
        let ch = chains.get(first).ok_or_else(|| {
            format!("GTF contains chromosome {first} not present in the chain file {chain_file}")
        })?;
        let Some((tokens, rest)) = tokens_and_rest(line, 5) else {
            // Fewer than 5 fields: not a liftable feature line.
            unlifted.push_str(line);
            unlifted.push('\n');
            continue;
        };
        let (str1, str2) = (tokens[1], tokens[2]);
        let c2_start = lift_coord(ch, parse_u64(tokens[3]), false);
        let c2_end = lift_coord(ch, parse_u64(tokens[4]), true);
        match (c2_start, c2_end) {
            (Some(s), Some(e)) if e >= s => {
                let _ = writeln!(lifted, "{}\t{str1}\t{str2}\t{s}\t{e}{rest}", ch.chr2);
            }
            _ => {
                unlifted.push_str(line);
                unlifted.push('\n');
            }
        }
    }
    Ok((lifted, unlifted))
}

/// `--runMode liftOver`: lift the `--sjdbGTFfile` coordinates through the first `--genomeChainFiles`
/// chain file, writing `<prefix>GTFliftOver_1.gtf` plus its `.unlifted` sidecar. As in STAR (whose
/// `exit(0)` sits inside the chain-file loop), only the first chain file is processed.
pub fn run(params: &Parameters) -> anyhow::Result<()> {
    let chain_file = params
        .genome_chain_files
        .first()
        .ok_or_else(|| anyhow::anyhow!("--runMode liftOver requires --genomeChainFiles <chain>"))?;
    let gtf_file = params
        .sjdb_gtf_file
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--runMode liftOver requires --sjdbGTFfile <gtf>"))?;

    let chains = chain_load(&std::fs::read_to_string(chain_file)?);
    let gtf = crate::junction::gtf::read_annotation_to_string(gtf_file)?;
    let (lifted, unlifted) =
        lift_over_gtf(&chains, &gtf, &chain_file.to_string_lossy()).map_err(anyhow::Error::msg)?;

    let out = params.output_path("GTFliftOver_1.gtf");
    if let Some(parent) = out.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(&out, lifted)?;
    std::fs::write(path_with_suffix(&out, ".unlifted"), unlifted)?;
    Ok(())
}

fn path_with_suffix(path: &Path, suffix: &str) -> std::path::PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    std::path::PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> BTreeMap<String, OneChain> {
        // Old chr1 -> new chrA: block0 [0,100), gap [100,110), block1 [110,210).
        chain_load("chain 1000 chr1 1000 + 0 1000 chrA 1000 + 0 1000 1\n100 10 10\n100\n\n")
    }

    #[test]
    fn binary_search_1a_boundaries() {
        let xs = [0u64, 110];
        assert_eq!(binary_search_1a(0, &xs), 0);
        assert_eq!(binary_search_1a(50, &xs), 0);
        assert_eq!(binary_search_1a(110, &xs), 1);
        assert_eq!(binary_search_1a(200, &xs), 1);
    }

    #[test]
    fn chain_load_builds_the_blocks() {
        let ch = &sample()["chr1"];
        assert_eq!(ch.chr2, "chrA");
        assert_eq!(ch.b_start1, [0, 110]);
        assert_eq!(ch.b_start2, [0, 110]);
        assert_eq!(ch.b_len, [100, 100]);
    }

    #[test]
    fn lift_inside_straddle_and_gap() {
        let chains = sample();
        let gtf = "#comment\n\
            chr1\tsrc\texon\t10\t50\t.\t+\t.\tgene_id \"g1\";\n\
            chr1\tsrc\texon\t90\t150\t.\t+\t.\tgene_id \"g2\";\n\
            chr1\tsrc\texon\t102\t105\t.\t+\t.\tgene_id \"g3\";\n";
        let (lifted, unlifted) = lift_over_gtf(&chains, gtf, "x.chain").unwrap();
        assert_eq!(
            lifted,
            "chrA\tsrc\texon\t10\t50\t.\t+\t.\tgene_id \"g1\";\n\
             chrA\tsrc\texon\t90\t150\t.\t+\t.\tgene_id \"g2\";\n"
        );
        assert_eq!(
            unlifted,
            "chr1\tsrc\texon\t102\t105\t.\t+\t.\tgene_id \"g3\";\n"
        );
    }

    #[test]
    fn missing_chromosome_errors() {
        let chains = sample();
        let err = lift_over_gtf(&chains, "chrZ\ts\te\t1\t2\t.\t+\t.\tx\n", "x.chain").unwrap_err();
        assert!(err.contains("chrZ"));
    }
}
