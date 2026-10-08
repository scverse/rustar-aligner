use std::fs::File;
use std::path::Path;

use byteorder::{LittleEndian, ReadBytesExt};

use crate::error::Error;
use crate::genome::Genome;
use crate::index::GenomeIndex;
use crate::index::packed_array::PackedArray;
use crate::index::sa_index::SaIndex;
use crate::index::suffix_array::SuffixArray;
use crate::junction::SpliceJunctionDb;
use crate::junction::sjdb_insert;
use crate::params::Parameters;
use crate::quant::transcriptome::TranscriptomeIndex;

impl GenomeIndex {
    /// Load a genome index from disk.
    ///
    /// Reads Genome, SA, and SAindex files from the specified directory.
    pub fn load(genome_dir: &Path, params: &Parameters) -> Result<Self, Error> {
        log::info!("Loading genome from {}...", genome_dir.display());

        // Refuse an index older than --versionGenome before reading a byte of
        // it: an older layout read with today's reader is silently wrong, not
        // loudly wrong.
        check_genome_version(genome_dir, &params.version_genome)?;

        // Load Genome file
        let genome = load_genome(genome_dir, params)?;
        log::info!(
            "Loaded genome: {} chromosomes, {} bytes",
            genome.n_chr_real,
            genome.n_genome
        );

        // Load SA file
        let suffix_array = load_suffix_array(genome_dir, &genome)?;
        log::info!("Loaded suffix array: {} entries", suffix_array.len());

        // Load SAindex file
        let sa_index = load_sa_index(genome_dir, suffix_array.gstrand_bit)?;
        log::info!(
            "Loaded SA index: nbases={}, {} indices",
            sa_index.nbases,
            sa_index.data.len()
        );

        // Load prepared junctions from the index (sjdbInfo.txt) if present.
        // STAR appends a Gsj flanking-sequence buffer to the genome at build
        // time; align-time code needs the parsed junctions to (a) decode SA hits
        // that land inside that buffer back to real (donor, acceptor) positions,
        // and (b) recognise annotated junctions when aligning against a
        // pre-built annotated index.
        let sjdb_info_path = genome_dir.join("sjdbInfo.txt");
        let (prepared_junctions, sjdb_overhang) = if sjdb_info_path.exists() {
            let tab = sjdb_insert::read_sjdb_info_tab(&sjdb_info_path, &genome)?;
            log::info!(
                "Loaded sjdbInfo.txt: {} junctions, sjdbOverhang={}",
                tab.junctions.len(),
                tab.sjdb_overhang,
            );
            (tab.junctions, tab.sjdb_overhang)
        } else {
            (Vec::new(), 0)
        };

        // Build the annotated-junction database consulted at stitch time.
        //   - If a GTF is supplied at align time, parse it (STAR's on-the-fly path).
        //   - Otherwise fall back to the junctions stored in the index
        //     (sjdbInfo.txt). Without this fallback, the standard workflow —
        //     build the index once with `--sjdbGTFfile`, then align with only
        //     `--genomeDir` — would treat every junction as novel (the runtime
        //     db would be empty), losing all `sjdbScore` bonuses and annotated
        //     junction recognition. Keyed on the stored (post-sjdbPrepare) donor/
        //     acceptor coordinates, matching what the stitch scan produces.
        let junction_db = if let Some(ref gtf_path) = params.sjdb_gtf_file {
            SpliceJunctionDb::from_gtf_configured(
                gtf_path,
                &genome,
                &params.sjdb_gtf_feature_exon,
                &params.sjdb_gtf_chr_prefix,
                &params.sjdb_gtf_tag_exon_parent_transcript,
            )?
        } else if !prepared_junctions.is_empty() {
            let raw: Vec<(usize, u64, u64, u8)> = prepared_junctions
                .iter()
                .map(|j| (j.chr_idx, j.stored_start(), j.stored_end(), j.strand))
                .collect();
            log::info!(
                "No GTF at align time; loaded {} annotated junctions from index sjdbInfo.txt",
                raw.len()
            );
            SpliceJunctionDb::from_raw_junctions(&raw)
        } else {
            log::info!("No GTF file provided, all junctions will be novel");
            SpliceJunctionDb::empty()
        };

        log::info!(
            "Junction database loaded: {} annotated junctions",
            junction_db.len()
        );

        // Prefer STAR-compatible transcriptInfo.tab / exonInfo.tab /
        // geneInfo.tab over re-parsing the GTF at align time. If the files
        // aren't present (legacy rustar-aligner index), fall back to on-the-fly
        // construction from the GTF when one is supplied — this matches
        // STAR's behavior in `sjdbInsertJunctions.cpp` (re-parse and regenerate).
        let transcriptome = if genome_dir.join("transcriptInfo.tab").exists() {
            log::info!(
                "Loading transcriptome index files from {}",
                genome_dir.display()
            );
            Some(TranscriptomeIndex::from_index_dir(genome_dir, &genome)?)
        } else if let Some(ref gtf_path) = params.sjdb_gtf_file {
            log::warn!(
                "transcriptInfo.tab not found in {}; re-parsing GTF at align time",
                genome_dir.display()
            );
            let exons = crate::junction::gtf::parse_gtf_configured(
                gtf_path,
                &params.sjdb_gtf_feature_exon,
                &params.sjdb_gtf_chr_prefix,
            )?;
            Some(TranscriptomeIndex::from_gtf_exons_configured(
                &exons,
                &genome,
                &params.sjdb_gtf_tag_exon_parent_transcript,
                &params.sjdb_gtf_tag_exon_parent_gene,
                &params.sjdb_gtf_tag_exon_parent_gene_name,
                &params.sjdb_gtf_tag_exon_parent_gene_type,
            )?)
        } else {
            None
        };

        if let Some(ref tr) = transcriptome {
            log::info!(
                "Transcriptome index ready: {} transcripts, {} genes",
                tr.n_transcripts(),
                tr.gene_ids.len()
            );
        }

        Ok(GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db,
            transcriptome,
            prepared_junctions,
            sjdb_overhang,
        })
    }
}

/// Parse a STAR genome version string (`2.7.4a`) into a comparable tuple:
/// the numeric components, then the trailing letter suffix (`a` -> 1).
///
/// STAR compares these versions as strings, which orders `2.7.10a` before
/// `2.7.4a`; comparing components avoids that trap.
fn parse_genome_version(v: &str) -> Option<(Vec<u32>, u32)> {
    let v = v.trim();
    if v.is_empty() {
        return None;
    }
    let digits_end = v.rfind(|c: char| c.is_ascii_digit()).map_or(0, |i| i + 1);
    let (numeric, suffix) = v.split_at(digits_end);
    let mut parts = Vec::new();
    for p in numeric.split('.') {
        parts.push(p.parse::<u32>().ok()?);
    }
    let suffix_rank = suffix.bytes().next().map_or(0, |b| {
        u32::from(b.to_ascii_lowercase().saturating_sub(b'a')) + 1
    });
    Some((parts, suffix_rank))
}

/// Read `versionGenome` from `genomeParameters.txt`.
fn read_genome_version(genome_dir: &Path) -> Result<Option<String>, Error> {
    let path = genome_dir.join("genomeParameters.txt");
    let contents = match std::fs::read_to_string(&path) {
        Ok(s) => s,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(Error::io(e, &path)),
    };
    for line in contents.lines() {
        if let Some(rest) = line.strip_prefix("versionGenome")
            && let Some(value) = rest.split_whitespace().next()
        {
            return Ok(Some(value.to_string()));
        }
    }
    Ok(None)
}

/// Fail when the index on disk is older than `required` (`--versionGenome`).
///
/// An index with no recorded version, or a version neither side can parse, is
/// accepted with a warning: refusing it would break directories that load
/// correctly today, and the version line is advisory metadata, not a checksum.
pub(crate) fn check_genome_version(genome_dir: &Path, required: &str) -> Result<(), Error> {
    let Some(found) = read_genome_version(genome_dir)? else {
        log::warn!(
            "{} has no versionGenome line; skipping the version check",
            genome_dir.join("genomeParameters.txt").display()
        );
        return Ok(());
    };
    let (Some(found_v), Some(required_v)) =
        (parse_genome_version(&found), parse_genome_version(required))
    else {
        log::warn!("could not compare genome version '{found}' against '{required}'");
        return Ok(());
    };
    if found_v < required_v {
        return Err(Error::Index(format!(
            "genome index in {} has versionGenome {found}, older than the required {required}.              Regenerate the index with --runMode genomeGenerate, or lower --versionGenome if              you know the layout is compatible",
            genome_dir.display()
        )));
    }
    Ok(())
}

/// Read `genomeFileSizes\t<n_genome> <sa_size>` from genomeParameters.txt
/// and return the first field (total genome byte count, including Gsj if
/// sjdb was baked in). Returns `Ok(None)` if the file or line is absent,
/// leaving the caller to fall back to the chr_start boundary.
fn read_genome_file_size(genome_dir: &Path) -> Result<Option<u64>, Error> {
    let path = genome_dir.join("genomeParameters.txt");
    let contents = match std::fs::read_to_string(&path) {
        Ok(s) => s,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(Error::io(e, &path)),
    };
    for line in contents.lines() {
        if let Some(rest) = line.strip_prefix("genomeFileSizes\t")
            && let Some(first) = rest.split_whitespace().next()
            && let Ok(v) = first.parse::<u64>()
        {
            return Ok(Some(v));
        }
    }
    Ok(None)
}

/// Load genome from disk.
fn load_genome(genome_dir: &Path, _params: &Parameters) -> Result<Genome, Error> {
    // Read chromosome metadata
    let chr_name_path = genome_dir.join("chrName.txt");
    let chr_name_contents =
        std::fs::read_to_string(&chr_name_path).map_err(|e| Error::io(e, &chr_name_path))?;
    let chr_name: Vec<String> = chr_name_contents.lines().map(ToString::to_string).collect();

    let chr_length_path = genome_dir.join("chrLength.txt");
    let chr_length_contents =
        std::fs::read_to_string(&chr_length_path).map_err(|e| Error::io(e, &chr_length_path))?;
    let chr_length: Vec<u64> = chr_length_contents
        .lines()
        .map(|s| s.parse().unwrap())
        .collect();

    let chr_start_path = genome_dir.join("chrStart.txt");
    let chr_start_contents =
        std::fs::read_to_string(&chr_start_path).map_err(|e| Error::io(e, &chr_start_path))?;
    let chr_start: Vec<u64> = chr_start_contents
        .lines()
        .map(|s| s.parse().unwrap())
        .collect();

    let n_chr_real = chr_name.len();

    // `chr_start[n_chr_real]` is the forward boundary of REAL chromosomes
    // only — it stays pinned at the pre-sjdb value in STAR (`chrStart.txt`).
    // When sjdb has been baked into the index, the total genome size
    // (real + Gsj) lives in `genomeParameters.txt` under `genomeFileSizes`.
    // Prefer that value; fall back to the chr_start boundary for indices
    // built without a GTF.
    let n_genome_real = chr_start[n_chr_real];
    let n_genome = read_genome_file_size(genome_dir)?.unwrap_or(n_genome_real);

    // Memory-map the Genome sequence file (forward strand only, `n_genome`
    // bytes). The reverse-complement half is computed on access by
    // `GenomeSeq::base`, so the ~`n_genome`-byte RC buffer is never
    // materialized and the forward bytes are reclaimable file-backed pages
    // rather than an anonymous `Vec`. The genome is accessed by single-byte
    // lookups during alignment, which `base` serves from the map.
    let genome_path = genome_dir.join("Genome");
    let file = File::open(&genome_path).map_err(|e| Error::io(e, &genome_path))?;
    // SAFETY: Genome is opened read-only and never mutated while loaded.
    let mmap = unsafe { memmap2::Mmap::map(&file).map_err(|e| Error::io(e, &genome_path))? };
    // Each `compare_seq_to_genome` touches only a read-length run of bytes (≪ one
    // page) at a genome position that is effectively random across reads, so kernel
    // readahead past that page is wasted I/O — same rationale as the SA/SAindex maps.
    advise_random(&mmap);

    if mmap.len() != n_genome as usize {
        return Err(Error::Index(format!(
            "Genome file size mismatch: expected {} bytes, got {}",
            n_genome,
            mmap.len()
        )));
    }

    let sequence = crate::genome::GenomeSeq::Mapped {
        fwd: std::sync::Arc::new(mmap),
        n_genome: n_genome as usize,
    };

    Ok(Genome {
        sequence,
        n_genome,
        n_genome_real,
        n_chr_real,
        chr_name,
        chr_length,
        chr_start,
        // Loading transformGenomeBlocks.tsv back is only needed for the
        // (not yet implemented) align-time back-transform.
        transform_blocks: None,
    })
}

/// Load suffix array from disk.
///
/// The `SA` file is **memory-mapped** rather than read into a `Vec`: it is the
/// largest index component (≈21 GB for mouse) and is accessed by random binary
/// search during alignment. mmap keeps it as reclaimable file-backed memory
/// (demand-loaded, dropped — not swapped — under pressure) instead of an
/// un-reclaimable anonymous allocation. `MADV_RANDOM` disables readahead, which
/// would waste I/O on the random access pattern.
/// Best-effort `MADV_RANDOM` on a read-only mmap. `madvise` (and `memmap2::Advice`)
/// is Unix-only, so this is a no-op on platforms without it (e.g. Windows).
#[cfg(unix)]
fn advise_random(mmap: &memmap2::Mmap) {
    let _ = mmap.advise(memmap2::Advice::Random); // best-effort; ignore if unsupported
}
#[cfg(not(unix))]
fn advise_random(_mmap: &memmap2::Mmap) {}

fn load_suffix_array(genome_dir: &Path, genome: &Genome) -> Result<SuffixArray, Error> {
    let sa_path = genome_dir.join("SA");
    let file = File::open(&sa_path).map_err(|e| Error::io(e, &sa_path))?;
    // SAFETY: the SA file is opened read-only and not mutated elsewhere while
    // the index is loaded; the mapping is only ever read.
    let mmap = unsafe { memmap2::Mmap::map(&file).map_err(|e| Error::io(e, &sa_path))? };
    advise_random(&mmap);

    let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
    let word_length = gstrand_bit + 1;
    let gstrand_mask = (1u64 << gstrand_bit) - 1;

    // Calculate expected length from file size
    // Formula from STAR: lengthByte = (length-1)*wordLength/8 + 8
    // We need to solve for length, accounting for integer division:
    // total_bits = (lengthByte - 8) * 8
    // length = (total_bits / wordLength) + 1
    // BUT we need ceiling division to account for partial entries
    let length_byte = mmap.len();
    let length = if length_byte < 8 {
        0
    } else {
        let total_bits = (length_byte - 8) * 8;
        let entries = total_bits.div_ceil(word_length as usize);
        entries + 1
    };

    let data = PackedArray::from_mmap(word_length, length, mmap);

    Ok(SuffixArray {
        data,
        gstrand_bit,
        gstrand_mask,
    })
}

/// Load SA index from disk.
///
/// The small fixed header (`nbases` + the `genomeSAindexStart` array) is read
/// normally; the packed-data region (≈1.8 GB for mouse) is **memory-mapped**
/// from its byte offset for the same reason as the SA — reclaimable, demand-
/// loaded file-backed memory instead of an anonymous `Vec`.
fn load_sa_index(genome_dir: &Path, gstrand_bit: u32) -> Result<SaIndex, Error> {
    let sai_path = genome_dir.join("SAindex");
    let mut file = File::open(&sai_path).map_err(|e| Error::io(e, &sai_path))?;

    // Read nbases (u64)
    let nbases = file
        .read_u64::<LittleEndian>()
        .map_err(|e| Error::io(e, &sai_path))? as u32;

    // Read genomeSAindexStart array (nbases + 1 entries)
    let mut genome_sa_index_start = Vec::with_capacity((nbases + 1) as usize);
    for _ in 0..=nbases {
        let val = file
            .read_u64::<LittleEndian>()
            .map_err(|e| Error::io(e, &sai_path))?;
        genome_sa_index_start.push(val);
    }

    // Map the packed-data region: header is `nbases` (8B) + (nbases+1)×8B.
    let header_len = 8 + 8 * (u64::from(nbases) + 1);
    // SAFETY: SAindex is opened read-only and never mutated while loaded.
    // memmap2 handles non-page-aligned offsets internally; the map runs from
    // `header_len` to EOF and is only ever read.
    let mmap = unsafe {
        memmap2::MmapOptions::new()
            .offset(header_len)
            .map(&file)
            .map_err(|e| Error::io(e, &sai_path))?
    };
    advise_random(&mmap);

    let word_length = gstrand_bit + 3;
    let num_indices = SaIndex::calculate_num_indices(nbases);

    let data = PackedArray::from_mmap(word_length, num_indices as usize, mmap);

    Ok(SaIndex {
        nbases,
        genome_sa_index_start,
        data,
        word_length,
        gstrand_bit,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn load_generated_index() {
        // Create a simple genome
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();

        let dir = tempfile::tempdir().unwrap();

        let args = vec![
            "rustar-aligner",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeDir",
            dir.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            "2",
            "--genomeSAindexNbases",
            "1",
        ];

        let params = Parameters::parse_from(args.clone());

        // Build index
        let index = GenomeIndex::build(&params).unwrap();
        index.write(dir.path(), &params).unwrap();

        // Load index back
        let loaded_index = GenomeIndex::load(dir.path(), &params).unwrap();

        // Verify
        assert_eq!(loaded_index.genome.n_genome, index.genome.n_genome);
        assert_eq!(loaded_index.genome.n_chr_real, index.genome.n_chr_real);
        assert_eq!(loaded_index.suffix_array.len(), index.suffix_array.len());
        assert_eq!(loaded_index.sa_index.nbases, index.sa_index.nbases);
        assert_eq!(loaded_index.sa_index.data.len(), index.sa_index.data.len());

        // Verify first few SA entries match
        for i in 0..loaded_index.suffix_array.len().min(5) {
            assert_eq!(loaded_index.suffix_array.get(i), index.suffix_array.get(i));
        }
    }

    // ── --versionGenome ──────────────────────────────────────────────────

    fn write_genome_params(dir: &std::path::Path, version_line: &str) {
        std::fs::write(dir.join("genomeParameters.txt"), version_line).unwrap();
    }

    #[test]
    fn version_genome_rejects_an_older_index() {
        let dir = tempfile::tempdir().unwrap();
        write_genome_params(dir.path(), "versionGenome\t2.7.1a\n");
        let err = check_genome_version(dir.path(), "2.7.4a").unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("2.7.1a"),
            "error should quote the found version: {msg}"
        );
        assert!(
            msg.contains("2.7.4a"),
            "error should quote the required version: {msg}"
        );
    }

    #[test]
    fn version_genome_accepts_equal_and_newer_indices() {
        let dir = tempfile::tempdir().unwrap();
        write_genome_params(dir.path(), "versionGenome\t2.7.4a\n");
        assert!(check_genome_version(dir.path(), "2.7.4a").is_ok());

        write_genome_params(dir.path(), "versionGenome\t2.7.10b\n");
        assert!(check_genome_version(dir.path(), "2.7.4a").is_ok());
    }

    #[test]
    fn version_genome_orders_by_component_not_lexically() {
        // A string comparison puts "2.7.10a" before "2.7.4a"; component
        // comparison must not.
        let older = parse_genome_version("2.7.4a").unwrap();
        let newer = parse_genome_version("2.7.10a").unwrap();
        assert!(newer > older, "2.7.10a must sort after 2.7.4a");
        assert!(parse_genome_version("2.7.4b").unwrap() > older);
    }

    #[test]
    fn version_genome_missing_line_is_accepted_with_a_warning() {
        let dir = tempfile::tempdir().unwrap();
        write_genome_params(dir.path(), "genomeType\tFull\n");
        assert!(check_genome_version(dir.path(), "2.7.4a").is_ok());
        // No genomeParameters.txt at all is also accepted.
        std::fs::remove_file(dir.path().join("genomeParameters.txt")).unwrap();
        assert!(check_genome_version(dir.path(), "2.7.4a").is_ok());
    }
}
