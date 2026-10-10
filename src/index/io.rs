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
        check_genome_version(genome_dir, params)?;

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

        // Align-time annotations (STAR's on-the-fly sjdbInsertJunctions):
        // collected here, inserted into the genome once the index is assembled.
        // STAR priorities: GTF 20, sjdbFileChrStartEnd 10.
        let mut otf_raw: Vec<(usize, u64, u64, u8)> = Vec::new();
        let mut otf_file_raw: Vec<(usize, u64, u64, u8)> = Vec::new();
        if let Some(ref gtf_path) = params.sjdb_gtf_file {
            let exons = crate::junction::gtf::parse_gtf_configured(
                gtf_path,
                &params.sjdb_gtf_feature_exon,
                &params.sjdb_gtf_chr_prefix,
            )?;
            otf_raw.extend(crate::junction::gtf::extract_junctions_configured(
                exons,
                &genome,
                &params.sjdb_gtf_tag_exon_parent_transcript,
            )?);
        }
        if !params.sjdb_file_chr_start_end.is_empty() {
            otf_file_raw.extend(crate::junction::chr_start_end::parse_sjdb_chr_start_end(
                &params.sjdb_file_chr_start_end,
                &genome,
            )?);
        }

        // The annotated-junction database consulted at stitch time comes from
        // the index's junctions (sjdbInfo.txt), keyed on STAR's stored
        // coordinates and carrying sjdbMotif/sjdbShift*/sjdbStrand. Without it
        // the standard workflow (annotated index, `--genomeDir` only at
        // mapping) would treat every junction as novel.
        let junction_db = if prepared_junctions.is_empty() {
            SpliceJunctionDb::empty()
        } else {
            log::info!(
                "Loaded {} annotated junctions from index sjdbInfo.txt",
                prepared_junctions.len()
            );
            SpliceJunctionDb::from_prepared(&prepared_junctions)
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

        let mut index = GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db,
            transcriptome,
            prepared_junctions,
            sjdb_overhang,
        };
        if !otf_raw.is_empty() || !otf_file_raw.is_empty() {
            log::info!(
                "Inserting {} align-time annotated junctions into the genome",
                otf_raw.len() + otf_file_raw.len()
            );
            let mut prepared = index.prepare_raw_junctions(&otf_raw, 20);
            prepared.extend(index.prepare_raw_junctions(&otf_file_raw, 10));
            crate::index::sjdb_otf::insert_junctions(&mut index, prepared, params.sjdb_overhang)?;
        }
        Ok(index)
    }
}

/// STAR's genome compatibility checks (`Genome_genomeLoad.cpp:56-101`): the
/// index's `versionGenome` must be the one this STAR release writes (2.7.4a),
/// and an annotated index from before `sjdbInsertSave` existed cannot take
/// junctions inserted at mapping time.
fn check_genome_version(genome_dir: &Path, params: &Parameters) -> Result<(), Error> {
    const VERSION_GENOME: &str = "2.7.4a";
    let path = genome_dir.join("genomeParameters.txt");
    let contents = std::fs::read_to_string(&path).map_err(|_| {
        Error::Index(format!(
            "EXITING because of FATAL ERROR: could not open genome file {}\n\
             SOLUTION: check that the path to genome files, specified in --genomeDir is \
             correct and the files are present, and have user read permsissions",
            path.display()
        ))
    })?;
    let value = |key: &str| {
        contents.lines().find_map(|l| {
            let mut f = l.split_whitespace();
            (f.next() == Some(key)).then(|| f.next().unwrap_or("").to_string())
        })
    };
    match value("versionGenome") {
        None => {
            return Err(Error::Index(
                "EXITING because of FATAL ERROR: read no value for the versionGenome parameter \
                 from genomeParameters.txt file\n\
                 SOLUTION: please re-generate genome from scratch with the latest version of STAR"
                    .to_string(),
            ));
        }
        Some(v) if v != VERSION_GENOME => {
            return Err(Error::Index(format!(
                "EXITING because of FATAL ERROR: Genome version: {v} is INCOMPATIBLE with \
                 running STAR version: 2.7.11b\n\
                 SOLUTION: please re-generate genome from scratch with running version of STAR, \
                 or with version: {VERSION_GENOME}"
            )));
        }
        Some(_) => {}
    }
    let sjdb_insert = params.sjdb_gtf_file.is_some()
        || !params.sjdb_file_chr_start_end.is_empty()
        || params.twopass_mode == crate::params::TwopassMode::Basic;
    if sjdb_insert && genome_dir.join("sjdbInfo.txt").exists() && value("sjdbInsertSave").is_none()
    {
        return Err(Error::Index(
            "EXITING because of FATAL ERROR: old Genome is INCOMPATIBLE with on the fly junction \
             insertion\n\
             SOLUTION: please re-generate genome from scratch with the latest version of STAR"
                .to_string(),
        ));
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

    /// STAR refuses an index written for another genome format version.
    #[test]
    fn rejects_incompatible_genome_version() {
        let dir = tempfile::tempdir().unwrap();
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "r.fq"]);
        std::fs::write(
            dir.path().join("genomeParameters.txt"),
            "versionGenome\t2.7.1a\n",
        )
        .unwrap();
        let err = check_genome_version(dir.path(), &params).unwrap_err();
        assert!(err.to_string().contains("INCOMPATIBLE"), "{err}");
        std::fs::write(
            dir.path().join("genomeParameters.txt"),
            "versionGenome\t2.7.4a\n",
        )
        .unwrap();
        assert!(check_genome_version(dir.path(), &params).is_ok());
    }

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
}
