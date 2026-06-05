pub mod io;
pub mod packed_array;
pub mod packed_stream;
pub mod sa_build;
pub mod sa_index;
pub mod suffix_array;

use std::fs;
use std::io::BufWriter;
use std::path::Path;

use crate::error::Error;
use crate::genome::Genome;
use crate::index::packed_array::PackedArray;
use crate::index::packed_stream::PackedStreamWriter;
use crate::junction::SpliceJunctionDb;
use crate::junction::sjdb_insert::{self, PreparedJunction};
use crate::params::Parameters;
use crate::quant::transcriptome::TranscriptomeIndex;
use sa_index::SaIndex;
use suffix_array::SuffixArray;

/// Complete genome index (genome + suffix array + SA index + junction database).
#[derive(Clone)]
pub struct GenomeIndex {
    pub genome: Genome,
    pub suffix_array: SuffixArray,
    pub sa_index: SaIndex,
    pub junction_db: SpliceJunctionDb,
    /// Populated when the index was built with a GTF (`--sjdbGTFfile`).
    /// Mirrors STAR's `Transcriptome` object and is written to disk as
    /// `transcriptInfo.tab` + friends at `genomeGenerate`, reloaded at
    /// `alignReads` from the same files.
    pub transcriptome: Option<TranscriptomeIndex>,
    /// Prepared splice junctions in their post-dedup order (the same
    /// order they occupy in the Gsj buffer appended to the genome).
    /// Populated on both build and load paths when sjdb is present;
    /// empty for indices built without a GTF. Used at align time to
    /// decode Gsj-region SA hits back to real-genome `(donor, acceptor)`
    /// pairs.
    pub prepared_junctions: Vec<PreparedJunction>,
    /// `sjdbOverhang` recorded in `sjdbInfo.txt`. Zero when no sjdb
    /// junctions are present.
    pub sjdb_overhang: u32,
}

/// Output of [`GenomeIndex::build_prep`] — the shared setup
/// stage between the in-memory [`GenomeIndex::build`] and the
/// streaming [`GenomeIndex::generate_streaming`]. Holds everything
/// needed to drive the SA construction and the per-file writes.
struct BuildPrep {
    genome: Genome,
    junction_db: SpliceJunctionDb,
    transcriptome: Option<TranscriptomeIndex>,
    prepared_junctions: Vec<PreparedJunction>,
}

impl GenomeIndex {
    /// Build a complete genome index from FASTA files.
    pub fn build(params: &Parameters) -> Result<Self, Error> {
        let BuildPrep {
            genome,
            junction_db,
            transcriptome,
            prepared_junctions,
        } = Self::build_prep(params)?;

        log::info!("Building suffix array...");
        let suffix_array = SuffixArray::build(&genome)?;
        log::info!("Suffix array built: {} entries", suffix_array.len());

        log::info!("Building SA index...");
        let sa_index = SaIndex::build(&genome, &suffix_array, params.genome_sa_index_nbases)?;
        log::info!(
            "SA index built: nbases={}, {} indices",
            sa_index.nbases,
            sa_index.data.len()
        );

        let sjdb_overhang = if prepared_junctions.is_empty() {
            0
        } else {
            params.sjdb_overhang
        };

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

    /// Streaming genome-index generation: writes every index file
    /// directly to `params.genome_dir` without materialising the
    /// 25 GB-class SA `PackedArray` in RAM. The flow:
    ///
    /// 1. [`build_prep`][Self::build_prep] — same as the in-memory path:
    ///    load FASTA, parse GTF, build junction database, append
    ///    `gsj` to the genome.
    /// 2. Write genome files (`Genome`, `chrInfo`, `genomeParameters.txt`,
    ///    etc.) immediately — no dependency on the SA.
    /// 3. Open `genome_dir/SA` through a [`PackedStreamWriter`]; build
    ///    a [`SaIndexBuilder`][sa_index::SaIndexBuilder]. The caps-sa
    ///    emit callback feeds each entry to **both**, so the SA file
    ///    grows as construction progresses and the SAindex is built
    ///    on the fly. Total peak RSS during this phase ≈
    ///    `genome.sequence` + caps-sa scratch + `SaIndex.data` —
    ///    ~5 GB on the human genome vs the ~47 GB the in-memory
    ///    path peaked at.
    /// 4. Finalise the SA writer (flush partial-byte + padding).
    /// 5. **Build the SAindex in parallel** from the on-disk SA via
    ///    mmap + [`SaIndex::build_parallel`]. caps-sa's phase-4
    ///    emit is single-threaded; doing the SAindex k-mer work
    ///    there sat the entire ~16 min of serial extraction
    ///    on top of caps-sa's parallel SA build. Deferring lets
    ///    the SAindex use all available rayon workers.
    /// 6. Write the SAindex file, patch `genomeParameters.txt` with
    ///    the SA size, write transcriptome + sjdb outputs.
    ///
    /// This is the path used by `--runMode genomeGenerate`. The
    /// in-memory [`build`][Self::build] + [`write`][Self::write]
    /// remains for tests and any caller that needs random access to
    /// the SA in RAM.
    pub fn generate_streaming(params: &Parameters) -> Result<(), Error> {
        let BuildPrep {
            genome,
            junction_db: _,
            transcriptome,
            prepared_junctions,
        } = Self::build_prep(params)?;

        let dir = &params.genome_dir;
        if !dir.exists() {
            fs::create_dir_all(dir).map_err(|e| Error::io(e, dir))?;
        }

        log::info!("Writing genome files to {}...", dir.display());
        genome.write_index_files(dir, params)?;

        let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
        let gstrand_mask = (1u64 << gstrand_bit) - 1;
        let word_length = gstrand_bit + 1;
        let nbases = params.genome_sa_index_nbases;

        log::info!(
            "Streaming SA to {} (gstrand_bit={gstrand_bit}, word_length={word_length}, nbases={nbases})",
            dir.join("SA").display()
        );

        let sa_path = dir.join("SA");
        let sa_file = fs::File::create(&sa_path).map_err(|e| Error::io(e, &sa_path))?;
        // Large buffer — at 33 bit / entry × 5.9 B entries the SA file
        // is ~24 GB; an 8 MB BufWriter keeps the write-side syscall
        // rate to ~3000/s, dominated by sequential bandwidth.
        let sa_buf = BufWriter::with_capacity(8 * 1024 * 1024, sa_file);
        let mut sa_writer = PackedStreamWriter::new(sa_buf, word_length);

        log::info!("Building suffix array...");
        let (got_gbit, got_gmask, n_entries) =
            sa_build::build_streaming(&genome, params.temp_dir.as_deref(), |packed_value| {
                // Emit is now lightweight: just bit-pack into the SA
                // file. caps-sa's phase-4 emit loop is single-threaded,
                // so anything we do here serialises the whole build.
                // The SAindex is built afterwards via a parallel pass
                // over the on-disk SA.
                sa_writer
                    .write_one(packed_value)
                    .map_err(|e| Error::io(e, &sa_path))?;
                Ok(())
            })?;
        debug_assert_eq!(got_gbit, gstrand_bit);
        debug_assert_eq!(got_gmask, gstrand_mask);
        log::info!("Suffix array streamed to disk: {n_entries} entries");

        let buf = sa_writer.finish().map_err(|e| Error::io(e, &sa_path))?;
        buf.into_inner()
            .map_err(|e| Error::io(e.into_error(), &sa_path))?
            .sync_all()
            .map_err(|e| Error::io(e, &sa_path))?;

        let sa_size = PackedArray::data_byte_len_for(word_length, n_entries);

        // Parallel SAindex build via chunked pread on the on-disk
        // SA. Avoiding `memmap2::Mmap` keeps the touched SA pages
        // out of process RSS (`read_at` hits the kernel page cache,
        // which is kernel-side memory — not counted in
        // `Maximum resident set size`). The only big new allocation
        // for this phase is the ~2.86 GB `Vec<AtomicU64>` inside
        // `build_parallel`.
        log::info!("Building SAindex in parallel from on-disk SA...");
        let sa_handle = fs::File::open(&sa_path).map_err(|e| Error::io(e, &sa_path))?;
        let sa_index = SaIndex::build_parallel(
            &genome,
            &sa_handle,
            word_length,
            gstrand_bit,
            gstrand_mask,
            n_entries,
            nbases,
        )?;
        drop(sa_handle);
        log::info!(
            "SA index built: nbases={}, {} indices",
            sa_index.nbases,
            sa_index.data.len()
        );

        // Write SAindex, then drop it — the transcriptome / sjdb
        // writers below don't read it, and the file on disk is the
        // canonical copy from this point on. On the human genome
        // this frees ~1.5 GB of resident `PackedArray` before the
        // last writes.
        write_sa_index_file(&dir.join("SAindex"), &sa_index)?;
        drop(sa_index);

        // Update genomeParameters.txt with the SA file size — matches
        // STAR's `genomeFileSizes\t<n_genome> <sa_size>\n` pattern.
        // Same edit as `GenomeIndex::write` does for the in-memory
        // path; factored into a helper.
        update_genome_params_sa_size(&dir.join("genomeParameters.txt"), genome.n_genome, sa_size)?;

        // Transcriptome + sjdb files. Matches the tail of
        // `GenomeIndex::write` byte-for-byte.
        if let Some(tr) = &transcriptome {
            tr.write_transcript_info(dir)?;
            tr.write_exon_info(dir)?;
            tr.write_gene_info(dir)?;
            tr.write_exon_ge_tr_info(dir)?;
            tr.write_sjdb_list_from_gtf(dir, &genome)?;
            log::info!(
                "Wrote transcriptome index files: {} transcripts, {} genes",
                tr.n_transcripts(),
                tr.gene_ids.len()
            );
        }
        if !prepared_junctions.is_empty() {
            sjdb_insert::write_sjdb_info_tab(
                &dir.join("sjdbInfo.txt"),
                &prepared_junctions,
                params.sjdb_overhang,
            )?;
            sjdb_insert::write_sjdb_list_out_tab(
                &dir.join("sjdbList.out.tab"),
                &prepared_junctions,
                &genome,
            )?;
            log::info!("Wrote sjdb files: {} junctions", prepared_junctions.len());
        }

        Ok(())
    }

    /// Shared setup: load FASTA, parse GTF, build junctions, append
    /// `gsj` to the genome. Used by both [`build`][Self::build] and
    /// [`generate_streaming`][Self::generate_streaming] so they
    /// share the same SA-input shape.
    fn build_prep(params: &Parameters) -> Result<BuildPrep, Error> {
        log::info!("Loading FASTA files...");
        let mut genome = Genome::from_fasta(params)?;

        log::info!(
            "Loaded {} chromosomes, total padded genome size: {} bytes",
            genome.n_chr_real,
            genome.n_genome
        );

        // Parse GTF once and share the result between the junction database,
        // the transcriptome index, and the sjdb insertion pipeline. Junction
        // preparation + Gsj append must happen BEFORE the suffix array is
        // built, because STAR indexes the flanking-sequence buffer alongside
        // the real genome in a single SA (`sjdbBuildIndex.cpp:293`).
        let (junction_db, transcriptome, prepared_junctions) = if let Some(ref gtf_path) =
            params.sjdb_gtf_file
        {
            let n_genome_real = genome.n_genome;

            let exons = crate::junction::gtf::parse_gtf_configured(
                gtf_path,
                &params.sjdb_gtf_feature_exon,
                &params.sjdb_gtf_chr_prefix,
            )?;
            log::debug!("Parsed {} exon features from GTF", exons.len());

            let tr = TranscriptomeIndex::from_gtf_exons_configured(
                &exons,
                &genome,
                &params.sjdb_gtf_tag_exon_parent_transcript,
                &params.sjdb_gtf_tag_exon_parent_gene,
            )?;
            log::info!(
                "Transcriptome index built from GTF: {} transcripts, {} genes",
                tr.n_transcripts(),
                tr.gene_ids.len()
            );

            let raw = crate::junction::gtf::extract_junctions_configured(
                exons,
                &genome,
                &params.sjdb_gtf_tag_exon_parent_transcript,
            )?;
            log::info!("Extracted {} annotated junctions from GTF", raw.len());
            let jdb = SpliceJunctionDb::from_raw_junctions(&raw);

            let prepared: Vec<PreparedJunction> = raw
                .iter()
                .map(|&(chr_idx, intron_start, intron_end, strand)| {
                    sjdb_insert::prepare_junction(
                        chr_idx,
                        intron_start,
                        intron_end,
                        strand,
                        &genome,
                        n_genome_real,
                    )
                })
                .collect();
            let prepared = sjdb_insert::sort_and_dedup(prepared);

            let gsj =
                sjdb_insert::build_gsj(&prepared, &genome, n_genome_real, params.sjdb_overhang)?;
            log::info!(
                "Built Gsj buffer: {} junctions × {} bytes = {} bytes",
                prepared.len(),
                2 * params.sjdb_overhang + 1,
                gsj.len()
            );
            genome.append_sjdb(&gsj);
            log::info!(
                "Extended genome with sjdb: n_genome = {} (pre-sjdb {})",
                genome.n_genome,
                n_genome_real
            );

            (jdb, Some(tr), prepared)
        } else {
            log::info!("No GTF file provided, all junctions will be novel");
            (SpliceJunctionDb::empty(), None, Vec::new())
        };

        log::info!(
            "Junction database initialized: {} annotated junctions",
            junction_db.len()
        );

        Ok(BuildPrep {
            genome,
            junction_db,
            transcriptome,
            prepared_junctions,
        })
    }

    /// Convert a raw SA position for a reverse-strand match to forward genome coordinates.
    ///
    /// The SA stores reverse-strand positions as offsets within the RC genome region.
    /// For chromosome identification and SAM output, we need the corresponding position
    /// in the forward genome (leftmost aligned base in forward coordinates).
    ///
    /// For forward strand: returns the position unchanged.
    /// For reverse strand: `forward_pos = n_genome - 1 - sa_pos - (match_length - 1)`
    ///
    /// The raw SA position is still needed for genome base access (add n_genome offset).
    pub fn sa_pos_to_forward(&self, sa_pos: u64, is_reverse: bool, match_length: usize) -> u64 {
        if is_reverse {
            self.genome
                .n_genome
                .saturating_sub(sa_pos)
                .saturating_sub(match_length as u64)
        } else {
            sa_pos
        }
    }

    /// Write index files to directory.
    pub fn write(&self, dir: &Path, params: &Parameters) -> Result<(), Error> {
        // Write genome files
        self.genome.write_index_files(dir, params)?;

        // Write SA file
        let sa_path = dir.join("SA");
        fs::write(&sa_path, self.suffix_array.data.data()).map_err(|e| Error::io(e, &sa_path))?;

        // Write SAindex file (factored helper, shared with the
        // streaming path).
        write_sa_index_file(&dir.join("SAindex"), &self.sa_index)?;

        // Update genomeParameters.txt with SA file size.
        let sa_size = self.suffix_array.data.data().len();
        update_genome_params_sa_size(
            &dir.join("genomeParameters.txt"),
            self.genome.n_genome,
            sa_size,
        )?;

        // Write transcriptome index files (STAR-compatible) when the GTF
        // was supplied. Matches STAR's `GTF_transcriptGeneSJ.cpp` outputs.
        if let Some(tr) = &self.transcriptome {
            tr.write_transcript_info(dir)?;
            tr.write_exon_info(dir)?;
            tr.write_gene_info(dir)?;
            tr.write_exon_ge_tr_info(dir)?;
            tr.write_sjdb_list_from_gtf(dir, &self.genome)?;
            log::info!(
                "Wrote transcriptome index files: {} transcripts, {} genes",
                tr.n_transcripts(),
                tr.gene_ids.len()
            );
        }

        // Write sjdbInfo.txt + sjdbList.out.tab — the sjdb-insertion outputs
        // STAR emits alongside the transcriptome files when junctions are
        // baked into the genome at `genomeGenerate` time.
        if !self.prepared_junctions.is_empty() {
            sjdb_insert::write_sjdb_info_tab(
                &dir.join("sjdbInfo.txt"),
                &self.prepared_junctions,
                params.sjdb_overhang,
            )?;
            sjdb_insert::write_sjdb_list_out_tab(
                &dir.join("sjdbList.out.tab"),
                &self.prepared_junctions,
                &self.genome,
            )?;
            log::info!(
                "Wrote sjdb files: {} junctions",
                self.prepared_junctions.len()
            );
        }

        Ok(())
    }
}

/// Write a [`SaIndex`] to `path` in STAR's `SAindex` file format
/// (8-byte little-endian `nbases` header, then the
/// `genome_sa_index_start[]` u64s, then the packed entries).
/// Shared between [`GenomeIndex::write`] (in-memory build) and
/// [`GenomeIndex::generate_streaming`] (streaming build).
fn write_sa_index_file(path: &Path, sa_index: &SaIndex) -> Result<(), Error> {
    use std::io::Write;
    let mut f = fs::File::create(path).map_err(|e| Error::io(e, path))?;
    f.write_all(&(sa_index.nbases as u64).to_le_bytes())
        .map_err(|e| Error::io(e, path))?;
    for &val in &sa_index.genome_sa_index_start {
        f.write_all(&val.to_le_bytes())
            .map_err(|e| Error::io(e, path))?;
    }
    f.write_all(sa_index.data.data())
        .map_err(|e| Error::io(e, path))?;
    Ok(())
}

/// Update `genomeParameters.txt` with the actual SA file size.
/// [`Genome::write_genome_parameters_txt`] emits the file with a `0`
/// placeholder for the SA size; this helper substitutes the real
/// number once it's known. Matches STAR's
/// `genomeFileSizes\t<n_genome> <sa_size>\n` line format (tab before
/// the first value, space before the second).
fn update_genome_params_sa_size(path: &Path, n_genome: u64, sa_size: usize) -> Result<(), Error> {
    let content = fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
    let updated = content.replace(
        &format!("genomeFileSizes\t{n_genome} 0"),
        &format!("genomeFileSizes\t{n_genome} {sa_size}"),
    );
    fs::write(path, updated).map_err(|e| Error::io(e, path))?;
    Ok(())
}
