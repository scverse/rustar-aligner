pub mod io;
pub mod packed_array;
pub mod sa_index;
pub mod suffix_array;

use std::fs;
use std::path::Path;

use crate::error::Error;
use crate::genome::Genome;
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
    /// Prepared splice junctions (sorted/deduped) used only on the build
    /// path to write `sjdbInfo.txt` + `sjdbList.out.tab`. Empty on the
    /// load path — those files have already been written and are not
    /// needed at align time.
    pub prepared_junctions: Vec<PreparedJunction>,
}

impl GenomeIndex {
    /// Build a complete genome index from FASTA files.
    pub fn build(params: &Parameters) -> Result<Self, Error> {
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

            // `extract_junctions_from_exons` returns chromosome-local 1-based
            // intron coordinates (matching STAR's `sjdbList.fromGTF.out.tab`).
            // `prepare_junction` + `sjdbInfo.txt` expect 0-based absolute
            // genome offsets (matching STAR's `sjdbPrepare.cpp` and
            // `detect_splice_motif`'s `genome.sequence[donor_pos]` access),
            // so convert here.
            let prepared: Vec<PreparedJunction> = raw
                .iter()
                .map(|&(chr_idx, start_local_1b, end_local_1b, strand)| {
                    let chr_off = genome.chr_start[chr_idx];
                    sjdb_insert::prepare_junction(
                        chr_idx,
                        chr_off + start_local_1b - 1,
                        chr_off + end_local_1b - 1,
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

        Ok(GenomeIndex {
            genome,
            suffix_array,
            sa_index,
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
        use std::io::Write;

        // Write genome files
        self.genome.write_index_files(dir, params)?;

        // Write SA file
        let sa_path = dir.join("SA");
        fs::write(&sa_path, self.suffix_array.data.data()).map_err(|e| Error::io(e, &sa_path))?;

        // Write SAindex file
        let sai_path = dir.join("SAindex");
        let mut sai_file = fs::File::create(&sai_path).map_err(|e| Error::io(e, &sai_path))?;

        // Write header: gSAindexNbases as u64
        sai_file
            .write_all(&(self.sa_index.nbases as u64).to_le_bytes())
            .map_err(|e| Error::io(e, &sai_path))?;

        // Write genomeSAindexStart array
        for &val in &self.sa_index.genome_sa_index_start {
            sai_file
                .write_all(&val.to_le_bytes())
                .map_err(|e| Error::io(e, &sai_path))?;
        }

        // Write packed SAindex data
        sai_file
            .write_all(self.sa_index.data.data())
            .map_err(|e| Error::io(e, &sai_path))?;

        // Update genomeParameters.txt with SA file size. Matches STAR's
        // `genomeFileSizes\t<n_genome> <sa_size>\n` pattern (tab before first
        // value, space between subsequent values) — written out in
        // Genome::write_genome_parameters_txt with `0` as the SA placeholder.
        let genome_params_path = dir.join("genomeParameters.txt");
        let sa_size = self.suffix_array.data.data().len();
        let content = fs::read_to_string(&genome_params_path)
            .map_err(|e| Error::io(e, &genome_params_path))?;
        let updated_content = content.replace(
            &format!("genomeFileSizes\t{} 0", self.genome.n_genome),
            &format!("genomeFileSizes\t{} {}", self.genome.n_genome, sa_size),
        );
        fs::write(&genome_params_path, updated_content)
            .map_err(|e| Error::io(e, &genome_params_path))?;

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
