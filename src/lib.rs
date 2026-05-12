#![allow(non_snake_case)]

pub mod error;
pub mod params;

pub mod align;
pub mod chimeric;
pub mod cpu;
pub mod genome;
pub mod index;
pub mod io;
pub mod junction;
pub mod mapq;
pub mod quant;
pub mod stats;

use log::info;

use crate::params::{Parameters, RunMode};

/// Top-level dispatcher. Called from `main()` after CLI parsing.
pub fn run(params: &Parameters) -> anyhow::Result<()> {
    params.validate()?;

    info!("rustar-aligner {}", env!("CARGO_PKG_VERSION"));
    info!("{}", env!("VERSION_BODY"));
    info!("{}", cpu::cpu_detected_line());
    if let Some(hint) = cpu::upgrade_hint() {
        info!("{hint}");
    }
    info!("runMode: {}", params.run_mode);
    info!("runThreadN: {}", params.run_thread_n);

    match params.run_mode {
        RunMode::GenomeGenerate => genome_generate(params),
        RunMode::AlignReads => align_reads(params),
    }
}

fn genome_generate(params: &Parameters) -> anyhow::Result<()> {
    use index::GenomeIndex;

    info!("genomeDir: {}", params.genome_dir.display());
    info!(
        "genomeFastaFiles: {:?}",
        params
            .genome_fasta_files
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
    );

    info!("Building genome index...");
    let index = GenomeIndex::build(params)?;

    info!("Writing index files to {}...", params.genome_dir.display());
    index.write(&params.genome_dir, params)?;

    info!("Genome generation complete!");
    Ok(())
}

/// Trait for alignment output writers (SAM or BAM).
/// `finish` flushes/sorts/closes the output; default is a no-op for streaming writers.
trait AlignmentWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error>;
    fn finish(&mut self) -> Result<(), error::Error> {
        Ok(())
    }
}

/// Null writer that discards all output (for two-pass mode pass 1)
struct NullWriter;

impl AlignmentWriter for NullWriter {
    fn write_batch(
        &mut self,
        _batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        Ok(()) // Discard all records
    }
}

impl AlignmentWriter for crate::io::sam::SamWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
}

impl AlignmentWriter for crate::io::bam::BamWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
    fn finish(&mut self) -> Result<(), error::Error> {
        self.finish()
    }
}

impl AlignmentWriter for crate::io::bam::SortedBamWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
    fn finish(&mut self) -> Result<(), error::Error> {
        self.finish()
    }
}

impl AlignmentWriter for crate::io::sam::SamStdoutWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
}

impl AlignmentWriter for crate::io::bam::BamStdoutWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
    fn finish(&mut self) -> Result<(), error::Error> {
        self.finish()
    }
}

impl AlignmentWriter for crate::io::bam::SortedBamStdoutWriter {
    fn write_batch(
        &mut self,
        batch: &[noodles::sam::alignment::record_buf::RecordBuf],
    ) -> Result<(), error::Error> {
        self.write_batch(batch)
    }
    fn finish(&mut self) -> Result<(), error::Error> {
        self.finish()
    }
}

fn align_reads(params: &Parameters) -> anyhow::Result<()> {
    use crate::index::GenomeIndex;

    use crate::params::TwopassMode;

    use std::sync::Arc;

    let time_start = chrono::Local::now();

    info!("Starting read alignment...");

    // Configure Rayon thread pool based on --runThreadN
    if params.run_thread_n > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(params.run_thread_n)
            .build_global()
            .map_err(|e| {
                error::Error::Parameter(format!("Failed to configure thread pool: {}", e))
            })?;
        info!("Using {} threads for alignment", params.run_thread_n);
    } else {
        info!("Using single-threaded mode");
    }

    // Validate read files
    if params.read_files_in.is_empty() {
        anyhow::bail!("No read files specified (--readFilesIn)");
    }

    // 1. Load genome index
    info!("Loading genome index from {}", params.genome_dir.display());
    let index = Arc::new(GenomeIndex::load(&params.genome_dir, params)?);
    info!(
        "Loaded {} chromosomes, {} bases",
        index.genome.n_chr_real, index.genome.n_genome
    );

    // Redefine window parameters based on genome size (STAR's Genome_genomeLoad.cpp)
    let mut params = params.clone();
    params.redefine_window_params(index.genome.n_genome);

    // Build gene-count context if --quantMode GeneCounts was requested.
    // GTF requirement is already validated in params.validate().
    let quant_ctx: Option<std::sync::Arc<crate::quant::QuantContext>> =
        if params.quant_gene_counts() {
            let gtf_path = params.sjdb_gtf_file.as_ref().unwrap();
            info!(
                "quantMode GeneCounts: building gene annotation from {}",
                gtf_path.display()
            );
            let ctx = crate::quant::QuantContext::build(
                gtf_path,
                &index.genome,
                &params.sjdb_gtf_feature_exon,
                &params.sjdb_gtf_chr_prefix,
                &params.sjdb_gtf_tag_exon_parent_gene,
            )?;
            Some(std::sync::Arc::new(ctx))
        } else {
            None
        };

    // Use the transcriptome index loaded alongside the genome (populated
    // from transcriptInfo.tab / exonInfo.tab / geneInfo.tab at load time
    // — see GenomeIndex::load). Only wire it through to the pipeline when
    // `--quantMode TranscriptomeSAM` is requested.
    let tr_idx: Option<std::sync::Arc<crate::quant::transcriptome::TranscriptomeIndex>> =
        if params.quant_transcriptome_sam() {
            let tr = index.transcriptome.as_ref().ok_or_else(|| {
                anyhow::anyhow!(
                    "--quantMode TranscriptomeSAM requires a GTF-aware index; \
                     re-run genomeGenerate with --sjdbGTFfile or pass --sjdbGTFfile \
                     at alignReads so transcriptInfo.tab can be (re)built"
                )
            })?;
            info!(
                "quantMode TranscriptomeSAM: using {} transcripts from genome index",
                tr.n_transcripts()
            );
            Some(std::sync::Arc::new(tr.clone()))
        } else {
            None
        };

    let time_map_start = chrono::Local::now();

    // 2. Dispatch based on two-pass mode
    let stats = match params.twopass_mode {
        TwopassMode::None => {
            info!("Running single-pass alignment");
            run_single_pass(&index, &params, quant_ctx.as_ref(), tr_idx.as_ref())?
        }
        TwopassMode::Basic => {
            info!("Running two-pass alignment mode");
            run_two_pass(&index, &params, quant_ctx.as_ref(), tr_idx.as_ref())?
        }
    };

    let time_finish = chrono::Local::now();

    // Write Log.final.out
    let log_path = params.out_file_name_prefix.join("Log.final.out");
    if let Some(parent) = log_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    stats.write_log_final(&log_path, time_start, time_map_start, time_finish)?;
    info!("Wrote {}", log_path.display());

    // Write ReadsPerGene.out.tab if quantMode GeneCounts was requested.
    if let Some(ref ctx) = quant_ctx {
        let quant_path = params.out_file_name_prefix.join("ReadsPerGene.out.tab");
        ctx.counts.write_output(&quant_path, &ctx.gene_ann)?;
        info!("Wrote {}", quant_path.display());
    }

    info!("Alignment complete!");
    Ok(())
}

/// Run single-pass alignment (original logic)
fn run_single_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
    quant_ctx: Option<&std::sync::Arc<crate::quant::QuantContext>>,
    tr_idx: Option<&std::sync::Arc<crate::quant::transcriptome::TranscriptomeIndex>>,
) -> anyhow::Result<std::sync::Arc<crate::stats::AlignmentStats>> {
    use crate::io::bam::{BamWriter, SortedBamWriter};
    use crate::io::sam::SamWriter;
    use crate::params::OutSamFormat;
    use std::sync::Arc;

    // Initialize statistics collectors
    let stats = Arc::new(crate::stats::AlignmentStats::new());
    let sj_stats = Arc::new(crate::junction::SpliceJunctionStats::new());

    // Clone the quant Arc so each dispatch call can own a reference.
    let quant = quant_ctx.map(Arc::clone);
    let tr = tr_idx.map(Arc::clone);

    // Open transcriptome BAM writer if requested.
    let mut tr_writer: Option<BamWriter> = if let Some(tidx) = tr.as_ref() {
        let path = params
            .out_file_name_prefix
            .join("Aligned.toTranscriptome.out.bam");
        info!("Writing transcriptome BAM to {}", path.display());
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        Some(BamWriter::create_transcriptome(&path, tidx, params)?)
    } else {
        None
    };

    // Create unmapped FASTQ writers if --outReadsUnmapped Fastx
    use crate::io::fastq::UnmappedFastqWriter;
    use crate::params::OutReadsUnmapped;

    let is_paired = params.read_files_in.len() == 2;
    let mut unmapped_w1: Option<UnmappedFastqWriter> =
        if params.out_reads_unmapped == OutReadsUnmapped::Fastx {
            let path = params.out_file_name_prefix.join("Unmapped.out.mate1");
            info!("Writing unmapped reads to {}", path.display());
            Some(UnmappedFastqWriter::create(&path)?)
        } else {
            None
        };
    let mut unmapped_w2: Option<UnmappedFastqWriter> =
        if params.out_reads_unmapped == OutReadsUnmapped::Fastx && is_paired {
            let path = params.out_file_name_prefix.join("Unmapped.out.mate2");
            info!("Writing unmapped mate2 reads to {}", path.display());
            Some(UnmappedFastqWriter::create(&path)?)
        } else {
            None
        };

    // 4. Route to SAM or BAM output based on --outSAMtype / --outStd
    use crate::params::{OutSamSortOrder, OutStd};

    let out_type = params
        .out_sam_type()
        .map_err(|e| anyhow::anyhow!("Invalid --outSAMtype: {}", e))?;

    // Build boxed writer — stdout takes precedence over file output.
    let mut writer: Box<dyn AlignmentWriter> = match params.out_std {
        OutStd::Sam => {
            info!("Writing SAM to stdout (--outStd SAM)");
            Box::new(crate::io::sam::SamStdoutWriter::create(
                &index.genome,
                params,
            )?)
        }
        OutStd::BamUnsorted => {
            info!("Writing unsorted BAM to stdout (--outStd BAM_Unsorted)");
            Box::new(crate::io::bam::BamStdoutWriter::create(
                &index.genome,
                params,
            )?)
        }
        OutStd::BamSortedByCoordinate => {
            info!("Writing coordinate-sorted BAM to stdout (--outStd BAM_SortedByCoordinate)");
            Box::new(crate::io::bam::SortedBamStdoutWriter::create(
                &index.genome,
                params,
            )?)
        }
        OutStd::None => match out_type.format {
            OutSamFormat::Sam => {
                let output_path = params.out_file_name_prefix.join("Aligned.out.sam");
                info!("Writing SAM to {}", output_path.display());
                if let Some(parent) = output_path.parent() {
                    std::fs::create_dir_all(parent)?;
                }
                Box::new(SamWriter::create(&output_path, &index.genome, params)?)
            }
            OutSamFormat::Bam => {
                let sorted = out_type.sort_order == Some(OutSamSortOrder::SortedByCoordinate);
                let output_path = if sorted {
                    params
                        .out_file_name_prefix
                        .join("Aligned.sortedByCoord.out.bam")
                } else {
                    params.out_file_name_prefix.join("Aligned.out.bam")
                };
                info!("Writing BAM to {}", output_path.display());
                if let Some(parent) = output_path.parent() {
                    std::fs::create_dir_all(parent)?;
                }
                if sorted {
                    Box::new(SortedBamWriter::create(
                        &output_path,
                        &index.genome,
                        params,
                    )?)
                } else {
                    Box::new(BamWriter::create(&output_path, &index.genome, params)?)
                }
            }
            OutSamFormat::None => {
                anyhow::bail!("Output format 'None' not yet implemented");
            }
        },
    };

    // Align reads through the boxed writer.
    match params.read_files_in.len() {
        1 => align_reads_single_end(
            params,
            index,
            writer.as_mut(),
            &stats,
            &sj_stats,
            quant.as_ref(),
            tr.as_ref(),
            tr_writer.as_mut(),
            unmapped_w1.as_mut(),
        ),
        2 => align_reads_paired_end(
            params,
            index,
            writer.as_mut(),
            &stats,
            &sj_stats,
            quant.as_ref(),
            tr.as_ref(),
            tr_writer.as_mut(),
            unmapped_w1.as_mut(),
            unmapped_w2.as_mut(),
        ),
        n => anyhow::bail!("Invalid number of read files: {} (expected 1 or 2)", n),
    }?;

    writer.finish()?;

    // Flush transcriptome BAM.
    if let Some(ref mut w) = tr_writer {
        w.finish()?;
    }

    // 5. Write SJ.out.tab file
    let sj_output_path = params.out_file_name_prefix.join("SJ.out.tab");
    if !sj_stats.is_empty() {
        info!(
            "Writing splice junction statistics to {}",
            sj_output_path.display()
        );
        sj_stats.write_output(&sj_output_path, &index.genome, params)?;
    }

    // 6. Print summary
    stats.print_summary();

    Ok(stats)
}

/// Run two-pass alignment mode
fn run_two_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
    quant_ctx: Option<&std::sync::Arc<crate::quant::QuantContext>>,
    tr_idx: Option<&std::sync::Arc<crate::quant::transcriptome::TranscriptomeIndex>>,
) -> anyhow::Result<std::sync::Arc<crate::stats::AlignmentStats>> {
    use std::sync::Arc;

    // PASS 1: Junction discovery (no quant counting in pass 1)
    info!("Two-pass mode: Pass 1 - Junction discovery");
    let (sj_stats_pass1, novel_junctions) = run_pass1(index, params)?;

    // Write SJ.pass1.out.tab
    let pass1_path = params.out_file_name_prefix.join("SJ.pass1.out.tab");

    // Create output directory if it doesn't exist
    if let Some(parent) = pass1_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    info!("Writing pass 1 junctions to {}", pass1_path.display());
    sj_stats_pass1.write_output(&pass1_path, &index.genome, params)?;
    info!(
        "Pass 1 discovered {} novel junctions",
        novel_junctions.len()
    );

    // Insert novel junctions into DB
    let mut merged_index = (**index).clone();
    merged_index
        .junction_db
        .insert_novel(novel_junctions.clone());
    info!(
        "Merged junction DB: {} total junctions",
        merged_index.junction_db.len()
    );

    // PASS 2: Re-alignment with merged DB (quant counts happen here)
    info!("Two-pass mode: Pass 2 - Re-alignment");
    let stats = run_single_pass(&Arc::new(merged_index), params, quant_ctx, tr_idx)?;

    Ok(stats)
}

/// Run pass 1 of two-pass mode (junction discovery)
fn run_pass1(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<(
    crate::junction::SpliceJunctionStats,
    Vec<(
        crate::junction::NovelJunctionKey,
        crate::junction::JunctionInfo,
    )>,
)> {
    use std::sync::Arc;

    let stats = Arc::new(crate::stats::AlignmentStats::new());
    let sj_stats = Arc::new(crate::junction::SpliceJunctionStats::new());

    // Modify params to limit reads for pass 1
    let mut params_pass1 = params.clone();
    if params.twopass1_reads_n >= 0 {
        params_pass1.read_map_number = params.twopass1_reads_n;
        info!("Pass 1 will align {} reads", params.twopass1_reads_n);
    } else {
        info!("Pass 1 will align all reads");
    }

    // Create NullWriter (discard SAM/BAM output in pass 1)
    let mut null_writer = NullWriter;

    // Align reads (single-end or paired-end); no quant counting in pass 1
    match params.read_files_in.len() {
        1 => align_reads_single_end(
            &params_pass1,
            index,
            &mut null_writer,
            &stats,
            &sj_stats,
            None,
            None,
            None,
            None,
        )?,
        2 => align_reads_paired_end(
            &params_pass1,
            index,
            &mut null_writer,
            &stats,
            &sj_stats,
            None,
            None,
            None,
            None,
            None,
        )?,
        n => anyhow::bail!("Invalid number of read files: {} (expected 1 or 2)", n),
    }

    info!("Pass 1 aligned {} reads", stats.total_reads());

    // Filter novel junctions
    let novel_junctions = crate::junction::filter_novel_junctions(&sj_stats, params);

    // Return ownership of sj_stats
    let sj_stats = Arc::try_unwrap(sj_stats).unwrap_or_else(|arc| (*arc).clone());

    Ok((sj_stats, novel_junctions))
}

/// Reverse-complement an encoded read (A=0,C=1,G=2,T=3,N=4).  Shared by the
/// SE and PE transcriptome builders for the STAR `Read1[2]` soft-clip path.
fn rc_encode(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| crate::io::fastq::complement_base(b))
        .collect()
}

/// Pick a random primary-hit index and compute the MAPQ for a set of
/// transcriptome projections.  Shared by the SE and PE builders.
fn pick_primary_and_mapq(
    n_alignments: usize,
    n_for_mapq: usize,
    read_name: &str,
    params: &Parameters,
) -> (usize, u8) {
    use crate::align::read_align::per_read_seed;
    use crate::mapq::calculate_mapq;
    use rand::{Rng, SeedableRng, rngs::StdRng};

    let mut rng = StdRng::seed_from_u64(per_read_seed(params.run_rng_seed, read_name));
    let primary_hit = rng.gen_range(0..n_alignments);
    let mapq = calculate_mapq(n_alignments.max(n_for_mapq), params.out_sam_mapq_unique);
    (primary_hit, mapq)
}

/// Per-read metadata for BySJout disk-buffered mode.
/// SAM records are written to a temp file; only small metadata stays in memory.
struct BySJReadMeta {
    /// Number of SAM records written to the temp file for this read.
    n_sam_records: u32,
    /// Junction keys from primary alignment. Empty if unmapped or no junctions.
    junction_keys: Vec<crate::junction::SjKey>,
    /// Chimeric alignments — kept in memory because they're rare (~0.1% of reads).
    chimeric_alns: Vec<crate::chimeric::ChimericAlignment>,
    /// Transcriptome SAM records (kept in memory — optional feature).
    transcriptome_records: Vec<noodles::sam::alignment::record_buf::RecordBuf>,
}

/// Helper struct to hold alignment results from parallel processing
struct AlignmentBatchResults {
    sam_records: crate::io::sam::BufferedSamRecords,
    chimeric_alns: Vec<crate::chimeric::ChimericAlignment>,
    /// Junction keys from the primary (best) alignment for BySJout filtering.
    /// Empty if unmapped or no junctions.
    primary_junction_keys: Vec<crate::junction::SjKey>,
    /// Transcriptome-space SAM records for `--quantMode TranscriptomeSAM`.
    /// Empty unless that mode is enabled.
    transcriptome_records: Vec<noodles::sam::alignment::record_buf::RecordBuf>,
    /// Unmapped reads for `--outReadsUnmapped Fastx` (name, encoded_seq, qual).
    /// mate1 file (also used for SE). Empty unless that mode is enabled.
    unmapped_mate1: Vec<(String, Vec<u8>, Vec<u8>)>,
    /// Unmapped mate2 reads (PE only). Empty unless outReadsUnmapped=Fastx.
    unmapped_mate2: Vec<(String, Vec<u8>, Vec<u8>)>,
}

/// Build transcriptome-space records for a single-end read.  Projects every
/// surviving genome-space alignment onto all compatible transcripts, picks one
/// projected alignment at random as the primary (seeded by `per_read_seed`),
/// and emits SAM records with the transcriptome header.
#[allow(clippy::too_many_arguments)]
fn build_transcriptome_records_se(
    transcripts: &[crate::align::transcript::Transcript],
    read_name: &str,
    read_seq: &[u8],
    read_qual: &[u8],
    genome: &crate::genome::Genome,
    tr_idx: &crate::quant::transcriptome::TranscriptomeIndex,
    params: &Parameters,
    n_for_mapq: usize,
) -> Result<Vec<noodles::sam::alignment::record_buf::RecordBuf>, error::Error> {
    use crate::io::sam::SamWriter;
    use crate::quant::transcriptome::filter_and_project;

    if transcripts.is_empty() || tr_idx.n_transcripts() == 0 {
        return Ok(Vec::new());
    }

    let mode = params.quant_transcriptome_sam_output;
    let lread = read_seq.len() as u32;
    // STAR passes the RC read to soft-clip extension on reverse-strand
    // alignments (`Read1[2]`); we mirror that here.
    let rc = rc_encode(read_seq);

    let mut projected_all: Vec<crate::align::transcript::Transcript> = Vec::new();
    for aln in transcripts {
        let bases: &[u8] = if aln.is_reverse { &rc } else { read_seq };
        projected_all.extend(filter_and_project(
            aln, bases, genome, tr_idx, lread, mode, params,
        ));
    }

    if projected_all.is_empty() {
        return Ok(Vec::new());
    }

    let (primary_hit, mapq) =
        pick_primary_and_mapq(projected_all.len(), n_for_mapq, read_name, params);

    SamWriter::build_transcriptome_records(
        read_name,
        read_seq,
        read_qual,
        &projected_all,
        mapq,
        params,
        primary_hit,
    )
}

/// Paired-end version of `build_transcriptome_records_se`.
///
/// For each `PairedAlignment`, project mate1 and mate2 onto all transcripts
/// and keep only transcripts where both mates project successfully.  Emit one
/// SAM record per projected pair per mate (2 records per projected hit, in
/// mate1-then-mate2 order).
#[allow(clippy::too_many_arguments)]
fn build_transcriptome_records_pe<'a, I>(
    both_mapped: I,
    read_name: &str,
    m1_seq: &[u8],
    m1_qual: &[u8],
    m2_seq: &[u8],
    m2_qual: &[u8],
    genome: &crate::genome::Genome,
    tr_idx: &crate::quant::transcriptome::TranscriptomeIndex,
    params: &Parameters,
    n_for_mapq: usize,
) -> Result<Vec<noodles::sam::alignment::record_buf::RecordBuf>, error::Error>
where
    I: IntoIterator<Item = &'a crate::align::read_align::PairedAlignment>,
{
    use crate::io::sam::SamWriter;
    use crate::quant::transcriptome::filter_and_project;
    use std::collections::HashMap;

    if tr_idx.n_transcripts() == 0 {
        return Ok(Vec::new());
    }

    let mode = params.quant_transcriptome_sam_output;
    let lread1 = m1_seq.len() as u32;
    let lread2 = m2_seq.len() as u32;
    let m1_rc = rc_encode(m1_seq);
    let m2_rc = rc_encode(m2_seq);

    // For each both-mapped pair, project each mate onto transcripts and pair
    // up projections that land on the same transcript.
    let mut all_projected: Vec<(
        crate::align::transcript::Transcript,
        crate::align::transcript::Transcript,
    )> = Vec::new();
    for pair in both_mapped {
        let m1 = &pair.mate1_transcript;
        let m2 = &pair.mate2_transcript;
        let m1_bases: &[u8] = if m1.is_reverse { &m1_rc } else { m1_seq };
        let m2_bases: &[u8] = if m2.is_reverse { &m2_rc } else { m2_seq };
        let proj_m1 = filter_and_project(m1, m1_bases, genome, tr_idx, lread1, mode, params);
        let proj_m2 = filter_and_project(m2, m2_bases, genome, tr_idx, lread2, mode, params);

        let mut by_tr1: HashMap<usize, Vec<&crate::align::transcript::Transcript>> = HashMap::new();
        for p in &proj_m1 {
            by_tr1.entry(p.chr_idx).or_default().push(p);
        }
        for p2 in &proj_m2 {
            if let Some(p1s) = by_tr1.get(&p2.chr_idx) {
                for p1 in p1s {
                    all_projected.push(((*p1).clone(), p2.clone()));
                }
            }
        }
    }

    if all_projected.is_empty() {
        return Ok(Vec::new());
    }

    let n_alignments = all_projected.len();
    let (primary_hit, mapq) = pick_primary_and_mapq(n_alignments, n_for_mapq, read_name, params);

    // Build one record per mate per projected pair in a single call each,
    // then stamp paired flags and interleave as mate1, mate2, mate1, mate2…
    let (p1s, p2s): (Vec<_>, Vec<_>) = all_projected.into_iter().unzip();
    let mut rec1s = SamWriter::build_transcriptome_records(
        read_name,
        m1_seq,
        m1_qual,
        &p1s,
        mapq,
        params,
        primary_hit,
    )?;
    let mut rec2s = SamWriter::build_transcriptome_records(
        read_name,
        m2_seq,
        m2_qual,
        &p2s,
        mapq,
        params,
        primary_hit,
    )?;

    use noodles::sam::alignment::record::Flags;
    for r in rec1s.iter_mut() {
        *r.flags_mut() |= Flags::SEGMENTED | Flags::FIRST_SEGMENT;
    }
    for r in rec2s.iter_mut() {
        *r.flags_mut() |= Flags::SEGMENTED | Flags::LAST_SEGMENT;
    }

    let mut out: Vec<noodles::sam::alignment::record_buf::RecordBuf> =
        Vec::with_capacity(n_alignments * 2);
    for (r1, r2) in rec1s.into_iter().zip(rec2s) {
        out.push(r1);
        out.push(r2);
    }
    Ok(out)
}

/// Extract SjKey junction identifiers from a transcript's CIGAR.
/// Used to check if a read's junctions survive outSJfilter* for BySJout mode.
fn extract_junction_keys(
    transcript: &crate::align::transcript::Transcript,
    index: &crate::index::GenomeIndex,
) -> Vec<crate::junction::SjKey> {
    use crate::align::score::AlignmentScorer;
    use crate::align::transcript::CigarOp;

    let scorer = AlignmentScorer::from_params_minimal();
    let mut keys = Vec::new();
    let mut genome_pos = transcript.genome_start;

    for op in &transcript.cigar {
        match op {
            CigarOp::RefSkip(len) => {
                let intron_len = *len;
                let intron_start = genome_pos + 1;
                let intron_end = genome_pos + intron_len as u64;

                let motif = scorer.detect_splice_motif(genome_pos, intron_len, &index.genome);
                let strand = match motif.implied_strand() {
                    Some('+') => 1u8,
                    Some('-') => 2u8,
                    _ => 0u8,
                };
                let encoded_motif = crate::junction::encode_motif(motif);

                keys.push(crate::junction::SjKey {
                    chr_idx: transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                    motif: encoded_motif,
                });

                genome_pos += intron_len as u64;
            }
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::Del(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
        }
    }

    keys
}

/// Align single-end reads
#[allow(clippy::too_many_arguments)]
fn align_reads_single_end<W: AlignmentWriter + ?Sized>(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut W,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
    quant_ctx: Option<&std::sync::Arc<crate::quant::QuantContext>>,
    tr_idx: Option<&std::sync::Arc<crate::quant::transcriptome::TranscriptomeIndex>>,
    mut tr_writer: Option<&mut crate::io::bam::BamWriter>,
    mut unmapped_writer: Option<&mut crate::io::fastq::UnmappedFastqWriter>,
) -> anyhow::Result<()> {
    use crate::align::read_align::align_read;
    use crate::io::fastq::{FastqReader, clip_read};
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::params::OutFilterType;
    use rayon::prelude::*;
    use std::sync::Arc;

    let quant = quant_ctx.map(Arc::clone);
    let tr = tr_idx.map(Arc::clone);

    let read_file = &params.read_files_in[0];
    info!("Reading single-end from {}", read_file.display());

    let mut reader = FastqReader::open(read_file, params.read_files_command.as_deref())?;

    // Create chimeric output writer if enabled
    let mut chimeric_writer = if params.chim_segment_min > 0 && params.chim_out_junctions() {
        use crate::chimeric::ChimericJunctionWriter;
        let prefix = params.out_file_name_prefix.to_str().unwrap_or(".");
        info!(
            "Chimeric detection enabled (chimSegmentMin={})",
            params.chim_segment_min
        );
        Some(ChimericJunctionWriter::new(prefix)?)
    } else {
        None
    };

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
    let mut read_count = 0u64;
    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };

    let batch_size = 10000;
    let clip5p = params.clip5p_nbases as usize;
    let clip3p = params.clip3p_nbases as usize;
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    let output_unmapped = params.out_sam_unmapped != params::OutSamUnmapped::None;
    let write_unmapped_fastq = params.out_reads_unmapped == params::OutReadsUnmapped::Fastx;
    let by_sjout = params.out_filter_type == OutFilterType::BySJout;
    let rg_id_owned = params.primary_rg_id()?;

    // BySJout disk buffer: SAM records written to a temp file; only compact metadata kept in RAM.
    // For 100M reads this avoids ~60 GB of Vec<RecordBuf> in memory.
    let bysj_temp = if by_sjout {
        info!("outFilterType=BySJout: disk-buffering reads for post-alignment junction filtering");
        let tf = tempfile::NamedTempFile::new()
            .map_err(|e| anyhow::anyhow!("BySJout: failed to create temp file: {}", e))?;
        Some(tf)
    } else {
        None
    };
    let (bysj_sam_header, mut bysj_temp_writer) = if let Some(ref tf) = bysj_temp {
        let write_file = tf
            .reopen()
            .map_err(|e| anyhow::anyhow!("BySJout: temp file reopen error: {}", e))?;
        let (hdr, w) = crate::io::sam::create_bysj_writer(write_file, &index.genome, params)?;
        (Some(hdr), Some(w))
    } else {
        (None, None)
    };
    let mut bysj_meta: Vec<BySJReadMeta> = Vec::new();

    info!("Aligning reads...");
    loop {
        // Sequential FASTQ reading (unavoidable bottleneck)
        let batch = reader.read_batch(batch_size)?;
        if batch.is_empty() {
            break;
        }

        // Check max reads limit
        let reads_to_process = if read_count + batch.len() as u64 > max_reads {
            (max_reads - read_count) as usize
        } else {
            batch.len()
        };

        let batch_to_process = &batch[..reads_to_process];

        // Parallel alignment processing
        let batch_results: Vec<Result<AlignmentBatchResults, error::Error>> = batch_to_process
            .par_iter()
            .map(|read| {
                #[allow(clippy::needless_borrow)]
                let index = Arc::clone(&index);
                #[allow(clippy::needless_borrow)]
                let stats = Arc::clone(&stats);
                #[allow(clippy::needless_borrow)]
                let sj_stats = Arc::clone(&sj_stats);
                let quant = quant.as_ref().map(Arc::clone);

                // Apply clipping
                let (clipped_seq, clipped_qual) =
                    clip_read(&read.sequence, &read.quality, clip5p, clip3p);

                let mut buffer = BufferedSamRecords::new();
                let mut chimeric_alns = Vec::new();
                let tr_local = tr.as_ref().map(Arc::clone);

                // Record read bases for Log.final.out
                stats.record_read_bases(clipped_seq.len() as u64);

                // Skip if read is too short after clipping
                if clipped_seq.is_empty() {
                    stats.record_alignment(0, max_multimaps);
                    stats.record_unmapped_reason(crate::stats::UnmappedReason::Other);
                    if let Some(ref q) = quant {
                        q.counts.count_se_read(&[], 0, &q.gene_ann);
                    }
                    if output_unmapped {
                        let record = SamWriter::build_unmapped_record(
                            &read.name,
                            &clipped_seq,
                            &clipped_qual,
                            rg_id_owned.as_deref(),
                        )?;
                        buffer.push(record);
                    }
                    let unmapped_m1 = if write_unmapped_fastq {
                        vec![(read.name.clone(), clipped_seq.clone(), clipped_qual.clone())]
                    } else {
                        Vec::new()
                    };
                    return Ok(AlignmentBatchResults {
                        sam_records: buffer,
                        chimeric_alns,
                        primary_junction_keys: Vec::new(),
                        transcriptome_records: Vec::new(),
                        unmapped_mate1: unmapped_m1,
                        unmapped_mate2: Vec::new(),
                    });
                }

                // Align read (CPU-intensive, pure function)
                let (transcripts, chimeric_results, n_for_mapq, unmapped_reason) =
                    align_read(&clipped_seq, &read.name, &index, params)?;

                // Collect chimeric alignments if enabled
                if params.chim_segment_min > 0 {
                    chimeric_alns.extend(chimeric_results);
                    if !chimeric_alns.is_empty() {
                        stats.record_chimeric();
                    }
                }

                // Record stats (atomic, lock-free)
                // For too-many-loci, n_for_mapq carries the true loci count
                // while transcripts is empty
                let n_for_stats = if transcripts.is_empty() && n_for_mapq > 0 {
                    n_for_mapq // too-many-loci: use true count for stats
                } else {
                    transcripts.len()
                };
                stats.record_alignment(n_for_stats, max_multimaps);
                if transcripts.is_empty() && unmapped_reason.is_some() {
                    stats.record_unmapped_reason(
                        unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
                    );
                } else if transcripts.len() == 1 {
                    stats.record_transcript_stats(&transcripts[0]);
                }

                // Gene-level quantification (lock-free atomic counts)
                if let Some(ref q) = quant {
                    q.counts
                        .count_se_read(&transcripts, n_for_mapq, &q.gene_ann);
                }

                // Record junction statistics
                let is_unique = transcripts.len() == 1;
                for transcript in &transcripts {
                    record_transcript_junctions(transcript, &index, &sj_stats, is_unique);
                }

                // Extract junction keys from primary alignment for BySJout filtering
                let primary_junction_keys =
                    if by_sjout && !transcripts.is_empty() && transcripts[0].n_junction > 0 {
                        extract_junction_keys(&transcripts[0], &index)
                    } else {
                        Vec::new()
                    };

                // Build SAM records (no I/O, just construction)
                let is_unmapped_se = transcripts.is_empty();
                if is_unmapped_se {
                    // Unmapped
                    if output_unmapped {
                        let record = SamWriter::build_unmapped_record(
                            &read.name,
                            &clipped_seq,
                            &clipped_qual,
                            rg_id_owned.as_deref(),
                        )?;
                        buffer.push(record);
                    }
                } else if transcripts.len() <= max_multimaps {
                    // Mapped (within multimap limit)
                    let records = SamWriter::build_alignment_records(
                        &read.name,
                        &clipped_seq,
                        &clipped_qual,
                        &transcripts,
                        &index.genome,
                        params,
                        n_for_mapq,
                    )?;
                    for record in records {
                        buffer.push(record);
                    }
                }
                // else: too many loci, skip output

                // Transcriptome SAM projection for --quantMode TranscriptomeSAM.
                let transcriptome_records: Vec<noodles::sam::alignment::record_buf::RecordBuf> =
                    if let Some(ref tidx) = tr_local {
                        build_transcriptome_records_se(
                            &transcripts,
                            &read.name,
                            &clipped_seq,
                            &clipped_qual,
                            &index.genome,
                            tidx,
                            params,
                            n_for_mapq,
                        )?
                    } else {
                        Vec::new()
                    };

                let unmapped_m1 = if write_unmapped_fastq && is_unmapped_se {
                    vec![(read.name.clone(), clipped_seq.clone(), clipped_qual.clone())]
                } else {
                    Vec::new()
                };

                Ok(AlignmentBatchResults {
                    sam_records: buffer,
                    chimeric_alns,
                    primary_junction_keys,
                    transcriptome_records,
                    unmapped_mate1: unmapped_m1,
                    unmapped_mate2: Vec::new(),
                })
            })
            .collect();

        if by_sjout {
            for result in batch_results {
                let batch = result?;
                // Write SAM records to temp file (disk, not RAM)
                let n_sam_records = batch.sam_records.records.len() as u32;
                if let (Some(tw), Some(hdr)) = (&mut bysj_temp_writer, &bysj_sam_header) {
                    crate::io::sam::bysj_write_records(tw, hdr, &batch.sam_records.records)?;
                }
                // Write unmapped reads immediately — they always pass BySJout (no junctions)
                if let Some(ref mut uw) = unmapped_writer {
                    for (name, seq, qual) in &batch.unmapped_mate1 {
                        uw.write_record(name, seq, qual)?;
                    }
                }
                // Store compact metadata in memory (chimeric kept here — rare, ~0.1%)
                bysj_meta.push(BySJReadMeta {
                    n_sam_records,
                    junction_keys: batch.primary_junction_keys,
                    chimeric_alns: batch.chimeric_alns,
                    transcriptome_records: batch.transcriptome_records,
                });
            }
        } else {
            // Normal mode: sequential writing (merge buffers in chunk order)
            for result in batch_results {
                let batch = result?;

                // Write SAM/BAM records
                writer.write_batch(&batch.sam_records.records)?;

                // Write transcriptome-space records (if enabled)
                if let Some(ref mut tw) = tr_writer {
                    tw.write_batch(&batch.transcriptome_records)?;
                }

                // Write chimeric alignments
                if let Some(ref mut chim_writer) = chimeric_writer {
                    for chim_aln in &batch.chimeric_alns {
                        chim_writer.write_alignment(
                            chim_aln,
                            &index.genome.chr_name,
                            &chim_aln.read_name,
                        )?;
                    }
                }
                if params.chim_out_within_bam() {
                    use crate::chimeric::build_within_bam_records;
                    for chim_aln in &batch.chimeric_alns {
                        let supp = build_within_bam_records(chim_aln, &index.genome, 255)?;
                        writer.write_batch(&supp)?;
                    }
                }

                // Write unmapped FASTQ records
                if let Some(ref mut uw) = unmapped_writer {
                    for (name, seq, qual) in &batch.unmapped_mate1 {
                        uw.write_record(name, seq, qual)?;
                    }
                }
            }
        }

        read_count += reads_to_process as u64;

        // Progress logging
        if read_count % 100000 < batch_size as u64 {
            info!("Processed {} reads...", read_count);
        }

        if read_count >= max_reads {
            break;
        }
    }

    // BySJout post-alignment filtering (disk-buffered reads)
    if by_sjout {
        let surviving_junctions = sj_stats.compute_surviving_junctions(params);
        info!(
            "BySJout filtering: {} surviving junctions from {} total",
            surviving_junctions.len(),
            sj_stats.len()
        );

        // Flush and close the temp writer before re-opening for reading
        drop(bysj_temp_writer);

        let mut filtered_count = 0u64;
        if let (Some(tf), Some(hdr)) = (&bysj_temp, &bysj_sam_header) {
            let read_file = tf
                .reopen()
                .map_err(|e| anyhow::anyhow!("BySJout: temp file reopen for reading: {}", e))?;
            let mut reader = noodles::sam::io::Reader::new(std::io::BufReader::new(read_file));
            reader.read_header()?;

            for meta in &bysj_meta {
                let all_survive = meta.junction_keys.is_empty()
                    || meta
                        .junction_keys
                        .iter()
                        .all(|key| surviving_junctions.contains(key));

                if all_survive {
                    let records = crate::io::sam::bysj_read_n_records(
                        &mut reader,
                        hdr,
                        meta.n_sam_records,
                        true,
                    )?;
                    writer.write_batch(&records)?;
                    if let Some(ref mut tw) = tr_writer {
                        tw.write_batch(&meta.transcriptome_records)?;
                    }
                    if let Some(ref mut chim_writer) = chimeric_writer {
                        for chim_aln in &meta.chimeric_alns {
                            chim_writer.write_alignment(
                                chim_aln,
                                &index.genome.chr_name,
                                &chim_aln.read_name,
                            )?;
                        }
                    }
                    if params.chim_out_within_bam() {
                        use crate::chimeric::build_within_bam_records;
                        for chim_aln in &meta.chimeric_alns {
                            let supp = build_within_bam_records(chim_aln, &index.genome, 255)?;
                            writer.write_batch(&supp)?;
                        }
                    }
                } else {
                    // Skip these records in the temp file (advance reader)
                    crate::io::sam::bysj_read_n_records(
                        &mut reader,
                        hdr,
                        meta.n_sam_records,
                        false,
                    )?;
                    filtered_count += 1;
                    stats.undo_mapped_record_bysj();
                }
            }
        }

        info!(
            "BySJout: filtered {} reads with non-surviving junctions",
            filtered_count
        );
    }

    // Flush chimeric output if enabled
    if let Some(ref mut chim_writer) = chimeric_writer {
        chim_writer.flush()?;
        info!("Chimeric junction output complete");
    }

    // Flush unmapped FASTQ writer
    if let Some(ref mut uw) = unmapped_writer {
        uw.flush()?;
    }

    Ok(())
}

/// Align paired-end reads
#[allow(clippy::too_many_arguments)]
fn align_reads_paired_end<W: AlignmentWriter + ?Sized>(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut W,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
    quant_ctx: Option<&std::sync::Arc<crate::quant::QuantContext>>,
    tr_idx: Option<&std::sync::Arc<crate::quant::transcriptome::TranscriptomeIndex>>,
    mut tr_writer: Option<&mut crate::io::bam::BamWriter>,
    mut unmapped_writer1: Option<&mut crate::io::fastq::UnmappedFastqWriter>,
    mut unmapped_writer2: Option<&mut crate::io::fastq::UnmappedFastqWriter>,
) -> anyhow::Result<()> {
    use crate::align::read_align::{PairedAlignment, PairedAlignmentResult, align_paired_read};
    use crate::io::fastq::{PairedFastqReader, clip_read};
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::params::OutFilterType;
    use rayon::prelude::*;
    use std::sync::Arc;

    let quant = quant_ctx.map(Arc::clone);
    let tr = tr_idx.map(Arc::clone);

    info!(
        "Reading paired-end from {} and {}",
        params.read_files_in[0].display(),
        params.read_files_in[1].display()
    );

    let mut reader = PairedFastqReader::open(
        &params.read_files_in[0],
        &params.read_files_in[1],
        params.read_files_command.as_deref(),
    )?;

    // Create chimeric output writer if enabled
    let mut chimeric_writer = if params.chim_segment_min > 0 && params.chim_out_junctions() {
        use crate::chimeric::ChimericJunctionWriter;
        let prefix = params.out_file_name_prefix.to_str().unwrap_or(".");
        info!(
            "Chimeric detection enabled (chimSegmentMin={})",
            params.chim_segment_min
        );
        Some(ChimericJunctionWriter::new(prefix)?)
    } else {
        None
    };

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
    let mut read_count = 0u64;
    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };

    let batch_size = 10000;
    let clip5p = params.clip5p_nbases as usize;
    let clip3p = params.clip3p_nbases as usize;
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    let output_unmapped = params.out_sam_unmapped != params::OutSamUnmapped::None;
    let write_unmapped_fastq = params.out_reads_unmapped == params::OutReadsUnmapped::Fastx;
    let by_sjout = params.out_filter_type == OutFilterType::BySJout;

    // BySJout disk buffer: SAM records to temp file, compact metadata in RAM.
    let bysj_temp = if by_sjout {
        info!("outFilterType=BySJout: disk-buffering pairs for post-alignment junction filtering");
        let tf = tempfile::NamedTempFile::new()
            .map_err(|e| anyhow::anyhow!("BySJout: failed to create temp file: {}", e))?;
        Some(tf)
    } else {
        None
    };
    let (bysj_sam_header, mut bysj_temp_writer) = if let Some(ref tf) = bysj_temp {
        let write_file = tf
            .reopen()
            .map_err(|e| anyhow::anyhow!("BySJout: temp file reopen error: {}", e))?;
        let (hdr, w) = crate::io::sam::create_bysj_writer(write_file, &index.genome, params)?;
        (Some(hdr), Some(w))
    } else {
        (None, None)
    };
    let mut bysj_meta: Vec<BySJReadMeta> = Vec::new();

    info!("Aligning paired-end reads...");
    loop {
        // Sequential FASTQ reading
        let batch = reader.read_paired_batch(batch_size)?;
        if batch.is_empty() {
            break;
        }

        // Check max reads limit (pairs, not individual reads)
        let pairs_to_process = if read_count + batch.len() as u64 > max_reads {
            (max_reads - read_count) as usize
        } else {
            batch.len()
        };

        let batch_to_process = &batch[..pairs_to_process];

        // Parallel alignment processing
        let batch_results: Vec<Result<AlignmentBatchResults, error::Error>> = batch_to_process
            .par_iter()
            .map(|paired_read| {
                #[allow(clippy::needless_borrow)]
                let index = Arc::clone(&index);
                #[allow(clippy::needless_borrow)]
                let stats = Arc::clone(&stats);
                #[allow(clippy::needless_borrow)]
                let sj_stats = Arc::clone(&sj_stats);
                let quant = quant.as_ref().map(Arc::clone);

                // Apply clipping to both mates
                let (m1_seq, m1_qual) = clip_read(
                    &paired_read.mate1.sequence,
                    &paired_read.mate1.quality,
                    clip5p,
                    clip3p,
                );
                let (m2_seq, m2_qual) = clip_read(
                    &paired_read.mate2.sequence,
                    &paired_read.mate2.quality,
                    clip5p,
                    clip3p,
                );

                let mut buffer = BufferedSamRecords::new();
                let tr_local = tr.as_ref().map(Arc::clone);

                // Record read bases for Log.final.out (both mates)
                stats.record_read_bases(m1_seq.len() as u64 + m2_seq.len() as u64);

                // Skip if either mate is too short after clipping
                if m1_seq.is_empty() || m2_seq.is_empty() {
                    stats.record_alignment(0, max_multimaps);
                    stats.record_unmapped_reason(crate::stats::UnmappedReason::Other);
                    if let Some(ref q) = quant {
                        q.counts.count_pe_read(&[], true, false, &q.gene_ann);
                    }
                    if output_unmapped {
                        let records = SamWriter::build_paired_unmapped_records(
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                            params,
                        )?;
                        for record in records {
                            buffer.push(record);
                        }
                    }
                    let (um1, um2) = if write_unmapped_fastq {
                        (
                            vec![(
                                paired_read.mate1.name.clone(),
                                m1_seq.clone(),
                                m1_qual.clone(),
                            )],
                            vec![(
                                paired_read.mate2.name.clone(),
                                m2_seq.clone(),
                                m2_qual.clone(),
                            )],
                        )
                    } else {
                        (Vec::new(), Vec::new())
                    };
                    return Ok(AlignmentBatchResults {
                        sam_records: buffer,
                        chimeric_alns: Vec::new(),
                        primary_junction_keys: Vec::new(),
                        transcriptome_records: Vec::new(),
                        unmapped_mate1: um1,
                        unmapped_mate2: um2,
                    });
                }

                // Align paired read (CPU-intensive)
                let (results, pe_chimeric, n_for_mapq, unmapped_reason) =
                    align_paired_read(&m1_seq, &m2_seq, &paired_read.name, &index, params)?;

                // Classify the result for stats and SAM output
                let has_half_mapped = results
                    .iter()
                    .any(|r| matches!(r, PairedAlignmentResult::HalfMapped { .. }));
                let both_mapped: Vec<_> = results
                    .iter()
                    .filter_map(|r| {
                        if let PairedAlignmentResult::BothMapped(pa) = r {
                            Some(pa)
                        } else {
                            None
                        }
                    })
                    .collect();

                if results.is_empty() {
                    // Both mates unmapped
                    stats.record_alignment(0, max_multimaps);
                    stats.record_unmapped_reason(
                        unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
                    );
                } else if has_half_mapped {
                    // Half-mapped: count as mapped for the mapped mate
                    stats.record_alignment(1, max_multimaps);
                    stats.record_half_mapped();
                    // Record transcript stats from the mapped mate only
                    if let Some(PairedAlignmentResult::HalfMapped {
                        mapped_transcript, ..
                    }) = results.first()
                    {
                        stats.record_transcript_stats(mapped_transcript);
                    }
                } else {
                    // Both-mapped pairs
                    let n = both_mapped.len();
                    stats.record_alignment(n, max_multimaps);
                    if n == 1 {
                        stats.record_transcript_stats(&both_mapped[0].mate1_transcript);
                        stats.record_transcript_stats(&both_mapped[0].mate2_transcript);
                    }
                }

                // Gene-level quantification (lock-free atomic counts)
                if let Some(ref q) = quant {
                    // Dereference Box<PairedAlignment> to get &PairedAlignment slice.
                    let bm_deref: Vec<&crate::align::read_align::PairedAlignment> =
                        both_mapped.iter().map(|b| b.as_ref()).collect();
                    q.counts.count_pe_read(
                        &bm_deref,
                        results.is_empty(),
                        has_half_mapped,
                        &q.gene_ann,
                    );
                }

                // Record junction statistics
                let is_unique = both_mapped.len() == 1 || (has_half_mapped && results.len() == 1);
                for result in &results {
                    match result {
                        PairedAlignmentResult::BothMapped(pair) => {
                            record_transcript_junctions(
                                &pair.mate1_transcript,
                                &index,
                                &sj_stats,
                                is_unique,
                            );
                            record_transcript_junctions(
                                &pair.mate2_transcript,
                                &index,
                                &sj_stats,
                                is_unique,
                            );
                        }
                        PairedAlignmentResult::HalfMapped {
                            mapped_transcript, ..
                        } => {
                            record_transcript_junctions(
                                mapped_transcript,
                                &index,
                                &sj_stats,
                                is_unique,
                            );
                        }
                    }
                }

                // Extract junction keys from primary alignment for BySJout
                let primary_junction_keys = if by_sjout && !results.is_empty() {
                    let mut keys = Vec::new();
                    match &results[0] {
                        PairedAlignmentResult::BothMapped(pair) => {
                            if pair.mate1_transcript.n_junction > 0 {
                                keys.extend(extract_junction_keys(&pair.mate1_transcript, &index));
                            }
                            if pair.mate2_transcript.n_junction > 0 {
                                keys.extend(extract_junction_keys(&pair.mate2_transcript, &index));
                            }
                        }
                        PairedAlignmentResult::HalfMapped {
                            mapped_transcript, ..
                        } => {
                            if mapped_transcript.n_junction > 0 {
                                keys.extend(extract_junction_keys(mapped_transcript, &index));
                            }
                        }
                    }
                    keys
                } else {
                    Vec::new()
                };

                // Build SAM records
                if results.is_empty() {
                    // Unmapped pair
                    if output_unmapped {
                        let records = SamWriter::build_paired_unmapped_records(
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                            params,
                        )?;
                        for record in records {
                            buffer.push(record);
                        }
                    }
                } else if has_half_mapped {
                    // Half-mapped pair
                    if let Some(PairedAlignmentResult::HalfMapped {
                        mapped_transcript,
                        mate1_is_mapped,
                    }) = results.first()
                    {
                        let records = SamWriter::build_half_mapped_records(
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                            mapped_transcript,
                            *mate1_is_mapped,
                            &index.genome,
                            params,
                            n_for_mapq,
                        )?;
                        for record in records {
                            buffer.push(record);
                        }
                    }
                } else if both_mapped.len() <= max_multimaps {
                    // Both-mapped pairs (within multimap limit)
                    // Extract PairedAlignments for the existing build_paired_records
                    let paired_alns: Vec<PairedAlignment> = both_mapped
                        .iter()
                        .map(|pa| PairedAlignment::clone(pa))
                        .collect();
                    let records = SamWriter::build_paired_records(
                        &paired_read.name,
                        &m1_seq,
                        &m1_qual,
                        &m2_seq,
                        &m2_qual,
                        &paired_alns,
                        &index.genome,
                        params,
                        n_for_mapq,
                    )?;
                    for record in records {
                        buffer.push(record);
                    }
                }
                // else: too many loci, skip output

                // Transcriptome SAM projection (both-mapped pairs only)
                let transcriptome_records: Vec<noodles::sam::alignment::record_buf::RecordBuf> =
                    if let Some(ref tidx) = tr_local {
                        build_transcriptome_records_pe(
                            both_mapped.iter().map(|b| b.as_ref()),
                            &paired_read.name,
                            &m1_seq,
                            &m1_qual,
                            &m2_seq,
                            &m2_qual,
                            &index.genome,
                            tidx,
                            params,
                            n_for_mapq,
                        )?
                    } else {
                        Vec::new()
                    };

                // Collect unmapped mates for --outReadsUnmapped Fastx.
                // Write both mates if: pair is fully unmapped OR half-mapped.
                // STAR writes both mates of half-mapped pairs to the unmapped files.
                let (unmapped_mate1, unmapped_mate2) = if write_unmapped_fastq {
                    let pair_unmapped = results.is_empty() || has_half_mapped;
                    if pair_unmapped {
                        (
                            vec![(
                                paired_read.mate1.name.clone(),
                                m1_seq.clone(),
                                m1_qual.clone(),
                            )],
                            vec![(
                                paired_read.mate2.name.clone(),
                                m2_seq.clone(),
                                m2_qual.clone(),
                            )],
                        )
                    } else {
                        (Vec::new(), Vec::new())
                    }
                } else {
                    (Vec::new(), Vec::new())
                };

                Ok(AlignmentBatchResults {
                    sam_records: buffer,
                    chimeric_alns: pe_chimeric,
                    primary_junction_keys,
                    transcriptome_records,
                    unmapped_mate1,
                    unmapped_mate2,
                })
            })
            .collect();

        if by_sjout {
            for result in batch_results {
                let batch = result?;
                let n_sam_records = batch.sam_records.records.len() as u32;
                if let (Some(tw), Some(hdr)) = (&mut bysj_temp_writer, &bysj_sam_header) {
                    crate::io::sam::bysj_write_records(tw, hdr, &batch.sam_records.records)?;
                }
                // Write unmapped reads immediately — they always pass BySJout
                if let Some(ref mut uw1) = unmapped_writer1 {
                    for (name, seq, qual) in &batch.unmapped_mate1 {
                        uw1.write_record(name, seq, qual)?;
                    }
                }
                if let Some(ref mut uw2) = unmapped_writer2 {
                    for (name, seq, qual) in &batch.unmapped_mate2 {
                        uw2.write_record(name, seq, qual)?;
                    }
                }
                bysj_meta.push(BySJReadMeta {
                    n_sam_records,
                    junction_keys: batch.primary_junction_keys,
                    chimeric_alns: batch.chimeric_alns,
                    transcriptome_records: batch.transcriptome_records,
                });
            }
        } else {
            // Normal mode: sequential SAM writing
            for result in batch_results {
                let batch = result?;
                writer.write_batch(&batch.sam_records.records)?;
                if let Some(ref mut tw) = tr_writer {
                    tw.write_batch(&batch.transcriptome_records)?;
                }
                if params.chim_out_within_bam() {
                    use crate::chimeric::build_within_bam_records;
                    for chim_aln in &batch.chimeric_alns {
                        let supp = build_within_bam_records(chim_aln, &index.genome, 255)?;
                        writer.write_batch(&supp)?;
                    }
                }
                if let Some(ref mut uw1) = unmapped_writer1 {
                    for (name, seq, qual) in &batch.unmapped_mate1 {
                        uw1.write_record(name, seq, qual)?;
                    }
                }
                if let Some(ref mut uw2) = unmapped_writer2 {
                    for (name, seq, qual) in &batch.unmapped_mate2 {
                        uw2.write_record(name, seq, qual)?;
                    }
                }
            }
        }

        read_count += pairs_to_process as u64;

        // Progress logging
        if read_count % 100000 < batch_size as u64 {
            info!("Processed {} pairs...", read_count);
        }

        if read_count >= max_reads {
            break;
        }
    }

    // BySJout post-alignment filtering (disk-buffered pairs)
    if by_sjout {
        let surviving_junctions = sj_stats.compute_surviving_junctions(params);
        info!(
            "BySJout filtering: {} surviving junctions from {} total",
            surviving_junctions.len(),
            sj_stats.len()
        );

        // Flush and close the temp writer before re-opening for reading
        drop(bysj_temp_writer);

        let mut filtered_count = 0u64;
        if let (Some(tf), Some(hdr)) = (&bysj_temp, &bysj_sam_header) {
            let read_file = tf
                .reopen()
                .map_err(|e| anyhow::anyhow!("BySJout: temp file reopen for reading: {}", e))?;
            let mut reader = noodles::sam::io::Reader::new(std::io::BufReader::new(read_file));
            reader.read_header()?;

            for meta in &bysj_meta {
                let all_survive = meta.junction_keys.is_empty()
                    || meta
                        .junction_keys
                        .iter()
                        .all(|key| surviving_junctions.contains(key));

                if all_survive {
                    let records = crate::io::sam::bysj_read_n_records(
                        &mut reader,
                        hdr,
                        meta.n_sam_records,
                        true,
                    )?;
                    writer.write_batch(&records)?;
                    if let Some(ref mut tw) = tr_writer {
                        tw.write_batch(&meta.transcriptome_records)?;
                    }
                    if params.chim_out_within_bam() {
                        use crate::chimeric::build_within_bam_records;
                        for chim_aln in &meta.chimeric_alns {
                            let supp = build_within_bam_records(chim_aln, &index.genome, 255)?;
                            writer.write_batch(&supp)?;
                        }
                    }
                } else {
                    crate::io::sam::bysj_read_n_records(
                        &mut reader,
                        hdr,
                        meta.n_sam_records,
                        false,
                    )?;
                    filtered_count += 1;
                    stats.undo_mapped_record_bysj();
                }
            }
        }

        info!(
            "BySJout: filtered {} pairs with non-surviving junctions",
            filtered_count
        );
    }

    // Flush chimeric output if enabled
    if let Some(ref mut chim_writer) = chimeric_writer {
        chim_writer.flush()?;
    }

    // Flush unmapped FASTQ writers
    if let Some(ref mut uw1) = unmapped_writer1 {
        uw1.flush()?;
    }
    if let Some(ref mut uw2) = unmapped_writer2 {
        uw2.flush()?;
    }

    Ok(())
}

/// Record junctions from a transcript into SJ statistics
fn record_transcript_junctions(
    transcript: &crate::align::transcript::Transcript,
    index: &crate::index::GenomeIndex,
    sj_stats: &crate::junction::SpliceJunctionStats,
    is_unique: bool,
) {
    use crate::align::score::AlignmentScorer;
    use crate::align::transcript::CigarOp;

    // First pass: compute exon segment lengths (query-consuming bases between N operations)
    // An "exon segment" is the query bases on each side of a splice junction.
    let mut exon_lengths: Vec<u32> = Vec::new();
    let mut current_exon_len = 0u32;

    for op in &transcript.cigar {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                current_exon_len += *len;
            }
            CigarOp::Ins(len) => {
                current_exon_len += *len;
            }
            CigarOp::RefSkip(_) => {
                exon_lengths.push(current_exon_len);
                current_exon_len = 0;
            }
            // Soft clips, deletions, hard clips do not contribute to overhang
            // STAR counts only matched/inserted bases (not soft-clipped bases)
            CigarOp::SoftClip(_) | CigarOp::Del(_) | CigarOp::HardClip(_) => {}
        }
    }
    exon_lengths.push(current_exon_len); // Final exon segment

    // Second pass: record junctions with computed overhangs
    let mut genome_pos = transcript.genome_start;
    let mut junction_idx = 0usize;

    let scorer = AlignmentScorer::from_params_minimal();

    for op in &transcript.cigar {
        match op {
            CigarOp::RefSkip(len) => {
                // This is a splice junction
                let intron_len = *len;
                let intron_start = genome_pos + 1; // 1-based, first intronic base
                let intron_end = genome_pos + intron_len as u64; // 1-based, last intronic base

                // Detect splice motif
                let motif = scorer.detect_splice_motif(genome_pos, intron_len, &index.genome);

                // Compute overhang: min(left_exon_length, right_exon_length)
                let left_exon = exon_lengths[junction_idx];
                let right_exon = exon_lengths[junction_idx + 1];
                let overhang = left_exon.min(right_exon);

                // Derive strand from splice motif (STAR convention)
                let strand = match motif.implied_strand() {
                    Some('+') => 1u8,
                    Some('-') => 2u8,
                    _ => 0u8, // non-canonical: unknown strand
                };
                let annotated = index.junction_db.is_annotated(
                    transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                );

                // Record junction
                sj_stats.record_junction(
                    transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                    motif,
                    is_unique,
                    overhang,
                    annotated,
                );

                // Advance genome position past the intron
                genome_pos += intron_len as u64;
                junction_idx += 1;
            }
            CigarOp::Match(len) | CigarOp::Equal(len) | CigarOp::Diff(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::Ins(_) => {}
            CigarOp::Del(len) => {
                genome_pos += *len as u64;
            }
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
        }
    }
}
