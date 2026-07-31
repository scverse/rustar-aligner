#![warn(clippy::pedantic)]
// TODO: enable these warnings eventually
#![allow(
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::doc_markdown,
    clippy::items_after_statements,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::must_use_candidate,
    clippy::similar_names,
    clippy::too_many_lines,
    // trailing comment because of https://github.com/rust-lang/rustfmt/issues/3277
)]
// These should stay disabled
#![allow(
    // we have a bunch of “`if !reverse`”
    clippy::if_not_else,
)]

pub mod error;
pub mod params;

pub mod align;
pub mod bam_dedup;
pub mod chimeric;
pub mod clip;
pub mod cpu;
pub mod genome;
pub mod index;
pub mod io;
pub mod junction;
pub mod liftover;
pub mod mapq;
pub mod quant;
pub mod rng;
pub mod signal;
pub mod solo;
pub mod stats;
pub mod wasp;

use log::info;
use noodles::sam::alignment::record::cigar;

use crate::params::{Parameters, RunMode};

/// Top-level dispatcher. Called from `main()` after CLI parsing.
pub fn run(params: &Parameters) -> anyhow::Result<()> {
    info!("rustar-aligner {}", env!("CARGO_PKG_VERSION"));
    info!("{}", env!("VERSION_BODY"));
    info!("{}", cpu::cpu_detected_line());
    if let Some(hint) = cpu::upgrade_hint() {
        info!("{hint}");
    }
    info!("runMode: {}", params.run_mode);
    info!("runThreadN: {}", params.run_thread_n);

    // Configure the rayon global pool from `--runThreadN` **before**
    // dispatching to either run-mode. Without this, rayon falls back
    // to `num_cpus::get()` (= logical cores), which on big servers
    // (e.g. 256-core machines) spawns hundreds of worker threads
    // independent of `--runThreadN`. With a thread-caching allocator
    // (we use mimalloc) each thread retains some MB of per-thread
    // heap, adding ~64 MB × n_threads of allocator overhead to peak
    // RSS for nothing. `build_global` errors if called twice; we
    // ignore the error so the in-process tests that already
    // initialised the pool still work.
    if usize::from(params.run_thread_n) > 1 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(params.run_thread_n.into())
            .build_global();
    }

    match params.run_mode {
        RunMode::GenomeGenerate => genome_generate(params),
        RunMode::AlignReads => align_reads(params),
        RunMode::InputAlignmentsFromBAM => bam_dedup::run(params),
        RunMode::LiftOver => liftover::run(params),
    }
}

fn genome_generate(params: &Parameters) -> anyhow::Result<()> {
    use index::GenomeIndex;

    info!("genomeDir: {}", params.genome_dir.display());
    if let Some(temp_dir) = &params.temp_dir {
        info!("tempDir: {}", temp_dir.display());
    }
    info!(
        "genomeFastaFiles: {:?}",
        params
            .genome_fasta_files
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
    );

    if params.limit_genome_generate_ram != 31_000_000_000 {
        log::warn!(
            "--limitGenomeGenerateRAM {} accepted but not enforced; rustar manages genome-generation memory independently",
            params.limit_genome_generate_ram
        );
    }

    info!("Building genome index (streaming SA + on-the-fly SAindex)...");
    // Streaming path: opens SA file early, packs each caps-sa emit
    // directly to disk + into the SAindex builder, never holding the
    // ~25 GB SA `PackedArray` in RAM. The in-memory `GenomeIndex::build`
    // + `write` remains for tests / random-access callers.
    GenomeIndex::generate_streaming(params)?;

    // --genomeTransformType Haploid/Diploid: also write the untransformed index to
    // `<genomeDir>/OriginalGenome/` (STAR's own layout), by re-running the same
    // build with the transform disabled and the output directory redirected.
    // `transformGenomeBlocks.tsv` (the reverse-conversion block map) was already
    // written into `genomeDir` above, as part of the transformed Genome's own
    // `write_index_files` call. One single original genome is written even for
    // Diploid (both haplotypes substitute the same shared reference).
    if params.genome_transform_type.eq_ignore_ascii_case("Haploid")
        || params.genome_transform_type.eq_ignore_ascii_case("Diploid")
    {
        let mut orig_params = params.clone();
        orig_params.genome_transform_type = "None".to_string();
        orig_params.genome_transform_vcf = None;
        orig_params.genome_dir = params.genome_dir.join("OriginalGenome");
        info!(
            "genomeTransformType {}: writing untransformed index to {}",
            params.genome_transform_type,
            orig_params.genome_dir.display()
        );
        GenomeIndex::generate_streaming(&orig_params)?;
    }

    info!("Genome generation complete!");
    Ok(())
}

/// Trait for alignment output writers (SAM or BAM).
/// `finish` flushes/sorts/closes the output; default is a no-op for streaming writers.
// `Send` supertrait so a `Box<dyn AlignmentWriter>` (and `&mut dyn AlignmentWriter`)
// can be moved onto the background writer thread in the align pipelines. Every impl
// wraps a `File`/`Stdout`/`Vec`, all of which are `Send`.
trait AlignmentWriter: Send {
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

    // Rayon thread pool was already configured by `run()` from
    // `--runThreadN`; just log the choice here for parity with the
    // previous behaviour.
    if usize::from(params.run_thread_n) > 1 {
        info!("Using {} threads for alignment", params.run_thread_n);
    } else {
        info!("Using single-threaded mode");
    }

    // Validate read files (SmartSeq supplies reads via --readFilesManifest).
    if params.read_files_in.is_empty() && params.solo_type != params::SoloType::SmartSeq {
        anyhow::bail!("No read files specified (--readFilesIn)");
    }

    // 1. Load genome index
    info!("Loading genome index from {}", params.genome_dir.display());
    let index = Arc::new(GenomeIndex::load(&params.genome_dir, params)?);
    let time_genome_loaded = chrono::Local::now();
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

    // SmartSeq has no barcodes/UMIs — a dedicated manifest-driven path.
    if params.solo_type == params::SoloType::SmartSeq {
        let stats = run_smartseq(&index, &params)?;
        let log_path = params.output_path("Log.final.out");
        if let Some(parent) = log_path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        stats.write_log_final(
            &log_path,
            time_start,
            chrono::Local::now(),
            chrono::Local::now(),
        )?;
        info!("Alignment complete!");
        return Ok(());
    }

    // Build the STARsolo context (whitelist + gene model) if a droplet solo run.
    let solo_ctx: Option<std::sync::Arc<crate::solo::SoloContext>> = if params.solo_enabled() {
        info!(
            "STARsolo: soloType={} — building barcode + gene context",
            params.solo_type
        );
        Some(std::sync::Arc::new(crate::solo::SoloContext::build(
            &params,
            &index.genome,
        )?))
    } else {
        None
    };

    let time_map_start = chrono::Local::now();

    // 2. Dispatch based on two-pass mode
    let stats = match params.twopass_mode {
        TwopassMode::None => {
            info!("Running single-pass alignment");
            run_single_pass(
                &index,
                &params,
                quant_ctx.as_ref(),
                tr_idx.as_ref(),
                solo_ctx.as_ref(),
            )?
        }
        TwopassMode::Basic => {
            info!("Running two-pass alignment mode");
            run_two_pass(
                &index,
                &params,
                quant_ctx.as_ref(),
                tr_idx.as_ref(),
                solo_ctx.as_ref(),
            )?
        }
    };

    let time_finish = chrono::Local::now();

    // Write Log.final.out
    let log_path = params.output_path("Log.final.out");
    if let Some(parent) = log_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    stats.write_log_final(&log_path, time_start, time_map_start, time_finish)?;
    info!("Wrote {}", log_path.display());

    // Write Log.out and Log.progress.out
    let log_out_path = params.output_path("Log.out");
    crate::io::log::write_log_out(
        &log_out_path,
        &params,
        &index.genome,
        time_start,
        time_genome_loaded,
        time_finish,
    )?;
    info!("Wrote {}", log_out_path.display());

    let log_progress_path = params.output_path("Log.progress.out");
    crate::io::log::write_log_progress_out(&log_progress_path, &stats, time_start, time_finish)?;
    info!("Wrote {}", log_progress_path.display());

    // Write ReadsPerGene.out.tab if quantMode GeneCounts was requested.
    if let Some(ref ctx) = quant_ctx {
        let quant_path = params.output_path("ReadsPerGene.out.tab");
        ctx.counts.write_output(&quant_path, &ctx.gene_ann)?;
        info!("Wrote {}", quant_path.display());
    }

    info!("Alignment complete!");
    Ok(())
}

/// Log STARsolo barcode/record stats and write the per-cell matrices (raw +
/// filtered), `Summary.csv`, and the SJ feature matrix. Called from the solo
/// branch of `run_single_pass`, where `sj_stats` is live.
fn write_solo_output(
    sctx: &std::sync::Arc<crate::solo::SoloContext>,
    params: &Parameters,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
) -> anyhow::Result<()> {
    use std::sync::atomic::Ordering;
    let s = &sctx.stats;
    info!(
        "STARsolo barcode stats: exact={} 1MM={} multiMM={} noMatch={} N-in-CB={} multReject={} N-in-UMI={} UMIhomopolymer={}",
        s.yes_exact.load(Ordering::Relaxed),
        s.yes_one_mm.load(Ordering::Relaxed),
        s.yes_mult_mm.load(Ordering::Relaxed),
        s.no_match.load(Ordering::Relaxed),
        s.n_in_cb.load(Ordering::Relaxed),
        s.mult_rejected.load(Ordering::Relaxed),
        s.n_in_umi.load(Ordering::Relaxed),
        s.umi_homopolymer.load(Ordering::Relaxed),
    );
    for (feature, recorder) in sctx.features.iter().zip(&sctx.recorders) {
        info!(
            "STARsolo {}: collected {} resolved (CB,UMI,gene) records ({} deferred 1MM_multi)",
            feature.dir_name(),
            recorder.n_records(),
            recorder.n_multi_records(),
        );
    }
    crate::solo::write_gene_matrix(sctx, params, stats, Some(&**sj_stats), &index.genome)?;
    Ok(())
}

/// `--soloType SmartSeq`: align each manifest cell's reads and count reads per
/// gene (no barcodes, no UMIs). Writes `Solo.out/Gene/raw/` (genes × cells) and
/// returns the alignment stats.
fn run_smartseq(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
) -> anyhow::Result<std::sync::Arc<crate::stats::AlignmentStats>> {
    use crate::align::read_align::{PairedAlignmentResult, align_paired_read, align_read};
    use crate::solo::{GeneAssignment, SoloStrand, classify_read};
    use rayon::prelude::*;
    use std::sync::Arc;

    let manifest = params
        .read_files_manifest
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--soloType SmartSeq requires --readFilesManifest"))?;
    let cells = crate::solo::smartseq::parse_manifest(manifest)?;
    info!(
        "STARsolo SmartSeq: {} cells from {}",
        cells.len(),
        manifest.display()
    );

    let gtf = params.sjdb_gtf_file.as_ref().ok_or_else(|| {
        anyhow::anyhow!("--soloType SmartSeq Gene counting requires --sjdbGTFfile")
    })?;
    let exons = crate::junction::gtf::parse_gtf_configured(
        gtf,
        &params.sjdb_gtf_feature_exon,
        &params.sjdb_gtf_chr_prefix,
    )?;
    let gene_ann = crate::quant::GeneAnnotation::from_gtf_exons_configured(
        &exons,
        &index.genome,
        &params.sjdb_gtf_tag_exon_parent_gene,
    );
    info!(
        "STARsolo SmartSeq: {} genes from {}",
        gene_ann.n_genes(),
        gtf.display()
    );
    let strand: SoloStrand = params.solo_strand.parse().unwrap_or_default();
    let max_multimaps = params.out_filter_multimap_nmax as usize;

    let stats = Arc::new(crate::stats::AlignmentStats::new());
    let cell_ids: Vec<String> = cells.iter().map(|c| c.cell_id.clone()).collect();
    let counts = crate::solo::smartseq::SmartSeqCounts::new(cell_ids, gene_ann.gene_ids.len());

    // Assign a (possibly multi-locus) read/fragment to a gene and count it.
    let assign_count = |ci: usize, transcripts: &[crate::align::transcript::Transcript]| {
        if let GeneAssignment::Gene(g) =
            classify_read(transcripts, &gene_ann, strand, true, false, false).gene
        {
            counts.add(ci, g);
        }
    };
    let cmd = params.read_files_command.as_deref();

    for (ci, cell) in cells.iter().enumerate() {
        match &cell.read2 {
            // Single-end: count reads.
            None => {
                let mut reader = crate::io::fastq::FastqReader::open(&cell.read1, cmd)?;
                loop {
                    let batch = reader.read_batch(10_000)?;
                    if batch.is_empty() {
                        break;
                    }
                    batch.par_iter().for_each(|read| {
                        stats.record_read_bases(read.sequence.len() as u64);
                        let Ok((transcripts, _chim, n_for_mapq, reason)) =
                            align_read(&read.sequence, &read.name, index, params)
                        else {
                            return;
                        };
                        let n = if transcripts.is_empty() && n_for_mapq > 0 {
                            n_for_mapq
                        } else {
                            transcripts.len()
                        };
                        stats.record_alignment(n, max_multimaps);
                        if transcripts.is_empty() {
                            stats.record_unmapped_reason(
                                reason.unwrap_or(crate::stats::UnmappedReason::Other),
                            );
                        } else if transcripts.len() == 1 {
                            stats.record_transcript_stats(&transcripts[0]);
                        }
                        assign_count(ci, &transcripts);
                    });
                }
            }
            // Paired-end: align both mates as a fragment, count the fragment once
            // (gene from the union of both mates' overlaps).
            Some(r2) => {
                let mut reader = crate::io::fastq::PairedFastqReader::open(&cell.read1, r2, cmd)?;
                loop {
                    let mut batch = Vec::with_capacity(10_000);
                    while batch.len() < 10_000 {
                        match reader.next_paired()? {
                            Some(p) => batch.push(p),
                            None => break,
                        }
                    }
                    if batch.is_empty() {
                        break;
                    }
                    batch.par_iter().for_each(|pr| {
                        stats.record_read_bases(
                            (pr.mate1.sequence.len() + pr.mate2.sequence.len()) as u64,
                        );
                        let Ok((results, _chim, n_for_mapq, reason)) = align_paired_read(
                            &pr.mate1.sequence,
                            &pr.mate2.sequence,
                            &pr.name,
                            index,
                            params,
                        ) else {
                            return;
                        };
                        let n_pairs = results.len();
                        let mut trs = Vec::with_capacity(n_pairs * 2);
                        for r in results {
                            match r {
                                PairedAlignmentResult::BothMapped(pa) => {
                                    trs.push(pa.mate1_transcript);
                                    trs.push(pa.mate2_transcript);
                                }
                                PairedAlignmentResult::HalfMapped {
                                    mapped_transcript, ..
                                } => trs.push(mapped_transcript),
                            }
                        }
                        let n = if trs.is_empty() && n_for_mapq > 0 {
                            n_for_mapq
                        } else {
                            n_pairs
                        };
                        stats.record_alignment(n, max_multimaps);
                        if trs.is_empty() {
                            stats.record_unmapped_reason(
                                reason.unwrap_or(crate::stats::UnmappedReason::Other),
                            );
                        }
                        assign_count(ci, &trs);
                    });
                }
            }
        }
    }

    let solo_dir = params
        .solo_out_file_names
        .first()
        .cloned()
        .unwrap_or_else(|| "Solo.out/".to_string());
    let raw_dir = params.output_path(&format!("{solo_dir}Gene/raw/"));
    let gzip = matches!(params.solo_out_gzip.as_str(), "yes" | "Yes" | "true");
    let nnz = counts.write_matrix(&raw_dir, &gene_ann.gene_ids, &gene_ann.gene_names, gzip)?;
    info!(
        "STARsolo SmartSeq: wrote Gene/raw matrix ({} genes × {} cells, {} entries)",
        gene_ann.n_genes(),
        cells.len(),
        nnz,
    );
    stats.print_summary();
    Ok(stats)
}

/// Run single-pass alignment (original logic)
fn run_single_pass(
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    params: &Parameters,
    quant_ctx: Option<&std::sync::Arc<crate::quant::QuantContext>>,
    tr_idx: Option<&std::sync::Arc<crate::quant::transcriptome::TranscriptomeIndex>>,
    solo_ctx: Option<&std::sync::Arc<crate::solo::SoloContext>>,
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
        let path = params.output_path("Aligned.toTranscriptome.out.bam");
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

    let is_paired = params.read_files_in.len() == 2 && !params.solo_enabled();
    let mut unmapped_w1: Option<UnmappedFastqWriter> =
        if params.out_reads_unmapped == OutReadsUnmapped::Fastx {
            let path = params.output_path("Unmapped.out.mate1");
            info!("Writing unmapped reads to {}", path.display());
            Some(UnmappedFastqWriter::create(&path)?)
        } else {
            None
        };
    let mut unmapped_w2: Option<UnmappedFastqWriter> =
        if params.out_reads_unmapped == OutReadsUnmapped::Fastx && is_paired {
            let path = params.output_path("Unmapped.out.mate2");
            info!("Writing unmapped mate2 reads to {}", path.display());
            Some(UnmappedFastqWriter::create(&path)?)
        } else {
            None
        };

    // 4. Route to SAM or BAM output based on --outSAMtype / --outStd
    use crate::params::{OutSamSortOrder, OutStd};

    let out_type = &params.out_sam_type;

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
                let output_path = params.output_path("Aligned.out.sam");
                info!("Writing SAM to {}", output_path.display());
                if let Some(parent) = output_path.parent() {
                    std::fs::create_dir_all(parent)?;
                }
                Box::new(SamWriter::create(&output_path, &index.genome, params)?)
            }
            OutSamFormat::Bam => {
                let sorted = out_type.sort_order == Some(OutSamSortOrder::SortedByCoordinate);
                let output_path = if sorted {
                    params.output_path("Aligned.sortedByCoord.out.bam")
                } else {
                    params.output_path("Aligned.out.bam")
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
                info!("--outSAMtype None: skipping alignment output (count/quant only)");
                Box::new(NullWriter)
            }
        },
    };

    // Align reads through the boxed writer.
    //
    // Solo runs supply two `--readFilesIn` files (cDNA read + barcode read) but
    // are single-end *alignment* runs: only the cDNA read (file 0) is aligned.
    // The dedicated solo loop reads the barcode read in lockstep, quantifies
    // per cell, and otherwise emits the cDNA alignments like the SE path.
    if let Some(sctx) = solo_ctx {
        // `--soloBarcodeMate 1` (5' 10x): barcode on mate 1, both mates aligned as
        // a pair. Otherwise the standard SE-solo path (barcode on a separate read).
        if params.solo_barcode_on_mate1() {
            align_reads_solo_pe(params, index, writer.as_mut(), &stats, &sj_stats, sctx)?;
        } else {
            align_reads_solo(params, index, writer.as_mut(), &stats, &sj_stats, sctx)?;
        }
        writer.finish()?;
        if let Some(ref mut w) = tr_writer {
            w.finish()?;
        }
        let sj_output_path = params.output_path("SJ.out.tab");
        if !sj_stats.is_empty() {
            sj_stats.write_output(&sj_output_path, &index.genome, params)?;
        }
        // Per-cell count matrices (raw + filtered), Summary.csv, and the SJ
        // feature matrix — written here where sj_stats is available.
        write_solo_output(sctx, params, &stats, &sj_stats, index)?;
        stats.print_summary();
        return Ok(stats);
    }

    let n_align_files = params.read_files_in.len();
    match n_align_files {
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
        n => anyhow::bail!("Invalid number of read files: {n} (expected 1 or 2)"),
    }?;

    writer.finish()?;

    // Flush transcriptome BAM.
    if let Some(ref mut w) = tr_writer {
        w.finish()?;
    }

    // 5. Write SJ.out.tab file
    let sj_output_path = params.output_path("SJ.out.tab");
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
    solo_ctx: Option<&std::sync::Arc<crate::solo::SoloContext>>,
) -> anyhow::Result<std::sync::Arc<crate::stats::AlignmentStats>> {
    use std::sync::Arc;

    // PASS 1: Junction discovery (no quant counting in pass 1)
    info!("Two-pass mode: Pass 1 - Junction discovery");
    let (sj_stats_pass1, novel_junctions) = run_pass1(index, params)?;

    let pass1_dir = params.output_path("_STARpass1");
    std::fs::create_dir_all(&pass1_dir)?;
    let pass1_path = pass1_dir.join("SJ.out.tab");

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
    let stats = run_single_pass(&Arc::new(merged_index), params, quant_ctx, tr_idx, solo_ctx)?;

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

    // Align reads (single-end or paired-end); no quant counting in pass 1.
    // Solo runs align only the cDNA read (file 0) — route to the SE path.
    let n_align_files = if params.solo_enabled() {
        1
    } else {
        params.read_files_in.len()
    };
    match n_align_files {
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
        n => anyhow::bail!("Invalid number of read files: {n} (expected 1 or 2)"),
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

    let mut rng = crate::rng::SplitMix64::seed(per_read_seed(params.run_rng_seed, read_name));
    let primary_hit = rng.below(n_alignments);
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
    /// `--outWigType bedGraph` signal contributions: (transcript, is_second_mate)
    /// per reported alignment. Empty unless enabled and the read/pair mapped within
    /// the multimap limit. `signal_n_tr` is the shared NH for all of them
    /// (transcripts.len() for SE, both_mapped.len() for PE).
    signal_contrib: Vec<(crate::align::transcript::Transcript, bool)>,
    signal_n_tr: usize,
}

/// One `Result` per read/pair in an aligned batch.
type BatchOut<T> = Vec<Result<T, error::Error>>;

/// Batches in flight beyond the one being consumed. 2 keeps the pool fed across
/// a batch's slow-tail read while bounding in-flight memory to a few batches.
const PIPELINE_DEPTH: usize = 2;

/// Drive an ordered, bounded cross-batch alignment pipeline (shared by the SE,
/// PE, and solo read loops).
///
/// Decoded input batches arrive on `read_rx`. Each is aligned on the rayon pool
/// via `align`, with up to `depth + 1` batches in flight at once, so when one
/// batch is down to a single slow (dense-repeat) read the pool's idle workers
/// steal reads from the next batch instead of stalling at a per-batch barrier.
/// Finished batches are handed to `consume` strictly in input order; `consume`
/// returns `Ok(false)` to stop early (e.g. the downstream writer has gone away).
/// A panic inside `align` becomes an error batch rather than aborting the process
/// (rayon's default). `progress` is called with the running read count per batch.
fn run_batch_pipeline<In, Out>(
    read_rx: std::sync::mpsc::Receiver<Result<Vec<In>, error::Error>>,
    max_reads: u64,
    depth: usize,
    progress: impl Fn(u64),
    align: impl Fn(u64, Vec<In>) -> BatchOut<Out> + Clone + Send + 'static,
    mut consume: impl FnMut(BatchOut<Out>) -> anyhow::Result<bool>,
) -> anyhow::Result<()>
where
    In: Send + 'static,
    Out: Send + 'static,
{
    let mut inflight: std::collections::VecDeque<std::sync::mpsc::Receiver<BatchOut<Out>>> =
        std::collections::VecDeque::new();
    let mut read_count = 0u64;
    for msg in &read_rx {
        let mut batch = msg?;
        if batch.is_empty() {
            break;
        }
        let take = ((max_reads - read_count) as usize).min(batch.len());
        batch.truncate(take);

        // Batch base offset in input (FASTQ) order — threaded to `align` so per-read
        // indices (e.g. --outSAMreadID Number) are deterministic regardless of pool
        // scheduling: dispatch is sequential here, so `read_count` is this batch's base.
        let base = read_count;
        let align = align.clone();
        let (tx, rx) = std::sync::mpsc::sync_channel::<BatchOut<Out>>(1);
        rayon::spawn(move || {
            let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| align(base, batch)))
                .unwrap_or_else(|_| {
                    vec![Err(error::Error::Alignment(
                        "alignment task panicked".to_string(),
                    ))]
                });
            let _ = tx.send(out);
        });
        inflight.push_back(rx);

        read_count += take as u64;
        progress(read_count);

        // Forward completed batches (oldest first) once enough are in flight.
        // Blocks only this dispatcher, never the pool.
        while inflight.len() > depth {
            match inflight.pop_front().unwrap().recv() {
                Ok(done) => {
                    if !consume(done)? {
                        return Ok(());
                    }
                }
                Err(_) => return Ok(()),
            }
        }
        if read_count >= max_reads {
            break;
        }
    }
    // Owned so it is dropped here: closes the decode channel, so a producer blocked
    // on `--readMapNumber` early-exit unblocks and the caller's scope can join.
    drop(read_rx);
    // Drain the remaining in-flight batches in input order.
    while let Some(rx) = inflight.pop_front() {
        match rx.recv() {
            Ok(done) => {
                if !consume(done)? {
                    break;
                }
            }
            Err(_) => break,
        }
    }
    Ok(())
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
    for r in &mut rec1s {
        *r.flags_mut() |= Flags::SEGMENTED | Flags::FIRST_SEGMENT;
    }
    for r in &mut rec2s {
        *r.flags_mut() |= Flags::SEGMENTED | Flags::LAST_SEGMENT;
    }

    for (((r1, r2), p1), p2) in rec1s
        .iter_mut()
        .zip(rec2s.iter_mut())
        .zip(p1s.iter())
        .zip(p2s.iter())
    {
        crate::io::sam::apply_pe_transcriptome_mate_fields(r1, r2, p1, p2)?;
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
    use cigar::op::Kind;

    let scorer = AlignmentScorer::from_params_minimal();
    let mut keys = Vec::new();
    let mut genome_pos = transcript.genome_start;

    for op in &transcript.cigar {
        match op.kind() {
            Kind::Skip => {
                let intron_len = op.len();
                let intron_start = genome_pos;
                let intron_end = genome_pos + intron_len as u64 - 1;

                let motif =
                    scorer.detect_splice_motif(genome_pos, intron_len as u32, &index.genome);
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
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                genome_pos += op.len() as u64;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }

    keys
}

/// Write the four `--outWigType bedGraph` stranded tracks (STAR's `Signal.{Unique,
/// UniqueMultiple}.str{1,2}.out.bg`). No-op if `signal` is `None` (feature disabled).
/// RPM scales the `Unique` tracks by `1e6/nUniq` and `UniqueMultiple` by
/// `1e6/(nUniq+nMult)` (STAR `signalFromBAM` `normFactor[0]` vs `[1]`); one norm
/// for both would double `UniqueMultiple` when multimappers exist.
fn write_signal_tracks(
    signal: Option<&crate::signal::Signal>,
    sig_n_uniq: u64,
    sig_n_mult: u64,
    params: &Parameters,
) -> anyhow::Result<()> {
    let Some(sig) = signal else {
        return Ok(());
    };
    let norm_unique = if sig_n_uniq > 0 {
        1.0e6 / sig_n_uniq as f64
    } else {
        0.0
    };
    let norm_mult = if sig_n_uniq + sig_n_mult > 0 {
        1.0e6 / (sig_n_uniq + sig_n_mult) as f64
    } else {
        0.0
    };
    for (unique, strand, file) in [
        (true, 0, "Signal.Unique.str1.out.bg"),
        (false, 0, "Signal.UniqueMultiple.str1.out.bg"),
        (true, 1, "Signal.Unique.str2.out.bg"),
        (false, 1, "Signal.UniqueMultiple.str2.out.bg"),
    ] {
        let body = if params.out_wig_rpm() {
            let norm = if unique { norm_unique } else { norm_mult };
            sig.bedgraph_rpm(unique, strand, norm)
        } else {
            sig.bedgraph(unique, strand)
        };
        std::fs::write(params.output_path(file), body)?;
    }
    Ok(())
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
    tr_writer: Option<&mut crate::io::bam::BamWriter>,
    unmapped_writer: Option<&mut crate::io::fastq::UnmappedFastqWriter>,
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

    let reader = FastqReader::open(read_file, params.read_files_command.as_deref())?;

    // Create chimeric output writer if enabled
    let chimeric_writer = if params.chim_segment_min > 0 && params.chim_out_junctions() {
        use crate::chimeric::ChimericJunctionWriter;
        info!(
            "Chimeric detection enabled (chimSegmentMin={})",
            params.chim_segment_min
        );
        Some(ChimericJunctionWriter::new(&params.out_file_name_prefix)?)
    } else {
        None
    };

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };

    let batch_size = 10000;
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    // `--outSAMtype None` (e.g. quant-only) skips building SAM records.
    let emit_sam = params.emits_alignments();
    let output_unmapped = emit_sam && params.out_sam_unmapped != params::OutSamUnmapped::None;
    let write_unmapped_fastq = params.out_reads_unmapped == params::OutReadsUnmapped::Fastx;
    let by_sjout = params.out_filter_type == OutFilterType::BySJout;

    // BySJout disk buffer: SAM records written to a temp file; only compact metadata kept in RAM.
    // For 100M reads this avoids ~60 GB of Vec<RecordBuf> in memory.
    let bysj_temp = if by_sjout {
        info!("outFilterType=BySJout: disk-buffering reads for post-alignment junction filtering");
        let tf = tempfile::NamedTempFile::new()
            .map_err(|e| anyhow::anyhow!("BySJout: failed to create temp file: {e}"))?;
        Some(tf)
    } else {
        None
    };
    let (bysj_sam_header, bysj_temp_writer) = if let Some(ref tf) = bysj_temp {
        let write_file = tf
            .reopen()
            .map_err(|e| anyhow::anyhow!("BySJout: temp file reopen error: {e}"))?;
        let (hdr, w) = crate::io::sam::create_bysj_writer(write_file, &index.genome, params)?;
        (Some(hdr), Some(w))
    } else {
        (None, None)
    };
    let bysj_meta: Vec<BySJReadMeta> = Vec::new();

    info!("Aligning reads...");
    // Three-stage pipeline: a producer thread decodes the next FASTQ batch, rayon
    // aligns the current batch in parallel, and a dedicated writer thread serializes
    // SAM/BAM output — so gzip inflate, alignment, and record encoding all overlap
    // instead of running one-after-another per batch. Bounded channels (depth 2) give
    // backpressure. Output order is preserved: aligned batches flow through the
    // channel in input order and the writer consumes them in that order.
    let stats_writer = Arc::clone(&stats);
    let sj_stats_writer = Arc::clone(&sj_stats);
    let index_writer = Arc::clone(index);
    // Shared, 'static parameters for the per-batch aligner tasks spawned below.
    let params_arc = Arc::new(params.clone());
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let (read_tx, read_rx) = std::sync::mpsc::sync_channel::<
            Result<Vec<crate::io::fastq::EncodedRead>, error::Error>,
        >(2);
        #[allow(clippy::type_complexity)]
        let (res_tx, res_rx) =
            std::sync::mpsc::sync_channel::<Vec<Result<AlignmentBatchResults, error::Error>>>(2);

        // Stage 1: decode.
        scope.spawn(move || {
            let mut reader = reader;
            loop {
                match reader.read_batch(batch_size) {
                    Ok(batch) => {
                        let last = batch.is_empty();
                        if read_tx.send(Ok(batch)).is_err() || last {
                            break;
                        }
                    }
                    Err(e) => {
                        let _ = read_tx.send(Err(e));
                        break;
                    }
                }
            }
        });

        // Stage 3: writer. Owns every output-side handle for its whole lifetime so
        // the main thread never touches the writers.
        let writer_handle = scope.spawn(move || -> anyhow::Result<()> {
            let stats = stats_writer;
            let sj_stats = sj_stats_writer;
            let index = index_writer;
            let writer = writer;
            let mut tr_writer = tr_writer;
            let mut unmapped_writer = unmapped_writer;
            let mut chimeric_writer = chimeric_writer;
            let mut bysj_temp_writer = bysj_temp_writer;
            let mut bysj_meta = bysj_meta;
            // --outWigType bedGraph: the coverage signal + RPM counters are folded
            // here (sequential, in chunk order) — the Vec<f64> accumulator isn't safe
            // to update from the parallel per-read closures. SE: 1 record/read.
            let mut signal = params.out_wig_bedgraph().then(|| {
                crate::signal::Signal::new(&index.genome.chr_name, &index.genome.chr_length)
            });
            let (mut sig_n_uniq, mut sig_n_mult) = (0u64, 0u64);
            for batch_results in &res_rx {
                if by_sjout {
                    for result in batch_results {
                        let batch = result?;
                        // --outWigType bedGraph: fold this read's contribution.
                        if let Some(sig) = signal.as_mut()
                            && !batch.signal_contrib.is_empty()
                        {
                            if batch.signal_n_tr == 1 {
                                sig_n_uniq += 1;
                            } else {
                                sig_n_mult += 1;
                            }
                            for (tr, second_mate) in &batch.signal_contrib {
                                sig.add_transcript(
                                    &index.genome,
                                    tr,
                                    batch.signal_n_tr,
                                    *second_mate,
                                );
                            }
                        }
                        // Write SAM records to temp file (disk, not RAM)
                        let n_sam_records = batch.sam_records.records.len() as u32;
                        if let (Some(tw), Some(hdr)) = (&mut bysj_temp_writer, &bysj_sam_header) {
                            crate::io::sam::bysj_write_records(
                                tw,
                                hdr,
                                &batch.sam_records.records,
                            )?;
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

                        // --outWigType bedGraph: fold this read's contribution.
                        if let Some(sig) = signal.as_mut()
                            && !batch.signal_contrib.is_empty()
                        {
                            if batch.signal_n_tr == 1 {
                                sig_n_uniq += 1;
                            } else {
                                sig_n_mult += 1;
                            }
                            for (tr, second_mate) in &batch.signal_contrib {
                                sig.add_transcript(
                                    &index.genome,
                                    tr,
                                    batch.signal_n_tr,
                                    *second_mate,
                                );
                            }
                        }

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
                    let read_file = tf.reopen().map_err(|e| {
                        anyhow::anyhow!("BySJout: temp file reopen for reading: {e}")
                    })?;
                    let mut reader =
                        noodles::sam::io::Reader::new(std::io::BufReader::new(read_file));
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
                                    let supp =
                                        build_within_bam_records(chim_aln, &index.genome, 255)?;
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

                info!("BySJout: filtered {filtered_count} reads with non-surviving junctions");
            }

            // --outWigType bedGraph: write the four Signal.*.out.bg tracks.
            write_signal_tracks(signal.as_ref(), sig_n_uniq, sig_n_mult, params)?;

            // Flush chimeric output if enabled
            if let Some(ref mut chim_writer) = chimeric_writer {
                // --chimOutJunctionFormat 1: STAR-Fusion comment trailer with read counts
                // (# Nreads <total>\tNreadsUnique <uniquely_mapped>\tNreadsMulti <multi_mapped>).
                if params.chim_out_junction_format == 1 {
                    use std::sync::atomic::Ordering;
                    let command_line = params.command_line.as_deref().unwrap_or("");
                    let n_reads = stats.total_reads.load(Ordering::Relaxed);
                    let n_unique = stats.uniquely_mapped.load(Ordering::Relaxed);
                    let n_multi = stats.multi_mapped.load(Ordering::Relaxed);
                    chim_writer.write_format1_trailer(command_line, n_reads, n_unique, n_multi)?;
                }
                chim_writer.flush()?;
                info!("Chimeric junction output complete");
            }

            // Flush unmapped FASTQ writer
            if let Some(ref mut uw) = unmapped_writer {
                uw.flush()?;
            }
            Ok(())
        });

        // WASP allele-specific filtering: load the VCF once (shared read-only across
        // the parallel per-read closures). `None` unless --waspOutputMode SAMtag.
        let wasp_ctx: Arc<Option<crate::wasp::WaspContext>> = Arc::new(
            if params.wasp_output_mode == params::WaspOutputMode::SAMtag {
                let vcf = params
                    .var_vcf_file
                    .as_ref()
                    .expect("validated: SAMtag requires --varVCFfile");
                let ctx = crate::wasp::WaspContext::load(
                    vcf,
                    &index.genome.chr_name,
                    &index.genome.chr_start,
                    params,
                )
                .map_err(|source| error::Error::Io {
                    source,
                    path: vcf.clone(),
                })?;
                info!(
                    "WASP: loaded {} heterozygous SNVs from {}",
                    ctx.snps.len(),
                    vcf.display()
                );
                Some(ctx)
            } else {
                None
            },
        );

        // Stage 2: align each decoded batch on the rayon pool via run_batch_pipeline,
        // forwarding finished batches to the writer thread in input order.
        let align_result = {
            let index = Arc::clone(index);
            let stats = Arc::clone(&stats);
            let sj_stats = Arc::clone(&sj_stats);
            let quant = quant.as_ref().map(Arc::clone);
            let tr = tr.as_ref().map(Arc::clone);
            let params_arc = Arc::clone(&params_arc);
            let wasp_ctx = Arc::clone(&wasp_ctx);
            let align = move |base: u64,
                              batch: Vec<crate::io::fastq::EncodedRead>|
                  -> BatchOut<AlignmentBatchResults> {
                let params: &Parameters = &params_arc;
                // Adapter-aware clip params (fixed 5'/3' Nbases + 3' adapter Hamming
                // scan); built once per batch, applied per read via clip_mate.
                let clip_params = crate::clip::clip_params_from(params, 0);
                batch
                    .par_iter()
                    .enumerate()
                    .map(|(read_idx, read)| {
                        // --outSAMreadID Number: replace the output QNAME with the
                        // read's 1-based input index (deterministic — from the FASTQ
                        // order via `base`, not parallel execution order). The seed
                        // name passed to align_read is left as the real read name.
                        let out_read_name =
                            if params.out_sam_read_id == crate::params::OutSamReadId::Number {
                                (base + read_idx as u64 + 1).to_string()
                            } else {
                                read.name.clone()
                            };
                        #[allow(clippy::needless_borrow)]
                        let index = Arc::clone(&index);
                        #[allow(clippy::needless_borrow)]
                        let stats = Arc::clone(&stats);
                        #[allow(clippy::needless_borrow)]
                        let sj_stats = Arc::clone(&sj_stats);
                        let quant = quant.as_ref().map(Arc::clone);

                        // Apply clipping: fixed Nbases + 3' adapter (clip_mate), then trim.
                        let (clip5p, clip3p) = crate::clip::clip_mate(&read.sequence, &clip_params);
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
                                // Unmapped reads keep the full original read (STAR: clipped
                                // bases are never removed from an unmapped record's SEQ).
                                let record = SamWriter::build_unmapped_record(
                                    &out_read_name,
                                    &read.sequence,
                                    &read.quality,
                                    params,
                                    crate::stats::UnmappedReason::Other,
                                )?;
                                buffer.push(record);
                            }
                            let unmapped_m1 = if write_unmapped_fastq {
                                vec![(
                                    out_read_name.clone(),
                                    read.sequence.clone(),
                                    read.quality.clone(),
                                )]
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
                                signal_contrib: Vec::new(),
                                signal_n_tr: 0,
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

                        // Record junction statistics (per-read dedup, fix A)
                        let is_unique = transcripts.len() == 1;
                        record_read_junctions(&transcripts, &index, &sj_stats, is_unique);

                        // Extract junction keys from primary alignment for BySJout filtering
                        let primary_junction_keys =
                            if by_sjout && !transcripts.is_empty() && transcripts[0].n_junction > 0
                            {
                                extract_junction_keys(&transcripts[0], &index)
                            } else {
                                Vec::new()
                            };

                        // --outWigType bedGraph: this read's coverage contribution is
                        // its reported alignments (same set written to SAM: mapped and
                        // within the multimap limit). signal_n_tr is the shared NH.
                        let signal_n_tr = transcripts.len();
                        let signal_contrib: Vec<(crate::align::transcript::Transcript, bool)> =
                            if params.out_wig_bedgraph()
                                && !transcripts.is_empty()
                                && transcripts.len() <= max_multimaps
                            {
                                transcripts.iter().cloned().map(|t| (t, false)).collect()
                            } else {
                                Vec::new()
                            };

                        // Build SAM records (no I/O, just construction).
                        // Skipped entirely under `--outSAMtype None`.
                        let is_unmapped_se = transcripts.is_empty();
                        if emit_sam {
                            if is_unmapped_se {
                                // Unmapped
                                if output_unmapped {
                                    // Full original read for unmapped (STAR convention).
                                    let record = SamWriter::build_unmapped_record(
                                        &out_read_name,
                                        &read.sequence,
                                        &read.quality,
                                        params,
                                        unmapped_reason
                                            .unwrap_or(crate::stats::UnmappedReason::Other),
                                    )?;
                                    buffer.push(record);
                                }
                            } else if transcripts.len() <= max_multimaps {
                                // Mapped (within multimap limit)
                                let mut records = SamWriter::build_alignment_records(
                                    &out_read_name,
                                    &read.sequence,
                                    &read.quality,
                                    clip5p,
                                    clip3p,
                                    &transcripts,
                                    &index.genome,
                                    params,
                                    n_for_mapq,
                                )?;
                                // WASP allele-specific filtering: stamp vW/vA/vG by
                                // re-mapping the allele-swapped read (STAR waspMap).
                                if let Some(ctx) = &*wasp_ctx {
                                    crate::wasp::annotate_records_se(
                                        &mut records,
                                        &transcripts,
                                        &clipped_seq,
                                        &read.name,
                                        &index,
                                        ctx,
                                        params.out_sam_attributes,
                                    )?;
                                }
                                for record in records {
                                    buffer.push(record);
                                }
                            }
                            // else: too many loci, skip output
                        }

                        // Transcriptome SAM projection for --quantMode TranscriptomeSAM.
                        let transcriptome_records: Vec<
                            noodles::sam::alignment::record_buf::RecordBuf,
                        > = if let Some(ref tidx) = tr_local {
                            build_transcriptome_records_se(
                                &transcripts,
                                &out_read_name,
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
                            vec![(
                                out_read_name.clone(),
                                read.sequence.clone(),
                                read.quality.clone(),
                            )]
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
                            signal_contrib,
                            signal_n_tr,
                        })
                    })
                    .collect()
            };
            run_batch_pipeline(
                read_rx,
                max_reads,
                PIPELINE_DEPTH,
                |n: u64| {
                    if n % 100_000 < batch_size as u64 {
                        info!("Processed {n} reads...");
                    }
                },
                align,
                |done| Ok(res_tx.send(done).is_ok()),
            )
        };
        // Disconnect the writer channel so the writer thread can finish and join.
        drop(res_tx);
        let writer_result = writer_handle
            .join()
            .map_err(|_| anyhow::anyhow!("SE writer thread panicked"))?;
        align_result?;
        writer_result?;
        Ok(())
    })?;

    Ok(())
}

/// Align a STARsolo single-cell run: the cDNA read (file 0) is aligned exactly
/// like the SE path, while the barcode read (file 1) is read in lockstep and
/// quantified per cell. Mapped cDNA alignments are written to the SAM/BAM output
/// just like a normal SE run; the per-cell (CB, UMI, gene) records are collected
/// into `solo_ctx.recorder` for the matrix output that follows in Phase 14.4.
///
/// Solo runs are single-pass and (for now) do not support BySJout / chimeric /
/// transcriptome-SAM side outputs — those are not part of the STARsolo MVP.
fn align_reads_solo<W: AlignmentWriter + ?Sized>(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut W,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
    solo_ctx: &std::sync::Arc<crate::solo::SoloContext>,
) -> anyhow::Result<()> {
    use crate::align::read_align::align_read;
    use crate::io::fastq::clip_read;
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::solo::{SoloCountRecord, SoloMultiRecord};
    use rayon::prelude::*;
    use std::sync::Arc;

    let cdna_file = &params.read_files_in[0];
    let barcode_file = &params.read_files_in[1];
    info!(
        "STARsolo: cDNA reads from {}, barcode reads from {}",
        cdna_file.display(),
        barcode_file.display()
    );
    let reader = crate::solo::open_reader(params)?;

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
    let solo = Arc::clone(solo_ctx);

    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };
    let batch_size = 10000;
    let clip5p = params.clip5p(0);
    let clip3p = params.clip3p(0);
    let cr4_clip = params.clip_adapter_type == "CellRanger4";
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    // With `--outSAMtype None` (count-only) we skip building SAM records entirely
    // — a large saving for solo runs that only need the count matrix.
    let emit_sam = params.emits_alignments();
    let output_unmapped = emit_sam && params.out_sam_unmapped != params::OutSamUnmapped::None;
    // Shared, 'static parameters for the per-batch aligner tasks spawned below.
    let params_arc = Arc::new(params.clone());

    /// Per-read result for the solo loop (one outcome per quantified feature).
    struct SoloReadProduct {
        sam_records: BufferedSamRecords,
        per_feature: Vec<crate::solo::FeatureOutcome>,
        sj: Vec<crate::solo::SjCountRecord>,
        velocyto: Option<crate::solo::VelocytoRecord>,
    }

    info!("STARsolo: aligning cDNA reads and quantifying barcodes...");
    // Decode the next batch on a background thread while the current batch is
    // aligning — overlaps the single-threaded gzip inflate (the per-batch serial
    // section that otherwise stalls every rayon worker) with the parallel work.
    // A bounded channel (depth 2) supplies backpressure so the reader stays at
    // most one batch ahead.
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let (tx, rx) =
            std::sync::mpsc::sync_channel::<Result<Vec<crate::solo::SoloRead>, error::Error>>(2);
        scope.spawn(move || {
            let mut reader = reader;
            loop {
                match reader.read_batch(batch_size) {
                    Ok(batch) => {
                        let last = batch.is_empty();
                        if tx.send(Ok(batch)).is_err() || last {
                            break;
                        }
                    }
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        break;
                    }
                }
            }
        });

        // Run the consumer in an IIFE so that on ANY exit path (EOF, a
        // `--readMapNumber` early break, or a `?` error) we still reach the
        // `drop(rx)` below. Disconnecting the receiver wakes a producer that is
        // blocked on the full bounded channel (its `send` returns Err, so it
        // breaks) — otherwise an early exit would deadlock joining the scope.
        // The in-order consumer (SAM write + recorder accumulation) runs on this
        // dispatcher thread — which owns `writer` — so SAM output and matrix order
        // are identical to the serial version.
        {
            let consume = |products: BatchOut<SoloReadProduct>| -> anyhow::Result<bool> {
                let n_feat = solo.features.len();
                let mut feat_records: Vec<Vec<SoloCountRecord>> =
                    (0..n_feat).map(|_| Vec::new()).collect();
                let mut feat_multi: Vec<Vec<SoloMultiRecord>> =
                    (0..n_feat).map(|_| Vec::new()).collect();
                let mut feat_multi_gene: Vec<Vec<crate::solo::MultiGeneRecord>> =
                    (0..n_feat).map(|_| Vec::new()).collect();
                let mut sj_batch: Vec<crate::solo::SjCountRecord> = Vec::new();
                let mut velo_batch: Vec<crate::solo::VelocytoRecord> = Vec::new();
                for result in products {
                    let product = result?;
                    writer.write_batch(&product.sam_records.records)?;
                    for (fi, fo) in product.per_feature.into_iter().enumerate() {
                        if let Some(r) = fo.record {
                            feat_records[fi].push(r);
                        }
                        if let Some(m) = fo.multi {
                            feat_multi[fi].push(m);
                        }
                        if let Some(mg) = fo.multi_gene {
                            feat_multi_gene[fi].push(mg);
                        }
                    }
                    sj_batch.extend(product.sj);
                    velo_batch.extend(product.velocyto);
                }
                for (fi, recorder) in solo.recorders.iter().enumerate() {
                    recorder.extend(
                        std::mem::take(&mut feat_records[fi]),
                        std::mem::take(&mut feat_multi[fi]),
                    );
                    let mg = std::mem::take(&mut feat_multi_gene[fi]);
                    if !mg.is_empty() {
                        recorder.multi_gene.lock().unwrap().extend(mg);
                    }
                }
                if !sj_batch.is_empty() {
                    solo.sj_records.lock().unwrap().extend(sj_batch);
                }
                if !velo_batch.is_empty() {
                    solo.velocyto_records.lock().unwrap().extend(velo_batch);
                }
                Ok(true)
            };
            // Cloned into a fresh nested scope so the outer `solo` stays available to
            // `consume` above; the align body clones from `&Arc` (fn-param style).
            let align = {
                let index = Arc::clone(index);
                let stats = Arc::clone(&stats);
                let sj_stats = Arc::clone(&sj_stats);
                let solo = Arc::clone(&solo);
                let params_arc = Arc::clone(&params_arc);
                move |base: u64, batch: Vec<crate::solo::SoloRead>| -> BatchOut<SoloReadProduct> {
                    let params: &Parameters = &params_arc;
                    let index = &index;
                    batch
                        .par_iter()
                        .enumerate()
                        .map(|(read_idx, sread)| {
                            let index = Arc::clone(index);
                            let stats = Arc::clone(&stats);
                            let sj_stats = Arc::clone(&sj_stats);
                            let solo = Arc::clone(&solo);

                            let read = &sread.cdna;
                            // --outSAMreadID Number: numeric output QNAME (1-based input order).
                            let out_read_name =
                                if params.out_sam_read_id == crate::params::OutSamReadId::Number {
                                    (base + read_idx as u64 + 1).to_string()
                                } else {
                                    read.name.clone()
                                };
                            // CellRanger4 adapter clipping (TSO 5' + polyA 3') runs before
                            // the fixed clip5p/clip3p Nbases trimming.
                            let (cr_seq, cr_qual, cr4_5p, cr4_3p) = if cr4_clip {
                                crate::solo::clip_adapter_cr4(&read.sequence, &read.quality)
                            } else {
                                (read.sequence.clone(), read.quality.clone(), 0, 0)
                            };
                            let (clipped_seq, clipped_qual) =
                                clip_read(&cr_seq, &cr_qual, clip5p, clip3p);
                            // Total clip against the ORIGINAL read = CR4 (TSO/polyA) + fixed
                            // Nbases; used to soft-clip all trimmed bases (STARsolo convention).
                            let total_clip5p = cr4_5p + clip5p;
                            let total_clip3p = cr4_3p + clip3p;
                            let mut buffer = BufferedSamRecords::new();
                            stats.record_read_bases(clipped_seq.len() as u64);

                            if clipped_seq.is_empty() {
                                stats.record_alignment(0, max_multimaps);
                                stats.record_unmapped_reason(crate::stats::UnmappedReason::Other);
                                // No alignment → barcode still counts toward stats (unmapped → no gene).
                                let outcome = solo.process_read(
                                    &[],
                                    0,
                                    sread.barcode.as_ref(),
                                    &[],
                                    &read.quality,
                                );
                                return Ok(SoloReadProduct {
                                    sam_records: buffer,
                                    per_feature: outcome.per_feature,
                                    sj: outcome.sj,
                                    velocyto: outcome.velocyto,
                                });
                            }

                            let (transcripts, _chimeric, n_for_mapq, unmapped_reason) =
                                align_read(&clipped_seq, &read.name, &index, params)?;

                            let n_for_stats = if transcripts.is_empty() && n_for_mapq > 0 {
                                n_for_mapq
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

                            let is_unique = transcripts.len() == 1;
                            record_read_junctions(&transcripts, &index, &sj_stats, is_unique);

                            // SJ feature: the junctions crossed by a uniquely-mapped read
                            // (absolute intron coords), mapped to SJ.out.tab rows at output.
                            let junctions: Vec<(u64, u64)> =
                                if solo.sj_enabled && is_unique && transcripts[0].n_junction > 0 {
                                    extract_junction_keys(&transcripts[0], &index)
                                        .into_iter()
                                        .map(|k| (k.intron_start, k.intron_end))
                                        .collect()
                                } else {
                                    Vec::new()
                                };

                            // Solo quantification (CB match + UMI check + gene assignment).
                            let outcome = solo.process_read(
                                &transcripts,
                                transcripts.len(),
                                sread.barcode.as_ref(),
                                &junctions,
                                &read.quality,
                            );

                            // Build SAM records for the cDNA alignment (same as SE path).
                            // Skipped entirely under `--outSAMtype None` (count-only).
                            if emit_sam {
                                if transcripts.is_empty() {
                                    if output_unmapped {
                                        let record = SamWriter::build_unmapped_record(
                                            &out_read_name,
                                            &clipped_seq,
                                            &clipped_qual,
                                            params,
                                            unmapped_reason
                                                .unwrap_or(crate::stats::UnmappedReason::Other),
                                        )?;
                                        buffer.push(record);
                                    }
                                } else if transcripts.len() <= max_multimaps {
                                    // Soft-clip ALL trimmed bases (CR4 TSO/polyA + fixed
                                    // clip5p/clip3p) against the original read, matching
                                    // STARsolo (e.g. 60M30S). Inert at default 10x.
                                    let mut records = SamWriter::build_alignment_records(
                                        &out_read_name,
                                        &read.sequence,
                                        &read.quality,
                                        total_clip5p,
                                        total_clip3p,
                                        &transcripts,
                                        &index.genome,
                                        params,
                                        n_for_mapq,
                                    )?;
                                    // STARsolo GX/GN gene tags (Gene-feature assignment).
                                    if params.out_sam_attributes.intersects(
                                        crate::params::SamAttributes::GX
                                            | crate::params::SamAttributes::GN,
                                    ) {
                                        let (gx, gn) = solo.gene_tags(&transcripts);
                                        crate::io::sam::add_gene_tags(
                                            &mut records,
                                            gx,
                                            gn,
                                            params.out_sam_attributes,
                                        );
                                    }
                                    for record in records {
                                        buffer.push(record);
                                    }
                                }
                            }

                            Ok(SoloReadProduct {
                                sam_records: buffer,
                                per_feature: outcome.per_feature,
                                sj: outcome.sj,
                                velocyto: outcome.velocyto,
                            })
                        })
                        .collect()
                }
            };
            run_batch_pipeline(
                rx,
                max_reads,
                PIPELINE_DEPTH,
                |n: u64| {
                    if n % 100_000 < batch_size as u64 {
                        info!("STARsolo: processed {n} reads...");
                    }
                },
                align,
                consume,
            )
        }
    })?;

    Ok(())
}

/// Align a 5' paired-end STARsolo run (`--soloBarcodeMate 1`): the cell barcode +
/// UMI are read from mate 1, and both mates (mate 1's cDNA after the clipped
/// barcode region + mate 2) are aligned as a pair, then quantified per cell.
/// Mirrors [`align_reads_solo`] (same cross-batch pipeline) but drives the paired
/// aligner and `SoloContext::process_read_pe`.
fn align_reads_solo_pe<W: AlignmentWriter + ?Sized>(
    params: &Parameters,
    index: &std::sync::Arc<crate::index::GenomeIndex>,
    writer: &mut W,
    stats: &std::sync::Arc<crate::stats::AlignmentStats>,
    sj_stats: &std::sync::Arc<crate::junction::SpliceJunctionStats>,
    solo_ctx: &std::sync::Arc<crate::solo::SoloContext>,
) -> anyhow::Result<()> {
    use crate::align::read_align::{PairedAlignment, PairedAlignmentResult, align_paired_read};
    use crate::io::fastq::clip_read;
    use crate::io::sam::{BufferedSamRecords, SamWriter};
    use crate::solo::{SoloCountRecord, SoloMultiRecord};
    use rayon::prelude::*;
    use std::sync::Arc;

    let (m1_file, m2_file) = params.solo_cdna_mate_files().ok_or_else(|| {
        anyhow::anyhow!("solo: --soloBarcodeMate 1 requires two --readFilesIn cDNA mate files")
    })?;
    info!(
        "STARsolo (5' paired-end): mate 1 (barcode + cDNA) from {}, mate 2 (cDNA) from {}",
        m1_file.display(),
        m2_file.display()
    );
    let reader = crate::solo::open_paired_reader(params)?;

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
    let solo = Arc::clone(solo_ctx);

    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };
    let batch_size = 10000;
    // Per-mate clip: mate 1 (--clip5pNbases[0], e.g. 39 to strip the 5' barcode
    // region) and mate 2 ([1], e.g. 0). CellRanger4 adapter clipping is not used
    // by the cellgeni 5' path (it uses clip5pNbases instead), so it is not applied.
    let (clip5p_m1, clip3p_m1) = (params.clip5p(0), params.clip3p(0));
    let (clip5p_m2, clip3p_m2) = (params.clip5p(1), params.clip3p(1));
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    let emit_sam = params.emits_alignments();
    let output_unmapped = emit_sam && params.out_sam_unmapped != params::OutSamUnmapped::None;
    let params_arc = Arc::new(params.clone());

    struct SoloReadProduct {
        sam_records: BufferedSamRecords,
        per_feature: Vec<crate::solo::FeatureOutcome>,
        sj: Vec<crate::solo::SjCountRecord>,
        velocyto: Option<crate::solo::VelocytoRecord>,
    }

    info!("STARsolo: aligning 5' paired-end reads and quantifying barcodes...");
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let (tx, rx) = std::sync::mpsc::sync_channel::<
            Result<Vec<crate::solo::SoloPairedRead>, error::Error>,
        >(2);
        scope.spawn(move || {
            let mut reader = reader;
            loop {
                match reader.read_batch(batch_size) {
                    Ok(batch) => {
                        let last = batch.is_empty();
                        if tx.send(Ok(batch)).is_err() || last {
                            break;
                        }
                    }
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        break;
                    }
                }
            }
        });

        {
            let consume = |products: BatchOut<SoloReadProduct>| -> anyhow::Result<bool> {
                let n_feat = solo.features.len();
                let mut feat_records: Vec<Vec<SoloCountRecord>> =
                    (0..n_feat).map(|_| Vec::new()).collect();
                let mut feat_multi: Vec<Vec<SoloMultiRecord>> =
                    (0..n_feat).map(|_| Vec::new()).collect();
                let mut feat_multi_gene: Vec<Vec<crate::solo::MultiGeneRecord>> =
                    (0..n_feat).map(|_| Vec::new()).collect();
                let mut sj_batch: Vec<crate::solo::SjCountRecord> = Vec::new();
                let mut velo_batch: Vec<crate::solo::VelocytoRecord> = Vec::new();
                for result in products {
                    let product = result?;
                    writer.write_batch(&product.sam_records.records)?;
                    for (fi, fo) in product.per_feature.into_iter().enumerate() {
                        if let Some(r) = fo.record {
                            feat_records[fi].push(r);
                        }
                        if let Some(m) = fo.multi {
                            feat_multi[fi].push(m);
                        }
                        if let Some(mg) = fo.multi_gene {
                            feat_multi_gene[fi].push(mg);
                        }
                    }
                    sj_batch.extend(product.sj);
                    velo_batch.extend(product.velocyto);
                }
                for (fi, recorder) in solo.recorders.iter().enumerate() {
                    recorder.extend(
                        std::mem::take(&mut feat_records[fi]),
                        std::mem::take(&mut feat_multi[fi]),
                    );
                    let mg = std::mem::take(&mut feat_multi_gene[fi]);
                    if !mg.is_empty() {
                        recorder.multi_gene.lock().unwrap().extend(mg);
                    }
                }
                if !sj_batch.is_empty() {
                    solo.sj_records.lock().unwrap().extend(sj_batch);
                }
                if !velo_batch.is_empty() {
                    solo.velocyto_records.lock().unwrap().extend(velo_batch);
                }
                Ok(true)
            };
            let align = {
                let index = Arc::clone(index);
                let stats = Arc::clone(&stats);
                let sj_stats = Arc::clone(&sj_stats);
                let solo = Arc::clone(&solo);
                let params_arc = Arc::clone(&params_arc);
                move |base: u64,
                      batch: Vec<crate::solo::SoloPairedRead>|
                      -> BatchOut<SoloReadProduct> {
                    let params: &Parameters = &params_arc;
                    let index = &index;
                    batch
                        .par_iter()
                        .enumerate()
                        .map(|(pair_idx, pread)| {
                            let index = Arc::clone(index);
                            let stats = Arc::clone(&stats);
                            let sj_stats = Arc::clone(&sj_stats);
                            let solo = Arc::clone(&solo);
                            // --outSAMreadID Number: numeric output QNAME (1-based input order).
                            let out_read_name =
                                if params.out_sam_read_id == crate::params::OutSamReadId::Number {
                                    (base + pair_idx as u64 + 1).to_string()
                                } else {
                                    pread.mate1.name.clone()
                                };

                            let (m1_seq, m1_qual) = clip_read(
                                &pread.mate1.sequence,
                                &pread.mate1.quality,
                                clip5p_m1,
                                clip3p_m1,
                            );
                            let (m2_seq, m2_qual) = clip_read(
                                &pread.mate2.sequence,
                                &pread.mate2.quality,
                                clip5p_m2,
                                clip3p_m2,
                            );
                            let mut buffer = BufferedSamRecords::new();
                            stats.record_read_bases((m1_seq.len() + m2_seq.len()) as u64);

                            let (results, _pe_chimeric, n_for_mapq, unmapped_reason) =
                                align_paired_read(
                                    &m1_seq,
                                    &m2_seq,
                                    &pread.mate1.name,
                                    &index,
                                    params,
                                )?;

                            let has_half_mapped = results
                                .iter()
                                .any(|r| matches!(r, PairedAlignmentResult::HalfMapped { .. }));
                            let both_mapped: Vec<_> = results
                                .iter()
                                .filter_map(|r| match r {
                                    PairedAlignmentResult::BothMapped(pa) => Some(pa),
                                    PairedAlignmentResult::HalfMapped { .. } => None,
                                })
                                .collect();

                            if results.is_empty() {
                                stats.record_alignment(n_for_mapq, max_multimaps);
                                stats.record_unmapped_reason(
                                    unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
                                );
                            } else if has_half_mapped {
                                stats.record_alignment(1, max_multimaps);
                                stats.record_half_mapped();
                                if let Some(PairedAlignmentResult::HalfMapped {
                                    mapped_transcript,
                                    ..
                                }) = results.first()
                                {
                                    stats.record_transcript_stats(mapped_transcript);
                                }
                            } else {
                                let n = both_mapped.len();
                                stats.record_alignment(n, max_multimaps);
                                if n == 1 {
                                    stats.record_transcript_stats(&both_mapped[0].mate1_transcript);
                                    stats.record_transcript_stats(&both_mapped[0].mate2_transcript);
                                }
                            }

                            let is_unique =
                                both_mapped.len() == 1 || (has_half_mapped && results.len() == 1);
                            // Per-read junction dedup (fix A): a junction crossed by
                            // both mates / several loci counts once.
                            let mut read_trs: Vec<&crate::align::transcript::Transcript> =
                                Vec::new();
                            for result in &results {
                                match result {
                                    PairedAlignmentResult::BothMapped(pair) => {
                                        read_trs.push(&pair.mate1_transcript);
                                        read_trs.push(&pair.mate2_transcript);
                                    }
                                    PairedAlignmentResult::HalfMapped {
                                        mapped_transcript, ..
                                    } => {
                                        read_trs.push(mapped_transcript);
                                    }
                                }
                            }
                            record_read_junctions(read_trs, &index, &sj_stats, is_unique);

                            // SJ feature: junctions crossed by the unique pair (both mates).
                            let junctions: Vec<(u64, u64)> =
                                if solo.sj_enabled && both_mapped.len() == 1 {
                                    let pair = both_mapped[0];
                                    let mut js = Vec::new();
                                    for tr in [&pair.mate1_transcript, &pair.mate2_transcript] {
                                        if tr.n_junction > 0 {
                                            js.extend(
                                                extract_junction_keys(tr, &index)
                                                    .into_iter()
                                                    .map(|k| (k.intron_start, k.intron_end)),
                                            );
                                        }
                                    }
                                    js
                                } else {
                                    Vec::new()
                                };

                            // Solo quantification: union both mates (strand from mate 1)
                            // for a both-mapped pair; fall back to the mapped mate for
                            // half-mapped.
                            let outcome = if !both_mapped.is_empty() {
                                let pairs: Vec<_> = both_mapped
                                    .iter()
                                    .map(|pa| (&pa.mate1_transcript, &pa.mate2_transcript))
                                    .collect();
                                solo.process_read_pe(
                                    &pairs,
                                    pread.barcode.as_ref(),
                                    &junctions,
                                    &pread.mate1.quality,
                                )
                            } else if let Some(PairedAlignmentResult::HalfMapped {
                                mapped_transcript,
                                ..
                            }) = results.first()
                            {
                                solo.process_read(
                                    std::slice::from_ref(mapped_transcript),
                                    1,
                                    pread.barcode.as_ref(),
                                    &junctions,
                                    &pread.mate1.quality,
                                )
                            } else {
                                solo.process_read(
                                    &[],
                                    0,
                                    pread.barcode.as_ref(),
                                    &[],
                                    &pread.mate1.quality,
                                )
                            };

                            // SAM records (skipped under `--outSAMtype None`).
                            if !emit_sam {
                            } else if results.is_empty() {
                                if output_unmapped {
                                    let records = SamWriter::build_paired_unmapped_records(
                                        &out_read_name,
                                        &m1_seq,
                                        &m1_qual,
                                        &m2_seq,
                                        &m2_qual,
                                        params,
                                        unmapped_reason
                                            .unwrap_or(crate::stats::UnmappedReason::Other),
                                    )?;
                                    for record in records {
                                        buffer.push(record);
                                    }
                                }
                            } else if has_half_mapped {
                                if let Some(PairedAlignmentResult::HalfMapped {
                                    mapped_transcript,
                                    mate1_is_mapped,
                                }) = results.first()
                                {
                                    let records = SamWriter::build_half_mapped_records(
                                        &out_read_name,
                                        &pread.mate1.sequence,
                                        &pread.mate1.quality,
                                        &pread.mate2.sequence,
                                        &pread.mate2.quality,
                                        clip5p_m1,
                                        clip3p_m1,
                                        clip5p_m2,
                                        clip3p_m2,
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
                                let paired_alns: Vec<PairedAlignment> = both_mapped
                                    .iter()
                                    .map(|pa| PairedAlignment::clone(pa))
                                    .collect();
                                // Soft-clip the fixed per-mate clips against the original
                                // mate reads (matching SE/PE non-solo). Inert at default 10x.
                                let records = SamWriter::build_paired_records(
                                    &out_read_name,
                                    &pread.mate1.sequence,
                                    &pread.mate1.quality,
                                    &pread.mate2.sequence,
                                    &pread.mate2.quality,
                                    clip5p_m1,
                                    clip3p_m1,
                                    clip5p_m2,
                                    clip3p_m2,
                                    &paired_alns,
                                    &index.genome,
                                    params,
                                    n_for_mapq,
                                )?;
                                for record in records {
                                    buffer.push(record);
                                }
                            }

                            Ok(SoloReadProduct {
                                sam_records: buffer,
                                per_feature: outcome.per_feature,
                                sj: outcome.sj,
                                velocyto: outcome.velocyto,
                            })
                        })
                        .collect()
                }
            };
            run_batch_pipeline(
                rx,
                max_reads,
                PIPELINE_DEPTH,
                |n: u64| {
                    if n % 100_000 < batch_size as u64 {
                        info!("STARsolo: processed {n} read pairs...");
                    }
                },
                align,
                consume,
            )
        }
    })?;

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
    tr_writer: Option<&mut crate::io::bam::BamWriter>,
    unmapped_writer1: Option<&mut crate::io::fastq::UnmappedFastqWriter>,
    unmapped_writer2: Option<&mut crate::io::fastq::UnmappedFastqWriter>,
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

    let reader = PairedFastqReader::open(
        &params.read_files_in[0],
        &params.read_files_in[1],
        params.read_files_command.as_deref(),
    )?;

    // Create chimeric output writer if enabled
    let chimeric_writer = if params.chim_segment_min > 0 && params.chim_out_junctions() {
        use crate::chimeric::ChimericJunctionWriter;
        info!(
            "Chimeric detection enabled (chimSegmentMin={})",
            params.chim_segment_min
        );
        Some(ChimericJunctionWriter::new(&params.out_file_name_prefix)?)
    } else {
        None
    };

    let stats = Arc::clone(stats);
    let sj_stats = Arc::clone(sj_stats);
    let max_reads = if params.read_map_number < 0 {
        u64::MAX
    } else {
        params.read_map_number as u64
    };

    let batch_size = 10000;
    let max_multimaps = params.out_filter_multimap_nmax as usize;
    // `--outSAMtype None` (e.g. quant-only) skips building SAM records.
    let emit_sam = params.emits_alignments();
    let output_unmapped = emit_sam && params.out_sam_unmapped != params::OutSamUnmapped::None;
    let write_unmapped_fastq = params.out_reads_unmapped == params::OutReadsUnmapped::Fastx;
    let by_sjout = params.out_filter_type == OutFilterType::BySJout;

    // BySJout disk buffer: SAM records to temp file, compact metadata in RAM.
    let bysj_temp = if by_sjout {
        info!("outFilterType=BySJout: disk-buffering pairs for post-alignment junction filtering");
        let tf = tempfile::NamedTempFile::new()
            .map_err(|e| anyhow::anyhow!("BySJout: failed to create temp file: {e}"))?;
        Some(tf)
    } else {
        None
    };
    let (bysj_sam_header, bysj_temp_writer) = if let Some(ref tf) = bysj_temp {
        let write_file = tf
            .reopen()
            .map_err(|e| anyhow::anyhow!("BySJout: temp file reopen error: {e}"))?;
        let (hdr, w) = crate::io::sam::create_bysj_writer(write_file, &index.genome, params)?;
        (Some(hdr), Some(w))
    } else {
        (None, None)
    };
    let bysj_meta: Vec<BySJReadMeta> = Vec::new();

    info!("Aligning paired-end reads...");
    // Three-stage pipeline: producer decodes the next pair batch, rayon aligns the
    // current batch, and a dedicated writer thread serializes output — overlapping
    // gzip inflate of both mate files, alignment, and record encoding. Bounded
    // channels (depth 2) give backpressure; output order is preserved.
    let stats_writer = Arc::clone(&stats);
    let sj_stats_writer = Arc::clone(&sj_stats);
    let index_writer = Arc::clone(index);
    // Shared, 'static parameters for the per-batch aligner tasks spawned below.
    let params_arc = Arc::new(params.clone());
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let (read_tx, read_rx) = std::sync::mpsc::sync_channel::<
            Result<Vec<crate::io::fastq::PairedRead>, error::Error>,
        >(2);
        #[allow(clippy::type_complexity)]
        let (res_tx, res_rx) =
            std::sync::mpsc::sync_channel::<Vec<Result<AlignmentBatchResults, error::Error>>>(2);

        // Stage 1: decode.
        scope.spawn(move || {
            let mut reader = reader;
            loop {
                match reader.read_paired_batch(batch_size) {
                    Ok(batch) => {
                        let last = batch.is_empty();
                        if read_tx.send(Ok(batch)).is_err() || last {
                            break;
                        }
                    }
                    Err(e) => {
                        let _ = read_tx.send(Err(e));
                        break;
                    }
                }
            }
        });

        // Stage 3: writer. Owns every output-side handle for its whole lifetime.
        let writer_handle = scope.spawn(move || -> anyhow::Result<()> {
            let stats = stats_writer;
            let sj_stats = sj_stats_writer;
            let index = index_writer;
            let writer = writer;
            let mut tr_writer = tr_writer;
            let mut unmapped_writer1 = unmapped_writer1;
            let mut unmapped_writer2 = unmapped_writer2;
            let mut chimeric_writer = chimeric_writer;
            let mut bysj_temp_writer = bysj_temp_writer;
            let mut bysj_meta = bysj_meta;
            // --outWigType bedGraph: coverage signal + RPM counters, folded here in
            // chunk order. PE: 2 mate records per pair, so a uniquely-mapped pair adds
            // 2 to nUniq and a multimapped pair adds 2 to nMult (STAR signalFromBAM
            // counts BAM records, not pairs — one per reported mate alignment).
            let mut signal = params.out_wig_bedgraph().then(|| {
                crate::signal::Signal::new(&index.genome.chr_name, &index.genome.chr_length)
            });
            let (mut sig_n_uniq, mut sig_n_mult) = (0u64, 0u64);
            for batch_results in &res_rx {
                if by_sjout {
                    for result in batch_results {
                        let batch = result?;
                        if let Some(sig) = signal.as_mut()
                            && !batch.signal_contrib.is_empty()
                        {
                            if batch.signal_n_tr == 1 {
                                sig_n_uniq += 2;
                            } else {
                                sig_n_mult += 2;
                            }
                            for (tr, second_mate) in &batch.signal_contrib {
                                sig.add_transcript(
                                    &index.genome,
                                    tr,
                                    batch.signal_n_tr,
                                    *second_mate,
                                );
                            }
                        }
                        let n_sam_records = batch.sam_records.records.len() as u32;
                        if let (Some(tw), Some(hdr)) = (&mut bysj_temp_writer, &bysj_sam_header) {
                            crate::io::sam::bysj_write_records(
                                tw,
                                hdr,
                                &batch.sam_records.records,
                            )?;
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
                        if let Some(sig) = signal.as_mut()
                            && !batch.signal_contrib.is_empty()
                        {
                            if batch.signal_n_tr == 1 {
                                sig_n_uniq += 2;
                            } else {
                                sig_n_mult += 2;
                            }
                            for (tr, second_mate) in &batch.signal_contrib {
                                sig.add_transcript(
                                    &index.genome,
                                    tr,
                                    batch.signal_n_tr,
                                    *second_mate,
                                );
                            }
                        }
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
                    let read_file = tf.reopen().map_err(|e| {
                        anyhow::anyhow!("BySJout: temp file reopen for reading: {e}")
                    })?;
                    let mut reader =
                        noodles::sam::io::Reader::new(std::io::BufReader::new(read_file));
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
                                    let supp =
                                        build_within_bam_records(chim_aln, &index.genome, 255)?;
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

                info!("BySJout: filtered {filtered_count} pairs with non-surviving junctions");
            }

            // --outWigType bedGraph: write the four Signal.*.out.bg tracks.
            write_signal_tracks(signal.as_ref(), sig_n_uniq, sig_n_mult, params)?;

            // Flush chimeric output if enabled
            if let Some(ref mut chim_writer) = chimeric_writer {
                // --chimOutJunctionFormat 1: STAR-Fusion comment trailer with read counts
                // (# Nreads <total>\tNreadsUnique <uniquely_mapped>\tNreadsMulti <multi_mapped>).
                if params.chim_out_junction_format == 1 {
                    use std::sync::atomic::Ordering;
                    let command_line = params.command_line.as_deref().unwrap_or("");
                    let n_reads = stats.total_reads.load(Ordering::Relaxed);
                    let n_unique = stats.uniquely_mapped.load(Ordering::Relaxed);
                    let n_multi = stats.multi_mapped.load(Ordering::Relaxed);
                    chim_writer.write_format1_trailer(command_line, n_reads, n_unique, n_multi)?;
                }
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
        });

        // WASP allele-specific filtering: load the VCF once (shared read-only).
        // `None` unless --waspOutputMode SAMtag.
        let wasp_ctx: Arc<Option<crate::wasp::WaspContext>> = Arc::new(
            if params.wasp_output_mode == params::WaspOutputMode::SAMtag {
                let vcf = params
                    .var_vcf_file
                    .as_ref()
                    .expect("validated: SAMtag requires --varVCFfile");
                let ctx = crate::wasp::WaspContext::load(
                    vcf,
                    &index.genome.chr_name,
                    &index.genome.chr_start,
                    params,
                )
                .map_err(|source| error::Error::Io {
                    source,
                    path: vcf.clone(),
                })?;
                info!(
                    "WASP: loaded {} heterozygous SNVs from {}",
                    ctx.snps.len(),
                    vcf.display()
                );
                Some(ctx)
            } else {
                None
            },
        );

        // Stage 2: align each decoded pair batch on the rayon pool via
        // run_batch_pipeline, forwarding finished batches to the writer in order.
        let align_result = {
            let index = Arc::clone(index);
            let stats = Arc::clone(&stats);
            let sj_stats = Arc::clone(&sj_stats);
            let quant = quant.as_ref().map(Arc::clone);
            let tr = tr.as_ref().map(Arc::clone);
            let params_arc = Arc::clone(&params_arc);
            let wasp_ctx = Arc::clone(&wasp_ctx);
            let align = move |base: u64,
                              batch: Vec<crate::io::fastq::PairedRead>|
                  -> BatchOut<AlignmentBatchResults> {
                let params: &Parameters = &params_arc;
                // Adapter-aware clip params per mate (fixed Nbases are per-mate; the
                // adapter is shared). Applied via clip_mate below.
                let clip_params_m1 = crate::clip::clip_params_from(params, 0);
                let clip_params_m2 = crate::clip::clip_params_from(params, 1);
                batch
                    .par_iter()
                    .enumerate()
                    .map(|(pair_idx, paired_read)| {
                        #[allow(clippy::needless_borrow)]
                        let index = Arc::clone(&index);
                        #[allow(clippy::needless_borrow)]
                        let stats = Arc::clone(&stats);
                        #[allow(clippy::needless_borrow)]
                        let sj_stats = Arc::clone(&sj_stats);
                        let quant = quant.as_ref().map(Arc::clone);

                        // --outSAMreadID Number: numeric output QNAME (1-based input order,
                        // deterministic via `base`). Applies to the SAM/transcriptome QNAME
                        // and — for the unmapped-FASTX files — replaces both mate names with
                        // the shared pair index; the align_paired_read seed name is unchanged.
                        let out_read_id_number =
                            params.out_sam_read_id == crate::params::OutSamReadId::Number;
                        let out_read_name = if out_read_id_number {
                            (base + pair_idx as u64 + 1).to_string()
                        } else {
                            paired_read.name.clone()
                        };
                        let (fastx_name1, fastx_name2) = if out_read_id_number {
                            (out_read_name.clone(), out_read_name.clone())
                        } else {
                            (
                                paired_read.mate1.name.clone(),
                                paired_read.mate2.name.clone(),
                            )
                        };

                        // Apply clipping to both mates: each mate's adapter position is
                        // scanned independently (fixed Nbases + 3' adapter Hamming), then trim.
                        let (m1_clip5p, m1_clip3p) =
                            crate::clip::clip_mate(&paired_read.mate1.sequence, &clip_params_m1);
                        let (m2_clip5p, m2_clip3p) =
                            crate::clip::clip_mate(&paired_read.mate2.sequence, &clip_params_m2);
                        let (m1_seq, m1_qual) = clip_read(
                            &paired_read.mate1.sequence,
                            &paired_read.mate1.quality,
                            m1_clip5p,
                            m1_clip3p,
                        );
                        let (m2_seq, m2_qual) = clip_read(
                            &paired_read.mate2.sequence,
                            &paired_read.mate2.quality,
                            m2_clip5p,
                            m2_clip3p,
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
                                // Full original mates for unmapped pairs (STAR convention).
                                let records = SamWriter::build_paired_unmapped_records(
                                    &out_read_name,
                                    &paired_read.mate1.sequence,
                                    &paired_read.mate1.quality,
                                    &paired_read.mate2.sequence,
                                    &paired_read.mate2.quality,
                                    params,
                                    crate::stats::UnmappedReason::Other,
                                )?;
                                for record in records {
                                    buffer.push(record);
                                }
                            }
                            let (um1, um2) = if write_unmapped_fastq {
                                (
                                    vec![(
                                        fastx_name1.clone(),
                                        paired_read.mate1.sequence.clone(),
                                        paired_read.mate1.quality.clone(),
                                    )],
                                    vec![(
                                        fastx_name2.clone(),
                                        paired_read.mate2.sequence.clone(),
                                        paired_read.mate2.quality.clone(),
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
                                signal_contrib: Vec::new(),
                                signal_n_tr: 0,
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
                            stats.record_alignment(n_for_mapq, max_multimaps);
                            stats.record_unmapped_reason(
                                unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
                            );
                        } else if has_half_mapped {
                            // Half-mapped: count as mapped for the mapped mate
                            stats.record_alignment(1, max_multimaps);
                            stats.record_half_mapped();
                            // Record transcript stats from the mapped mate only
                            if let Some(PairedAlignmentResult::HalfMapped {
                                mapped_transcript,
                                ..
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

                        // Chimeric stats
                        if params.chim_segment_min > 0 && !pe_chimeric.is_empty() {
                            stats.record_chimeric();
                        }

                        // Gene-level quantification (lock-free atomic counts)
                        if let Some(ref q) = quant {
                            // Dereference Box<PairedAlignment> to get &PairedAlignment slice.
                            let bm_deref: Vec<&crate::align::read_align::PairedAlignment> =
                                both_mapped.iter().map(AsRef::as_ref).collect();
                            q.counts.count_pe_read(
                                &bm_deref,
                                results.is_empty(),
                                has_half_mapped,
                                &q.gene_ann,
                            );
                        }

                        // Record junction statistics
                        let is_unique =
                            both_mapped.len() == 1 || (has_half_mapped && results.len() == 1);
                        // Per-read junction dedup (fix A): a junction crossed by
                        // both mates / several loci counts once.
                        let mut read_trs: Vec<&crate::align::transcript::Transcript> = Vec::new();
                        for result in &results {
                            match result {
                                PairedAlignmentResult::BothMapped(pair) => {
                                    read_trs.push(&pair.mate1_transcript);
                                    read_trs.push(&pair.mate2_transcript);
                                }
                                PairedAlignmentResult::HalfMapped {
                                    mapped_transcript, ..
                                } => {
                                    read_trs.push(mapped_transcript);
                                }
                            }
                        }
                        record_read_junctions(read_trs, &index, &sj_stats, is_unique);

                        // --outWigType bedGraph: both-mapped pairs' reported alignments,
                        // both mates per pair (mate1 sense, mate2 flipped). signal_n_tr is
                        // the pair NH; each mate record is counted in the writer thread.
                        let signal_n_tr = both_mapped.len();
                        let signal_contrib: Vec<(crate::align::transcript::Transcript, bool)> =
                            if params.out_wig_bedgraph()
                                && !both_mapped.is_empty()
                                && both_mapped.len() <= max_multimaps
                            {
                                both_mapped
                                    .iter()
                                    .flat_map(|pa| {
                                        let pa = pa.as_ref();
                                        [
                                            (pa.mate1_transcript.clone(), false),
                                            (pa.mate2_transcript.clone(), true),
                                        ]
                                    })
                                    .collect()
                            } else {
                                Vec::new()
                            };

                        // Extract junction keys from primary alignment for BySJout
                        let primary_junction_keys = if by_sjout && !results.is_empty() {
                            let mut keys = Vec::new();
                            match &results[0] {
                                PairedAlignmentResult::BothMapped(pair) => {
                                    if pair.mate1_transcript.n_junction > 0 {
                                        keys.extend(extract_junction_keys(
                                            &pair.mate1_transcript,
                                            &index,
                                        ));
                                    }
                                    if pair.mate2_transcript.n_junction > 0 {
                                        keys.extend(extract_junction_keys(
                                            &pair.mate2_transcript,
                                            &index,
                                        ));
                                    }
                                }
                                PairedAlignmentResult::HalfMapped {
                                    mapped_transcript, ..
                                } => {
                                    if mapped_transcript.n_junction > 0 {
                                        keys.extend(extract_junction_keys(
                                            mapped_transcript,
                                            &index,
                                        ));
                                    }
                                }
                            }
                            keys
                        } else {
                            Vec::new()
                        };

                        // Build SAM records (skipped entirely under `--outSAMtype None`).
                        if !emit_sam {
                            // count/quant-only: no SAM record construction
                        } else if results.is_empty() {
                            // Unmapped pair
                            if output_unmapped {
                                // Full original mates for unmapped pairs (STAR convention).
                                let records = SamWriter::build_paired_unmapped_records(
                                    &out_read_name,
                                    &paired_read.mate1.sequence,
                                    &paired_read.mate1.quality,
                                    &paired_read.mate2.sequence,
                                    &paired_read.mate2.quality,
                                    params,
                                    unmapped_reason.unwrap_or(crate::stats::UnmappedReason::Other),
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
                                    &out_read_name,
                                    &paired_read.mate1.sequence,
                                    &paired_read.mate1.quality,
                                    &paired_read.mate2.sequence,
                                    &paired_read.mate2.quality,
                                    m1_clip5p,
                                    m1_clip3p,
                                    m2_clip5p,
                                    m2_clip3p,
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
                            let mut records = SamWriter::build_paired_records(
                                &out_read_name,
                                &paired_read.mate1.sequence,
                                &paired_read.mate1.quality,
                                &paired_read.mate2.sequence,
                                &paired_read.mate2.quality,
                                m1_clip5p,
                                m1_clip3p,
                                m2_clip5p,
                                m2_clip3p,
                                &paired_alns,
                                &index.genome,
                                params,
                                n_for_mapq,
                            )?;
                            // WASP allele-specific filtering: stamp vW/vA/vG by
                            // re-mapping the allele-swapped pair (STAR waspMap).
                            if let Some(ctx) = &*wasp_ctx {
                                crate::wasp::annotate_records_pe(
                                    &mut records,
                                    &paired_alns,
                                    &m1_seq,
                                    &m2_seq,
                                    &paired_read.name,
                                    &index,
                                    ctx,
                                    params.out_sam_attributes,
                                )?;
                            }
                            for record in records {
                                buffer.push(record);
                            }
                        }
                        // else: too many loci, skip output

                        // Transcriptome SAM projection (both-mapped pairs only)
                        let transcriptome_records: Vec<
                            noodles::sam::alignment::record_buf::RecordBuf,
                        > = if let Some(ref tidx) = tr_local {
                            build_transcriptome_records_pe(
                                both_mapped.iter().map(AsRef::as_ref),
                                &out_read_name,
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
                                    vec![(fastx_name1.clone(), m1_seq.clone(), m1_qual.clone())],
                                    vec![(fastx_name2.clone(), m2_seq.clone(), m2_qual.clone())],
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
                            signal_contrib,
                            signal_n_tr,
                        })
                    })
                    .collect()
            };
            run_batch_pipeline(
                read_rx,
                max_reads,
                PIPELINE_DEPTH,
                |n: u64| {
                    if n % 100_000 < batch_size as u64 {
                        info!("Processed {n} pairs...");
                    }
                },
                align,
                |done| Ok(res_tx.send(done).is_ok()),
            )
        };
        // Disconnect the writer channel so the writer thread can finish and join.
        drop(res_tx);
        let writer_result = writer_handle
            .join()
            .map_err(|_| anyhow::anyhow!("PE writer thread panicked"))?;
        align_result?;
        writer_result?;
        Ok(())
    })?;

    Ok(())
}

/// Record junctions from a transcript into SJ statistics
/// One splice junction extracted from a transcript's CIGAR (before per-read
/// deduplication / recording).
struct ReadJunction {
    chr_idx: usize,
    intron_start: u64,
    intron_end: u64,
    strand: u8,
    motif: crate::align::score::SpliceMotif,
    overhang: u32,
    annotated: bool,
}

/// Record a read's splice junctions into `sj_stats`, **deduplicated per read**:
/// a junction crossed by several of the read's alignments/mates counts ONCE,
/// taking the max overhang across occurrences (STAR `outputTranscriptSJ`;
/// confirmed against the byte-faithful STAR-rs `star_sj.rs::record_read`).
/// `is_unique` reflects the read's overall mapping multiplicity (n_loci == 1),
/// not per-locus. Recording per-locus/per-mate instead double-counts junctions
/// crossed by both mates of a pair, inflating the SJ.out.tab multi counts.
fn record_read_junctions<'a>(
    transcripts: impl IntoIterator<Item = &'a crate::align::transcript::Transcript>,
    index: &crate::index::GenomeIndex,
    sj_stats: &crate::junction::SpliceJunctionStats,
    is_unique: bool,
) {
    use std::collections::HashMap;
    let mut per_read: HashMap<(usize, u64, u64), ReadJunction> = HashMap::new();
    for t in transcripts {
        for j in extract_transcript_junctions(t, index) {
            per_read
                .entry((j.chr_idx, j.intron_start, j.intron_end))
                .and_modify(|e| {
                    e.overhang = e.overhang.max(j.overhang);
                    e.annotated |= j.annotated;
                })
                .or_insert(j);
        }
    }
    for j in per_read.values() {
        sj_stats.record_junction(
            j.chr_idx,
            j.intron_start,
            j.intron_end,
            j.strand,
            j.motif,
            is_unique,
            j.overhang,
            j.annotated,
        );
    }
}

/// Extract all splice junctions (`N` / Skip ops) from one transcript's CIGAR,
/// with per-junction overhang, motif, strand and annotation — without recording
/// them (recording is done per-read by `record_read_junctions`).
fn extract_transcript_junctions(
    transcript: &crate::align::transcript::Transcript,
    index: &crate::index::GenomeIndex,
) -> Vec<ReadJunction> {
    use crate::align::score::AlignmentScorer;
    use cigar::op::Kind;

    let mut out: Vec<ReadJunction> = Vec::new();

    // First pass: compute exon segment lengths (query-consuming bases between N operations)
    // An "exon segment" is the query bases on each side of a splice junction.
    let mut exon_lengths: Vec<u32> = Vec::new();
    let mut current_exon_len = 0u32;

    for op in &transcript.cigar {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion => {
                current_exon_len += op.len() as u32;
            }
            Kind::Skip => {
                exon_lengths.push(current_exon_len);
                current_exon_len = 0;
            }
            // Soft clips, deletions, hard clips do not contribute to overhang
            // STAR counts only matched/inserted bases (not soft-clipped bases)
            Kind::Deletion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    exon_lengths.push(current_exon_len); // Final exon segment

    // Second pass: record junctions with computed overhangs
    let mut genome_pos = transcript.genome_start;
    let mut junction_idx = 0usize;

    let scorer = AlignmentScorer::from_params_minimal();

    for op in &transcript.cigar {
        match op.kind() {
            Kind::Skip => {
                // This is a splice junction
                let intron_len = op.len();
                let intron_start = genome_pos;
                let intron_end = genome_pos + intron_len as u64 - 1;

                // Detect splice motif
                let motif =
                    scorer.detect_splice_motif(genome_pos, intron_len as u32, &index.genome);

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

                out.push(ReadJunction {
                    chr_idx: transcript.chr_idx,
                    intron_start,
                    intron_end,
                    strand,
                    motif,
                    overhang,
                    annotated,
                });

                // Advance genome position past the intron
                genome_pos += intron_len as u64;
                junction_idx += 1;
            }
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                genome_pos += op.len() as u64;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }

    out
}
