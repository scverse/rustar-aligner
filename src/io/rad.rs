//! RAD output of transcriptome-space alignments (`--quantTranscriptomeRAD Yes`).
//!
//! RAD ("Reduced Alignment Data") is the binary format of the COMBINE-lab tools.
//! A fragment is stored as the list of its placements (transcript, strands,
//! position, fragment length, score), with no read name, sequence, qualities or
//! CIGAR, which is all `salmon quant --rad` reads.
//!
//! The layout written here is the `salmon_bulk_sa` profile of salmon 2.x
//! (`crates/salmon-rad`; tested with 2.8.0), which salmon itself writes from a
//! transcriptome BAM (BSD 3-clause, see `LICENSES/SALMON_LICENSE`). All integers
//! are little-endian.
//!
//! ```text
//! prelude   is_paired: u8, ref_count: u64, ref_names: (len: u16, bytes)*,
//!           num_chunks: u64
//!           file tags, read tags, alignment tags, each
//!             count: u16, (name_len: u16, name, type_id: u8,
//!                          [len_type: u8, value_type: u8 for arrays])*
//! values    the file tag values, in the order they were declared
//! chunks    nbytes: u32 (including this 8-byte header), nrec: u32, records
//! record    na: u32, frag_map_type: u8,
//!           (compressed_ori_ref: u32, pos: u32, frag_len: u16,
//!            alignment_score: i32) * na
//! ```

use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::align::read_align::calculate_insert_size;
use crate::align::transcript::Transcript;
use crate::error::Error;

/// `rad_type` value that identifies the profile to salmon. Reads as "salmon's
/// format", not "written by salmon": `mapping_type` and `source_programs` say
/// where the alignments came from.
const RAD_TYPE: &str = "salmon_bulk_sa";
/// `mapping_type` value for fragments from an external aligner.
const MAPPING_TYPE_ALIGNMENT: &str = "alignment";

const TYPE_U8: u8 = 1;
const TYPE_U16: u8 = 2;
const TYPE_U32: u8 = 3;
const TYPE_U64: u8 = 4;
const TYPE_I32: u8 = 12;
const TYPE_ARRAY: u8 = 7;
const TYPE_STRING: u8 = 8;

/// `baked_flags` bit: `num_processed` holds a real count.
const BAKED_MAP_COUNTERS: u8 = 0x20;
/// `baked_flags` bit: the writer reached the end of the run.
const WRITE_COMPLETE: u8 = 0x10;

/// `frag_map_type` of a single-end fragment.
pub const FRAG_SINGLE: u8 = 1;
/// `frag_map_type` of a fragment mapped as a proper pair.
pub const FRAG_MAPPED_PAIR: u8 = 4;

/// `compressed_ori_ref` bit: read 1 is on the forward strand of the transcript.
const ORI_FW: u32 = 0x8000_0000;
/// `compressed_ori_ref` bit: read 2 is on the forward strand of the transcript.
const MATE_FW: u32 = 0x4000_0000;
/// `compressed_ori_ref` mask of the transcript index.
const TID_MASK: u32 = 0x3FFF_FFFF;
/// A `frag_len` of this value reads as "no fragment length".
const FRAG_LEN_UNPAIRED: u16 = u16::MAX;

/// Chunks are flushed at this many record bytes (salmon's own value).
const CHUNK_FLUSH_BYTES: usize = 64 * 1024;

/// One placement of a fragment on a transcript.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RadHit {
    pub tid: u32,
    pub is_fw: bool,
    pub mate_fw: bool,
    /// 0-based leftmost transcript position of the fragment.
    pub pos: u32,
    /// Fragment length of a pair; aligned length on the transcript of a
    /// single-end read, which salmon uses to bound its fragment length.
    pub frag_len: u16,
    /// STAR's alignment score of the whole read or pair (what `AS` holds in the
    /// genomic BAM). `Aligned.toTranscriptome.out.bam` has no `AS` tag, so this is
    /// the one thing here the BAM does not carry.
    pub score: i32,
}

impl RadHit {
    fn compressed_ori_ref(&self) -> u32 {
        (self.tid & TID_MASK)
            | if self.is_fw { ORI_FW } else { 0 }
            | if self.mate_fw { MATE_FW } else { 0 }
    }
}

/// One fragment: all of its transcript-space placements.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RadFragment {
    pub frag_map_type: u8,
    pub hits: Vec<RadHit>,
}

impl RadFragment {
    /// A single-end read's projections, as returned by
    /// [`filter_and_project`](crate::quant::transcriptome::filter_and_project).
    pub fn single_end(projected: &[Transcript]) -> Self {
        let hits = projected
            .iter()
            .map(|t| RadHit {
                tid: t.chr_idx as u32,
                is_fw: !t.is_reverse,
                mate_fw: false,
                pos: t.genome_start as u32,
                // MAX - 1: a length must not be mistaken for the "unpaired" sentinel.
                frag_len: (t.genome_end - t.genome_start).min(u64::from(FRAG_LEN_UNPAIRED - 1))
                    as u16,
                score: t.score,
            })
            .collect();
        Self {
            frag_map_type: FRAG_SINGLE,
            hits,
        }
    }

    /// A pair's projections, as returned by
    /// [`filter_and_project_pair`](crate::quant::transcriptome::filter_and_project_pair),
    /// with the score of the genomic alignment each came from (one per pair).
    pub fn paired(projected: &[(Transcript, Transcript)], scores: &[i32]) -> Self {
        debug_assert_eq!(projected.len(), scores.len());
        let hits = projected
            .iter()
            .zip(scores)
            .map(|((m1, m2), &score)| RadHit {
                tid: m1.chr_idx as u32,
                is_fw: !m1.is_reverse,
                mate_fw: !m2.is_reverse,
                pos: m1.genome_start.min(m2.genome_start) as u32,
                // Saturates at the sentinel: an absurd length must not become a
                // small one.
                frag_len: calculate_insert_size(m1, m2)
                    .unsigned_abs()
                    .min(u32::from(FRAG_LEN_UNPAIRED)) as u16,
                score,
            })
            .collect();
        Self {
            frag_map_type: FRAG_MAPPED_PAIR,
            hits,
        }
    }

    fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&(self.hits.len() as u32).to_le_bytes());
        out.push(self.frag_map_type);
        for h in &self.hits {
            out.extend_from_slice(&h.compressed_ori_ref().to_le_bytes());
            out.extend_from_slice(&h.pos.to_le_bytes());
            out.extend_from_slice(&h.frag_len.to_le_bytes());
            out.extend_from_slice(&h.score.to_le_bytes());
        }
    }
}

/// Sequential RAD writer. The file appears under its final name only once
/// [`finish`](Self::finish) has run, so an interrupted run cannot leave a truncated
/// file that looks complete.
pub struct RadWriter {
    file: BufWriter<File>,
    final_path: PathBuf,
    partial_path: PathBuf,
    num_chunks_offset: u64,
    baked_flags_offset: u64,
    num_processed_offset: u64,
    chunk: Vec<u8>,
    chunk_records: u32,
    num_chunks: u64,
}

impl RadWriter {
    /// Create the file and write its prelude. `ref_names` and `ref_lengths` define
    /// the transcript indices that hits refer to. `source_program` is a SAM `@PG`
    /// line, which salmon copies to `meta_info.json`.
    pub fn create(
        path: &Path,
        ref_names: &[String],
        ref_lengths: &[u32],
        is_paired: bool,
        source_program: &str,
    ) -> Result<Self, Error> {
        let mut partial_name = path.file_name().unwrap_or_default().to_os_string();
        partial_name.push(".partial");
        let partial_path = path.with_file_name(partial_name);

        let prelude = Prelude::encode(ref_names, ref_lengths, is_paired, source_program)?;
        let file = File::create(&partial_path).map_err(|e| Error::io(e, &partial_path))?;
        let mut file = BufWriter::with_capacity(1 << 20, file);
        file.write_all(&prelude.bytes)
            .map_err(|e| Error::io(e, &partial_path))?;
        Ok(Self {
            file,
            final_path: path.to_path_buf(),
            partial_path,
            num_chunks_offset: prelude.num_chunks_offset,
            baked_flags_offset: prelude.baked_flags_offset,
            num_processed_offset: prelude.num_processed_offset,
            chunk: Vec::with_capacity(CHUNK_FLUSH_BYTES + 4096),
            chunk_records: 0,
            num_chunks: 0,
        })
    }

    /// Append one fragment. Fragments without a placement are not written:
    /// salmon learns their number from `num_processed`.
    pub fn write(&mut self, fragment: &RadFragment) -> Result<(), Error> {
        if fragment.hits.is_empty() {
            return Ok(());
        }
        fragment.encode(&mut self.chunk);
        self.chunk_records += 1;
        if self.chunk.len() >= CHUNK_FLUSH_BYTES {
            self.flush_chunk()?;
        }
        Ok(())
    }

    fn flush_chunk(&mut self) -> Result<(), Error> {
        if self.chunk_records == 0 {
            return Ok(());
        }
        let nbytes = self.chunk.len() as u32 + 8;
        let result = self
            .file
            .write_all(&nbytes.to_le_bytes())
            .and_then(|()| self.file.write_all(&self.chunk_records.to_le_bytes()))
            .and_then(|()| self.file.write_all(&self.chunk));
        result.map_err(|e| Error::io(e, &self.partial_path))?;
        self.chunk.clear();
        self.chunk_records = 0;
        self.num_chunks += 1;
        Ok(())
    }

    /// Write the last chunk, record the number of input fragments (mapped or not)
    /// in the header, and give the file its final name.
    pub fn finish(mut self, num_processed: u64) -> Result<(), Error> {
        self.flush_chunk()?;
        let patches: [(u64, &[u8]); 3] = [
            (self.num_chunks_offset, &self.num_chunks.to_le_bytes()),
            (
                self.baked_flags_offset,
                &[BAKED_MAP_COUNTERS | WRITE_COMPLETE],
            ),
            (self.num_processed_offset, &num_processed.to_le_bytes()),
        ];
        for (offset, bytes) in patches {
            self.file
                .seek(SeekFrom::Start(offset))
                .and_then(|_| self.file.write_all(bytes))
                .map_err(|e| Error::io(e, &self.partial_path))?;
        }
        self.file
            .flush()
            .map_err(|e| Error::io(e, &self.partial_path))?;
        drop(self.file);
        std::fs::rename(&self.partial_path, &self.final_path)
            .map_err(|e| Error::io(e, &self.final_path))
    }
}

/// Prelude and file tag values, with the offsets [`RadWriter::finish`] patches.
struct Prelude {
    bytes: Vec<u8>,
    num_chunks_offset: u64,
    baked_flags_offset: u64,
    num_processed_offset: u64,
}

impl Prelude {
    fn encode(
        ref_names: &[String],
        ref_lengths: &[u32],
        is_paired: bool,
        source_program: &str,
    ) -> Result<Self, Error> {
        if ref_names.len() != ref_lengths.len() {
            return Err(Error::Alignment(format!(
                "RAD output: {} transcript names but {} lengths",
                ref_names.len(),
                ref_lengths.len()
            )));
        }
        if ref_names.len() > TID_MASK as usize {
            return Err(Error::Alignment(format!(
                "RAD output: {} transcripts do not fit in 30 bits",
                ref_names.len()
            )));
        }
        let mut b = Vec::with_capacity(ref_names.iter().map(|n| n.len() + 2).sum::<usize>() + 4096);
        b.push(u8::from(is_paired));
        b.extend_from_slice(&(ref_names.len() as u64).to_le_bytes());
        for name in ref_names {
            put_str(&mut b, name)?;
        }
        let num_chunks_offset = b.len() as u64;
        b.extend_from_slice(&0u64.to_le_bytes());

        // File tags, read tags, alignment tags. The `Array(u32, u32)` of
        // `ref_lengths` is what salmon requires; everything else is the minimum
        // that a RAD without a baked fragment-length distribution needs.
        put_tag_section(
            &mut b,
            &[
                ("rad_type", &[TYPE_STRING]),
                ("ref_lengths", &[TYPE_ARRAY, TYPE_U32, TYPE_U32]),
                ("baked_flags", &[TYPE_U8]),
                ("num_processed", &[TYPE_U64]),
                ("mapping_type", &[TYPE_STRING]),
                ("source_programs", &[TYPE_STRING]),
            ],
        )?;
        put_tag_section(&mut b, &[("frag_map_type", &[TYPE_U8])])?;
        put_tag_section(
            &mut b,
            &[
                ("compressed_ori_ref", &[TYPE_U32]),
                ("pos", &[TYPE_U32]),
                ("frag_len", &[TYPE_U16]),
                ("alignment_score", &[TYPE_I32]),
            ],
        )?;

        // File tag values, in declaration order.
        put_str(&mut b, RAD_TYPE)?;
        b.extend_from_slice(&(ref_lengths.len() as u32).to_le_bytes());
        for len in ref_lengths {
            b.extend_from_slice(&len.to_le_bytes());
        }
        let baked_flags_offset = b.len() as u64;
        b.push(0);
        let num_processed_offset = b.len() as u64;
        b.extend_from_slice(&0u64.to_le_bytes());
        put_str(&mut b, MAPPING_TYPE_ALIGNMENT)?;
        // Provenance, not worth failing a run over: past a string tag's 64 KiB
        // limit, record nothing rather than a cut command line.
        put_str(
            &mut b,
            if u16::try_from(source_program.len()).is_ok() {
                source_program
            } else {
                ""
            },
        )?;
        Ok(Self {
            bytes: b,
            num_chunks_offset,
            baked_flags_offset,
            num_processed_offset,
        })
    }
}

fn put_str(b: &mut Vec<u8>, s: &str) -> Result<(), Error> {
    let len = u16::try_from(s.len()).map_err(|_| {
        Error::Alignment(format!(
            "RAD output: string of {} bytes does not fit in a u16 length",
            s.len()
        ))
    })?;
    b.extend_from_slice(&len.to_le_bytes());
    b.extend_from_slice(s.as_bytes());
    Ok(())
}

fn put_tag_section(b: &mut Vec<u8>, tags: &[(&str, &[u8])]) -> Result<(), Error> {
    b.extend_from_slice(&(tags.len() as u16).to_le_bytes());
    for (name, type_ids) in tags {
        put_str(b, name)?;
        b.extend_from_slice(type_ids);
    }
    Ok(())
}

/// A file tag value, as [`read_rad`] returns it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RadTagValue {
    U8(u8),
    U64(u64),
    ArrayU32(Vec<u32>),
    String(String),
}

/// The contents of a file written by [`RadWriter`].
#[derive(Debug)]
pub struct RadFile {
    pub is_paired: bool,
    pub ref_names: Vec<String>,
    /// `num_chunks` in the header.
    pub num_chunks: u64,
    pub file_tags: Vec<(String, RadTagValue)>,
    /// Names of the read tags, in order.
    pub read_tags: Vec<String>,
    /// Names of the alignment tags, in order.
    pub aln_tags: Vec<String>,
    /// Number of records of each chunk, as its header states.
    pub chunk_records: Vec<u32>,
    pub fragments: Vec<RadFragment>,
}

impl RadFile {
    /// The value of file tag `name`.
    pub fn file_tag(&self, name: &str) -> Option<&RadTagValue> {
        self.file_tags
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v)
    }
}

/// Read back a file written by [`RadWriter`]. Other RAD profiles are refused
/// rather than misread; quantifying is salmon's job.
pub fn read_rad(path: &Path) -> Result<RadFile, Error> {
    let bytes = std::fs::read(path).map_err(|e| Error::io(e, path))?;
    parse_rad(&bytes).ok_or_else(|| {
        Error::Alignment(format!(
            "{} is truncated or not a RAD file written by rustar-aligner",
            path.display()
        ))
    })
}

fn parse_rad(bytes: &[u8]) -> Option<RadFile> {
    let mut c = Cursor { bytes, pos: 0 };
    let is_paired = c.u8()? != 0;
    let ref_count = c.u64()?;
    let mut ref_names = Vec::new();
    for _ in 0..ref_count {
        ref_names.push(c.string()?);
    }
    let num_chunks = c.u64()?;
    let file_desc = c.tag_section()?;
    let read_tags = c.tag_section()?;
    let aln_tags = c.tag_section()?;
    let mut file_tags = Vec::with_capacity(file_desc.len());
    for (name, type_id) in file_desc {
        let value = match type_id {
            TYPE_U8 => RadTagValue::U8(c.u8()?),
            TYPE_U64 => RadTagValue::U64(c.u64()?),
            TYPE_ARRAY => {
                let n = c.u32()?;
                let mut v = Vec::new();
                for _ in 0..n {
                    v.push(c.u32()?);
                }
                RadTagValue::ArrayU32(v)
            }
            TYPE_STRING => RadTagValue::String(c.string()?),
            _ => return None,
        };
        file_tags.push((name, value));
    }

    let mut chunk_records = Vec::new();
    let mut fragments = Vec::new();
    while c.pos < bytes.len() {
        let chunk_start = c.pos;
        let nbytes = c.u32()? as usize;
        let nrec = c.u32()?;
        chunk_records.push(nrec);
        for _ in 0..nrec {
            let na = c.u32()?;
            let frag_map_type = c.u8()?;
            let mut hits = Vec::new();
            for _ in 0..na {
                let ori = c.u32()?;
                let pos = c.u32()?;
                let frag_len = c.u16()?;
                let score = c.u32()? as i32;
                hits.push(RadHit {
                    tid: ori & TID_MASK,
                    is_fw: ori & ORI_FW != 0,
                    mate_fw: ori & MATE_FW != 0,
                    pos,
                    frag_len,
                    score,
                });
            }
            fragments.push(RadFragment {
                frag_map_type,
                hits,
            });
        }
        if c.pos - chunk_start != nbytes {
            return None;
        }
    }
    Some(RadFile {
        is_paired,
        ref_names,
        num_chunks,
        file_tags,
        read_tags: read_tags.into_iter().map(|(n, _)| n).collect(),
        aln_tags: aln_tags.into_iter().map(|(n, _)| n).collect(),
        chunk_records,
        fragments,
    })
}

struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Cursor<'_> {
    fn take(&mut self, n: usize) -> Option<&[u8]> {
        let s = self.bytes.get(self.pos..self.pos.checked_add(n)?)?;
        self.pos += n;
        Some(s)
    }

    fn u8(&mut self) -> Option<u8> {
        Some(self.take(1)?[0])
    }

    fn u16(&mut self) -> Option<u16> {
        Some(u16::from_le_bytes(self.take(2)?.try_into().ok()?))
    }

    fn u32(&mut self) -> Option<u32> {
        Some(u32::from_le_bytes(self.take(4)?.try_into().ok()?))
    }

    fn u64(&mut self) -> Option<u64> {
        Some(u64::from_le_bytes(self.take(8)?.try_into().ok()?))
    }

    fn string(&mut self) -> Option<String> {
        let n = usize::from(self.u16()?);
        String::from_utf8(self.take(n)?.to_vec()).ok()
    }

    /// Name and type id of each tag of a section.
    fn tag_section(&mut self) -> Option<Vec<(String, u8)>> {
        let n = self.u16()?;
        let mut tags = Vec::with_capacity(usize::from(n));
        for _ in 0..n {
            let name = self.string()?;
            let type_id = self.u8()?;
            if type_id == TYPE_ARRAY && (self.u8()?, self.u8()?) != (TYPE_U32, TYPE_U32) {
                return None;
            }
            tags.push((name, type_id));
        }
        Some(tags)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::Exon;
    use tempfile::TempDir;

    fn projected(tid: usize, start: u64, end: u64, is_reverse: bool, score: i32) -> Transcript {
        Transcript {
            output_rank: 0,
            chr_idx: tid,
            genome_start: start,
            genome_end: end,
            is_reverse,
            exons: vec![Exon {
                genome_start: start,
                genome_end: end,
                read_start: 0,
                read_end: (end - start) as usize,
                i_frag: 0,
            }],
            cigar: Vec::new(),
            score,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: Vec::new(),
            junction_annotated: Vec::new(),
            junction_strand: Vec::new(),
        }
    }

    fn hit(tid: u32, is_fw: bool, mate_fw: bool, pos: u32, frag_len: u16, score: i32) -> RadHit {
        RadHit {
            tid,
            is_fw,
            mate_fw,
            pos,
            frag_len,
            score,
        }
    }

    fn names(n: &[&str]) -> Vec<String> {
        n.iter().map(ToString::to_string).collect()
    }

    fn create(dir: &TempDir, is_paired: bool) -> (RadWriter, PathBuf) {
        let path = dir.path().join("t.rad");
        let writer = RadWriter::create(
            &path,
            &names(&["ENST1", "ENST2"]),
            &[1000, 2000],
            is_paired,
            "@PG\tID:rustar-aligner",
        )
        .unwrap();
        (writer, path)
    }

    #[test]
    fn single_end_hit_carries_the_aligned_length() {
        let f = RadFragment::single_end(&[projected(3, 10, 110, true, 98)]);
        assert_eq!(f.frag_map_type, FRAG_SINGLE);
        assert_eq!(f.hits, vec![hit(3, false, false, 10, 100, 98)]);
    }

    #[test]
    fn single_end_length_stays_below_the_unpaired_sentinel() {
        let f = RadFragment::single_end(&[projected(0, 0, 1 << 20, false, 1)]);
        assert_eq!(f.hits[0].frag_len, u16::MAX - 1);
    }

    #[test]
    fn pair_on_forward_strand_uses_read1_start_and_template_length() {
        let pair = (
            projected(7, 100, 200, false, 0),
            projected(7, 250, 350, true, 0),
        );
        let f = RadFragment::paired(&[pair], &[198]);
        assert_eq!(f.frag_map_type, FRAG_MAPPED_PAIR);
        assert_eq!(f.hits, vec![hit(7, true, false, 100, 250, 198)]);
    }

    #[test]
    fn pair_on_reverse_strand_starts_at_mate2() {
        let pair = (
            projected(7, 250, 350, true, 0),
            projected(7, 100, 200, false, 0),
        );
        let f = RadFragment::paired(&[pair], &[198]);
        assert_eq!(f.hits, vec![hit(7, false, true, 100, 250, 198)]);
    }

    #[test]
    fn absurd_pair_length_saturates_at_the_sentinel() {
        let pair = (
            projected(0, 0, 100, false, 0),
            projected(0, 1 << 20, (1 << 20) + 100, true, 0),
        );
        let f = RadFragment::paired(&[pair], &[0]);
        assert_eq!(f.hits[0].frag_len, FRAG_LEN_UNPAIRED);
    }

    #[test]
    fn orientation_bits_match_the_piscem_layout() {
        assert_eq!(
            hit(5, true, false, 0, 0, 0).compressed_ori_ref(),
            0x8000_0005
        );
        assert_eq!(
            hit(5, false, true, 0, 0, 0).compressed_ori_ref(),
            0x4000_0005
        );
        assert_eq!(hit(5, false, false, 0, 0, 0).compressed_ori_ref(), 5);
    }

    #[test]
    fn prelude_round_trips() {
        let dir = TempDir::new().unwrap();
        let (writer, path) = create(&dir, true);
        writer.finish(42).unwrap();

        let rad = read_rad(&path).unwrap();
        assert!(rad.is_paired);
        assert_eq!(rad.ref_names, names(&["ENST1", "ENST2"]));
        assert_eq!(rad.num_chunks, 0);
        assert_eq!(rad.read_tags, names(&["frag_map_type"]));
        assert_eq!(
            rad.aln_tags,
            names(&["compressed_ori_ref", "pos", "frag_len", "alignment_score"])
        );
        assert_eq!(
            rad.file_tag("rad_type"),
            Some(&RadTagValue::String("salmon_bulk_sa".into()))
        );
        assert_eq!(
            rad.file_tag("ref_lengths"),
            Some(&RadTagValue::ArrayU32(vec![1000, 2000]))
        );
        assert_eq!(
            rad.file_tag("mapping_type"),
            Some(&RadTagValue::String("alignment".into()))
        );
        assert_eq!(
            rad.file_tag("source_programs"),
            Some(&RadTagValue::String("@PG\tID:rustar-aligner".into()))
        );
        assert_eq!(rad.file_tag("num_processed"), Some(&RadTagValue::U64(42)));
        assert_eq!(rad.file_tag("baked_flags"), Some(&RadTagValue::U8(0x30)));
    }

    #[test]
    fn prelude_bytes_match_the_documented_layout() {
        // Fixed bytes: a reader in another language sees exactly this.
        let p = Prelude::encode(&names(&["T"]), &[7], false, "").unwrap();
        let mut expected: Vec<u8> = vec![0]; // is_paired
        expected.extend_from_slice(&1u64.to_le_bytes()); // ref_count
        expected.extend_from_slice(&[1, 0, b'T']);
        assert_eq!(&p.bytes[..expected.len()], &expected[..]);
        assert_eq!(p.num_chunks_offset, expected.len() as u64);
        assert_eq!(
            &p.bytes[expected.len() + 8..expected.len() + 10],
            &6u16.to_le_bytes(),
            "six file tags"
        );
    }

    #[test]
    fn records_round_trip_across_chunks() {
        let dir = TempDir::new().unwrap();
        let (mut writer, path) = create(&dir, true);
        // Enough to overflow one chunk: 5 + 14 * 2 bytes a record.
        let n = CHUNK_FLUSH_BYTES / 33 * 2 + 17;
        let fragments: Vec<RadFragment> = (0..n)
            .map(|i| RadFragment {
                frag_map_type: FRAG_MAPPED_PAIR,
                hits: vec![
                    hit(1, true, false, i as u32, 300, 198),
                    hit(0, false, true, 3, u16::MAX, -7),
                ],
            })
            .collect();
        for f in &fragments {
            writer.write(f).unwrap();
        }
        writer.finish(n as u64 + 5).unwrap();

        let rad = read_rad(&path).unwrap();
        assert_eq!(rad.fragments, fragments);
        assert_eq!(rad.num_chunks, rad.chunk_records.len() as u64);
        assert!(rad.num_chunks >= 2, "got {} chunks", rad.num_chunks);
        assert_eq!(rad.chunk_records.iter().sum::<u32>() as usize, n);
    }

    #[test]
    fn record_bytes_match_the_documented_layout() {
        let mut b = Vec::new();
        RadFragment {
            frag_map_type: FRAG_MAPPED_PAIR,
            hits: vec![hit(2, true, false, 16, 300, -1)],
        }
        .encode(&mut b);
        let mut expected = 1u32.to_le_bytes().to_vec();
        expected.push(4);
        expected.extend_from_slice(&0x8000_0002u32.to_le_bytes());
        expected.extend_from_slice(&16u32.to_le_bytes());
        expected.extend_from_slice(&300u16.to_le_bytes());
        expected.extend_from_slice(&(-1i32).to_le_bytes());
        assert_eq!(b, expected);
    }

    #[test]
    fn fragments_without_placements_are_not_written() {
        let dir = TempDir::new().unwrap();
        let (mut writer, path) = create(&dir, false);
        writer
            .write(&RadFragment {
                frag_map_type: FRAG_SINGLE,
                hits: Vec::new(),
            })
            .unwrap();
        writer.finish(1).unwrap();
        let rad = read_rad(&path).unwrap();
        assert_eq!(rad.fragments.len(), 0);
        assert_eq!(rad.num_chunks, 0);
    }

    #[test]
    fn file_appears_only_when_finished() {
        let dir = TempDir::new().unwrap();
        let (writer, path) = create(&dir, false);
        let partial = dir.path().join("t.rad.partial");
        assert!(partial.exists());
        assert!(!path.exists());
        writer.finish(0).unwrap();
        assert!(path.exists());
        assert!(!partial.exists());
    }

    #[test]
    fn oversized_command_line_is_dropped_whole() {
        let p = Prelude::encode(&names(&["T"]), &[7], false, &"x".repeat(70_000)).unwrap();
        assert!(p.bytes.len() < 1000);
    }

    #[test]
    fn mismatched_reference_table_is_refused() {
        assert!(Prelude::encode(&names(&["T"]), &[], false, "").is_err());
    }

    #[test]
    fn truncated_file_is_refused() {
        let dir = TempDir::new().unwrap();
        let (mut writer, path) = create(&dir, false);
        writer
            .write(&RadFragment::single_end(&[projected(0, 0, 10, false, 1)]))
            .unwrap();
        writer.finish(1).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        std::fs::write(&path, &bytes[..bytes.len() - 1]).unwrap();
        assert!(read_rad(&path).is_err());
    }
}
