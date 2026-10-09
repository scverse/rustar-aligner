//! Input decompression with format detection by magic bytes (#218).
//!
//! A FASTQ input is opened, its first bytes are peeked, and the matching
//! decoder is stacked on top. The file name plays no part: a gzip file called
//! `reads.fq` is decoded, and a plain file called `reads.fq.gz` is read as
//! plain text. Detection reads from any `Read`, so pipes and `/dev/stdin` are
//! handled the same way as regular files.
//!
//! gzip is always available (`flate2`, multi-member aware). bzip2, Zstandard
//! and XZ are behind the optional Cargo features `bz2`, `zstd` and `xz`
//! (`compressed-input` enables all three); all three backends are pure Rust.
//! A build without the feature for a detected format fails with an error
//! naming the feature and the `--readFilesCommand` alternative, instead of
//! handing compressed bytes to the FASTQ parser.
//!
//! STAR itself does no detection and relies on `--readFilesCommand`; that
//! option still takes precedence over everything here (see `DIVERGENCE.md`).

use flate2::read::MultiGzDecoder;
use std::fmt;
use std::io::{self, BufRead, BufReader, Cursor, Read};

/// Buffer size for both the compressed side and the decoded side. Larger than
/// the 8 KiB default to cut read syscalls on the decode hot path.
pub const DECODE_BUF: usize = 1 << 19; // 512 KiB

/// Number of leading bytes needed to tell the supported formats apart.
const MAGIC_LEN: usize = 6;

/// Compression format of an input stream, as detected from its first bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Not a recognised compressed format; read as-is.
    Plain,
    /// gzip, including multi-member files and BGZF (`1f 8b`).
    Gzip,
    /// bzip2 (`42 5a 68`, "BZh").
    Bzip2,
    /// Zstandard, a regular frame (`28 b5 2f fd`) or a skippable frame
    /// (`5? 2a 4d 18`, which `pzstd` writes first).
    Zstd,
    /// XZ (`fd 37 7a 58 5a 00`).
    Xz,
}

impl Format {
    /// Identify the format from the first bytes of a stream. Fewer than
    /// [`MAGIC_LEN`] bytes are fine: a short prefix simply matches fewer
    /// signatures.
    #[must_use]
    pub fn detect(prefix: &[u8]) -> Self {
        match prefix {
            [0x1f, 0x8b, ..] => Format::Gzip,
            [b'B', b'Z', b'h', ..] => Format::Bzip2,
            [0x28, 0xb5, 0x2f, 0xfd, ..] => Format::Zstd,
            [m, 0x2a, 0x4d, 0x18, ..] if m & 0xf0 == 0x50 => Format::Zstd,
            [0xfd, b'7', b'z', b'X', b'Z', 0x00, ..] => Format::Xz,
            _ => Format::Plain,
        }
    }

    /// The Cargo feature that enables decoding this format, if any.
    fn feature(self) -> Option<&'static str> {
        match self {
            Format::Plain | Format::Gzip => None,
            Format::Bzip2 => Some("bz2"),
            Format::Zstd => Some("zstd"),
            Format::Xz => Some("xz"),
        }
    }

    /// A decompressor that `--readFilesCommand` can use instead.
    fn command_hint(self) -> &'static str {
        match self {
            Format::Plain | Format::Gzip => "zcat",
            Format::Bzip2 => "bzcat",
            Format::Zstd => "zstdcat",
            Format::Xz => "xzcat",
        }
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Format::Plain => "uncompressed",
            Format::Gzip => "gzip",
            Format::Bzip2 => "bzip2",
            Format::Zstd => "zstd",
            Format::Xz => "xz",
        })
    }
}

/// Peek the stream's magic bytes and return a buffered reader over its
/// decoded contents, together with the detected format.
///
/// The peeked bytes are replayed in front of the stream, so nothing is lost.
///
/// # Errors
/// An I/O error while peeking, or `ErrorKind::Unsupported` when the format is
/// recognised but its Cargo feature is not compiled in.
pub fn open_decoded<R>(mut source: R) -> io::Result<(Box<dyn BufRead + Send>, Format)>
where
    R: Read + Send + 'static,
{
    let mut magic = [0u8; MAGIC_LEN];
    let n = read_up_to(&mut source, &mut magic)?;
    let format = Format::detect(&magic[..n]);
    let replayed = Cursor::new(magic).take(n as u64).chain(source);
    let compressed = BufReader::with_capacity(DECODE_BUF, replayed);

    let reader: Box<dyn BufRead + Send> = match format {
        Format::Plain => Box::new(compressed),
        Format::Gzip => Box::new(BufReader::with_capacity(
            DECODE_BUF,
            MultiGzDecoder::new(compressed),
        )),
        #[cfg(feature = "bz2")]
        Format::Bzip2 => Box::new(BufReader::with_capacity(
            DECODE_BUF,
            bzip2::read::MultiBzDecoder::new(compressed),
        )),
        #[cfg(feature = "zstd")]
        Format::Zstd => Box::new(BufReader::with_capacity(
            DECODE_BUF,
            zstd_frames::ZstdFrames::new(compressed),
        )),
        #[cfg(feature = "xz")]
        Format::Xz => Box::new(BufReader::with_capacity(
            DECODE_BUF,
            // `true`: keep going across concatenated XZ streams.
            lzma_rust2::XzReader::new(compressed, true),
        )),
        #[allow(unreachable_patterns)]
        other => return Err(unsupported(other)),
    };
    Ok((reader, format))
}

/// Error for a format whose decoder was not compiled in.
fn unsupported(format: Format) -> io::Error {
    io::Error::new(
        io::ErrorKind::Unsupported,
        format!(
            "input is {format}-compressed, but this build does not include the `{}` \
             feature; rebuild with `--features {}` or pass `--readFilesCommand {}`",
            format.feature().unwrap_or("?"),
            format.feature().unwrap_or("?"),
            format.command_hint(),
        ),
    )
}

/// Fill `buf` from `r` until it is full or `r` is exhausted, returning the
/// number of bytes read. Unlike `read_exact`, a short stream is not an error,
/// and unlike a single `read`, a pipe delivering bytes piecemeal is handled.
fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(k) => filled += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// Multi-frame Zstandard decoding on top of `ruzstd`.
///
/// `ruzstd::decoding::StreamingDecoder` decodes exactly one frame. A `.zst`
/// file is a sequence of frames: `zstd -T0` and `pzstd` write several, `pzstd`
/// adds skippable frames, and `cat a.zst b.zst` is valid input to `zstd -d`.
/// Stopping after the first frame would silently truncate the reads, the same
/// bug multi-member gzip had, so this adapter walks every frame, skips
/// skippable ones, and verifies each frame's content checksum when present.
#[cfg(feature = "zstd")]
mod zstd_frames {
    use ruzstd::decoding::{FrameDecoder, StreamingDecoder};
    use std::io::{self, Chain, Cursor, Read};

    /// The frame's magic number, already consumed, replayed ahead of the rest.
    type Source<R> = Chain<Cursor<[u8; 4]>, R>;

    enum State<R: Read> {
        /// Between frames. The `FrameDecoder` is kept to reuse its buffers.
        Idle(R, Box<FrameDecoder>),
        /// Inside a frame.
        Frame(Box<StreamingDecoder<Source<R>, FrameDecoder>>),
        /// End of input, or a previous error.
        Done,
    }

    pub struct ZstdFrames<R: Read> {
        state: State<R>,
    }

    impl<R: Read> ZstdFrames<R> {
        pub fn new(source: R) -> Self {
            Self {
                state: State::Idle(source, Box::new(FrameDecoder::new())),
            }
        }
    }

    fn truncated() -> io::Error {
        io::Error::new(io::ErrorKind::UnexpectedEof, "truncated zstd frame")
    }

    impl<R: Read> Read for ZstdFrames<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if buf.is_empty() {
                return Ok(0);
            }
            loop {
                // Any early return below leaves `Done`, so an error is sticky.
                match std::mem::replace(&mut self.state, State::Done) {
                    State::Done => return Ok(0),
                    State::Frame(mut dec) => {
                        let n = dec.read(buf)?;
                        if n > 0 {
                            self.state = State::Frame(dec);
                            return Ok(n);
                        }
                        let (source, fd) = dec.into_parts();
                        if let (Some(stored), Some(computed)) =
                            (fd.get_checksum_from_data(), fd.get_calculated_checksum())
                            && stored != computed
                        {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "zstd frame content checksum mismatch",
                            ));
                        }
                        let (_, inner) = source.into_inner();
                        self.state = State::Idle(inner, Box::new(fd));
                    }
                    State::Idle(mut inner, fd) => {
                        let mut magic = [0u8; 4];
                        match super::read_up_to(&mut inner, &mut magic)? {
                            0 => return Ok(0),
                            4 => {}
                            _ => return Err(truncated()),
                        }
                        let number = u32::from_le_bytes(magic);
                        if (0x184D_2A50..=0x184D_2A5F).contains(&number) {
                            let mut len = [0u8; 4];
                            inner.read_exact(&mut len)?;
                            let len = u64::from(u32::from_le_bytes(len));
                            let skipped = io::copy(&mut (&mut inner).take(len), &mut io::sink())?;
                            if skipped != len {
                                return Err(truncated());
                            }
                            self.state = State::Idle(inner, fd);
                            continue;
                        }
                        let dec = StreamingDecoder::new_with_decoder(
                            Cursor::new(magic).chain(inner),
                            *fd,
                        )
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                        self.state = State::Frame(Box::new(dec));
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const FASTQ: &[u8] = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nHHHH\n";

    fn decode(bytes: Vec<u8>) -> io::Result<(Vec<u8>, Format)> {
        let (mut r, f) = open_decoded(Cursor::new(bytes))?;
        let mut out = Vec::new();
        r.read_to_end(&mut out)?;
        Ok((out, f))
    }

    /// A reader that returns at most one byte per call, like a slow pipe.
    struct Trickle(Cursor<Vec<u8>>);
    impl Read for Trickle {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = buf.len().min(1);
            self.0.read(&mut buf[..n])
        }
    }

    #[test]
    fn detects_magic_bytes() {
        assert_eq!(Format::detect(b"@r1\nACGT"), Format::Plain);
        assert_eq!(Format::detect(b""), Format::Plain);
        assert_eq!(Format::detect(&[0x1f, 0x8b, 0x08]), Format::Gzip);
        assert_eq!(Format::detect(b"BZh91AY"), Format::Bzip2);
        assert_eq!(Format::detect(&[0x28, 0xb5, 0x2f, 0xfd, 0]), Format::Zstd);
        assert_eq!(Format::detect(&[0x5a, 0x2a, 0x4d, 0x18]), Format::Zstd);
        assert_eq!(
            Format::detect(&[0xfd, b'7', b'z', b'X', b'Z', 0x00]),
            Format::Xz
        );
        // A truncated signature is not a match.
        assert_eq!(Format::detect(&[0xfd, b'7', b'z']), Format::Plain);
    }

    #[test]
    fn plain_passes_through_including_short_input() {
        assert_eq!(
            decode(FASTQ.to_vec()).unwrap(),
            (FASTQ.to_vec(), Format::Plain)
        );
        assert_eq!(
            decode(b"@r".to_vec()).unwrap(),
            (b"@r".to_vec(), Format::Plain)
        );
        assert_eq!(decode(Vec::new()).unwrap(), (Vec::new(), Format::Plain));
    }

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn gzip_multi_member_and_trickled_source() {
        let (a, b) = FASTQ.split_at(16);
        let mut bytes = gzip(a);
        bytes.extend(gzip(b));

        let (out, f) = decode(bytes.clone()).unwrap();
        assert_eq!((out.as_slice(), f), (FASTQ, Format::Gzip));

        // Magic bytes arriving one at a time must still be detected.
        let (mut r, f) = open_decoded(Trickle(Cursor::new(bytes))).unwrap();
        let mut out = Vec::new();
        r.read_to_end(&mut out).unwrap();
        assert_eq!((out.as_slice(), f), (FASTQ, Format::Gzip));
    }

    #[cfg(feature = "bz2")]
    fn bzip2(data: &[u8]) -> Vec<u8> {
        let mut e = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    #[cfg(feature = "bz2")]
    #[test]
    fn bzip2_round_trip_multi_stream() {
        let (a, b) = FASTQ.split_at(16);
        let mut bytes = bzip2(a);
        bytes.extend(bzip2(b));
        let (out, f) = decode(bytes).unwrap();
        assert_eq!((out.as_slice(), f), (FASTQ, Format::Bzip2));
    }

    #[cfg(feature = "zstd")]
    fn zstd(data: &[u8]) -> Vec<u8> {
        ruzstd::encoding::compress_to_vec(data, ruzstd::encoding::CompressionLevel::Fastest)
    }

    #[cfg(feature = "zstd")]
    #[test]
    fn zstd_round_trip_multi_frame_with_skippable() {
        let (a, b) = FASTQ.split_at(16);
        // A skippable frame first, as pzstd writes, then two data frames.
        let mut bytes = vec![0x50, 0x2a, 0x4d, 0x18, 3, 0, 0, 0, 9, 9, 9];
        bytes.extend(zstd(a));
        bytes.extend(zstd(b));
        let (out, f) = decode(bytes).unwrap();
        assert_eq!((out.as_slice(), f), (FASTQ, Format::Zstd));

        // Large enough to span several blocks.
        let big: Vec<u8> = FASTQ.iter().copied().cycle().take(1 << 20).collect();
        let (out, _) = decode(zstd(&big)).unwrap();
        assert_eq!(out, big);
    }

    #[cfg(feature = "zstd")]
    #[test]
    fn zstd_truncated_is_an_error() {
        let mut bytes = zstd(FASTQ);
        bytes.extend([0x28, 0xb5]);
        assert!(decode(bytes).is_err());
    }

    fn xz(data: &[u8]) -> Vec<u8> {
        let mut w =
            lzma_rust2::XzWriter::new(Vec::new(), lzma_rust2::XzOptions::with_preset(6)).unwrap();
        w.write_all(data).unwrap();
        w.finish().unwrap()
    }

    #[cfg(feature = "xz")]
    #[test]
    fn xz_round_trip_multi_stream() {
        let (a, b) = FASTQ.split_at(16);
        let mut bytes = xz(a);
        bytes.extend(xz(b));
        let (out, f) = decode(bytes).unwrap();
        assert_eq!((out.as_slice(), f), (FASTQ, Format::Xz));
    }

    /// Without the feature, a detected format is a clear error rather than
    /// compressed bytes fed to the FASTQ parser.
    #[cfg(not(feature = "xz"))]
    #[test]
    fn xz_without_feature_is_unsupported() {
        let err = decode(xz(FASTQ)).expect_err("must not decode");
        assert_eq!(err.kind(), io::ErrorKind::Unsupported);
        assert!(err.to_string().contains("--features xz"));
        assert!(err.to_string().contains("xzcat"));
    }
}
