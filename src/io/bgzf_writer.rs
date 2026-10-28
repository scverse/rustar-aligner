//! Order-preserving BGZF writer with optional parallel block compression.
//!
//! The align pipelines hand finished records to a single writer thread
//! (#223). With BAM output that thread used to encode every record *and*
//! DEFLATE every 64 KiB block itself, so compression was a serial stage that
//! capped scaling at roughly 4-8 `--runThreadN`.
//!
//! [`BgzfWriter`] keeps the calling thread's job down to filling 64 KiB
//! staging blocks. Full blocks are compressed by a small pool of dedicated
//! worker threads and written, strictly in submission order, by one I/O
//! thread. The pool is private to the writer (plain `std::thread`s, not
//! rayon), so it never queues behind alignment work on rayon's global pool;
//! idle workers are parked on a channel and cost nothing.
//!
//! # Output identity
//!
//! The byte stream is identical to `noodles_bgzf::io::Writer` built with the
//! same compression level (with the `libdeflate` feature this crate enables):
//!
//! * block boundaries depend only on the uncompressed byte stream: a block is
//!   emitted when the staging buffer reaches [`MAX_BUF_SIZE`] or on an
//!   explicit `flush()`, exactly as noodles does, never on thread timing;
//! * each block is compressed with a fresh-state libdeflate compressor at the
//!   same level, which is deterministic;
//! * the gzip member header/trailer and the EOF marker are the ones noodles
//!   writes (unit tests below compare against noodles byte for byte).
//!
//! So `--outSAMtype BAM ...` output is byte-identical for every thread count.

use std::io::{self, Write};
use std::mem;
use std::sync::mpsc::{self, Receiver, Sender, SyncSender};
use std::sync::{Arc, Mutex, PoisonError};
use std::thread::JoinHandle;

use libdeflater::{CompressionLvl, Compressor, Crc};
use noodles::bgzf::io::writer::CompressionLevel;

/// gzip header (10) + XLEN (2) + BGZF `BC` subfield (6).
const BGZF_HEADER_SIZE: usize = 18;
/// gzip trailer: CRC32 + ISIZE.
const GZ_TRAILER_SIZE: usize = 8;
/// Worst-case DEFLATE overhead of a stored (level 0) 64 KiB block, as
/// budgeted by noodles-bgzf.
const COMPRESSION_LEVEL_0_OVERHEAD: usize = 15;
/// Largest uncompressed payload per BGZF block. Must equal noodles-bgzf's
/// `MAX_BUF_SIZE` so that block boundaries (and therefore bytes) match.
pub const MAX_BUF_SIZE: usize =
    (1 << 16) - BGZF_HEADER_SIZE - GZ_TRAILER_SIZE - COMPRESSION_LEVEL_0_OVERHEAD;

/// SAM spec 4.1.2 end-of-file marker block.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Blocks allowed in flight per compression worker before the producer
/// blocks. Bounds memory at roughly `workers * 4 * 128 KiB`.
const IN_FLIGHT_PER_WORKER: usize = 4;

/// A compressed block: DEFLATE payload, CRC32 and uncompressed size.
struct Frame {
    data: Vec<u8>,
    crc32: u32,
    isize: usize,
}

struct Job {
    data: Vec<u8>,
    reply: SyncSender<io::Result<Frame>>,
}

enum Mode<W: Write + Send + 'static> {
    /// Compress and write on the calling thread (`--runThreadN 1`).
    Inline {
        inner: W,
        compressor: Box<Compressor>,
    },
    /// Compress on `workers`, write in order on `io_thread`.
    Threaded {
        job_tx: Sender<Job>,
        order_tx: SyncSender<Receiver<io::Result<Frame>>>,
        workers: Vec<JoinHandle<()>>,
        io_thread: JoinHandle<io::Result<W>>,
    },
    /// `finish()` has run; the inner writer is gone.
    Done,
}

/// BGZF writer; see the module docs.
pub struct BgzfWriter<W: Write + Send + 'static> {
    staging: Vec<u8>,
    mode: Mode<W>,
}

impl<W: Write + Send + 'static> BgzfWriter<W> {
    /// Create a writer that compresses at `level`.
    ///
    /// `threads <= 1` compresses inline on the calling thread. Otherwise
    /// `threads` compression workers plus one I/O thread are spawned.
    pub fn new(inner: W, level: CompressionLevel, threads: usize) -> io::Result<Self> {
        let lvl = CompressionLvl::from(level);
        let mode = if threads <= 1 {
            Mode::Inline {
                inner,
                compressor: Box::new(Compressor::new(lvl)),
            }
        } else {
            Self::spawn(inner, lvl, threads)?
        };
        Ok(Self {
            staging: Vec::with_capacity(MAX_BUF_SIZE),
            mode,
        })
    }

    fn spawn(inner: W, lvl: CompressionLvl, threads: usize) -> io::Result<Mode<W>> {
        let (job_tx, job_rx) = mpsc::channel::<Job>();
        let job_rx = Arc::new(Mutex::new(job_rx));
        let mut workers = Vec::with_capacity(threads);
        for i in 0..threads {
            let job_rx = Arc::clone(&job_rx);
            workers.push(
                std::thread::Builder::new()
                    .name(format!("bgzf-deflate-{i}"))
                    .spawn(move || {
                        let mut compressor = Compressor::new(lvl);
                        loop {
                            // Hold the lock only for the dequeue, not the work.
                            let job = job_rx.lock().unwrap_or_else(PoisonError::into_inner).recv();
                            let Ok(job) = job else { break };
                            let frame = compress_block(&mut compressor, &job.data);
                            // The I/O thread may have quit on an error; the
                            // producer surfaces that, so dropping is fine.
                            let _ = job.reply.send(frame);
                        }
                    })?,
            );
        }

        let (order_tx, order_rx) =
            mpsc::sync_channel::<Receiver<io::Result<Frame>>>(threads * IN_FLIGHT_PER_WORKER);
        let io_thread = std::thread::Builder::new()
            .name("bgzf-write".into())
            .spawn(move || -> io::Result<W> {
                let mut inner = inner;
                // Receivers arrive in submission order; waiting on each in
                // turn is what keeps the output order deterministic.
                for reply in order_rx {
                    let frame = reply.recv().map_err(|_| {
                        io::Error::other("BGZF compression worker exited unexpectedly")
                    })??;
                    write_frame(&mut inner, &frame)?;
                }
                inner.write_all(&BGZF_EOF)?;
                inner.flush()?;
                Ok(inner)
            })?;

        Ok(Mode::Threaded {
            job_tx,
            order_tx,
            workers,
            io_thread,
        })
    }

    /// Emit the staged bytes as one BGZF block.
    fn emit_block(&mut self) -> io::Result<()> {
        match &mut self.mode {
            Mode::Inline { inner, compressor } => {
                let frame = compress_block(compressor, &self.staging)?;
                self.staging.clear();
                write_frame(inner, &frame)
            }
            Mode::Threaded {
                job_tx, order_tx, ..
            } => {
                let data = mem::replace(&mut self.staging, Vec::with_capacity(MAX_BUF_SIZE));
                let (reply_tx, reply_rx) = mpsc::sync_channel(1);
                if order_tx.send(reply_rx).is_err()
                    || job_tx
                        .send(Job {
                            data,
                            reply: reply_tx,
                        })
                        .is_err()
                {
                    // The I/O thread stopped early: report its error.
                    return match self.finish_inner() {
                        Err(e) => Err(e),
                        Ok(_) => Err(io::Error::other("BGZF writer stopped unexpectedly")),
                    };
                }
                Ok(())
            }
            Mode::Done => Err(io::Error::other("BGZF writer already finished")),
        }
    }

    /// Flush staged data, write the EOF marker, and return the inner writer.
    ///
    /// Blocks until every in-flight block has been compressed and written.
    pub fn finish(&mut self) -> io::Result<W> {
        if !self.staging.is_empty() {
            self.emit_block()?;
        }
        self.finish_inner()
    }

    fn finish_inner(&mut self) -> io::Result<W> {
        match mem::replace(&mut self.mode, Mode::Done) {
            Mode::Inline { mut inner, .. } => {
                inner.write_all(&BGZF_EOF)?;
                inner.flush()?;
                Ok(inner)
            }
            Mode::Threaded {
                job_tx,
                order_tx,
                workers,
                io_thread,
            } => {
                drop(order_tx);
                drop(job_tx);
                let result = io_thread
                    .join()
                    .unwrap_or_else(|_| Err(io::Error::other("BGZF I/O thread panicked")));
                for w in workers {
                    let _ = w.join();
                }
                result
            }
            Mode::Done => Err(io::Error::other("BGZF writer already finished")),
        }
    }
}

impl<W: Write + Send + 'static> Write for BgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let amt = (MAX_BUF_SIZE - self.staging.len()).min(buf.len());
        self.staging.extend_from_slice(&buf[..amt]);
        if self.staging.len() >= MAX_BUF_SIZE {
            self.emit_block()?;
        }
        Ok(amt)
    }

    /// Like noodles: closes the current block, does not flush the inner
    /// writer (that happens in [`BgzfWriter::finish`]).
    fn flush(&mut self) -> io::Result<()> {
        if self.staging.is_empty() {
            Ok(())
        } else {
            self.emit_block()
        }
    }
}

impl<W: Write + Send + 'static> Drop for BgzfWriter<W> {
    /// Same contract as noodles' writer: an unfinished stream is finished
    /// (EOF marker included) on drop, with errors discarded.
    fn drop(&mut self) {
        if !matches!(self.mode, Mode::Done) {
            let _ = self.finish();
        }
    }
}

fn compress_block(compressor: &mut Compressor, src: &[u8]) -> io::Result<Frame> {
    let mut dst = vec![0u8; compressor.deflate_compress_bound(src.len())];
    let len = compressor
        .deflate_compress(src, &mut dst)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.truncate(len);
    let mut crc = Crc::new();
    crc.update(src);
    Ok(Frame {
        data: dst,
        crc32: crc.sum(),
        isize: src.len(),
    })
}

fn write_frame<W: Write>(writer: &mut W, frame: &Frame) -> io::Result<()> {
    let block_size = BGZF_HEADER_SIZE + frame.data.len() + GZ_TRAILER_SIZE;
    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let isize =
        u32::try_from(frame.isize).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut header = [0u8; BGZF_HEADER_SIZE];
    header[..16].copy_from_slice(&[
        0x1f, 0x8b, // gzip magic
        0x08, // CM = DEFLATE
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME = 0
        0x00, // XFL
        0xff, // OS = unknown
        0x06, 0x00, // XLEN = 6
        b'B', b'C', // BGZF subfield id
        0x02, 0x00, // SLEN = 2
    ]);
    header[16..].copy_from_slice(&bsize.to_le_bytes());

    writer.write_all(&header)?;
    writer.write_all(&frame.data)?;
    writer.write_all(&frame.crc32.to_le_bytes())?;
    writer.write_all(&isize.to_le_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::bgzf;

    /// Pseudo-random but compressible bytes (BAM-like: repetitive text).
    fn payload(n: usize) -> Vec<u8> {
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        (0..n)
            .map(|i| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                if i % 7 == 0 {
                    b"ACGT"[(x & 3) as usize]
                } else {
                    (x >> 56) as u8 & 0x1f
                }
            })
            .collect()
    }

    /// Write `data` in irregular chunks with a few explicit flushes, the
    /// same way through both writers.
    fn drive<Wr: Write>(w: &mut Wr, data: &[u8]) {
        let mut off = 0;
        let mut step = 1usize;
        while off < data.len() {
            let end = (off + step).min(data.len());
            w.write_all(&data[off..end]).unwrap();
            if step.is_multiple_of(11) {
                w.flush().unwrap();
            }
            off = end;
            step = step * 7 % 20_011 + 1;
        }
    }

    fn noodles_bytes(data: &[u8], level: CompressionLevel) -> Vec<u8> {
        let mut w = bgzf::io::writer::Builder::default()
            .set_compression_level(level)
            .build_from_writer(Vec::new());
        drive(&mut w, data);
        w.finish().unwrap()
    }

    fn ours_bytes(data: &[u8], level: CompressionLevel, threads: usize) -> Vec<u8> {
        let mut w = BgzfWriter::new(Vec::new(), level, threads).unwrap();
        drive(&mut w, data);
        w.finish().unwrap()
    }

    #[test]
    fn max_buf_size_matches_noodles() {
        // noodles: BGZF_MAX_ISIZE(65536) - 18 - 8 - 15
        assert_eq!(MAX_BUF_SIZE, 65_495);
    }

    #[test]
    fn byte_identical_to_noodles() {
        let sizes = [
            0,
            1,
            1000,
            MAX_BUF_SIZE,
            MAX_BUF_SIZE + 1,
            3 * MAX_BUF_SIZE + 17,
            700_001,
        ];
        let levels = [0u8, 1, 6, 12].map(|l| CompressionLevel::try_from(l).unwrap());
        for &n in &sizes {
            let data = payload(n);
            for &level in &levels {
                let expected = noodles_bytes(&data, level);
                for threads in [1, 2, 5] {
                    assert_eq!(
                        ours_bytes(&data, level, threads),
                        expected,
                        "size={n} level={} threads={threads}",
                        level.get()
                    );
                }
            }
        }
    }

    #[test]
    fn roundtrips_through_noodles_reader() {
        use std::io::Read;
        let data = payload(1_000_003);
        let bytes = ours_bytes(&data, CompressionLevel::FAST, 4);
        let mut out = Vec::new();
        bgzf::io::Reader::new(&bytes[..])
            .read_to_end(&mut out)
            .unwrap();
        assert_eq!(out, data);
    }

    #[test]
    fn drop_without_finish_writes_eof() {
        struct Shared(Arc<Mutex<Vec<u8>>>);
        impl Write for Shared {
            fn write(&mut self, b: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(b);
                Ok(b.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        for threads in [1, 3] {
            let sink = Arc::new(Mutex::new(Vec::new()));
            {
                let mut w =
                    BgzfWriter::new(Shared(Arc::clone(&sink)), CompressionLevel::FAST, threads)
                        .unwrap();
                w.write_all(b"noodles").unwrap();
            }
            let bytes = sink.lock().unwrap().clone();
            assert!(bytes.ends_with(&BGZF_EOF), "threads={threads}");
            assert_eq!(bytes, noodles_bytes(b"noodles", CompressionLevel::FAST));
        }
    }

    #[test]
    fn io_error_is_reported() {
        struct Failing;
        impl Write for Failing {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::other("disk full"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        for threads in [1, 3] {
            let mut w = BgzfWriter::new(Failing, CompressionLevel::FAST, threads).unwrap();
            let data = payload(10 * MAX_BUF_SIZE);
            let write_res = w.write_all(&data);
            let finish_res = w.finish();
            assert!(
                write_res.is_err() || finish_res.is_err(),
                "threads={threads}: error swallowed"
            );
        }
    }
}
