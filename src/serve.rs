//! Resident-database server.
//!
//! Loads the database once into an [`Engine`] and then answers many sequence batches over either
//! a Unix domain socket or stdin/stdout. Every connection holds only an `Arc<Engine>` clone, and
//! each batch is dispatched across the engine's shared rayon pool, so concurrent connections and
//! single large batches both keep all cores busy.
//!
//! # Wire format
//!
//! All integers are little-endian `u32`. A connection may send many requests back to back.
//!
//! ```text
//! Request:   n_seqs, then n_seqs * ( len, len bytes of sequence )
//! Response:  n_seqs, then n_seqs * ( n_kmers, n_kmers * taxid )
//! Error:     0xFFFF_FFFF, then msg_len, msg_len bytes of UTF-8
//! ```
//!
//! Binary rather than JSON because the payload is per-k-mer `u32` vectors on the hot path; a
//! client can decode a response body with a single `numpy.frombuffer`.

use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::engine::{Engine, EngineConfig};

/// Sentinel in the `n_seqs` slot marking an error response.
const ERROR_TAG: u32 = u32::MAX;

/// Guards against a malformed or hostile frame causing a huge allocation.
const MAX_SEQS_PER_BATCH: u32 = 10_000_000;
const MAX_SEQ_LEN: u32 = 1 << 28; // 256 MiB

pub struct ServeConfig {
    pub db_dir: String,
    pub threads: usize,
    pub use_accessions: bool,
    pub background: bool,
    /// Unix socket path. Ignored when `stdio` is set.
    pub socket: Option<String>,
    /// Serve a single session over stdin/stdout instead of a socket.
    pub stdio: bool,
}

pub fn run_serve(config: ServeConfig) -> Result<()> {
    let engine = Engine::load(&EngineConfig {
        db_dir: config.db_dir.clone(),
        accessions: config.use_accessions,
        background: config.background,
        threads: config.threads,
    })?;

    eprintln!(
        "Database loaded in {:.3}s (k={}, l={}), {} worker threads",
        engine.load_time().as_secs_f64(),
        engine.k(),
        engine.l(),
        config.threads
    );

    if config.stdio {
        eprintln!("READY stdio");
        let stdin = std::io::stdin();
        let stdout = std::io::stdout();
        return serve_connection(&engine, stdin.lock(), stdout.lock());
    }

    let socket_path = config
        .socket
        .as_deref()
        .context("serve needs either --socket <path> or --stdio")?;
    serve_socket(engine, socket_path)
}

#[cfg(unix)]
fn serve_socket(engine: Arc<Engine>, socket_path: &str) -> Result<()> {
    use std::os::unix::net::UnixListener;

    // A leftover socket file from a previous run would make bind() fail with EADDRINUSE.
    match std::fs::remove_file(socket_path) {
        Ok(()) => eprintln!("Removed stale socket {}", socket_path),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).context(format!("Can't clear existing socket {}", socket_path)),
    }

    let listener = UnixListener::bind(socket_path)
        .with_context(|| format!("Can't bind socket {}", socket_path))?;
    eprintln!("READY {}", socket_path);

    for stream in listener.incoming() {
        let stream = match stream {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Rejected connection: {}", e);
                continue;
            }
        };

        // Arc clone only; the database itself is never copied.
        let engine = Arc::clone(&engine);
        std::thread::spawn(move || {
            let write_half = match stream.try_clone() {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("Connection dropped: {}", e);
                    return;
                }
            };
            if let Err(e) = serve_connection(&engine, stream, write_half) {
                eprintln!("Connection ended: {}", e);
            }
        });
    }

    Ok(())
}

#[cfg(not(unix))]
fn serve_socket(_engine: Arc<Engine>, _socket_path: &str) -> Result<()> {
    anyhow::bail!("Unix socket serving is only available on unix; use --stdio instead")
}

/// Handles requests on one connection until the peer closes it.
fn serve_connection<R: Read, W: Write>(engine: &Engine, read: R, write: W) -> Result<()> {
    let mut reader = BufReader::with_capacity(1 << 20, read);
    let mut writer = BufWriter::with_capacity(1 << 20, write);

    loop {
        let seqs = match read_request(&mut reader)? {
            Some(seqs) => seqs,
            None => return Ok(()), // clean EOF between requests
        };

        let start = std::time::Instant::now();
        let results = engine.query_batch(&seqs);
        eprintln!(
            "Answered {} sequences in {:.3}ms",
            results.len(),
            start.elapsed().as_secs_f64() * 1000.0
        );

        write_response(&mut writer, &results)?;
        writer.flush()?;
    }
}

/// Reads one request frame. `Ok(None)` means the peer closed cleanly at a frame boundary.
fn read_request<R: Read>(reader: &mut R) -> Result<Option<Vec<Vec<u8>>>> {
    let n_seqs = match read_u32_opt(reader)? {
        Some(n) => n,
        None => return Ok(None),
    };

    if n_seqs > MAX_SEQS_PER_BATCH {
        anyhow::bail!("Batch of {} sequences exceeds the {} limit", n_seqs, MAX_SEQS_PER_BATCH);
    }

    let mut seqs = Vec::with_capacity(n_seqs as usize);
    for i in 0..n_seqs {
        let len = read_u32(reader).with_context(|| format!("Truncated length for sequence {}", i))?;
        if len > MAX_SEQ_LEN {
            anyhow::bail!("Sequence {} of {} bytes exceeds the {} limit", i, len, MAX_SEQ_LEN);
        }
        let mut seq = vec![0u8; len as usize];
        reader
            .read_exact(&mut seq)
            .with_context(|| format!("Truncated body for sequence {}", i))?;
        seqs.push(seq);
    }

    Ok(Some(seqs))
}

fn write_response<W: Write>(writer: &mut W, results: &[Vec<u32>]) -> Result<()> {
    writer.write_all(&(results.len() as u32).to_le_bytes())?;
    for taxids in results {
        writer.write_all(&(taxids.len() as u32).to_le_bytes())?;
        for &taxid in taxids {
            writer.write_all(&taxid.to_le_bytes())?;
        }
    }
    Ok(())
}

/// Reports a failure to the peer in-band so a client blocked on a read doesn't hang.
#[allow(dead_code)]
fn write_error<W: Write>(writer: &mut W, msg: &str) -> Result<()> {
    let bytes = msg.as_bytes();
    writer.write_all(&ERROR_TAG.to_le_bytes())?;
    writer.write_all(&(bytes.len() as u32).to_le_bytes())?;
    writer.write_all(bytes)?;
    writer.flush()?;
    Ok(())
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

/// Like [`read_u32`] but distinguishes a clean EOF at a frame boundary from a truncated frame.
fn read_u32_opt<R: Read>(reader: &mut R) -> Result<Option<u32>> {
    let mut buf = [0u8; 4];
    let mut filled = 0;
    while filled < 4 {
        match reader.read(&mut buf[filled..]) {
            Ok(0) if filled == 0 => return Ok(None),
            Ok(0) => anyhow::bail!("Truncated frame header: got {} of 4 bytes", filled),
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    Ok(Some(u32::from_le_bytes(buf)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(seqs: &[&[u8]]) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&(seqs.len() as u32).to_le_bytes());
        for s in seqs {
            buf.extend_from_slice(&(s.len() as u32).to_le_bytes());
            buf.extend_from_slice(s);
        }
        buf
    }

    #[test]
    fn roundtrips_a_request_frame() {
        let bytes = frame(&[b"ACGT", b"TTTTTT"]);
        let parsed = read_request(&mut bytes.as_slice()).unwrap().unwrap();
        assert_eq!(parsed, vec![b"ACGT".to_vec(), b"TTTTTT".to_vec()]);
    }

    #[test]
    fn clean_eof_ends_the_session() {
        assert!(read_request(&mut [].as_slice()).unwrap().is_none());
    }

    #[test]
    fn truncated_header_is_an_error() {
        assert!(read_request(&mut [1u8, 2].as_slice()).is_err());
    }

    #[test]
    fn truncated_body_is_an_error() {
        let mut bytes = frame(&[b"ACGT"]);
        bytes.truncate(bytes.len() - 2);
        assert!(read_request(&mut bytes.as_slice()).is_err());
    }

    #[test]
    fn oversized_batch_is_rejected_without_allocating() {
        let bytes = (MAX_SEQS_PER_BATCH + 1).to_le_bytes();
        assert!(read_request(&mut bytes.as_slice()).is_err());
    }

    #[test]
    fn writes_the_expected_response_layout() {
        let mut out = Vec::new();
        write_response(&mut out, &[vec![0u32, 9606], vec![]]).unwrap();

        let mut expected = Vec::new();
        expected.extend_from_slice(&2u32.to_le_bytes()); // n_seqs
        expected.extend_from_slice(&2u32.to_le_bytes()); // n_kmers seq 0
        expected.extend_from_slice(&0u32.to_le_bytes());
        expected.extend_from_slice(&9606u32.to_le_bytes());
        expected.extend_from_slice(&0u32.to_le_bytes()); // n_kmers seq 1
        assert_eq!(out, expected);
    }
}
