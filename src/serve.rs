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
//! With -a:  n_seqs, then n_seqs * ( n_kmers, n_classes,
//!                                  n_classes * ( class_id, hit_kmers ) )
//! Table request: 0xFFFF_FFFE (instead of n_seqs, no body)
//! Table response: n_accessions, n_accessions * ( name_len, name bytes ),
//!                 n_classes, n_classes * ( n_ids, n_ids * accession_id )
//! Error:     0xFFFF_FFFF, then msg_len, msg_len bytes of UTF-8
//! ```
//!
//! Classes are sorted by ID within each sequence; unclassed k-mers still count in n_kmers.
//! Table positions are zero-based accession/class IDs, including the empty class at ID 0.
//! Accessions are per minimizer and therefore candidates, not exact full-k-mer matches.
//! Class IDs are valid only for the database that issued them. Table requests require -a.
//!
//! Binary rather than JSON because the payload is per-k-mer `u32` vectors on the hot path; a
//! client can decode a response body with a single `numpy.frombuffer`.

use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::engine::{Engine, EngineConfig};

mod progress;
use progress::ServerProgress;

/// Sentinel in the `n_seqs` slot marking an error response.
const ERROR_TAG: u32 = u32::MAX;
const CLASS_TABLE_TAG: u32 = u32::MAX - 1;

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
    if config.use_accessions {
        anyhow::ensure!(
            engine.db().accessions.is_some() && engine.acc_registry().is_some(),
            "serve -a requires an accession index and registry; build the database with -a"
        );
    }

    eprintln!(
        "Database loaded in {:.3}s (k={}, l={}), {} worker threads",
        engine.load_time().as_secs_f64(),
        engine.k(),
        engine.l(),
        config.threads
    );

    let progress = Arc::new(ServerProgress::new());
    if config.stdio {
        progress.message("READY stdio");
        let stdin = std::io::stdin();
        let stdout = std::io::stdout();
        let result = serve_connection(&engine, stdin.lock(), stdout.lock(), &progress);
        if result.is_ok() {
            progress.finish();
        } else {
            progress.end_line();
        }
        return result;
    }

    let socket_path = config
        .socket
        .as_deref()
        .context("serve needs either --socket <path> or --stdio")?;
    serve_socket(engine, socket_path, progress)
}

#[cfg(unix)]
fn serve_socket(
    engine: Arc<Engine>,
    socket_path: &str,
    progress: Arc<ServerProgress>,
) -> Result<()> {
    use std::os::unix::net::UnixListener;

    // A leftover socket file from a previous run would make bind() fail with EADDRINUSE.
    match std::fs::remove_file(socket_path) {
        Ok(()) => progress.message(&format!("Removed stale socket {}", socket_path)),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).context(format!("Can't clear existing socket {}", socket_path)),
    }

    let listener = UnixListener::bind(socket_path)
        .with_context(|| format!("Can't bind socket {}", socket_path))?;
    progress.message(&format!("READY {}", socket_path));

    for stream in listener.incoming() {
        let stream = match stream {
            Ok(s) => s,
            Err(e) => {
                progress.message(&format!("Rejected connection: {}", e));
                continue;
            }
        };

        // Arc clone only; the database itself is never copied.
        let engine = Arc::clone(&engine);
        let progress = Arc::clone(&progress);
        std::thread::spawn(move || {
            let write_half = match stream.try_clone() {
                Ok(s) => s,
                Err(e) => {
                    progress.message(&format!("Connection dropped: {}", e));
                    return;
                }
            };
            if let Err(e) = serve_connection(&engine, stream, write_half, &progress) {
                progress.message(&format!("Connection ended: {}", e));
            }
        });
    }

    Ok(())
}

#[cfg(not(unix))]
fn serve_socket(
    _engine: Arc<Engine>,
    _socket_path: &str,
    _progress: Arc<ServerProgress>,
) -> Result<()> {
    anyhow::bail!("Unix socket serving is only available on unix; use --stdio instead")
}

/// Handles requests on one connection until the peer closes it.
fn serve_connection<R: Read, W: Write>(
    engine: &Engine,
    read: R,
    write: W,
    progress: &ServerProgress,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(1 << 20, read);
    // Keep block buffering for stdio's underlying line-buffered writer as well as sockets.
    let mut writer = BufWriter::with_capacity(1 << 20, write);

    let mut batch = RequestBatch::default();
    let mut response = Vec::new();
    let mut scratch = Vec::new();
    loop {
        match read_request_kind(&mut reader, &mut batch)? {
            None => return Ok(()),
            Some(RequestKind::ClassTable) => {
                if engine.db().accessions.is_some() && engine.acc_registry().is_some() {
                    write_class_table(engine, &mut writer)?;
                } else {
                    write_error(&mut writer, "Class table requires serve -a")?;
                }
                writer.flush()?;
                continue;
            }
            Some(RequestKind::Batch) => {}
        }

        progress.start();
        let sequences = batch.ranges.len();
        let bases = batch.ranges.last().map_or(0, |r| r.end as u64);
        if engine.db().accessions.is_some() {
            engine.query_classes_wire(&batch.data, &batch.ranges, &mut response);
        } else {
            engine.query_batch_wire(&batch.data, &batch.ranges, &mut response, &mut scratch)?;
        }

        writer.write_all(&response)?;
        writer.flush()?;
        progress.completed(sequences, bases);
    }
}

fn write_class_table<W: Write>(engine: &Engine, writer: &mut W) -> Result<()> {
    let registry = engine
        .acc_registry()
        .as_ref()
        .context("Missing accession registry")?;
    let index = engine
        .db()
        .accessions
        .as_ref()
        .context("Missing accession index")?;
    writer.write_all(&(registry.len() as u32).to_le_bytes())?;
    for id in 0..registry.len() {
        let name = registry.get_name(id as u32).as_bytes();
        writer.write_all(&(name.len() as u32).to_le_bytes())?;
        writer.write_all(name)?;
    }
    writer.write_all(&(index.num_classes() as u32).to_le_bytes())?;
    let mut ids = Vec::new();
    for id in 0..index.num_classes() {
        index.decode_class_into(id as u32, &mut ids);
        writer.write_all(&(ids.len() as u32).to_le_bytes())?;
        for accession in &ids {
            writer.write_all(&accession.to_le_bytes())?;
        }
    }
    Ok(())
}

#[derive(Default)]
struct RequestBatch {
    data: Vec<u8>,
    ranges: Vec<std::ops::Range<usize>>,
}

#[derive(PartialEq)]
enum RequestKind {
    Batch,
    ClassTable,
}

/// Reuse contiguous storage; None means clean EOF at a frame boundary.
fn read_request_kind<R: Read>(
    reader: &mut R,
    batch: &mut RequestBatch,
) -> Result<Option<RequestKind>> {
    let n_seqs = match read_u32_opt(reader)? {
        Some(n) => n,
        None => return Ok(None),
    };
    if n_seqs == CLASS_TABLE_TAG {
        return Ok(Some(RequestKind::ClassTable));
    }
    if n_seqs > MAX_SEQS_PER_BATCH {
        anyhow::bail!(
            "Batch of {} sequences exceeds the {} limit",
            n_seqs,
            MAX_SEQS_PER_BATCH
        );
    }
    // Keep the initialized allocation between requests as well as its capacity.
    batch.ranges.clear();
    batch.ranges.reserve(n_seqs as usize);
    let mut end = 0usize;
    for i in 0..n_seqs {
        let len =
            read_u32(reader).with_context(|| format!("Truncated length for sequence {}", i))?;
        if len > MAX_SEQ_LEN {
            anyhow::bail!(
                "Sequence {} of {} bytes exceeds the {} limit",
                i,
                len,
                MAX_SEQ_LEN
            );
        }
        let start = end;
        end = start
            .checked_add(len as usize)
            .context("Request size overflow")?;
        if batch.data.len() < end {
            batch.data.resize(end, 0);
        }
        reader
            .read_exact(&mut batch.data[start..end])
            .with_context(|| format!("Truncated body for sequence {}", i))?;
        batch.ranges.push(start..end);
    }
    Ok(Some(RequestKind::Batch))
}

#[cfg(test)]
fn read_request_into<R: Read>(reader: &mut R, batch: &mut RequestBatch) -> Result<bool> {
    Ok(read_request_kind(reader, batch)? == Some(RequestKind::Batch))
}

#[cfg(test)]
fn read_request<R: Read>(reader: &mut R) -> Result<Option<Vec<Vec<u8>>>> {
    let mut batch = RequestBatch::default();
    if !read_request_into(reader, &mut batch)? {
        return Ok(None);
    }
    Ok(Some(
        batch
            .ranges
            .iter()
            .map(|r| batch.data[r.clone()].to_vec())
            .collect(),
    ))
}

#[cfg(test)]
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

    fn test_engine() -> Arc<Engine> {
        Engine::load(&EngineConfig {
            db_dir: "data/influenza_test_db".to_string(),
            accessions: false,
            background: false,
            threads: 1,
        })
        .unwrap()
    }

    #[test]
    fn session_responses_are_unchanged_and_count_only_completed_work() {
        let engine = test_engine();
        let seqs: Vec<Vec<u8>> = vec![b"ACGT".to_vec(), vec![b'A'; 50]];
        let request = frame(&[&seqs[0], &seqs[1]]);
        let mut requests = request.clone();
        requests.extend_from_slice(&request);
        let mut expected = Vec::new();
        write_response(&mut expected, &engine.query_batch(&seqs)).unwrap();
        let expected = expected.repeat(2);
        let mut output = Vec::new();
        let progress = ServerProgress::new();
        serve_connection(&engine, requests.as_slice(), &mut output, &progress).unwrap();
        assert_eq!(output, expected);
        assert_eq!(progress.totals(), (2, 4, 108));
    }

    #[test]
    fn reused_batches_preserve_every_kmer_and_response_order() {
        let fixture = std::fs::read_to_string("data/test_reads/mock.fastq").unwrap();
        let hit = fixture.lines().nth(1).unwrap().as_bytes().to_vec();
        let cases = [
            Vec::new(),
            vec![b'A'; 34],
            vec![b'A'; 35],
            hit,
            b"ACGTNacgt".repeat(40),
            b"ACGT".repeat(12000),
        ];
        for threads in [1, 4] {
            let engine = Engine::load(&EngineConfig {
                db_dir: "data/influenza_test_db".to_string(),
                accessions: false,
                background: false,
                threads,
            })
            .unwrap();
            let mut requests = Vec::new();
            let mut expected = Vec::new();
            let mut sequences = 0;
            let mut bases = 0;
            // Grow, shrink, and clear connection buffers, crossing parallel chunk boundaries.
            for size in [0, 1, 100, 256, 257, 600, 2, 0, 300] {
                let batch: Vec<_> = (0..size)
                    .map(|i| cases[i % cases.len()].as_slice())
                    .collect();
                requests.extend(frame(&batch));
                let results: Vec<_> = batch.iter().map(|seq| engine.query_one(seq)).collect();
                write_response(&mut expected, &results).unwrap();
                sequences += size as u64;
                bases += batch.iter().map(|seq| seq.len() as u64).sum::<u64>();
            }
            let mut output = Vec::new();
            let progress = ServerProgress::new();
            serve_connection(&engine, requests.as_slice(), &mut output, &progress).unwrap();
            assert_eq!(output, expected);
            assert_eq!(progress.totals(), (9, sequences, bases));
        }
    }

    #[test]
    fn request_buffer_reuse_ignores_old_tail_and_checks_sequence_limit() {
        let mut batch = RequestBatch::default();
        let large = frame(&[&[b'A'; 1000]]);
        read_request_into(&mut large.as_slice(), &mut batch).unwrap();
        let ptr = batch.data.as_ptr();
        let small = frame(&[b"ACGT", b""]);
        read_request_into(&mut small.as_slice(), &mut batch).unwrap();
        assert_eq!(batch.data.as_ptr(), ptr);
        assert_eq!(batch.ranges, [0..4, 4..4]);
        assert_eq!(&batch.data[batch.ranges[0].clone()], b"ACGT");
        let mut oversized = 1u32.to_le_bytes().to_vec();
        oversized.extend_from_slice(&(MAX_SEQ_LEN + 1).to_le_bytes());
        assert!(read_request_into(&mut oversized.as_slice(), &mut batch).is_err());
    }

    #[test]
    fn accession_table_and_sorted_counts_roundtrip_across_batches() {
        // Use three overlapping classes and unclassed slots; the bundled index has one class.
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("target/serve-class-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for entry in std::fs::read_dir("data/influenza_test_db").unwrap() {
            let entry = entry.unwrap();
            if entry.file_name().to_string_lossy().starts_with("db.") {
                std::fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
            }
        }
        let slots = test_engine().db().num_minimizers;
        let mut index = crate::database::EqClassAccessions::new_empty(slots);
        let mut hashes = rustc_hash::FxHashMap::default();
        let mut build_map = rustc_hash::FxHashMap::default();
        for slot in 0..slots {
            let ids: &[u32] = match slot % 4 {
                0 => &[0],
                1 => &[1],
                2 => &[0, 1],
                _ => &[],
            };
            index.add_accessions(slot, ids, &mut hashes, &mut build_map);
        }
        index.finalize_for_save(build_map);
        index
            .save(dir.join("db.accession").to_str().unwrap())
            .unwrap();
        let mut registry = crate::database::AccessionRegistry::new();
        registry.get_or_create("reference_A.1");
        registry.get_or_create("reference_β.2");
        registry
            .save(dir.join("db.accessions").to_str().unwrap())
            .unwrap();
        let fixture = std::fs::read_to_string("data/test_reads/mock.fastq").unwrap();
        let hit = fixture.lines().nth(1).unwrap().as_bytes().to_vec();
        let target: Vec<u8> = std::fs::read_to_string("data/influenza_test_db/target.fasta")
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.as_bytes().to_vec())
            .collect();
        for threads in [1, 4] {
            let engine = Engine::load(&EngineConfig {
                db_dir: dir.to_string_lossy().into_owned(),
                accessions: true,
                background: false,
                threads,
            })
            .unwrap();
            let cases = [
                Vec::new(),
                vec![b'A'; engine.k() - 1],
                vec![b'N'; engine.k() + 10],
                hit.clone(),
                hit.repeat(3),
                target.clone(),
            ];
            let sequences: Vec<_> = (0..601)
                .map(|i| cases[i % cases.len()].as_slice())
                .collect();
            let mut request = CLASS_TABLE_TAG.to_le_bytes().to_vec();
            request.extend(frame(&sequences));
            request.extend(frame(&[]));
            request.extend(frame(&[&hit]));
            let mut output = Vec::new();
            let progress = ServerProgress::new();
            serve_connection(&engine, request.as_slice(), &mut output, &progress).unwrap();
            let mut reader = output.as_slice();
            let registry = engine.acc_registry().as_ref().unwrap();
            assert_eq!(read_u32(&mut reader).unwrap() as usize, registry.len());
            for id in 0..registry.len() {
                let len = read_u32(&mut reader).unwrap() as usize;
                assert_eq!(&reader[..len], registry.get_name(id as u32).as_bytes());
                reader = &reader[len..];
            }
            let index = engine.db().accessions.as_ref().unwrap();
            assert_eq!(read_u32(&mut reader).unwrap() as usize, index.num_classes());
            for class in 0..index.num_classes() {
                let n = read_u32(&mut reader).unwrap();
                let ids: Vec<_> = (0..n).map(|_| read_u32(&mut reader).unwrap()).collect();
                let mut expected = Vec::new();
                index.decode_class_into(class as u32, &mut expected);
                assert_eq!(ids, expected);
                assert!(ids.iter().all(|id| (*id as usize) < registry.len()));
            }
            let mut total_classed = 0;
            let mut observed = std::collections::BTreeSet::new();
            for batch in [sequences, vec![], vec![hit.as_slice()]] {
                assert_eq!(read_u32(&mut reader).unwrap() as usize, batch.len());
                for seq in batch {
                    let mut expected = std::collections::BTreeMap::<u32, u32>::new();
                    let n = seq.len().checked_sub(engine.k()).map_or(0, |n| n + 1);
                    // Independent one-window lookups verify repeated minimizer expansion too.
                    for i in 0..n {
                        let mut minimizers = Vec::new();
                        engine
                            .scanner()
                            .scan_into(&seq[i..i + engine.k()], &mut minimizers);
                        if minimizers[0] == u64::MAX {
                            continue;
                        }
                        let mut hits = vec![crate::querydb::Lookup::new(minimizers[0])];
                        engine.db().lookup_batch(&mut hits, &mut Vec::new(), &[]);
                        let h = hits[0];
                        if h.is_hit && h.class_id != u32::MAX && h.class_id != 0 {
                            *expected.entry(h.class_id).or_default() += 1;
                        }
                    }
                    assert_eq!(read_u32(&mut reader).unwrap() as usize, n);
                    assert_eq!(read_u32(&mut reader).unwrap() as usize, expected.len());
                    for (id, count) in expected {
                        observed.insert(id);
                        assert_eq!(read_u32(&mut reader).unwrap(), id);
                        assert_eq!(read_u32(&mut reader).unwrap(), count);
                        total_classed += count;
                    }
                }
            }
            assert!(total_classed > 0);
            assert_eq!(observed.len(), 3);
            assert!(reader.is_empty());
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn table_without_accessions_returns_error_and_session_continues() {
        let engine = test_engine();
        let mut request = CLASS_TABLE_TAG.to_le_bytes().to_vec();
        request.extend(frame(&[b"ACGT"]));
        let mut output = Vec::new();
        serve_connection(
            &engine,
            request.as_slice(),
            &mut output,
            &ServerProgress::new(),
        )
        .unwrap();
        let mut reader = output.as_slice();
        assert_eq!(read_u32(&mut reader).unwrap(), ERROR_TAG);
        let len = read_u32(&mut reader).unwrap() as usize;
        assert_eq!(&reader[..len], b"Class table requires serve -a");
        reader = &reader[len..];
        assert_eq!(read_u32(&mut reader).unwrap(), 1);
        assert_eq!(read_u32(&mut reader).unwrap(), 0);
        assert!(reader.is_empty());
    }

    #[test]
    fn failed_response_flush_is_not_counted_as_completed_work() {
        struct FailingFlush;
        impl Write for FailingFlush {
            fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
                Ok(data.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Err(std::io::Error::new(
                    std::io::ErrorKind::BrokenPipe,
                    "peer closed",
                ))
            }
        }
        let engine = test_engine();
        let bytes = frame(&[b"ACGT"]);
        let progress = ServerProgress::new();
        assert!(serve_connection(&engine, bytes.as_slice(), FailingFlush, &progress).is_err());
        assert_eq!(progress.totals(), (0, 0, 0));
    }
}
