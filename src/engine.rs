//! Resident, thread-safe query engine.
//!
//! [`Engine`] owns a fully loaded [`QueryDb`] plus its optional domain bloom filters and
//! accession registry, so the expensive one-time load (MPHF deserialization, accession offset
//! walk, rank-cache rebuild, bloom filter deserialization) is paid once per process instead of
//! once per invocation.
//!
//! `QueryDb` has no interior mutability and is already shared by `&`-reference across rayon
//! workers in the CLI path, so it is `Sync` and an `Arc<Engine>` can be cloned freely across
//! connections or threads without copying the database.

use std::sync::Arc;
use std::time::Duration;

use anyhow::Result;
use rayon::prelude::*;

use crate::database::{AccessionRegistry, DomainBloomFilter};
use crate::querydb::{Lookup as Hit, QueryDb};
use crate::minimizer::MinimizerScanner;

/// How to load a resident database.
#[derive(Clone, Debug)]
pub struct EngineConfig {
    /// Database directory (the one passed to `-d/--db`); the file prefix is `<db_dir>/db`.
    pub db_dir: String,
    /// Load the accession equivalence-class structure (`-a`).
    pub accessions: bool,
    /// Load the per-domain background bloom filters (`-b`).
    pub background: bool,
    /// Worker threads for the engine's own rayon pool.
    pub threads: usize,
}

impl EngineConfig {
    /// Minimal config: database only, one thread per core.
    pub fn new(db_dir: impl Into<String>) -> Self {
        Self {
            db_dir: db_dir.into(),
            accessions: false,
            background: false,
            threads: num_cpus::get(),
        }
    }
}

/// A loaded database plus everything needed to answer queries against it.
pub struct Engine {
    db: QueryDb,
    scanner: MinimizerScanner,
    bg_filters: Vec<(DomainBloomFilter, u32)>,
    acc_registry: Option<AccessionRegistry>,
    pool: rayon::ThreadPool,
    load_time: Duration,
}

// Coarse jobs amortize scheduling and keep scratch allocations across requests.
const QUERY_CHUNK_SIZE: usize = 256;

// Neighboring jobs must not contend over cache lines containing mutable Vec metadata.
#[repr(align(128))]
#[derive(Default)]
pub(crate) struct BatchScratch {
    minimizers: Vec<u64>,
    hits: Vec<Hit>,
    probes: Vec<u32>,
}

impl Engine {
    /// Load a database into memory.
    ///
    /// This performs all of the work that `run_query` used to do on every invocation. The
    /// resulting `Arc` is cheap to clone and safe to share across threads.
    pub fn load(cfg: &EngineConfig) -> Result<Arc<Engine>> {
        let start = std::time::Instant::now();

        // Built before the filters load so the parallel read below has a pool to run on.
        // Deliberately a private pool built with .build(), not rayon's global pool. build_global()
        // can only succeed once per process and silently swallows later failures, which would make
        // a long-lived server ignore its own thread setting and would collide with whatever pool a
        // host process (e.g. Python) has already configured.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(cfg.threads.max(1))
            .thread_name(|i| format!("tdkc-worker-{i}"))
            .build()?;

        let db_prefix = format!("{}/db", cfg.db_dir);
        let mut db = pool.install(|| QueryDb::load(&db_prefix, cfg.accessions))?;

        let mut bg_filters = Vec::new();
        if cfg.background {
            let domains = crate::custom::load_domains(std::path::Path::new(&cfg.db_dir))?;
            let mut paths: Vec<(String, u32)> =
                domains.iter()
                    .map(|d| (format!("{}.{}.bloom", db_prefix, d.identity), d.taxid))
                    .collect();
            if domains.is_empty() {
                for label in ["bacteria", "archaea", "viral", "fungi"] {
                    paths.push((
                        format!("{}.{}.bloom", db_prefix, label),
                        crate::custom::default_domain_taxid(label).unwrap(),
                    ));
                }
            }
            paths.retain(|(path, _)| std::path::Path::new(path).exists());

            // Load the filters concurrently -- they are independent files, and each individual
            // load is itself internally parallel. Results are collected in registry order so
            // background LCA resolution stays deterministic regardless of completion order.
            let loaded: Vec<Option<(DomainBloomFilter, u32)>> = pool.install(|| {
                paths
                    .par_iter()
                    .map(|(path, taxid)| {
                        eprintln!("Loading background filter: {}", path);
                        match DomainBloomFilter::load(path) {
                            Ok(filter) => Some((filter, *taxid)),
                            Err(e) => {
                                eprintln!("Failed to load {}: {}", path, e);
                                None
                            }
                        }
                    })
                    .collect()
            });
            bg_filters.extend(loaded.into_iter().flatten());
        }

        if !bg_filters.is_empty() {
            db.add_domain_nodes(&bg_filters.iter().map(|(_, taxid)| *taxid).collect::<Vec<_>>());
        }

        let acc_registry = if cfg.accessions {
            let acc_path = format!("{}.accessions", db_prefix);
            match AccessionRegistry::load_for_query(&acc_path) {
                Ok(reg) => Some(reg),
                Err(_) => {
                    eprintln!("No accession registry found, skipping accession output");
                    None
                }
            }
        } else {
            None
        };

        let scanner = MinimizerScanner::new(db.k, db.l, db.spaced_seed_mask, db.toggle_mask);

        Ok(Arc::new(Engine {
            db,
            scanner,
            bg_filters,
            acc_registry,
            pool,
            load_time: start.elapsed(),
        }))
    }

    /// Per-k-mer taxids for a batch of sequences, dispatched across the engine's thread pool.
    ///
    /// The returned outer `Vec` is index-aligned with `seqs`. Each inner `Vec` has one entry per
    /// k-mer in sequence order: `0` for a miss, the resolved taxid for a database hit, and the
    /// domain taxid for a background bloom-filter hit. Sequences shorter than `k` yield an empty
    /// vector.
    pub fn query_batch<S: AsRef<[u8]> + Sync>(&self, seqs: &[S]) -> Vec<Vec<u32>> {
        self.pool.install(|| {
            seqs.par_iter()
                .map_init(
                    // Per-worker scratch, reused across every sequence that worker handles.
                    // Allocating these per sequence would dominate at short read lengths.
                    || (Vec::<u64>::new(), Vec::<Hit>::new(), Vec::<u32>::new()),
                    |(minimizer_buf, hit_buf, probes), seq| {
                        hit_buf.clear();
                        self.db.query_into(
                            &self.scanner,
                            seq.as_ref(),
                            minimizer_buf,
                            hit_buf,
                            probes,
                            &self.bg_filters,
                        );
                        hit_buf.iter().map(|h| self.hit_to_taxid(h)).collect()
                    },
                )
                .collect()
        })
    }

    /// Encode the existing per-k-mer wire response directly into a reusable flat buffer.
    /// Each parallel job owns a disjoint output slice; no per-sequence result allocation.
    pub(crate) fn query_batch_wire(
        &self,
        data: &[u8],
        ranges: &[std::ops::Range<usize>],
        output: &mut Vec<u8>,
        scratch: &mut Vec<BatchScratch>,
    ) -> Result<()> {
        let encoded_len = |range: &std::ops::Range<usize>| {
            let kmers = range.len().checked_sub(self.k()).map_or(0, |n| n + 1);
            4 + 4 * kmers
        };
        let total = ranges.iter().try_fold(4usize, |total, range| {
            total
                .checked_add(encoded_len(range))
                .ok_or_else(|| anyhow::anyhow!("Response size overflow"))
        })?;
        output.resize(total, 0);
        output[..4].copy_from_slice(&(ranges.len() as u32).to_le_bytes());
        // Keep long-read batches parallel too: cap a job at roughly 32 KiB of input.
        let bases: usize = ranges.iter().map(|r| r.len()).sum();
        let average_len = bases / ranges.len().max(1);
        // Several jobs per worker balance variable lengths and keep small batches parallel.
        let balanced = ranges.len().div_ceil(self.pool.current_num_threads() * 4);
        let chunk_size = (32768 / average_len.max(1))
            .clamp(1, QUERY_CHUNK_SIZE)
            .min(balanced.max(1));
        let chunks = ranges.len().div_ceil(chunk_size);
        if scratch.len() < chunks {
            scratch.resize_with(chunks, BatchScratch::default);
        }
        let mut remaining = &mut output[4..];
        let mut jobs = Vec::with_capacity(chunks);
        for (ranges, scratch) in ranges.chunks(chunk_size).zip(scratch.iter_mut()) {
            let len = ranges.iter().map(encoded_len).sum();
            let (chunk, rest) = remaining.split_at_mut(len);
            remaining = rest;
            jobs.push((ranges, chunk, scratch));
        }
        let run = |(ranges, output, scratch): (
            &[std::ops::Range<usize>],
            &mut [u8],
            &mut BatchScratch,
        )| {
            let mut remaining = output;
            for range in ranges {
                scratch.hits.clear();
                self.db.query_into(
                    &self.scanner,
                    &data[range.clone()],
                    &mut scratch.minimizers,
                    &mut scratch.hits,
                    &mut scratch.probes,
                    &self.bg_filters,
                );
                let size = 4 + 4 * scratch.hits.len();
                let (record, rest) = remaining.split_at_mut(size);
                remaining = rest;
                record[..4].copy_from_slice(&(scratch.hits.len() as u32).to_le_bytes());
                for (bytes, hit) in record[4..].chunks_exact_mut(4).zip(&scratch.hits) {
                    bytes.copy_from_slice(&self.hit_to_taxid(hit).to_le_bytes());
                }
            }
            debug_assert!(remaining.is_empty());
        };
        // All connections share the configured worker limit, including small requests.
        self.pool.install(|| {
            if jobs.len() <= 1 {
                jobs.into_iter().for_each(run);
            } else {
                jobs.into_par_iter().for_each(run);
            }
        });
        Ok(())
    }

    /// Per-k-mer taxids for a single sequence. Runs on the calling thread.
    pub fn query_one(&self, seq: &[u8]) -> Vec<u32> {
        let mut minimizer_buf = Vec::new();
        let mut hit_buf = Vec::new();
        let mut probes = Vec::new();
        self.db.query_into(
            &self.scanner,
            seq,
            &mut minimizer_buf,
            &mut hit_buf,
            &mut probes,
            &self.bg_filters,
        );
        hit_buf.iter().map(|h| self.hit_to_taxid(h)).collect()
    }

    #[inline]
    fn hit_to_taxid(&self, hit: &Hit) -> u32 {
        if hit.is_hit {
            self.db.true_taxid(hit.taxid_idx)
        } else {
            // 0 when the background filters missed too, or when they were not loaded at all.
            hit.bg_taxid
        }
    }

    /// Run `f` inside the engine's thread pool.
    ///
    /// Used by the CLI path so its existing reader/`par_bridge`/writer pipeline shares this pool
    /// rather than rayon's global one.
    pub fn install<R: Send>(&self, f: impl FnOnce() -> R + Send) -> R {
        self.pool.install(f)
    }

    pub fn db(&self) -> &QueryDb {
        &self.db
    }

    pub fn scanner(&self) -> &MinimizerScanner {
        &self.scanner
    }

    pub fn bg_filters(&self) -> &[(DomainBloomFilter, u32)] {
        &self.bg_filters
    }

    pub fn acc_registry(&self) -> &Option<AccessionRegistry> {
        &self.acc_registry
    }

    pub fn k(&self) -> usize {
        self.db.k
    }

    pub fn l(&self) -> usize {
        self.db.l
    }

    /// How long [`Engine::load`] took — the cost a resident process pays once instead of per run.
    pub fn load_time(&self) -> Duration {
        self.load_time
    }

    pub fn taxid_name(&self, taxid: u32) -> String {
        self.db.get_taxid_name(taxid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_DB: &str = "data/influenza_test_db";

    fn test_engine() -> Arc<Engine> {
        Engine::load(&EngineConfig {
            db_dir: TEST_DB.to_string(),
            accessions: false,
            background: false,
            threads: 2,
        })
        .expect("test database should load")
    }

    /// Reads the first record's sequence out of the bundled mock FASTQ.
    fn first_mock_read() -> Vec<u8> {
        let text = std::fs::read_to_string("data/test_reads/mock.fastq").expect("mock reads");
        text.lines()
            .nth(1)
            .expect("sequence line")
            .as_bytes()
            .to_vec()
    }

    #[test]
    fn emits_one_taxid_per_kmer() {
        let engine = test_engine();
        let seq = first_mock_read();
        let taxids = engine.query_one(&seq);
        assert_eq!(
            taxids.len(),
            seq.len() - engine.k() + 1,
            "expected exactly one entry per k-mer, in sequence order"
        );
    }

    /// Pins the Hit -> taxid flattening against the CLI's own run-length output. The golden
    /// harness reports `C read_1_influenza_mock 114727 50 114727:16` for this read, i.e. 16
    /// consecutive k-mers all hitting taxid 114727.
    #[test]
    fn flattening_matches_cli_hit_pattern() {
        let engine = test_engine();
        let taxids = engine.query_one(&first_mock_read());
        assert_eq!(taxids, vec![114727u32; 16]);
    }

    /// The second mock read is random sequence; the CLI reports `0:16` for it.
    #[test]
    fn misses_flatten_to_zero() {
        let engine = test_engine();
        let text = std::fs::read_to_string("data/test_reads/mock.fastq").expect("mock reads");
        let seq = text
            .lines()
            .nth(5)
            .expect("second record sequence")
            .as_bytes();
        assert_eq!(engine.query_one(seq), vec![0u32; 16]);
    }

    #[test]
    fn short_sequence_yields_no_kmers() {
        let engine = test_engine();
        let short = vec![b'A'; engine.k() - 1];
        assert!(engine.query_one(&short).is_empty());
    }

    #[test]
    fn batch_matches_single_and_preserves_order() {
        let engine = test_engine();
        let seq = first_mock_read();
        let other = vec![b'A'; engine.k() + 10];

        let batch = engine.query_batch(&[seq.clone(), other.clone()]);
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0], engine.query_one(&seq));
        assert_eq!(batch[1], engine.query_one(&other));
    }

    #[test]
    fn repeated_queries_are_stable() {
        // Nothing may leak between requests: the engine carries no per-request state.
        let engine = test_engine();
        let seq = first_mock_read();
        let first = engine.query_batch(&[seq.clone()]);
        let second = engine.query_batch(&[seq]);
        assert_eq!(first, second);
    }

    #[test]
    fn wide_taxids_and_multiple_custom_background_filters() -> Result<()> {
        use crate::database::KmerDatabase;
        use crate::minimizer::{create_spaced_seed_mask, SPACED_PATTERN, TOGGLE_MASK};

        let dir = std::env::temp_dir().join(format!("tdkc-wide-engine-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("genome"))?;
        let result = (|| -> Result<()> {
            let seed = create_spaced_seed_mask(SPACED_PATTERN);
            let scanner = MinimizerScanner::new(35, 31, seed, TOGGLE_MASK);
            let seq = first_mock_read();
            let mut keys = Vec::new();
            scanner.scan_into(&seq, &mut keys);
            keys.retain(|&m| m != u64::MAX);
            keys.sort_unstable();
            keys.dedup();
            let mphf = boomphf::Mphf::new(2.0, &keys);
            let mut fingerprints = vec![0; keys.len()];
            for &key in &keys {
                fingerprints[mphf.hash(&key) as usize] = crate::hash::compute_fingerprint(key);
            }
            let mut taxids: Vec<u32> = (2_000_000..2_000_300).collect();
            taxids[0] = crate::taxonomy::TAXID_ROOT;
            let mut ancestors = vec![0; 300 * 300];
            ancestors[..300].fill(1);
            for i in 0..300 { ancestors[i * 300 + i] = 1; }
            // Five custom domains share this parent, exercising LCA beyond four filters.
            ancestors[302..307].fill(1);
            KmerDatabase::new(
                35, 31, seed, TOGGLE_MASK, keys.len(), mphf, fingerprints,
                vec![299; keys.len()], taxids, ancestors, None, None,
            ).save(&dir.join("db").to_string_lossy())?;
            let mut domains = Vec::new();
            for i in 2..7 {
                let identity = format!("custom_{i}");
                let filter = DomainBloomFilter::new(vec![255], 1, 8);
                std::fs::write(dir.join(format!("db.{identity}.bloom")), bincode::serialize(&filter)?)?;
                domains.push(crate::custom::CustomDomain {
                    identity, label: "custom".into(), taxid: 2_000_000 + i,
                    fasta: String::new(), sha256: String::new(),
                });
            }
            crate::custom::save_domains(&dir, &domains)?;
            let engine = Engine::load(&EngineConfig {
                db_dir: dir.to_string_lossy().into_owned(), accessions: false,
                background: true, threads: 2,
            })?;
            assert_eq!(engine.query_one(&seq), vec![2_000_299; seq.len() - 34]);
            let text = std::fs::read_to_string("data/test_reads/mock.fastq")?;
            assert_eq!(engine.query_one(text.lines().nth(5).unwrap().as_bytes()), vec![2_000_001; 16]);
            let mut interrupted = seq.clone();
            interrupted.extend_from_slice(&[b'N'; 35]);
            interrupted.extend_from_slice(&seq);
            let hits = engine.query_one(&interrupted);
            assert_eq!(&hits[..16], &[2_000_299; 16]);
            assert!(hits[16..seq.len() + 1].iter().all(|&t| t == 0));
            assert_eq!(&hits[hits.len() - 16..], &[2_000_299; 16]);
            assert_eq!(engine.query_batch(&[interrupted.clone()])[0], hits);
            // Legacy upstream databases discover standard filters without a registry.
            std::fs::remove_file(dir.join("genome/domains.tsv"))?;
            for label in ["archaea", "viral"] {
                let filter = DomainBloomFilter::new(vec![255], 1, 8);
                std::fs::write(dir.join(format!("db.{label}.bloom")), bincode::serialize(&filter)?)?;
            }
            let legacy = Engine::load(&EngineConfig {
                db_dir: dir.to_string_lossy().into_owned(), accessions: false,
                background: true, threads: 2,
            })?;
            assert_eq!(legacy.query_one(text.lines().nth(5).unwrap().as_bytes()), vec![1; 16]);
            Ok(())
        })();
        std::fs::remove_dir_all(dir)?;
        result
    }
}
