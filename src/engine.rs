//! Resident, thread-safe query engine.
//!
//! [`Engine`] owns a fully loaded [`KmerDatabase`] plus its optional domain bloom filters and
//! accession registry, so the expensive one-time load (MPHF deserialization, accession offset
//! walk, rank-cache rebuild, bloom filter deserialization) is paid once per process instead of
//! once per invocation.
//!
//! `KmerDatabase` has no interior mutability and is already shared by `&`-reference across rayon
//! workers in the CLI path, so it is `Sync` and an `Arc<Engine>` can be cloned freely across
//! connections or threads without copying the database.

use std::sync::Arc;
use std::time::Duration;

use anyhow::Result;
use rayon::prelude::*;

use crate::database::{AccessionRegistry, DomainBloomFilter, Hit, KmerDatabase};
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
    db: KmerDatabase,
    scanner: MinimizerScanner,
    bg_filters: Vec<(DomainBloomFilter, u32)>,
    acc_registry: Option<AccessionRegistry>,
    pool: rayon::ThreadPool,
    load_time: Duration,
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
        let db = KmerDatabase::load(&db_prefix, cfg.accessions)?;

        let mut bg_filters = Vec::new();
        if cfg.background {
            let paths: Vec<(String, u32)> = crate::custom::load_domains(std::path::Path::new(&cfg.db_dir))?
                .into_iter()
                .map(|d| (format!("{}.{}.bloom", db_prefix, d.identity), d.taxid))
                .filter(|(p, _)| std::path::Path::new(p).exists())
                .collect();

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

        let scanner = MinimizerScanner::new(db.k(), db.l(), db.spaced_seed_mask, db.toggle_mask);

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
                    || (Vec::<u64>::new(), Vec::<Hit>::new()),
                    |(minimizer_buf, hit_buf), seq| {
                        hit_buf.clear();
                        self.db.query_into(
                            &self.scanner,
                            seq.as_ref(),
                            minimizer_buf,
                            hit_buf,
                            &self.bg_filters,
                        );
                        hit_buf.iter().map(|h| self.hit_to_taxid(h)).collect()
                    },
                )
                .collect()
        })
    }

    /// Per-k-mer taxids for a single sequence. Runs on the calling thread.
    pub fn query_one(&self, seq: &[u8]) -> Vec<u32> {
        let mut minimizer_buf = Vec::new();
        let mut hit_buf = Vec::new();
        self.db
            .query_into(&self.scanner, seq, &mut minimizer_buf, &mut hit_buf, &self.bg_filters);
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

    pub fn db(&self) -> &KmerDatabase {
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
        self.db.k()
    }

    pub fn l(&self) -> usize {
        self.db.l()
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
        text.lines().nth(1).expect("sequence line").as_bytes().to_vec()
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
        let seq = text.lines().nth(5).expect("second record sequence").as_bytes();
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
}
