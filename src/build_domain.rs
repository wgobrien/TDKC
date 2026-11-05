use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::Arc;
use std::time::Instant;
use anyhow::{Context, Result};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use crate::database::DomainBloomFilter;
use crate::minimizer::MinimizerScanner;
use crate::taxonomy::{TaxonomyTree, TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_FUNGI, TAXID_VIRAL};
use crate::utils::{extract_accession, init_thread_pool};
use crate::hash::{BuildFmix64Hasher, bloom_hash_pair};

const RAW_FPR: f64 = 0.0001;

const BACTERIA: usize = 0;
const OTHER: usize = 4;
const GROUP_NAMES: [&str; 5] = ["bacteria", "archaea", "viral", "fungi", "other"];

struct AtomicBloom {
    data: Vec<AtomicU8>,
    num_hashes: u32,
    num_bits: u64,
}

impl AtomicBloom {
    fn new(expected_items: usize, fpr: f64) -> Self {
        let n = (expected_items.max(1)) as f64;
        let ln2 = std::f64::consts::LN_2;
        let m = (-n * fpr.ln() / (ln2 * ln2)).ceil() as u64;
        let k = ((m as f64 / n) * ln2).round() as u32;
        let k = k.max(1);

        let num_bytes = ((m + 7) / 8) as usize;
        let mut data = Vec::with_capacity(num_bytes);
        for _ in 0..num_bytes {
            data.push(AtomicU8::new(0));
        }

        eprintln!("  {} bits ({:.2} GiB), {} hashes for {} items @ {}% FPR", m, num_bytes as f64 / 1_073_741_824.0, k, expected_items, fpr * 100.0);

        Self { data, num_hashes: k, num_bits: m }
    }

    #[inline]
    fn insert(&self, item: u64) {
        let (h1, h2) = bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            self.data[byte_idx].fetch_or(bit_mask, Ordering::Relaxed);
        }
    }

    #[inline]
    fn contains(&self, item: u64) -> bool {
        let (h1, h2) = bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            if self.data[byte_idx].load(Ordering::Relaxed) & bit_mask == 0 {
                return false;
            }
        }
        true
    }

    fn into_filter(self) -> DomainBloomFilter {
        let (k, m) = (self.num_hashes, self.num_bits);
        DomainBloomFilter::new(self.data.into_iter().map(|a| a.into_inner()).collect(), k, m)
    }
}

pub struct BuildDomainConfig {
    pub db_dir: String,
    pub threads: usize,
    pub fpr: f64,
    pub bacteria: bool,
    pub archaea: bool,
    pub viral: bool,
    pub fungi: bool,
}

struct RefFile {
    path: String,
    fixed: Option<usize>,
}

type Hll = HyperLogLogPlus<u64, BuildFmix64Hasher>;

fn new_hll() -> Hll {
    HyperLogLogPlus::<u64, BuildFmix64Hasher>::new(16, BuildFmix64Hasher::default()).unwrap()
}

pub fn run_build_domain(config: BuildDomainConfig) -> Result<()> {
    init_thread_pool(config.threads);

    let total_start = Instant::now();
    let db_path = std::path::PathBuf::from(&config.db_dir);
    let db_prefix = db_path.join("db").to_string_lossy().into_owned();
    let genome_dir = db_path.join("genome");

    let (k, l, spaced_seed_mask, toggle_mask) = load_db_meta(&db_prefix)?;
    let scanner = MinimizerScanner::new(k, l, spaced_seed_mask, toggle_mask);

    let none_specified = !config.bacteria && !config.archaea && !config.viral && !config.fungi;
    let mut registered = crate::custom::load_domains(&db_path)?;
    // Older upstream databases have no registry.
    if registered.is_empty() {
        for name in &GROUP_NAMES[..OTHER] {
            let fasta = format!("{}.fna", name);
            if genome_dir.join(&fasta).exists() {
                registered.push(crate::custom::CustomDomain {
                    identity: name.to_string(), label: name.to_string(),
                    taxid: crate::custom::default_domain_taxid(name).unwrap(),
                    fasta, sha256: String::new(),
                });
            }
        }
    }
    registered.retain(|domain| {
        let path = genome_dir.join(&domain.fasta);
        if path.exists() { return true; }
        eprintln!("Can't find registered domain file {}", path.display());
        false
    });
    let mut group_taxids = vec![TAXID_BACTERIA, TAXID_ARCHAEA, TAXID_VIRAL, TAXID_FUNGI, 0];
    let mut group_names: Vec<String> = GROUP_NAMES.iter().map(|n| n.to_string()).collect();
    let mut registered_groups = Vec::new();
    let mut requested = Vec::new();
    for domain in &registered {
        let g = match group_taxids.iter().position(|&t| t == domain.taxid) {
            Some(g) => g,
            None => {
                group_taxids.push(domain.taxid);
                group_names.push(domain.label.clone());
                group_taxids.len() - 1
            }
        };
        registered_groups.push(g);
        let selected = none_specified || match domain.label.as_str() {
            "bacteria" => config.bacteria,
            "archaea" => config.archaea,
            "viral" | "viruses" => config.viral,
            "fungi" => config.fungi,
            _ => false,
        };
        if selected && genome_dir.join(&domain.fasta).exists() {
            requested.push((g, domain));
        }
    }
    if requested.is_empty() {
        eprintln!("No registered domains selected.");
        return Ok(());
    }
    let ngroups = group_taxids.len();
    let mut ref_files = Vec::new();
    let mut names: Vec<String> = std::fs::read_dir(&genome_dir)
        .with_context(|| format!("Can't read {}", genome_dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".fna") && (!n.starts_with('_') || n.starts_with("_target.")))
        .collect();
    names.extend(registered.iter().map(|d| d.fasta.clone()));
    names.sort();
    names.dedup();
    for name in names {
        let fixed = registered.iter().position(|d| d.fasta == name)
            .map(|i| registered_groups[i])
            .or_else(|| (0..OTHER).find(|&g| GROUP_NAMES[g] == name.trim_end_matches(".fna")));
        ref_files.push(RefFile { path: genome_dir.join(&name).to_string_lossy().into_owned(), fixed });
    }
    eprintln!("Reference files:");
    for f in &ref_files {
        eprintln!("  {}{}", f.path, f.fixed.map(|g| format!("  (all records = {})", group_names[g])).unwrap_or_default());
    }

    let t = Instant::now();
    let taxonomy = TaxonomyTree::load(&db_path.join("taxonomy/nodes.dmp").to_string_lossy())?;
    let prelim = crate::prep::load_prelim_map(&db_path.join("prelim_map.txt").to_string_lossy())?;
    let group_of = Arc::new(GroupResolver::new(taxonomy, prelim, group_taxids));
    eprintln!("  [taxonomy] {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("\nPass A: estimating unique minimizers per domain group over the reference...");
    let t = Instant::now();
    let rec_counts: Vec<AtomicU64> = (0..ngroups).map(|_| AtomicU64::new(0)).collect();
    let mut hlls = stream_reference(&ref_files, &group_of, &rec_counts, |batch| {
        let mut local: Vec<Hll> = (0..ngroups).map(|_| new_hll()).collect();
        let mut mins = Vec::new();
        for (g, seq) in batch {
            scanner.scan_into(seq, &mut mins);
            for &m in &mins {
                if m != u64::MAX {
                    local[*g].insert(&m);
                }
            }
        }
        local
    })?;
    let mut raw_counts = vec![0usize; ngroups];
    for g in 0..ngroups {
        raw_counts[g] = hlls[g].count().round() as usize;
        eprintln!("  {:9} {:>9} records  ~{} unique minimizers", group_names[g], rec_counts[g].load(Ordering::Relaxed), raw_counts[g]);
    }
    drop(hlls);
    eprintln!("  [pass A] {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("\nPass B: building raw per-group filters @ {}% FPR...", RAW_FPR * 100.0);
    let t = Instant::now();
    let mut raw: Vec<Option<AtomicBloom>> = (0..ngroups)
        .map(|g| {
            if raw_counts[g] == 0 {
                None
            } else {
                eprintln!(" {}:", group_names[g]);
                Some(AtomicBloom::new(raw_counts[g], RAW_FPR))
            }
        })
        .collect();
    {
        let raw_ref = &raw;
        let dummy: Vec<AtomicU64> = (0..ngroups).map(|_| AtomicU64::new(0)).collect();
        stream_reference(&ref_files, &group_of, &dummy, |batch| {
            let mut mins = Vec::new();
            for (g, seq) in batch {
                if let Some(b) = &raw_ref[*g] {
                    scanner.scan_into(seq, &mut mins);
                    for &m in &mins {
                        if m != u64::MAX {
                            b.insert(m);
                        }
                    }
                }
            }
            Vec::new()
        })?;
    }
    eprintln!("  [pass B] {:.1}s", t.elapsed().as_secs_f64());

    // bacteria last so its raw filter, the largest, is freed before its output filter is built
    requested.sort_by_key(|&(g, _)| (g == BACTERIA, g));
    for &(d, domain) in &requested {
        if d == BACTERIA {
            raw[BACTERIA] = None;
        }
        build_distilled(d, domain, &raw, &scanner, &genome_dir.join(&domain.fasta).to_string_lossy(), &db_prefix, config.fpr)?;
    }

    eprintln!("\nDone. It took {:.2}s", total_start.elapsed().as_secs_f64());
    Ok(())
}

fn build_distilled(
    d: usize,
    domain: &crate::custom::CustomDomain,
    raw: &[Option<AtomicBloom>],
    scanner: &MinimizerScanner,
    path: &str,
    db_prefix: &str,
    fpr: f64,
) -> Result<()> {
    let name = &domain.identity;
    eprintln!("\nDomain {}:", name);
    let others: Vec<&AtomicBloom> = (0..raw.len()).filter(|&g| g != d).filter_map(|g| raw[g].as_ref()).collect();
    let exclusive = |m: u64| !others.iter().any(|b| b.contains(m));

    let t = Instant::now();
    let total = AtomicU64::new(0);
    let kept = AtomicU64::new(0);
    let hll = stream_file(&path, |batch| {
        let mut h = new_hll();
        let mut mins = Vec::new();
        let (mut tot, mut kp) = (0u64, 0u64);
        for seq in batch {
            scanner.scan_into(seq, &mut mins);
            for &m in &mins {
                if m == u64::MAX {
                    continue;
                }
                tot += 1;
                if exclusive(m) {
                    kp += 1;
                    h.insert(&m);
                }
            }
        }
        total.fetch_add(tot, Ordering::Relaxed);
        kept.fetch_add(kp, Ordering::Relaxed);
        h
    })?;
    let mut hll = hll;
    let unique = hll.count().round() as usize;
    let (tot, kp) = (total.load(Ordering::Relaxed), kept.load(Ordering::Relaxed));
    eprintln!("  Pass C: {} minimizers (incl. duplicates), {} domain-exclusive ({:.4}% removed), ~{} unique exclusive  [{:.1}s]",
        tot, kp, 100.0 * (tot - kp) as f64 / tot.max(1) as f64, unique, t.elapsed().as_secs_f64());
    if unique == 0 {
        eprintln!("WARNING: no domain-exclusive minimizers for {}", name);
    }

    let t = Instant::now();
    let bloom = AtomicBloom::new(unique, fpr);
    {
        let bloom_ref = &bloom;
        stream_file(&path, |batch| {
            let mut mins = Vec::new();
            for seq in batch {
                scanner.scan_into(seq, &mut mins);
                for &m in &mins {
                    if m != u64::MAX && exclusive(m) {
                        bloom_ref.insert(m);
                    }
                }
            }
            new_hll()
        })?;
    }
    eprintln!("  Pass D: inserted [{:.1}s]", t.elapsed().as_secs_f64());

    let out_path = format!("{}.{}.bloom", db_prefix, name);
    let f = File::create(&out_path)?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, f);
    bincode::serialize_into(&mut writer, &bloom.into_filter())?;
    writer.flush()?;
    eprintln!("  Saved {}", out_path);
    Ok(())
}

struct GroupResolver {
    taxonomy: TaxonomyTree,
    prelim: FxHashMap<String, u32>,
    group_taxids: Vec<u32>,
}

impl GroupResolver {
    fn new(taxonomy: TaxonomyTree, prelim: FxHashMap<String, u32>, group_taxids: Vec<u32>) -> Self {
        Self { taxonomy, prelim, group_taxids }
    }

    fn group_of_taxid(&self, taxid: u32) -> usize {
        for t in self.taxonomy.lineage_path(taxid).into_iter().rev() {
            if let Some(g) = self.group_taxids.iter().position(|&candidate| candidate == t) {
                return g;
            }
        }
        OTHER
    }

    fn group_of_header(&self, id: &[u8]) -> usize {
        let acc = extract_accession(id);
        let taxid = self.prelim.get(acc).copied().or_else(|| {
            acc.strip_prefix("kraken:taxid|").and_then(|r| r.split('|').next()).and_then(|s| s.parse().ok())
        });
        match taxid {
            Some(t) => self.group_of_taxid(t),
            None => OTHER,
        }
    }
}

fn stream_reference<F>(files: &[RefFile], groups: &Arc<GroupResolver>, rec_counts: &[AtomicU64], work: F) -> Result<Vec<Hll>>
where
    F: Fn(&[(usize, Vec<u8>)]) -> Vec<Hll> + Sync,
{
    let ngroups = groups.group_taxids.len();
    let (tx, rx) = crossbeam_channel::bounded::<Vec<(usize, Vec<u8>)>>(32);
    let jobs: Vec<(String, Option<usize>)> = files.iter().map(|f| (f.path.clone(), f.fixed)).collect();
    let groups = Arc::clone(groups);
    let counts: Arc<Vec<AtomicU64>> = Arc::new((0..ngroups).map(|_| AtomicU64::new(0)).collect());
    let counts_reader = Arc::clone(&counts);
    let reader = std::thread::spawn(move || -> Result<()> {
        for (path, fixed) in jobs {
            let mut r = needletail::parse_fastx_file(&path).map_err(|e| anyhow::anyhow!("Cannot open FASTA {}: {}", path, e))?;
            let mut batch = Vec::with_capacity(1000);
            let mut bytes = 0usize;
            while let Some(rec) = r.next() {
                let rec = rec.map_err(|e| anyhow::anyhow!("FASTA error in {}: {}", path, e))?;
                let g = fixed.unwrap_or_else(|| groups.group_of_header(rec.id()));
                counts_reader[g].fetch_add(1, Ordering::Relaxed);
                let seq = rec.seq().into_owned();
                bytes += seq.len();
                batch.push((g, seq));
                if bytes >= 8 * 1024 * 1024 {
                    if tx.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).is_err() {
                        return Ok(());
                    }
                    bytes = 0;
                }
            }
            if !batch.is_empty() {
                let _ = tx.send(batch);
            }
        }
        Ok(())
    });

    let merged = rx
        .into_iter()
        .par_bridge()
        .map(|batch| work(&batch))
        .reduce(Vec::new, merge_hll_vecs);

    match reader.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }
    for g in 0..ngroups {
        rec_counts[g].fetch_add(counts[g].load(Ordering::Relaxed), Ordering::Relaxed);
    }
    if merged.is_empty() {
        return Ok((0..ngroups).map(|_| new_hll()).collect());
    }
    Ok(merged)
}

fn merge_hll_vecs(mut a: Vec<Hll>, b: Vec<Hll>) -> Vec<Hll> {
    if a.is_empty() {
        return b;
    }
    for (x, y) in a.iter_mut().zip(b.iter()) {
        x.merge(y).unwrap();
    }
    a
}

fn stream_file<F>(path: &str, work: F) -> Result<Hll>
where
    F: Fn(&[Vec<u8>]) -> Hll + Sync,
{
    let (tx, rx) = crossbeam_channel::bounded::<Vec<Vec<u8>>>(32);
    let p = path.to_string();
    let reader = std::thread::spawn(move || -> Result<()> {
        let mut r = needletail::parse_fastx_file(&p).map_err(|e| anyhow::anyhow!("Cannot open FASTA {}: {}", p, e))?;
        let mut batch = Vec::with_capacity(1000);
        let mut bytes = 0usize;
        while let Some(rec) = r.next() {
            let rec = rec.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let seq = rec.seq().into_owned();
            bytes += seq.len();
            batch.push(seq);
            if bytes >= 8 * 1024 * 1024 {
                if tx.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).is_err() {
                    return Ok(());
                }
                bytes = 0;
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        Ok(())
    });
    let h = rx
        .into_iter()
        .par_bridge()
        .map(|batch| work(&batch))
        .reduce(new_hll, |mut a, b| {
            a.merge(&b).unwrap();
            a
        });
    match reader.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }
    Ok(h)
}

fn load_db_meta(prefix: &str) -> Result<(usize, usize, u64, u64)> {
    let mut f = File::open(format!("{}.meta", prefix)).context("Cannot open .meta file to read scanner parameters")?;
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    f.read_exact(&mut buf4)?; let k = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf4)?; let l = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf8)?; let spaced_seed_mask = u64::from_le_bytes(buf8);
    f.read_exact(&mut buf8)?; let toggle_mask = u64::from_le_bytes(buf8);

    Ok((k, l, spaced_seed_mask, toggle_mask))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distills_registered_custom_domains_against_each_other() -> Result<()> {
        let dir = std::env::temp_dir().join(format!("tdkc-domain-merge-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("genome"))?;
        std::fs::create_dir_all(dir.join("taxonomy"))?;
        let result = (|| -> Result<()> {
            std::fs::copy("data/influenza_test_db/db.meta", dir.join("db.meta"))?;
            std::fs::write(dir.join("taxonomy/nodes.dmp"), "1\t|\t1\t|\tno rank\t|\n")?;
            std::fs::write(dir.join("prelim_map.txt"), "")?;
            let (k, l, seed, toggle) = load_db_meta(&dir.join("db").to_string_lossy())?;
            let scanner = MinimizerScanner::new(k, l, seed, toggle);
            let mut state = 42u64;
            let mut sequence = || -> Vec<u8> {
                (0..k).map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize]
                }).collect()
            };
            let shared = sequence();
            let unique = [sequence(), sequence(), sequence()];
            let mut domains = Vec::new();
            for (i, (label, taxid)) in [("bacteria", 2), ("viral", 10239), ("human", 9606)].into_iter().enumerate() {
                let identity = format!("{label}_custom");
                let fasta = format!("_domain.{identity}.fna");
                std::fs::write(dir.join("genome").join(&fasta), format!(
                    ">shared\n{}\n>unique\n{}\n",
                    String::from_utf8_lossy(&shared), String::from_utf8_lossy(&unique[i]),
                ))?;
                domains.push(crate::custom::CustomDomain {
                    identity, label: label.into(), taxid, fasta, sha256: String::new(),
                });
            }
            crate::custom::save_domains(&dir, &domains)?;
            run_build_domain(BuildDomainConfig {
                db_dir: dir.to_string_lossy().into_owned(), threads: 2, fpr: 1e-8,
                bacteria: false, archaea: false, viral: false, fungi: false,
            })?;
            let mut mins = Vec::new();
            for (i, domain) in domains.iter().enumerate() {
                let filter = DomainBloomFilter::load(&dir.join(format!("db.{}.bloom", domain.identity)).to_string_lossy())?;
                scanner.scan_into(&unique[i], &mut mins);
                assert!(mins.iter().any(|&m| m != u64::MAX));
                assert!(mins.iter().all(|&m| m == u64::MAX || filter.contains(m)));
                scanner.scan_into(&shared, &mut mins);
                assert!(mins.iter().all(|&m| m == u64::MAX || !filter.contains(m)));
            }
            Ok(())
        })();
        std::fs::remove_dir_all(dir)?;
        result
    }
}
