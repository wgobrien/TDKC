mod download;
mod pipeline;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use anyhow::{Context, Result};
use crate::taxonomy::{load_target_taxids, TargetTaxIDManager, TaxonomyTree};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use download::*;
use pipeline::*;

pub struct PrepConfig {
    pub domains: Option<String>,
    pub targets_file: String,
    pub db_dir: String,
    pub backend: String,
    pub concurrent_downloads: usize,
    pub in_flight_chunks: usize,
    pub threads: usize,
    pub custom_fasta: Option<String>,
    pub custom_domain: bool,
    pub domain_taxid: Option<u32>,
    pub no_mask: bool,
}

fn ensure_dustmasker_available() -> Result<()> {
    match std::process::Command::new("dustmasker")
        .arg("-version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
    {
        Ok(_) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(anyhow::anyhow!(
            "`dustmasker` not found on PATH. Install NCBI BLAST+ (e.g. `brew install blast` \
             or `conda install -c bioconda blast`), or re-run with --no-mask to skip masking."
        )),
        Err(e) => Err(e).context("failed to invoke dustmasker"),
    }
}

pub async fn run_prep(cfg: PrepConfig) -> Result<()> {
    if !cfg.no_mask {
        ensure_dustmasker_available()?;
    }

    let out_dir = PathBuf::from(&cfg.db_dir);
    let genome_dir = out_dir.join("genome");

    tokio::fs::create_dir_all(&genome_dir)
        .await
        .context("Failed to create output genome directory")?;

    let targets = load_target_taxids(&cfg.targets_file)?;
    let stored_targets = out_dir.join("targets.txt");
    if stored_targets.exists() {
        let existing = load_target_taxids(stored_targets.to_string_lossy().as_ref())?;
        if existing != targets {
            anyhow::bail!("The supplied targets do not match this database's existing targets.txt");
        }
    } else {
        tokio::fs::copy(&cfg.targets_file, &stored_targets)
            .await.context("Failed to copy targets file into db dir")?;
    }

    let backend = make_backend(BackendKind::Http)?;
    let tax_dir = out_dir.join("taxonomy");
    let taxonomy_existed = tax_dir.exists();
    let tax_paths = if cfg.custom_fasta.is_some() && taxonomy_existed {
        crate::custom::validate_existing_taxonomy(&tax_dir)?;
        download::TaxdumpPaths { nodes_dmp: tax_dir.join("nodes.dmp"), names_dmp: tax_dir.join("names.dmp") }
    } else {
        download_taxdump(&backend, &out_dir).await?
    };
    if cfg.custom_fasta.is_some() && !taxonomy_existed {
        prefetch_acc2taxid(&backend, &tax_dir).await?;
    }


    let nodes_path = tax_paths.nodes_dmp.to_string_lossy();
    let taxonomy = TaxonomyTree::load(nodes_path.as_ref())?;
    
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);
    
    let relevant_taxids = Arc::new(taxid_manager.all_relevant_taxids());
    
    eprintln!("Found {} target taxids in the tree.", relevant_taxids.len());

    let mut domain_list: Vec<String> = if cfg.custom_fasta.is_some() {
        Vec::new()
    } else {
        cfg.domains.as_deref().unwrap_or("bacteria,viral,archaea,human,plasmid").split(',')
            .map(|s| s.trim().to_lowercase())
            .filter(|s| !s.is_empty() && s != "none")
            .collect()
    };

    if cfg.custom_fasta.is_none() && !domain_list.iter().any(|d| d == "univec" || d == "univec_core") {
        domain_list.push("univec_core".to_string());
    }

    if domain_list.iter().any(|d| d == "plasmid") {
        prefetch_acc2taxid(&backend, &tax_dir).await?;
    }

    let mut total_stats = PipelineStats::default();

    for domain in &domain_list {
        let done_marker = genome_dir.join(format!(".{}.done", domain));
        if done_marker.exists() {
            eprintln!("\nDomain '{}' already exists, skipping.", domain);
            continue;
        }

        eprintln!("\nProcessing Domain: {}", domain);

        let mut sources = Vec::new();
        let mut domain_custom_map: Option<Arc<FxHashMap<String, u32>>> = None;

        //univec is always downloaded
        if domain == "univec" || domain == "univec_core" {
            let meta = resolve_univec(domain);

            sources.push(GenomeSource {
                fetch_path: meta.fetch_path,
                gzipped: false,
                taxid: meta.taxid,
                assembly_accession: meta.filename,
            });

        } else if domain == "plasmid" {
            // plasmid has no assembly_summary.txt; use accession2taxid
            let cache_dir = genome_dir.join(".plasmid_cache");
            tokio::fs::create_dir_all(&cache_dir).await?;

            let files = list_refseq_catalog_files(&backend, "plasmid").await?;

            let mut wanted = std::collections::HashSet::new();
            for fname in files {
                let remote = format!("/genomes/refseq/plasmid/{}", fname);
                let local = cache_dir.join(fname.trim_end_matches(".gz"));
                let accs =
                    fetch_and_decompress_to(&backend, &remote, &local).await?;
                for a in accs {
                    wanted.insert(a);
                }
                sources.push(GenomeSource {
                    fetch_path: local.to_string_lossy().into_owned(),
                    gzipped: false,
                    taxid: 0, 
                    assembly_accession: fname,
                });
            }

            let map = fetch_acc2taxid_filtered(&backend, &wanted, &tax_dir, cfg.threads).await?;
            domain_custom_map = Some(Arc::new(map));

        } else {
            // Standard RefSeq 
            let entries = fetch_assembly_summary(&backend, domain).await?;
            
            for entry in entries {
                match entry.fna_gz_path() {
                    Ok(path) => {
                        sources.push(GenomeSource {
                            fetch_path: path,
                            gzipped: true, // RefSeq is always gzipped
                            taxid: entry.taxid,
                            assembly_accession: entry.assembly_accession,
                        });
                    }
                    Err(e) => eprintln!("Warning: Skipping genome due to path error: {}", e),
                }
            }
        }

        let (shared_writers, shared_handles) = SharedWriters::new(&genome_dir, domain)?;

        let pipeline_cfg = PipelineConfig {
            domain_name: domain.clone(),
            source_label: format!("refseq:{}", domain),
            output_genome_dir: genome_dir.clone(),
            backend: Arc::clone(&backend),
            concurrent_downloads: cfg.concurrent_downloads,
            max_in_flight_chunks: cfg.in_flight_chunks,
            batch_threshold_bytes: 64 * 1024 * 1024,
            relevant_taxids: Arc::clone(&relevant_taxids),
            shared: Arc::clone(&shared_writers),
            custom_map: domain_custom_map,
            background_output: BackgroundOutput::Challenger,
            no_mask: cfg.no_mask,
            // plasmid sequences with no accession2taxid hit are dropped
            drop_unmapped: domain == "plasmid",
        };

        let stats = pipeline::run_pipeline(pipeline_cfg, sources).await?;
        finalize_shared(shared_writers, shared_handles).await?;

        if domain == "plasmid" {
            let _ = tokio::fs::remove_dir_all(genome_dir.join(".plasmid_cache")).await;
        }

        total_stats.genomes_processed += stats.genomes_processed;
        total_stats.records_processed += stats.records_processed;
        total_stats.records_in_target += stats.records_in_target;
        total_stats.chunks_masked += stats.chunks_masked;

        tokio::fs::write(&done_marker, b"")
            .await
            .with_context(|| format!("write done marker {}", done_marker.display()))?;
    }

    if cfg.custom_fasta.is_none() {
        let mut registered = crate::custom::load_domains(&out_dir)?;
        let mut changed = false;
        for domain in &domain_list {
            let fasta = format!("{}.fna", domain);
            if !genome_dir.join(&fasta).is_file() {
                continue;
            }
            let taxid = crate::custom::default_domain_taxid(domain)
                .or_else(|| crate::custom::resolve_taxid(domain, None, &tax_paths.names_dmp).ok())
                .ok_or_else(|| anyhow::anyhow!("Cannot infer reporting taxid for downloaded domain '{}'", domain))?;
            if !taxonomy.contains(taxid) {
                anyhow::bail!("Reporting taxid {} for downloaded domain '{}' is not present in nodes.dmp", taxid, domain);
            }
            if let Some(existing) = registered.iter().find(|entry| entry.identity == *domain) {
                if existing.taxid != taxid || existing.fasta != fasta {
                    anyhow::bail!("Domain registry entry '{}' conflicts with the prepared domain", domain);
                }
            } else {
                registered.push(crate::custom::CustomDomain {
                    identity: domain.clone(),
                    label: domain.clone(),
                    taxid,
                    fasta,
                    sha256: "-".to_string(),
                });
                changed = true;
            }
        }
        if changed {
            registered.sort_by(|a, b| a.identity.cmp(&b.identity));
            crate::custom::save_domains(&out_dir, &registered)?;
        }
    }

    if let Some(custom_path) = &cfg.custom_fasta {
        let custom_path_obj = Path::new(custom_path);
        eprintln!("\nInspecting custom FASTA: {}", custom_path);
        let inspect_started = Instant::now();
        let inspect_path = custom_path_obj.to_path_buf();
        let info = tokio::task::spawn_blocking(move || crate::custom::inspect(&inspect_path))
            .await
            .context("custom FASTA inspection task panicked")??;
        eprintln!(
            "Inspected {} accessions in {:.2?}; identity is {}.",
            info.accessions.len(),
            inspect_started.elapsed(),
            info.identity
        );
        if info.unparseable_headers > 0 {
            eprintln!(
                "WARNING: {} custom FASTA headers have no valid accession token and will be ignored.",
                info.unparseable_headers
            );
        }
        let done_marker = genome_dir.join(format!(".custom.{}.done", info.identity));
        let partition_marker = genome_dir.join(format!(".custom.{}.partitioned", info.identity));
        let custom_domain_path = genome_dir.join(format!("_domain.{}.fna", info.identity));
        let manifest_shard = genome_dir.join(format!("_manifest.{}.tsv", info.identity));
        if manifest_shard.exists() {
            let manifest = std::fs::read_to_string(&manifest_shard)?;
            if let Some(previous_hash) = manifest.lines().last().and_then(|line| line.split('\t').nth(4)) {
                if previous_hash != info.sha256 {
                    anyhow::bail!("Short hash collision for custom identity {}", info.identity);
                }
            }
        }
        let mut registered = crate::custom::load_domains(&out_dir)?;
        if let Some(existing) = registered.iter().find(|d| d.identity == info.identity) {
            if existing.sha256 != info.sha256 {
                anyhow::bail!("Short hash collision for custom identity {}", info.identity);
            }
        }

        let requested_domain_taxid = if cfg.custom_domain {
            let taxid = crate::custom::resolve_taxid(&info.label, cfg.domain_taxid, &tax_paths.names_dmp)?;
            if !taxonomy.contains(taxid) { anyhow::bail!("Domain taxid {} is not present in nodes.dmp", taxid); }
            Some(taxid)
        } else { None };

        let needs_domain_partition = cfg.custom_domain
            && partition_marker.exists()
            && !custom_domain_path.exists();
        if done_marker.exists() && !needs_domain_partition {
            eprintln!("\nCustom FASTA {} already added, skipping sequence processing.", info.identity);
        } else {
            eprintln!("\nProcessing custom FASTA: {} as {}", custom_path, info.identity);
            let custom_map = fetch_acc2taxid_filtered(
                &backend,
                &info.accessions,
                &tax_dir,
                cfg.threads,
            )
            .await?;
            let mut unresolved: Vec<&String> = info
                .accessions
                .iter()
                .filter(|accession| !custom_map.contains_key(*accession))
                .collect();
            unresolved.sort_unstable();
            if !unresolved.is_empty() {
                let examples = unresolved
                    .iter()
                    .take(5)
                    .map(|accession| accession.as_str())
                    .collect::<Vec<_>>()
                    .join(", ");
                eprintln!(
                    "WARNING: {} custom accessions could not be mapped and will be ignored. Examples: {}",
                    unresolved.len(),
                    examples
                );
            }
            let custom_source = GenomeSource {
                fetch_path: custom_path.clone(),
                gzipped: custom_path.to_ascii_lowercase().ends_with(".gz"),
                taxid: 0,
                assembly_accession: info.identity.clone(),
            };
            let (shared_writers, shared_handles) = SharedWriters::new(&genome_dir, &info.identity)?;
            let pipeline_cfg = PipelineConfig {
                domain_name: info.identity.clone(),
                source_label: "custom_file".to_string(),
                output_genome_dir: genome_dir.clone(),
                backend: Arc::clone(&backend),
                concurrent_downloads: 1,
                max_in_flight_chunks: cfg.in_flight_chunks,
                batch_threshold_bytes: 64 * 1024 * 1024,
                relevant_taxids: Arc::clone(&relevant_taxids),
                shared: Arc::clone(&shared_writers),
                custom_map: Some(Arc::new(custom_map)),
                background_output: if cfg.custom_domain {
                    BackgroundOutput::Domain
                } else {
                    BackgroundOutput::Discard
                },
                no_mask: cfg.no_mask,
                drop_unmapped: true,
            };
            let stats = run_pipeline(pipeline_cfg, vec![custom_source]).await?;
            finalize_shared(shared_writers, shared_handles).await?;
            total_stats.genomes_processed += stats.genomes_processed;
            total_stats.records_processed += stats.records_processed;
            total_stats.records_in_target += stats.records_in_target;
            total_stats.records_in_domain += stats.records_in_domain;
            total_stats.records_ignored_unmapped += stats.records_ignored_unmapped;
            total_stats.chunks_masked += stats.chunks_masked;
            if stats.records_in_target == 0 {
                let _ = tokio::fs::remove_file(genome_dir.join(format!("_target.{}.fna", info.identity))).await;
            }
            tokio::fs::write(&done_marker, b"").await
                .with_context(|| format!("write done marker {}", done_marker.display()))?;
            tokio::fs::write(&partition_marker, b"").await
                .with_context(|| format!("write partition marker {}", partition_marker.display()))?;
        }

        if let Some(taxid) = requested_domain_taxid {
            if let Some(existing) = registered.iter().find(|d| d.identity == info.identity) {
                if existing.taxid != taxid { anyhow::bail!("Custom domain {} is already registered with taxid {}", info.identity, existing.taxid); }
            } else {
                let domain_fasta = if partition_marker.exists() {
                    format!("_domain.{}.fna", info.identity)
                } else {
                    // Preserve the layout of custom entries prepared before partitioning.
                    format!("{}.fna", info.identity)
                };
                registered.push(crate::custom::CustomDomain {
                    identity: info.identity.clone(), label: info.label.clone(), taxid,
                    fasta: domain_fasta, sha256: info.sha256.clone(),
                });
                registered.sort_by(|a, b| a.identity.cmp(&b.identity));
                crate::custom::save_domains(&out_dir, &registered)?;
            }
        }
    }

    finalize_outputs(&out_dir, &genome_dir).await?;

    eprintln!("Prep is done!");
    eprintln!("Total target accessions: {}", total_stats.records_in_target);
    Ok(())
}

pub fn load_prelim_map(path: &str) -> Result<FxHashMap<String, u32>> {
    eprintln!("Loading prelim map from {}...", path);
    let file = File::open(path).with_context(|| format!("Failed to open {}", path))?;
    let reader = BufReader::new(file);
    let mut map = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with("accession") {
            continue;
        }
        
        let mut parts = line.split('\t');
        if let (Some(acc), Some(taxid_str)) = (parts.next(), parts.next()) {
            if let Ok(taxid) = taxid_str.trim().parse::<u32>() {
                map.insert(acc.to_string(), taxid);
            }
        }
    }
    Ok(map)
}
