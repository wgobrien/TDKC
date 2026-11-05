use tdkc::{build, build_domain, custom, inspect, prep, query, serve};
use clap::{Parser, Subcommand};
use anyhow::Context;

fn parse_positive_usize(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("'{value}' is not a valid positive integer"))?;
    if parsed == 0 {
        return Err("value must be at least 1".to_string());
    }
    Ok(parsed)
}

#[derive(Parser)]
#[command(name = "tdkc")]
#[command(about = "TDKC (Target Distilled K-mer Classifier) - Ultrafast and Memory-Efficient Metagenomic Classification for Target Pathogen Diagnostics")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Prep {
        /// You can download archaea, bacteria, viral, fungi, invertebrate, plant, protozoa, human, plasmid, univec_core, or none
        #[arg(long, conflicts_with = "custom")]
        domains: Option<String>,

        #[arg(short = 't', long)]
        targets: String,

        #[arg(short = 'd', long, default_value = "tdkc_db")]
        db: String,

        #[arg(long, default_value = "http")]
        backend: String,

        #[arg(long, default_value_t = 6)]
        concurrent_downloads: usize,

        #[arg(long, default_value_t = 2)]
        in_flight_chunks: usize,

        /// CPU workers used for accession-to-taxid lookup
        #[arg(
            short = 'j',
            long,
            default_value_t = num_cpus::get(),
            value_parser = parse_positive_usize
        )]
        threads: usize,

        /// Local FASTA to process without downloading RefSeq or UniVec sequences
        #[arg(long)]
        custom: Option<String>,

        /// Register the custom FASTA as a query-time domain
        #[arg(long, requires = "custom")]
        domain: bool,

        /// Reporting taxid when a custom domain name cannot be inferred
        #[arg(long, requires = "domain")]
        domain_taxid: Option<u32>,

        /// Skip dustmasker low-complexity masking
        #[arg(long)]
        no_mask: bool,
    },

    Build {

        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        #[arg(short = 'k', long, default_value_t = 35)]
        window_size: usize,

        #[arg(short = 'l', long, default_value_t = 31)]
        minimizer_size: usize,

        #[arg(short = 'a', long)]
        accession: bool,

    },

    Query {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = '1', long, required_unless_present = "input_dir")]
        read1: Option<String>,

        #[arg(short = '2', long, requires = "read1")]
        read2: Option<String>,

        /// Read adjacent, matching records as paired mates in each input FASTQ.
        #[arg(long, conflicts_with = "read2")]
        interleaved: bool,

        #[arg(short = 'i', long, required_unless_present = "read1")]
        input_dir: Option<String>,

        #[arg(short = 'j', long, default_value_t = 1)]
        threads: usize,

        #[arg(short = 'a', long = "accession-tracking")]
        accession: bool,

        #[arg(short = 'g', long, default_value_t = 2)]
        minimum_hit_groups: usize,

        #[arg(short = 'o', long)]
        output_prefix: String,

        #[arg(short = 'b', long)]
        background: bool,

        #[arg(long, help = "only generates .report")]
        no_output: bool,

        #[arg(short = 'z', long = "gzip-output", conflicts_with = "no_output",
              help = "write the per-read output gzip-compressed (.output.gz, multi-member gzip)")]
        gzip_output: bool,

        #[arg(short = 'f', long, requires = "no_output")]
        fast_gzip_decompression: bool,

        #[arg(short = 'p', long, default_value_t = 1, requires = "background")]
        domain_penalty: usize,

        #[arg(short = 'c', long = "counts", requires = "accession")]
        counts: bool,
    },

    BuildDomain {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        #[arg(short = 'p', long, default_value_t = 0.0001)]
        fpr: f64,

        #[arg(long)]
        bacteria: bool,

        #[arg(long)]
        archaea: bool,

        #[arg(long)]
        viral: bool,

        #[arg(long)]
        fungi: bool,
    },

    Inspect {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = 'o', long, default_value = "inspect.txt")]
        output: String,
    },

    /// Keep the database resident in memory and answer sequence batches over a socket or stdio
    Serve {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        /// Unix socket to listen on
        #[arg(long, required_unless_present = "stdio", conflicts_with = "stdio")]
        socket: Option<String>,

        /// Serve a single session over stdin/stdout instead of a socket
        #[arg(long)]
        stdio: bool,

        #[arg(short = 'a', long = "accession-tracking")]
        accession: bool,

        #[arg(short = 'b', long)]
        background: bool,
    },

}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Prep {
            domains,
            targets,
            db,
            backend,
            concurrent_downloads,
            in_flight_chunks,
            threads,
            custom,
            domain,
            domain_taxid,
            no_mask,
        } => {
            let rt = tokio::runtime::Runtime::new()?;
            rt.block_on(prep::run_prep(prep::PrepConfig {
                domains,
                targets_file: targets,
                db_dir: db,
                backend,
                concurrent_downloads,
                in_flight_chunks,
                threads,
                custom_fasta: custom,
                custom_domain: domain,
                domain_taxid,
                no_mask,
            }))?;
        }

        Commands::Build {
            db,
            threads,
            accession,
            window_size,
            minimizer_size,
        } => {
            let db_path = std::path::PathBuf::from(&db);
            let genome_dir = db_path.join("genome");

            let mut fasta_files: Vec<String> = Vec::new();
            for entry in std::fs::read_dir(&genome_dir).with_context(|| format!("Can't read {}", genome_dir.display()))? {
                let entry = entry?;
                let path = entry.path();
                if path.extension().and_then(|s| s.to_str()) != Some("fna") {
                    continue;
                }
                let fname = path.file_name().unwrap_or_default().to_string_lossy();
                if fname.starts_with('_') {
                    continue;
                }
                fasta_files.push(path.to_string_lossy().into_owned());
            }

            fasta_files.sort();

            build::run_build(build::BuildConfig {
                fasta_files,
                target_fasta_file: db_path.join("target.fasta").to_string_lossy().into_owned(),
                prelim_map_file: db_path.join("prelim_map.txt").to_string_lossy().into_owned(),
                targets_file: db_path.join("targets.txt").to_string_lossy().into_owned(),
                nodes_file: db_path.join("taxonomy/nodes.dmp").to_string_lossy().into_owned(),
                names_dmp_path: db_path.join("taxonomy/names.dmp").to_string_lossy().into_owned(),
                domain_taxids: custom::load_domains(&db_path)?.into_iter().map(|d| d.taxid).collect(),
                db_prefix: db_path.join("db").to_string_lossy().into_owned(),
                threads,
                track_accessions: accession,
                k: window_size,
                l: minimizer_size,
            })?;
        }

        Commands::Query {
            db,
            read1,
            read2,
            interleaved,
            input_dir,
            threads,
            accession,
            minimum_hit_groups,
            output_prefix,
            background,
            no_output,
            gzip_output,
            fast_gzip_decompression,
            domain_penalty,
            counts,
        } => {
            if accession && no_output && !counts {
                anyhow::bail!("With --no-output, accession tracking (-a) needs at least the count summary: add -c to write <sample>.strain.txt");
            }

            let mut samples = Vec::new();
            if let Some(dir_path) = input_dir {
                std::fs::create_dir_all(&output_prefix)?;

                samples = if interleaved {
                    query::read_finder::discover_interleaved_samples(std::path::Path::new(&dir_path))?
                } else {
                    query::read_finder::discover_samples(std::path::Path::new(&dir_path))?
                };
                let paired = samples.iter().filter(|s| s.interleaved || s.r2.is_some()).count();
                let single = samples.len() - paired;
                eprintln!("Detected {} paired-end reads and {} single-end reads", paired, single);
                
            } else if let Some(r1) = read1 {
                samples.push(query::read_finder::Sample {
                    name: "read".to_string(),
                    r1: std::path::PathBuf::from(r1),
                    r2: read2.map(std::path::PathBuf::from),
                    interleaved,
                });
            }

            if samples.is_empty() {
                anyhow::bail!("No valid FASTQ files to query, check your files!");
            }

            query::run_query(query::QueryConfig {
                db_dir: db,
                samples,
                threads,
                use_accessions: accession,
                minimum_hit_groups,
                output_prefix,
                background,
                no_output,
                gzip_output,
                gzip_multithreaded: fast_gzip_decompression,
                domain_penalty,
                counts,
            })?;
        }

        Commands::BuildDomain {
            db,
            threads,
            fpr, 
            bacteria,
            archaea,
            viral,
            fungi,
        } => {
            build_domain::run_build_domain(build_domain::BuildDomainConfig {
                db_dir: db,
                threads,
                fpr, 
                bacteria,
                archaea,
                viral,
                fungi,
            })?;
        }

        Commands::Inspect { db, output } => {
            inspect::run_inspect(inspect::InspectConfig {
                db_dir: db,
                output,
            })?;
        }

        Commands::Serve {
            db,
            threads,
            socket,
            stdio,
            accession,
            background,
        } => {
            serve::run_serve(serve::ServeConfig {
                db_dir: db,
                threads,
                use_accessions: accession,
                background,
                socket,
                stdio,
            })?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod interleaved_cli_tests {
    use super::*;

    #[test]
    fn interleaved_inputs_and_conflicts() {
        for input in ["-1", "-i"] {
            assert!(Cli::try_parse_from([
                "tdkc",
                "query",
                "-d",
                "db",
                input,
                "reads",
                "--interleaved",
                "-o",
                "out"
            ])
            .is_ok());
        }
        assert!(Cli::try_parse_from([
            "tdkc",
            "query",
            "-d",
            "db",
            "-1",
            "r1",
            "-2",
            "r2",
            "--interleaved",
            "-o",
            "out"
        ])
        .is_err());
        assert!(
            Cli::try_parse_from(["tdkc", "query", "-d", "db", "--interleaved", "-o", "out"])
                .is_err()
        );
        assert!(Cli::try_parse_from([
            "tdkc", "query", "-d", "db", "-1", "r1", "-2", "r2", "-o", "out"
        ])
        .is_ok());
    }
}
