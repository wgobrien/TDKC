use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CustomDomain {
    pub identity: String,
    pub label: String,
    pub taxid: u32,
    pub fasta: String,
    pub sha256: String,
}

pub struct CustomIdentity {
    pub label: String,
    pub identity: String,
    pub sha256: String,
    pub accessions: HashSet<String>,
    pub unparseable_headers: u64,
}

struct HashingReader<R> {
    inner: R,
    hasher: Sha256,
}

impl<R> HashingReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    fn finalize(self) -> String {
        hex::encode(self.hasher.finalize())
    }
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.hasher.update(&buf[..n]);
        Ok(n)
    }
}

fn label_from_path(path: &Path) -> Result<String> {
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow::anyhow!("Custom FASTA has no valid UTF-8 filename"))?;
    let mut stem = name.to_string();
    if stem.to_ascii_lowercase().ends_with(".gz") {
        stem.truncate(stem.len() - 3);
    }
    for suffix in [".fasta", ".fna", ".fa"] {
        if stem.to_ascii_lowercase().ends_with(suffix) {
            stem.truncate(stem.len() - suffix.len());
            break;
        }
    }
    let mut label = String::new();
    let mut separator = false;
    for ch in stem.chars() {
        if ch.is_ascii_alphanumeric() {
            label.push(ch.to_ascii_lowercase());
            separator = false;
        } else if !separator && !label.is_empty() {
            label.push('_');
            separator = true;
        }
    }
    while label.ends_with('_') {
        label.pop();
    }
    if label.is_empty() {
        anyhow::bail!("Custom FASTA filename has no usable name");
    }
    Ok(label)
}

fn collect_accessions<R: BufRead>(reader: &mut R) -> Result<(HashSet<String>, u64)> {
    let mut wanted = HashSet::new();
    let mut unparseable_headers = 0u64;
    let mut line = Vec::with_capacity(256);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.first() != Some(&b'>') {
            continue;
        }
        let accession = line[1..]
            .split(|b| b.is_ascii_whitespace())
            .next()
            .unwrap_or_default();
        if accession.is_empty() {
            unparseable_headers += 1;
            continue;
        }
        match String::from_utf8(accession.to_vec()) {
            Ok(accession) => {
                wanted.insert(accession);
            }
            Err(_) => unparseable_headers += 1,
        }
    }
    Ok((wanted, unparseable_headers))
}

pub fn inspect(path: &Path) -> Result<CustomIdentity> {
    let label = label_from_path(path)?;
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let gzipped = path
        .to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".gz");

    let (accessions, unparseable_headers, sha256) = if gzipped {
        let hashing = HashingReader::new(file);
        let decoder = GzDecoder::new(hashing);
        let mut reader = BufReader::with_capacity(1024 * 1024, decoder);
        let (accessions, unparseable_headers) = collect_accessions(&mut reader)?;
        let decoder = reader.into_inner();
        let mut hashing = decoder.into_inner();
        std::io::copy(&mut hashing, &mut std::io::sink())?;
        (accessions, unparseable_headers, hashing.finalize())
    } else {
        let hashing = HashingReader::new(file);
        let mut reader = BufReader::with_capacity(1024 * 1024, hashing);
        let (accessions, unparseable_headers) = collect_accessions(&mut reader)?;
        let mut hashing = reader.into_inner();
        std::io::copy(&mut hashing, &mut std::io::sink())?;
        (accessions, unparseable_headers, hashing.finalize())
    };

    let identity = format!("{}_{}", label, &sha256[..12]);
    Ok(CustomIdentity {
        label,
        identity,
        sha256,
        accessions,
        unparseable_headers,
    })
}

pub fn validate_existing_taxonomy(tax_dir: &Path) -> Result<()> {
    let required = [
        "nodes.dmp",
        "names.dmp",
        "nucl_gb.accession2taxid",
        "nucl_wgs.accession2taxid",
    ];
    let missing: Vec<&str> = required
        .into_iter()
        .filter(|n| !tax_dir.join(n).is_file())
        .collect();
    if !missing.is_empty() {
        anyhow::bail!("Existing taxonomy directory {} is incomplete; missing: {}. Custom mode will not download into an existing taxonomy directory.", tax_dir.display(), missing.join(", "));
    }
    Ok(())
}

pub fn default_domain_taxid(label: &str) -> Option<u32> {
    match label {
        "bacteria" => Some(crate::taxonomy::TAXID_BACTERIA),
        "archaea" => Some(crate::taxonomy::TAXID_ARCHAEA),
        "viral" | "viruses" => Some(crate::taxonomy::TAXID_VIRAL),
        "fungi" => Some(crate::taxonomy::TAXID_FUNGI),
        "human" => Some(9606),
        "invertebrate" => Some(33208),
        "plant" => Some(33090),
        "protozoa" => Some(2759),
        "plasmid" => Some(36549),
        "univec" | "univec_core" => Some(28384),
        _ => None,
    }
}

pub fn resolve_taxid(label: &str, supplied: Option<u32>, names_path: &Path) -> Result<u32> {
    if let Some(taxid) = supplied {
        return Ok(taxid);
    }
    if let Some(taxid) = default_domain_taxid(label) {
        return Ok(taxid);
    }
    let wanted = label.replace('_', " ");
    let mut matches = HashSet::new();
    for line in BufReader::new(File::open(names_path)?).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split("\t|\t").collect();
        if fields.len() >= 2 && fields[1].trim().eq_ignore_ascii_case(&wanted) {
            if let Ok(taxid) = fields[0].trim().parse::<u32>() {
                matches.insert(taxid);
            }
        }
    }
    match matches.len() {
        1 => Ok(*matches.iter().next().unwrap()),
        0 => anyhow::bail!(
            "Cannot infer a taxid for custom domain '{}'; pass --domain-taxid",
            label
        ),
        _ => anyhow::bail!(
            "Custom domain '{}' is ambiguous; pass --domain-taxid",
            label
        ),
    }
}

fn registry_path(db_dir: &Path) -> PathBuf {
    db_dir.join("genome/domains.tsv")
}

pub fn load_domains(db_dir: &Path) -> Result<Vec<CustomDomain>> {
    let path = registry_path(db_dir);
    let file = match File::open(&path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
        Err(e) => return Err(e.into()),
    };
    let mut domains = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with("identity\t") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 5 {
            anyhow::bail!("Malformed domain registry {}", path.display());
        }
        domains.push(CustomDomain {
            identity: f[0].into(),
            label: f[1].into(),
            taxid: f[2].parse()?,
            fasta: f[3].into(),
            sha256: f[4].into(),
        });
    }
    Ok(domains)
}

pub fn save_domains(db_dir: &Path, domains: &[CustomDomain]) -> Result<()> {
    let path = registry_path(db_dir);
    let tmp = path.with_extension("tsv.tmp");
    let mut out = File::create(&tmp)?;
    writeln!(out, "identity\tlabel\ttaxid\tfasta\tsha256")?;
    for d in domains {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            d.identity, d.label, d.taxid, d.fasta, d.sha256
        )?;
    }
    out.flush()?;
    std::fs::rename(tmp, path)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn known_domain_taxids() {
        let missing = Path::new("not-needed");
        assert_eq!(resolve_taxid("human", None, missing).unwrap(), 9606);
        assert_eq!(resolve_taxid("bacteria", None, missing).unwrap(), 2);
        assert_eq!(default_domain_taxid("plasmid"), Some(36549));
        assert_eq!(default_domain_taxid("univec_core"), Some(28384));
    }

    #[test]
    fn sanitizes_fasta_names() {
        assert_eq!(
            label_from_path(Path::new("My Samples.fasta.gz")).unwrap(),
            "my_samples"
        );
        assert_eq!(label_from_path(Path::new("viral.fna")).unwrap(), "viral");
    }

    #[test]
    fn content_hash_changes_identity() {
        let base = std::env::temp_dir().join(format!("tdkc-custom-test-{}", std::process::id()));
        std::fs::create_dir_all(&base).unwrap();
        let path = base.join("Sample.fa");
        std::fs::write(&path, b">A\nAAAA\n").unwrap();
        let first = inspect(&path).unwrap();
        let again = inspect(&path).unwrap();
        assert_eq!(first.accessions, HashSet::from(["A".to_string()]));
        assert_eq!(first.identity, again.identity);
        std::fs::write(&path, b">A\nCCCC\n").unwrap();
        assert_ne!(first.identity, inspect(&path).unwrap().identity);
        std::fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn inspects_gzip_in_one_pass() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let base = std::env::temp_dir().join(format!(
            "tdkc-custom-gzip-test-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&base).unwrap();
        let path = base.join("Sample.fa.gz");
        let file = File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder
            .write_all(b">NC_000001.1 description\nACGT\n>ABCD01000001.2\nNN\n")
            .unwrap();
        encoder.finish().unwrap();

        let raw = std::fs::read(&path).unwrap();
        let expected_hash = hex::encode(Sha256::digest(&raw));
        let info = inspect(&path).unwrap();
        assert_eq!(info.sha256, expected_hash);
        assert_eq!(
            info.accessions,
            HashSet::from([
                "NC_000001.1".to_string(),
                "ABCD01000001.2".to_string()
            ])
        );
        std::fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn ignores_headers_without_valid_accession_tokens() {
        let input = b">\nAAAA\n>GOOD.1 description\nCCCC\n>\xff\nGGGG\n";
        let (accessions, ignored) = collect_accessions(&mut &input[..]).unwrap();
        assert_eq!(accessions, HashSet::from(["GOOD.1".to_string()]));
        assert_eq!(ignored, 2);
    }
}
