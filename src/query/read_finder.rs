use std::collections::HashMap;
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};

pub struct Sample {
    pub name: String,
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
    pub interleaved: bool,
}

// This will automatically detect reads, and run them sequentially
pub fn discover_samples(dir: &Path) -> Result<Vec<Sample>> {
    let mut sample_map: HashMap<String, Sample> = HashMap::new();

    for entry in std::fs::read_dir(dir).context("Failed to read input dir")? {
        let entry = entry?;
        let path = entry.path();

        if !path.is_file() {
            continue;
        }

        let file_name = path.file_name().unwrap_or_default().to_string_lossy().to_string();
        
        let is_fastq = file_name.ends_with(".fastq") 
            || file_name.ends_with(".fq") 
            || file_name.ends_with(".fastq.gz") 
            || file_name.ends_with(".fq.gz");

        if !is_fastq {
            continue;
        }

        let (sample_name, is_r1, is_r2) = if let Some(idx) = file_name.find("_R1") {
            (file_name[..idx].to_string(), true, false)
        } else if let Some(idx) = file_name.find("_R2") {
            (file_name[..idx].to_string(), false, true)
        } else {
            // now just assume this is single ended
            let name = file_name.split(".f").next().unwrap_or(&file_name).to_string();
            (name, true, false)
        };

        let entry = sample_map.entry(sample_name.clone()).or_insert(Sample {
            name: sample_name,
            r1: PathBuf::new(),
            r2: None,
            interleaved: false,
        });

        if is_r1 {
            entry.r1 = path;
        } else if is_r2 {
            entry.r2 = Some(path);
        }
    }

    let mut valid_samples: Vec<Sample> = sample_map.into_values()
        .filter(|s| s.r1.exists())
        .collect();

    // sort so it processes from A01 through H12 cleanly
    valid_samples.sort_by(|a, b| a.name.cmp(&b.name));
    
    Ok(valid_samples)
}

/// Every FASTQ is an independent interleaved sample, regardless of mate-like filenames.
pub fn discover_interleaved_samples(dir: &Path) -> Result<Vec<Sample>> {
    let mut samples = Vec::new();
    let mut names = std::collections::HashSet::new();
    for entry in std::fs::read_dir(dir).context("Failed to read input dir")? {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let filename = path.file_name().unwrap_or_default().to_string_lossy();
        let Some(name) = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
            .iter()
            .find_map(|ext| filename.strip_suffix(ext))
        else {
            continue;
        };
        if !names.insert(name.to_owned()) {
            anyhow::bail!(
                "Duplicate interleaved sample name '{}' in {}",
                name,
                dir.display()
            );
        }
        samples.push(Sample {
            name: name.to_owned(),
            r1: path,
            r2: None,
            interleaved: true,
        });
    }
    samples.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interleaved_discovery_preserves_filenames_and_sorting() -> Result<()> {
        let dir = std::env::temp_dir().join(format!("tdkc-discovery-{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let result = (|| -> Result<()> {
            for name in [
                "z_R2.fq.gz",
                "z_R1.fastq",
                "a.v1.fastq.gz",
                "b.fq",
                "ignored.txt",
            ] {
                std::fs::write(dir.join(name), b"")?;
            }
            std::fs::create_dir_all(dir.join("ignored.fastq"))?;
            let samples = discover_interleaved_samples(&dir)?;
            assert_eq!(
                samples.iter().map(|s| s.name.as_str()).collect::<Vec<_>>(),
                ["a.v1", "b", "z_R1", "z_R2"]
            );
            assert!(samples.iter().all(|s| s.interleaved && s.r2.is_none()));
            let normal = discover_samples(&dir)?;
            assert!(normal.iter().all(|s| !s.interleaved));
            assert!(normal.iter().any(|s| s.name == "z" && s.r2.is_some()));
            std::fs::write(dir.join("b.fastq.gz"), b"")?;
            assert!(discover_interleaved_samples(&dir)
                .err()
                .unwrap()
                .to_string()
                .contains("Duplicate"));
            Ok(())
        })();
        std::fs::remove_dir_all(dir)?;
        result
    }
}
