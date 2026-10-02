<div align="center">

<img src="./tdkc_logo.png" width="500" alt="TDKC Logo"/>

# TDKC

### Target Distilled K-mer Classifier: Ultrafast and Memory-Efficient Metagenomic Sequence Classification for Target Pathogen Diagnostics

[![Build](https://github.com/543090lee/TDKC/actions/workflows/rust.yml/badge.svg)](https://github.com/543090lee/TDKC/actions/workflows/rust.yml)
[![Tests](https://github.com/543090lee/TDKC/actions/workflows/query-tests.yml/badge.svg)](https://github.com/543090lee/TDKC/actions/workflows/query-tests.yml)
[![Release](https://github.com/543090lee/TDKC/actions/workflows/release.yml/badge.svg)](https://github.com/543090lee/TDKC/actions/workflows/release.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/built%20with-Rust-orange?logo=rust)](https://www.rust-lang.org/)
[![Paper](https://img.shields.io/badge/paper-preprint-purple)](https://www.biorxiv.org/content/10.64898/2026.06.05.730319v1)
[![version](https://img.shields.io/github/v/release/543090lee/TDKC?style=flat-square&labelColor=21262d&color=30363d)](https://github.com/543090lee/TDKC/releases/latest)

[Paper](https://www.biorxiv.org/content/10.64898/2026.06.05.730319v1) · [Getting Started](#quick-start) · [Citation](#citation)

</div>

> [TIP]
> **Just want to run TDKC?** A prebuilt, target-distilled index is available on Zenodo. See [Prebuilt Database](#prebuilt-database).

---
## Installation

```bash
git clone https://github.com/543090lee/TDKC.git
cd TDKC
mamba env create -f environment.yml
mamba activate tdkc
cargo install --path .

```

This installs the tdkc binary to `~/.cargo/bin` (on your `PATH`)

---

## Quick Start

TDKC has three main steps: `prep` → `build` → `query`.

### 1. Prep

Download reference sequences from NCBI via HTTPS concurrently (no FTP or rsync), filter to your target taxa, and produce a local accession2taxid map. The output directory is the database, and every later command points at it.

```bash
tdkc prep \
  --domains bacteria,viral,human,archaea,plasmid \
  -t data/targets.txt \
  --db tdkc_db/
```

> **Tip:** When making your own target list, it's recommended to put genus-level target taxIDs along with species-level. Many k-mers get pushed up due to conserved regions.

**Adding a user-supplied FASTA (GenBank / WGS)**

Use `--custom` to process a local FASTA without downloading RefSeq or UniVec sequences. If the database has no `taxonomy/` directory, TDKC downloads the taxonomy and accession maps once; an existing taxonomy directory is always used offline.

```bash
tdkc prep \
  -t data/targets.txt \
  --db tdkc_db/ \
  --custom /data/my_genbank_sequences.fasta
```

Custom files are identified from the sanitized input basename plus a short SHA-256. Headers must use an NCBI accession (including its version when available) as the first whitespace-delimited token. TDKC resolves those accessions through the GenBank and WGS accession maps and stages only target-subtree records as `genome/_target.my_genbank_sequences_a1b2c3d4e5f6.fna`. Unmapped records are reported and ignored. Repeating `prep` adds more hashed FASTAs to the same database.

Pass `--domain` to retain mapped non-target records as `genome/_domain.my_genbank_sequences_a1b2c3d4e5f6.fna` and register that disjoint shard for `build-domain`. Without `--domain`, non-target records are discarded. Custom non-target records never participate in the main build's minimizer challenge, so they cannot remove target minimizers contributed by the same custom FASTA. Known names such as `human`, `bacteria`, `archaea`, `viral`, and `fungi` infer their reporting taxid; otherwise pass `--domain-taxid`.

Custom entries prepared by earlier TDKC versions keep their existing full-FASTA layout. Re-add or rebuild those entries to opt into partitioning.

```bash
tdkc prep -t data/targets.txt --db tdkc_db/ \
  --custom /data/human.fna --domain

tdkc prep -t data/targets.txt --db tdkc_db/ \
  --custom /data/wastewater.fna --domain --domain-taxid 12345
```

| Flag | Default | Description |
|------|---------|-------------|
| `--domains` | `bacteria,viral,archaea,human,plasmid` | Comma-separated RefSeq domains to download. Valid values: `bacteria`, `viral`, `archaea`, `human`, `fungi`, `invertebrate`, `plant`, `plasmid`, `protozoa`. UniVec_Core is always included automatically. |
| `-t` / `--targets` | — | Path to targets file (one NCBI taxid per line, any rank). A copy is saved into the db dir as `targets.txt`. |
| `-d` / `--db` | `tdkc_db` | Database output directory. This becomes the input for every subsequent command. |
| `--custom` | — | Process a local FASTA without downloading sequence datasets. |
| `--domain` | off | Retain mapped custom non-target records for a query-time domain Bloom filter. |
| `--domain-taxid` | inferred | Reporting taxid for an unrecognized custom domain name. |
| `--concurrent-downloads` | `6` | Number of genome files being streamed from NCBI in parallel. I don't recommend going over 6, NCBI server might complain... |
| `-j` / `--threads` | all cores | CPU workers used to scan accession-to-taxid tables. |
| `--in-flight-chunks` | `2` | Number of dust-masking jobs to run in parallel. Increase this separately from `-j` when memory allows. |
| `--no-mask` | - | Skips dust-masking |



### 2. Build

Distill the target k-mer index. 
```bash
tdkc build --db tdkc_db/ -j 32
```

Add `-a` to enable per-minimizer accession tracking (TDKC-A mode).

You can find an example `targets.txt` in the `data/` directory of this repo. The target taxa list covers the most common human respiratory and enteric viruses — feel free to use it as-is or as a starting point.

| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--db` | — | Database directory (output of `prep`) |
| `-j` | all cores | Threads |
| `-k` | `35` | Window size k |
| `-l` | `31` | Minimizer length l |
| `-a` | off | Enable accession tracking |

---

### 3. Query

#### Default Classification

Query a single sample:

```bash
tdkc query \
  --db tdkc_db/ \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -j 32 \
  -o results/sample
```

Or query an entire directory of FASTQ files at once:

```bash
tdkc query \
  --db tdkc_db/ \
  -i /data/fastq_dir/ \
  -j 32 \
  -o results/
```

TDKC auto-detects reads in the directory and processes them sequentially. Output files are written into the directory specified by `-o`.

Outputs: `results/sample.output` (per-read) and `results/sample.report`.

#### Subtype Classification

If you would like to also classify which specific subtype/strain the sample has,  make sure you built the database with accession tracking (-a).

Query with `-a`, `-c` flags:

```bash
tdkc query \
  --db tdkc_db/ \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -j 32 \
  -o results/sample \
  -a -c
```

Using these flags provides detailed strain-level information:
* `-a`: Includes the per k-mer set of accession hits in .output file
* `-c`: Generate `results/sample.strain.txt` summary file that includes the per-accession total and unique hit counts for the entire sample


| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--db` | — | Database directory |
| `-1` | — | R1 FASTQ (required unless `-i` is used) |
| `-2` | — | R2 FASTQ (optional, for paired-end) |
| `-i` | — | Input directory of FASTQ files (alternative to `-1`/`-2`) |
| `-g` | `2` | Min distinct minimizer hit groups |
| `-a` | off | Output per k-mer accession hits |
| `-b` | off | Enable domain-level detection (requires built bloom filters) |
| `-p` | `1` | Penalty applied to domain hit group prevent FP read classifications |
| `-f` | off | Enable multithreaded gzip decompression. Recommended for short-read (50-100bp) |
| `-c` | off | Output per-accession total and unique hit counts to .strain.txt |
| `--no-output` | off | only generates .report |



---

### (Optional) Build Domain Bloom Filters

Domain bloom filters are driven by `<db>/genome/domains.tsv`. `prep` registers every downloaded domain and the disjoint non-target shard of every custom FASTA passed with `--domain`. With no selection flags, `build-domain` builds every registered entry; built-in flags such as `--bacteria` filter registry entries by label, including hashed custom files with that label.

```bash
# Build all available domains
tdkc build-domain --db tdkc_db/ -j 32

# Or pick specific domains
tdkc build-domain \
  --db tdkc_db/ \
  --bacteria \
  --viral \
  --archaea \
  -j 32
```

Activate at query time with `-b` to classify against broad domains too.

| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--db` | — | Database directory |
| `-p` | `0.0001` | Bloom filter false positive rate (e.g. `0.001` for 0.1%) |
| `-j` | all cores | Threads |
| `--bacteria` / `--archaea` / `--viral` / `--fungi` | off | Select domains you want to include only. If none are passed, all available are built. |

> **Note:** Lower FPR values produce more accurate domain classification but require more memory. The default 0.01% FPR is recommended for most use cases — raising it significantly (e.g. above 0.1%) may introduce spurious hits...

---

### Inspect TDKC Database

Inspect pre-built TDKC database to analyze k-mer compositions of targets

```bash
tdkc inspect --db tdkc_db/ 
```

This will output inspect.txt with number of minimizers each taxID has and its ratio.

---

## Resident Database Server (`serve`)

`tdkc query` reloads the database on every invocation. `tdkc serve` loads it **once**, keeps it
resident in memory, and answers many sequence batches with low latency — useful for interactive
or high-throughput querying where the load cost would otherwise dominate.

```bash
# Unix socket (accepts many concurrent clients)
tdkc serve --db tdkc_db/ -j 32 --socket /tmp/tdkc.sock

# Or a single session over stdin/stdout
tdkc serve --db tdkc_db/ -j 32 --stdio
```

It accepts the same `-a` (accession tracking) and `-b` (background domain filters) flags as
`query`. On startup it prints the database load time and then a `READY <target>` line on stderr.

Each connection is handled on its own thread and shares one `Arc<Engine>`, so the database is
never copied. Batches are dispatched across a single shared rayon pool, meaning both one large
batch and many concurrent small batches keep all `-j` threads busy.

### Wire format

All integers are little-endian `u32`. A connection may send many requests back to back.

```text
Request:   n_seqs, then n_seqs * ( len, len bytes of sequence )
Response:  n_seqs, then n_seqs * ( n_kmers, n_kmers * taxid )
Error:     0xFFFFFFFF, then msg_len, msg_len bytes of UTF-8
```

Responses are index-aligned with the request. Each sequence yields one taxid **per k-mer, in
sequence order**: `0` for a miss, the resolved taxid for a database hit, and the domain taxid for
a background bloom-filter hit. Sequences shorter than `k` yield zero k-mers. Callers get this raw
per-k-mer signal and can apply their own LCA or scoring logic.

The format is binary rather than JSON so a client can decode a response body with a single
`numpy.frombuffer`:

```python
import socket, struct, numpy as np

def query(sock_path, seqs):
    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    s.connect(sock_path)
    f = s.makefile("rwb")
    f.write(struct.pack("<I", len(seqs)))
    for q in seqs:
        f.write(struct.pack("<I", len(q)) + q.encode())
    f.flush()

    (n,) = struct.unpack("<I", f.read(4))
    out = []
    for _ in range(n):
        (k,) = struct.unpack("<I", f.read(4))
        out.append(np.frombuffer(f.read(4 * k), dtype="<u4"))
    return out

print(query("/tmp/tdkc.sock", ["TAGTGTGTGCCACTTGTGAACAGATTGCTGATTCACAGCATCGGTCTCAC"]))
```

### Embedding in Rust

TDKC also builds as a library. `tdkc::engine::Engine` is the same resident core `serve` uses:

```rust
use tdkc::engine::{Engine, EngineConfig};

let engine = Engine::load(&EngineConfig::new("tdkc_db"))?;   // pay the load once
let hits = engine.query_batch(&[b"ACGT...".to_vec()]);       // Vec<Vec<u32>>, one taxid per k-mer
```

`Engine` is `Sync` and holds no per-request state, so an `Arc<Engine>` can be shared freely.

---


## Prebuilt Database
  
**https://doi.org/10.5281/zenodo.21343246** &nbsp;·&nbsp; CC-BY-4.0
 
```bash
wget "https://zenodo.org/records/21343246/files/tdkc_db_v1.0.0.tar.gz"
tar xzf tdkc_db_v1.0.0.tar.gz
 
# Verify integrity
cd tdkc_db && sha256sum -c checksums.sha256 && cd ..
 
# Classify
tdkc query \
  --db tdkc_db/ \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -j 32 \
  -o results
```

> **Coming soon:** We are planning to build a dedicated platform serving periodically updated TDKC databases

---

## Citation

```bibtex
@article{lee2026tdkc,
  title   = {TDKC: Ultrafast and Memory-Efficient Sequence Classification
             for Target Pathogen Diagnostics},
  author  = {Lee, Seungmo and Eskin, Eleazar},
  year    = {2026},
  url     = {https://github.com/543090lee/TDKC}
}
```

---

## License

MIT © Seungmo Lee
