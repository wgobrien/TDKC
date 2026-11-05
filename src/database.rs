use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use anyhow::{Context, Result};
use boomphf::Mphf;
use rustc_hash::FxHashMap;
use crate::compression::{delta_vbyte_encode, delta_vbyte_skip, vbyte_decode};

pub struct AccessionRegistry {
    name_to_id: HashMap<String, u32>,
    id_to_name: Vec<String>,
}

impl AccessionRegistry {
    pub fn new() -> Self {
        Self {
            name_to_id: HashMap::new(),
            id_to_name: Vec::new(),
        }
    }

    pub fn get_or_create(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.name_to_id.get(name) {
            return id;
        }
        let id = self.id_to_name.len() as u32;
        self.id_to_name.push(name.to_string());
        self.name_to_id.insert(name.to_string(), id);
        id
    }

    pub fn get_name(&self, id: u32) -> &str {
        self.id_to_name.get(id as usize).map(|s| s.as_str()).unwrap_or("Unknown")
    }

    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }

    pub fn save(&self, path: &str) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "ID\tAccession")?;
        for (i, name) in self.id_to_name.iter().enumerate() {
            writeln!(w, "{}\t{}", i, name)?;
        }
        w.flush().map_err(Into::into)
    }

    pub fn load_for_query(path: &str) -> Result<Self> {
        let mut id_to_name = Vec::new();

        let mut lines = BufReader::new(
            File::open(path).context("Cannot open accession registry")?
        ).lines();
        lines.next();

        for line in lines {
            let line = line?;
            if let Some(tab_pos) = line.find('\t') {
                id_to_name.push(line[tab_pos + 1..].to_string());
            }
        }
        Ok(Self { name_to_id: HashMap::new(), id_to_name })
    }
}

pub struct EqClassAccessions {
    bitset: Vec<u64>,
    dense_ids: Vec<u32>,
    num_slots: usize,
    offsets: Vec<u32>,
    blob: Vec<u8>,
    pub num_classes: usize,
}

impl EqClassAccessions {
    pub fn new_empty(num_slots: usize) -> Self {
        let num_words = (num_slots + 63) / 64;
        let encoded_empty = delta_vbyte_encode(&[]);
        Self {
            bitset: vec![0u64; num_words],
            dense_ids: Vec::new(),
            num_slots,
            offsets: vec![0u32],
            blob: encoded_empty,
            num_classes: 1,
        }
    }

    #[inline]
    fn set_bit(&mut self, idx: usize) {
        let word = idx / 64;
        let bit = idx % 64;
        self.bitset[word] |= 1u64 << bit;
    }

    fn matches(&self, class_id: u32, acc_list: &[u32]) -> bool {
        let mut pos = self.offsets[class_id as usize] as usize;
        let count = vbyte_decode(&self.blob, &mut pos) as usize;
        if count != acc_list.len() {
            return false;
        }
        let mut prev = 0u32;
        for &expected in acc_list {
            let delta = vbyte_decode(&self.blob, &mut pos);
            prev += delta;
            if prev != expected {
                return false;
            }
        }
        true
    }

    pub fn add_accessions(
        &mut self,
        idx: usize,
        acc_list: &[u32],
        hash_to_class: &mut FxHashMap<u64, u32>,
        // This build map is temp and will be dropped before save. This is just for fast look up
        build_map: &mut FxHashMap<u32, u32>,
    ) {
        if acc_list.is_empty() {
            return;
        }

        let mut h = crate::hash::hash_slice(acc_list);
        let class_id = loop {
            match hash_to_class.get(&h) {
                Some(&id) => {
                    if self.matches(id, acc_list) {
                        break id;
                    }
                    // Hash collision... linear probe
                    h = h.wrapping_add(0x9E3779B97F4A7C15);
                }
                None => {
                    let id = self.num_classes as u32;
                    self.offsets.push(self.blob.len() as u32);
                    self.blob.extend_from_slice(&delta_vbyte_encode(acc_list));
                    hash_to_class.insert(h, id);
                    self.num_classes += 1;
                    break id;
                }
            }
        };

        self.set_bit(idx);
        build_map.insert(idx as u32, class_id);
    }

    pub fn finalize_for_save(&mut self, build_map: FxHashMap<u32, u32>) {
        let total_set: usize = self.bitset.iter().map(|w| w.count_ones() as usize).sum();
        self.dense_ids = vec![0u32; total_set];

        let mut dense_pos = 0usize;
        for word_idx in 0..self.bitset.len() {
            let mut word = self.bitset[word_idx];
            while word != 0 {
                let bit = word.trailing_zeros() as usize;
                let idx = word_idx * 64 + bit;
                if let Some(&class_id) = build_map.get(&(idx as u32)) {
                    self.dense_ids[dense_pos] = class_id;
                }
                dense_pos += 1;
                word &= word - 1;
            }
        }
    }

    pub fn into_parts(self) -> (Vec<u64>, Vec<u32>, Vec<u32>, Vec<u8>, usize) {
        (self.bitset, self.dense_ids, self.offsets, self.blob, self.num_classes)
    }

    pub fn save(&self, path: &str) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        let num_entries = self.dense_ids.len() as u64;
        let num_classes = self.num_classes as u64;
        let num_slots = self.num_slots as u64;
        w.write_all(&num_entries.to_le_bytes())?;
        w.write_all(&num_classes.to_le_bytes())?;
        w.write_all(&num_slots.to_le_bytes())?;

        let width: u8 = if num_classes <= 256 { 1 } else if num_classes <= 65536 { 2 } else { 4 };
        w.write_all(&[width])?;

        let num_words = self.bitset.len() as u64;
        w.write_all(&num_words.to_le_bytes())?;
        let bitset_bytes = unsafe {
            std::slice::from_raw_parts(
                self.bitset.as_ptr() as *const u8,
                self.bitset.len() * 8,
            )
        };
        w.write_all(bitset_bytes)?;

        match width {
            1 => w.write_all(&self.dense_ids.iter().map(|&c| c as u8).collect::<Vec<_>>())?,
            2 => {
                for &c in &self.dense_ids {
                    w.write_all(&(c as u16).to_le_bytes())?;
                }
            }
            _ => {
                let id_bytes = unsafe {
                    std::slice::from_raw_parts(
                        self.dense_ids.as_ptr() as *const u8,
                        self.dense_ids.len() * 4,
                    )
                };
                w.write_all(id_bytes)?;
            }
        }

        w.write_all(&(self.blob.len() as u64).to_le_bytes())?;
        w.write_all(&self.blob)?;
        w.flush().map_err(Into::into)
    }

    pub fn load(path: &str) -> Result<Self> {
        let mut f = BufReader::new(
            File::open(path).context("Cant open eq class accession file")?
        );
        let mut buf8 = [0u8; 8];
        let mut buf1 = [0u8; 1];

        f.read_exact(&mut buf8)?;
        let num_entries = u64::from_le_bytes(buf8) as usize;
        f.read_exact(&mut buf8)?;
        let num_classes = u64::from_le_bytes(buf8) as usize;
        f.read_exact(&mut buf8)?;
        let num_slots = u64::from_le_bytes(buf8) as usize;
        f.read_exact(&mut buf1)?;
        let width = buf1[0];

        f.read_exact(&mut buf8)?;
        let num_words = u64::from_le_bytes(buf8) as usize;
        let mut bitset = vec![0u64; num_words];
        let bitset_bytes = unsafe {
            std::slice::from_raw_parts_mut(
                bitset.as_mut_ptr() as *mut u8,
                num_words * 8,
            )
        };
        f.read_exact(bitset_bytes)?;

        let dense_ids: Vec<u32> = match width {
            1 => {
                let mut bytes = vec![0u8; num_entries];
                f.read_exact(&mut bytes)?;
                bytes.into_iter().map(|b| b as u32).collect()
            }
            2 => {
                let mut bytes = vec![0u8; num_entries * 2];
                f.read_exact(&mut bytes)?;
                bytes.chunks_exact(2)
                    .map(|c| u16::from_le_bytes([c[0], c[1]]) as u32)
                    .collect()
            }
            4 => {
                let mut raw = vec![0u32; num_entries];
                let byte_slice = unsafe {
                    std::slice::from_raw_parts_mut(
                        raw.as_mut_ptr() as *mut u8,
                        num_entries * 4,
                    )
                };
                f.read_exact(byte_slice)?;
                raw
            }
            _ => anyhow::bail!("Invalid class_id width: {}", width),
        };

        f.read_exact(&mut buf8)?;
        let blob_len = u64::from_le_bytes(buf8) as usize;
        let mut blob = vec![0u8; blob_len];
        f.read_exact(&mut blob)?;

        let mut offsets = Vec::with_capacity(num_classes);
        let mut bpos: usize = 0;
        for _ in 0..num_classes {
            offsets.push(bpos as u32);
            delta_vbyte_skip(&blob, &mut bpos);
        }

        Ok(Self {
            bitset,
            dense_ids,
            num_slots,
            offsets,
            blob,
            num_classes,
        })
    }
}

pub struct KmerDatabase {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    num_minimizers: usize,
    mphf: Mphf<u64>,
    fingerprints: Vec<u16>,
    taxid_indices: Vec<TaxIdx>,
    pub index_to_taxid: Vec<u32>,
    pub ancestor_matrix: Vec<u8>,
    accessions: Option<EqClassAccessions>,
    pub names_map: Option<rustc_hash::FxHashMap<u32, String>>,
}

impl KmerDatabase {

    // This is called by build.rs 
    pub fn new(
        k: usize,
        l: usize,
        spaced_seed_mask: u64,
        toggle_mask: u64,
        num_minimizers: usize,
        mphf: Mphf<u64>,
        fingerprints: Vec<u16>,
        taxid_indices: Vec<TaxIdx>,
        index_to_taxid: Vec<u32>,
        ancestor_matrix: Vec<u8>,
        accessions: Option<EqClassAccessions>,
        names_map: Option<rustc_hash::FxHashMap<u32, String>>,
    ) -> Self {
        Self {
            k, l, spaced_seed_mask, toggle_mask,
            num_minimizers, mphf, fingerprints, taxid_indices,
            index_to_taxid, ancestor_matrix, accessions, names_map
        }
    }

    pub fn true_taxid(&self, idx: TaxIdx) -> u32 {
        self.index_to_taxid.get(idx as usize).copied().unwrap_or(0)
    }

    pub fn get_taxid_name(&self, taxid: u32) -> String {
        if taxid == u32::MAX {
            return "\t".to_string();
        }

        if let Some(names_map) = &self.names_map {
            if let Some(name) = names_map.get(&taxid) {
                return name.clone();
            }
        }
        "Unknown".to_string()
    }
    
    pub fn taxid_indices(&self) -> &[TaxIdx] {
        &self.taxid_indices
    }

    pub fn save(&self, prefix: &str) -> Result<()> {
        // .meta
        {
            let mut w = BufWriter::new(File::create(format!("{}.meta", prefix))?);
            w.write_all(&(self.k as u32).to_le_bytes())?;
            w.write_all(&(self.l as u32).to_le_bytes())?;
            w.write_all(&self.spaced_seed_mask.to_le_bytes())?;
            w.write_all(&self.toggle_mask.to_le_bytes())?;
            w.write_all(&(self.num_minimizers as u64).to_le_bytes())?;
            w.write_all(&[self.accessions.is_some() as u8])?;
            w.write_all(&[taxid_width(self.index_to_taxid.len()) as u8])?;
            w.flush()?;
        }
        // .mphf
        {
            let encoded = bincode::serialize(&self.mphf)?;
            std::fs::write(format!("{}.mphf", prefix), &encoded)?;
        }
        // .fp
        {
            let mut w = BufWriter::new(File::create(format!("{}.fp", prefix))?);
            let fp_bytes = unsafe {
                std::slice::from_raw_parts(
                    self.fingerprints.as_ptr() as *const u8,
                    self.fingerprints.len() * 2,
                )
            };
            w.write_all(fp_bytes)?;
            w.flush()?;
        }
        {
            let mut w = BufWriter::new(File::create(format!("{}.taxid", prefix))?);
            if taxid_width(self.index_to_taxid.len()) == 1 {
                let bytes: Vec<u8> = self.taxid_indices.iter().map(|&t| t as u8).collect();

                w.write_all(&bytes)?;
            } else {
                let bytes = unsafe {
                    std::slice::from_raw_parts(
                        self.taxid_indices.as_ptr() as *const u8,
                        self.taxid_indices.len() * 2,
                    )
                };
                w.write_all(bytes)?;
            }
            w.flush()?;
        }
        // .taxmap - binary
        {
            let mut w = BufWriter::new(File::create(format!("{}.taxmap", prefix))?);
            w.write_all(&(self.index_to_taxid.len() as u64).to_le_bytes())?;
            let taxmap_bytes = unsafe {
                std::slice::from_raw_parts(
                    self.index_to_taxid.as_ptr() as *const u8,
                    self.index_to_taxid.len() * 4,
                )
            };
            w.write_all(taxmap_bytes)?;
            w.write_all(&self.ancestor_matrix)?;
            w.flush()?;
        }
        
        //names.dmp is saved separately and directly in build.rs

        if let Some(ref acc) = self.accessions {
            acc.save(&format!("{}.accession", prefix))?;
        }
        Ok(())
    }

    pub fn load(prefix: &str, load_accessions: bool) -> Result<Self> {
        eprintln!("Loading database from {}", prefix);

        let MetaHeader { k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc, taxid_width: meta_width } = read_meta(prefix)?;

        // .mphf
       let mphf: Mphf<u64> = {
            let f = File::open(format!("{}.mphf", prefix)).context("Cannot open .mphf file")?;
            let reader = BufReader::new(f);
            bincode::deserialize_from(reader)?
        };

        // .fp
        let fingerprints = {
            let mut f = File::open(format!("{}.fp", prefix)).context("Cannot open .fp file")?;
            let mut fp = vec![0u16; num_minimizers];
            let byte_slice = unsafe {
                std::slice::from_raw_parts_mut(fp.as_mut_ptr() as *mut u8, num_minimizers * 2)
            };
            f.read_exact(byte_slice)?;
            fp
        };

        let taxid_indices = {
            let mut f = File::open(format!("{}.taxid", prefix)).context("Cannot open .taxid file")?;
            let width = infer_taxid_width(f.metadata()?.len(), num_minimizers)?;
            check_taxid_width(meta_width, width)?;
            let mut raw = vec![0u8; num_minimizers * width];
            f.read_exact(&mut raw)?;
            if width == 1 {
                raw.into_iter().map(|b| b as TaxIdx).collect::<Vec<_>>()
            } else {
                raw.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()
            }
        };

        let (index_to_taxid, ancestor_matrix) = read_taxmap(prefix)?;
        let names_map = read_names(prefix);

        let acc_path = format!("{}.accession", prefix);
        let accessions = if load_accessions && has_acc && Path::new(&acc_path).exists() {
            Some(EqClassAccessions::load(&acc_path)?)
        } else {
            None
        };
        Ok(Self {
            k, l, spaced_seed_mask, toggle_mask, num_minimizers,
            mphf, fingerprints, taxid_indices, index_to_taxid, ancestor_matrix, accessions, names_map, 
        })
    }
}

pub struct MetaHeader {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub num_minimizers: usize,
    pub has_acc: bool,
    pub taxid_width: usize,
}

pub fn read_meta(prefix: &str) -> Result<MetaHeader> {
    let mut f = File::open(format!("{}.meta", prefix)).context("Cannot open .meta file")?;
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];
    let mut buf1 = [0u8; 1];
    f.read_exact(&mut buf4)?; let k = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf4)?; let l = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf8)?; let spaced_seed_mask = u64::from_le_bytes(buf8);
    f.read_exact(&mut buf8)?; let toggle_mask = u64::from_le_bytes(buf8);
    f.read_exact(&mut buf8)?; let num_minimizers = u64::from_le_bytes(buf8) as usize;
    f.read_exact(&mut buf1)?; let has_acc = buf1[0] == 1;
    let taxid_width = read_meta_taxid_width(&mut f)?;
    Ok(MetaHeader { k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc, taxid_width })
}

pub fn read_taxmap(prefix: &str) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut f = File::open(format!("{}.taxmap", prefix)).context("Cannot open .taxmap file")?;
    let mut all = Vec::new();
    f.read_to_end(&mut all)?;
    if all.len() < 8 {
        anyhow::bail!(".taxmap is truncated");
    }
    let sz = u64::from_le_bytes(all[0..8].try_into().unwrap()) as usize;
    if all.len() < 8 + sz * 4 {
        anyhow::bail!(".taxmap is truncated");
    }
    let v: Vec<u32> = all[8..8 + sz * 4].chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
    let mut anc = vec![0u8; sz * sz];
    let rest = &all[8 + sz * 4..];
    if rest.len() >= sz * sz {
        anc.copy_from_slice(&rest[..sz * sz]);
    }
    Ok((v, anc))
}

pub fn read_names(prefix: &str) -> Option<FxHashMap<u32, String>> {
    let file = File::open(format!("{}.names", prefix)).ok()?;
    let mut map = FxHashMap::default();
    for line in BufReader::new(file).lines().flatten() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 2 {
            if let Ok(taxid) = parts[0].parse::<u32>() {
                map.insert(taxid, parts[1].to_string());
            }
        }
    }
    Some(map)
}

pub type TaxIdx = u16;
pub const MAX_TAXA: usize = 65_536;

#[inline]
pub fn taxid_width(num_taxa: usize) -> usize {
    if num_taxa <= 256 { 1 } else { 2 }
}

pub fn read_meta_taxid_width<R: Read>(f: &mut R) -> Result<usize> {
    let mut b = [0u8; 1];
    match f.read(&mut b)? {
        0 => Ok(1),
        _ => match b[0] {
            1 => Ok(1),
            2 => Ok(2),
            x => anyhow::bail!(".meta taxid width byte is {}, expected 1 or 2", x),
        },
    }
}

pub fn check_taxid_width(meta_width: usize, file_width: usize) -> Result<()> {
    if meta_width != file_width {
        anyhow::bail!(".meta says {} byte(s) per taxid but .taxid size implies {}", meta_width, file_width);
    }
    Ok(())
}

pub fn infer_taxid_width(file_len: u64, num_minimizers: usize) -> Result<usize> {
    let n = num_minimizers as u64;
    if file_len == n { Ok(1) }
    else if file_len == n * 2 { Ok(2) }
    else { anyhow::bail!(".taxid length {} is neither {} nor {} bytes (num_minimizers={})", file_len, n, n * 2, n) }
}
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct DomainBloomFilter {
    pub data: Vec<u8>,
    pub num_hashes: u32,
    pub num_bits: u64,
}

impl DomainBloomFilter {
    pub fn new(data: Vec<u8>, num_hashes: u32, num_bits: u64) -> Self {
        Self { data, num_hashes, num_bits }
    }

    #[inline]
    pub fn contains(&self, item: u64) -> bool {
        let (h1, h2) = crate::hash::bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            if self.data[byte_idx] & bit_mask == 0 {
                return false;
            }
        }
        true
    }

    pub fn load(path: &str) -> Result<Self> {
        // `data` is a plain Vec<u8>, and serde has no specialization for byte vectors, so
        // bincode::deserialize_from walks it as a *sequence* -- one deserialize_u8 call per byte.
        // On a 42 GB filter that is ~42 billion calls through a BufReader and takes minutes.
        //
        // bincode encodes a seq as a fixed-width u64 length followed by the elements, and for u8
        // elements those are just raw bytes, so the on-disk layout is:
        //     u64 data_len | data_len raw bytes | u32 num_hashes | u64 num_bits
        // which we can pull in with bulk reads. Byte-identical to what bincode writes, so existing
        // .bloom files load unchanged and files written by older builds stay compatible.
        let f = File::open(path)?;
        let file_len = f.metadata()?.len();

        let mut header = [0u8; 8];
        (&f).read_exact(&mut header)?;
        let data_len = u64::from_le_bytes(header);

        // Validate the layout before trusting it; fall back to bincode if anything is off (e.g. a
        // file written with non-default bincode options).
        let trailer_len = 4 + 8; // num_hashes + num_bits
        if data_len.checked_add(8 + trailer_len) != Some(file_len) {
            let mut reader = BufReader::new(File::open(path)?);
            return bincode::deserialize_from(&mut reader).map_err(Into::into);
        }

        let mut data = vec![0u8; data_len as usize];
        read_parallel_at(&f, &mut data, 8)?;

        let mut trailer = [0u8; 12];
        read_exact_at_offset(&f, &mut trailer, 8 + data_len)?;
        let num_hashes = u32::from_le_bytes(trailer[0..4].try_into().unwrap());
        let num_bits = u64::from_le_bytes(trailer[4..12].try_into().unwrap());

        Ok(Self { data, num_hashes, num_bits })
    }
}

/// Size of each positioned read when filling a large buffer in parallel.
const BLOOM_READ_CHUNK: usize = 64 * 1024 * 1024;

/// Fills `buf` from `file` starting at `offset`, using parallel positioned reads when the buffer
/// is large enough to be worth it. Chunks are disjoint slices of `buf`, so no locking is needed.
#[cfg(unix)]
fn read_parallel_at(file: &File, buf: &mut [u8], offset: u64) -> Result<()> {
    use rayon::prelude::*;
    use std::os::unix::fs::FileExt;

    if buf.len() <= BLOOM_READ_CHUNK {
        return file.read_exact_at(buf, offset).map_err(Into::into);
    }

    buf.par_chunks_mut(BLOOM_READ_CHUNK)
        .enumerate()
        .try_for_each(|(i, chunk)| {
            file.read_exact_at(chunk, offset + (i * BLOOM_READ_CHUNK) as u64)
        })?;
    Ok(())
}

#[cfg(not(unix))]
fn read_parallel_at(file: &File, buf: &mut [u8], offset: u64) -> Result<()> {
    read_exact_at_offset(file, buf, offset)
}

fn read_exact_at_offset(file: &File, buf: &mut [u8], offset: u64) -> Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileExt;
        file.read_exact_at(buf, offset)?;
        return Ok(());
    }
    #[cfg(not(unix))]
    {
        let mut f = file.try_clone()?;
        f.seek(SeekFrom::Start(offset))?;
        f.read_exact(buf)?;
        Ok(())
    }
}

#[cfg(test)]
mod bloom_load_tests {
    use super::*;

    fn write_with_bincode(filter: &DomainBloomFilter, path: &std::path::Path) {
        let f = File::create(path).unwrap();
        let mut w = BufWriter::new(f);
        bincode::serialize_into(&mut w, filter).unwrap();
        w.flush().unwrap();
    }

    fn roundtrip(data: Vec<u8>, num_hashes: u32, num_bits: u64) {
        let dir = std::env::temp_dir().join(format!(
            "tdkc-bloom-{}-{}",
            std::process::id(),
            data.len()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.bloom");

        let original = DomainBloomFilter::new(data.clone(), num_hashes, num_bits);
        write_with_bincode(&original, &path);

        let loaded = DomainBloomFilter::load(path.to_str().unwrap()).unwrap();
        assert_eq!(loaded.data, data, "data mismatch for len {}", data.len());
        assert_eq!(loaded.num_hashes, num_hashes);
        assert_eq!(loaded.num_bits, num_bits);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn matches_bincode_for_small_filter() {
        roundtrip((0u8..255).cycle().take(4096).collect(), 13, 32768);
    }

    #[test]
    fn matches_bincode_for_empty_filter() {
        roundtrip(Vec::new(), 1, 0);
    }

    /// Crosses the parallel-read threshold so the chunked path and its final short chunk are
    /// both exercised.
    #[test]
    fn matches_bincode_across_chunk_boundary() {
        let len = BLOOM_READ_CHUNK + 12345;
        let data: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        let num_bits = (len as u64) * 8;
        roundtrip(data, 13, num_bits);
    }

    /// A filter loaded by the fast path must answer membership exactly as the in-memory one did.
    #[test]
    fn membership_survives_the_fast_load() {
        let mut data = vec![0u8; 1 << 16];
        let num_bits = (data.len() as u64) * 8;
        let num_hashes = 13;
        // Insert a few items the same way AtomicBloom does.
        let items: Vec<u64> = (0..64).map(|i| i * 2_654_435_761).collect();
        for &item in &items {
            let (h1, h2) = crate::hash::bloom_hash_pair(item);
            for i in 0..num_hashes {
                let bit = h1.wrapping_add((i as u64).wrapping_mul(h2)) % num_bits;
                data[(bit / 8) as usize] |= 1u8 << (bit % 8);
            }
        }

        let dir = std::env::temp_dir().join(format!("tdkc-bloom-mem-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("m.bloom");
        let original = DomainBloomFilter::new(data, num_hashes, num_bits);
        write_with_bincode(&original, &path);

        let loaded = DomainBloomFilter::load(path.to_str().unwrap()).unwrap();
        for &item in &items {
            assert!(loaded.contains(item), "inserted item {} went missing", item);
        }
        // And agree with the original on non-inserted probes, whatever the answer.
        for probe in 1_000_000u64..1_000_500 {
            assert_eq!(loaded.contains(probe), original.contains(probe));
        }

        std::fs::remove_dir_all(&dir).ok();
    }
}
