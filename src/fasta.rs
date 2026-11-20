// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Minimal indexed FASTA access (.fai) with uppercase output.

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};

#[derive(Clone, Debug)]
pub struct FaiEntry {
    pub name: String,
    pub len: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

pub struct Fasta {
    pub path: PathBuf,
    pub entries: Vec<FaiEntry>,
    idx: HashMap<String, usize>,
}

impl Fasta {
    pub fn open(path: &Path) -> Result<Self> {
        let mut fai = path.as_os_str().to_owned();
        fai.push(".fai");
        let fai = PathBuf::from(fai);
        let text = std::fs::read_to_string(&fai)
            .with_context(|| format!("FASTA index not found: {} (run `samtools faidx`)", fai.display()))?;
        let mut entries = Vec::new();
        let mut idx = HashMap::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            let e = FaiEntry {
                name: f[0].to_string(),
                len: f[1].parse()?,
                offset: f[2].parse()?,
                line_bases: f[3].parse()?,
                line_width: f[4].parse()?,
            };
            idx.insert(e.name.clone(), entries.len());
            entries.push(e);
        }
        if entries.is_empty() {
            bail!("empty FASTA index {}", fai.display());
        }
        Ok(Self { path: path.to_path_buf(), entries, idx })
    }

    pub fn entry(&self, chrom: &str) -> Option<&FaiEntry> {
        self.idx.get(chrom).map(|&i| &self.entries[i])
    }

    /// Resolve `chr9` vs `9` naming differences.
    pub fn resolve(&self, chrom: &str) -> Option<String> {
        if self.idx.contains_key(chrom) {
            return Some(chrom.to_string());
        }
        let alt = match chrom.strip_prefix("chr") {
            Some(s) => s.to_string(),
            None => format!("chr{chrom}"),
        };
        self.idx.contains_key(&alt).then_some(alt)
    }

    /// Fetch [start, end) of `chrom` as uppercase bytes (clamped to the contig).
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let e = self.entry(chrom).ok_or_else(|| anyhow!("contig {chrom} not in reference"))?;
        let mut f = File::open(&self.path)?;
        read_range(&mut f, e, start, end)
    }
}

pub fn read_range(f: &mut File, e: &FaiEntry, start: u64, end: u64) -> Result<Vec<u8>> {
    let end = end.min(e.len);
    if start >= end {
        return Ok(Vec::new());
    }
    let off = |p: u64| e.offset + (p / e.line_bases) * e.line_width + (p % e.line_bases);
    let a = off(start);
    let b = off(end - 1) + 1;
    let mut buf = vec![0u8; (b - a) as usize];
    f.seek(SeekFrom::Start(a))?;
    f.read_exact(&mut buf)?;
    let mut out = Vec::with_capacity((end - start) as usize);
    for &c in &buf {
        if c != b'\n' && c != b'\r' {
            out.push(c.to_ascii_uppercase());
        }
    }
    Ok(out)
}
