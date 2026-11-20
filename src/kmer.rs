// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Canonical k-mer encoding and the flank k-mer panel.
//!
//! The panel holds k-mers from the unique flanks of each locus. A k-mer is kept only if it
//! occurs exactly once (either strand) in the primary reference contigs, so a hit places a
//! read on that flank without alignment. This is what lets unaligned fail BAMs be scanned.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::Result;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::catalog::Locus;
use crate::fasta::{Fasta, read_range};

/// ASCII -> 2-bit (A0 C1 G2 T3), 4 = invalid.
pub static ENC: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

/// BAM 4-bit code (=ACMGRSVTWYHKDBN) -> 2-bit, 4 = invalid.
pub static NIB: [u8; 16] = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];
pub static NIB_ASCII: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

#[inline]
fn mask(k: usize) -> u64 {
    if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 }
}

/// Call `f(pos, canonical, fwd_is_canonical)` for every valid k-mer of an ASCII sequence.
#[inline]
pub fn for_each_kmer<F: FnMut(usize, u64, bool)>(seq: &[u8], k: usize, mut f: F) {
    let m = mask(k);
    let shift = 2 * (k as u64 - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = ENC[b as usize];
        if c > 3 {
            valid = 0;
            continue;
        }
        let c = c as u64;
        fw = ((fw << 2) | c) & m;
        rv = (rv >> 2) | ((3 - c) << shift);
        valid += 1;
        if valid >= k {
            if fw < rv { f(i + 1 - k, fw, true) } else { f(i + 1 - k, rv, false) }
        }
    }
}

/// Same as [`for_each_kmer`] over a BAM 4-bit packed sequence; stops when `f` returns true.
#[inline]
pub fn any_kmer_packed<F: FnMut(u64) -> bool>(packed: &[u8], len: usize, k: usize, mut f: F) -> bool {
    let m = mask(k);
    let shift = 2 * (k as u64 - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for i in 0..len {
        let byte = packed[i >> 1];
        let code = if i & 1 == 0 { byte >> 4 } else { byte & 0x0f };
        let c = NIB[code as usize];
        if c > 3 {
            valid = 0;
            continue;
        }
        let c = c as u64;
        fw = ((fw << 2) | c) & m;
        rv = (rv >> 2) | ((3 - c) << shift);
        valid += 1;
        if valid >= k && f(fw.min(rv)) {
            return true;
        }
    }
    false
}

/// True if some locus collects at least `need` panel hits in a BAM 4-bit packed sequence.
pub fn locus_hits_packed(packed: &[u8], len: usize, panel: &Panel, need: u16) -> bool {
    let k = panel.k;
    let m = mask(k);
    let shift = 2 * (k as u64 - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut counts: Vec<(u32, u16)> = Vec::new();
    for i in 0..len {
        let byte = packed[i >> 1];
        let code = if i & 1 == 0 { byte >> 4 } else { byte & 0x0f };
        let c = NIB[code as usize];
        if c > 3 {
            valid = 0;
            continue;
        }
        let c = c as u64;
        fw = ((fw << 2) | c) & m;
        rv = (rv >> 2) | ((3 - c) << shift);
        valid += 1;
        if valid >= k {
            for e in panel.get(fw.min(rv)) {
                match counts.iter_mut().find(|x| x.0 == e.locus) {
                    Some(x) => {
                        x.1 += 1;
                        if x.1 >= need {
                            return true;
                        }
                    }
                    None => counts.push((e.locus, 1)),
                }
            }
        }
    }
    false
}

pub fn decode_packed(packed: &[u8], len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| {
            let byte = packed[i >> 1];
            let code = if i & 1 == 0 { byte >> 4 } else { byte & 0x0f };
            NIB_ASCII[code as usize]
        })
        .collect()
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            x => x,
        })
        .collect()
}

#[derive(Clone, Copy, Debug)]
pub struct Entry {
    pub locus: u32,
    /// 0 = left flank, 1 = right flank
    pub side: u8,
    /// 0-based reference position of the k-mer's first base
    pub refpos: u32,
    /// true if the reference-forward k-mer is the canonical one
    pub canon_is_fwd: bool,
}

#[derive(Clone, Debug, Default)]
pub struct SideStat {
    pub locus: usize,
    pub side: u8,
    pub region_start: u64,
    pub region_end: u64,
    pub kmers: usize,
    pub unique: usize,
    /// unique (primary) k-mers that also occur on ALT/fix contigs
    pub alt_shared: usize,
}

const BITS_LOG2: u32 = 27;

/// k-mer -> entries. A k-mer can belong to several loci when their flank windows overlap
/// (e.g. neighbouring loci in one gene); it is still one reference position.
pub struct Panel {
    pub k: usize,
    pub flank: u64,
    pub pad: u64,
    map: FxHashMap<u64, (u32, u16)>,
    entries: Vec<Entry>,
    bits: Vec<u64>,
    pub stats: Vec<SideStat>,
}

#[inline]
fn bit_idx(x: u64) -> usize {
    (x.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - BITS_LOG2)) as usize
}

pub fn flank_regions(l: &Locus, flank: u64, pad: u64) -> [(u64, u64); 2] {
    let ls = l.start.saturating_sub(pad + flank);
    let le = l.start.saturating_sub(pad);
    let rs = l.end + pad;
    let re = l.end + pad + flank;
    [(ls, le), (rs, re)]
}

pub fn is_primary_contig(name: &str) -> bool {
    !(name.contains("_alt") || name.contains("_fix") || name.contains("_hap"))
}

impl Panel {
    #[inline]
    pub fn get(&self, canon: u64) -> &[Entry] {
        let i = bit_idx(canon);
        if self.bits[i >> 6] & (1u64 << (i & 63)) == 0 {
            return &[];
        }
        match self.map.get(&canon) {
            Some(&(s, n)) => &self.entries[s as usize..s as usize + n as usize],
            None => &[],
        }
    }

    #[inline]
    pub fn contains(&self, canon: u64) -> bool {
        !self.get(canon).is_empty()
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    fn finish(k: usize, flank: u64, pad: u64, mut pairs: Vec<(u64, Entry)>, stats: Vec<SideStat>) -> Self {
        pairs.sort_by_key(|p| (p.0, p.1.locus, p.1.side));
        let mut map = FxHashMap::default();
        let mut entries = Vec::with_capacity(pairs.len());
        let mut bits = vec![0u64; (1usize << BITS_LOG2) / 64];
        let mut i = 0;
        while i < pairs.len() {
            let key = pairs[i].0;
            let start = entries.len();
            while i < pairs.len() && pairs[i].0 == key {
                entries.push(pairs[i].1);
                i += 1;
            }
            map.insert(key, (start as u32, (entries.len() - start) as u16));
            let b = bit_idx(key);
            bits[b >> 6] |= 1u64 << (b & 63);
        }
        Panel { k, flank, pad, map, entries, bits, stats }
    }

    /// Build from the reference: collect flank k-mers, then drop any k-mer that occurs at more
    /// reference positions (either strand, primary contigs) than the flank positions it came from.
    pub fn build(fa: &Fasta, loci: &[Locus], k: usize, flank: u64, pad: u64, threads: usize) -> Result<Self> {
        let mut pairs: Vec<(u64, Entry)> = Vec::new();
        let mut stats = Vec::new();
        for (li, l) in loci.iter().enumerate() {
            let chrom = fa.resolve(&l.chrom).unwrap_or_else(|| l.chrom.clone());
            for (side, (s, e)) in flank_regions(l, flank, pad).into_iter().enumerate() {
                let seq = fa.fetch(&chrom, s, e)?;
                let mut local: FxHashMap<u64, Option<Entry>> = FxHashMap::default();
                let mut n = 0;
                for_each_kmer(&seq, k, |off, canon, fwd| {
                    n += 1;
                    let e = Entry { locus: li as u32, side: side as u8, refpos: (s + off as u64) as u32, canon_is_fwd: fwd };
                    local.entry(canon).and_modify(|v| *v = None).or_insert(Some(e));
                });
                // k-mers repeated within one flank are not informative
                pairs.extend(local.into_iter().filter_map(|(key, v)| v.map(|e| (key, e))));
                stats.push(SideStat {
                    locus: li,
                    side: side as u8,
                    region_start: s,
                    region_end: e,
                    kmers: n,
                    unique: 0,
                    alt_shared: 0,
                });
            }
        }
        // expected reference occurrences = distinct (chrom, refpos) among a k-mer's entries
        let chrom_ids: Vec<u32> = loci
            .iter()
            .map(|l| {
                let name = fa.resolve(&l.chrom).unwrap_or_else(|| l.chrom.clone());
                fa.entries.iter().position(|x| x.name == name).unwrap_or(0) as u32
            })
            .collect();
        // distinct reference positions per k-mer (almost always 1)
        pairs.sort_unstable_by_key(|(key, e)| (*key, chrom_ids[e.locus as usize], e.refpos));
        let mut keys: FxHashSet<u64> = FxHashSet::default();
        let mut multi: FxHashMap<u64, u32> = FxHashMap::default();
        let mut i = 0;
        while i < pairs.len() {
            let key = pairs[i].0;
            let mut distinct = 0u32;
            let mut last: Option<(u32, u32)> = None;
            while i < pairs.len() && pairs[i].0 == key {
                let pos = (chrom_ids[pairs[i].1.locus as usize], pairs[i].1.refpos);
                if last != Some(pos) {
                    distinct += 1;
                    last = Some(pos);
                }
                i += 1;
            }
            keys.insert(key);
            if distinct > 1 {
                multi.insert(key, distinct);
            }
        }
        let (counts, alt) = count_in_reference(fa, &keys, k, threads)?;
        drop(keys);
        pairs.retain(|(key, _)| counts.get(key).copied().unwrap_or(0) <= multi.get(key).copied().unwrap_or(1));
        let mut acc: FxHashMap<(usize, u8), (usize, usize)> = FxHashMap::default();
        for (key, e) in &pairs {
            let a = acc.entry((e.locus as usize, e.side)).or_default();
            a.0 += 1;
            if alt.contains_key(key) {
                a.1 += 1;
            }
        }
        for st in stats.iter_mut() {
            let (u, a) = acc.get(&(st.locus, st.side)).copied().unwrap_or((0, 0));
            st.unique = u;
            st.alt_shared = a;
        }
        Ok(Self::finish(k, flank, pad, pairs, stats))
    }

    pub fn save(&self, path: &Path, meta: &str) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "#meta\t{meta}")?;
        for s in &self.stats {
            writeln!(
                w,
                "#side\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.locus, s.side, s.region_start, s.region_end, s.kmers, s.unique, s.alt_shared
            )?;
        }
        for (key, &(st, n)) in &self.map {
            for e in &self.entries[st as usize..st as usize + n as usize] {
                writeln!(w, "{}\t{}\t{}\t{}\t{}", key, e.locus, e.side, e.refpos, e.canon_is_fwd as u8)?;
            }
        }
        Ok(())
    }

    pub fn load(path: &Path, meta: &str, k: usize, flank: u64, pad: u64) -> Option<Self> {
        let f = File::open(path).ok()?;
        let mut pairs = Vec::new();
        let mut stats = Vec::new();
        let mut ok = false;
        for line in BufReader::new(f).lines() {
            let line = line.ok()?;
            let f: Vec<&str> = line.split('\t').collect();
            if f[0] == "#meta" {
                ok = f.get(1).map(|m| *m == meta).unwrap_or(false);
                if !ok {
                    return None;
                }
            } else if f[0] == "#side" {
                stats.push(SideStat {
                    locus: f[1].parse().ok()?,
                    side: f[2].parse().ok()?,
                    region_start: f[3].parse().ok()?,
                    region_end: f[4].parse().ok()?,
                    kmers: f[5].parse().ok()?,
                    unique: f[6].parse().ok()?,
                    alt_shared: f.get(7).and_then(|x| x.parse().ok()).unwrap_or(0),
                });
            } else if f.len() == 5 {
                pairs.push((
                    f[0].parse().ok()?,
                    Entry { locus: f[1].parse().ok()?, side: f[2].parse().ok()?, refpos: f[3].parse().ok()?, canon_is_fwd: f[4] == "1" },
                ));
            }
        }
        ok.then(|| Self::finish(k, flank, pad, pairs, stats))
    }
}

type Counts = FxHashMap<u64, u32>;

/// Count occurrences of panel k-mers in parallel chunks: (primary contigs, ALT/fix contigs).
fn count_in_reference(fa: &Fasta, raw: &FxHashSet<u64>, k: usize, threads: usize) -> Result<(Counts, Counts)> {
    const CHUNK: u64 = 8_000_000;
    let mut jobs: Vec<(usize, u64, u64)> = Vec::new();
    for (ei, e) in fa.entries.iter().enumerate() {
        let mut s = 0;
        while s < e.len {
            let t = (s + CHUNK).min(e.len);
            jobs.push((ei, s, t));
            s = t;
        }
    }
    // A cheap membership filter in front of the hash map.
    let mut bits = vec![0u64; (1usize << BITS_LOG2) / 64];
    for &key in raw.iter() {
        let i = bit_idx(key);
        bits[i >> 6] |= 1u64 << (i & 63);
    }
    let next = AtomicUsize::new(0);
    let total: Mutex<(Counts, Counts)> = Mutex::new((FxHashMap::default(), FxHashMap::default()));
    let err: Mutex<Option<anyhow::Error>> = Mutex::new(None);
    std::thread::scope(|sc| {
        for _ in 0..threads.max(1) {
            sc.spawn(|| {
                let mut local: Counts = FxHashMap::default();
                let mut local_alt: Counts = FxHashMap::default();
                let mut f = match File::open(&fa.path) {
                    Ok(f) => f,
                    Err(e) => {
                        *err.lock().unwrap() = Some(e.into());
                        return;
                    }
                };
                loop {
                    let j = next.fetch_add(1, Ordering::Relaxed);
                    if j >= jobs.len() {
                        break;
                    }
                    let (ei, s, t) = jobs[j];
                    let e = &fa.entries[ei];
                    let primary = is_primary_contig(&e.name);
                    // extend by k-1 so k-mers starting in [s, t) are all seen exactly once
                    let seq = match read_range(&mut f, e, s, t + k as u64 - 1) {
                        Ok(v) => v,
                        Err(x) => {
                            *err.lock().unwrap() = Some(x);
                            return;
                        }
                    };
                    for_each_kmer(&seq, k, |_, canon, _| {
                        let i = bit_idx(canon);
                        if bits[i >> 6] & (1u64 << (i & 63)) != 0 && raw.contains(&canon) {
                            let m = if primary { &mut local } else { &mut local_alt };
                            *m.entry(canon).or_insert(0) += 1;
                        }
                    });
                }
                let mut t = total.lock().unwrap();
                for (key, c) in local {
                    *t.0.entry(key).or_insert(0) += c;
                }
                for (key, c) in local_alt {
                    *t.1.entry(key).or_insert(0) += c;
                }
            });
        }
    });
    if let Some(e) = err.into_inner().unwrap() {
        return Err(e);
    }
    Ok(total.into_inner().unwrap())
}
