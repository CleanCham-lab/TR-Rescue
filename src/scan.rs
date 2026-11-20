// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Streaming scan of unaligned PacBio BAMs (fail_reads, hifi_reads) for reads that carry
//! panel k-mers. Decompression and k-mer checks run on separate threads.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::Read;
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use anyhow::{Context, Result};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::anchor::{Call, MotifSet, analyze};
use crate::catalog::Locus;
use crate::kmer::{Panel, decode_packed, locus_hits_packed};

pub const TAG_FF: Tag = Tag::new(b'f', b'f');
pub const TAG_RQ: Tag = Tag::new(b'r', b'q');
pub const TAG_NP: Tag = Tag::new(b'n', b'p');

/// Primary fail category, same priority as failscope.
pub fn ff_category(ff: Option<i64>) -> &'static str {
    match ff {
        None | Some(0) => "fail_unflagged",
        Some(x) if x & 0x02 != 0 => "fail_control",
        Some(x) if x & (0x10 | 0x20 | 0x40) != 0 => "fail_adapter",
        Some(x) if x & 0x80 != 0 => "fail_segment",
        Some(x) if x & 0x08 != 0 => "fail_subread",
        Some(x) if x & 0x04 != 0 => "fail_single_strand",
        Some(x) if x & 0x01 != 0 => "fail_lowqv",
        _ => "fail_unflagged",
    }
}

#[derive(Clone, Debug, Default)]
pub struct CatStat {
    pub reads: u64,
    pub bases: u64,
}

#[derive(Clone, Debug, Default)]
pub struct ScanStats {
    pub reads: u64,
    pub bases: u64,
    /// reads passing the k-mer prefilter
    pub prefilter: u64,
    /// reads placed on at least one locus
    pub candidates: u64,
    pub seconds: f64,
    pub by_cat: BTreeMap<String, CatStat>,
}

/// A read that was placed on at least one locus (sequence not kept).
#[derive(Clone, Debug)]
pub struct Candidate {
    pub name: String,
    pub len: usize,
    pub ff: Option<i64>,
    pub rq: Option<f32>,
    pub np: Option<i64>,
    pub calls: Vec<Call>,
}

pub struct ScanCtx<'a> {
    pub panel: &'a Panel,
    pub loci: &'a [Locus],
    pub motif_sets: &'a [MotifSet],
    pub min_hits: usize,
}

/// Counts compressed bytes read, for progress reporting.
struct Counting<R> {
    inner: R,
    n: Arc<AtomicU64>,
}

impl<R: Read> Read for Counting<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let k = self.inner.read(buf)?;
        self.n.fetch_add(k as u64, Ordering::Relaxed);
        Ok(k)
    }
}

fn hms(s: f64) -> String {
    let s = s.max(0.0) as u64;
    if s >= 3600 {
        format!("{}h{:02}m", s / 3600, (s % 3600) / 60)
    } else if s >= 60 {
        format!("{}m{:02}s", s / 60, s % 60)
    } else {
        format!("{s}s")
    }
}

struct Raw {
    name: String,
    packed: Vec<u8>,
    len: usize,
    ff: Option<i64>,
    rq: Option<f32>,
    np: Option<i64>,
}

pub fn int_tag(data: &bam::record::Data<'_>, tag: &Tag) -> Option<i64> {
    data.get(tag).and_then(|r| r.ok()).and_then(|v| v.as_int())
}

pub fn float_tag(data: &bam::record::Data<'_>, tag: &Tag) -> Option<f32> {
    match data.get(tag) {
        Some(Ok(Value::Float(x))) => Some(x),
        _ => None,
    }
}

/// Scan one unaligned BAM. `is_fail` controls how reads are categorised in the stats.
pub fn scan_bam(path: &Path, ctx: &ScanCtx<'_>, threads: usize, is_fail: bool) -> Result<(ScanStats, Vec<Candidate>)> {
    let panel = ctx.panel;
    let t0 = Instant::now();
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let total_bytes = file.metadata().map(|m| m.len()).unwrap_or(0);
    let counter = Arc::new(AtomicU64::new(0));
    let file = Counting { inner: file, n: counter.clone() };
    let workers = NonZero::new(threads.max(1)).unwrap();
    let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(workers, file));
    eprintln!("    {:.1} GB compressed; progress every 20 s", total_bytes as f64 / 1e9);
    let mut last_report = Instant::now();
    reader.read_header().with_context(|| format!("reading header of {}", path.display()))?;

    let (tx, rx) = crossbeam_channel::bounded::<Vec<Raw>>(threads.max(1) * 4);
    let ctx_ref = ctx;
    let pre = AtomicU64::new(0);
    let pre_ref = &pre;
    let (ctx, crx) = crossbeam_channel::unbounded::<Candidate>();
    let k = panel.k;
    let mut stats = ScanStats::default();

    let res: Result<()> = std::thread::scope(|sc| {
        for _ in 0..threads.max(1) {
            let rx = rx.clone();
            let ctx = ctx.clone();
            let sctx = &*ctx_ref;
            sc.spawn(move || {
                for batch in rx {
                    for r in batch {
                        if r.len >= k && locus_hits_packed(&r.packed, r.len, panel, sctx.min_hits as u16) {
                            pre_ref.fetch_add(1, Ordering::Relaxed);
                            let seq = decode_packed(&r.packed, r.len);
                            let calls = analyze(&seq, panel, sctx.loci, sctx.motif_sets, sctx.min_hits);
                            if calls.is_empty() {
                                continue;
                            }
                            let _ = ctx.send(Candidate { len: r.len, name: r.name, ff: r.ff, rq: r.rq, np: r.np, calls });
                        }
                    }
                }
            });
        }
        drop(ctx);
        drop(rx);

        let mut rec = bam::Record::default();
        let mut batch: Vec<Raw> = Vec::with_capacity(512);
        let mut batch_bases = 0usize;
        loop {
            let n = reader.read_record(&mut rec).with_context(|| format!("reading {}", path.display()))?;
            if n == 0 {
                break;
            }
            let fl = rec.flags();
            if fl.is_secondary() || fl.is_supplementary() {
                continue;
            }
            let seq = rec.sequence();
            let len = seq.len();
            let data = rec.data();
            let ff = int_tag(&data, &TAG_FF);
            let cat = if is_fail { ff_category(ff) } else { "hifi" };
            stats.reads += 1;
            stats.bases += len as u64;
            let cs = stats.by_cat.entry(cat.to_string()).or_default();
            cs.reads += 1;
            cs.bases += len as u64;
            batch.push(Raw {
                name: rec.name().map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()).unwrap_or_default(),
                packed: seq.as_bytes().to_vec(),
                len,
                ff,
                rq: float_tag(&data, &TAG_RQ),
                np: int_tag(&data, &TAG_NP),
            });
            batch_bases += len;
            if last_report.elapsed().as_secs() >= 20 {
                last_report = Instant::now();
                let done = counter.load(Ordering::Relaxed) as f64;
                let el = t0.elapsed().as_secs_f64();
                let rate = done / el.max(1e-3);
                let pct = if total_bytes > 0 { 100.0 * done / total_bytes as f64 } else { 0.0 };
                let eta = if rate > 0.0 && total_bytes > 0 { (total_bytes as f64 - done) / rate } else { 0.0 };
                eprintln!(
                    "    {:5.1}%  {:.1}/{:.1} GB  {:.0} MB/s  {:.2}M reads  ETA {}",
                    pct,
                    done / 1e9,
                    total_bytes as f64 / 1e9,
                    rate / 1e6,
                    stats.reads as f64 / 1e6,
                    hms(eta)
                );
            }
            if batch_bases >= 4_000_000 {
                let b = std::mem::replace(&mut batch, Vec::with_capacity(512));
                batch_bases = 0;
                if tx.send(b).is_err() {
                    break;
                }
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        drop(tx);
        Ok(())
    });
    res?;
    let cands: Vec<Candidate> = crx.into_iter().collect();
    stats.candidates = cands.len() as u64;
    stats.prefilter = pre.load(Ordering::Relaxed);
    stats.seconds = t0.elapsed().as_secs_f64();
    Ok((stats, cands))
}
