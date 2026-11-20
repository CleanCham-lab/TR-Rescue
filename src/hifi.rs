// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Indexed access to the aligned (ideally HiPhase-haplotagged) HiFi BAM:
//! depth profiles around a locus and the reads overlapping its flanks.

use std::fs::File;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow};
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;

use crate::scan::int_tag;

pub const TAG_HP: Tag = Tag::new(b'H', b'P');
pub const TAG_PS: Tag = Tag::new(b'P', b'S');

pub struct AlnBam {
    pub path: PathBuf,
    reader: bam::io::Reader<noodles::bgzf::io::Reader<File>>,
    pub header: sam::Header,
    index: bam::Index,
}

#[derive(Clone, Debug)]
pub struct AlnRead {
    pub name: String,
    pub seq: Vec<u8>,
    pub hp: Option<i64>,
    pub mapq: u8,
    pub ref_start: u64,
    pub ref_end: u64,
}

#[derive(Clone, Debug, Default)]
pub struct Depth {
    pub win_start: u64,
    pub bin: u64,
    pub total: Vec<f64>,
    pub hp1: Vec<f64>,
    pub hp2: Vec<f64>,
    pub unphased: Vec<f64>,
    pub read_lens: Vec<u64>,
    /// (ref start, ref end, hp 1/2/0, phase set or -1) of every counted read
    pub spans: Vec<(u64, u64, u8, i64)>,
    pub n_reads: usize,
    pub n_phased: usize,
}

fn find_index(path: &Path) -> Result<bam::Index> {
    match bam::fs::read_associated_index(path) {
        Ok(i) => Ok(i),
        Err(_) => {
            // samtools/pbmm2 sometimes write x.bai instead of x.bam.bai
            let alt = path.with_extension("bai");
            let i = bam::bai::fs::read(&alt)
                .with_context(|| format!("no BAM index for {} (.bam.bai / .bai / .csi)", path.display()))?;
            Ok(bam::Index::Bai(i))
        }
    }
}

fn ref_span(rec: &bam::Record) -> u64 {
    rec.cigar()
        .iter()
        .filter_map(|r| r.ok())
        .filter(|op| matches!(op.kind(), Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch))
        .map(|op| op.len() as u64)
        .sum()
}

impl AlnBam {
    pub fn open(path: &Path) -> Result<Self> {
        let mut reader = File::open(path)
            .map(bam::io::Reader::new)
            .with_context(|| format!("opening {}", path.display()))?;
        let header = reader.read_header()?;
        let index = find_index(path)?;
        Ok(Self { path: path.to_path_buf(), reader, header, index })
    }

    /// `chr9` vs `9` naming.
    pub fn resolve(&self, chrom: &str) -> Option<String> {
        let refs = self.header.reference_sequences();
        if refs.contains_key(chrom.as_bytes()) {
            return Some(chrom.to_string());
        }
        let alt = match chrom.strip_prefix("chr") {
            Some(s) => s.to_string(),
            None => format!("chr{chrom}"),
        };
        refs.contains_key(alt.as_bytes()).then_some(alt)
    }

    fn query<F: FnMut(&bam::Record)>(&mut self, chrom: &str, start: u64, end: u64, mut f: F) -> Result<()> {
        let region: noodles::core::Region = format!("{}:{}-{}", chrom, start.max(0) + 1, end)
            .parse()
            .map_err(|e| anyhow!("bad region {chrom}:{start}-{end}: {e:?}"))?;
        let q = self.reader.query(&self.header, &self.index, &region)?;
        for r in q.records() {
            let rec = r?;
            f(&rec);
        }
        Ok(())
    }

    /// Mean base depth (alignment span) per bin over [start, end), split by HP.
    pub fn depth(&mut self, chrom: &str, start: u64, end: u64, bin: u64, min_mapq: u8) -> Result<Depth> {
        let nb = ((end - start).div_ceil(bin)) as usize;
        let mut d = Depth {
            win_start: start,
            bin,
            total: vec![0.0; nb],
            hp1: vec![0.0; nb],
            hp2: vec![0.0; nb],
            unphased: vec![0.0; nb],
            ..Default::default()
        };
        self.query(chrom, start, end, |rec| {
            let fl = rec.flags();
            if fl.is_unmapped() || fl.is_secondary() || fl.is_supplementary() || fl.is_duplicate() || fl.is_qc_fail() {
                return;
            }
            let mq = rec.mapping_quality().map(|m| m.get()).unwrap_or(255);
            if mq != 255 && mq < min_mapq {
                return;
            }
            let Some(Ok(pos)) = rec.alignment_start() else { return };
            let s = pos.get() as u64 - 1;
            let e = s + ref_span(rec);
            let hp = int_tag(&rec.data(), &TAG_HP);
            d.n_reads += 1;
            if hp.is_some() {
                d.n_phased += 1;
            }
            d.read_lens.push(rec.sequence().len() as u64);
            let ps = int_tag(&rec.data(), &TAG_PS).unwrap_or(-1);
            d.spans.push((s, e, match hp { Some(1) => 1, Some(2) => 2, _ => 0 }, ps));
            let a = s.max(start);
            let b = e.min(end);
            if a >= b {
                return;
            }
            let mut p = a;
            while p < b {
                let bi = ((p - start) / bin) as usize;
                let bend = start + (bi as u64 + 1) * bin;
                let q = b.min(bend);
                let w = (q - p) as f64 / bin as f64;
                d.total[bi] += w;
                match hp {
                    Some(1) => d.hp1[bi] += w,
                    Some(2) => d.hp2[bi] += w,
                    _ => d.unphased[bi] += w,
                }
                p = q;
            }
        })?;
        Ok(d)
    }

    /// Primary reads overlapping [start, end) with their sequences.
    pub fn reads(&mut self, chrom: &str, start: u64, end: u64) -> Result<Vec<AlnRead>> {
        let mut out = Vec::new();
        self.query(chrom, start, end, |rec| {
            let fl = rec.flags();
            if fl.is_unmapped() || fl.is_secondary() || fl.is_supplementary() || fl.is_duplicate() || fl.is_qc_fail() {
                return;
            }
            let Some(Ok(pos)) = rec.alignment_start() else { return };
            let s = pos.get() as u64 - 1;
            out.push(AlnRead {
                name: rec.name().map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()).unwrap_or_default(),
                seq: rec.sequence().iter().collect(),
                hp: int_tag(&rec.data(), &TAG_HP),
                mapq: rec.mapping_quality().map(|m| m.get()).unwrap_or(255),
                ref_start: s,
                ref_end: s + ref_span(rec),
            });
        })?;
        Ok(out)
    }
}
