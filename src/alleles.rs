// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! HiFi allele summary per sample x locus: up to two length clusters from spanning reads.

use serde::Serialize;

use crate::anchor::{Bound, Cls};
use crate::catalog::Locus;
use crate::i18n::Tr;
use crate::report::ReadRow;

#[derive(Clone, Debug, Serialize)]
pub struct Cluster {
    pub n: usize,
    pub median_bp: i64,
    pub min_bp: i64,
    pub max_bp: i64,
    pub units: f64,
    pub min_units: f64,
    pub max_units: f64,
    pub hp1: usize,
    pub hp2: usize,
    pub motif: String,
    /// normal / intermediate / pathogenic / benign_motif
    pub class: &'static str,
    /// wide size spread (somatic mosaicism or noise)
    pub spread: bool,
}

#[derive(Clone, Debug, Serialize, Default)]
pub struct AlleleSummary {
    pub clusters: Vec<Cluster>,
    /// largest HiFi lower bound from reads ending inside the repeat (bp)
    pub hifi_lb_max: Option<i64>,
    /// fail-read expanded calls (spanning or lower bound) with pure repeat
    pub fail_expanded: usize,
    pub fail_expanded_max_bp: Option<i64>,
    /// most severe cluster class
    pub range: &'static str,
}

fn class_of(l: &Locus, units: f64, motif_pathogenic: bool) -> &'static str {
    let units = units.round();
    if l.contraction {
        return if l.path_max_units.is_some_and(|p| units <= p) { "pathogenic" } else { "normal" };
    }
    if units <= l.normal_max_units {
        "normal"
    } else if !motif_pathogenic {
        "benign_motif"
    } else if units >= l.path_min_units {
        "pathogenic"
    } else {
        "intermediate"
    }
}

fn severity(c: &str) -> u8 {
    match c {
        "pathogenic" => 4,
        "intermediate" => 3,
        "benign_motif" => 2,
        "normal" => 1,
        _ => 0,
    }
}

fn median(v: &[i64]) -> i64 {
    let mut v = v.to_vec();
    v.sort_unstable();
    v[v.len() / 2]
}

pub fn summarize(rows: &[&ReadRow], l: &Locus) -> AlleleSummary {
    let unit = l.unit_len() as f64;
    let mut sp: Vec<&ReadRow> = rows
        .iter()
        .copied()
        .filter(|r| r.source == "hifi_aligned" && r.call.cls == Cls::Spanning && r.call.bound == Bound::Exact && r.call.state != "NONREPEAT")
        .collect();
    sp.sort_by_key(|r| r.call.allele_bp.unwrap_or(0));
    let mut out = AlleleSummary::default();
    if !sp.is_empty() {
        let bps: Vec<i64> = sp.iter().map(|r| r.call.allele_bp.unwrap_or(0)).collect();
        // split at the largest relative gap, if it is clearly a second allele
        let mut cut = None;
        let mut best = 0.0;
        for i in 1..bps.len() {
            let gap = (bps[i] - bps[i - 1]) as f64;
            let rel = gap / (bps[i - 1].max(1) as f64);
            if gap >= 2.0 * unit && rel >= 0.12 && gap > best {
                best = gap;
                cut = Some(i);
            }
        }
        let groups: Vec<&[&ReadRow]> = match cut {
            Some(i) => vec![&sp[..i], &sp[i..]],
            None => vec![&sp[..]],
        };
        for g in groups {
            let v: Vec<i64> = g.iter().map(|r| r.call.allele_bp.unwrap_or(0)).collect();
            let med = median(&v);
            let (mn, mx) = (*v.iter().min().unwrap(), *v.iter().max().unwrap());
            let mut motifs: std::collections::BTreeMap<&str, usize> = Default::default();
            for r in g.iter() {
                *motifs.entry(r.call.dom_motif.as_str()).or_insert(0) += 1;
            }
            let motif = motifs.into_iter().max_by_key(|(_, c)| *c).map(|(m, _)| m.to_string()).unwrap_or_default();
            let path = motif == "-" || l.motif_is_pathogenic(&motif);
            let units = l.units(med);
            out.clusters.push(Cluster {
                n: g.len(),
                median_bp: med,
                min_bp: mn,
                max_bp: mx,
                units,
                min_units: l.units(mn),
                max_units: l.units(mx),
                hp1: g.iter().filter(|r| r.hp == Some(1)).count(),
                hp2: g.iter().filter(|r| r.hp == Some(2)).count(),
                motif,
                class: class_of(l, units, path),
                spread: (mx - mn) as f64 > (0.2 * med as f64).max(10.0 * unit),
            });
        }
    }
    // HiFi reads that end inside an expanded, pure, pathogenic-motif repeat (state already checks motif and purity)
    let mut lbs: Vec<i64> = rows
        .iter()
        .filter(|r| r.source == "hifi_aligned" && r.call.bound == Bound::Lower && r.call.state == "EXPANDED_LB")
        .filter_map(|r| r.call.allele_bp)
        .collect();
    lbs.sort_unstable_by(|a, b| b.cmp(a));
    out.hifi_lb_max = lbs.first().copied();
    // expansions seen only outside the aligned HiFi BAM: raw subreads carry ~10% indels, so a read counts only
    // if it is clearly beyond both the pathogenic threshold and every HiFi allele of this sample
    let hifi_max = out.clusters.iter().map(|c| c.max_bp).max().unwrap_or(0) as f64;
    let min_bp = (l.path_min_bp() + 3 * l.unit_len() as i64).max((1.5 * hifi_max).ceil() as i64);
    let fe: Vec<i64> = rows
        .iter()
        .filter(|r| (r.source.starts_with("fail") || r.source.starts_with("unassigned")) && matches!(r.call.state, "EXPANDED" | "EXPANDED_LB"))
        .filter_map(|r| r.call.allele_bp)
        .filter(|&b| b >= min_bp)
        .collect();
    out.fail_expanded = fe.len();
    out.fail_expanded_max_bp = fe.iter().copied().max();
    let mut range = out.clusters.iter().map(|c| c.class).max_by_key(|c| severity(c)).unwrap_or("NA");
    // HiFi lower bounds beyond the pathogenic threshold count too, when at least two reads agree
    if let Some(&lb2) = lbs.get(1) {
        if !l.contraction && l.units(lb2).round() >= l.path_min_units && severity(range) < 4 {
            range = "pathogenic";
        }
    }
    out.range = range;
    out
}

/// "84 bp (28회, n=23, HP2) / 430 bp (143회, n=16, HP1)"
pub fn text(a: &AlleleSummary) -> Tr {
    if a.clusters.is_empty() {
        return Tr::new("HiFi spanning read 없음", "no spanning HiFi read");
    }
    let one = |en: bool| {
        a.clusters
            .iter()
            .map(|c| {
                let hp = match (c.hp1, c.hp2) {
                    (0, 0) => String::new(),
                    (x, y) if x >= 2 * y.max(1) => ", HP1".into(),
                    (x, y) if y >= 2 * x.max(1) => ", HP2".into(),
                    _ => if en { ", mixed HP".into() } else { ", HP 혼합".into() },
                };
                let u = if en { "" } else { "회" };
                if c.spread {
                    format!("{}–{} bp ({:.0}–{:.0}{u}, n={}{})", c.min_bp, c.max_bp, c.min_units, c.max_units, c.n, hp)
                } else {
                    format!("{} bp ({:.0}{u}, n={}{})", c.median_bp, c.units, c.n, hp)
                }
            })
            .collect::<Vec<_>>()
            .join(" / ")
    };
    Tr::new(one(false), one(true))
}
