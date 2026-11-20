// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Place a read on a locus from flank k-mer hits and measure the repeat it carries.
//!
//! Works the same for aligned HiFi reads and unaligned (noisy) fail reads: only the read
//! sequence and the unique flank k-mers are used. Hits are grouped per locus and strand;
//! each flank gives an "inner edge" position in read coordinates, and the repeat sits
//! between the edges (spanning) or runs from one edge (anchored).

use rustc_hash::{FxHashMap, FxHashSet};

use crate::catalog::Locus;
use crate::kmer::{ENC, Panel, for_each_kmer, revcomp};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cls {
    Spanning,
    LeftAnchor,
    RightAnchor,
    FlankLeft,
    FlankRight,
    Ambiguous,
}

impl Cls {
    pub fn parse(s: &str) -> Cls {
        match s {
            "SPANNING" => Cls::Spanning,
            "LEFT_ANCHOR" => Cls::LeftAnchor,
            "RIGHT_ANCHOR" => Cls::RightAnchor,
            "FLANK_LEFT" => Cls::FlankLeft,
            "FLANK_RIGHT" => Cls::FlankRight,
            _ => Cls::Ambiguous,
        }
    }
    pub fn as_str(&self) -> &'static str {
        match self {
            Cls::Spanning => "SPANNING",
            Cls::LeftAnchor => "LEFT_ANCHOR",
            Cls::RightAnchor => "RIGHT_ANCHOR",
            Cls::FlankLeft => "FLANK_LEFT",
            Cls::FlankRight => "FLANK_RIGHT",
            Cls::Ambiguous => "AMBIGUOUS",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Bound {
    Exact,
    Lower,
    Estimate,
    NA,
}

impl Bound {
    pub fn parse(s: &str) -> Bound {
        match s {
            "exact" => Bound::Exact,
            "lower_bound" => Bound::Lower,
            "estimate" => Bound::Estimate,
            _ => Bound::NA,
        }
    }
    pub fn as_str(&self) -> &'static str {
        match self {
            Bound::Exact => "exact",
            Bound::Lower => "lower_bound",
            Bound::Estimate => "estimate",
            Bound::NA => "NA",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Call {
    pub locus: usize,
    pub strand: char,
    pub cls: Cls,
    pub allele_bp: Option<i64>,
    pub bound: Bound,
    pub state: &'static str,
    pub dom_motif: String,
    pub dom_frac: f64,
    pub purity: f64,
    pub path_motif: bool,
    pub hits_l: usize,
    pub hits_r: usize,
    pub flags: Vec<&'static str>,
}

/// Motif unit sets per locus (all rotations, both strands), encoded 2-bit by length.
pub struct MotifSet {
    pub names: Vec<String>,
    pub lens: Vec<usize>,
    pub units: Vec<FxHashSet<u64>>,
}

fn encode(s: &[u8]) -> Option<u64> {
    let mut x = 0u64;
    for &b in s {
        let c = ENC[b as usize];
        if c > 3 {
            return None;
        }
        x = (x << 2) | c as u64;
    }
    Some(x)
}

impl MotifSet {
    pub fn new(motifs: &[String]) -> Self {
        let mut names = Vec::new();
        let mut lens = Vec::new();
        let mut units = Vec::new();
        for name in motifs {
            let name = name.to_ascii_uppercase();
            if name.is_empty() || name.len() > 31 {
                continue;
            }
            // degenerate motifs (GCN for polyalanine loci) are expanded to all concrete motifs
            let concrete = crate::catalog::expand_motif(&name);
            let mut set = FxHashSet::default();
            for c in &concrete {
                let m = c.as_bytes().to_vec();
                for src in [m.clone(), revcomp(&m)] {
                    for r in 0..src.len() {
                        let rot: Vec<u8> = src[r..].iter().chain(src[..r].iter()).copied().collect();
                        if let Some(x) = encode(&rot) {
                            set.insert(x);
                        }
                    }
                }
            }
            if set.is_empty() {
                continue;
            }
            lens.push(name.len());
            names.push(name);
            units.push(set);
        }
        Self { names, lens, units }
    }

    /// Per-position motif index + 1 (0 = no motif unit starts here).
    pub fn indicator(&self, seq: &[u8]) -> Vec<u8> {
        let mut ind = vec![0u8; seq.len()];
        let mut by_len: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
        for (i, &l) in self.lens.iter().enumerate() {
            by_len.entry(l).or_default().push(i);
        }
        for (&l, idxs) in &by_len {
            let m = if l == 32 { u64::MAX } else { (1u64 << (2 * l)) - 1 };
            let (mut x, mut valid) = (0u64, 0usize);
            for (i, &b) in seq.iter().enumerate() {
                let c = ENC[b as usize];
                if c > 3 {
                    valid = 0;
                    continue;
                }
                x = ((x << 2) | c as u64) & m;
                valid += 1;
                if valid >= l {
                    let p = i + 1 - l;
                    if ind[p] == 0 {
                        for &mi in idxs {
                            if self.units[mi].contains(&x) {
                                ind[p] = (mi + 1) as u8;
                                break;
                            }
                        }
                    }
                }
            }
        }
        ind
    }
}

const WIN: usize = 50;
const STEP: usize = 5;
const DENS: f64 = 0.30;

fn density(prefix: &[u32], a: usize, b: usize) -> f64 {
    (prefix[b] - prefix[a]) as f64 / (b - a).max(1) as f64
}

/// Repeat run forward from p0. Returns (run length, reached read end).
fn scan_fwd(prefix: &[u32], p0: usize) -> (usize, bool) {
    let n = prefix.len() - 1;
    if p0 >= n {
        return (0, true);
    }
    let mut j = p0;
    while j + WIN <= n {
        if density(prefix, j, j + WIN) < DENS {
            // repeat ended inside this window; credit the covered fraction
            let extra = (density(prefix, j, j + WIN) * WIN as f64) as usize;
            return (j - p0 + extra, false);
        }
        j += STEP;
    }
    let tail_ok = n - j < WIN / 2 || density(prefix, j, n) >= DENS;
    if tail_ok { (n - p0, true) } else { (j - p0, false) }
}

/// Repeat run backward ending at p1 (exclusive).
fn scan_bwd(prefix: &[u32], p1: usize) -> (usize, bool) {
    let n = prefix.len() - 1;
    let p1 = p1.min(n);
    if p1 == 0 {
        return (0, true);
    }
    let mut j = p1;
    while j >= WIN {
        if density(prefix, j - WIN, j) < DENS {
            let extra = (density(prefix, j - WIN, j) * WIN as f64) as usize;
            return (p1 - j + extra, false);
        }
        j -= STEP.min(j);
        if j < WIN {
            break;
        }
    }
    let tail_ok = j < WIN / 2 || density(prefix, 0, j) >= DENS;
    if tail_ok { (p1, true) } else { (p1 - j, false) }
}

struct SideAnchor {
    edge: i64,
    n: usize,
}

/// Inner-edge position (read coordinates, reference orientation) from consistent hits.
fn side_anchor(hits: &[(i64, i64)], edge_ref: i64, left: bool, min_hits: usize) -> Option<SideAnchor> {
    if hits.len() < min_hits {
        return None;
    }
    let mut diag: Vec<i64> = hits.iter().map(|(r, f)| f - r).collect();
    diag.sort_unstable();
    let med = diag[diag.len() / 2];
    let mut kept: Vec<(i64, i64)> = hits.iter().copied().filter(|(r, f)| ((f - r) - med).abs() <= 300).collect();
    if kept.len() < min_hits {
        return None;
    }
    if left {
        kept.sort_by_key(|&(_, f)| std::cmp::Reverse(f));
    } else {
        kept.sort_by_key(|&(_, f)| f);
    }
    let mut est: Vec<i64> = kept.iter().take(5).map(|(r, f)| r + (edge_ref - f)).collect();
    est.sort_unstable();
    Some(SideAnchor { edge: est[est.len() / 2], n: kept.len() })
}

fn composition(seq: &[u8], ind: &[u8], ms: &MotifSet, a: usize, b: usize) -> (String, f64, f64) {
    let b = b.min(seq.len());
    if a >= b {
        return ("-".into(), 0.0, 0.0);
    }
    let mut cnt = vec![0usize; ms.names.len()];
    let mut any = 0usize;
    for &v in &ind[a..b] {
        if v > 0 {
            cnt[(v - 1) as usize] += 1;
            any += 1;
        }
    }
    let len = (b - a) as f64;
    let (bi, bc) = cnt.iter().enumerate().max_by_key(|(_, c)| **c).map(|(i, c)| (i, *c)).unwrap_or((0, 0));
    if bc == 0 {
        return ("-".into(), 0.0, any as f64 / len);
    }
    (ms.names[bi].clone(), bc as f64 / len, any as f64 / len)
}

pub fn state_of(l: &Locus, bp: i64, bound: Bound) -> &'static str {
    if l.is_control() {
        return "CONTROL";
    }
    // compare whole repeat units (nearest integer), so a partial unit in the region does not cross a threshold
    let u = l.units(bp).round();
    if l.contraction {
        return match bound {
            Bound::Exact | Bound::Estimate if l.path_max_units.is_some_and(|p| u <= p) => "CONTRACTED",
            Bound::Exact | Bound::Estimate => "NORMAL",
            Bound::Lower => "UNRESOLVED",
            Bound::NA => "FLANK_ONLY",
        };
    }
    let (nm, pm) = (l.normal_max_units, l.path_min_units);
    match bound {
        Bound::Exact | Bound::Estimate => {
            if u <= nm {
                "NORMAL"
            } else if u >= pm {
                "EXPANDED"
            } else {
                "INTERMEDIATE"
            }
        }
        Bound::Lower => {
            if u >= pm {
                "EXPANDED_LB"
            } else if u > nm {
                "INTERMEDIATE_LB"
            } else {
                "UNRESOLVED"
            }
        }
        Bound::NA => "FLANK_ONLY",
    }
}

/// State of a placed read from its class, length, purity and dominant motif (also used when re-reporting).
pub fn finalize_state(call: &mut Call, l: &Locus) {
    call.state = match call.cls {
        Cls::Ambiguous => "AMBIGUOUS",
        Cls::FlankLeft | Cls::FlankRight => {
            if l.is_control() { "CONTROL" } else { "FLANK_ONLY" }
        }
        _ => state_of(l, call.allele_bp.unwrap_or(0), call.bound),
    };
    if matches!(call.state, "EXPANDED" | "EXPANDED_LB" | "INTERMEDIATE" | "INTERMEDIATE_LB") {
        if call.purity < MIN_PURITY {
            call.state = "NONREPEAT";
            if !call.flags.contains(&"LOW_PURITY_INSERT") {
                call.flags.push("LOW_PURITY_INSERT");
            }
        } else if !call.path_motif {
            call.state = "BENIGN_MOTIF";
        }
    }
}

/// Minimum repeat purity (fraction of positions starting a motif unit) for a long allele to
/// count as an expansion; long low-purity inserts are artefacts or non-repeat insertions.
pub const MIN_PURITY: f64 = 0.55;

/// Analyse one read against every locus it hits. `min_hits`: consistent flank k-mer hits
/// required on a side (real reads give tens to thousands; chance matches give 1–5).
pub fn analyze(seq: &[u8], panel: &Panel, loci: &[Locus], motif_sets: &[MotifSet], min_hits: usize) -> Vec<Call> {
    let k = panel.k as i64;
    let len = seq.len() as i64;
    let pad = panel.pad as i64;
    // (strand, side, rpos in reference orientation, refpos)
    let mut by_locus: FxHashMap<u32, Vec<(bool, u8, i64, i64)>> = FxHashMap::default();
    for_each_kmer(seq, panel.k, |pos, canon, fwd| {
        for e in panel.get(canon) {
            let plus = fwd == e.canon_is_fwd;
            let rp = if plus { pos as i64 } else { len - pos as i64 - k };
            by_locus.entry(e.locus).or_default().push((plus, e.side, rp, e.refpos as i64));
        }
    });
    let mut calls = Vec::new();
    let mut rc_cache: Option<Vec<u8>> = None;
    for (li, hits) in by_locus {
        let li = li as usize;
        let l = &loci[li];
        let n_plus = hits.iter().filter(|h| h.0).count();
        let n_minus = hits.len() - n_plus;
        let plus = n_plus >= n_minus;
        let mut flags = Vec::new();
        if n_plus.min(n_minus) >= 5 && (n_plus.min(n_minus) as f64) >= 0.3 * n_plus.max(n_minus) as f64 {
            flags.push("BOTH_STRANDS");
        }
        let hl: Vec<(i64, i64)> = hits.iter().filter(|h| h.0 == plus && h.1 == 0).map(|h| (h.2, h.3)).collect();
        let hr: Vec<(i64, i64)> = hits.iter().filter(|h| h.0 == plus && h.1 == 1).map(|h| (h.2, h.3)).collect();
        let e_l = l.start as i64 - pad;
        let e_r = l.end as i64 + pad;
        let a = side_anchor(&hl, e_l, true, min_hits);
        let b = side_anchor(&hr, e_r, false, min_hits);
        if a.is_none() && b.is_none() {
            continue; // chance k-mer matches only
        }
        let reach_l = a.as_ref().map(|x| x.edge + pad < len).unwrap_or(false);
        let reach_r = b.as_ref().map(|x| x.edge - pad > 0).unwrap_or(false);

        let oseq: &[u8] = if plus {
            seq
        } else {
            rc_cache.get_or_insert_with(|| revcomp(seq)).as_slice()
        };
        let ms = &motif_sets[li];

        let mut call = Call {
            locus: li,
            strand: if plus { '+' } else { '-' },
            cls: Cls::Ambiguous,
            allele_bp: None,
            bound: Bound::NA,
            state: "AMBIGUOUS",
            dom_motif: "-".into(),
            dom_frac: 0.0,
            purity: 0.0,
            path_motif: false,
            hits_l: a.as_ref().map(|x| x.n).unwrap_or(0),
            hits_r: b.as_ref().map(|x| x.n).unwrap_or(0),
            flags,
        };

        let need_ind = (a.is_some() && reach_l) || (b.is_some() && reach_r);
        let ind = if need_ind { ms.indicator(oseq) } else { Vec::new() };
        let prefix: Vec<u32> = if need_ind {
            let mut p = Vec::with_capacity(ind.len() + 1);
            p.push(0u32);
            let mut s = 0u32;
            for &v in &ind {
                s += (v > 0) as u32;
                p.push(s);
            }
            p
        } else {
            Vec::new()
        };

        let region: Option<(usize, usize)>;
        match (&a, &b) {
            (Some(x), Some(y)) if reach_l && reach_r => {
                let inner = y.edge - x.edge;
                if inner < -50 {
                    call.cls = Cls::Ambiguous;
                    call.flags.push("FLANK_ORDER");
                    region = None;
                } else {
                    call.cls = Cls::Spanning;
                    call.allele_bp = Some((inner - 2 * pad).max(0));
                    call.bound = Bound::Exact;
                    let s = (x.edge + pad).clamp(0, len) as usize;
                    let e = (y.edge - pad).clamp(0, len) as usize;
                    region = Some((s, e.max(s)));
                }
            }
            (Some(x), None) if reach_l => {
                call.cls = Cls::LeftAnchor;
                let p0 = (x.edge + pad).clamp(0, len) as usize;
                let (run, end) = scan_fwd(&prefix, p0);
                call.allele_bp = Some(run as i64);
                if end {
                    call.bound = Bound::Lower;
                } else {
                    call.bound = Bound::Estimate;
                    call.flags.push("NO_OPPOSITE_FLANK");
                }
                region = Some((p0, p0 + run));
            }
            (None, Some(y)) if reach_r => {
                call.cls = Cls::RightAnchor;
                let p1 = (y.edge - pad).clamp(0, len) as usize;
                let (run, end) = scan_bwd(&prefix, p1);
                call.allele_bp = Some(run as i64);
                if end {
                    call.bound = Bound::Lower;
                } else {
                    call.bound = Bound::Estimate;
                    call.flags.push("NO_OPPOSITE_FLANK");
                }
                region = Some((p1.saturating_sub(run), p1));
            }
            (Some(_), None) => {
                call.cls = Cls::FlankLeft;
                region = None;
            }
            (None, Some(_)) => {
                call.cls = Cls::FlankRight;
                region = None;
            }
            (Some(_), Some(_)) => {
                call.cls = Cls::Ambiguous;
                call.flags.push("FLANK_REACH");
                region = None;
            }
            (None, None) => {
                call.cls = Cls::Ambiguous;
                call.flags.push("FEW_HITS");
                region = None;
            }
        }
        if let Some((s, e)) = region {
            let (m, f, p) = composition(oseq, &ind, ms, s, e);
            call.path_motif = m != "-" && l.motif_is_pathogenic(&m);
            call.dom_motif = m;
            call.dom_frac = f;
            call.purity = p;
        }
        finalize_state(&mut call, l);
        calls.push(call);
    }
    calls.sort_by_key(|c| c.locus);
    calls
}
