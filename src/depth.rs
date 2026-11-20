// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Read-length-wide depth canyon: expected shape under whole-read loss and fits.
//!
//! Total depth: if a fraction f of the molecules that contain the repeat is lost (whole read,
//! flanks included), the relative depth at distance x from the repeat is 1 - f * g(x), where
//! g(x) = sum(max(0, L - x)) / sum(L) over the sample's read lengths L. f = 0.5 is a complete
//! one-allele loss, f = 1 both alleles.
//!
//! Haplotypes: phase blocks rarely cover the whole ±window, so each haplotype's baseline is
//! taken inside the phase block that covers the repeat, from bins where the expected canyon
//! is shallow (g <= 0.1). f_hp = 1 means that haplotype lost every molecule containing the
//! repeat; one ≈ 1 and the other ≈ 0 is the one-allele pattern, both high the two-allele one.

use crate::hifi::Depth;

#[derive(Default)]
pub struct Canyon {
    /// bin centre offsets from the repeat (bp; negative = left)
    pub offsets: Vec<i64>,
    pub g: Vec<f64>,
    pub rel_total: Vec<f64>,
    /// HP_k / its in-block baseline (NaN outside the phase block)
    pub rel_hp1: Vec<f64>,
    pub rel_hp2: Vec<f64>,
    pub base_total: f64,
    /// per-haplotype baseline inside the phase block (NaN when not estimable)
    pub base_hp1: f64,
    pub base_hp2: f64,
    pub centre_total: f64,
    pub f_total: Option<f64>,
    /// per-haplotype loss fraction
    pub f_hp1: Option<f64>,
    pub f_hp2: Option<f64>,
    pub se_total: Option<f64>,
    pub se_hp1: Option<f64>,
    pub se_hp2: Option<f64>,
    pub phase_set: Option<i64>,
    pub read_n50: u64,
}

fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}

pub fn n50(lens: &[u64]) -> u64 {
    let mut v = lens.to_vec();
    v.sort_unstable_by(|a, b| b.cmp(a));
    let tot: u64 = v.iter().sum();
    let mut acc = 0;
    for x in v {
        acc += x;
        if acc * 2 >= tot {
            return x;
        }
    }
    0
}

pub fn g_of(lens: &[u64], x: f64) -> f64 {
    let tot: f64 = lens.iter().map(|&l| l as f64).sum();
    if tot <= 0.0 {
        return 0.0;
    }
    lens.iter().map(|&l| (l as f64 - x).max(0.0)).sum::<f64>() / tot
}

type Span = (u64, u64, u8, i64);

/// Bin spans into total / HP1 / HP2 mean depth; HP only for reads in phase set `ps`.
fn bin_spans<'a, I: Iterator<Item = &'a Span>>(spans: I, ps: Option<i64>, win_start: u64, bin: u64, nb: usize) -> [Vec<f64>; 3] {
    let mut t = [vec![0.0; nb], vec![0.0; nb], vec![0.0; nb]];
    let win_end = win_start + bin * nb as u64;
    for &(s, e, hp, sps) in spans {
        let hp = if ps.is_some() && Some(sps) == ps { hp } else { 0 };
        let a = s.max(win_start);
        let b = e.min(win_end);
        let mut p = a;
        while p < b {
            let bi = ((p - win_start) / bin) as usize;
            let q = b.min(win_start + (bi as u64 + 1) * bin);
            let w = (q - p) as f64 / bin as f64;
            t[0][bi] += w;
            if hp == 1 {
                t[1][bi] += w;
            } else if hp == 2 {
                t[2][bi] += w;
            }
            p = q;
        }
    }
    t
}

/// Total-depth fit against the distal baseline.
fn fit_total(v: &[f64], g: &[f64], offsets: &[i64], distal: i64) -> Option<f64> {
    let mut b: Vec<f64> = (0..v.len()).filter(|&i| offsets[i].abs() >= distal).map(|i| v[i]).collect();
    let base = median(&mut b);
    if !(base.is_finite() && base > 0.5) {
        return None;
    }
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..v.len() {
        if offsets[i].abs() >= distal {
            continue;
        }
        let rel = v[i] / base;
        num += g[i] * (1.0 - rel);
        den += g[i] * g[i];
    }
    (den > 1e-6).then(|| num / den)
}

/// Haplotype baseline: median depth over in-block bins with a shallow expected canyon.
fn hap_base(h: &[f64], g: &[f64], in_block: &[bool]) -> Option<f64> {
    let mut v: Vec<f64> = (0..h.len()).filter(|&i| in_block[i] && g[i] <= 0.1).map(|i| h[i]).collect();
    if v.len() < 8 {
        return None;
    }
    let m = median(&mut v);
    // per-haplotype depth below ~4x is too noisy to say anything
    (m.is_finite() && m >= 4.0).then_some(m)
}

/// Per-haplotype canyon fit inside the phase block.
fn fit_hap(h: &[f64], g: &[f64], in_block: &[bool]) -> Option<(f64, f64)> {
    let base = hap_base(h, g, in_block)?;
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..h.len() {
        if !in_block[i] || g[i] < 0.02 {
            continue;
        }
        num += g[i] * (1.0 - h[i] / base);
        den += g[i] * g[i];
    }
    // the block must reach into the canyon (sum of g^2 over >= ~3 kb of strong signal)
    (den > 2.0).then(|| (num / den, base))
}

fn sd(v: &[f64]) -> Option<f64> {
    if v.len() < 20 {
        return None;
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    Some((v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() - 1) as f64).sqrt())
}

/// Phase set carried by most phased reads overlapping the repeat ±5 kb.
fn dominant_ps(spans: &[Span], rep_start: u64, rep_end: u64) -> Option<i64> {
    let (a, b) = (rep_start.saturating_sub(5000), rep_end + 5000);
    let mut counts: std::collections::BTreeMap<i64, usize> = Default::default();
    for &(s, e, hp, ps) in spans {
        if hp != 0 && ps >= 0 && s < b && e > a {
            *counts.entry(ps).or_insert(0) += 1;
        }
    }
    counts.into_iter().max_by_key(|(_, c)| *c).map(|(p, _)| p)
}

/// `rep_start`/`rep_end`: repeat interval (0-based, half-open). `distal`: |offset| beyond
/// which bins form the total-depth baseline and are ignored by the fits.
pub fn canyon(d: &Depth, rep_start: u64, rep_end: u64, distal: i64, centre_half: i64) -> Canyon {
    let nb = d.total.len();
    let mut offsets = Vec::with_capacity(nb);
    for i in 0..nb {
        let c = d.win_start as i64 + (i as i64) * d.bin as i64 + d.bin as i64 / 2;
        let off = if c < rep_start as i64 {
            c - rep_start as i64
        } else if c >= rep_end as i64 {
            c - rep_end as i64
        } else {
            0
        };
        offsets.push(off);
    }
    let g: Vec<f64> = offsets.iter().map(|&o| g_of(&d.read_lens, o.unsigned_abs() as f64)).collect();
    let ps = dominant_ps(&d.spans, rep_start, rep_end);
    let [tot, h1, h2] = bin_spans(d.spans.iter(), ps, d.win_start, d.bin, nb);
    // phase block extent (reads of the dominant phase set)
    let (mut lo, mut hi) = (u64::MAX, 0u64);
    for &(s0, e0, hp, sps) in &d.spans {
        if hp != 0 && Some(sps) == ps {
            lo = lo.min(s0);
            hi = hi.max(e0);
        }
    }
    let in_block: Vec<bool> = (0..nb)
        .map(|i| {
            let c = d.win_start + i as u64 * d.bin + d.bin / 2;
            ps.is_some() && c >= lo && c < hi
        })
        .collect();

    let mut b: Vec<f64> = (0..nb).filter(|&i| offsets[i].abs() >= distal).map(|i| tot[i]).collect();
    let base_total = median(&mut b);
    let rel_total: Vec<f64> = tot.iter().map(|&x| if base_total > 0.0 { x / base_total } else { f64::NAN }).collect();
    let sel: Vec<f64> = (0..nb).filter(|&i| offsets[i].abs() <= centre_half).map(|i| tot[i]).collect();
    let centre_total = if sel.is_empty() { f64::NAN } else { sel.iter().sum::<f64>() / sel.len() as f64 };
    let f_total = fit_total(&tot, &g, &offsets, distal);
    let r1 = fit_hap(&h1, &g, &in_block);
    let r2 = fit_hap(&h2, &g, &in_block);
    let rel = |h: &[f64], r: Option<(f64, f64)>| -> Vec<f64> {
        (0..nb)
            .map(|i| match r {
                Some((_, b)) if in_block[i] => h[i] / b,
                _ => f64::NAN,
            })
            .collect()
    };

    // bootstrap SE: resample reads with replacement
    let (mut bt, mut b1, mut b2) = (Vec::new(), Vec::new(), Vec::new());
    let n = d.spans.len();
    if n >= 10 {
        let mut state: u64 = 0x2545_F491_4F6C_DD1D ^ n as u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..200 {
            let idx: Vec<usize> = (0..n).map(|_| (next() % n as u64) as usize).collect();
            let [t, a1, a2] = bin_spans(idx.iter().map(|&i| &d.spans[i]), ps, d.win_start, d.bin, nb);
            if let Some(f) = fit_total(&t, &g, &offsets, distal) {
                bt.push(f);
            }
            if let Some((f, _)) = fit_hap(&a1, &g, &in_block) {
                b1.push(f);
            }
            if let Some((f, _)) = fit_hap(&a2, &g, &in_block) {
                b2.push(f);
            }
        }
    }

    Canyon {
        rel_hp1: rel(&h1, r1),
        rel_hp2: rel(&h2, r2),
        base_hp1: r1.map(|x| x.1).unwrap_or(f64::NAN),
        base_hp2: r2.map(|x| x.1).unwrap_or(f64::NAN),
        se_total: if f_total.is_some() { sd(&bt) } else { None },
        se_hp1: if r1.is_some() { sd(&b1) } else { None },
        se_hp2: if r2.is_some() { sd(&b2) } else { None },
        f_total,
        f_hp1: r1.map(|x| x.0),
        f_hp2: r2.map(|x| x.0),
        base_total,
        centre_total,
        rel_total,
        phase_set: ps,
        read_n50: n50(&d.read_lens),
        offsets,
        g,
    }
}
