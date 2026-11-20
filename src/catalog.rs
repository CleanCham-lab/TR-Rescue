// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Tandem-repeat disease locus catalog.
//!
//! Built-in: an STRchive snapshot (hg38, 0-based half-open, thresholds in repeat units),
//! embedded at compile time from `data/strchive_hg38.tsv`. Extra loci: TRGT-style BED.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::kmer::revcomp;

const STRCHIVE_TSV: &str = include_str!("../data/strchive_hg38.tsv");

#[derive(Clone, Debug, Default)]
pub struct Locus {
    pub id: String,
    pub gene: String,
    pub disease_id: String,
    pub disease: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// all motifs used for repeat detection (reference, pathogenic, benign, unknown)
    pub motifs: Vec<String>,
    pub ref_motifs: Vec<String>,
    /// motifs whose expansion is pathogenic (empty = any)
    pub path_motifs: Vec<String>,
    pub unit: usize,
    pub normal_max_units: f64,
    pub path_min_units: f64,
    pub path_max_units: Option<f64>,
    pub thresholds_known: bool,
    pub inheritance: String,
    pub evidence: String,
    pub location: String,
    pub source: String,
    pub stripy: String,
    pub gnomad: String,
    pub omim: String,
    pub genereviews: String,
    pub control_of: Option<usize>,
    pub control_offset: i64,
    /// bp of non-repeat sequence inside the catalog region (repeat count = (allele bp - this) / unit)
    pub count_offset_bp: i64,
    /// why the region or count was adjusted (shown in the report)
    pub count_note: crate::i18n::Tr,
    /// non-motif bp at the left/right edge of the reference region (QC, filled at run time)
    pub ref_edge_bp: Option<(i64, i64)>,
    pub benign_min_units: Option<f64>,
    /// disease is caused by a contraction (pathogenic_max < benign_min)
    pub contraction: bool,
    /// the reference allele itself is above the normal maximum (region and thresholds do not match)
    pub ref_over_normal: bool,
}

impl Locus {
    pub fn unit_len(&self) -> usize {
        self.unit.max(1)
    }
    pub fn normal_max_bp(&self) -> i64 {
        (self.normal_max_units * self.unit_len() as f64).round() as i64 + self.count_offset_bp
    }
    pub fn path_min_bp(&self) -> i64 {
        (self.path_min_units * self.unit_len() as f64).round() as i64 + self.count_offset_bp
    }
    /// Repeat count (clinical convention) for an allele length in bp.
    pub fn units(&self, bp: i64) -> f64 {
        (bp - self.count_offset_bp) as f64 / self.unit_len() as f64
    }
    pub fn is_control(&self) -> bool {
        self.control_of.is_some()
    }
    pub fn centre(&self) -> u64 {
        (self.start + self.end) / 2
    }
    pub fn region_str(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start + 1, self.end)
    }
    pub fn is_sex_chrom(&self) -> bool {
        let c = self.chrom.trim_start_matches("chr");
        c == "X" || c == "Y"
    }
    /// Is the given (dominant) motif pathogenic at this locus? Rotation/strand insensitive.
    pub fn motif_is_pathogenic(&self, m: &str) -> bool {
        self.path_motifs.is_empty() || self.path_motifs.iter().any(|p| motif_eq(p, m))
    }
}

const IUPAC: [(u8, &[u8]); 11] = [
    (b'N', b"ACGT"),
    (b'R', b"AG"),
    (b'Y', b"CT"),
    (b'S', b"CG"),
    (b'W', b"AT"),
    (b'K', b"GT"),
    (b'M', b"AC"),
    (b'B', b"CGT"),
    (b'D', b"AGT"),
    (b'H', b"ACT"),
    (b'V', b"ACG"),
];

fn bases_of(b: u8) -> &'static [u8] {
    match b.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        x => IUPAC.iter().find(|(k, _)| *k == x).map(|(_, v)| *v).unwrap_or(b""),
    }
}

/// Concrete motifs for a possibly degenerate one (GCN -> GCA, GCC, GCG, GCT). Capped at 64.
pub fn expand_motif(m: &str) -> Vec<String> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    for &b in m.as_bytes() {
        let opts = bases_of(b);
        if opts.is_empty() {
            return vec![];
        }
        out = out.iter().flat_map(|p| opts.iter().map(move |&o| [p.as_slice(), &[o]].concat())).collect();
        if out.len() > 64 {
            return vec![];
        }
    }
    out.into_iter().map(|v| String::from_utf8(v).unwrap()).collect()
}

/// Same motif up to rotation / strand, with IUPAC codes in either one matching their bases.
pub fn motif_eq(a: &str, b: &str) -> bool {
    if a.len() != b.len() || a.is_empty() {
        return false;
    }
    let a = a.to_ascii_uppercase().into_bytes();
    let b = b.to_ascii_uppercase().into_bytes();
    let rc = |s: &[u8]| -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&x| match x {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'R' => b'Y',
                b'Y' => b'R',
                b'K' => b'M',
                b'M' => b'K',
                b'B' => b'V',
                b'V' => b'B',
                b'D' => b'H',
                b'H' => b'D',
                o => o,
            })
            .collect()
    };
    let compat = |x: u8, y: u8| bases_of(x).iter().any(|c| bases_of(y).contains(c));
    for src in [b.clone(), rc(&b)] {
        for r in 0..src.len() {
            if (0..src.len()).all(|i| compat(a[i], src[(i + r) % src.len()])) {
                return true;
            }
        }
    }
    false
}

/// Loci where STRchive's region and the clinical repeat count differ (checked against the ExpansionHunter
/// catalog and the clinical convention): (id, new start, new end, count offset bp, extra motifs, note).
const COUNT_FIXES: [(&str, Option<u64>, Option<u64>, i64, &str, &str, &str); 3] = [
    (
        "SBMA_AR",
        None,
        None,
        34,
        "",
        "STRchive 구간 103 bp 중 CAG 반복은 앞 69 bp(ExpansionHunter 구간)이고 뒤 34 bp는 반복이 아님 → 반복 수 = (길이 − 34) ÷ 3. bp는 TRGT와 같음",
        "Of the 103 bp STRchive region, the CAG repeat is the first 69 bp (ExpansionHunter region) and the last 34 bp are not repeat → count = (length − 34) / 3. bp values match TRGT",
    ),
    (
        "FRAXE_AFF2",
        None,
        None,
        89,
        "",
        "STRchive 구간 149 bp 중 GCC 반복은 60 bp(ExpansionHunter 구간), 앞 27 bp·뒤 62 bp는 반복이 아님 → 반복 수 = (길이 − 89) ÷ 3. bp는 TRGT와 같음",
        "Of the 149 bp STRchive region, the GCC repeat is 60 bp (ExpansionHunter region); 27 bp before and 62 bp after are not repeat → count = (length − 89) / 3. bp values match TRGT",
    ),
    (
        "SCA8_ATXN8OS",
        Some(70_139_353),
        None,
        0,
        "CTA",
        "임상 관례(GeneReviews)대로 CTA+CTG 합산: 구간을 앞쪽 CTA 반복까지 넓힘(TRGT STRchive 카탈로그와 같은 구간)",
        "Combined CTA+CTG count as in clinical practice (GeneReviews): region extended over the preceding CTA repeat (same region as TRGT's STRchive catalog)",
    ),
];

/// Apply COUNT_FIXES to built-in loci.
fn apply_fixes(loci: &mut [Locus]) {
    for l in loci.iter_mut() {
        if let Some((_, s, e, off, extra, note, note_en)) = COUNT_FIXES.iter().find(|f| f.0 == l.id) {
            if let Some(s) = s {
                l.start = *s;
            }
            if let Some(e) = e {
                l.end = *e;
            }
            l.count_offset_bp = *off;
            for m in extra.split(',').filter(|m| !m.is_empty()) {
                if !l.motifs.iter().any(|x| x == m) {
                    l.motifs.push(m.to_string());
                }
            }
            l.count_note = crate::i18n::Tr::new(*note, *note_en);
        }
    }
}

/// Make thresholds self-consistent: contraction loci, pathogenic minimum at or below the normal maximum
/// (e.g. FAME loci where the pathogenic count refers to the inserted motif only), reference above normal.
pub fn finalize_thresholds(l: &mut Locus) {
    let ref_u = l.units((l.end - l.start) as i64).round();
    if let (Some(pmax), Some(bmin)) = (l.path_max_units, l.benign_min_units) {
        if pmax < bmin {
            l.contraction = true;
            return;
        }
    }
    if l.path_min_units <= l.normal_max_units {
        if l.path_min_units <= ref_u {
            // pathogenic count is of an inserted motif: require that many units on top of the reference allele
            l.normal_max_units = l.normal_max_units.max(ref_u);
            l.path_min_units += l.normal_max_units;
        } else {
            l.path_min_units = l.normal_max_units + 1.0;
        }
        l.thresholds_known = false;
    }
    l.ref_over_normal = ref_u > l.normal_max_units;
}

/// Non-motif bp at the left and right edge of a reference repeat region (best phase).
pub fn edge_nonrepeat(seq: &[u8], motifs: &[String], unit: usize) -> Option<(i64, i64)> {
    if unit == 0 || unit > 6 || seq.len() < 2 * unit {
        return None;
    }
    let mut set: Vec<String> = Vec::new();
    for m in motifs {
        if m.len() == unit {
            set.extend(expand_motif(m));
        }
    }
    if set.is_empty() {
        return None;
    }
    let is_rep = |t: &[u8]| {
        let t = String::from_utf8_lossy(t).to_ascii_uppercase();
        set.iter().any(|m| motif_eq(m, &t))
    };
    // main repeat window per phase = maximum-scoring run of tiles (+1 motif tile, -2 other tile)
    let mut best: Option<(i64, i64, i64)> = None;
    for p in 0..unit {
        let tiles: Vec<bool> = (p..seq.len().saturating_sub(unit - 1)).step_by(unit).map(|i| is_rep(&seq[i..i + unit])).collect();
        let (mut cur, mut cur_start) = (0i64, 0usize);
        let mut top: Option<(i64, usize, usize)> = None;
        for (i, &t) in tiles.iter().enumerate() {
            if cur <= 0 {
                cur = 0;
                cur_start = i;
            }
            cur += if t { 1 } else { -2 };
            if cur > 0 && top.map(|x| cur > x.0).unwrap_or(true) {
                top = Some((cur, cur_start, i));
            }
        }
        if let Some((sc, i, j)) = top {
            if sc < 2 {
                continue;
            }
            let left = (p + i * unit) as i64;
            let right = seq.len() as i64 - (p + (j + 1) * unit) as i64;
            if best.map(|b| sc > b.0).unwrap_or(true) {
                best = Some((sc, left, right));
            }
        }
    }
    best.map(|b| (b.1, b.2))
}

/// Rotation- and strand-independent motif key.
pub fn canonical(m: &str) -> String {
    let m = m.to_ascii_uppercase().into_bytes();
    if m.is_empty() {
        return String::new();
    }
    let mut best: Option<Vec<u8>> = None;
    for src in [m.clone(), revcomp(&m)] {
        for r in 0..src.len() {
            let rot: Vec<u8> = src[r..].iter().chain(src[..r].iter()).copied().collect();
            if best.as_ref().map(|b| rot < *b).unwrap_or(true) {
                best = Some(rot);
            }
        }
    }
    String::from_utf8(best.unwrap()).unwrap_or_default()
}

fn list(s: &str) -> Vec<String> {
    s.split(',').map(|x| x.trim().to_ascii_uppercase()).filter(|x| !x.is_empty()).collect()
}

fn num(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.is_empty() { None } else { s.parse().ok() }
}

/// Threshold fill-in when the source leaves a bound empty.
fn thresholds(ref_units: f64, benign_max: Option<f64>, int_min: Option<f64>, path_min: Option<f64>) -> (f64, f64, bool) {
    match (benign_max, path_min) {
        (Some(b), Some(p)) => (b, p, true),
        (None, Some(p)) => (int_min.map(|i| i - 1.0).unwrap_or((p * 0.5).max(ref_units)), p, false),
        (Some(b), None) => (b, (b * 2.0).max(b + 20.0), false),
        (None, None) => {
            let n = (ref_units * 1.5).max(ref_units + 10.0);
            (n, (n * 2.0).max(n + 50.0), false)
        }
    }
}

/// STRchive snapshot embedded in the binary.
pub fn builtin() -> Vec<Locus> {
    let mut out = Vec::new();
    let mut header: Option<Vec<&str>> = None;
    for line in STRCHIVE_TSV.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let Some(h) = &header else {
            header = Some(f);
            continue;
        };
        let get = |k: &str| h.iter().position(|c| *c == k).and_then(|i| f.get(i)).copied().unwrap_or("");
        let (Ok(start), Ok(end)) = (get("start_hg38").parse::<u64>(), get("stop_hg38").parse::<u64>()) else { continue };
        let ref_m = list(get("reference_motif_reference_orientation"));
        let path_m = list(get("pathogenic_motif_reference_orientation"));
        let mut motifs: Vec<String> = Vec::new();
        for m in ref_m
            .iter()
            .chain(path_m.iter())
            .chain(list(get("benign_motif_reference_orientation")).iter())
            .chain(list(get("unknown_motif_reference_orientation")).iter())
        {
            if !m.contains('N') && !motifs.iter().any(|x| canonical(x) == canonical(m)) {
                motifs.push(m.clone());
            }
        }
        let unit = get("motif_len").parse::<usize>().ok().or_else(|| ref_m.first().map(|m| m.len())).unwrap_or(3);
        let ref_units = (end - start) as f64 / unit.max(1) as f64;
        let (nm, pm, known) = thresholds(ref_units, num(get("benign_max")), num(get("intermediate_min")), num(get("pathogenic_min")));
        // a pathogenic motif list identical to the reference motif adds nothing
        let path_motifs = if !path_m.is_empty() && path_m.iter().all(|p| ref_m.iter().any(|r| canonical(r) == canonical(p))) {
            vec![]
        } else {
            path_m
        };
        out.push(Locus {
            id: get("id").to_string(),
            gene: get("gene").to_string(),
            disease_id: get("disease_id").to_string(),
            disease: get("disease").to_string(),
            chrom: get("chrom").to_string(),
            start,
            end,
            motifs,
            ref_motifs: ref_m,
            path_motifs,
            unit,
            normal_max_units: nm,
            path_min_units: pm,
            path_max_units: num(get("pathogenic_max")),
            benign_min_units: num(get("benign_min")),
            thresholds_known: known,
            inheritance: get("inheritance").to_string(),
            evidence: get("evidence").to_string(),
            location: get("location_in_gene").to_string(),
            source: "STRchive".into(),
            stripy: get("stripy").to_string(),
            gnomad: get("gnomad").to_string(),
            omim: get("omim").to_string(),
            genereviews: get("genereviews").to_string(),
            control_of: None,
            control_offset: 0,
            ..Default::default()
        });
    }
    apply_fixes(&mut out);
    out.iter_mut().for_each(finalize_thresholds);
    out
}

fn info_get<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| {
        let mut it = kv.splitn(2, '=');
        let k = it.next()?.trim();
        let v = it.next()?.trim();
        if k.eq_ignore_ascii_case(key) { Some(v) } else { None }
    })
}

/// TRGT-style catalog BED: `chrom start end ID=..;MOTIFS=A,B;...`.
/// Optional keys: NORMAL_MAX, PATH_MIN (repeat units), PATH_MOTIFS, INHERITANCE, GENE, DISEASE.
pub fn parse_bed(path: &Path) -> Result<Vec<Locus>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading catalog {}", path.display()))?;
    let mut out = Vec::new();
    for (ln, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            bail!("{}:{}: expected >=4 tab-separated columns", path.display(), ln + 1);
        }
        let start: u64 = f[1].parse().with_context(|| format!("line {}: start", ln + 1))?;
        let end: u64 = f[2].parse().with_context(|| format!("line {}: end", ln + 1))?;
        let info = f[3];
        let id = info_get(info, "ID").unwrap_or(f[3]).to_string();
        let motifs = info_get(info, "MOTIFS").map(list).unwrap_or_default();
        if motifs.is_empty() {
            bail!("{}:{}: no MOTIFS for {}", path.display(), ln + 1, id);
        }
        let unit = motifs[0].len().max(1);
        let ref_units = (end.saturating_sub(start)) as f64 / unit as f64;
        let (nm, pm, known) = thresholds(
            ref_units,
            info_get(info, "NORMAL_MAX").and_then(|v| v.parse().ok()),
            None,
            info_get(info, "PATH_MIN").and_then(|v| v.parse().ok()),
        );
        out.push(Locus {
            gene: info_get(info, "GENE").unwrap_or(&id).to_string(),
            disease: info_get(info, "DISEASE").unwrap_or("").to_string(),
            id,
            chrom: f[0].to_string(),
            start,
            end,
            ref_motifs: motifs.clone(),
            motifs,
            path_motifs: info_get(info, "PATH_MOTIFS").map(list).unwrap_or_default(),
            unit,
            normal_max_units: nm,
            path_min_units: pm,
            thresholds_known: known,
            inheritance: info_get(info, "INHERITANCE").unwrap_or("").to_string(),
            source: format!("BED:{}", path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default()),
            ..Default::default()
        });
    }
    Ok(out)
}

/// Find a locus by STRchive id, disease id or (unique) gene symbol, case-insensitive.
pub fn find(cat: &[Locus], key: &str) -> Result<Locus> {
    if let Some(l) = cat.iter().find(|l| l.id.eq_ignore_ascii_case(key)) {
        return Ok(l.clone());
    }
    let by_gene: Vec<&Locus> = cat.iter().filter(|l| l.gene.eq_ignore_ascii_case(key)).collect();
    if by_gene.len() == 1 {
        return Ok(by_gene[0].clone());
    }
    if by_gene.len() > 1 {
        bail!(
            "{key} matches several loci ({}); use the full ID",
            by_gene.iter().map(|l| l.id.as_str()).collect::<Vec<_>>().join(", ")
        );
    }
    let by_dis: Vec<&Locus> = cat.iter().filter(|l| l.disease_id.eq_ignore_ascii_case(key)).collect();
    if by_dis.len() == 1 {
        return Ok(by_dis[0].clone());
    }
    bail!("locus {key} not in catalog (see `tr-rescue catalog`)")
}

/// Neighbourhood control pseudo-loci: zero-length points outside the read-length canyon.
pub const CONTROL_OFFSETS: [i64; 6] = [-160_000, -130_000, -100_000, 100_000, 130_000, 160_000];

pub fn add_controls(targets: &[Locus]) -> Vec<Locus> {
    let mut all: Vec<Locus> = targets.to_vec();
    for (ti, t) in targets.iter().enumerate() {
        for off in CONTROL_OFFSETS {
            let c = t.centre() as i64 + off;
            if c < 10_000 {
                continue;
            }
            let c = c as u64;
            all.push(Locus {
                id: format!("{}_ctl{:+}k", t.id, off / 1000),
                gene: t.gene.clone(),
                chrom: t.chrom.clone(),
                start: c,
                end: c,
                motifs: t.motifs.clone(),
                ref_motifs: t.ref_motifs.clone(),
                unit: t.unit,
                normal_max_units: t.normal_max_units,
                path_min_units: t.path_min_units,
                inheritance: "control".into(),
                source: "control".into(),
                control_of: Some(ti),
                control_offset: off,
                ..Default::default()
            });
        }
    }
    all
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degenerate_motifs() {
        assert_eq!(expand_motif("GCN").len(), 4);
        assert!(expand_motif("GCN").contains(&"GCG".to_string()));
        assert!(motif_eq("GCN", "GCG") && motif_eq("NGC", "CGC") && motif_eq("GCN", "CGC"));
        assert!(!motif_eq("GCN", "CAT"));
        // a CAG repeat read in another frame is (GCA)n, which is a GCN repeat
        assert!(motif_eq("GCN", "CAG"));
        assert!(motif_eq("CAG", "CTG") && motif_eq("AAGGG", "CCCTT") && !motif_eq("AAGGG", "AAAAG"));
    }

    #[test]
    fn count_fixes() {
        let cat = builtin();
        let ar = cat.iter().find(|l| l.id == "SBMA_AR").unwrap();
        assert_eq!(ar.count_offset_bp, 34);
        assert!((ar.units(112) - 26.0).abs() < 0.01 && (ar.units(166) - 44.0).abs() < 0.01);
        assert_eq!(ar.normal_max_bp(), 34 * 3 + 34);
        let a8 = cat.iter().find(|l| l.id == "SCA8_ATXN8OS").unwrap();
        assert_eq!(a8.start, 70_139_353);
        assert!(a8.motifs.iter().any(|m| m == "CTA"));
        let fe = cat.iter().find(|l| l.id == "FRAXE_AFF2").unwrap();
        assert!((fe.units(149) - 20.0).abs() < 0.01);
        // polyalanine loci have usable motifs
        for l in cat.iter().filter(|l| l.motifs.iter().any(|m| m.contains('N'))) {
            assert!(l.motifs.iter().all(|m| !expand_motif(m).is_empty()), "{}", l.id);
            assert!(l.motif_is_pathogenic("GCG") || l.motif_is_pathogenic("CGC") || l.path_motifs.is_empty(), "{}", l.id);
        }
    }

    #[test]
    fn edges() {
        // (CAG)23 + 34 bp of non-repeat: right edge flagged
        let mut s = "CAG".repeat(23);
        s.push_str("CAAGAGACTAGCCCCAGGCAGCAGCAGCAGCAGC");
        let (l, r) = edge_nonrepeat(s.as_bytes(), &["CAG".into()], 3).unwrap();
        assert_eq!(l, 0);
        assert!(r >= 18, "right edge {r}");
        // interrupted CGG (AGG) is not flagged at the edges
        let f = format!("{}AGG{}AGG{}", "CGG".repeat(10), "CGG".repeat(9), "CGG".repeat(9));
        let (l, r) = edge_nonrepeat(f.as_bytes(), &["CGG".into()], 3).unwrap();
        assert_eq!((l, r), (0, 0));
    }
}
