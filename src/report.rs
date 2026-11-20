// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! TSV outputs and the interactive HTML report (data as JSON + embedded template).

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

use crate::alleles::{self, AlleleSummary};
use crate::anchor::Call;
use crate::catalog::Locus;
use crate::i18n::Tr;
use crate::depth::Canyon;
use crate::scan::ScanStats;

const TEMPLATE: &str = include_str!("report_template.html");
const ICON: &[u8] = include_bytes!("../assets/icon.webp");
const BANNER: &[u8] = include_bytes!("../assets/banner.webp");
const EGG: &[u8] = include_bytes!("../assets/egg.webp");
const FAVICON: &[u8] = include_bytes!("../assets/favicon.png");

/// Standard base64 (for embedding images as data URIs).
fn b64(data: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for c in data.chunks(3) {
        let n = (c[0] as u32) << 16 | (*c.get(1).unwrap_or(&0) as u32) << 8 | *c.get(2).unwrap_or(&0) as u32;
        out.push(T[(n >> 18) as usize & 63] as char);
        out.push(T[(n >> 12) as usize & 63] as char);
        out.push(if c.len() > 1 { T[(n >> 6) as usize & 63] as char } else { '=' });
        out.push(if c.len() > 2 { T[n as usize & 63] as char } else { '=' });
    }
    out
}

pub struct ReadRow {
    pub sample: usize,
    pub source: String,
    pub bam: String,
    pub name: String,
    pub read_len: usize,
    pub ff: Option<i64>,
    pub rq: Option<f32>,
    pub np: Option<i64>,
    pub hp: Option<i64>,
    pub call: Call,
}

pub struct SampleInfo {
    pub id: String,
    pub group: String,
    /// positive / control / test
    pub role: &'static str,
    /// from the manifest ("M"/"F")
    pub sex_declared: Option<&'static str>,
    /// from chrX / autosome baseline depth ("M"/"F")
    pub sex_inferred: Option<&'static str>,
    /// median chrX / autosome baseline ratio over chrX loci
    pub x_ratio: Option<f64>,
    /// sample-level warnings (sex mismatch, missing locus for a positive, low depth)
    pub qc: Vec<Tr>,
    pub expected: Vec<usize>,
    pub expected_raw: String,
    pub note: String,
    pub hifi_bam: String,
    pub fail_bams: Vec<String>,
    pub hifi_raw_bams: Vec<String>,
    pub locate_notes: Vec<String>,
    /// unassigned BAMs found but not used
    pub unassigned: Vec<String>,
    /// unassigned BAMs scanned (single-plex cells, or listed in the manifest)
    pub unassigned_used: Vec<String>,
}

pub struct BamStat {
    pub sample: usize,
    pub kind: String,
    pub path: String,
    pub stats: ScanStats,
}

pub struct LocusSummary {
    pub sample: usize,
    pub locus: usize,
    pub canyon: Option<Canyon>,
    pub phased_frac: f64,
    pub error: Option<String>,
    pub x_ratio: Option<f64>,
    /// set on chrX/Y loci: how the copy number was decided (shown as a hint)
    pub sex_note: Option<Tr>,
}

pub const SRC_GROUPS: [&str; 7] =
    ["hifi_aligned", "hifi_not_in_aln", "fail_subread", "fail_lowqv", "fail_other", "unassigned_hifi", "unassigned_fail"];
pub const STATES: [&str; 11] = [
    "NORMAL",
    "CONTRACTED",
    "INTERMEDIATE",
    "INTERMEDIATE_LB",
    "EXPANDED",
    "EXPANDED_LB",
    "BENIGN_MOTIF",
    "NONREPEAT",
    "UNRESOLVED",
    "FLANK_ONLY",
    "AMBIGUOUS",
];

pub fn src_group(s: &str) -> &'static str {
    match s {
        "hifi_aligned" => "hifi_aligned",
        "hifi_not_in_aln" => "hifi_not_in_aln",
        "fail_subread" => "fail_subread",
        "fail_lowqv" => "fail_lowqv",
        "unassigned_hifi" => "unassigned_hifi",
        "unassigned_fail" => "unassigned_fail",
        _ => "fail_other",
    }
}

pub fn src_label(s: &str) -> (&'static str, &'static str) {
    match s {
        "hifi_aligned" => ("HiFi", "HiFi"),
        "hifi_not_in_aln" => ("HiFi(정렬본 밖)", "HiFi (not in aligned BAM)"),
        "fail_subread" => ("fail 0x08", "fail 0x08"),
        "fail_lowqv" => ("fail 0x01", "fail 0x01"),
        "unassigned_hifi" => ("unassigned HiFi", "unassigned HiFi"),
        "unassigned_fail" => ("unassigned fail", "unassigned fail"),
        _ => ("fail 기타", "fail other"),
    }
}

type Tally = BTreeMap<&'static str, BTreeMap<&'static str, usize>>;

fn f2(v: f64) -> String {
    if v.is_finite() { format!("{v:.2}") } else { "NA".into() }
}
fn of2(v: Option<f64>) -> String {
    v.map(f2).unwrap_or_else(|| "NA".into())
}
fn opt<T: std::fmt::Display>(v: &Option<T>) -> String {
    v.as_ref().map(|x| x.to_string()).unwrap_or_else(|| "NA".into())
}
pub fn fse(f: Option<f64>, se: Option<f64>) -> String {
    match (f, se) {
        (Some(f), Some(se)) => format!("{f:.2}±{se:.2}"),
        (Some(f), None) => format!("{f:.2}"),
        _ => "NA".into(),
    }
}
fn z(f: Option<f64>, se: Option<f64>) -> f64 {
    match (f, se) {
        (Some(f), Some(se)) if se > 0.0 => f / se,
        _ => f64::NAN,
    }
}

/// All per (sample, target) derived numbers, computed once.
pub struct Derived {
    pub tally: HashMap<(usize, usize), Tally>,
    pub ctl_mean: HashMap<(usize, usize), BTreeMap<&'static str, f64>>,
    pub n_ctl: HashMap<usize, usize>,
    pub alleles: HashMap<(usize, usize), AlleleSummary>,
    pub cohort: HashMap<usize, LocusCohort>,
    /// robust SD of depth loss f around each locus' cohort median (samples without a diagnosis there)
    pub null_sd: Option<f64>,
}

fn touches(t: Option<&Tally>, src: &str) -> usize {
    t.and_then(|m| m.get(src)).map(|x| x.iter().filter(|(k, _)| **k != "AMBIGUOUS").map(|(_, v)| *v).sum()).unwrap_or(0)
}
fn n_state(t: Option<&Tally>, src: &str, st: &str) -> usize {
    t.and_then(|m| m.get(src)).and_then(|x| x.get(st)).copied().unwrap_or(0)
}

pub fn derive(rows: &[ReadRow], sums: &[LocusSummary], loci: &[Locus], usable: &[bool]) -> Derived {
    let mut tally: HashMap<(usize, usize), Tally> = HashMap::new();
    for r in rows {
        *tally
            .entry((r.sample, r.call.locus))
            .or_default()
            .entry(src_group(&r.source))
            .or_default()
            .entry(r.call.state)
            .or_insert(0) += 1;
    }
    let mut ctls: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, l) in loci.iter().enumerate() {
        if let Some(t) = l.control_of {
            if usable[i] {
                ctls.entry(t).or_default().push(i);
            }
        }
    }
    let mut ctl_mean = HashMap::new();
    let mut by_sl: HashMap<(usize, usize), Vec<&ReadRow>> = HashMap::new();
    for r in rows {
        if !loci[r.call.locus].is_control() {
            by_sl.entry((r.sample, r.call.locus)).or_default().push(r);
        }
    }
    let mut alleles_m = HashMap::new();
    for s in sums {
        let cl = ctls.get(&s.locus).cloned().unwrap_or_default();
        let mut m = BTreeMap::new();
        if !cl.is_empty() {
            for src in SRC_GROUPS {
                let tot: usize = cl.iter().map(|c| touches(tally.get(&(s.sample, *c)), src)).sum();
                m.insert(src, tot as f64 / cl.len() as f64);
            }
        }
        ctl_mean.insert((s.sample, s.locus), m);
        let rr = by_sl.get(&(s.sample, s.locus)).cloned().unwrap_or_default();
        alleles_m.insert((s.sample, s.locus), alleles::summarize(&rr, &loci[s.locus]));
    }
    let n_ctl = ctls.iter().map(|(k, v)| (*k, v.len())).collect();
    Derived { tally, ctl_mean, n_ctl, alleles: alleles_m, cohort: HashMap::new(), null_sd: None }
}

/// Like `derive`, with control means taken from an earlier summary (control reads are not stored).
pub fn derive_from(
    rows: &[ReadRow],
    sums: &[LocusSummary],
    loci: &[Locus],
    ctl_mean: HashMap<(usize, usize), BTreeMap<&'static str, f64>>,
    n_ctl: HashMap<usize, usize>,
) -> Derived {
    let mut tally: HashMap<(usize, usize), Tally> = HashMap::new();
    let mut by_sl: HashMap<(usize, usize), Vec<&ReadRow>> = HashMap::new();
    for r in rows {
        *tally
            .entry((r.sample, r.call.locus))
            .or_default()
            .entry(src_group(&r.source))
            .or_default()
            .entry(r.call.state)
            .or_insert(0) += 1;
        by_sl.entry((r.sample, r.call.locus)).or_default().push(r);
    }
    let mut alleles_m = HashMap::new();
    for s in sums {
        let rr = by_sl.get(&(s.sample, s.locus)).cloned().unwrap_or_default();
        alleles_m.insert((s.sample, s.locus), alleles::summarize(&rr, &loci[s.locus]));
    }
    Derived { tally, ctl_mean, n_ctl, alleles: alleles_m, cohort: HashMap::new(), null_sd: None }
}

/// Evidence counts per source: (expanded exact, expanded lower bound, intermediate).
fn evidence_text(t: Option<&Tally>) -> Tr {
    let mut ko = Vec::new();
    let mut en = Vec::new();
    for src in SRC_GROUPS {
        let e = n_state(t, src, "EXPANDED");
        let lb = n_state(t, src, "EXPANDED_LB");
        let im = n_state(t, src, "INTERMEDIATE") + n_state(t, src, "INTERMEDIATE_LB");
        if e + lb + im > 0 {
            let (lk, le) = src_label(src);
            let mut pk = format!("{lk} 확장 {e}");
            let mut pe = format!("{le} expanded {e}");
            if lb > 0 {
                pk.push_str(&format!(" · 하한 {lb}"));
                pe.push_str(&format!(" · lower bound {lb}"));
            }
            if im > 0 {
                pk.push_str(&format!(" · 중간 {im}"));
                pe.push_str(&format!(" · intermediate {im}"));
            }
            ko.push(pk);
            en.push(pe);
        }
    }
    if ko.is_empty() { Tr::new("없음", "none") } else { Tr::new(ko.join(" / "), en.join(" / ")) }
}

/// Depth below which a sample's loci cannot be judged.
pub const LOW_DEPTH: f64 = 5.0;

/// Sample warning when the aligned BAM has almost no reads at the analysed loci.
pub fn low_depth_qc(sums: &[LocusSummary], infos: &mut [SampleInfo]) {
    for (si, info) in infos.iter_mut().enumerate() {
        let mut b: Vec<f64> = sums
            .iter()
            .filter(|s| s.sample == si)
            .filter_map(|s| s.canyon.as_ref().map(|c| if c.base_total.is_finite() { c.base_total } else { 0.0 }))
            .collect();
        if let Some(m) = med(&mut b) {
            if m < LOW_DEPTH && !info.qc.iter().any(|q| q.en.starts_with("Insufficient data")) {
                info.qc.push(Tr::new(
                    format!("데이터 부족: 정렬 BAM depth 중앙값 {m:.1}x — 이 샘플은 판정할 수 없음 ('이상 없음'이 음성을 뜻하지 않음)"),
                    format!("Insufficient data: median depth of the aligned BAM {m:.1}x — this sample cannot be judged ('no finding' does not mean negative)"),
                ));
            }
        }
    }
}

/// Warning for a positive sample whose known-diagnosis locus is missing or not analysed.
pub fn qc_positive_without_locus(group: &str) -> Tr {
    Tr::new(
        format!("group '{group}' = 양성인데 locus(기존 진단 좌위)가 비어 있거나 분석 좌위에 없음"),
        format!("group '{group}' = positive, but the locus column (known diagnosis) is empty or not among the analysed loci"),
    )
}

/// Sample warning when the declared sex disagrees with the chrX depth.
pub fn qc_sex_mismatch(declared: &str, inferred: &str, ratio: f64) -> Tr {
    Tr::new(
        format!(
            "성별 불일치: manifest {} / chrX depth 비 {ratio:.2} → {} 추정 (샘플 바뀜·성염색체 이수성 확인)",
            if declared == "M" { "남" } else { "여" },
            if inferred == "M" { "X 1개" } else { "X 2개" }
        ),
        format!(
            "Sex mismatch: manifest {declared} / chrX depth ratio {ratio:.2} → {} (check for a sample swap or sex-chromosome aneuploidy)",
            if inferred == "M" { "one X" } else { "two X" }
        ),
    )
}

/// Hint for chrX/Y loci: how the X copy number was decided.
pub fn sex_note(declared: Option<&str>, inferred: Option<&str>, ratio: Option<f64>) -> Tr {
    let r = ratio.unwrap_or(f64::NAN);
    let (ko, en) = match (declared, inferred) {
        (Some("M"), _) => ("남성 (manifest)".to_string(), "male (manifest)".to_string()),
        (Some("F"), Some("M")) => (format!("manifest 여성이나 chrX depth 비 {r:.2}"), format!("manifest female but chrX depth ratio {r:.2}")),
        (Some("F"), _) => ("여성 (manifest)".to_string(), "female (manifest)".to_string()),
        (_, Some("M")) => (format!("남성 추정 (chrX depth 비 {r:.2})"), format!("male by chrX depth (ratio {r:.2})")),
        (_, Some("F")) => (format!("여성 추정 (chrX depth 비 {r:.2})"), format!("female by chrX depth (ratio {r:.2})")),
        _ => ("성별 미상".to_string(), "sex unknown".to_string()),
    };
    if declared == Some("M") || inferred == Some("M") {
        Tr::new(format!("{ko}: X 1개 → HP 비교 없음"), format!("{en}: one X → no haplotype comparison"))
    } else {
        Tr::new(ko, en)
    }
}

pub struct Verdict {
    pub level: &'static str,
    pub hints: Vec<Tr>,
}

/// One piece of evidence behind a verdict. `kind` groups signals across samples for the cohort check.
#[derive(Clone)]
struct Signal {
    kind: &'static str,
    text: Tr,
    level: &'static str,
}

/// Per-locus cohort statistics over samples without a known diagnosis at that locus.
#[derive(Default, Clone)]
pub struct LocusCohort {
    pub n_eval: usize,
    pub kind_n: HashMap<&'static str, usize>,
    pub f_med: Option<f64>,
    pub hifi_ratio_med: Option<f64>,
}

/// Cohort comparison is used when at least this many samples (without a known diagnosis there) were analysed.
pub const COHORT_MIN: usize = 4;

fn rank(x: &str) -> u8 {
    match x {
        "high" => 3,
        "watch" => 2,
        "info" => 1,
        _ => 0,
    }
}

/// Upper tail P(X >= k) for X ~ Poisson(lambda).
fn pois_tail(k: usize, lambda: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    let mut term = (-lambda).exp();
    let mut cdf = term;
    for i in 1..k {
        term *= lambda / i as f64;
        cdf += term;
    }
    (1.0 - cdf).max(0.0)
}

fn med(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(v[v.len() / 2])
}

fn signals(s: &LocusSummary, d: &Derived, l: &Locus, co: Option<&LocusCohort>) -> Vec<Signal> {
    let key = (s.sample, s.locus);
    let t = d.tally.get(&key);
    let a = d.alleles.get(&key);
    let cm = d.ctl_mean.get(&key);
    let co = co.filter(|c| c.n_eval >= COHORT_MIN);
    let mut out: Vec<Signal> = Vec::new();
    let mut push = |kind: &'static str, text: Tr, level: &'static str| out.push(Signal { kind, text, level });
    if let Some(a) = a {
        match a.range {
            "pathogenic" => push(
                "range_path",
                if l.contraction {
                    Tr::new("HiFi에서 병적 범위 allele (축소)", "pathogenic-range allele in HiFi (contraction)")
                } else {
                    Tr::new("HiFi에서 병적 범위 allele", "pathogenic-range allele in HiFi")
                },
                "high",
            ),
            "intermediate" if l.ref_over_normal => push(
                "range_int",
                Tr::new(
                    format!(
                        "HiFi에서 중간 범위 allele이지만 참조 allele({:.0}회)부터 정상 상한({:.0}회)을 넘음 → 구간·기준 불일치 가능",
                        l.units((l.end - l.start) as i64),
                        l.normal_max_units
                    ),
                    format!(
                        "intermediate-range allele in HiFi, but the reference allele ({:.0} units) already exceeds the normal maximum ({:.0}) → region/threshold mismatch likely",
                        l.units((l.end - l.start) as i64),
                        l.normal_max_units
                    ),
                ),
                "info",
            ),
            "intermediate" => push("range_int", Tr::new("HiFi에서 중간 범위 allele", "intermediate-range allele in HiFi"), "watch"),
            "benign_motif" => push("benign_motif", Tr::new("긴 allele이지만 병적 motif 아님", "long allele, but not the pathogenic motif"), "info"),
            _ => {}
        }
        let hifi_exp = a.clusters.iter().any(|c| matches!(c.class, "pathogenic" | "intermediate"));
        if a.fail_expanded > 0 && !hifi_exp {
            let lv = match a.fail_expanded {
                0..=2 => "info",
                3..=4 => "watch",
                _ => "high",
            };
            push(
                "fail_only",
                Tr::new(
                    format!("확장 근거가 정렬본 밖 read(fail·unassigned)에만 있음 ({}개)", a.fail_expanded),
                    format!("expansion seen only in reads outside the aligned BAM (fail/unassigned: {})", a.fail_expanded),
                ),
                lv,
            );
        }
    }
    // read counts relative to neighbourhood controls
    let n_ctl = d.n_ctl.get(&s.locus).copied().unwrap_or(0);
    if let (Some(cm), true) = (cm, n_ctl >= 3) {
        for (src, kind) in [
            ("fail_subread", "enrich_fail_subread"),
            ("fail_lowqv", "enrich_fail_lowqv"),
            ("fail_other", "enrich_fail_other"),
            ("unassigned_hifi", "enrich_unassigned_hifi"),
            ("unassigned_fail", "enrich_unassigned_fail"),
        ] {
            let tt = touches(t, src);
            let m = cm.get(src).copied().unwrap_or(0.0);
            if tt >= 4 && tt as f64 >= 3.0 * m.max(0.5) && pois_tail(tt, m.max(0.5)) < 1e-3 {
                let (lk, le) = src_label(src);
                push(
                    kind,
                    Tr::new(format!("{lk} read 농축 ({tt} vs 이웃 평균 {m:.1})"), format!("{le} reads enriched ({tt} vs neighbour mean {m:.1})")),
                    "watch",
                );
            }
        }
        let th = touches(t, "hifi_aligned");
        let mh = cm.get("hifi_aligned").copied().unwrap_or(0.0);
        let ratio = if mh > 0.0 { th as f64 / mh } else { f64::NAN };
        let thr = 0.6 * co.and_then(|c| c.hifi_ratio_med).unwrap_or(1.0).min(1.0);
        if mh >= 10.0 && ratio < thr {
            push("hifi_reduced", Tr::new(format!("HiFi read 감소 ({th} vs 이웃 평균 {mh:.0})"), format!("fewer HiFi reads ({th} vs neighbour mean {mh:.0})")), "watch");
        }
    }
    if let Some(c) = &s.canyon {
        let (f1, f2v) = (c.f_hp1.unwrap_or(f64::NAN), c.f_hp2.unwrap_or(f64::NAN));
        let (z1, z2, zt) = (z(c.f_hp1, c.se_hp1), z(c.f_hp2, c.se_hp2), z(c.f_total, c.se_total));
        let ft = c.f_total.unwrap_or(f64::NAN);
        // Depth loss must stand out from the other samples at this locus (GC-rich dips are shared) and from the
        // cohort's own spread: bootstrap SEs miss correlated depth noise, so z alone over-calls (in-house cohort:
        // 1% of loci reach z 3.6 with nothing there). Without a cohort, fall back to f and z only.
        let fmed = co.and_then(|c| c.f_med).unwrap_or(0.0);
        let delta = ft - fmed;
        let (strong, weak) = match d.null_sd {
            Some(sd) if co.is_some() => (delta >= (4.0 * sd).max(0.7) && zt >= 3.0, delta >= (3.0 * sd).max(0.5) && zt >= 2.5),
            _ => (ft >= 0.7 && zt >= 4.0, ft >= 0.5 && zt >= 3.0),
        };
        let asym = (weak || delta >= 0.3)
            && ((f1 >= 0.8 && z1 >= 3.5 && f2v < 0.25) || (f2v >= 0.8 && z2 >= 3.5 && f1 < 0.25));
        if asym {
            push("depth_asym", Tr::new("depth: 한 haplotype만 함몰 (한쪽 확장 패턴)", "depth: one haplotype lost (one-allele expansion pattern)"), "watch");
        } else if strong && f1 >= 0.5 && f2v >= 0.5 && z1 >= 2.0 && z2 >= 2.0 {
            push("depth_both", Tr::new("depth: 두 haplotype 함께 함몰 (양쪽 확장 패턴)", "depth: both haplotypes lost (biallelic expansion pattern)"), "watch");
        } else if strong {
            push("depth_total", Tr::new("depth: 총 depth 함몰", "depth: total depth loss"), "watch");
        } else if weak {
            push("depth_total", Tr::new("depth: 약한 함몰 (30x에서는 우연히도 나오는 크기)", "depth: weak dip (within chance at 30x)"), "info");
        }
        let any_exp = a.map(|a| a.fail_expanded > 0 || a.clusters.iter().any(|c| c.class != "normal")).unwrap_or(false);
        if strong && !any_exp {
            push("absence", Tr::new("확장 read 없이 depth만 크게 함몰 (부재 패턴 후보)", "strong depth loss without expanded reads (candidate absence pattern)"), "high");
        }
    }
    out
}

/// Cohort statistics per target locus (samples with a known diagnosis at that locus are left out).
pub fn add_cohort(d: &mut Derived, sums: &[LocusSummary], loci: &[Locus], samples: &[SampleInfo]) {
    let mut by_locus: HashMap<usize, Vec<&LocusSummary>> = HashMap::new();
    for s in sums.iter().filter(|s| s.error.is_none() && s.canyon.as_ref().is_some_and(|c| c.base_total >= LOW_DEPTH)) {
        if !samples[s.sample].expected.contains(&s.locus) {
            by_locus.entry(s.locus).or_default().push(s);
        }
    }
    let mut out: HashMap<usize, LocusCohort> = HashMap::new();
    for (li, ss) in by_locus {
        let l = &loci[li];
        let mut c = LocusCohort { n_eval: ss.len(), ..Default::default() };
        let mut f: Vec<f64> = ss.iter().filter_map(|s| s.canyon.as_ref().and_then(|c| c.f_total)).filter(|v| v.is_finite()).collect();
        c.f_med = med(&mut f);
        let mut r: Vec<f64> = ss
            .iter()
            .filter_map(|s| {
                let m = d.ctl_mean.get(&(s.sample, s.locus))?.get("hifi_aligned").copied()?;
                (m >= 10.0).then(|| touches(d.tally.get(&(s.sample, s.locus)), "hifi_aligned") as f64 / m)
            })
            .collect();
        c.hifi_ratio_med = med(&mut r);
        for s in &ss {
            let mut kinds: Vec<&'static str> = signals(s, d, l, None).into_iter().map(|g| g.kind).collect();
            kinds.sort_unstable();
            kinds.dedup();
            for k in kinds {
                *c.kind_n.entry(k).or_default() += 1;
            }
        }
        out.insert(li, c);
    }
    // cohort-wide spread of f around the per-locus median
    let mut res: Vec<f64> = Vec::new();
    for s in sums.iter().filter(|s| s.error.is_none() && !samples[s.sample].expected.contains(&s.locus)) {
        if let (Some(f), Some(c)) = (s.canyon.as_ref().and_then(|c| c.f_total), out.get(&s.locus)) {
            if let (Some(m), true) = (c.f_med, c.n_eval >= COHORT_MIN && f.is_finite()) {
                res.push((f - m).abs());
            }
        }
    }
    d.null_sd = (res.len() >= 40).then(|| 1.4826 * med(&mut res).unwrap_or(0.0)).filter(|v| *v > 0.0);
    d.cohort = out;
}

pub fn verdict(s: &LocusSummary, d: &Derived, l: &Locus) -> Verdict {
    let co = d.cohort.get(&s.locus);
    let mut h: Vec<Tr> = Vec::new();
    let mut level = "ok";
    if let Some(n) = &s.sex_note {
        h.push(n.clone());
    }
    if let Some(c) = &s.canyon {
        if !(c.base_total >= LOW_DEPTH) {
            let x = if c.base_total.is_finite() { c.base_total } else { 0.0 };
            h.push(Tr::new(format!("depth {x:.1}x — 데이터 부족으로 판정 불가"), format!("depth {x:.1}x — not enough data to judge")));
            return Verdict { level: "info", hints: h };
        }
    }
    for g in signals(s, d, l, co) {
        // a signal seen in at least half of the other samples at this locus is a locus property, not this sample's
        let common = co.filter(|c| c.n_eval >= COHORT_MIN).and_then(|c| {
            let k = c.kind_n.get(g.kind).copied().unwrap_or(0);
            (k * 2 >= c.n_eval).then_some((k, c.n_eval))
        });
        let lv = if common.is_some() && rank(g.level) > 1 { "info" } else { g.level };
        if rank(lv) > rank(level) {
            level = lv;
        }
        match common {
            Some((k, n)) => h.push(Tr::new(
                format!("{} [다른 샘플 {k}/{n}명도 같음 → 좌위 특성 가능성]", g.text.ko),
                format!("{} [also in {k}/{n} other samples → likely a property of the locus]", g.text.en),
            )),
            None => h.push(g.text),
        }
    }
    if !l.thresholds_known {
        h.push(Tr::new("이 좌위는 출처에 정상/병적 기준 일부가 없음", "the source gives only part of the normal/pathogenic thresholds for this locus"));
    }
    Verdict { level, hints: h }
}

// ------------------------------------------------------------------ TSV

pub fn write_reads(path: &Path, rows: &[ReadRow], samples: &[SampleInfo], loci: &[Locus]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "sample_id\tlocus\tgene\tsource\tbam\tread_name\tread_len\tff\trq\tnp\thp\tstrand\tclass\tallele_bp\tbound\trepeat_units\tstate\tdominant_motif\tmotif_frac\tpurity\tpathogenic_motif\thits_left\thits_right\tflags"
    )?;
    for r in rows.iter().filter(|r| !loci[r.call.locus].is_control()) {
        let l = &loci[r.call.locus];
        let units = r.call.allele_bp.map(|b| format!("{:.1}", l.units(b))).unwrap_or_else(|| "NA".into());
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}",
            samples[r.sample].id,
            l.id,
            l.gene,
            r.source,
            r.bam,
            r.name,
            r.read_len,
            opt(&r.ff),
            r.rq.map(|x| format!("{x:.4}")).unwrap_or_else(|| "NA".into()),
            opt(&r.np),
            opt(&r.hp),
            r.call.strand,
            r.call.cls.as_str(),
            opt(&r.call.allele_bp),
            r.call.bound.as_str(),
            units,
            r.call.state,
            r.call.dom_motif,
            r.call.dom_frac,
            r.call.purity,
            if r.call.dom_motif == "-" { "NA" } else if r.call.path_motif { "yes" } else { "no" },
            r.call.hits_l,
            r.call.hits_r,
            if r.call.flags.is_empty() { "-".to_string() } else { r.call.flags.join(",") },
        )?;
    }
    Ok(())
}

pub fn write_depth(path: &Path, sums: &[LocusSummary], samples: &[SampleInfo], loci: &[Locus]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "sample_id\tlocus\toffset_bp\trel_total\trel_hp1\trel_hp2\ttemplate_g\tif_one_allele_lost\tif_both_alleles_lost")?;
    for s in sums {
        let Some(c) = &s.canyon else { continue };
        for i in 0..c.offsets.len() {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                samples[s.sample].id,
                loci[s.locus].id,
                c.offsets[i],
                f2(c.rel_total[i]),
                f2(c.rel_hp1[i]),
                f2(c.rel_hp2[i]),
                c.g[i],
                1.0 - 0.5 * c.g[i],
                1.0 - c.g[i]
            )?;
        }
    }
    Ok(())
}

pub fn write_bam_stats(path: &Path, stats: &[BamStat], samples: &[SampleInfo]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "sample_id\tkind\tbam\treads\tbases\tprefilter_pass\tplaced\tseconds\tcategory_reads")?;
    for b in stats {
        let cats: Vec<String> = b.stats.by_cat.iter().map(|(k, v)| format!("{k}:{}", v.reads)).collect();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
            samples[b.sample].id,
            b.kind,
            b.path,
            b.stats.reads,
            b.stats.bases,
            b.stats.prefilter,
            b.stats.candidates,
            b.stats.seconds,
            cats.join(",")
        )?;
    }
    Ok(())
}

pub fn write_summary(path: &Path, sums: &[LocusSummary], samples: &[SampleInfo], loci: &[Locus], d: &Derived, lang: &str) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    let mut head = String::from(
        "sample_id\tgroup\trole\tsex\tsex_inferred\texpected\tlocus\tgene\tdisease\tlevel\thifi_alleles\trange\thifi_lb_max_bp\tfail_expanded\tchrX_auto_ratio\tbaseline_depth\tcentre_rel\tf_total\tf_hp1\tf_hp2\tse_total\tse_hp1\tse_hp2\tphased_frac\tread_n50\tusable_controls",
    );
    for src in SRC_GROUPS {
        for st in ["NORMAL", "INTERMEDIATE", "EXPANDED", "EXPANDED_LB", "BENIGN_MOTIF", "NONREPEAT", "FLANK_ONLY"] {
            head.push_str(&format!("\t{src}.{st}"));
        }
        head.push_str(&format!("\t{src}.touch\t{src}.ctl_mean\t{src}.rel_to_ctl"));
    }
    head.push_str("\tevidence\thint\terror");
    writeln!(w, "{head}")?;
    for s in sums {
        let si = &samples[s.sample];
        let l = &loci[s.locus];
        let key = (s.sample, s.locus);
        let t = d.tally.get(&key);
        let a = d.alleles.get(&key).cloned().unwrap_or_default();
        let v = verdict(s, d, l);
        let c = s.canyon.as_ref();
        let mut line = vec![
            si.id.clone(),
            si.group.clone(),
            si.role.into(),
            si.sex_declared.unwrap_or("NA").into(),
            si.sex_inferred.unwrap_or("NA").into(),
            if si.expected.contains(&s.locus) { "yes".into() } else { "no".into() },
            l.id.clone(),
            l.gene.clone(),
            l.disease.clone(),
            v.level.into(),
            alleles::text(&a).get(lang).to_string(),
            a.range.into(),
            opt(&a.hifi_lb_max),
            a.fail_expanded.to_string(),
            of2(s.x_ratio),
            c.map(|c| f2(c.base_total)).unwrap_or("NA".into()),
            c.map(|c| f2(c.centre_total / c.base_total)).unwrap_or("NA".into()),
            of2(c.and_then(|c| c.f_total)),
            of2(c.and_then(|c| c.f_hp1)),
            of2(c.and_then(|c| c.f_hp2)),
            of2(c.and_then(|c| c.se_total)),
            of2(c.and_then(|c| c.se_hp1)),
            of2(c.and_then(|c| c.se_hp2)),
            f2(s.phased_frac),
            c.map(|c| c.read_n50.to_string()).unwrap_or("NA".into()),
            d.n_ctl.get(&s.locus).copied().unwrap_or(0).to_string(),
        ];
        let cm = d.ctl_mean.get(&key);
        for src in SRC_GROUPS {
            let im = n_state(t, src, "INTERMEDIATE") + n_state(t, src, "INTERMEDIATE_LB");
            for st in ["NORMAL", "INTERMEDIATE", "EXPANDED", "EXPANDED_LB", "BENIGN_MOTIF", "NONREPEAT", "FLANK_ONLY"] {
                line.push(if st == "INTERMEDIATE" { im.to_string() } else { n_state(t, src, st).to_string() });
            }
            let tt = touches(t, src);
            let m = cm.and_then(|m| m.get(src)).copied().unwrap_or(f64::NAN);
            line.push(tt.to_string());
            line.push(f2(m));
            line.push(if m > 0.0 { f2(tt as f64 / m) } else { "NA".into() });
        }
        line.push(evidence_text(t).get(lang).to_string());
        line.push(if v.hints.is_empty() { "-".into() } else { v.hints.iter().map(|h| h.get(lang)).collect::<Vec<_>>().join("; ") });
        line.push(s.error.clone().unwrap_or_else(|| "-".into()));
        writeln!(w, "{}", line.join("\t"))?;
    }
    Ok(())
}

// ------------------------------------------------------------------ JSON for the HTML report

#[derive(Serialize)]
struct JMeta<'a> {
    version: &'a str,
    created: String,
    reference: String,
    params: BTreeMap<&'a str, String>,
    catalog: &'a str,
    /// default language of the page (ko / en); the reader can switch
    lang: &'a str,
}

#[derive(Serialize)]
struct JSample<'a> {
    id: &'a str,
    group: &'a str,
    role: &'a str,
    sex: Option<&'a str>,
    sex_inferred: Option<&'a str>,
    x_ratio: Option<f64>,
    qc: &'a [Tr],
    expected: &'a [usize],
    expected_raw: &'a str,
    note: &'a str,
    hifi_bam: &'a str,
    fail_bams: &'a [String],
    hifi_raw_bams: &'a [String],
    notes: &'a [String],
    unassigned: &'a [String],
    unassigned_used: &'a [String],
}

#[derive(Serialize)]
struct JLocus<'a> {
    id: &'a str,
    gene: &'a str,
    disease: &'a str,
    disease_id: &'a str,
    region: String,
    motifs: &'a [String],
    path_motifs: &'a [String],
    unit: usize,
    /// non-repeat bp inside the region: units = (bp - offset) / unit
    offset: i64,
    count_note: &'a Tr,
    ref_edge: Option<(i64, i64)>,
    normal_max: f64,
    path_min: f64,
    path_max: Option<f64>,
    known: bool,
    inheritance: &'a str,
    evidence: &'a str,
    location: &'a str,
    source: &'a str,
    stripy: &'a str,
    gnomad: &'a str,
    omim: &'a str,
    genereviews: &'a str,
    uniq_l: f64,
    uniq_r: f64,
    alt: f64,
    controls: usize,
}

#[derive(Serialize)]
struct JCanyon {
    step: i64,
    x0: i64,
    tot: Vec<Option<f32>>,
    h1: Vec<Option<f32>>,
    h2: Vec<Option<f32>>,
    g: Vec<f32>,
    base: f64,
    base_h1: Option<f64>,
    base_h2: Option<f64>,
    centre_rel: Option<f64>,
    f: Option<f64>,
    f1: Option<f64>,
    f2: Option<f64>,
    se: Option<f64>,
    se1: Option<f64>,
    se2: Option<f64>,
    phased: f64,
    n50: u64,
}

#[derive(Serialize)]
struct JResult<'a> {
    s: usize,
    l: usize,
    expected: bool,
    level: &'a str,
    hints: Vec<Tr>,
    alleles: AlleleSummary,
    allele_text: Tr,
    evidence: Tr,
    x_ratio: Option<f64>,
    error: Option<String>,
    canyon: Option<JCanyon>,
    counts: BTreeMap<&'static str, BTreeMap<&'static str, usize>>,
    ctl: BTreeMap<&'static str, f64>,
    n_ctl: usize,
    /// [source, class, bp, bound, hp, purity%, motif, state, name?, read_len]
    reads: Vec<serde_json::Value>,
}

#[derive(Serialize)]
struct JBam<'a> {
    sample: usize,
    kind: &'a str,
    path: &'a str,
    reads: u64,
    gb: f64,
    prefilter: u64,
    placed: u64,
    seconds: f64,
    cats: BTreeMap<String, u64>,
}

#[derive(Serialize)]
struct JData<'a> {
    meta: JMeta<'a>,
    samples: Vec<JSample<'a>>,
    loci: Vec<JLocus<'a>>,
    bams: Vec<JBam<'a>>,
    results: Vec<JResult<'a>>,
}

fn round(v: f64) -> Option<f32> {
    v.is_finite().then(|| (v * 100.0).round() as f32 / 100.0)
}

/// Average pairs of bins to a coarser display grid.
fn coarse(c: &Canyon, factor: usize) -> JCanyon {
    let n = c.offsets.len() / factor;
    let avg = |v: &[f64]| -> Vec<Option<f32>> {
        (0..n)
            .map(|i| {
                let s: Vec<f64> = v[i * factor..(i + 1) * factor].iter().copied().filter(|x| x.is_finite()).collect();
                if s.is_empty() { None } else { round(s.iter().sum::<f64>() / s.len() as f64) }
            })
            .collect()
    };
    let step = if c.offsets.len() > 1 { (c.offsets.last().unwrap() - c.offsets[0]) / (c.offsets.len() as i64 - 1) } else { 500 };
    JCanyon {
        step: step * factor as i64,
        x0: c.offsets[0] + step * (factor as i64 - 1) / 2,
        tot: avg(&c.rel_total),
        h1: avg(&c.rel_hp1),
        h2: avg(&c.rel_hp2),
        g: (0..n).map(|i| (c.g[i * factor..(i + 1) * factor].iter().sum::<f64>() / factor as f64) as f32).collect(),
        base: c.base_total,
        base_h1: c.base_hp1.is_finite().then_some(c.base_hp1),
        base_h2: c.base_hp2.is_finite().then_some(c.base_hp2),
        centre_rel: (c.base_total > 0.0).then(|| c.centre_total / c.base_total),
        f: c.f_total,
        f1: c.f_hp1,
        f2: c.f_hp2,
        se: c.se_total,
        se1: c.se_hp1,
        se2: c.se_hp2,
        phased: 0.0,
        n50: c.read_n50,
    }
}

#[allow(clippy::too_many_arguments)]
pub fn write_html(
    path: &Path,
    sums: &[LocusSummary],
    rows: &[ReadRow],
    samples: &[SampleInfo],
    loci: &[Locus],
    usable: &[bool],
    panel_frac: &[(f64, f64)],
    alt_frac: &[f64],
    bam_stats: &[BamStat],
    d: &Derived,
    version: &str,
    reference: &str,
    params: BTreeMap<&str, String>,
    lang: &str,
) -> Result<()> {
    let targets: Vec<usize> = (0..loci.len()).filter(|&i| !loci[i].is_control()).collect();
    let tindex: HashMap<usize, usize> = targets.iter().enumerate().map(|(j, &i)| (i, j)).collect();
    let jl: Vec<JLocus> = targets
        .iter()
        .map(|&i| {
            let l = &loci[i];
            JLocus {
                id: &l.id,
                gene: &l.gene,
                disease: &l.disease,
                disease_id: &l.disease_id,
                region: l.region_str(),
                motifs: &l.motifs,
                path_motifs: &l.path_motifs,
                unit: l.unit_len(),
                offset: l.count_offset_bp,
                count_note: &l.count_note,
                ref_edge: l.ref_edge_bp,
                normal_max: l.normal_max_units,
                path_min: l.path_min_units,
                path_max: l.path_max_units,
                known: l.thresholds_known,
                inheritance: &l.inheritance,
                evidence: &l.evidence,
                location: &l.location,
                source: &l.source,
                stripy: &l.stripy,
                gnomad: &l.gnomad,
                omim: &l.omim,
                genereviews: &l.genereviews,
                uniq_l: panel_frac[i].0,
                uniq_r: panel_frac[i].1,
                alt: alt_frac[i],
                controls: (0..loci.len()).filter(|&c| loci[c].control_of == Some(i) && usable[c]).count(),
            }
        })
        .collect();
    let mut by_sl: HashMap<(usize, usize), Vec<&ReadRow>> = HashMap::new();
    for r in rows {
        if !loci[r.call.locus].is_control() && r.call.allele_bp.is_some() {
            by_sl.entry((r.sample, r.call.locus)).or_default().push(r);
        }
    }
    let mut results = Vec::new();
    for s in sums {
        let key = (s.sample, s.locus);
        let l = &loci[s.locus];
        let v = verdict(s, d, l);
        let a = d.alleles.get(&key).cloned().unwrap_or_default();
        let reads: Vec<serde_json::Value> = by_sl
            .get(&key)
            .map(|rs| {
                rs.iter()
                    .map(|r| {
                        let notable = !matches!(r.call.state, "NORMAL" | "FLANK_ONLY" | "AMBIGUOUS");
                        serde_json::json!([
                            r.source,
                            r.call.cls.as_str(),
                            r.call.allele_bp,
                            r.call.bound.as_str(),
                            r.hp,
                            (r.call.purity * 100.0).round() as i64,
                            r.call.dom_motif,
                            r.call.state,
                            if notable { r.name.as_str() } else { "" },
                            r.read_len
                        ])
                    })
                    .collect()
            })
            .unwrap_or_default();
        let mut canyon = s.canyon.as_ref().map(|c| coarse(c, 2));
        if let Some(c) = canyon.as_mut() {
            c.phased = s.phased_frac;
        }
        results.push(JResult {
            s: s.sample,
            l: tindex[&s.locus],
            expected: samples[s.sample].expected.contains(&s.locus),
            level: v.level,
            hints: v.hints,
            allele_text: alleles::text(&a),
            alleles: a,
            evidence: evidence_text(d.tally.get(&key)),
            x_ratio: s.x_ratio,
            error: s.error.clone(),
            canyon,
            counts: d.tally.get(&key).cloned().unwrap_or_default(),
            ctl: d.ctl_mean.get(&key).cloned().unwrap_or_default(),
            n_ctl: d.n_ctl.get(&s.locus).copied().unwrap_or(0),
            reads,
        });
    }
    let expected_t: Vec<Vec<usize>> =
        samples.iter().map(|si| si.expected.iter().filter_map(|i| tindex.get(i).copied()).collect()).collect();
    let data = JData {
        meta: JMeta {
            lang,
            version,
            created: crate::now_string(),
            reference: reference.to_string(),
            params,
            catalog: "STRchive snapshot 2026-11-17 (CC BY 4.0)",
        },
        samples: samples
            .iter()
            .enumerate()
            .map(|(i, si)| JSample {
                id: &si.id,
                group: &si.group,
                role: si.role,
                sex: si.sex_declared,
                sex_inferred: si.sex_inferred,
                x_ratio: si.x_ratio.filter(|v| v.is_finite()),
                qc: &si.qc,
                expected: &expected_t[i],
                expected_raw: &si.expected_raw,
                note: &si.note,
                hifi_bam: &si.hifi_bam,
                fail_bams: &si.fail_bams,
                hifi_raw_bams: &si.hifi_raw_bams,
                notes: &si.locate_notes,
                unassigned: &si.unassigned,
                unassigned_used: &si.unassigned_used,
            })
            .collect(),
        loci: jl,
        bams: bam_stats
            .iter()
            .map(|b| JBam {
                sample: b.sample,
                kind: &b.kind,
                path: &b.path,
                reads: b.stats.reads,
                gb: b.stats.bases as f64 / 1e9,
                prefilter: b.stats.prefilter,
                placed: b.stats.candidates,
                seconds: b.stats.seconds,
                cats: b.stats.by_cat.iter().map(|(k, v)| (k.clone(), v.reads)).collect(),
            })
            .collect(),
        results,
    };
    let json = serde_json::to_string(&data)?.replace("</", "<\\/");
    let html = TEMPLATE
        .replace("__ICON__", &b64(ICON))
        .replace("__BANNER__", &b64(BANNER))
        .replace("__EGG__", &b64(EGG))
        .replace("__FAVICON__", &b64(FAVICON))
        .replace("__TR_RESCUE_DATA__", &json);
    std::fs::write(path, html)?;
    Ok(())
}
