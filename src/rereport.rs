// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! `tr-rescue report`: rebuild summary / reads / report.html from an existing `trace` output folder
//! with the current interpretation rules (thresholds, repeat-count corrections, cohort checks),
//! without scanning any BAM again. Read placement and lengths are taken as they were measured.

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use serde_json::Value;

use crate::anchor::{Bound, Call, Cls, finalize_state};
use crate::catalog::{self, Locus, add_controls, builtin, parse_bed};
use crate::depth::Canyon;
use crate::fasta::Fasta;
use crate::report::{self, BamStat, Derived, LocusSummary, ReadRow, SRC_GROUPS, SampleInfo};
use crate::scan::{CatStat, ScanStats};

fn tsv(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = text.lines();
    let head: Vec<String> = lines.next().ok_or_else(|| anyhow!("{} is empty", path.display()))?.split('\t').map(String::from).collect();
    let rows = lines.filter(|l| !l.is_empty()).map(|l| l.split('\t').map(String::from).collect()).collect();
    Ok((head, rows))
}

fn col(head: &[String], name: &str) -> Option<usize> {
    head.iter().position(|h| h == name)
}

fn num(s: &str) -> Option<f64> {
    s.parse::<f64>().ok().filter(|v| v.is_finite())
}

fn st(s: &str) -> &'static str {
    match s {
        "positive" => "positive",
        "control" => "control",
        "M" => "M",
        "F" => "F",
        _ => "test",
    }
}

fn flag(s: &str) -> Option<&'static str> {
    ["BOTH_STRANDS", "FEW_HITS", "FLANK_ORDER", "FLANK_REACH", "LOW_PURITY_INSERT", "NO_OPPOSITE_FLANK"].into_iter().find(|f| *f == s)
}

fn strs(v: &Value) -> Vec<String> {
    v.as_array().map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default()
}

/// The JSON block embedded in report.html.
fn report_json(dir: &Path) -> Result<Value> {
    let p = dir.join("report.html");
    let html = fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?;
    let tag = "<script id=\"data\" type=\"application/json\">";
    let a = html.find(tag).ok_or_else(|| anyhow!("{}: no data block (not a tr-rescue report?)", p.display()))? + tag.len();
    let b = html[a..].find("</script>").ok_or_else(|| anyhow!("{}: data block not closed", p.display()))? + a;
    serde_json::from_str(&html[a..b]).with_context(|| format!("parsing data in {}", p.display()))
}

pub fn run(from: &Path, out: &Path, catalog_bed: Option<&Path>, reference: Option<&Path>, version: &str, lang: &str) -> Result<()> {
    if from == out {
        bail!("--out must differ from --from (the original results are kept)");
    }
    let d = report_json(from)?;
    let meta = &d["meta"];
    let version_from = meta["version"].as_str().unwrap_or("?").to_string();
    eprintln!("[report] rebuilding {} (trace v{version_from}) with v{version}", from.display());

    // ---- loci: current catalog definitions, in the order of the old report
    let mut cat = builtin();
    if let Some(p) = catalog_bed {
        let extra = parse_bed(p)?;
        cat.retain(|b| !extra.iter().any(|e| e.id.eq_ignore_ascii_case(&b.id)));
        cat.extend(extra);
    }
    let fa = match reference {
        Some(p) => Some(Fasta::open(p)?),
        None => None,
    };
    let jl = d["loci"].as_array().ok_or_else(|| anyhow!("report has no loci"))?;
    let mut targets: Vec<Locus> = Vec::new();
    let mut panel_frac = Vec::new();
    let mut alt_frac = Vec::new();
    let mut n_ctl_json = Vec::new();
    for x in jl {
        let id = x["id"].as_str().unwrap_or("");
        let Some(mut l) = cat.iter().find(|c| c.id == id).cloned() else {
            bail!("locus {id} of the old report is not in the current catalog (give the same --catalog)");
        };
        // keep the contig name used in that run; warn when the region itself changed
        let region = x["region"].as_str().unwrap_or("");
        if let Some((c, rng)) = region.split_once(':') {
            l.chrom = c.to_string();
            if let Some((a, b)) = rng.split_once('-') {
                let (a, b): (u64, u64) = (a.parse().unwrap_or(0), b.parse().unwrap_or(0));
                if a > 0 && (a - 1 != l.start || b != l.end) {
                    eprintln!(
                        "[report] WARNING {id}: region was {region}, current catalog has {}; lengths stay as measured, counts may be off — run trace again for this locus",
                        l.region_str()
                    );
                }
            }
        }
        if let Some(fa) = &fa {
            if let Some(c) = fa.resolve(&l.chrom) {
                if let Ok(seq) = fa.fetch(&c, l.start, l.end) {
                    l.ref_edge_bp = catalog::edge_nonrepeat(&seq, &l.motifs, l.unit_len());
                }
            }
        }
        panel_frac.push((x["uniq_l"].as_f64().unwrap_or(f64::NAN), x["uniq_r"].as_f64().unwrap_or(f64::NAN)));
        alt_frac.push(x["alt"].as_f64().unwrap_or(0.0));
        n_ctl_json.push(x["controls"].as_u64().unwrap_or(0) as usize);
        targets.push(l);
    }
    let loci = add_controls(&targets);
    let nt = targets.len();
    // usable flags: only the number per target matters for the report
    let mut usable = vec![true; loci.len()];
    let mut seen: HashMap<usize, usize> = HashMap::new();
    for (i, l) in loci.iter().enumerate() {
        if let Some(t) = l.control_of {
            let k = seen.entry(t).or_default();
            usable[i] = *k < n_ctl_json[t];
            *k += 1;
        }
    }
    let lidx: HashMap<String, usize> = targets.iter().enumerate().map(|(i, l)| (l.id.clone(), i)).collect();
    let mut pf = panel_frac.clone();
    pf.resize(loci.len(), (f64::NAN, f64::NAN));
    let mut af = alt_frac.clone();
    af.resize(loci.len(), 0.0);

    // ---- samples
    let mut infos: Vec<SampleInfo> = Vec::new();
    for x in d["samples"].as_array().ok_or_else(|| anyhow!("report has no samples"))? {
        infos.push(SampleInfo {
            id: x["id"].as_str().unwrap_or("").to_string(),
            group: x["group"].as_str().unwrap_or("-").to_string(),
            role: st(x["role"].as_str().unwrap_or("test")),
            sex_declared: x["sex"].as_str().map(st),
            sex_inferred: x["sex_inferred"].as_str().map(st),
            x_ratio: x["x_ratio"].as_f64(),
            // warnings are recomputed below with the current wording
            qc: Vec::new(),
            expected: x["expected"].as_array().map(|a| a.iter().filter_map(|v| v.as_u64().map(|v| v as usize)).filter(|&v| v < nt).collect()).unwrap_or_default(),
            expected_raw: x["expected_raw"].as_str().unwrap_or("").to_string(),
            note: x["note"].as_str().unwrap_or("").to_string(),
            hifi_bam: x["hifi_bam"].as_str().unwrap_or("").to_string(),
            fail_bams: strs(&x["fail_bams"]),
            hifi_raw_bams: strs(&x["hifi_raw_bams"]),
            locate_notes: strs(&x["notes"]),
            unassigned: strs(&x["unassigned"]),
            unassigned_used: strs(&x["unassigned_used"]),
        });
    }
    for info in infos.iter_mut() {
        if info.role == "positive" && info.expected.is_empty() {
            let g = info.group.clone();
            info.qc.push(report::qc_positive_without_locus(&g));
        }
        if let (Some(d), Some(i), Some(r)) = (info.sex_declared, info.sex_inferred, info.x_ratio) {
            if d != i {
                info.qc.push(report::qc_sex_mismatch(d, i, r));
            }
        }
    }
    let sidx: HashMap<String, usize> = infos.iter().enumerate().map(|(i, s)| (s.id.clone(), i)).collect();

    // ---- reads (target loci); states recomputed with the current thresholds
    let (h, rr) = tsv(&from.join("reads.tsv"))?;
    let c = |n: &str| col(&h, n).ok_or_else(|| anyhow!("reads.tsv: no column {n}"));
    let (ci_s, ci_l, ci_src, ci_bam, ci_name, ci_len) = (c("sample_id")?, c("locus")?, c("source")?, c("bam")?, c("read_name")?, c("read_len")?);
    let (ci_ff, ci_rq, ci_np, ci_hp, ci_strand, ci_cls) = (c("ff")?, c("rq")?, c("np")?, c("hp")?, c("strand")?, c("class")?);
    let (ci_bp, ci_bound, ci_dm, ci_df, ci_pur) = (c("allele_bp")?, c("bound")?, c("dominant_motif")?, c("motif_frac")?, c("purity")?);
    let (ci_hl, ci_hr, ci_fl) = (c("hits_left")?, c("hits_right")?, c("flags")?);
    let mut rows: Vec<ReadRow> = Vec::with_capacity(rr.len());
    for f in &rr {
        let (Some(&si), Some(&li)) = (sidx.get(&f[ci_s]), lidx.get(&f[ci_l])) else { continue };
        let l = &loci[li];
        let dm = f[ci_dm].clone();
        let mut call = Call {
            locus: li,
            strand: f[ci_strand].chars().next().unwrap_or('+'),
            cls: Cls::parse(&f[ci_cls]),
            allele_bp: f[ci_bp].parse().ok(),
            bound: Bound::parse(&f[ci_bound]),
            state: "",
            path_motif: dm != "-" && !dm.is_empty() && l.motif_is_pathogenic(&dm),
            dom_motif: dm,
            dom_frac: num(&f[ci_df]).unwrap_or(0.0),
            purity: num(&f[ci_pur]).unwrap_or(0.0),
            hits_l: f[ci_hl].parse().unwrap_or(0),
            hits_r: f[ci_hr].parse().unwrap_or(0),
            flags: f[ci_fl].split(',').filter_map(flag).collect(),
        };
        finalize_state(&mut call, l);
        rows.push(ReadRow {
            sample: si,
            source: f[ci_src].clone(),
            bam: f[ci_bam].clone(),
            name: f[ci_name].clone(),
            read_len: f[ci_len].parse().unwrap_or(0),
            ff: f[ci_ff].parse().ok(),
            rq: f[ci_rq].parse().ok(),
            np: f[ci_np].parse().ok(),
            hp: f[ci_hp].parse().ok(),
            call,
        });
    }
    eprintln!("[report] {} reads, {} samples, {} loci", rows.len(), infos.len(), nt);

    // ---- depth curves
    let (h, dr) = tsv(&from.join("depth.tsv"))?;
    let c = |n: &str| col(&h, n).ok_or_else(|| anyhow!("depth.tsv: no column {n}"));
    let (ds, dl, doff, dt, d1, d2, dg) = (c("sample_id")?, c("locus")?, c("offset_bp")?, c("rel_total")?, c("rel_hp1")?, c("rel_hp2")?, c("template_g")?);
    let mut curves: HashMap<(usize, usize), Canyon> = HashMap::new();
    for f in &dr {
        let (Some(&si), Some(&li)) = (sidx.get(&f[ds]), lidx.get(&f[dl])) else { continue };
        let cy = curves.entry((si, li)).or_default();
        cy.offsets.push(f[doff].parse().unwrap_or(0));
        cy.rel_total.push(num(&f[dt]).unwrap_or(f64::NAN));
        cy.rel_hp1.push(num(&f[d1]).unwrap_or(f64::NAN));
        cy.rel_hp2.push(num(&f[d2]).unwrap_or(f64::NAN));
        cy.g.push(num(&f[dg]).unwrap_or(0.0));
    }

    // ---- per sample x locus numbers
    let (h, sr) = tsv(&from.join("summary.tsv"))?;
    let c = |n: &str| col(&h, n).ok_or_else(|| anyhow!("summary.tsv: no column {n}"));
    let (ss, sl, serr, sx, sb, sc) = (c("sample_id")?, c("locus")?, c("error")?, c("chrX_auto_ratio")?, c("baseline_depth")?, c("centre_rel")?);
    let (sft, sf1, sf2, sst, ss1, ss2, sph, sn50, sctl) =
        (c("f_total")?, c("f_hp1")?, c("f_hp2")?, c("se_total")?, c("se_hp1")?, c("se_hp2")?, c("phased_frac")?, c("read_n50")?, c("usable_controls")?);
    let mut sums: Vec<LocusSummary> = Vec::new();
    let mut ctl_mean: HashMap<(usize, usize), BTreeMap<&'static str, f64>> = HashMap::new();
    let mut n_ctl: HashMap<usize, usize> = HashMap::new();
    for f in &sr {
        let (Some(&si), Some(&li)) = (sidx.get(&f[ss]), lidx.get(&f[sl])) else { continue };
        let err = f[serr].trim();
        let error = (!err.is_empty() && err != "-" && err != "NA").then(|| err.to_string());
        let canyon = curves.remove(&(si, li)).map(|mut cy| {
            cy.base_total = num(&f[sb]).unwrap_or(f64::NAN);
            cy.centre_total = num(&f[sc]).unwrap_or(f64::NAN) * cy.base_total;
            cy.base_hp1 = f64::NAN;
            cy.base_hp2 = f64::NAN;
            cy.f_total = num(&f[sft]);
            cy.f_hp1 = num(&f[sf1]);
            cy.f_hp2 = num(&f[sf2]);
            cy.se_total = num(&f[sst]);
            cy.se_hp1 = num(&f[ss1]);
            cy.se_hp2 = num(&f[ss2]);
            cy.read_n50 = f[sn50].parse().unwrap_or(0);
            cy
        });
        let info = &infos[si];
        let sex_note = targets[li].is_sex_chrom().then(|| report::sex_note(info.sex_declared, info.sex_inferred, info.x_ratio));
        sums.push(LocusSummary {
            sample: si,
            locus: li,
            canyon: if error.is_some() { None } else { canyon },
            phased_frac: num(&f[sph]).unwrap_or(0.0),
            error,
            x_ratio: num(&f[sx]),
            sex_note,
        });
        let mut m = BTreeMap::new();
        for src in SRC_GROUPS {
            if let Some(i) = col(&h, &format!("{src}.ctl_mean")) {
                if let Some(v) = num(&f[i]) {
                    m.insert(src, v);
                }
            }
        }
        ctl_mean.insert((si, li), m);
        n_ctl.insert(li, f[sctl].parse().unwrap_or(0));
    }
    sums.sort_by_key(|s| (s.sample, s.locus));

    // ---- BAM stats
    let mut bam_stats = Vec::new();
    if let Ok((h, br)) = tsv(&from.join("bam_stats.tsv")) {
        let g = |f: &Vec<String>, n: &str| col(&h, n).and_then(|i| f.get(i)).cloned().unwrap_or_default();
        for f in &br {
            let Some(&si) = sidx.get(&g(f, "sample_id")) else { continue };
            let mut by_cat = BTreeMap::new();
            for kv in g(f, "category_reads").split(',') {
                if let Some((k, v)) = kv.split_once(':') {
                    by_cat.insert(k.to_string(), CatStat { reads: v.parse().unwrap_or(0), bases: 0 });
                }
            }
            bam_stats.push(BamStat {
                sample: si,
                kind: g(f, "kind"),
                path: g(f, "bam"),
                stats: ScanStats {
                    reads: g(f, "reads").parse().unwrap_or(0),
                    bases: g(f, "bases").parse().unwrap_or(0),
                    prefilter: g(f, "prefilter_pass").parse().unwrap_or(0),
                    candidates: g(f, "placed").parse().unwrap_or(0),
                    seconds: g(f, "seconds").parse().unwrap_or(0.0),
                    by_cat,
                },
            });
        }
    }

    // ---- derived numbers and outputs
    report::low_depth_qc(&sums, &mut infos);
    let mut dv: Derived = report::derive_from(&rows, &sums, &loci, ctl_mean, n_ctl);
    report::add_cohort(&mut dv, &sums, &loci, &infos);
    fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    report::write_reads(&out.join("reads.tsv"), &rows, &infos, &loci)?;
    report::write_depth(&out.join("depth.tsv"), &sums, &infos, &loci)?;
    report::write_bam_stats(&out.join("bam_stats.tsv"), &bam_stats, &infos)?;
    report::write_summary(&out.join("summary.tsv"), &sums, &infos, &loci, &dv, lang)?;
    let mut params: BTreeMap<&str, String> = BTreeMap::new();
    if let Some(p) = meta["params"].as_object() {
        for (k, v) in p {
            let k: &'static str = Box::leak(k.clone().into_boxed_str());
            params.insert(k, v.as_str().map(String::from).unwrap_or_else(|| v.to_string()));
        }
    }
    params.insert("rebuilt_from", format!("{} (trace v{version_from}, {})", from.display(), meta["created"].as_str().unwrap_or("")));
    report::write_html(
        &out.join("report.html"),
        &sums,
        &rows,
        &infos,
        &loci,
        &usable,
        &pf,
        &af,
        &bam_stats,
        &dv,
        version,
        meta["reference"].as_str().unwrap_or(""),
        params,
        lang,
    )?;
    let run = serde_json::json!({
        "tool": "tr-rescue", "version": version, "command": "report",
        "created": crate::now_string(), "from": from.display().to_string(), "trace_version": version_from,
    });
    fs::write(out.join("run.json"), serde_json::to_string_pretty(&run)?)?;
    eprintln!("[report] done → {}", out.join("report.html").display());
    Ok(())
}
