// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! TR-Rescue: dropout-aware evidence for tandem-repeat expansions in PacBio HiFi WGS.
//!
//! `trace` (step 0): for every catalog locus and sample, measure HiFi alleles, the
//! read-length-wide depth canyon, and where repeat-containing molecules went (HiFi,
//! fail 0x01 / 0x08, or nowhere), using unique flank k-mers so unaligned fail reads count.

#![allow(dead_code)] // fields kept for the M1–M3 modules that build on step 0

mod alleles;
mod anchor;
mod catalog;
mod depth;
mod fasta;
mod hifi;
mod i18n;
mod kmer;
mod manifest;
mod report;
mod rereport;
mod scan;

use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::{Args, Parser, Subcommand};

use crate::anchor::{MotifSet, analyze};
use crate::catalog::{Locus, add_controls, builtin, find, parse_bed};
use crate::fasta::Fasta;
use crate::hifi::AlnBam;
use crate::kmer::Panel;
use crate::manifest::{index_run_root, locate, read_manifest, rg_fields};
use crate::report::{BamStat, LocusSummary, ReadRow, SampleInfo};
use crate::scan::{ScanCtx, ff_category, scan_bam};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[command(
    name = "tr-rescue",
    version,
    about = "Dropout-aware tandem repeat expansion evidence for PacBio HiFi WGS",
    after_help = "Copyright (c) 2026 Seung Won Chae · MIT License · https://github.com/CleanCham-lab/TR-Rescue\nLocus catalog: STRchive (CC BY 4.0, Hiatt et al. Genome Medicine 2025). Research use only — not a diagnosis."
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Step 0: alleles, depth canyon and read whereabouts at every catalog locus
    Trace(TraceArgs),
    /// Show which fail/hifi BAMs match an aligned BAM's @RG (run before `trace`)
    Locate {
        /// Revio run root(s), e.g. /data/revio/revio_run (several allowed)
        #[arg(long, required = true, num_args = 1..)]
        run_root: Vec<PathBuf>,
        /// Aligned BAM(s)
        #[arg(long, required = true, num_args = 1..)]
        bam: Vec<PathBuf>,
    },
    /// Check a manifest (.xlsx/.csv/.tsv) before `trace`: how each row is read, files, loci
    Check {
        /// Sample manifest (.xlsx / .csv / .tsv)
        #[arg(long)]
        manifest: PathBuf,
        /// Custom catalog BED (default: built-in STRchive)
        #[arg(long)]
        catalog: Option<PathBuf>,
        /// Revio run root(s), to resolve fail_bams=auto
        #[arg(long, num_args = 1..)]
        run_root: Vec<PathBuf>,
        /// Output language: ko or en
        #[arg(long, default_value = "ko", value_parser = ["ko", "en"])]
        lang: String,
    },
    /// Rebuild summary / reads / report.html from an earlier `trace` folder with this version's rules (no BAM access)
    Report {
        /// Existing trace output folder (report.html, summary.tsv, reads.tsv, depth.tsv)
        #[arg(long)]
        from: PathBuf,
        /// New output folder (must differ from --from)
        #[arg(long)]
        out: PathBuf,
        /// Custom catalog BED used in that run, if any
        #[arg(long)]
        catalog: Option<PathBuf>,
        /// Reference FASTA (optional; only for the reference-region QC column)
        #[arg(long)]
        reference: Option<PathBuf>,
        /// Language of summary.tsv hints and the report's starting language: ko or en
        #[arg(long, default_value = "ko", value_parser = ["ko", "en"])]
        lang: String,
    },
    /// Print the built-in locus catalog (STRchive snapshot, hg38)
    Catalog,
}

#[derive(Args)]
struct TraceArgs {
    /// Reference FASTA used for alignment (needs .fai)
    #[arg(long)]
    reference: PathBuf,
    /// Sample manifest: .xlsx / .csv / .tsv (sample_id, group, sex, locus, hifi_bam, fail_bams, hifi_raw_bams, note)
    #[arg(long)]
    manifest: PathBuf,
    /// Output directory
    #[arg(long)]
    out: PathBuf,
    /// Restrict to these loci (comma-separated STRchive IDs or gene symbols). Default: all
    #[arg(long)]
    loci: Option<String>,
    /// Extra TRGT-style catalog BED (ID=..;MOTIFS=..; optional NORMAL_MAX, PATH_MIN in units)
    #[arg(long)]
    catalog: Option<PathBuf>,
    /// Revio run root(s) for `auto` fail_bams / hifi_raw_bams (several allowed)
    #[arg(long, num_args = 1..)]
    run_root: Vec<PathBuf>,
    /// Threads
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// k-mer size for flank anchoring (odd, <= 31)
    #[arg(short = 'k', long, default_value_t = 21)]
    k: usize,
    /// Minimum consistent flank k-mer hits on a side to place a read
    #[arg(long, default_value_t = 10)]
    min_hits: usize,
    /// Flank length on each side used for anchoring (bp)
    #[arg(long, default_value_t = 3000)]
    flank: u64,
    /// Gap between repeat edge and flank (bp)
    #[arg(long, default_value_t = 30)]
    pad: u64,
    /// Half-width of the depth profile window (bp)
    #[arg(long, default_value_t = 100_000)]
    window: u64,
    /// Depth bin size (bp)
    #[arg(long, default_value_t = 500)]
    bin: u64,
    /// |offset| beyond which depth bins form the baseline (bp)
    #[arg(long, default_value_t = 30_000)]
    distal: u64,
    /// Minimum MAPQ for depth (0 keeps reads split with ALT contigs)
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,
    /// Language of summary.tsv hints and the report's starting language (the report can switch): ko or en
    #[arg(long, default_value = "ko", value_parser = ["ko", "en"])]
    lang: String,
    /// Only these samples from the manifest (comma-separated sample_id). Default: all
    #[arg(long)]
    samples: Option<String>,
    /// Unassigned BAMs (<movie>.fail_reads/hifi_reads.unassigned.bam) with `auto`: single-plex cells only (default),
    /// never, or always (also multiplexed cells; reads of other samples get mixed in)
    #[arg(long, default_value = "auto", value_parser = ["auto", "never", "always"])]
    unassigned: String,
    /// Reuse / save the flank k-mer panel here (~100 MB for the full catalog; default: rebuild each run)
    #[arg(long)]
    panel_cache: Option<PathBuf>,
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(16)
}

static START: std::sync::OnceLock<Instant> = std::sync::OnceLock::new();

macro_rules! log {
    ($($arg:tt)*) => {{
        let t = START.get_or_init(Instant::now).elapsed().as_secs_f64();
        eprintln!("[{:>7.1}s] {}", t, format!($($arg)*));
    }};
}

/// UTC "YYYY-MM-DD HH:MM" without external crates.
pub fn now_string() -> String {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0) as i64;
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    // civil from days (Howard Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!("{y:04}-{m:02}-{d:02} {:02}:{:02} UTC", rem / 3600, (rem % 3600) / 60)
}

fn main() {
    START.get_or_init(Instant::now);
    let cli = Cli::parse();
    let r = match cli.cmd {
        Cmd::Trace(a) => trace(a),
        Cmd::Locate { run_root, bam } => cmd_locate(&run_root, &bam),
        Cmd::Check { manifest, catalog, run_root, lang } => cmd_check(&manifest, catalog.as_deref(), &run_root, &lang),
        Cmd::Report { from, out, catalog, reference, lang } => rereport::run(&from, &out, catalog.as_deref(), reference.as_deref(), VERSION, &lang),
        Cmd::Catalog => {
            let mut out = std::io::stdout().lock();
            let _ = writeln!(out, "#id\tgene\tregion_hg38\tmotifs\tpathogenic_motifs\tnormal_max\tpathogenic_min\tinheritance\tevidence\tdisease");
            for l in builtin() {
                if writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    l.id,
                    l.gene,
                    l.region_str(),
                    l.motifs.join(","),
                    l.path_motifs.join(","),
                    l.normal_max_units,
                    l.path_min_units,
                    l.inheritance,
                    l.evidence,
                    l.disease
                )
                .is_err()
                {
                    break;
                }
            }
            Ok(())
        }
    };
    if let Err(e) = r {
        eprintln!("error: {e:#}");
        std::process::exit(1);
    }
}

fn index_roots(roots: &[PathBuf], kind: &str) -> BTreeMap<String, Vec<PathBuf>> {
    let mut all: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    for r in roots {
        if !r.is_dir() {
            log!("WARNING run root not found: {}", r.display());
            continue;
        }
        for (k, v) in index_run_root(r, kind) {
            all.entry(k).or_default().extend(v);
        }
    }
    all
}

fn cmd_locate(roots: &[PathBuf], bams: &[PathBuf]) -> Result<()> {
    log!("indexing {}", roots.iter().map(|r| r.display().to_string()).collect::<Vec<_>>().join(", "));
    let fail = index_roots(roots, "fail_reads");
    let hifi = index_roots(roots, "hifi_reads");
    log!("{} movies with fail_reads, {} with hifi_reads", fail.len(), hifi.len());
    for b in bams {
        let rgs = rg_fields(b)?;
        println!("# {}", b.display());
        for (pu, sm, bc) in &rgs {
            println!("#   RG PU={pu} SM={sm} BC={bc}");
        }
        let lf = locate(&rgs, &fail);
        let lh = locate(&rgs, &hifi);
        for p in &lf.matched {
            println!("fail_reads\t{}", p.display());
        }
        for p in &lh.matched {
            println!("hifi_reads\t{}", p.display());
        }
        let lhu = locate(&rgs, &hifi);
        for p in lf.unassigned_safe.iter().chain(lhu.unassigned_safe.iter()) {
            println!("unassigned(single-plex, used)\t{}", p.display());
        }
        for p in lf.unassigned.iter().chain(lhu.unassigned.iter()) {
            println!("unassigned(multiplexed, not used)\t{}", p.display());
        }
        for n in lf.notes.iter().chain(lh.notes.iter()) {
            println!("note\t{n}");
        }
    }
    Ok(())
}

fn is_unassigned(p: &Path) -> bool {
    p.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.contains(".unassigned."))
}

fn bam_index(p: &Path) -> Option<PathBuf> {
    let s = p.display().to_string();
    [format!("{s}.bai"), format!("{s}.csi"), s.trim_end_matches(".bam").to_string() + ".bai"]
        .into_iter()
        .map(PathBuf::from)
        .find(|x| x.is_file())
}

fn cmd_check(manifest: &Path, catalog: Option<&Path>, roots: &[PathBuf], lang: &str) -> Result<()> {
    let t = |ko: &'static str, en: &'static str| i18n::t(lang, ko, en);
    let samples = read_manifest(manifest)?;
    let cat = match catalog {
        Some(p) => parse_bed(p)?,
        None => builtin(),
    };
    let (fail_idx, hifi_idx) = if roots.is_empty() {
        (BTreeMap::new(), BTreeMap::new())
    } else {
        (index_roots(roots, "fail_reads"), index_roots(roots, "hifi_reads"))
    };
    let mut errors = 0usize;
    let mut warns = 0usize;
    println!("{} {} ({})", samples.len(), t("샘플", "samples"), manifest.display());
    for s in &samples {
        println!("\n# {}", s.id);
        println!(
            "  group  {:<14} → {} {}",
            s.group,
            t("역할", "role"),
            match s.role {
                "positive" => t("양성", "positive"),
                "control" => t("대조", "control"),
                _ => t("검사", "test"),
            }
        );
        println!(
            "  sex    {:<14} → {}",
            if s.sex_raw.is_empty() { "-" } else { &s.sex_raw },
            match s.sex {
                Some("M") => t("남", "male"),
                Some("F") => t("여", "female"),
                _ => t("미기재 → chrX depth로 추정", "not given → inferred from chrX depth"),
            }
        );
        let keys: Vec<&str> = s.locus.split(',').map(|x| x.trim()).filter(|x| !x.is_empty() && *x != "-" && !x.eq_ignore_ascii_case("NA")).collect();
        if keys.is_empty() {
            println!("  locus  -");
            if s.role == "positive" {
                println!("  WARN   {}", t("양성인데 기존 진단 좌위(locus)가 비어 있음", "positive, but the locus column (known diagnosis) is empty"));
                warns += 1;
            }
        }
        for k in keys {
            match find(&cat, k) {
                Ok(l) => println!("  locus  {k:<14} → {} ({} {})", l.id, l.chrom, l.start),
                Err(e) => {
                    println!("  ERROR  locus {k}: {e:#}");
                    errors += 1;
                }
            }
        }
        let ok = s.hifi_bam.is_file();
        println!("  bam    {} {}", if ok { "OK " } else { t("없음", "MISSING") }, s.hifi_bam.display());
        if !ok {
            errors += 1;
        } else if bam_index(&s.hifi_bam).is_none() {
            println!("  ERROR  {}", t("BAM index(.bai/.csi) 없음", "no BAM index (.bai/.csi)"));
            errors += 1;
        }
        let mut list = |name: &str, paths: &[PathBuf], auto: bool, idx: &BTreeMap<String, Vec<PathBuf>>| {
            if auto {
                if roots.is_empty() {
                    println!("  {name:<6} auto ({})", t("trace에 --run-root 필요", "needs --run-root in trace"));
                } else if ok {
                    match rg_fields(&s.hifi_bam) {
                        Ok(rgs) => {
                            let l = locate(&rgs, idx);
                            if l.matched.is_empty() {
                                println!("  WARN   {name} auto: {} ({})", t("맞는 파일 없음", "no matching file"), l.notes.join("; "));
                                warns += 1;
                            }
                            for p in &l.matched {
                                println!("  {name:<6} auto → {}", p.display());
                            }
                            if name == "fail" {
                                for p in &l.unassigned_safe {
                                    println!("  unasg  single-plex → {} {}", t("사용", "used"), p.display());
                                }
                                for p in &l.unassigned {
                                    println!("  unasg  {} {}", t("여러 barcode 셀 → 사용 안 함", "multiplexed cell → not used"), p.display());
                                }
                            }
                        }
                        Err(e) => {
                            println!("  WARN   {name} auto: {} {e:#}", t("@RG 읽기 실패", "cannot read @RG"));
                            warns += 1;
                        }
                    }
                }
            } else if paths.is_empty() {
                println!("  {name:<6} -");
            }
            for p in paths {
                let ok = p.is_file();
                println!("  {name:<6} {} {}", if ok { "OK " } else { t("없음", "MISSING") }, p.display());
                if !ok {
                    errors += 1;
                }
            }
        };
        list("fail", &s.fail_bams, s.fail_auto, &fail_idx);
        if s.fail_auto && ok && !roots.is_empty() {
            if let Ok(rgs) = rg_fields(&s.hifi_bam) {
                let l = locate(&rgs, &hifi_idx);
                for p in &l.unassigned_safe {
                    println!("  unasg  single-plex → {} {}", t("사용", "used"), p.display());
                }
                for p in &l.unassigned {
                    println!("  unasg  {} {}", t("여러 barcode 셀 → 사용 안 함", "multiplexed cell → not used"), p.display());
                }
            }
        }
        list("raw", &s.hifi_raw_bams, s.hifi_raw_auto, &hifi_idx);
        if !s.note.is_empty() {
            println!("  note   {}", s.note);
        }
    }
    println!("\n{} {errors} · {} {warns}", t("오류", "errors"), t("경고", "warnings"));
    if errors > 0 {
        bail!("{errors} problem(s) in the manifest");
    }
    Ok(())
}

/// Results of the aligned-BAM pass for one sample.
struct AlnOut {
    sums: Vec<LocusSummary>,
    rows: Vec<ReadRow>,
    names: HashSet<String>,
}

/// Process a set of loci on one aligned BAM (one reader per thread).
fn aligned_pass(bam_path: &Path, si: usize, loci: &[Locus], idx: &[usize], ctx: &ScanCtx<'_>, a: &TraceArgs) -> Result<AlnOut> {
    let mut bam = AlnBam::open(bam_path)?;
    let mut out = AlnOut { sums: Vec::new(), rows: Vec::new(), names: HashSet::new() };
    let bam_name = file_name(bam_path);
    for &li in idx {
        let l = &loci[li];
        let Some(chrom) = bam.resolve(&l.chrom) else {
            if !l.is_control() {
                out.sums.push(LocusSummary {
                    sample: si,
                    locus: li,
                    canyon: None,
                    phased_frac: 0.0,
                    error: Some(format!("contig {} not in BAM header", l.chrom)),
                    x_ratio: None, sex_note: None,
                });
            }
            continue;
        };
        if !l.is_control() {
            let c = l.centre();
            let ws = c.saturating_sub(a.window);
            match bam.depth(&chrom, ws, c + a.window, a.bin, a.min_mapq) {
                Ok(d) => {
                    let phased = if d.n_reads > 0 { d.n_phased as f64 / d.n_reads as f64 } else { 0.0 };
                    let cy = depth::canyon(&d, l.start, l.end, a.distal as i64, 1000);
                    out.sums.push(LocusSummary { sample: si, locus: li, canyon: Some(cy), phased_frac: phased, error: None, x_ratio: None, sex_note: None });
                }
                Err(e) => out.sums.push(LocusSummary {
                    sample: si,
                    locus: li,
                    canyon: None,
                    phased_frac: 0.0,
                    error: Some(format!("{e:#}")),
                    x_ratio: None, sex_note: None,
                }),
            }
        }
        let rs = l.start.saturating_sub(a.pad + a.flank);
        let re = l.end + a.pad + a.flank;
        for r in bam.reads(&chrom, rs, re)? {
            out.names.insert(r.name.clone());
            for call in analyze(&r.seq, ctx.panel, ctx.loci, ctx.motif_sets, ctx.min_hits) {
                if call.locus != li {
                    continue;
                }
                out.rows.push(ReadRow {
                    sample: si,
                    source: "hifi_aligned".into(),
                    bam: bam_name.clone(),
                    name: r.name.clone(),
                    read_len: r.seq.len(),
                    ff: None,
                    rq: None,
                    np: None,
                    hp: r.hp,
                    call,
                });
            }
        }
    }
    Ok(out)
}

fn trace(a: TraceArgs) -> Result<()> {
    if a.k % 2 == 0 || a.k < 15 || a.k > 31 {
        bail!("-k must be odd and in 15..=31");
    }
    fs::create_dir_all(&a.out).with_context(|| format!("creating {}", a.out.display()))?;
    let mut samples = read_manifest(&a.manifest)?;
    if let Some(sel) = &a.samples {
        let want: Vec<&str> = sel.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()).collect();
        for w in &want {
            if !samples.iter().any(|s| s.id == *w) {
                bail!("--samples: {w} is not in the manifest");
            }
        }
        samples.retain(|s| want.contains(&s.id.as_str()));
        log!("--samples: {} of the manifest", samples.iter().map(|s| s.id.as_str()).collect::<Vec<_>>().join(", "));
    }
    let fa = Fasta::open(&a.reference)?;

    // ---- loci
    let mut cat = builtin();
    if let Some(p) = &a.catalog {
        let extra = parse_bed(p)?;
        cat.retain(|b| !extra.iter().any(|e| e.id.eq_ignore_ascii_case(&b.id)));
        cat.extend(extra);
    }
    let chosen: Vec<Locus> = match &a.loci {
        Some(s) => {
            let mut v = Vec::new();
            for key in s.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()) {
                v.push(find(&cat, key)?);
            }
            v
        }
        None => cat.clone(),
    };
    let mut targets = Vec::new();
    for mut l in chosen {
        match fa.resolve(&l.chrom) {
            Some(c) => {
                l.chrom = c;
                if let Ok(seq) = fa.fetch(&l.chrom, l.start, l.end) {
                    l.ref_edge_bp = catalog::edge_nonrepeat(&seq, &l.motifs, l.unit_len());
                }
                targets.push(l);
            }
            None => log!("skipping {} (contig {} not in reference)", l.id, l.chrom),
        }
    }
    for l in &targets {
        if !l.count_note.is_empty() {
            log!("{}: {}", l.id, l.count_note.en);
        } else if let Some((le, ri)) = l.ref_edge_bp {
            if le + ri >= 2 * l.unit_len() as i64 {
                log!(
                    "QC {}: reference region has {} bp (left {le}, right {ri}) not matching its motifs; repeat counts may be over by ~{:.0}",
                    l.id,
                    le + ri,
                    (le + ri) as f64 / l.unit_len() as f64
                );
            }
        }
    }
    if targets.is_empty() {
        bail!("no target loci");
    }
    let loci = add_controls(&targets);
    log!("{} target loci (+{} neighbourhood controls)", targets.len(), loci.len() - targets.len());

    // expected (confirmed) loci per sample, from the manifest
    let mut infos: Vec<SampleInfo> = Vec::new();
    for s in &samples {
        let mut exp = Vec::new();
        for key in s.locus.split(',').map(|x| x.trim()).filter(|x| !x.is_empty() && *x != "-" && !x.eq_ignore_ascii_case("NA")) {
            match find(&targets, key) {
                Ok(l) => {
                    if let Some(i) = targets.iter().position(|t| t.id == l.id) {
                        exp.push(i);
                    }
                }
                Err(e) => log!("WARNING {}: {e:#}", s.id),
            }
        }
        let mut qc = Vec::new();
        if s.role == "positive" && exp.is_empty() {
            qc.push(report::qc_positive_without_locus(&s.group));
            log!("WARNING {}: positive sample without a usable locus", s.id);
        }
        infos.push(SampleInfo {
            id: s.id.clone(),
            group: s.group.clone(),
            role: s.role,
            sex_declared: s.sex,
            sex_inferred: None,
            x_ratio: None,
            qc,
            expected: exp,
            expected_raw: s.locus.clone(),
            note: s.note.clone(),
            hifi_bam: s.hifi_bam.display().to_string(),
            fail_bams: vec![],
            hifi_raw_bams: vec![],
            locate_notes: vec![],
            unassigned: vec![],
            unassigned_used: vec![],
        });
    }

    // ---- panel
    let ref_size = fs::metadata(&a.reference).map(|m| m.len()).unwrap_or(0);
    let meta = format!(
        "v3;k={};flank={};pad={};ref={}:{};loci={}",
        a.k,
        a.flank,
        a.pad,
        a.reference.display(),
        ref_size,
        loci.iter().map(|l| format!("{}:{}:{}:{}", l.id, l.chrom, l.start, l.end)).collect::<Vec<_>>().join(",")
    );
    let cached = a.panel_cache.as_ref().and_then(|c| Panel::load(c, &meta, a.k, a.flank, a.pad));
    let panel = match cached {
        Some(p) => {
            log!("panel loaded ({} k-mers)", p.len());
            p
        }
        None => {
            log!("building flank k-mer panel (k={}, flank={} bp) and counting in reference ...", a.k, a.flank);
            let p = Panel::build(&fa, &loci, a.k, a.flank, a.pad, a.threads)?;
            log!("panel: {} unique k-mers", p.len());
            if let Some(c) = &a.panel_cache {
                p.save(c, &meta)?;
                log!("panel saved to {}", c.display());
            }
            p
        }
    };
    let mut panel_frac = vec![(0.0, 0.0); loci.len()];
    let mut alt_frac = vec![0.0f64; loci.len()];
    for st in &panel.stats {
        let f = if st.kmers > 0 { st.unique as f64 / st.kmers as f64 } else { 0.0 };
        if st.side == 0 {
            panel_frac[st.locus].0 = f;
        } else {
            panel_frac[st.locus].1 = f;
        }
        if st.unique > 0 {
            alt_frac[st.locus] = alt_frac[st.locus].max(st.alt_shared as f64 / st.unique as f64);
        }
    }
    let usable: Vec<bool> = loci
        .iter()
        .enumerate()
        .map(|(i, l)| !l.is_control() || (panel_frac[i].0 >= 0.4 && panel_frac[i].1 >= 0.4))
        .collect();
    let low: Vec<&str> = loci
        .iter()
        .enumerate()
        .filter(|(i, l)| !l.is_control() && (panel_frac[*i].0 < 0.2 || panel_frac[*i].1 < 0.2))
        .map(|(_, l)| l.id.as_str())
        .collect();
    if !low.is_empty() {
        log!("few unique flank k-mers (<20% on a side): {}", low.join(", "));
    }
    let motif_sets: Vec<MotifSet> = loci.iter().map(|l| MotifSet::new(&l.motifs)).collect();
    let ctx = ScanCtx { panel: &panel, loci: &loci, motif_sets: &motif_sets, min_hits: a.min_hits };

    // ---- run root index
    let need_root = samples.iter().any(|s| s.fail_auto || s.hifi_raw_auto);
    let (fail_idx, hifi_idx) = if need_root {
        if a.run_root.is_empty() {
            bail!("manifest uses `auto` but --run-root was not given");
        }
        log!("indexing run root(s) {}", a.run_root.iter().map(|r| r.display().to_string()).collect::<Vec<_>>().join(", "));
        let f = index_roots(&a.run_root, "fail_reads");
        let h = if samples.iter().any(|s| s.hifi_raw_auto) || (a.unassigned != "never" && samples.iter().any(|s| s.fail_auto)) {
            index_roots(&a.run_root, "hifi_reads")
        } else {
            BTreeMap::new()
        };
        log!("found fail_reads for {} movies", f.len());
        (f, h)
    } else {
        (BTreeMap::new(), BTreeMap::new())
    };

    let mut rows: Vec<ReadRow> = Vec::new();
    let mut sums: Vec<LocusSummary> = Vec::new();
    let mut bam_stats: Vec<BamStat> = Vec::new();

    for (si, s) in samples.iter().enumerate() {
        log!("== sample {} ({}) ==", s.id, s.group);
        let info = &mut infos[si];
        let mut fail_bams = s.fail_bams.clone();
        let mut raw_bams = s.hifi_raw_bams.clone();
        let mut un_bams: Vec<PathBuf> = Vec::new();
        if s.fail_auto || s.hifi_raw_auto {
            let rgs = rg_fields(&s.hifi_bam).with_context(|| format!("reading @RG of {}", s.hifi_bam.display()))?;
            let mut not_used: Vec<PathBuf> = Vec::new();
            let take = |l: &manifest::Located, un_bams: &mut Vec<PathBuf>, not_used: &mut Vec<PathBuf>| match a.unassigned.as_str() {
                "never" => {
                    not_used.extend(l.unassigned_safe.iter().cloned());
                    not_used.extend(l.unassigned.iter().cloned());
                }
                "always" => {
                    un_bams.extend(l.unassigned_safe.iter().cloned());
                    un_bams.extend(l.unassigned.iter().cloned());
                }
                _ => {
                    un_bams.extend(l.unassigned_safe.iter().cloned());
                    not_used.extend(l.unassigned.iter().cloned());
                }
            };
            if s.fail_auto {
                let l = locate(&rgs, &fail_idx);
                fail_bams.extend(l.matched.iter().cloned());
                info.locate_notes.extend(l.notes.iter().cloned());
                take(&l, &mut un_bams, &mut not_used);
                if a.unassigned != "never" {
                    let lh = locate(&rgs, &hifi_idx);
                    take(&lh, &mut un_bams, &mut not_used);
                }
            }
            if s.hifi_raw_auto {
                raw_bams.extend(locate(&rgs, &hifi_idx).matched);
            }
            info.unassigned = not_used.iter().map(|p| p.display().to_string()).collect();
            if !un_bams.is_empty() {
                log!("unassigned BAMs ({}): {}", a.unassigned, un_bams.iter().map(|p| file_name(p)).collect::<Vec<_>>().join(", "));
            }
        }
        // explicitly listed unassigned files are treated as unassigned too
        for v in [&mut fail_bams, &mut raw_bams] {
            let (u, k): (Vec<PathBuf>, Vec<PathBuf>) = v.drain(..).partition(|p| is_unassigned(p));
            un_bams.extend(u);
            *v = k;
        }
        un_bams.sort();
        un_bams.dedup();
        info.unassigned_used = un_bams.iter().map(|p| p.display().to_string()).collect();
        info.fail_bams = fail_bams.iter().map(|p| p.display().to_string()).collect();
        info.hifi_raw_bams = raw_bams.iter().map(|p| p.display().to_string()).collect();
        log!("fail BAMs: {}", if info.fail_bams.is_empty() { "none".to_string() } else { info.fail_bams.join(", ") });

        // ---- aligned HiFi, loci split across threads
        let t0 = Instant::now();
        let nthr = a.threads.max(1).min(loci.len());
        let chunks: Vec<Vec<usize>> = (0..nthr).map(|t| (t..loci.len()).step_by(nthr).collect()).collect();
        let results: Vec<Result<AlnOut>> = std::thread::scope(|sc| {
            let hs: Vec<_> = chunks.iter().map(|c| sc.spawn(|| aligned_pass(&s.hifi_bam, si, &loci, c, &ctx, &a))).collect();
            hs.into_iter().map(|h| h.join().unwrap_or_else(|_| Err(anyhow::anyhow!("thread panicked")))).collect()
        });
        let mut aligned_names: HashSet<String> = HashSet::new();
        let mut sample_sums = Vec::new();
        let mut failed = None;
        for r in results {
            match r {
                Ok(o) => {
                    sample_sums.extend(o.sums);
                    rows.extend(o.rows);
                    aligned_names.extend(o.names);
                }
                Err(e) => failed = Some(format!("{e:#}")),
            }
        }
        if let Some(e) = failed {
            log!("ERROR aligned BAM {}: {e}", s.hifi_bam.display());
            for (li, _) in loci.iter().enumerate().filter(|(_, l)| !l.is_control()) {
                if !sample_sums.iter().any(|x: &LocusSummary| x.locus == li) {
                    sample_sums.push(LocusSummary { sample: si, locus: li, canyon: None, phased_frac: 0.0, error: Some(e.clone()), x_ratio: None, sex_note: None });
                }
            }
        }
        // X copy number: chrX baseline / median autosomal baseline.
        // Declared sex (manifest) wins; depth is used when it is missing and to flag mismatches.
        let mut auto: Vec<f64> = sample_sums
            .iter()
            .filter(|x| !loci[x.locus].is_sex_chrom())
            .filter_map(|x| x.canyon.as_ref().map(|c| c.base_total))
            .filter(|b| b.is_finite() && *b > 0.0)
            .collect();
        auto.sort_by(|p, q| p.partial_cmp(q).unwrap());
        let mut xr: Vec<f64> = Vec::new();
        if !auto.is_empty() {
            let med = auto[auto.len() / 2];
            for x in sample_sums.iter_mut().filter(|x| loci[x.locus].chrom.trim_start_matches("chr") == "X") {
                if let Some(c) = x.canyon.as_ref() {
                    let r = c.base_total / med;
                    if r.is_finite() {
                        x.x_ratio = Some(r);
                        xr.push(r);
                    }
                }
            }
        }
        xr.sort_by(|p, q| p.partial_cmp(q).unwrap());
        let ratio = (!xr.is_empty()).then(|| xr[xr.len() / 2]);
        let inferred = ratio.and_then(|r| if r < 0.7 { Some("M") } else if r >= 0.8 { Some("F") } else { None });
        {
            info.x_ratio = ratio;
            info.sex_inferred = inferred;
            if let (Some(d), Some(i), Some(r)) = (info.sex_declared, inferred, ratio) {
                if d != i {
                    let msg = report::qc_sex_mismatch(d, i, r);
                    log!("WARNING {}: {}", info.id, msg.en);
                    info.qc.push(msg);
                }
            }
        }
        let declared = info.sex_declared;
        // one X if declared male, or (no declaration and) depth says one X; also when depth clearly says one X
        let one_x = declared == Some("M") || inferred == Some("M");
        for x in sample_sums.iter_mut().filter(|x| loci[x.locus].is_sex_chrom()) {
            x.sex_note = Some(report::sex_note(declared, inferred, ratio));
            if one_x {
                if let Some(c) = x.canyon.as_mut() {
                    c.f_hp1 = None;
                    c.f_hp2 = None;
                    c.se_hp1 = None;
                    c.se_hp2 = None;
                    c.rel_hp1.iter_mut().for_each(|v| *v = f64::NAN);
                    c.rel_hp2.iter_mut().for_each(|v| *v = f64::NAN);
                }
            }
        }
        log!("aligned HiFi: {} loci in {:.0}s", sample_sums.len(), t0.elapsed().as_secs_f64());
        sums.extend(sample_sums);

        // ---- unaligned BAMs
        let mut jobs: Vec<(PathBuf, bool)> = fail_bams.iter().map(|p| (p.clone(), true)).collect();
        jobs.extend(raw_bams.iter().map(|p| (p.clone(), false)));
        jobs.extend(un_bams.iter().map(|p| (p.clone(), p.display().to_string().contains(".fail_reads."))));
        for (p, is_fail) in jobs {
            let unassigned = is_unassigned(&p);
            log!("scanning {} ...", p.display());
            match scan_bam(&p, &ctx, a.threads, is_fail) {
                Err(e) => {
                    log!("ERROR scanning {}: {e:#}", p.display());
                    info.locate_notes.push(format!("{}: {e:#}", file_name(&p)));
                }
                Ok((st, cands)) => {
                    log!(
                        "  {} reads, {:.1} Gb, prefilter {}, placed {}, {:.0}s",
                        st.reads,
                        st.bases as f64 / 1e9,
                        st.prefilter,
                        cands.len(),
                        st.seconds
                    );
                    for c in cands {
                        let source = if aligned_names.contains(&c.name) && !is_fail {
                            continue;
                        } else if unassigned {
                            if is_fail { "unassigned_fail".to_string() } else { "unassigned_hifi".to_string() }
                        } else if is_fail {
                            ff_category(c.ff).to_string()
                        } else {
                            "hifi_not_in_aln".to_string()
                        };
                        for call in c.calls {
                            rows.push(ReadRow {
                                sample: si,
                                source: source.clone(),
                                bam: file_name(&p),
                                name: c.name.clone(),
                                read_len: c.len,
                                ff: c.ff,
                                rq: c.rq,
                                np: c.np,
                                hp: None,
                                call,
                            });
                        }
                    }
                    bam_stats.push(BamStat {
                        sample: si,
                        kind: match (is_fail, unassigned) {
                            (true, false) => "fail_reads".into(),
                            (false, false) => "hifi_reads".into(),
                            (true, true) => "fail_reads (unassigned)".into(),
                            (false, true) => "hifi_reads (unassigned)".into(),
                        },
                        path: p.display().to_string(),
                        stats: st,
                    });
                }
            }
        }
    }

    // ---- outputs
    sums.sort_by_key(|s| (s.sample, s.locus));
    report::low_depth_qc(&sums, &mut infos);
    let mut d = report::derive(&rows, &sums, &loci, &usable);
    report::add_cohort(&mut d, &sums, &loci, &infos);
    report::write_reads(&a.out.join("reads.tsv"), &rows, &infos, &loci)?;
    report::write_depth(&a.out.join("depth.tsv"), &sums, &infos, &loci)?;
    report::write_bam_stats(&a.out.join("bam_stats.tsv"), &bam_stats, &infos)?;
    report::write_summary(&a.out.join("summary.tsv"), &sums, &infos, &loci, &d, &a.lang)?;
    let mut params = BTreeMap::new();
    params.insert("k", a.k.to_string());
    params.insert("min_hits", a.min_hits.to_string());
    params.insert("flank_bp", a.flank.to_string());
    params.insert("pad_bp", a.pad.to_string());
    params.insert("window_bp", format!("±{}", a.window));
    params.insert("bin_bp", a.bin.to_string());
    params.insert("distal_bp", a.distal.to_string());
    params.insert("min_mapq", a.min_mapq.to_string());
    params.insert("min_purity", format!("{}", anchor::MIN_PURITY));
    params.insert("threads", a.threads.to_string());
    params.insert("manifest", a.manifest.display().to_string());
    params.insert("run_root", a.run_root.iter().map(|r| r.display().to_string()).collect::<Vec<_>>().join(", "));
    report::write_html(
        &a.out.join("report.html"),
        &sums,
        &rows,
        &infos,
        &loci,
        &usable,
        &panel_frac,
        &alt_frac,
        &bam_stats,
        &d,
        VERSION,
        &a.reference.display().to_string(),
        params,
        &a.lang,
    )?;
    write_run_json(&a, &infos, &loci, &panel_frac)?;
    log!("done → {}", a.out.join("report.html").display());
    Ok(())
}

fn file_name(p: &Path) -> String {
    p.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default()
}

fn write_run_json(a: &TraceArgs, infos: &[SampleInfo], loci: &[Locus], frac: &[(f64, f64)]) -> Result<()> {
    let v = serde_json::json!({
        "tool": "tr-rescue",
        "version": VERSION,
        "command": "trace",
        "created": now_string(),
        "catalog": "STRchive snapshot 2026-11-17 (CC BY 4.0)",
        "args": {
            "reference": a.reference.display().to_string(),
            "manifest": a.manifest.display().to_string(),
            "k": a.k, "min_hits": a.min_hits, "flank": a.flank, "pad": a.pad,
            "window": a.window, "bin": a.bin, "distal": a.distal, "min_mapq": a.min_mapq, "threads": a.threads
        },
        "samples": infos.iter().map(|i| serde_json::json!({
            "sample_id": i.id, "group": i.group, "role": i.role, "sex": i.sex_declared, "sex_inferred": i.sex_inferred,
            "x_ratio": i.x_ratio, "qc": i.qc, "hifi_bam": i.hifi_bam,
            "fail_bams": i.fail_bams, "hifi_raw_bams": i.hifi_raw_bams, "unassigned_bams": i.unassigned_used, "notes": i.locate_notes
        })).collect::<Vec<_>>(),
        "loci": loci.iter().enumerate().filter(|(_, l)| !l.is_control()).map(|(i, l)| serde_json::json!({
            "id": l.id, "gene": l.gene, "region": l.region_str(),
            "unique_left": frac[i].0, "unique_right": frac[i].1
        })).collect::<Vec<_>>(),
    });
    let mut f = fs::File::create(a.out.join("run.json"))?;
    writeln!(f, "{}", serde_json::to_string_pretty(&v)?)?;
    Ok(())
}
