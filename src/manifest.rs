// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Sample manifest and automatic discovery of Revio fail/hifi BAMs from @RG headers.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use noodles::bam;

#[derive(Clone, Debug)]
pub struct Sample {
    pub id: String,
    /// free text as written in the manifest
    pub group: String,
    /// positive / control / test, derived from `group`
    pub role: &'static str,
    /// "M" / "F" if given in the manifest
    pub sex: Option<&'static str>,
    pub sex_raw: String,
    pub locus: String,
    pub hifi_bam: PathBuf,
    pub fail_bams: Vec<PathBuf>,
    pub fail_auto: bool,
    pub hifi_raw_bams: Vec<PathBuf>,
    pub hifi_raw_auto: bool,
    pub note: String,
}

fn split_paths(s: &str) -> (Vec<PathBuf>, bool) {
    let s = s.trim();
    if s.is_empty() || s == "-" || s == "." || s.eq_ignore_ascii_case("na") || s.eq_ignore_ascii_case("none") {
        return (vec![], false);
    }
    if s.eq_ignore_ascii_case("auto") {
        return (vec![], true);
    }
    (
        s.split([',', ';']).map(|p| PathBuf::from(p.trim())).filter(|p| !p.as_os_str().is_empty()).collect(),
        false,
    )
}

/// Map the free-text group to a role used by the analysis.
pub fn role_of(group: &str) -> &'static str {
    let g = group.trim().to_lowercase();
    const POS: [&str; 14] =
        ["positive", "pos", "case", "affected", "patient", "proband", "carrier", "양성", "환자", "확진", "발단자", "보인자", "이환", "affect"];
    const CTL: [&str; 15] = [
        "negative", "neg", "control", "ctrl", "normal", "healthy", "reference", "ref", "unaffected", "음성", "대조", "대조군", "정상",
        "참조", "비이환",
    ];
    if POS.iter().any(|k| g == *k || g.starts_with(k)) {
        "positive"
    } else if CTL.iter().any(|k| g == *k || g.starts_with(k)) {
        "control"
    } else {
        "test"
    }
}

/// Parse a declared sex: M / F / unknown.
pub fn sex_of(raw: &str) -> Option<&'static str> {
    let g = raw.trim().to_lowercase();
    match g.as_str() {
        "m" | "male" | "man" | "boy" | "xy" | "1" | "남" | "남자" | "남성" | "남아" => Some("M"),
        "f" | "female" | "woman" | "girl" | "xx" | "2" | "여" | "여자" | "여성" | "여아" => Some("F"),
        _ => None,
    }
}

fn norm(h: &str) -> String {
    h.trim_start_matches('\u{feff}')
        .trim_start_matches('#')
        .trim()
        .to_lowercase()
        .chars()
        .filter(|c| !matches!(c, ' ' | '_' | '-' | '.' | '(' | ')'))
        .collect()
}

/// Canonical column for a header cell (English and Korean aliases).
fn column_of(h: &str) -> Option<&'static str> {
    let n = norm(h);
    let is = |xs: &[&str]| xs.iter().any(|x| n == *x);
    if is(&["sampleid", "sample", "id", "samplename", "샘플", "샘플id", "샘플명", "검체", "검체id"]) {
        Some("sample_id")
    } else if is(&["group", "role", "type", "category", "그룹", "구분", "분류", "역할"]) {
        Some("group")
    } else if is(&["sex", "gender", "성별", "성"]) {
        Some("sex")
    } else if is(&["locus", "loci", "expectedlocus", "confirmedlocus", "gene", "좌위", "확진좌위", "유전자"]) {
        Some("locus")
    } else if is(&["hifibam", "bam", "alignedbam", "hiphasebam", "alignedhifibam", "정렬bam", "bam경로"]) {
        Some("hifi_bam")
    } else if is(&["failbams", "failbam", "fail", "failreads", "failreadsbam", "failbam경로"]) {
        Some("fail_bams")
    } else if is(&["hifirawbams", "hifirawbam", "hifiraw", "hifireads", "rawbam", "rawhifibam", "hifireadsbam"]) {
        Some("hifi_raw_bams")
    } else if is(&["note", "notes", "comment", "comments", "memo", "메모", "비고", "설명"]) {
        Some("note")
    } else {
        None
    }
}

/// Decode text that may be UTF-8 (with or without BOM) or CP949/EUC-KR (Korean Excel CSV).
fn decode(bytes: &[u8]) -> String {
    let b = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes);
    match std::str::from_utf8(b) {
        Ok(s) => s.to_string(),
        Err(_) => encoding_rs::EUC_KR.decode(b).0.into_owned(),
    }
}

/// Minimal RFC 4180 reader (quoted fields, doubled quotes, CRLF).
fn parse_delimited(text: &str, sep: char) -> Vec<Vec<String>> {
    let mut rows = Vec::new();
    let mut row = Vec::new();
    let mut field = String::new();
    let mut inq = false;
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if inq {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    field.push('"');
                    chars.next();
                } else {
                    inq = false;
                }
            } else {
                field.push(c);
            }
        } else if c == '"' && field.is_empty() {
            inq = true;
        } else if c == sep {
            row.push(std::mem::take(&mut field));
        } else if c == '\n' || c == '\r' {
            if c == '\r' && chars.peek() == Some(&'\n') {
                chars.next();
            }
            row.push(std::mem::take(&mut field));
            rows.push(std::mem::take(&mut row));
        } else {
            field.push(c);
        }
    }
    if !field.is_empty() || !row.is_empty() {
        row.push(field);
        rows.push(row);
    }
    rows
}

fn read_sheet(path: &Path) -> Result<Vec<Vec<String>>> {
    use calamine::{Data, Reader, open_workbook_auto};
    let mut wb = open_workbook_auto(path).with_context(|| format!("opening workbook {}", path.display()))?;
    let names = wb.sheet_names();
    // prefer a sheet named like the manifest, else the first sheet with a sample_id header
    let mut order: Vec<String> = names
        .iter()
        .filter(|n| matches!(norm(n).as_str(), "manifest" | "samples" | "sample" | "샘플" | "샘플목록"))
        .cloned()
        .collect();
    order.extend(names.iter().cloned());
    for name in order {
        let Ok(range) = wb.worksheet_range(&name) else { continue };
        let rows: Vec<Vec<String>> = range
            .rows()
            .map(|r| {
                r.iter()
                    .map(|c| match c {
                        Data::Float(f) if f.fract() == 0.0 && f.abs() < 1e15 => format!("{}", *f as i64),
                        Data::Empty => String::new(),
                        other => other.to_string(),
                    })
                    .collect()
            })
            .collect();
        if rows.iter().take(20).any(|r| r.iter().any(|c| column_of(c) == Some("sample_id"))) {
            return Ok(rows);
        }
    }
    bail!("no sheet in {} has a sample_id header", path.display())
}

/// Read a sample manifest: .tsv/.txt (tab), .csv (comma; UTF-8 or CP949), .xlsx/.xls/.ods.
pub fn read_manifest(path: &Path) -> Result<Vec<Sample>> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    let rows: Vec<Vec<String>> = match ext.as_str() {
        "xlsx" | "xlsm" | "xlsb" | "xls" | "ods" => read_sheet(path)?,
        _ => {
            let bytes = fs::read(path).with_context(|| format!("reading manifest {}", path.display()))?;
            let text = decode(&bytes);
            let first = text.lines().find(|l| !l.trim().is_empty()).unwrap_or("");
            let sep = if ext == "csv" || (!first.contains('\t') && first.contains(',')) { ',' } else { '\t' };
            parse_delimited(&text, sep)
        }
    };
    // header = first row naming a sample_id column
    let Some(hi) = rows.iter().position(|r| r.iter().any(|c| column_of(c) == Some("sample_id"))) else {
        bail!("{}: no header row with sample_id (see examples/manifest_template.*)", path.display());
    };
    let cols: Vec<Option<&'static str>> = rows[hi].iter().map(|c| column_of(c)).collect();
    if !cols.contains(&Some("hifi_bam")) {
        bail!("{}: no hifi_bam column", path.display());
    }
    let mut out = Vec::new();
    for (ri, r) in rows.iter().enumerate().skip(hi + 1) {
        let get = |name: &str| -> String {
            cols.iter().position(|c| *c == Some(name)).and_then(|i| r.get(i)).map(|s| s.trim().to_string()).unwrap_or_default()
        };
        let id = get("sample_id");
        if id.is_empty() || id.starts_with('#') {
            continue;
        }
        let hifi = get("hifi_bam");
        if hifi.is_empty() {
            bail!("{} row {}: empty hifi_bam for {}", path.display(), ri + 1, id);
        }
        let (fail_bams, fail_auto) = split_paths(&get("fail_bams"));
        let (hifi_raw_bams, hifi_raw_auto) = split_paths(&get("hifi_raw_bams"));
        let group = get("group");
        let sex_raw = get("sex");
        let sex = sex_of(&sex_raw);
        if sex.is_none() && !sex_raw.is_empty() && sex_raw != "-" && !sex_raw.eq_ignore_ascii_case("na") {
            eprintln!("[manifest] WARNING {id}: sex '{sex_raw}' not recognised (use M/F); will infer from chrX depth");
        }
        out.push(Sample {
            sex,
            sex_raw,
            role: role_of(&group),
            group: if group.is_empty() { "-".into() } else { group },
            id,
            locus: get("locus"),
            hifi_bam: PathBuf::from(hifi),
            fail_bams,
            fail_auto,
            hifi_raw_bams,
            hifi_raw_auto,
            note: get("note"),
        });
    }
    if out.is_empty() {
        bail!("manifest {} has no samples", path.display());
    }
    let mut seen = std::collections::HashSet::new();
    for s in &out {
        if !seen.insert(s.id.clone()) {
            bail!("duplicate sample_id {} in {}", s.id, path.display());
        }
    }
    Ok(out)
}

/// Files under `root` (depth <= 5) named `<movie>.<kind>.*.bam`, keyed by movie.
pub fn index_run_root(root: &Path, kind: &str) -> BTreeMap<String, Vec<PathBuf>> {
    let mut out: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    let pat = format!(".{kind}.");
    fn walk(dir: &Path, depth: usize, pat: &str, out: &mut BTreeMap<String, Vec<PathBuf>>) {
        let Ok(rd) = fs::read_dir(dir) else { return };
        for e in rd.flatten() {
            let p = e.path();
            let Ok(ft) = e.file_type() else { continue };
            if ft.is_dir() {
                if depth < 5 {
                    walk(&p, depth + 1, pat, out);
                }
            } else if let Some(name) = p.file_name().and_then(|n| n.to_str()) {
                if name.ends_with(".bam") {
                    if let Some(i) = name.find(pat) {
                        out.entry(name[..i].to_string()).or_default().push(p.clone());
                    }
                }
            }
        }
    }
    walk(root, 0, &pat, &mut out);
    out
}

pub fn rg_fields(path: &Path) -> Result<Vec<(String, String, String)>> {
    let mut r = File::open(path).map(bam::io::Reader::new)?;
    let h = r.read_header()?;
    let mut out = Vec::new();
    for (_id, rg) in h.read_groups() {
        let get = |t: [u8; 2]| {
            rg.other_fields()
                .iter()
                .find(|(k, _)| AsRef::<[u8; 2]>::as_ref(*k) == &t)
                .map(|(_, v)| v.to_string())
                .unwrap_or_default()
        };
        out.push((get(*b"PU"), get(*b"SM"), get(*b"BC")));
    }
    Ok(out)
}

#[derive(Clone, Debug, Default)]
pub struct Located {
    pub matched: Vec<PathBuf>,
    /// unassigned BAMs of single-plex cells (one barcode file, and it is this sample's)
    pub unassigned_safe: Vec<PathBuf>,
    /// unassigned BAMs of cells with several barcodes (other samples' reads mixed in)
    pub unassigned: Vec<PathBuf>,
    pub notes: Vec<String>,
}

/// Match `<movie>.<kind>.*.bam` files to the aligned BAM's read groups by SM, then BC.
pub fn locate(aligned_rgs: &[(String, String, String)], index: &BTreeMap<String, Vec<PathBuf>>) -> Located {
    let mut loc = Located::default();
    for (pu, sm, bc) in aligned_rgs {
        if pu.is_empty() {
            loc.notes.push("@RG without PU (movie name)".into());
            continue;
        }
        let Some(files) = index.get(pu) else {
            loc.notes.push(format!("{pu}: no files under run root"));
            continue;
        };
        let mut hit = false;
        let mut n_barcoded = 0usize;
        let mut un: Vec<PathBuf> = Vec::new();
        for f in files {
            let name = f.file_name().and_then(|n| n.to_str()).unwrap_or("");
            if name.contains(".unassigned.") {
                un.push(f.clone());
                continue;
            }
            n_barcoded += 1;
            let Ok(rgs) = rg_fields(f) else {
                loc.notes.push(format!("{}: unreadable header", f.display()));
                continue;
            };
            let ok = rgs.iter().any(|(_, s2, b2)| (!sm.is_empty() && s2 == sm) || (!bc.is_empty() && b2 == bc));
            if ok && !loc.matched.contains(f) {
                loc.matched.push(f.clone());
                hit = true;
            }
        }
        if !hit {
            loc.notes.push(format!("{pu}: no barcode file matched SM/BC"));
        }
        if hit && n_barcoded == 1 {
            loc.unassigned_safe.extend(un);
        } else {
            if !un.is_empty() && hit {
                loc.notes.push(format!("{pu}: {n_barcoded} barcodes in this cell, unassigned reads may be other samples' (not used)"));
            }
            loc.unassigned.extend(un);
        }
    }
    for v in [&mut loc.unassigned, &mut loc.unassigned_safe] {
        v.sort();
        v.dedup();
    }
    loc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roles() {
        for g in ["positive", "Positive", "patient", "case", "양성", "환자", "확진"] {
            assert_eq!(role_of(g), "positive", "{g}");
        }
        for g in ["reference", "control", "normal_range", "Normal", "음성", "정상", "대조군", "healthy"] {
            assert_eq!(role_of(g), "control", "{g}");
        }
        for g in ["sample", "test", "", "-", "unknown", "검사", "cohort"] {
            assert_eq!(role_of(g), "test", "{g}");
        }
    }

    #[test]
    fn sexes() {
        for x in ["M", "m", "male", "남", "남성", "1"] {
            assert_eq!(sex_of(x), Some("M"), "{x}");
        }
        for x in ["F", "female", "여", "여성", "2"] {
            assert_eq!(sex_of(x), Some("F"), "{x}");
        }
        for x in ["", "-", "unknown", "?"] {
            assert_eq!(sex_of(x), None, "{x}");
        }
    }

    #[test]
    fn headers() {
        assert_eq!(column_of("\u{feff}sample_id"), Some("sample_id"));
        assert_eq!(column_of("#Sample ID"), Some("sample_id"));
        assert_eq!(column_of("샘플"), Some("sample_id"));
        assert_eq!(column_of("성별"), Some("sex"));
        assert_eq!(column_of("Gender"), Some("sex"));
        assert_eq!(column_of("확진 좌위"), Some("locus"));
        assert_eq!(column_of("HiFi BAM"), Some("hifi_bam"));
        assert_eq!(column_of("fail_bams"), Some("fail_bams"));
        assert_eq!(column_of("비고"), Some("note"));
        assert_eq!(column_of("whatever"), None);
    }

    #[test]
    fn csv_quotes() {
        let rows = parse_delimited("a,b,c\r\n1,\"x, y\",\"he said \"\"hi\"\"\"\r\n2,,\n", ',');
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[1], vec!["1", "x, y", "he said \"hi\""]);
        assert_eq!(rows[2], vec!["2", "", ""]);
    }
}
