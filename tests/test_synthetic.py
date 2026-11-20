#!/usr/bin/env python3
"""End-to-end self-test: synthetic data -> tr-rescue trace -> expected calls.

Usage: python3 tests/test_synthetic.py [path/to/tr-rescue]   (needs pysam)
"""
import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "bin", "tr-rescue")


def main():
    tmp = tempfile.mkdtemp(prefix="trrescue_test_")
    subprocess.run([sys.executable, os.path.join(HERE, "make_synthetic.py"), tmp], check=True, stdout=subprocess.DEVNULL)
    out = os.path.join(tmp, "out")
    subprocess.run(
        [BIN, "trace", "--reference", f"{tmp}/ref.fa", "--manifest", f"{tmp}/manifest.tsv", "--catalog", f"{tmp}/catalog.bed",
         "--run-root", f"{tmp}/run", "--out", out, "--threads", "2", "--loci", "TGAA,TRFC,TX"],
        check=True, stderr=subprocess.DEVNULL,
    )
    summ = {(r["sample_id"], r["locus"]): r for r in csv.DictReader(open(f"{out}/summary.tsv"), delimiter="\t")}
    reads = list(csv.DictReader(open(f"{out}/reads.tsv"), delimiter="\t"))
    fails = []

    def check(cond, msg):
        print(("PASS " if cond else "FAIL ") + msg)
        if not cond:
            fails.append(msg)

    def n(sample, locus, source, states):
        return sum(1 for r in reads if r["sample_id"] == sample and r["locus"] == locus and r["source"] == source and r["state"] in states)

    exp = ("EXPANDED", "EXPANDED_LB")
    check(n("P1", "TGAA", "fail_subread", exp) >= 10, "P1 GAA expansion found in fail 0x08 reads")
    check(n("P1", "TGAA", "hifi_aligned", exp) >= 1, "P1 GAA expansion found in HiFi reads")
    spans = sorted(int(r["allele_bp"]) for r in reads if r["sample_id"] == "P1" and r["locus"] == "TGAA"
                   and r["source"] == "hifi_aligned" and r["class"] == "SPANNING" and r["state"] == "EXPANDED")
    check(bool(spans) and all(b in (1800, 2400) for b in spans), f"P1 HiFi spanning lengths exact (1800/2400): {spans}")
    check("두 haplotype" in summ[("P1", "TGAA")]["hint"], "P1 GAA: both-haplotype canyon")
    check(float(summ[("P1", "TGAA")]["fail_subread.rel_to_ctl"]) > 2.0, "P1 GAA: fail 0x08 enriched vs neighbourhood controls")
    check(all(n("N1", l, s, exp) == 0 for l in ("TGAA", "TRFC") for s in ("hifi_aligned", "fail_subread", "fail_lowqv")),
          "N1: no expanded reads")
    check(summ[("N1", "TGAA")]["hint"] == "-" and summ[("N1", "TRFC")]["hint"] == "-", "N1: no canyon hints")
    check("두 haplotype" in summ[("P2", "TRFC")]["hint"], "P2 AAGGG: both-haplotype canyon (molecules mostly lost)")
    p2 = [r for r in reads if r["sample_id"] == "P2" and r["locus"] == "TRFC" and r["state"] == "EXPANDED"]
    check(all(r["dominant_motif"] == "AAGGG" for r in p2), "P2: dominant motif AAGGG")
    f31, f32 = float(summ[("P3", "TGAA")]["f_hp1"]), float(summ[("P3", "TGAA")]["f_hp2"])
    check(f31 > 0.5 and f31 > f32 + 0.5, f"P3 one-allele, narrow phase block: HP1 lost, HP2 kept (f_hp1 {f31:.2}, f_hp2 {f32:.2})")
    check("두 haplotype" not in summ[("P3", "TGAA")]["hint"], "P3: not called a two-allele pattern")
    check(n("P3", "TGAA", "fail_subread", exp) + n("P3", "TGAA", "hifi_aligned", exp) >= 2, "P3: expanded reads found")
    check(all("한 haplotype" not in summ[(p, "TRFC")]["hint"] for p in ("P1", "N1", "P3")), "no single-haplotype hint at unexpanded TRFC")
    check(all(float(summ[("P1", "TGAA")][k]) > 0.5 for k in ("f_hp1", "f_hp2")), "P1 biallelic: both haplotypes lost")
    sx = {p: summ[(p, "TX")]["sex_inferred"] for p in ("P1", "N1", "P2", "P3")}
    check(sx == {"P1": "M", "N1": "F", "P2": "F", "P3": "M"}, f"sex inferred from chrX depth: {sx}")
    check(summ[("P1", "TX")]["f_hp1"] == "NA" and "HP 비교 없음" in summ[("P1", "TX")]["hint"], "male chrX: no HP fit, hint says so")
    check(summ[("N1", "TX")]["level"] == "ok" and summ[("P1", "TX")]["level"] == "ok", "chrX normal locus not flagged for sex alone")
    check(summ[("P1", "TGAA")]["role"] == "positive" and summ[("N1", "TGAA")]["role"] == "control", "group -> role (positive/negative)")

    # ---- manifest formats: CSV (UTF-8 BOM, Korean headers), CSV (CP949), XLSX (guide sheet first, header not on row 1)
    rows = list(csv.reader(open(f"{tmp}/manifest.tsv"), delimiter="\t"))
    hdr, body = rows[0], rows[1:]
    ko = {"sample_id": "샘플", "group": "그룹", "locus": "확진 좌위", "hifi_bam": "HiFi BAM", "fail_bams": "fail_bams", "hifi_raw_bams": "hifi_raw_bams", "note": "비고"}
    groups = {"P1": "환자", "N1": "reference", "P2": "Positive", "P3": "sample"}
    sexes = {"P1": "M", "N1": "여", "P2": "", "P3": "female"}
    out_rows = [[ko[h] for h in hdr[:2]] + ["성별"] + [ko[h] for h in hdr[2:]]]
    for r in body:
        r = list(r)
        r[1] = groups[r[0]]
        r[-1] = r[-1] + ", 쉼표 \"따옴표\" 포함"
        out_rows.append(r[:2] + [sexes[r[0]]] + r[2:])
    with open(f"{tmp}/manifest_utf8.csv", "w", newline="", encoding="utf-8-sig") as f:
        csv.writer(f).writerows(out_rows)
    with open(f"{tmp}/manifest_cp949.csv", "w", newline="", encoding="cp949") as f:
        csv.writer(f).writerows(out_rows)
    try:
        import openpyxl
        wb = openpyxl.Workbook()
        g = wb.active
        g.title = "설명"
        g.append(["이 시트는 설명입니다"])
        m = wb.create_sheet("manifest")
        m.append(["TR-Rescue manifest"])
        m.append([])
        for r in out_rows:
            m.append(r)
        wb.save(f"{tmp}/manifest.xlsx")
        have_xlsx = True
    except ImportError:
        have_xlsx = False
        print("SKIP xlsx (no openpyxl)")
    for name in ["manifest_utf8.csv", "manifest_cp949.csv"] + (["manifest.xlsx"] if have_xlsx else []):
        r = subprocess.run([BIN, "check", "--manifest", f"{tmp}/{name}", "--catalog", f"{tmp}/catalog.bed", "--run-root", f"{tmp}/run"],
                           capture_output=True, text=True)
        o = r.stdout
        good = (r.returncode == 0 and o.count("역할 양성") == 2 and o.count("역할 대조") == 1 and o.count("역할 검사") == 1
                and o.count("→ 남") == 1 and o.count("→ 여") == 2 and "쉼표 \"따옴표\" 포함" in o and "fail   auto →" in o)
        check(good, f"check reads {name} (roles, sex, quoted note, auto fail BAM)")
        if not good:
            print(o, r.stderr)
    mf = f"{tmp}/manifest.xlsx" if have_xlsx else f"{tmp}/manifest_cp949.csv"
    out2 = os.path.join(tmp, "out2")
    subprocess.run(
        [BIN, "trace", "--reference", f"{tmp}/ref.fa", "--manifest", mf, "--catalog", f"{tmp}/catalog.bed",
         "--run-root", f"{tmp}/run", "--out", out2, "--threads", "2", "--loci", "TGAA,TRFC,TX"],
        check=True, stderr=subprocess.DEVNULL,
    )
    summ2 = {(r["sample_id"], r["locus"]): r for r in csv.DictReader(open(f"{out2}/summary.tsv"), delimiter="\t")}
    same = all(summ2[k][c] == summ[k][c] for k in summ for c in ("level", "hifi_alleles", "range", "f_total") + (() if k[1] == "TX" else ("hint",)))
    check(same, f"{os.path.basename(mf)} manifest gives the same results as TSV")
    check(summ2[("P1", "TGAA")]["sex"] == "M" and summ2[("N1", "TGAA")]["sex"] == "F" and summ2[("P2", "TRFC")]["sex"] == "NA",
          "sex column carried to summary.tsv")
    run2 = open(f"{out2}/run.json").read()
    check("성별 불일치" in run2 and summ2[("P3", "TX")]["sex"] == "F" and summ2[("P3", "TX")]["f_hp1"] == "NA",
          "declared F but chrX says one X -> sex mismatch warning, HP off")
    check("남성 (manifest)" in summ2[("P1", "TX")]["hint"], "declared male used for chrX hint")
    check(summ2[("P3", "TGAA")]["role"] == "test" and summ2[("P2", "TRFC")]["role"] == "positive", "free-text group mapped (sample -> test)")
    # report rebuild, English
    out_en = os.path.join(tmp, "out_en")
    subprocess.run([BIN, "report", "--from", out, "--out", out_en, "--catalog", f"{tmp}/catalog.bed", "--lang", "en"],
                   check=True, stderr=subprocess.DEVNULL)
    se = {(r["sample_id"], r["locus"]): r for r in csv.DictReader(open(f"{out_en}/summary.tsv"), delimiter="\t")}
    hangul = any("\uac00" <= ch <= "\ud7a3" for r in se.values() for ch in r["hint"] + r["hifi_alleles"] + r["evidence"])
    check("pathogenic-range allele in HiFi" in se[("P1", "TGAA")]["hint"] and not hangul, "report --lang en: English hints, no Korean left")
    check(all(se[k]["level"] == summ[k]["level"] for k in summ), "report rebuild gives the same calls as trace")
    html = open(f"{out_en}/report.html", encoding="utf-8").read()
    check('"lang":"en"' in html and "확인 필요" in html and "Review" in html, "report.html carries both languages, opens in English")
    stats = list(csv.DictReader(open(f"{out}/bam_stats.tsv"), delimiter="\t"))
    un = {(r["sample_id"], r["kind"]) for r in stats if "unassigned" in r["kind"]}
    check(un == {("P2", "fail_reads (unassigned)"), ("P2", "hifi_reads (unassigned)")},
          f"unassigned used only for the single-plex cell (P2), not the 4-sample cell: {sorted(un)}")
    uh = n("P2", "TRFC", "unassigned_hifi", exp)
    uf = sum(1 for r in reads if r["sample_id"] == "P2" and r["locus"] == "TRFC" and r["source"] == "unassigned_fail")
    check(uh == 3 and uf > 0, f"P2 unassigned reads counted separately (HiFi expanded {uh}, fail {uf})")
    check("unassigned_hifi.EXPANDED" in summ[("P2", "TRFC")], "summary.tsv has unassigned columns")
    check(os.path.getsize(f"{out}/report.html") > 10000, "report.html written")
    print("\nall passed" if not fails else f"\n{len(fails)} failed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
