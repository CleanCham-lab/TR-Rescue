#!/usr/bin/env python3
"""Write manifest templates (examples/manifest_template.{xlsx,csv,tsv}).

Usage:
  python3 scripts/make_manifest_template.py            # templates into examples/
  python3 scripts/make_manifest_template.py in.tsv out.xlsx   # convert an existing TSV manifest to xlsx
"""
import csv
import os
import sys

from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.worksheet.datavalidation import DataValidation

COLS = ["sample_id", "group", "sex", "locus", "hifi_bam", "fail_bams", "hifi_raw_bams", "note"]
WIDTH = {"sample_id": 12, "group": 12, "sex": 6, "locus": 14, "hifi_bam": 70, "fail_bams": 14, "hifi_raw_bams": 14, "note": 40}
HELP = [
    ("sample_id", "필수", "샘플 이름. 보고서에 그대로 나오므로 환자명·등록번호 대신 가명(S01 등)을 쓰세요. 중복 불가."),
    ("group", "선택", "자유 입력. positive·patient·case·양성·환자·확진 → 양성 / control·negative·normal·reference·음성·정상·대조 → 대조 / "
                    "그 밖(sample·test·빈칸) → 검사. 적은 글자는 보고서에 그대로 표시되고 역할 필터로도 고를 수 있습니다."),
    ("sex", "권장", "M/F (남/여, male/female 도 됨). 비우면 chrX depth로 추정합니다. 적은 값과 depth가 다르면 ⚠ 경고가 뜹니다 "
                  "(샘플 바뀜·성염색체 이수성 확인용). 남성이면 X 좌위에서 haplotype 비교를 하지 않습니다."),
    ("locus", "선택", "이미 알고 있는 확진 좌위. 유전자명(FMR1)이나 STRchive ID(FXS_FMR1). 여러 개면 쉼표로. 없으면 - 또는 빈칸. "
                    "분석 대상 좌위를 제한하지는 않습니다(전 좌위 분석). 보고서에 '확진' 표시와 재현율 집계에만 씁니다."),
    ("hifi_bam", "필수", "정렬 BAM 경로 (HiPhase haplotag 본 권장, .bai 필요). WSL 경로(/mnt/...)로 적으세요."),
    ("fail_bams", "선택", "auto = 정렬 BAM의 @RG(PU 무비 이름, SM, BC)로 --run-root 아래 <movie>.fail_reads.*.bam 을 찾음. "
                        "직접 적을 땐 경로를 쉼표로. 안 쓰면 -."),
    ("hifi_raw_bams", "선택", "원본 hifi_reads BAM (정렬본에서 빠진 HiFi read 확인용). 느려서 기본은 -. auto 또는 경로."),
    ("note", "선택", "자유 메모. 보고서에 표시만 되고 분석에는 쓰지 않습니다."),
]
EXAMPLE = ["SAMPLE_01", "positive", "F", "DMPK", "/path/to/SAMPLE_01.hiphase.haplotagged.bam", "auto", "-", "예시 행 — 지우고 쓰세요 / example row, overwrite it"]

FONT = "Arial"
thin = Side(style="thin", color="BBBBBB")


def style_header(ws, ncol):
    for c in range(1, ncol + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = Font(name=FONT, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="1F4E79")
        cell.alignment = Alignment(vertical="center")
        cell.border = Border(bottom=thin)
    ws.freeze_panes = "B2"


def add_validation(ws, nrows):
    g = DataValidation(type="list", formula1='"positive,control,test"', allow_blank=True, errorStyle="information",
                       showErrorMessage=True, errorTitle="group", error="목록에 없는 값도 됩니다 (자유 입력). 계속하려면 확인.",
                       promptTitle="group", prompt="positive / control / test 또는 자유 입력", showInputMessage=True)
    s = DataValidation(type="list", formula1='"M,F"', allow_blank=True, errorStyle="information", showErrorMessage=True,
                       errorTitle="sex", error="M 또는 F (남/여도 인식). 계속하려면 확인.")
    f = DataValidation(type="list", formula1='"auto,-"', allow_blank=True, errorStyle="information", showErrorMessage=True,
                       errorTitle="fail_bams", error="auto, - 또는 경로(쉼표 구분). 계속하려면 확인.")
    for dv in (g, s, f):
        ws.add_data_validation(dv)
    g.add(f"B2:B{nrows}")
    s.add(f"C2:C{nrows}")
    f.add(f"F2:G{nrows}")


def help_sheet(wb):
    h = wb.create_sheet("설명")
    h.append(["열", "필수", "설명"])
    for r in HELP:
        h.append(list(r))
    h.append([])
    h.append(["", "", "manifest 시트의 흰 칸을 채우면 됩니다. 2행은 예시이니 지우거나 덮어쓰세요. '#'로 시작하는 sample_id 행은 건너뜁니다."])
    h.append(["", "", "같은 내용을 .csv(UTF-8 또는 Excel CSV)나 .tsv로 저장해도 됩니다. 열 순서는 상관없고, 열 이름은 한글(샘플, 그룹, 성별, 확진 좌위, 비고)도 인식합니다."])
    h.append(["", "", "실행 전 확인: tr-rescue check --manifest 이파일.xlsx --run-root <Revio run 폴더>"])
    for c in range(1, 4):
        h.cell(row=1, column=c).font = Font(name=FONT, bold=True, color="FFFFFF")
        h.cell(row=1, column=c).fill = PatternFill("solid", fgColor="1F4E79")
    for row in h.iter_rows(min_row=2):
        for cell in row:
            cell.font = Font(name=FONT, bold=cell.column == 1)
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    h.column_dimensions["A"].width = 14
    h.column_dimensions["B"].width = 6
    h.column_dimensions["C"].width = 110


HELP_EN = [
    ("sample_id", "required", "Sample name, shown as written in the report. Use a pseudonym (S01 …), not a patient name or record number. Must be unique."),
    ("group", "optional", "Free text. positive/patient/case → positive; control/negative/normal/reference → control; anything else (sample, test, empty) → test. "
                         "The text itself is shown in the report and can be filtered on."),
    ("sex", "recommended", "M/F (male/female also accepted). If empty, inferred from chrX depth; if it disagrees with the depth a ⚠ warning appears "
                            "(sample swap / sex-chromosome aneuploidy). Male → no haplotype comparison at X-linked loci."),
    ("locus", "optional", "Diagnosis already known: gene (FMR1) or STRchive ID (FXS_FMR1); comma-separate several; '-' or empty if none. "
                         "Does not restrict the analysis (all loci are analysed); used for the 'Known dx' tag and the reproduction count."),
    ("hifi_bam", "required", "Aligned BAM (HiPhase-tagged recommended, .bai needed), as a Linux/WSL path."),
    ("fail_bams", "optional", "auto = find <movie>.fail_reads.*.bam under --run-root from the aligned BAM's @RG (PU movie, SM, BC). "
                             "Or comma-separated paths. '-' = none."),
    ("hifi_raw_bams", "optional", "Raw hifi_reads BAM (to find HiFi reads missing from the aligned BAM). Slow, so '-' by default. auto or paths."),
    ("note", "optional", "Free text. Shown in the report only, not used in the analysis."),
]


def help_sheet_en(wb):
    h = wb.create_sheet("Help (EN)")
    h.append(["Column", "Need", "Description"])
    for r in HELP_EN:
        h.append(list(r))
    h.append([])
    h.append(["", "", "Fill in the white cells of the 'manifest' sheet; row 2 is an example — overwrite or delete it. Rows whose sample_id starts with '#' are skipped."])
    h.append(["", "", "The same table can be saved as .csv (UTF-8 or Excel CSV) or .tsv. Column order does not matter; Korean column names are recognised too."])
    h.append(["", "", "Check before running: tr-rescue check --manifest this.xlsx --run-root <Revio run folder> --lang en"])
    for c in range(1, 4):
        h.cell(row=1, column=c).font = Font(name=FONT, bold=True, color="FFFFFF")
        h.cell(row=1, column=c).fill = PatternFill("solid", fgColor="1F4E79")
    for row in h.iter_rows(min_row=2):
        for cell in row:
            cell.font = Font(name=FONT, bold=cell.column == 1)
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    h.column_dimensions["A"].width = 14
    h.column_dimensions["B"].width = 12
    h.column_dimensions["C"].width = 110


def write_xlsx(path, header, rows, example_row=False):
    wb = Workbook()
    ws = wb.active
    ws.title = "manifest"
    ws.append(header)
    for r in rows:
        ws.append(r)
    style_header(ws, len(header))
    for i, name in enumerate(header, 1):
        ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = WIDTH.get(name, 14)
        ws.cell(row=1, column=i).comment = Comment(next((h[2] for h in HELP if h[0] == name), ""), "TR-Rescue", width=320, height=140)
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.font = Font(name=FONT, color="808080" if example_row and cell.row == 2 else "000000",
                             italic=example_row and cell.row == 2)
    if "sex" in header:  # highlight cells still to fill
        c = header.index("sex") + 1
        for r in range(2, len(rows) + 2):
            if not ws.cell(row=r, column=c).value:
                ws.cell(row=r, column=c).fill = PatternFill("solid", fgColor="FFFF00")
    add_validation(ws, max(200, len(rows) + 50))
    help_sheet(wb)
    help_sheet_en(wb)
    wb.save(path)


def main():
    if len(sys.argv) == 3:
        src, dst = sys.argv[1:]
        rows = list(csv.reader(open(src, encoding="utf-8-sig"), delimiter="\t"))
        header = rows[0]
        if "sex" not in header:
            gi = header.index("group") + 1 if "group" in header else 1
            header = header[:gi] + ["sex"] + header[gi:]
            n = len(rows[0])
            rows = [rows[0]] + [r[:gi] + [""] + r[gi:] if len(r) == n else r for r in rows[1:]]
        write_xlsx(dst, header, [r for r in rows[1:] if r and r[0]])
        print("wrote", dst)
        return
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples")
    os.makedirs(out, exist_ok=True)
    write_xlsx(os.path.join(out, "manifest_template.xlsx"), COLS, [EXAMPLE], example_row=True)
    with open(os.path.join(out, "manifest_template.csv"), "w", newline="", encoding="utf-8-sig") as f:
        csv.writer(f).writerows([COLS, EXAMPLE])
    with open(os.path.join(out, "manifest_template.tsv"), "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows([COLS, EXAMPLE])
    print("wrote templates to", os.path.normpath(out))


if __name__ == "__main__":
    main()
