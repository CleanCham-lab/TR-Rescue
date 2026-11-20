# Changelog

## 1.0.0 — 2026-11-20
First public release.
- Report: browser back / forward (and mouse side buttons) move between tabs and loci — every view has its own
  address (`#summary`, `#analysis/<sample>/<locus>` …); the overview fills wide screens (banner with an
  at-a-glance panel: calls, known-diagnosis reproduction, sample warnings, all clickable); text in overview cards
  uses the card width.
- Copyright, affiliation (Seoul National University Cancer Research Institute; Molecular Diagnostics Laboratory, Department of
  Laboratory Medicine, Seoul National University Hospital), license and citation in the report (sidebar and overview),
  `--help`, README, NOTICE and source headers.
- README in Korean and English without site-specific paths.

## 0.3.0 — 2026-11-19
- Korean / English: the report has a KR/EN switch (remembered per browser); every hint, warning, allele and evidence
  text is carried in both languages. `--lang ko|en` on `trace`, `report` and `check` sets the summary.tsv hint language,
  the report's starting language and the check output. README.en.md; English help sheet in the xlsx manifest template.
- `report` recomputes sample warnings (positive without locus, sex mismatch, low depth) in both languages.

## 0.2.3 — 2026-11-19
- New `report` command: rebuild summary/reads/report.html from an earlier trace folder with the current rules (seconds,
  no BAM access). `trace --samples` to run a subset of the manifest.
- Verdicts compare samples: with >= 4 samples, a signal present in at least half of the other samples at a locus is
  marked as a locus property and demoted to info; depth loss must exceed the cohort's robust spread (4 SD, >= 0.7 for
  watch/high; 3 SD, >= 0.5 for info); HiFi read reduction is relative to the cohort median ratio.
- Fail/unassigned-only expansions need reads >= pathogenic threshold + 3 units and >= 1.5x the longest HiFi allele;
  1-2 reads info, 3-4 watch, >= 5 high. Enrichment needs >= 3 usable controls, >= 4 reads and Poisson p < 0.001.
- Thresholds: whole-unit (rounded) comparison; contraction loci (MIR7-2); pathogenic minimum below normal maximum fixed
  (FAME loci); reference allele above the normal maximum flagged (POLG). HiFi lower bounds raise the range only when
  two reads agree and the read state is EXPANDED_LB (motif and purity checked) — fixes FAME1 TTTTA long alleles called
  pathogenic.
- Samples with median depth < 5x get a '데이터 부족' warning and their loci are not judged.
- Logo, banner, favicon in the report (and a small easter egg); README banner. Summary opens with high+watch.

## 0.2.2 — 2026-11-17
- Unassigned BAMs: with `auto`, `<movie>.hifi_reads/fail_reads.unassigned.bam` of single-plex cells (one barcode file,
  matched to the sample) are scanned and counted as separate sources `unassigned_hifi` / `unassigned_fail`
  (strip plot rows, whereabouts, summary columns, enrichment vs neighbourhood controls). Multiplexed cells are skipped.
  `--unassigned auto|never|always`. Unassigned files listed in the manifest are also counted as unassigned.
- Repeat counts: per-locus offset for non-repeat bp inside the catalog region; thresholds and units use it, bp unchanged.
  Fixed against the ExpansionHunter catalog / clinical convention: SBMA_AR (−34 bp), FRAXE_AFF2 (−89 bp; v0.2.0–0.2.1
  over-counted by ~30 and could call normal alleles intermediate), SCA8_ATXN8OS (region extended to CTA, CTA+CTG combined).
  Run-time QC lists other loci whose reference region has non-motif edges (log + loci tab).
- Degenerate motifs (GCN/NGC, polyalanine loci incl. ARX, PHOX2B, PABPN1, HOXA13, SOX3, ZIC2/3, RUNX2, FOXL2, TBX1,
  HOXD13) are expanded to concrete codons; v0.2.0–0.2.1 could not measure purity there and could miss expansions.
- Report: '확진' renamed '기존 진단' (manifest locus, not a call) with tooltip; '재현 안 된 것만' filter; status row
  '기존 진단 재현' opens the not-reproduced rows; analysis header flags '이 데이터에서 재현 안 됨'.
- Run tab: every file name has the full path on hover and copies it on click (BAM scan table had no tooltip).

## 0.2.1 — 2026-11-17
- Manifest: .xlsx / .xls / .ods (sheet named `manifest`, else the first sheet with a sample_id header), .csv (UTF-8 with or
  without BOM, CP949/EUC-KR from Korean Excel; quoted fields) and .tsv. Columns found by name in any order; Korean
  header aliases (샘플, 그룹, 성별, 확진 좌위, 비고 …); header need not be on row 1; `#` rows skipped; duplicate IDs rejected.
- `group` is free text, mapped to a role (positive / control / test) shown as a badge and filter; the text is kept.
  Positive samples without a usable `locus` get a warning.
- New `sex` column (M/F, 남/여, male/female, 1/2). Declared male → one X; otherwise chrX/autosome depth ratio decides
  (<0.7 one X, ≥0.8 two X). Declared vs depth mismatch → sample warning (⚠ in report, log, run.json). summary.tsv gains
  `role`, `sex`, `sex_inferred`.
- New `check` command: shows how each manifest row is read (role, sex, locus mapping, BAM/index presence, auto fail BAMs).
- Report: level chips now filter strictly (confirmed loci no longer bypass the filter; '확진 좌위만' turns all levels on);
  role filter with manifest group values; role/sex badges; sample warnings in status panel, summary and analysis.
- Templates: `examples/manifest_template.{xlsx,csv,tsv}` (xlsx has a 설명 sheet and dropdowns).
- Synthetic test: chrX locus with male/female samples; manifest format round-trip.

## 0.2.0 — 2026-11-17
- Locus catalog: STRchive snapshot (82 disease loci, hg38, CC BY 4.0) embedded; default analyses all loci;
  loci addressable by gene symbol or STRchive ID; same coordinates as TRGT's STRchive catalog.
- Read placement: at least 10 consistent flank k-mer hits per side (`--min-hits`); prefilter requires hits for one
  locus. Removes chance placements of noisy subreads that inflated fail-read counts in v0.1.x.
- Repeat purity >= 0.55 required for expanded/intermediate; low-purity long inserts -> NONREPEAT.
- BENIGN_MOTIF state for long alleles whose motif is not pathogenic (e.g. RFC1 AAAAG/AAAGG).
- HiFi allele summary (up to two clusters, length/units/HP/range, size-spread flag) and per-result level
  (high / watch / info / ok) with reasons.
- Haplotype fit needs >= 4x per haplotype in the phase block; HP analysis off for single-X samples.
- Panel k-mers may be shared by neighbouring loci; panel cache only with `--panel-cache`.
- Read analysis runs inside the scan threads (sequences not kept); aligned BAM processed in parallel.
- New report: full-screen layout with left tab bar (개요/요약/분석/좌위 정보/실행 정보), status panel, summary
  filters, one sample×locus per analysis view, plain-language help, hover tooltips, STRchive/STRipy/gnomAD/OMIM links;
  works on landscape, portrait and narrow screens.

## 0.1.3 — 2026-11-16
- Scan progress every 20 s (compressed GB read, MB/s, reads, ETA).
- Haplotype canyon: restricted to the phase set covering the repeat; each haplotype's baseline is taken
  inside that phase block where the expected canyon is shallow (g <= 0.1), so narrow HiPhase blocks work.
- Single-haplotype hint also requires some total-depth loss (z >= 1.5).
- Synthetic test: one-allele expansion with a narrow phase block.

## 0.1.2 — 2026-11-16
- Built-in AR (SBMA) locus.

## 0.1.1 — 2026-11-16
- Built-in loci: DMPK, FMR1, ATXN8OS, HTT, NOTCH2NLC (hg38, PanelApp thresholds) in addition to FXN, RFC1.
- FXN and RFC1 are always analysed (unless `--loci` is given), so every sample serves as their negative control.
- `--run-root` accepts several roots (e.g. revio/revio_run and revio_30/revio_run).
- chrX/autosome baseline ratio for X-linked loci (FMR1) with a single-copy hint.

## 0.1.0 — 2026-11-16
- `trace` (step 0): flank k-mer panel with reference-uniqueness and ALT-sharing checks; aligned HiFi depth canyon
  (total/HP1/HP2) with whole-read-loss templates, f fit and bootstrap SE; streaming scan of unaligned fail/hifi BAMs;
  sequence-based read placement (spanning / anchored / flank-only) with allele length and motif composition;
  neighbourhood control pseudo-loci; TSV/JSON/HTML outputs.
- `locate`: match fail_reads / hifi_reads BAMs to an aligned BAM by @RG (PU, SM, BC); unassigned BAMs listed, not used.
- `catalog`: built-in FXN, RFC1 (hg38).
- Static x86_64 Linux binary; synthetic end-to-end test.
