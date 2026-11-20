## TR-Rescue v1.0.0

TR-Rescue finds tandem repeat expansions in PacBio HiFi WGS that genotypers can miss, because the expanded molecules ended up in the fail BAM or left no read at all. For each sample × disease locus (STRchive, 82 loci) it reports:

- the HiFi alleles;
- the depth loss around the repeat;
- where the reads went: aligned HiFi, fail 0x01, fail 0x08, unassigned, or nowhere.

The report is a single HTML page with a Korean/English switch.

TR-Rescue는 확장된 분자가 fail BAM으로 빠지거나 read를 아예 남기지 않아 genotyper가 놓칠 수 있는 tandem repeat 확장을 PacBio HiFi WGS에서 찾습니다. 샘플 × 질환 좌위(STRchive 82개)마다 다음을 보여 줍니다.

- HiFi allele
- 반복 주변 depth 함몰
- read가 어디로 갔는지: 정렬 HiFi, fail 0x01, fail 0x08, unassigned, 또는 어디에도 없음

보고서는 HTML 한 장이며, 한국어와 영어로 전환할 수 있습니다.

### Download

`tr-rescue_v1.0.0_linux_x86_64.zip` contains:

- a static binary for Linux x86_64 and WSL (`bin/tr-rescue`; no Rust or other dependencies needed);
- the source code;
- manifest templates (.xlsx / .csv / .tsv);
- an example report built from synthetic data.

```bash
unzip tr-rescue_v1.0.0_linux_x86_64.zip && cd tr-rescue_v1.0.0
./bin/tr-rescue --help
```

### Commands

- `check`: check a manifest before running.
- `trace`: scan aligned, fail and unassigned BAMs and write the report.
- `report`: rebuild the report from an earlier trace folder with the current rules, in seconds, without opening any BAM.
- `locate`: find a sample's fail and unassigned BAMs from its @RG.
- `catalog`: list the built-in loci.

See README.md (Korean) or README.en.md (English), and CHANGELOG.md.

Research software — not a diagnosis. Code: MIT © 2026 Seung Won Chae. Locus catalog: STRchive, CC BY 4.0 (Hiatt et al., Genome Medicine 2025).
