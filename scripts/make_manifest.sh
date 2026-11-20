#!/usr/bin/env bash
# Draft a tr-rescue manifest from a Revio WGS project folder that has HiPhase_QC/.
# sample_id comes from a <PREFIX>_<NNN> sample folder name if there is one, else S01, S02 …
# (keep hospital IDs out of sample_id);
# fill in group / locus / note by hand afterwards.
#
# usage: scripts/make_manifest.sh <project_dir> > manifest.tsv
set -euo pipefail
P=${1:?usage: make_manifest.sh <project_dir>}
printf "sample_id\tgroup\tsex\tlocus\thifi_bam\tfail_bams\thifi_raw_bams\tnote\n"
n=0
for bam in "$P"/HiPhase_QC/*.hiphase.haplotagged.snp.vaf30.bam; do
  [ -e "$bam" ] || continue
  stem=$(basename "$bam" .hiphase.haplotagged.snp.vaf30.bam)
  n=$((n + 1))
  id=""
  for d in "$P"/*"$stem"*; do
    [ -d "$d" ] || continue
    id=$(basename "$d" | grep -oE '^[A-Za-z]+_[0-9]+' || true)
    [ -n "$id" ] && break
  done
  [ -z "$id" ] && id=$(printf "S%02d" "$n")
  printf "%s\tsample\t\t-\t%s\tauto\t-\t\n" "$id" "$bam"
done
