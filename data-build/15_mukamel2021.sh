#!/bin/bash
# mukamel2021/: the association statistics of the PER3 coding VNTR (VNTR_HG38_chr1_7829874_7830126) for the 787
# UK Biobank phenotypes of Mukamel et al. 2021 (Science 373:1499), from the authors' public archive, and the mean
# and SD of its rescaled dosage.  Run:  bash 15_mukamel2021.sh OUTDIR      (made on 21 October 2026)
set -euo pipefail
OUT=${1:-.}; T=$(mktemp -d)
B=https://data.broadinstitute.org/lohlab/UKB_codingVNTR_sumstats
curl -sS -o "$T/a.tar.gz" "$B/coding_vntr_stats.tar.gz"; curl -sS -o "$OUT/readme_broad.txt" "$B/readme.txt"
mkdir "$T/x"; tar xzf "$T/a.tar.gz" -C "$T/x"
( printf "phenotype\tA1FREQ\tF_MISS\tCHISQ_LINREG\tP_LINREG\tBETA\tSE\n"
  for f in "$T"/x/*.vntr_stats.gz; do
    zcat "$f" | awk -v f="$(basename "$f" .vntr_stats.gz)" -F'\t' '$1 == "VNTR_HG38_chr1_7829874_7830126" {print f "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}'
  done | sort ) > "$OUT/mukamel2021_PER3_VNTR_787_phenotypes.tsv"
grep "chr1_7829874" "$T"/x/vntr_list118.mean_std.txt > "$OUT/mukamel2021_PER3_VNTR_dosage_mean_sd.txt"
rm -rf "$T"
