#!/bin/bash
# hprc/hprc_pggb_per3_7.78-7.86Mb.vcf.gz: the HPRC release 2 PGGB graph (GRCh38-based VCF of 2025-09-30, p98-k311,
# 'laced'), chr1:7,780,000-7,860,000, all samples (CHM13 included). The VCF on S3 has no index, so the file is
# streamed from its start and cut once past the region (chr1 comes first; records are sorted).
# Needs: curl, bgzip, tabix.   Run:  bash 13_hprc_pggb.sh OUTDIR      (made on 21 October 2026)
set -euo pipefail
OUT=${1:-.}
U=https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/pggb/vcfs/whole-genome/20250930_hprc25272.GRCh38.p98-k311.laced.vcf.gz
curl -sS "$U" | zcat 2>/dev/null \
  | awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} {if ($2 > 7860000) exit; if ($2 >= 7780000) print}' \
  | bgzip > "$OUT/hprc_pggb_per3_7.78-7.86Mb.vcf.gz" || true
tabix -f -p vcf "$OUT/hprc_pggb_per3_7.78-7.86Mb.vcf.gz"
echo "$(zcat "$OUT/hprc_pggb_per3_7.78-7.86Mb.vcf.gz" | grep -vc '^#') records"      # 1115
