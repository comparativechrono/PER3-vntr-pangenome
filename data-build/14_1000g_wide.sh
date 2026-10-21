#!/bin/bash
# 1000g/1kgp2504_chr1_7.58-8.08Mb.vcf.gz: the 2,504 unrelated people of the 1000 Genomes high-coverage phased panel
# (NYGC, 2022), chr1:7,580,000-8,080,000, variants present in them, without INFO - the LD reference of the
# conditional analyses. Needs: bcftools, tabix (network), and 1000g/1kgp_samples.tsv.
# Run:  bash 14_1000g_wide.sh OUTDIR SAMPLES_TSV      (made on 21 October 2026)
set -euo pipefail
OUT=${1:-.}; SMP=${2:-../data/1000g/1kgp_samples.tsv}
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
awk -F'\t' 'NR > 1 && $6 == "yes" {print $1}' "$SMP" > "$OUT/unrelated2504.txt"
bcftools view -r chr1:7580000-8080000 -S "$OUT/unrelated2504.txt" "$URL" -Ou | bcftools view -c 1 -Ou \
  | bcftools annotate -x INFO -Oz -o "$OUT/1kgp2504_chr1_7.58-8.08Mb.vcf.gz"
tabix -f -p vcf "$OUT/1kgp2504_chr1_7.58-8.08Mb.vcf.gz"
echo "$(bcftools view -H "$OUT/1kgp2504_chr1_7.58-8.08Mb.vcf.gz" | wc -l) records"      # 13204
