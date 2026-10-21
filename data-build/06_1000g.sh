#!/bin/bash
# 1000g/1kgp_per3_gene.vcf.gz (+ .tbi): the phased SNVs, indels and SVs of the 3,202 people of the high-coverage
# 1000 Genomes Project (NYGC phased panel, 2022) across PER3, chr1:7,780,000-7,850,000, without INFO.
# Needs: bcftools and tabix with network access.   Run:  bash 06_1000g.sh OUTDIR
set -euo pipefail
OUT=${1:-.}
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
bcftools view -r chr1:7780000-7850000 "$URL" -Ou | bcftools annotate -x INFO -Oz -o "$OUT/1kgp_per3_gene.vcf.gz"
tabix -f -p vcf "$OUT/1kgp_per3_gene.vcf.gz"
echo "$(bcftools view -H "$OUT/1kgp_per3_gene.vcf.gz" | wc -l) records"      # 1775
