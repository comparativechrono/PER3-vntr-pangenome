#!/bin/bash
# reference/dbsnp157_per3_500kb.vcf.gz (+ .tbi): the variants of dbSNP build 157 between chr1:7,580,000 and
# 8,080,000 (GRCh38) - positions, rsIDs and alleles, without INFO; the chromosome renamed from NC_000001.11 to chr1.
# Needs: bcftools, bgzip, tabix (htslib with network access).   Run:  bash 03_dbsnp.sh OUTDIR
set -euo pipefail
OUT=${1:-.}
URL=https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz
{
  printf '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##source=dbSNP build 157 (NCBI), GRCh38; positions chr1:7,580,000-8,080,000; INFO removed\n'
  printf '##contig=<ID=chr1,length=248956422,assembly=GRCh38>\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  bcftools query -r NC_000001.11:7580000-8080000 -f 'chr1\t%POS\t%ID\t%REF\t%ALT\t.\t.\t.\n' "$URL"
} | bgzip > "$OUT/dbsnp157_per3_500kb.vcf.gz"
tabix -f -p vcf "$OUT/dbsnp157_per3_500kb.vcf.gz"
bcftools view -H "$OUT/dbsnp157_per3_500kb.vcf.gz" | wc -l      # 209029
