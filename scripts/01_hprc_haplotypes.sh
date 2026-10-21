#!/usr/bin/env bash
# 01. Rebuild the 600 bp around the PER3 VNTR (GRCh38 chr1:7,829,701-7,830,300) in each of the 462 haplotypes of
#     the HPRC release 2 pangenome (Minigraph-Cactus v2.1, GRCh38-based VCF of whole bubbles), with bcftools consensus;
#     and in the 1000 Genomes haplotypes with a deletion call in the repeat (short-read structural-variant calls).
# Output: outputs/work/hprc_haplotypes.fa, outputs/work/kg_sv_haplotypes.fa, outputs/work/grch38_region.fa,
#         outputs/logs/01_consensus.log
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
VCF=$ROOT/data/hprc/hprc_per3_vntr.vcf.gz
REF=$ROOT/data/reference/GRCh38_chr1_8Mb.fa.gz
W=$ROOT/outputs/work; mkdir -p "$W" "$ROOT/outputs/logs"
LOG=$ROOT/outputs/logs/01_consensus.log; : > "$LOG"

samtools faidx "$REF" chr1:7829701-7830300 > "$W/grch38_region.fa"
# only the records that start in the region (the VCF covers 2 kb; records outside would be reported as overlapping)
bcftools view -t chr1:7829701-7830300 "$VCF" -Oz -o "$W/hprc_region.vcf.gz"
tabix -f -p vcf "$W/hprc_region.vcf.gz"
: > "$W/hprc_haplotypes.fa"
for s in $(bcftools query -l "$W/hprc_region.vcf.gz"); do
  for h in 1 2; do
    echo "== $s hap$h" >> "$LOG"
    bcftools consensus -f "$W/grch38_region.fa" -s "$s" -H "$h" "$W/hprc_region.vcf.gz" 2>> "$LOG" \
      | sed "s/^>.*/>${s}_hap${h}/" >> "$W/hprc_haplotypes.fa"
  done
done
n=$(grep -c '^>' "$W/hprc_haplotypes.fa")
echo "01: $n haplotypes rebuilt"

# The 1000 Genomes haplotypes in which the short-read panel called a deletion of 50 bp or more in exon 19: the same
# 600 bp as the panel describes them (all of its variants on that haplotype), for step 02 to read their units.
KG=$ROOT/data/1000g/1kgp_per3_gene.vcf.gz
bcftools view -t chr1:7829701-7830300 "$KG" -Oz -o "$W/kg_region.vcf.gz"
tabix -f -p vcf "$W/kg_region.vcf.gz"
bcftools query -f '%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' "$W/kg_region.vcf.gz" \
  | awk -F'\t' '$1 >= 7829834 && $1 <= 7830161 && length($2) - length($3) >= 50 {
        for (i = 4; i <= NF; i++) { split($i, a, "="); split(a[2], g, "|"); if (g[1] == "1") print a[1] "\t1"; if (g[2] == "1") print a[1] "\t2" } }' \
  | sort -u > "$W/kg_sv_carriers.txt"
: > "$W/kg_sv_haplotypes.fa"
while IFS=$'\t' read -r s h; do
  echo "== 1000G $s hap$h" >> "$LOG"
  bcftools consensus -f "$W/grch38_region.fa" -s "$s" -H "$h" "$W/kg_region.vcf.gz" 2>> "$LOG" \
    | sed "s/^>.*/>${s}_hap${h}/" >> "$W/kg_sv_haplotypes.fa"
done < "$W/kg_sv_carriers.txt"
echo "01: 1000 Genomes haplotypes with a deletion call in the repeat, rebuilt: $(grep -c '^>' "$W/kg_sv_haplotypes.fa")"
# bcftools consensus reports only 'Applied N variants' when no variant had to be skipped; anything else is a warning
problems=$(grep -c -v -E '^(== |Applied [0-9]+ variants)' "$LOG" || true)
if [ "$problems" -gt 0 ]; then
  echo "01: WARNING: bcftools consensus reported $problems problems, see $LOG"; exit 1
fi
echo "01: records used: $(bcftools view -H "$W/hprc_region.vcf.gz" | wc -l); no warnings from bcftools consensus"
