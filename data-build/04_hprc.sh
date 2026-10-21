#!/bin/bash
# hprc/hprc_per3_vntr.vcf.gz and hprc/hprc_per3_500kb.vcf.gz (+ .tbi): the HPRC release 2 pangenome
# (Minigraph-Cactus v2.1, GRCh38-based VCFs) around PER3, without the reference haplotype CHM13 and without the
# INFO fields AF, NS and AT and the graph's IDs.
#   - hprc-v2.1-mc-grch38.vcf.gz       the graph's bubbles, one record each, with whole alleles: chr1:7,829,001-7,831,000
#   - hprc-v2.1-mc-grch38.wave.vcf.gz  the same, decomposed into SNVs, indels and larger variants: chr1:7,580,000-8,080,000
# Needs: bcftools and tabix with network access.   Run:  bash 04_hprc.sh OUTDIR
set -euo pipefail
OUT=${1:-.}
B=https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38
bcftools view -r chr1:7829001-7831000 "$B/hprc-v2.1-mc-grch38.vcf.gz" -Ou \
  | bcftools view -s ^CHM13 -Ou \
  | bcftools annotate -x ID,INFO/AT,INFO/AF,INFO/NS -Oz -o "$OUT/hprc_per3_vntr.vcf.gz"
bcftools view -r chr1:7580000-8080000 "$B/hprc-v2.1-mc-grch38.wave.vcf.gz" -Ou \
  | bcftools view -s ^CHM13 -Ou \
  | bcftools annotate -x ID,INFO/AF,INFO/NS -Oz -o "$OUT/hprc_per3_500kb.vcf.gz"
for f in hprc_per3_vntr hprc_per3_500kb; do tabix -f -p vcf "$OUT/$f.vcf.gz"; echo "$f: $(bcftools view -H "$OUT/$f.vcf.gz" | wc -l) records, $(bcftools query -l "$OUT/$f.vcf.gz" | wc -l) samples"; done
