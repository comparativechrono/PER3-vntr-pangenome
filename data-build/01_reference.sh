#!/bin/bash
# reference/GRCh38_chr1_8Mb.fa.gz (+ .fai, .gzi): the first 8,100,000 bases of GRCh38 chromosome 1 (NCBI RefSeq
# NC_000001.11), named chr1, 60 letters a line, compressed with bgzip and indexed.
# Needs: curl, seqtk, bgzip and samtools.   Run:  bash 01_reference.sh OUTDIR
set -euo pipefail
OUT=${1:-.}
curl -sS --retry 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000001.11&rettype=fasta&seq_start=1&seq_stop=8100000" > "$OUT/chr1_part.fa"
{
  echo ">chr1 GRCh38 chromosome 1, bases 1-8,100,000 (NCBI NC_000001.11); positions are those of the whole chromosome"
  seqtk seq -l 60 "$OUT/chr1_part.fa" | grep -v '^>' | tr 'acgtn' 'ACGTN'
} > "$OUT/GRCh38_chr1_8Mb.fa"
rm "$OUT/chr1_part.fa"
bgzip -f -@ 2 -i "$OUT/GRCh38_chr1_8Mb.fa"        # -i writes the .gzi index too
samtools faidx "$OUT/GRCh38_chr1_8Mb.fa.gz"
cut -f 1,2 "$OUT/GRCh38_chr1_8Mb.fa.gz.fai"        # chr1  8100000
