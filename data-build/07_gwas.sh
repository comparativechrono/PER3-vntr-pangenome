#!/bin/bash
# gwas/*.tsv: the summary statistics of six GWAS of sleep in UK Biobank (GWAS Catalog FTP), cut to the variants of
# chr1:7,600,000-8,200,000 (GRCh37, as published). The full files are large (0.5-2.5 GB each, 2.5 GB for the
# accelerometer one): this takes a while, and the full files are deleted once cut.
# Needs: curl, gzip, unzip, awk.   Run:  bash 07_gwas.sh OUTDIR
set -u
OUT=${1:-.}
cd "$OUT"
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000
win() { awk -F'\t' 'NR==1 || ($2==1 && $3>=7600000 && $3<=8200000)'; }     # the header, and column 2 = 1, column 3 in the window
get() { curl -sS --retry 3 -o "$2" "$1" && echo "downloaded $2"; }

# chronotype (Jones et al. 2019, GCST007576) and morning person (GCST007565)
get "$B/GCST007576/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz" chrono.gz
zcat chrono.gz | win > chronotype_Jones2019.tsv && rm chrono.gz
get "$B/GCST007565/morning_person_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3_logORs.txt.gz" morning.gz
zcat morning.gz | win > morning_person_Jones2019.tsv && rm morning.gz
# sleep duration (Dashti et al. 2019, GCST007561) and insomnia (Lane et al. 2019, GCST007387): space-separated
get "$B/GCST007561/Dashti_30846698_sleepdurationsumstats.txt.zip" sleepdur.zip
unzip -p sleepdur.zip | tr ' ' '\t' | win > sleep_duration_Dashti2019.tsv && rm sleepdur.zip
get "$B/GCST007387/Lane_30804566_fullUKBB_Insomnia_summary_stats.zip" insomnia.zip
unzip -p insomnia.zip | tr ' ' '\t' | win > insomnia_Lane2019.tsv && rm insomnia.zip
# daytime napping (Dashti et al. 2021, GCST011494; GRCh37 file)
get "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST011001-GCST012000/GCST011494/GCST011494_buildGRCh37.tsv" nap.tsv
win < nap.tsv > daytime_nap_Dashti2021.tsv && rm nap.tsv
# accelerometer-derived sleep (Jones et al. 2019, GCST007801-007808, all eight measures in one file)
get "$B/GCST007805/accel_GWAS_all_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.zip.gz" accel.zip.gz
gunzip -c accel.zip.gz > accel.zip && unzip -p accel.zip | win > accelerometer_Jones2019.tsv && rm accel.zip accel.zip.gz
wc -l *.tsv
