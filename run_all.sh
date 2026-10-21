#!/usr/bin/env bash
# Reproduce every table, figure and number of the paper from the files in data/.
#   bash run_all.sh            (about two minutes on a laptop; no network needed)
# Needs bcftools, samtools and tabix (htslib), kalign (3.x), Python 3 (standard library only) and R (base; the R
# package susieR for the optional fine-mapping step). GCTA (gcta64) is optional: when it is on the PATH, or at
# tools/gcta64, the conditional analysis is also run with it as a check. See environment.yml.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p outputs/tables outputs/figures outputs/logs outputs/work

echo "== checking the input data"
sha256sum -c CHECKSUMS_data.sha256 --quiet && echo "   all $(wc -l < CHECKSUMS_data.sha256) data files as published"

run() {                      # run NAME COMMAND...: log to outputs/logs/NAME.log, show the last lines
  local name=$1; shift
  echo "== $name"
  if ! "$@" > "outputs/logs/$name.log" 2>&1; then
    echo "   FAILED - see outputs/logs/$name.log"; tail -n 20 "outputs/logs/$name.log"; exit 1
  fi
  grep -E '^[0-9]{2}: ' "outputs/logs/$name.log" | head -n 4 | sed 's/^/   /' || true
}
run 01_hprc_haplotypes        bash scripts/01_hprc_haplotypes.sh
run 02_repeat_alleles         python3 scripts/02_repeat_alleles.py
run 03_graph_concordance      python3 scripts/03_graph_concordance.py
run 04_populations            Rscript scripts/04_populations.R
run 05_ld_hprc                Rscript scripts/05_ld_hprc.R
run 06_1000g                  Rscript scripts/06_1000g.R
run 07_lineage_snps           Rscript scripts/07_lineage_snps.R
run 08_gwas_tag               Rscript scripts/08_gwas_tag.R
run 09_conditional            Rscript scripts/09_conditional.R
run 10_finemap_susie          Rscript scripts/10_finemap_susie.R
run 11_gtex                   Rscript scripts/11_gtex.R
run 12_primates               python3 scripts/12_primates.py
run 13_mukamel_comparison     Rscript scripts/13_mukamel_comparison.R
run 14_implied_length_effects Rscript scripts/14_implied_length_effects.R
run 15_figures                Rscript scripts/15_figures.R
run 16_numbers                Rscript scripts/16_numbers.R
run 17_manuscript_tables      Rscript scripts/17_manuscript_tables.R

echo "== software versions"
set +o pipefail              # 'head' and 'grep -m1' close their pipes early
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "uname: $(uname -srm)"
  bcftools --version | head -n 2
  samtools --version | head -n 1
  tabix --version 2>&1 | head -n 1
  echo "kalign: $(kalign --version 2>&1 | head -n 1)"
  python3 --version
  R --version | head -n 1
  if command -v gcta64 > /dev/null; then gcta64 2>&1 | grep -m1 version; elif [ -x tools/gcta64 ]; then tools/gcta64 2>&1 | grep -m1 version; else echo "GCTA: not found (check skipped)"; fi
} > outputs/logs/versions.txt 2>&1
set -o pipefail
Rscript -e 'suppressWarnings(suppressMessages(requireNamespace("susieR", quietly = TRUE))); sessionInfo()' > outputs/logs/R_sessionInfo.txt 2>&1
sed 's/^/   /' outputs/logs/versions.txt

echo "== comparing the numbers with those of the paper (expected/manuscript_numbers.tsv)"
if cmp -s <(cut -f1,2 outputs/manuscript_numbers.tsv) <(cut -f1,2 expected/manuscript_numbers.tsv); then
  echo "   identical: all $(($(wc -l < outputs/manuscript_numbers.tsv) - 1)) numbers"
else
  echo "   DIFFERENT - lines that differ (expected first):"
  diff <(cut -f1,2 expected/manuscript_numbers.tsv) <(cut -f1,2 outputs/manuscript_numbers.tsv) | head -n 40 || true
fi
(cd outputs && find tables figures -type f | sort | xargs sha256sum > CHECKSUMS_outputs.sha256)
echo "== done: outputs/tables, outputs/figures, outputs/manuscript_numbers.tsv, outputs/logs"
