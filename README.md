# The PER3 coding VNTR in a human pangenome – reproducible analysis

Data, code and results behind the manuscript *"Long-read haplotypes of the PER3 coding VNTR: chronotype associations
track sub-lineages of one four-unit allele, not repeat length"* (Chronomic Medicine Group, draft of 21 October 2026). Everything in the paper – each number, table and figure – is rebuilt from
the files in `data/` by one command, in about two minutes, without network access.

```bash
git clone https://github.com/comparativechrono/PER3-vntr-pangenome
cd PER3-vntr-pangenome
conda env create -f environment.yml      # bcftools, samtools/htslib, kalign, Python 3, R (+ susieR, GCTA: optional)
conda activate per3-vntr
bash run_all.sh                          # writes outputs/; ends by comparing every number with expected/
```

`run_all.sh` checks the input data against `CHECKSUMS_data.sha256`, runs the 17 steps in `scripts/` (a log of each
in `outputs/logs/`), records the software versions, and compares the 474 numbers of
`outputs/manuscript_numbers.tsv` with `expected/manuscript_numbers.tsv` (the run behind the manuscript). On
21 October 2026 two runs from a clean copy gave identical numbers.

## The analysis, step by step

| Step | Script | What it does |
| --- | --- | --- |
| 01 | `01_hprc_haplotypes.sh` | Rebuilds the 600 bp around the repeat (GRCh38 chr1:7,829,701-7,830,300) in each of the 462 haplotypes of the HPRC release 2 pangenome (Minigraph-Cactus v2.1, GRCh38 VCF of whole bubbles), with `bcftools consensus`; and the same 600 bp in the 1000 Genomes haplotypes in which the short-read panel called a deletion in the repeat. |
| 02 | `02_repeat_alleles.py` | Splits each haplotype's repeat into its 54-bp units (two independent counts: motif chain and length), names each unit after the five units of GRCh38 (R1-R5), and assigns lineages 4A, 4B, 5A, 5B; translates the units; reads the units of the 1000 Genomes haplotypes with a deletion call. |
| 03 | `03_graph_concordance.py` | Compares the units of every haplotype with the second graph of HPRC release 2 (PGGB). |
| 04 | `04_populations.R` | Allele and lineage frequencies by continental group (Wilson 95% CIs), heterogeneity (chi-squared), Hardy-Weinberg (exact test). |
| 05 | `05_ld_hprc.R` | LD (r², D′) between repeat length and every SNV within 250 kb in the 462 phased long-read haplotypes, by group; the best tag of each lineage. |
| 06 | `06_1000g.R` | The repeat in the 3,202 short-read genomes of 1000 Genomes through its tags: frequencies in 26 populations, lineages, validation against the long reads of the 222 people in both data sets, and what the short-read SV calls see, haplotype by haplotype. |
| 07 | `07_lineage_snps.R` | On which lineage each variant of the PER3 chronotype GWAS signal lies (HPRC and 1000 Genomes); the 4B sub-lineages 4B.1 and 4B.2. |
| 08 | `08_gwas_tag.R` | The tag (rs1773135-A = five units) in 13 GWAS of UK Biobank: chronotype, morning person, sleep duration, insomnia, napping and eight accelerometer measures; 95% CIs in SD units; detectable effects. |
| 09 | `09_conditional.R` | A joint model of markers of the lineages ('lineage model': 5A, 5B, 4B.1, 4B.2 and a rare branch of 4B.1, each against 4A, and five-unit haplotypes as a whole against 4A), stepwise conditional-joint analysis (COJO; this R implementation and GCTA 1.94.1), the single-SNP conditional test; three LD reference panels; the R code checked against GCTA. |
| 10 | `10_finemap_susie.R` | Fine-mapping with SuSiE-RSS, with each of the three LD panels (needs the R package susieR; skipped otherwise). |
| 11 | `11_gtex.R` | PER3 eQTL results (GTEx v8) for the tags and the lineage markers; the associations with plasma PER3 protein in the GWAS Catalog. |
| 12 | `12_primates.py` | Units of the repeat in 35 genome sequences of 28 primate species, human included (both haplotypes of six T2T apes), and in the PER3 proteins of 53 species (35 primates); the assembly behind each; a protein alignment (kalign). |
| 13 | `13_mukamel_comparison.R` | The same repeat in the phenome-wide scan of Mukamel et al. (2021), whose VNTR dosages were imputed in UK Biobank. |
| 14 | `14_implied_length_effects.R` | The association of repeat length that the lineage effects alone would produce in each 1000 Genomes population (an illustration). |
| 15 | `15_figures.R` | Figures 1-5 (PDF, 300-dpi PNG). |
| 16 | `16_numbers.R` | Every number quoted in the manuscript, with its source table: `outputs/manuscript_numbers.tsv`. |
| 17 | `17_manuscript_tables.R` | Tables 1-3 of the manuscript. |

`per3lib.py`, `lib.R` and `cojo.R` are shared code; `cojo.R` implements the conditional and joint analysis of Yang
et al. (2012) and reproduces GCTA's joint estimates to a relative difference below 1e-5 (`gcta_crosscheck.tsv`).

## Where each result of the paper is

| Manuscript | Output |
| --- | --- |
| Figure 1 | `outputs/figures/Figure1.*` (from `hprc_haplotype_alleles.tsv`, `hprc_lineages_by_group.tsv`, `kg_lineages_by_group.tsv`, `kg_freq_by_population.tsv`) |
| Figure 2 | `Figure2.*` (`hprc_ld_vntr_snvs.tsv`, `hprc_tags_top.tsv`) |
| Figure 3 | `Figure3.*` (`gwas_tag_results.tsv`, `lineage_model.tsv`) |
| Figure 4 | `Figure4.*` (1000 Genomes EUR haplotypes; chronotype GWAS) |
| Figure 5 | `Figure5.*` (`primate_units_dna.tsv`) |
| Tables 1-3 | `outputs/tables/Table1.tsv`, `Table2.tsv`, `Table3.tsv` |
| Table S1 | `data-build/data_sources.tsv` – the data sources, versions and the script that made each data file |
| Table S2 | `hprc_haplotype_alleles.tsv`, `graph_concordance.tsv` – every haplotype |
| Table S3 | `unit_catalogue.tsv`, `allele_structures.tsv` – the units and allele structures |
| Table S4 | `hprc_freq_by_population.tsv`, `hprc_lineages_by_group.tsv`, `hprc_tests.tsv`, `kg_freq_by_population.tsv`, `kg_lineages_by_population.tsv`, `sublineages_by_group.tsv` |
| Table S5 | `hprc_ld_vntr_snvs.tsv`, `hprc_tags_top.tsv`, `hprc_ld_extent.tsv`, `hprc_lineage_tags.tsv`, `kg_tag_ld.tsv` |
| Table S6 | `kg_validation.tsv`, `kg_validation_discordant.tsv`, `kg_sv_calls.tsv`, `kg_sv_by_group.tsv`, `kg_sv_carriers.tsv`, `kg_sv_alleles.tsv` |
| Table S7 | `mukamel_sleep_traits.tsv`, `mukamel_top_phenotypes.tsv` |
| Table S8 | `gwas_tag_results.tsv`, `gwas_tag_accelerometer_raw_units.tsv`, `gwas_lineage_marker_results.tsv` |
| Table S9 | `lineage_snps_hprc.tsv`, `lineage_snps_1kg.tsv` |
| Table S10 | `lineage_model.tsv`, `cojo_selected.tsv`, `cojo_tag_conditional.tsv`, `single_snp_conditional.tsv`, `cojo_reference_panels.tsv`, `gcta_crosscheck.tsv` |
| Table S11 | `susie_credible_sets.tsv`, `susie_tags.tsv`, `susie_diagnostics.tsv` |
| Table S12 | `gtex_summary.tsv`, `gtex_all.tsv`, `gwas_catalog_per3_protein.tsv` |
| Table S13 | `primate_units_dna.tsv`, `primate_units_protein.tsv`, `primate_units_by_species.tsv`, `primate_aligned_region.tsv`, `data/primates/assemblies.tsv` |
| Table S14 | `implied_length_effects.tsv` |
| Tables S1-S14 in Excel | `supplementary_tables/Table_S1.xlsx` to `Table_S14.xlsx` – each opens on an About sheet, then one sheet per file of the table as listed in the manuscript, values unchanged (made by `supplementary_tables/make_supplementary_xlsx.py`) |

## The data

All public; small extracts of large data sets (about 10 MB in all). `data/README_bundle.md` describes the files that
are shared with the student project; `data-build/` holds the scripts that made every file, from the sources below
(scripts 01-06, 08, 09 and 13-15 were re-run and gave identical files).

| Folder | Content | Source |
| --- | --- | --- |
| `reference/` | GRCh38 chr1:1-8,100,000; PER3 MANE exons; dbSNP build 157 rsIDs of chr1:7.58-8.08 Mb | NCBI; Ensembl REST |
| `hprc/` | HPRC release 2 (231 people, 462 haplotypes): Minigraph-Cactus v2.1 VCFs (repeat; 500 kb), PGGB VCF (80 kb); sample table | Human Pangenome Reference Consortium, AWS open data |
| `1000g/` | 1000 Genomes high-coverage phased panel: all 3,202 people across PER3 (70 kb); the 2,504 unrelated across 500 kb; sample table | IGSR / NYGC (2022) |
| `gwas/` | UK Biobank GWAS summary statistics for chr1:7.6-8.2 Mb (GRCh37) of chronotype, morning person, sleep duration, insomnia, napping and eight accelerometer measures; GWAS Catalog associations mapped to PER3 | GWAS Catalog (GCST007576, 007565, 007561, 007387, 011494, 007801-007808) |
| `gtex/` | GTEx v8 eQTL (dyneqtl) results for PER3 and the tags and lineage markers, 49 tissues | GTEx Portal API v2 |
| `mukamel2021/` | The PER3 VNTR's results for 787 phenotypes from Mukamel et al. 2021 | Broad Institute (authors' archive) |
| `primates/` | PER3 exon 19 ± 150 bp in 35 primate genome sequences; PER3 proteins of 53 species; species; assemblies | NCBI Datasets, E-utilities, BLAST+; T2T ape genomes |

Terms: the data keep the terms of their sources (open access for HPRC, 1000 Genomes, GWAS Catalog summary statistics,
GTEx eQTL summaries, NCBI and Ensembl; UK Biobank summary statistics as published by the original studies, which
should be cited). No individual-level data are included or were used.

## Things to know

- **Optional tools.** Without GCTA or the R package susieR the run still completes; the comparison at the end then
  lists the numbers that come from them (GCTA's stepwise selection, the fine-mapping) as missing. Everything else
  needs only bcftools, samtools/tabix, kalign, Python 3 and base R.
- **UK Biobank appears twice.** The GWAS used here and the phenome-wide scan of Mukamel et al. (2021) both analysed
  UK Biobank; the comparison (step 13) is between two ways of genotyping the repeat, not two samples.
- **LD reference.** The conditional and fine-mapping analyses take LD from 1000 Genomes Europeans (503; or 404, or
  190) rather than from UK Biobank; step 09 shows how much the answers depend on the panel.
- **The terms of the lineage model.** With markers of five units (rs1773135) and of 5A (rs6680585) in the same
  model, the five-unit marker's own term compares 5B with 4A. `lineage_model.tsv` therefore also gives 5A vs 4A,
  five-unit haplotypes as a whole vs 4A (5A and 5B weighted by their frequencies) and a 2-df test of 5A = 5B = 4A
  (rows of type `derived` and `test`), with standard errors from the covariance of the joint effects.
- **Stepwise selection is not unique.** COJO's choice among SNPs in near-perfect LD depends on ties and small LD
  differences (R and GCTA pick different but equivalent SNPs); the joint model of the lineage markers avoids this.
- **The accelerometer association** of repeat length with sleep midpoint (P = 9e-4, ~0.5 min per allele) is
  nominal, in ~85,000 people, and not mirrored by self-reported chronotype in 450,000. In the lineage model it is
  weaker for five-unit haplotypes as a whole (P = 0.053) and lies with the uncommon 5B allele; it needs replication.
- **Primate assemblies** differ in quality: `primate_units_dna.tsv` marks those without long reads, in which a
  tandem repeat can be collapsed.

## Licence

Code (`scripts/`, `data-build/`, `run_all.sh`): MIT licence. Figures, tables and text: CC BY 4.0. Data: the terms of
their sources (above). Please cite the manuscript and the data sources listed in it.
