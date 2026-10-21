# The data of the PER3 VNTR project

Small pieces of large public data sets – the parts around the gene *PER3* on chromosome 1 – so that the whole
project runs in a web browser. Read-only: copy a file into your home folder to change it. The guide's page
*The project data* says the same as this file. How each file was made is in `data-build/` of the site's repository.

Positions are GRCh38 unless a file says otherwise; the GWAS files are GRCh37.

## `reference/` – the human genome around *PER3*

| File | What it is |
| --- | --- |
| `GRCh38_chr1_8Mb.fa.gz` (+ `.fai`, `.gzi`) | The first 8,100,000 bases of chromosome 1 of GRCh38 (NCBI `NC_000001.11`), compressed with bgzip and indexed. The sequence is called `chr1`, and its positions are those of the whole chromosome. |
| `genes_GRCh38.gff3` | The seven protein-coding genes of chr1:7.4–8.1 Mb – *CAMTA1*, *VAMP3*, *PER3*, *UTS2*, *TNFRSF9*, *PARK7*, *ERRFI1* – with the exons and CDS of their MANE Select transcript (Ensembl REST API, October 2026). |
| `PER3_exons.bed` | The 22 exons of *PER3*'s MANE Select transcript, ENST00000377532, numbered (BED: counted from 0). |
| `dbsnp157_per3_500kb.vcf.gz` (+ `.tbi`) | The 209,029 variants of dbSNP build 157 between chr1:7,580,000 and 8,080,000: positions, alleles and rsIDs. |

## `hprc/` – 231 long-read genomes

| File | What it is |
| --- | --- |
| `hprc_per3_vntr.vcf.gz` (+ `.tbi`) | The variants of the 462 haplotypes of 231 people at chr1:7,829,001–7,831,000 – the repeat and 2 kb around it – from the pangenome of the Human Pangenome Reference Consortium (release 2, Minigraph-Cactus graph v2.1, GRCh38-based VCF before decomposition: one record per bubble of the graph, with whole alleles). Phased genotypes. 26 records. |
| `hprc_per3_500kb.vcf.gz` (+ `.tbi`) | The same 231 people, chr1:7,580,000–8,080,000, from the decomposed VCF of the same graph (`wave`): each record an SNV, an indel or a larger variant. 7,441 records. |
| `hprc_samples.tsv` | The 231 people: population (1000 Genomes codes, with the 1000 Genomes Project's names for them), continental group (`AFR`, `AMR`, `EAS`, `EUR`, `SAS`), sex, and whether the person is also in the high-coverage 1000 Genomes data. HG002 and HG005 come from Genome in a Bottle. |

The reference haplotype of the pangenome, CHM13, was left out; so were the INFO fields `AF`, `NS` and `AT` and the
graph's names of the variants.

## `1000g/` – 3,202 short-read genomes

| File | What it is |
| --- | --- |
| `1kgp_per3_gene.vcf.gz` (+ `.tbi`) | The phased SNVs, indels and structural variants of the 3,202 people of the high-coverage 1000 Genomes Project across *PER3*, chr1:7,780,000–7,850,000 (New York Genome Center, phased panel of 2022), without their INFO fields. 1,775 records. IDs are `chromosome:position:REF:ALT`. |
| `1kgp_samples.tsv` | Their population, its name, continental group and sex; their parents if they are in the data; and whether they are one of the 2,504 unrelated people of the original project (`unrelated_2504`). |

## `gwas/` – genome-wide association studies of sleep (UK Biobank, European ancestry)

The results for every variant between chr1:7,600,000 and 8,200,000, in **GRCh37** positions, from the GWAS
Catalog's summary statistics. The effect belongs to `ALLELE1` (`effect_allele` in the napping file).

| File | Trait | People | Study |
| --- | --- | --- | --- |
| `chronotype_Jones2019.tsv` | chronotype score (higher: more of a morning person) | 449,734 | Jones et al. 2019, GCST007576 |
| `morning_person_Jones2019.tsv` | morning person or not (LOGOR) | 403,195 | Jones et al. 2019, GCST007565 |
| `sleep_duration_Dashti2019.tsv` | usual hours of sleep | 446,118 | Dashti et al. 2019, GCST007561 |
| `insomnia_Lane2019.tsv` | insomnia symptoms | 453,379 | Lane et al. 2019, GCST007387 |
| `daytime_nap_Dashti2021.tsv` | daytime napping | 452,633 | Dashti et al. 2021, GCST011494 |
| `accelerometer_Jones2019.tsv` | eight measures of sleep from a wrist accelerometer worn for a week (the midpoint of sleep, L5 and M10 timing, duration, efficiency, number of episodes, variability of duration, daytime inactivity), each as `_RAW` (untransformed: hours, a proportion for efficiency) and `_RAW_SIN` (rank-normalised: standard deviations); missing values are `-` | about 85,000 | Jones et al. 2019, GCST007801–007808 |
| `gwas_catalog_PER3_associations.tsv` | the 28 associations the GWAS Catalog maps to *PER3* (October 2026) | | GWAS Catalog |

## `gtex/` – gene expression

| File | What it is |
| --- | --- |
| `GTEx_v8_PER3_rs228701.tsv` | For 49 tissues of GTEx v8 (13 of them parts of the brain): the effect of rs228701 on the amount of *PER3* RNA – the normalised effect size (NES) of its alternative allele, G – its p-value, and the number of donors of each genotype (GTEx Portal API, `dyneqtl`). |

## `primates/` – *PER3* in other species

| File | What it is |
| --- | --- |
| `PER3_repeat_exon_DNA.fa` | The repeat's exon and about 150 bases on each side, in the direction of the gene, from 35 genomes of 28 primate species: NCBI RefSeq reference assemblies, and for the apes both haplotypes of their T2T genomes (Yoo et al. 2025). Found with BLAST of the human exon. |
| `PER3_proteins.fa` | The PER3 protein (the longest RefSeq isoform) of 53 species: primates, and colugo (a fragment), rodents, carnivores, cattle, opossum and chicken. |
| `species.tsv` | The species, their common names and groups, how many DNA sequences each has here, and the accession of its protein. |

## To cite

- **HPRC pangenome**: Liao et al. 2023, *Nature* 617: 312–324; Hickey et al. 2024, *Nature Biotechnology* 42: 663–673; Human Pangenome Reference Consortium, Release 2.
- **1000 Genomes, high coverage**: Byrska-Bishop et al. 2022, *Cell* 185: 3426–3440.
- **GWAS**: Jones et al. 2019 (*Nature Communications* 10: 343 and 10: 1585), Dashti et al. 2019 (*Nature Communications* 10: 1100) and 2021 (12: 900), Lane et al. 2019 (*Nature Genetics* 51: 387–393); the GWAS Catalog: Cerezo et al. 2025, *Nucleic Acids Research* 53: D998–D1005. This research used data of UK Biobank, through those studies.
- **GTEx**: GTEx Consortium 2020, *Science* 369: 1318–1330.
- **Ape genomes**: Yoo et al. 2025, *Nature* 641: 401–418. **Other species**: NCBI RefSeq (O'Leary et al. 2016, *Nucleic Acids Research* 44: D733–D745).
- **GRCh38**: Genome Reference Consortium; **dbSNP**: Sherry et al. 2001, *Nucleic Acids Research* 29: 308–311; **Ensembl**: Dyer et al. 2025, *Nucleic Acids Research* 53: D948–D957.
