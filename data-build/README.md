# How the project's data were made

The files in `data/examples/per3/` are small pieces of large public data sets. These scripts made them, on
21 October 2026, from the sources named in each. They are here so that the data can be checked and remade (for
instance when the sources are updated) – students do not need them.

They need a Linux or macOS computer with `bcftools`, `samtools`, `tabix`/`bgzip` (htslib, with network access),
`seqtk`, `curl`, Python 3 and, for the primates, BLAST+. Each script writes into the folder it is given:

```bash
mkdir out
bash 01_reference.sh out                         # reference/GRCh38_chr1_8Mb.fa.gz (+ .fai, .gzi)
python3 02_genes_gff3.py out                     # reference/genes_GRCh38.gff3, reference/PER3_exons.bed
bash 03_dbsnp.sh out                             # reference/dbsnp157_per3_500kb.vcf.gz
bash 04_hprc.sh out                              # hprc/hprc_per3_vntr.vcf.gz, hprc/hprc_per3_500kb.vcf.gz
python3 05_samples.py out/hprc_per3_vntr.vcf.gz out   # hprc/hprc_samples.tsv, 1000g/1kgp_samples.tsv
bash 06_1000g.sh out                             # 1000g/1kgp_per3_gene.vcf.gz
bash 07_gwas.sh out                              # gwas/*_Jones2019.tsv and the other GWAS (several GB of downloads)
python3 08_gwas_catalog.py out                   # gwas/gwas_catalog_PER3_associations.tsv
python3 09_gtex.py out                           # gtex/GTEx_v8_PER3_rs228701.tsv
python3 11_gtex_lineage_markers.py out           # gtex/GTEx_v8_PER3_lineage_markers.tsv
python3 12_primate_assemblies.py ../data/primates/PER3_repeat_exon_DNA.fa out   # primates/assemblies.tsv
bash 13_hprc_pggb.sh out                         # hprc/hprc_pggb_per3_7.78-7.86Mb.vcf.gz (streams ~100 MB)
bash 14_1000g_wide.sh out ../data/1000g/1kgp_samples.tsv   # 1000g/1kgp2504_chr1_7.58-8.08Mb.vcf.gz
bash 15_mukamel2021.sh out                       # mukamel2021/*
```

**Checked.** Scripts 01–06, 08, 09 and 13–15 were run again on 21 October 2026 and gave the same files as those in
`data/` (the same sequence, the same records, the same rows; the GWAS Catalog's rows may come in another order).
`07_gwas.sh` downloads several gigabytes, and 10–12 query services whose answers change over time; they are the
scripts that were used, run once. (In this package the data are in `data/`, not `data/examples/per3/` as in the
student project for which scripts 01–10 were first written.)

## Sources

The same list, with the script that made each file, is `data_sources.tsv` (Table S1 of the manuscript).

| Files | Source | Version / date |
| --- | --- | --- |
| `reference/GRCh38_chr1_8Mb.fa.gz` | NCBI E-utilities, NC_000001.11 bases 1–8,100,000 | GRCh38.p14 |
| `reference/genes_GRCh38.gff3`, `PER3_exons.bed` | Ensembl REST API, `overlap/region` (GFF3), protein-coding genes and their MANE Select transcripts | October 2026 |
| `reference/dbsnp157_per3_500kb.vcf.gz` | NCBI dbSNP build 157, `snp/archive/b157/VCF/GCF_000001405.40.gz` | build 157 |
| `hprc/*.vcf.gz` | Human Pangenome Reference Consortium, release 2, Minigraph-Cactus v2.1, GRCh38-based VCFs (`hprc-v2.1-mc-grch38.vcf.gz`, `.wave.vcf.gz`) on AWS S3 | release 2 (v2.1) |
| `hprc/hprc_samples.tsv` | HPRC release 2 sample table (GitHub: human-pangenomics/hprc_intermediate_assembly); 1000 Genomes population names | October 2026 |
| `1000g/1kgp_per3_gene.vcf.gz` | 1000 Genomes high-coverage phased panel (NYGC, 20220422) | 2022 |
| `1000g/1kgp_samples.tsv` | 1000 Genomes: the 3,202-sample pedigree file, the phase 3 panel of 2,504, the phase 3 table of populations | 2013–2022 |
| `gwas/*_Jones2019.tsv`, `*_Dashti2019.tsv`, `*_Lane2019.tsv`, `*_Dashti2021.tsv` | GWAS Catalog summary statistics FTP (GCST007576, 007565, 007561, 007387, 011494, 007801–007808) | as deposited (2019–2021) |
| `gwas/gwas_catalog_PER3_associations.tsv` | GWAS Catalog REST API v2, `associations?mapped_gene=PER3` | 21 October 2026 |
| `gtex/GTEx_v8_PER3_rs228701.tsv` | GTEx Portal API v2, `association/dyneqtl` for PER3 (ENSG00000049246.14) and chr1_7829744_C_G_b38 | GTEx v8 |
| `primates/*` | NCBI Datasets (orthologues of human PER3, gene 8863: gene regions and proteins); NCBI E-utilities and BLAST+ for the T2T apes | October 2026 |
| `gtex/GTEx_v8_PER3_lineage_markers.tsv` | GTEx Portal API v2, `association/dyneqtl` for PER3 and rs1773135, rs228663, rs61773390, rs228697 | GTEx v8 |
| `primates/assemblies.tsv` | NCBI E-utilities (elink nuccore→assembly) and NCBI Datasets v2 (genome dataset_report) | 21 October 2026 |
| `hprc/hprc_pggb_per3_7.78-7.86Mb.vcf.gz` | HPRC release 2 PGGB graph, `20250930_hprc25272.GRCh38.p98-k311.laced.vcf.gz` (AWS S3) | release 2 |
| `1000g/1kgp2504_chr1_7.58-8.08Mb.vcf.gz` | 1000 Genomes high-coverage phased panel (NYGC, 20220422), the 2,504 unrelated | 2022 |
| `mukamel2021/*` | `data.broadinstitute.org/lohlab/UKB_codingVNTR_sumstats/coding_vntr_stats.tar.gz` (Mukamel et al. 2021) | 2021-11-05 |

## The primates

These files were made in steps, not by one script:

1. **Orthologues.** NCBI Datasets: the orthologues of human *PER3* (gene ID 8863; `api.ncbi.nlm.nih.gov/datasets/v2/gene/id/8863/orthologs`, 510 species), and a data package of their gene regions (`gene.fna`) and proteins (`protein.faa`).
2. **The exon in the primates' genomes.** The human exon 19 (GRCh38 chr1:7,829,834–7,830,161) was searched against the gene regions of the primates with `blastn -task blastn` (E-value 10⁻²⁰); around the best hit of each, the exon and about 150 bases on each side were cut out, in the direction of the gene. These are the sequences of `PER3_repeat_exon_DNA.fa` named after a RefSeq chromosome (`reference_assembly`), with GRCh38's and CHM13's for human.
3. **The second haplotype of each ape.** `10_ape_haplotypes.py`: for the T2T assemblies of chimpanzee, bonobo, gorilla, both orangutans and the siamang (Yoo et al. 2025), the chromosome of the other haplotype that carries *PER3* is found by name, a 6-Mb window around the gene's expected position is fetched with E-utilities, and the human exon is found in it with BLAST. These are the sequences named `second_haplotype`.
4. **Proteins.** For each species chosen – the primates of step 2, the other primates and, for comparison, colugo, rodents, carnivores, cattle, opossum and chicken – the longest RefSeq protein of *PER3* in `protein.faa`: `PER3_proteins.fa`. `species.tsv` lists the species, their common names and groups, and the accession of each protein.
