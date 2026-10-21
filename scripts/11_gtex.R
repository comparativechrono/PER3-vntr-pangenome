# 11. PER3 expression (GTEx v8 eQTL, every tissue with eQTL results): the repeat's tags (rs228701-G and
#     rs1773135-G = four units) and the markers of 4B (rs228663-A), 4B.1 (rs61773390-T) and 4B.2 (rs228697-G).
#     NES is the normalised effect of the ALT allele; Bonferroni over tissues: 0.05 / number of tissues.
#     Also: the associations with plasma PER3 protein levels in the GWAS Catalog extract (data/gwas), and the
#     lineages on which their variants lie (1000 Genomes, step 07).
# Outputs (outputs/tables): gtex_all.tsv, gtex_summary.tsv, gwas_catalog_per3_protein.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

a <- read_tsv(file.path(DATA, "gtex", "GTEx_v8_PER3_rs228701.tsv"))
a <- data.frame(rsid = "rs228701", variant_id = "chr1_7829744_C_G_b38", tissue = a$tissue, n_hom_ref = a$n_CC, n_het = a$n_CG,
                n_hom_alt = a$n_GG, NES = a$NES, p_value = a$p_value)
b <- read_tsv(file.path(DATA, "gtex", "GTEx_v8_PER3_lineage_markers.tsv"))
x <- rbind(a, b)
x$marks <- c(rs228701 = "four units (ALT G)", rs1773135 = "four units (ALT G)", rs228663 = "4B (ALT A)",
             rs61773390 = "4B.1 (ALT T)", rs228697 = "4B.2 (ALT G)")[x$rsid]
x$brain <- grepl("^Brain", x$tissue)
write_tsv(x, "gtex_all.tsv")
s <- do.call(rbind, lapply(split(x, x$rsid), function(d) {
  nt <- nrow(d); thr <- 0.05 / nt; i <- which.min(d$p_value); bb <- d[d$brain, ]; j <- which.min(bb$p_value)
  data.frame(rsid = d$rsid[1], marks = d$marks[1], tissues = nt, bonferroni = signif(thr, 2),
             tissues_p05 = sum(d$p_value < 0.05), tissues_bonferroni = sum(d$p_value < thr),
             top_tissue = d$tissue[i], top_NES = d$NES[i], top_p = d$p_value[i],
             brain_tissues = nrow(bb), brain_p05 = sum(bb$p_value < 0.05), top_brain = bb$tissue[j], top_brain_NES = bb$NES[j],
             top_brain_p = bb$p_value[j],
             hypothalamus_NES = d$NES[d$tissue == "Brain_Hypothalamus"], hypothalamus_p = d$p_value[d$tissue == "Brain_Hypothalamus"])
}))
s <- s[match(c("rs228701", "rs1773135", "rs228663", "rs61773390", "rs228697"), s$rsid), ]
write_tsv(s, "gtex_summary.tsv")
cat("11: GTEx v8, PER3:\n"); print(s, row.names = FALSE)
cat("11: all results with P below the Bonferroni threshold:\n")
print(x[x$p_value < 0.05 / 49, c("rsid", "marks", "tissue", "NES", "p_value")], row.names = FALSE)

# plasma PER3 protein levels in the GWAS Catalog (associations mapped to PER3), with the lineages of their variants
gc <- read_tsv(file.path(DATA, "gwas", "gwas_catalog_PER3_associations.tsv"))
pq <- gc[grepl("period circadian protein homolog 3", gc$trait, ignore.case = TRUE), ]
pq <- pq[order(as.numeric(pq$p_value)), ]
ls1 <- read_tsv(file.path(TAB, "lineage_snps_1kg.tsv"))
pick <- function(r, d, col) { v <- ls1[[col]][ls1$rsid == r & grepl(d, ls1$data)]; if (length(v)) v[1] else NA }
pq$carriers_1000G_all <- sapply(pq$rsid, pick, d = "All", col = "alt_haplotypes")
pq$on_4B_1000G_all <- sapply(pq$rsid, pick, d = "All", col = "alt_on_4B")
write_tsv(pq, "gwas_catalog_per3_protein.tsv")
cat("11: plasma PER3 protein associations in the GWAS Catalog:\n"); print(pq, row.names = FALSE)
