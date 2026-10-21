# 07. On which repeat lineages do the PER3 variants of the chronotype GWAS sit?
#     Long-read haplotypes (HPRC, lineages read directly) and 1000 Genomes (lineages from tags, step 06).
# Outputs (outputs/tables): lineage_snps_hprc.tsv, lineage_snps_1kg.tsv, lineage_snps_r2.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

snps <- data.frame(
  rsid = c("rs150812083", "rs10462020", "rs61773390", "rs228697", "rs228726", "rs228654", "rs2640909", "rs228663",
           "rs2859387", "rs12023156", "rs6680585", "rs7545893", "rs2087947", "rs148515378", "rs10462021"),
  pos = c(7809893, 7820623, 7824465, 7827519, 7788883, 7837168, 7830057, 7832049, 7827188, 7830004, 7805003, 7796207,
          7791033, 7814731, 7837073),
  ref = c("C", "T", "G", "C", "G", "G", "T", "T", "G", "A", "G", "C", "C", "C", "A"),
  alt = c("G", "G", "T", "G", "A", "A", "C", "A", "A", "G", "A", "A", "T", "A", "G"),
  note = c("rare; chronotype lead variant (Jones 2019)", "p.Val647Gly (V647G)", "chronotype lead variant (Jones 2019)",
           "p.Pro864Ala", "chronotype lead variant (Jones 2019)", "", "p.Met1037Thr, inside unit R4",
           "best tag of 4B", "best tag of 4A", "best tag of 5A, inside unit R3", "tag of 5A outside the repeat",
           "lead variant of Lane 2016", "proxy of the VNTR used by Webb 2022", "PER3 protein QTL (Loya 2025)",
           "p.His1158Arg (H1158R)"),
  stringsAsFactors = FALSE)

lineages <- c("4A", "4B", "5A", "5B", "other")
tabulate_snps <- function(h, info, lineage, label) {
  do.call(rbind, lapply(seq_len(nrow(snps)), function(i) {
    j <- which(info$pos == snps$pos[i] & info$ref == snps$ref[i] & info$alt == snps$alt[i])
    if (length(j) != 1) return(NULL)
    x <- h[j, ]
    cnt <- sapply(lineages, function(l) sum(x[lineage == l] == 1, na.rm = TRUE))
    tot <- sapply(lineages, function(l) sum(lineage == l & !is.na(x)))
    data.frame(data = label, rsid = snps$rsid[i], pos = snps$pos[i], alt = snps$alt[i], note = snps$note[i],
               alt_haplotypes = sum(x == 1, na.rm = TRUE), haplotypes = sum(!is.na(x)),
               t(setNames(paste0(cnt, "/", tot), paste0("alt_on_", lineages))),
               r2_with_five = NA,   # filled in below, from the number of units of each haplotype
               r2_with_4A = round(r2(x, as.integer(lineage == "4A")), 3),
               r2_with_4B = round(r2(x, as.integer(lineage == "4B")), 3),
               r2_with_5A = round(r2(x, as.integer(lineage == "5A")), 3),
               r2_with_5B = round(r2(x, as.integer(lineage == "5B")), 3), check.names = FALSE)
  }))
}

# HPRC: lineages from the long reads; SNVs from the decomposed VCF of the same graph
a <- read_tsv(file.path(TAB, "hprc_haplotype_alleles.tsv"))
vcf_h <- file.path(DATA, "hprc", "hprc_per3_500kb.vcf.gz")
s_h <- system2("bcftools", c("query", "-l", vcf_h), stdout = TRUE)
mh <- haplotype_matrix(vcf_h, "chr1:7780000-7850000", s_h)
a <- a[match(colnames(mh$h), a$haplotype), ]
five_h <- as.integer(a$units == 5)
th <- tabulate_snps(mh$h, mh$info, a$lineage, "HPRC (462 haplotypes)")
th$r2_with_five <- sapply(seq_len(nrow(th)), function(i) {
  j <- which(mh$info$pos == th$pos[i] & mh$info$alt == th$alt[i]); round(r2(mh$h[j, ], five_h), 3) })
write_tsv(th, "lineage_snps_hprc.tsv")

# 1000 Genomes: lineages from the tags, as in step 06
vcf_k <- file.path(DATA, "1000g", "1kgp_per3_gene.vcf.gz")
smp <- read_tsv(file.path(DATA, "1000g", "1kgp_samples.tsv"))
s_k <- system2("bcftools", c("query", "-l", vcf_k), stdout = TRUE)
smp <- smp[match(s_k, smp$sample), ]
mk <- haplotype_matrix(vcf_k, "chr1:7780000-7850000", s_k)
g1 <- function(pos, ref, alt) mk$h[which(mk$info$pos == pos & mk$info$ref == ref & mk$info$alt == alt), ]
five_k <- as.integer(g1(7830204, "A", "G") == 0)
lin_k <- ifelse(five_k == 1, ifelse(g1(7830004, "A", "G") == 1, "5A", "5B"), ifelse(g1(7832049, "T", "A") == 1, "4B", "4A"))
grp_k <- rep(smp$superpopulation, each = 2); unrel <- rep(smp$unrelated_2504 == "yes", each = 2)
sub_k <- ifelse(lin_k == "4B", ifelse(g1(7824465, "G", "T") == 1, "4B.1", ifelse(g1(7827519, "C", "G") == 1, "4B.2", "4B.0")), lin_k)
tk <- do.call(rbind, lapply(c("EUR", "All"), function(g) {
  k <- unrel & (g == "All" | grp_k == g)
  t <- tabulate_snps(mk$h[, k], mk$info, lin_k[k], paste0("1000G ", g, " (", sum(k), " haplotypes)"))
  t$r2_with_five <- sapply(seq_len(nrow(t)), function(i) {
    j <- which(mk$info$pos == t$pos[i] & mk$info$alt == t$alt[i]); round(r2(mk$h[j, k], five_k[k]), 3) })
  for (sl_ in c("4B.1", "4B.2", "4B.0")) t[[paste0("alt_on_", sl_)]] <- sapply(seq_len(nrow(t)), function(i) {
    j <- which(mk$info$pos == t$pos[i] & mk$info$alt == t$alt[i]); x <- mk$h[j, k]
    paste0(sum(x[sub_k[k] == sl_] == 1, na.rm = TRUE), "/", sum(sub_k[k] == sl_)) })
  t
}))
write_tsv(tk, "lineage_snps_1kg.tsv")

# 4B sub-lineages (4B.1: rs61773390-T; 4B.2: rs228697-G; 4B.0: neither) by group, HPRC (long reads) and 1000G
sub_h <- function(h, info, lineage, grp, label) {
  g <- function(pos, alt) h[which(info$pos == pos & info$alt == alt), ]
  s1 <- g(7824465, "T"); s2 <- g(7827519, "G")
  sub <- ifelse(lineage == "4B", ifelse(s1 == 1, ifelse(s2 == 1, "4B.1+4B.2", "4B.1"), ifelse(s2 == 1, "4B.2", "4B.0")), lineage)
  do.call(rbind, lapply(c(GROUPS, "All"), function(gg) { k <- if (gg == "All") rep(TRUE, length(sub)) else grp == gg
    t <- table(factor(sub[k], c("4A", "4B.0", "4B.1", "4B.2", "4B.1+4B.2", "5A", "5B", "other")))
    data.frame(data = label, group = gg, haplotypes = sum(k), t(as.matrix(t)), check.names = FALSE) }))
}
sl <- rbind(sub_h(mh$h, mh$info, a$lineage, a$superpopulation, "HPRC"),
            sub_h(mk$h[, unrel], mk$info, lin_k[unrel], grp_k[unrel], "1000G (2,504 unrelated)"))
write_tsv(sl, "sublineages_by_group.tsv")
cat("07: 4B sub-lineages:\n"); print(sl, row.names = FALSE)
show <- c("data", "rsid", "alt", "alt_haplotypes", "haplotypes", paste0("alt_on_", lineages), "r2_with_five", "r2_with_4B", "r2_with_5A")
cat("07: HPRC\n"); print(th[, show], row.names = FALSE)
cat("07: 1000 Genomes EUR (unrelated)\n"); print(tk[grepl("EUR", tk$data), show], row.names = FALSE)
