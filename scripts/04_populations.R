# 04. Who carries which allele, in the 231 HPRC people (long-read haplotypes, step 02).
# Outputs (outputs/tables): hprc_freq_by_group.tsv, hprc_lineages_by_group.tsv, hprc_genotypes.tsv,
#                           hprc_freq_by_population.tsv, hprc_tests.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

a <- read_tsv(file.path(TAB, "hprc_haplotype_alleles.tsv"))
stopifnot(nrow(a) == 462)
a$five <- as.integer(a$units == 5)

# allele frequency by continental group, with Wilson 95% CIs
grp <- do.call(rbind, lapply(c(GROUPS, "All"), function(g) {
  x <- if (g == "All") a else a[a$superpopulation == g, ]
  ci <- wilson(sum(x$five), nrow(x))
  data.frame(group = g, people = length(unique(x$sample)), haplotypes = nrow(x), four = sum(x$units == 4),
             five = sum(x$five), freq_five = round(mean(x$five), 3), ci_lower = round(ci[1], 3), ci_upper = round(ci[2], 3))
}))
write_tsv(grp, "hprc_freq_by_group.tsv")

# heterogeneity among the five groups
tab <- table(factor(a$superpopulation, GROUPS), a$units)
chi <- chisq.test(tab)
# genotypes and Hardy-Weinberg
gt <- tapply(a$units, a$sample, function(x) paste(sort(x), collapse = "/"))
gtab <- table(factor(gt, c("4/4", "4/5", "5/5")))
hwe_p <- hwe_exact(gtab[["4/4"]], gtab[["4/5"]], gtab[["5/5"]])
hwe_group <- sapply(GROUPS, function(g) {
  s <- unique(a$sample[a$superpopulation == g]); t <- table(factor(gt[s], c("4/4", "4/5", "5/5")))
  hwe_exact(t[["4/4"]], t[["4/5"]], t[["5/5"]])
})
write_tsv(data.frame(genotype = names(gtab), people = as.integer(gtab)), "hprc_genotypes.tsv")
tests <- data.frame(test = c("chi-squared, five groups x two alleles", "Hardy-Weinberg exact, all 231",
                             paste0("Hardy-Weinberg exact, ", GROUPS)),
                    statistic = c(round(unname(chi$statistic), 2), rep(NA, 6)),
                    df = c(unname(chi$parameter), rep(NA, 6)),
                    p = signif(c(chi$p.value, hwe_p, hwe_group), 3))
write_tsv(tests, "hprc_tests.tsv")

# lineages by group
lin <- as.data.frame.matrix(table(factor(a$superpopulation, GROUPS), factor(a$lineage, names(LINEAGE_COL))))
lin <- rbind(lin, All = colSums(lin))
lin$haplotypes <- rowSums(lin)
lin$share_4B_of_four <- round(lin[["4B"]] / (lin[["4A"]] + lin[["4B"]]), 3)
lin$share_5B_of_five <- round(lin[["5B"]] / (lin[["5A"]] + lin[["5B"]]), 3)
write_tsv(cbind(group = rownames(lin), lin), "hprc_lineages_by_group.tsv")

# by population
pp <- aggregate(cbind(haplotypes = 1, five = a$five) ~ population + superpopulation, data = a, FUN = sum)
pp$freq_five <- round(pp$five / pp$haplotypes, 3)
write_tsv(pp[order(pp$superpopulation, pp$population), ], "hprc_freq_by_population.tsv")

cat("04: five-unit allele frequency by group (HPRC):\n"); print(grp, row.names = FALSE)
cat(sprintf("04: heterogeneity chi2 = %.1f, df = %d, P = %.2g; genotypes 4/4 %d, 4/5 %d, 5/5 %d; HWE P = %.2f\n",
            chi$statistic, chi$parameter, chi$p.value, gtab[["4/4"]], gtab[["4/5"]], gtab[["5/5"]], hwe_p))
cat("04: lineages by group:\n"); print(lin)
