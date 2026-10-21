# 17. Tables 1-3 of the manuscript, formatted (outputs/tables/Table1.tsv, Table2.tsv, Table3.tsv)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))
t <- function(f) read_tsv(file.path(TAB, f))
ci <- function(f, lo, hi) sprintf("%.3f (%.3f-%.3f)", f, lo, hi)

# Table 1: the repeat by continental group: long reads (HPRC) and tags (1000 Genomes, 2,504 unrelated)
hf <- t("hprc_freq_by_group.tsv"); kf <- t("kg_freq_by_group.tsv"); kl <- t("kg_lineages_by_group.tsv"); sl <- t("sublineages_by_group.tsv")
sl <- sl[grepl("1000G", sl$data), ]
tab1 <- data.frame(Group = c(GROUP_NAMES[GROUPS], All = "All")[hf$group],
                   `HPRC people` = hf$people, `HPRC five-unit frequency (95% CI)` = ci(hf$freq_five, hf$ci_lower, hf$ci_upper),
                   `1000G haplotypes` = kf$haplotypes[match(hf$group, kf$group)],
                   `1000G five-unit frequency (95% CI)` = with(kf[match(hf$group, kf$group), ], ci(freq_five, ci_lower, ci_upper)),
                   `4A` = sprintf("%.2f", kl[["4A"]][match(hf$group, kl$group)]), `4B` = sprintf("%.2f", kl[["4B"]][match(hf$group, kl$group)]),
                   `5A` = sprintf("%.2f", kl[["5A"]][match(hf$group, kl$group)]), `5B` = sprintf("%.2f", kl[["5B"]][match(hf$group, kl$group)]),
                   `4B.1 : 4B.2` = with(sl[match(hf$group, sl$group), ], paste0(`4B.1`, " : ", `4B.2`)), check.names = FALSE)
write_tsv(tab1, "Table1.tsv")

# Table 2: the five-unit allele (rs1773135-A) and 13 traits
tr <- t("gwas_tag_results.tsv"); tr <- tr[tr$rsid == "rs1773135", ]; tt <- t("gwas_traits.tsv")
tab2 <- data.frame(Trait = tr$label, N = formatC(tt$N[match(tr$trait, tt$trait)], big.mark = ",", format = "d"), Units = tr$scale,
                   `Beta (SE)` = sprintf("%.4f (%.4f)", tr$beta, tr$se),
                   `SD per allele (95% CI)` = sprintf("%.4f (%.4f to %.4f)", tr$beta_sd, tr$ci_lower_sd, tr$ci_upper_sd),
                   P = fmt_p(tr$p), `Detectable with 80% power at P < 0.05 (SD)` = sprintf("%.4f", tr$mde80_p05_sd),
                   Study = tr$study, check.names = FALSE)
write_tsv(tab2, "Table2.tsv")

# Table 3: the lineage model (joint effects; reference = 4A); EUR LD panel, with the range over the three panels
lm_ <- t("lineage_model.tsv")
rows3 <- data.frame(contrast = c("five units (5A and 5B by frequency) vs 4A", "5A vs 4A", "5B vs 4A", "4B.1 vs 4A", "4B.2 vs 4A",
                                 "rare branch of 4B.1 (vs rest of 4B.1)"),
                    label = c("Five units (5A and 5B) vs 4A", "5A vs 4A", "5B vs 4A", "4B.1 vs 4A", "4B.2 vs 4A",
                              "Rare branch of 4B.1 vs other 4B.1"), stringsAsFactors = FALSE)
traits3 <- c(chronotype = "Chronotype", morning = "Morning person", acc_midpoint = "Sleep midpoint (accelerometer)")
tab3 <- do.call(rbind, lapply(names(traits3), function(tr) do.call(rbind, lapply(seq_len(nrow(rows3)), function(i) {
  x <- lm_[lm_$trait == tr & lm_$contrast == rows3$contrast[i], ]; e <- x[x$panel == "EUR", ]
  data.frame(Trait = traits3[[tr]], `Haplotypes compared` = rows3$label[i],
             `Joint effect (SE)` = sprintf("%.4f (%.4f)", e$bJ, e$bJ_se), `Joint P` = fmt_p(e$pJ),
             `Joint P, three LD panels` = paste(fmt_p(range(x$pJ)), collapse = " to "), check.names = FALSE)
}))))
write_tsv(tab3, "Table3.tsv")
cat("17: Tables 1-3 written\n"); print(tab1, row.names = FALSE)
