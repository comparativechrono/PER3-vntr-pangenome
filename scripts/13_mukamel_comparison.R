# 13. The repeat in an earlier phenome-wide scan: Mukamel et al. (2021, Science 373:1499) imputed the length of
#     118 protein-coding VNTRs, this one among them (VNTR_HG38_chr1_7829874_7830126), in up to 337,466 unrelated
#     White British UK Biobank participants and tested 787 phenotypes. Their dosage is the summed length of the two
#     alleles rescaled to 0-2, so betas are not per allele; z-scores are compared. Effect allele: 'longer'.
#     Both analyses use UK Biobank: they are not independent samples, but independent ways of genotyping the repeat.
# Outputs (outputs/tables): mukamel_sleep_traits.tsv, mukamel_top_phenotypes.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

mk <- read_tsv(file.path(DATA, "mukamel2021", "mukamel2021_PER3_VNTR_787_phenotypes.tsv"))
mk$z <- mk$BETA / mk$SE
ours <- read_tsv(file.path(TAB, "gwas_tag_results.tsv")); ours <- ours[ours$rsid == "rs1773135", ]
pairs <- data.frame(
  mukamel = c("sleep_EVENING_PERSON", "sleep_HOURS", "sleep_INSOMNIA", "sleep_NAP_DURING_DAY", "sleep_GET_UP_EASE", "sleep_NARCOLEPSY", "sleep_SNORING"),
  ours = c("chronotype", "sleepdur", "insomnia", "nap", NA, NA, NA),
  # +1: same direction of coding as our trait; -1: opposite (EVENING_PERSON: higher = evening; our chronotype: higher = morning)
  sign = c(-1, 1, 1, 1, NA, NA, NA), stringsAsFactors = FALSE)
out <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  m <- mk[mk$phenotype == pairs$mukamel[i], ]
  o <- if (is.na(pairs$ours[i])) NULL else ours[ours$trait == pairs$ours[i], ]
  data.frame(mukamel_phenotype = pairs$mukamel[i], mukamel_z_longer = round(m$z, 2), mukamel_p = m$P_LINREG,
             our_trait = ifelse(is.null(o), "", pairs$ours[i]),
             our_z_five_aligned = if (is.null(o)) NA else round(pairs$sign[i] * o$z, 2),
             our_p = if (is.null(o)) NA else signif(o$p, 2),
             same_direction = if (is.null(o)) NA else sign(m$z) == sign(pairs$sign[i] * o$z))
}))
write_tsv(out, "mukamel_sleep_traits.tsv")
top <- mk[order(mk$P_LINREG), ][1:10, c("phenotype", "BETA", "SE", "P_LINREG", "z")]
write_tsv(top, "mukamel_top_phenotypes.tsv")
ms <- readLines(file.path(DATA, "mukamel2021", "mukamel2021_PER3_VNTR_dosage_mean_sd.txt"))
cat("13: Mukamel et al. 2021, PER3 VNTR dosage:", trimws(ms), "\n")
cat(sprintf("13: phenotypes tested: %d; smallest P: %s (%s); phenotypes with P < 0.05: %d (expected by chance: %.0f)\n",
            nrow(mk), format(min(mk$P_LINREG)), mk$phenotype[which.min(mk$P_LINREG)], sum(mk$P_LINREG < 0.05), 0.05 * nrow(mk)))
cat(sprintf("13: Bonferroni threshold for 787 x 118 tests: %.1e\n", 0.05 / (787 * 118)))
print(out, row.names = FALSE)
