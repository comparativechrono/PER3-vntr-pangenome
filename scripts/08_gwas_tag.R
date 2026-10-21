# 08. The repeat, through its tag, in GWAS of sleep and circadian traits of UK Biobank (European ancestry):
#     self-reported chronotype, morning person, sleep duration, insomnia and napping (~400,000-450,000 people), and
#     eight measures from a wrist accelerometer (~85,000). Effects are given for the five-unit allele
#     (rs1773135-A; rs228701-C, a strand-ambiguous C/G SNP, is shown as a check).
#     Also the marginal associations of the SNPs that mark the repeat's lineages and sub-lineages (step 07).
# Outputs (outputs/tables): gwas_traits.tsv, gwas_tag_results.tsv, gwas_lineage_marker_results.tsv,
#                           gwas_harmonised_<trait>.tsv (outputs/work), regional data for figures
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

G <- file.path(DATA, "gwas")
traits <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, text = "
trait|label|file|prefix|beta|se|p|freq|N|scale|study
chronotype|Chronotype (self-report; higher = morning)|chronotype_Jones2019.tsv||BETA|SE|P_BOLT_LMM|A1FREQ|449734|score units|Jones 2019, GCST007576
morning|Morning person (self-report)|morning_person_Jones2019.tsv||LOGOR|LOGOR_SE|P_BOLT_LMM|A1FREQ|403195|log odds ratio|Jones 2019, GCST007565
sleepdur|Sleep duration (self-report)|sleep_duration_Dashti2019.tsv||BETA_SLEEPDURATION|SE_SLEEPDURATION|P_SLEEPDURATION|A1FREQ|446118|hours|Dashti 2019, GCST007561
insomnia|Insomnia symptoms (self-report)|insomnia_Lane2019.tsv||BETA_INSOMNIA|SE_INSOMNIA|P_INSOMNIA|A1FREQ|453379|linear-model beta (case = 1)|Lane 2019, GCST007387
nap|Daytime napping (self-report)|daytime_nap_Dashti2021.tsv||beta|standard_error|p_value|effect_allele_frequency|452633|napping category|Dashti 2021, GCST011494
acc_midpoint|Sleep midpoint (accelerometer)|accelerometer_Jones2019.tsv|ACC_SLEEP_MIDP_RAW_SIN|BETA|SE|P|A1FREQ|84810|SD|Jones 2019b, GCST007807
acc_l5|L5 timing (accelerometer)|accelerometer_Jones2019.tsv|ACC_L5_TIME_RAW_SIN|BETA|SE|P|A1FREQ|85205|SD|Jones 2019b, GCST007805
acc_m10|M10 timing (accelerometer)|accelerometer_Jones2019.tsv|ACC_M10_TIME_RAW_SIN|BETA|SE|P|A1FREQ|85670|SD|Jones 2019b, GCST007806
acc_duration|Sleep duration (accelerometer)|accelerometer_Jones2019.tsv|ACC_SLEEP_DUR_RAW_SIN|BETA|SE|P|A1FREQ|85449|SD|Jones 2019b, GCST007801
acc_efficiency|Sleep efficiency (accelerometer)|accelerometer_Jones2019.tsv|ACC_SLEEP_EFF_RAW_SIN|BETA|SE|P|A1FREQ|84810|SD|Jones 2019b, GCST007803
acc_episodes|Number of sleep episodes (accelerometer)|accelerometer_Jones2019.tsv|ACC_N_SLEEP_EPISODES_RAW_SIN|BETA|SE|P|A1FREQ|84810|SD|Jones 2019b, GCST007804
acc_duration_sd|Variability of sleep duration (accelerometer)|accelerometer_Jones2019.tsv|ACC_SLEEP_DUR_SD_RAW_SIN|BETA|SE|P|A1FREQ|84441|SD|Jones 2019b, GCST007802
acc_inactivity|Diurnal inactivity (accelerometer)|accelerometer_Jones2019.tsv|ACC_DIURNAL_INACT_RAW_SIN|BETA|SE|P|A1FREQ|84757|SD|Jones 2019b, GCST007808
")
write_tsv(traits, "gwas_traits.tsv")
n_tests <- nrow(traits); alpha_bonf <- 0.05 / n_tests

# read one trait into a common format: rsid, chr, bp37, a1 (effect), a0, freq1, beta, se, p, n, info
load_trait <- function(i) {
  t <- traits[i, ]
  d <- read.delim(file.path(G, t$file), stringsAsFactors = FALSE, na.strings = c("NA", "-", ""), check.names = FALSE)
  col <- function(x) if (t$prefix == "") x else paste0(t$prefix, "_", x)
  if (t$trait == "nap") {
    out <- data.frame(rsid = d$variant_id, chr = d$chromosome, bp37 = d$base_pair_location, a1 = d$effect_allele,
                      a0 = d$other_allele, freq1 = d$effect_allele_frequency, beta = d$beta, se = d$standard_error,
                      p = d$p_value, info = d$INFO)
  } else {
    out <- data.frame(rsid = d$SNP, chr = d$CHR, bp37 = d$BP, a1 = d$ALLELE1, a0 = d$ALLELE0,
                      freq1 = d[[col(t$freq)]], beta = d[[col(t$beta)]], se = d[[col(t$se)]], p = d[[col(t$p)]],
                      info = d$INFO)
  }
  out$n <- t$N
  out <- out[!is.na(out$rsid) & !is.na(out$beta) & !is.na(out$se) & out$se > 0, ]
  out$z <- out$beta / out$se
  out
}
gw <- lapply(seq_len(nrow(traits)), load_trait); names(gw) <- traits$trait
for (k in names(gw)) write.table(gw[[k]], file.path(WORK, paste0("gwas_harmonised_", k, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# phenotypic SD implied by the summary statistics (median over SNPs of sqrt(2pq (N SE^2 + beta^2)))
sd_y <- sapply(gw, function(d) { q <- d$freq1 * (1 - d$freq1); k <- d$freq1 > 0.01 & d$freq1 < 0.99
  median(sqrt(2 * q[k] * (d$n[k] * d$se[k]^2 + d$beta[k]^2))) })

lookup <- function(d, rsid, allele) {
  x <- d[d$rsid == rsid, ]
  if (nrow(x) != 1) return(NULL)
  s <- if (x$a1 == allele) 1 else if (x$a0 == allele) -1 else stop("allele mismatch ", rsid)
  data.frame(rsid = rsid, allele = allele, freq = round(if (s == 1) x$freq1 else 1 - x$freq1, 4),
             beta = s * x$beta, se = x$se, z = s * x$z, p = x$p, n = x$n, info = x$info)
}

# the tag, for every trait
res <- do.call(rbind, lapply(seq_len(nrow(traits)), function(i) {
  d <- gw[[i]]; tr <- traits[i, ]
  do.call(rbind, lapply(list(c("rs1773135", "A"), c("rs228701", "C")), function(s) {
    x <- lookup(d, s[1], s[2]); if (is.null(x)) return(NULL)
    p2 <- 2 * pnorm(-abs(x$z))
    data.frame(trait = tr$trait, label = tr$label, scale = tr$scale, x,
               p_from_z = signif(p2, 3),
               ci_lower = x$beta - 1.959964 * x$se, ci_upper = x$beta + 1.959964 * x$se,
               sd_trait = sd_y[[tr$trait]],
               beta_sd = x$beta / sd_y[[tr$trait]], ci_lower_sd = (x$beta - 1.959964 * x$se) / sd_y[[tr$trait]],
               ci_upper_sd = (x$beta + 1.959964 * x$se) / sd_y[[tr$trait]],
               mde80_p05 = (qnorm(0.975) + qnorm(0.8)) * x$se, mde80_p5e8 = (qnorm(1 - 2.5e-8) + qnorm(0.8)) * x$se,
               mde80_p05_sd = (qnorm(0.975) + qnorm(0.8)) * x$se / sd_y[[tr$trait]],
               bonferroni_p = min(1, x$p * n_tests), study = tr$study)
  }))
}))
write_tsv(res, "gwas_tag_results.tsv")

# lineage markers (marginal associations; alleles that mark the lineage or sub-lineage)
markers <- data.frame(rsid = c("rs1773135", "rs2859387", "rs228663", "rs61773390", "rs10462020", "rs228697", "rs228654",
                               "rs228726", "rs150812083", "rs6680585"),
                      allele = c("A", "A", "A", "T", "G", "G", "A", "A", "G", "A"),
                      marks = c("five units (all five-unit lineages)", "4A", "4B", "4B.1 (with V647G)", "4B.1 (p.Val647Gly)",
                                "4B.2 (p.Pro864Ala)", "4B.2", "4B.2", "rare branch of 4B.1", "5A"), stringsAsFactors = FALSE)
lm_res <- do.call(rbind, lapply(seq_len(nrow(traits)), function(i) {
  do.call(rbind, lapply(seq_len(nrow(markers)), function(j) {
    x <- lookup(gw[[i]], markers$rsid[j], markers$allele[j]); if (is.null(x)) return(NULL)
    data.frame(trait = traits$trait[i], marks = markers$marks[j], x, beta_sd = x$beta / sd_y[[traits$trait[i]]])
  }))
}))
write_tsv(lm_res, "gwas_lineage_marker_results.tsv")

# accelerometer measures in their own units (the untransformed _RAW results): hours -> minutes for timing and duration
acc <- read.delim(file.path(G, "accelerometer_Jones2019.tsv"), stringsAsFactors = FALSE, na.strings = c("NA", "-", ""))
raw <- do.call(rbind, lapply(c("ACC_SLEEP_MIDP", "ACC_L5_TIME", "ACC_M10_TIME", "ACC_SLEEP_DUR"), function(pfx) {
  x <- acc[acc$SNP == "rs1773135", ]; sgn <- if (x$ALLELE1 == "A") 1 else -1
  b <- sgn * x[[paste0(pfx, "_RAW_BETA")]]; se <- x[[paste0(pfx, "_RAW_SE")]]
  data.frame(measure = pfx, beta_hours = b, se_hours = se, beta_minutes = 60 * b, ci_minutes = paste(round(60 * (b - 1.96 * se), 2), "to", round(60 * (b + 1.96 * se), 2)),
             p = x[[paste0(pfx, "_RAW_P")]])
}))
write_tsv(raw, "gwas_tag_accelerometer_raw_units.tsv")
cat("08: the five-unit allele and accelerometer measures in minutes:\n"); print(raw, row.names = FALSE)
cat(sprintf("08: %d traits; Bonferroni threshold %.4f\n", n_tests, alpha_bonf))
r1 <- res[res$rsid == "rs1773135", ]
print(data.frame(trait = r1$trait, freq = r1$freq, beta = signif(r1$beta, 3), se = signif(r1$se, 3), p = signif(r1$p, 2),
                 beta_sd = round(r1$beta_sd, 4), ci_sd = paste0(round(r1$ci_lower_sd, 4), " to ", round(r1$ci_upper_sd, 4)),
                 mde80_sd = round(r1$mde80_p05_sd, 4)), row.names = FALSE)
chk <- merge(res[res$rsid == "rs1773135", c("trait", "z")], res[res$rsid == "rs228701", c("trait", "z")], by = "trait")
cat("08: largest difference in z between the two tags:", round(max(abs(chk$z.x - chk$z.y)), 3), "\n")
cat("08: chronotype, lineage markers:\n")
print(lm_res[lm_res$trait == "chronotype", c("marks", "rsid", "allele", "freq", "beta", "se", "p")], row.names = FALSE)
