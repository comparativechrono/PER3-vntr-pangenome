# 10. Fine-mapping of the PER3 window with SuSiE-RSS (Zou et al. 2022): which variants are in the 95% credible sets
#     of chronotype, morning person and sleep midpoint, which lineage does each set mark, and what is the posterior
#     inclusion probability of the repeat's tags? LD from 1000 Genomes EUR (503; the main analysis), EUR without FIN
#     (404) and CEU+GBR (190), as in step 09.
#     Needs the R package susieR; skipped (with a note) if it is not installed.
# Outputs (outputs/tables): susie_credible_sets.tsv, susie_tags.tsv, susie_diagnostics.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))
if (!requireNamespace("susieR", quietly = TRUE)) {
  cat("10: susieR is not installed; fine-mapping skipped\n"); quit(save = "no", status = 0)
}
suppressPackageStartupMessages(library(susieR))
set.seed(1)

vcf <- file.path(DATA, "1000g", "1kgp2504_chr1_7.58-8.08Mb.vcf.gz")
smp <- read_tsv(file.path(DATA, "1000g", "1kgp_samples.tsv"))
panels <- list(
  EUR = smp$sample[smp$unrelated_2504 == "yes" & smp$superpopulation == "EUR"],
  EUR_noFIN = smp$sample[smp$unrelated_2504 == "yes" & smp$superpopulation == "EUR" & smp$population != "FIN"],
  CEU_GBR = smp$sample[smp$unrelated_2504 == "yes" & smp$population %in% c("CEU", "GBR")])

db <- read.delim(text = bcf_query(file.path(DATA, "reference", "dbsnp157_per3_500kb.vcf.gz"), "%POS\t%REF\t%ALT\t%ID\n"),
                 header = FALSE, col.names = c("pos", "ref", "alt", "rsid"), stringsAsFactors = FALSE)
alts <- strsplit(db$alt, ",")
db <- data.frame(key = paste(rep(db$pos, lengths(alts)), rep(db$ref, lengths(alts)), unlist(alts)), rsid = rep(db$rsid, lengths(alts)))

cs_rows <- list(); tag_rows <- list(); diag_rows <- list()
for (pn in names(panels)) {
m <- haplotype_matrix(vcf, "chr1:7580000-8080000", panels[[pn]])
keep <- nchar(m$info$ref) <= 50 & nchar(m$info$alt) <= 50 & !grepl(",", m$info$alt) & rowSums(is.na(m$h)) == 0
info <- m$info[keep, ]; h <- m$h[keep, ]
dos <- h[, seq(1, ncol(h), 2)] + h[, seq(2, ncol(h), 2)]
gh <- function(pos, ref, alt) h[which(info$pos == pos & info$ref == ref & info$alt == alt), ]
five <- as.integer(gh(7830204, "A", "G") == 0)
lin <- ifelse(five == 1, ifelse(gh(7830004, "A", "G") == 1, "5A", "5B"),
              ifelse(gh(7832049, "T", "A") == 1, ifelse(gh(7824465, "G", "T") == 1, "4B.1", ifelse(gh(7827519, "C", "G") == 1, "4B.2", "4B.0")), "4A"))
rs <- db$rsid[match(paste(info$pos, info$ref, info$alt), db$key)]
for (tr in c("chronotype", "morning", "acc_midpoint")) {
  set.seed(1)
  key <- paste(pn, tr)
  g <- read.delim(file.path(WORK, paste0("gwas_harmonised_", tr, ".tsv")), stringsAsFactors = FALSE)
  mm <- match(rs, g$rsid); k <- which(!is.na(mm)); gg <- g[mm[k], ]
  d <- dos[k, , drop = FALSE]; hh <- h[k, , drop = FALSE]; inf <- info[k, ]
  same <- gg$a1 == inf$alt & gg$a0 == inf$ref; swap <- gg$a1 == inf$ref & gg$a0 == inf$alt
  d[swap, ] <- 2 - d[swap, ]; hh[swap, ] <- 1 - hh[swap, ]
  fr <- rowMeans(d) / 2
  pal <- paste0(gg$a1, gg$a0) %in% c("AT", "TA", "CG", "GC")
  ok <- (same | swap) & abs(fr - gg$freq1) < 0.2 & !(pal & (pmin(gg$freq1, 1 - gg$freq1) > 0.4 | abs(fr - gg$freq1) > 0.1)) &
        pmin(fr, 1 - fr) * 2 * ncol(d) >= 2 & (is.na(gg$info) | gg$info >= 0.8) & apply(d, 1, sd) > 0 & !duplicated(gg$rsid)
  gg <- gg[ok, ]; d <- d[ok, , drop = FALSE]; hh <- hh[ok, , drop = FALSE]; inf <- inf[ok, ]
  R <- cor(t(d)); z <- gg$beta / gg$se; n <- gg$n[1]
  lam <- estimate_s_rss(z, R, n = n)
  fit <- susie_rss(z = z, R = R, n = n, L = 10, estimate_residual_variance = FALSE, max_iter = 500, coverage = 0.95,
                   min_abs_corr = 0.5)
  pip <- fit$pip
  diag_rows[[key]] <- data.frame(panel = pn, trait = tr, variants = length(z), n = n, ld_mismatch_s = round(lam, 4),
                                credible_sets = length(fit$sets$cs), converged = fit$converged)
  sets <- fit$sets$cs
  for (ci in seq_along(sets)) {
    idx <- sets[[ci]]
    lead <- idx[which.max(pip[idx])]
    r2l <- sapply(c("4A", "4B.1", "4B.2", "5A", "5B"), function(l) r2(hh[lead, ], as.integer(lin == l)))
    r2l <- c(r2l, `4B` = r2(hh[lead, ], as.integer(substr(lin, 1, 2) == "4B")), five = r2(hh[lead, ], five))
    best <- names(which.max(r2l))
    # on which lineages is the minor allele of the lead variant? (haplotype counts)
    minor <- if (mean(hh[lead, ]) < 0.5) 1 else 0
    carr <- table(factor(lin[hh[lead, ] == minor], c("4A", "4B.0", "4B.1", "4B.2", "5A", "5B")))
    minor_allele <- if (minor == 1) gg$a1[lead] else gg$a0[lead]
    cs_rows[[paste(key, ci)]] <- data.frame(panel = pn, trait = tr, credible_set = ci, size = length(idx),
                                           lead_minor_allele = minor_allele,
                                           minor_allele_haplotypes_by_lineage = paste(paste0(names(carr), ":", as.integer(carr)), collapse = " "),
                                           purity_min_abs_r = round(fit$sets$purity[ci, "min.abs.corr"], 3),
                                           lead = gg$rsid[lead], lead_pip = round(pip[lead], 3), lead_a1 = gg$a1[lead],
                                           lead_freq = round(gg$freq1[lead], 4), lead_z = round(z[lead], 2),
                                           members = paste(gg$rsid[idx][order(-pip[idx])], collapse = ","),
                                           best_lineage = best, r2_best_lineage = round(r2l[[best]], 3),
                                           r2_with_five = round(r2l[["five"]], 3),
                                           in_set_rs1773135 = "rs1773135" %in% gg$rsid[idx], in_set_rs228701 = "rs228701" %in% gg$rsid[idx])
  }
  for (t in c("rs1773135", "rs228701")) {
    i <- match(t, gg$rsid)
    tag_rows[[paste(key, t)]] <- data.frame(panel = pn, trait = tr, tag = t, z = if (is.na(i)) NA else round(z[i], 2),
                                           pip = if (is.na(i)) NA else signif(pip[i], 3), in_any_credible_set = !is.na(i) && any(sapply(sets, function(s) i %in% s)))
  }
}
}
write_tsv(do.call(rbind, cs_rows), "susie_credible_sets.tsv")
write_tsv(do.call(rbind, tag_rows), "susie_tags.tsv")
write_tsv(do.call(rbind, diag_rows), "susie_diagnostics.tsv")
cat("10: SuSiE-RSS diagnostics:\n"); print(do.call(rbind, diag_rows), row.names = FALSE)
cat("10: credible sets:\n"); x <- do.call(rbind, cs_rows); print(x[, setdiff(names(x), "members")], row.names = FALSE)
cat("10: the repeat's tags:\n"); print(do.call(rbind, tag_rows), row.names = FALSE)
cat("10: susieR version", as.character(packageVersion("susieR")), "\n")
