# 09. Is there any association of repeat length once the other signals at PER3 are accounted for?
#   (a) a 'lineage model': the joint effects of markers of the repeat's lineages and sub-lineages,
#       rs1773135-A (five units: 5A+5B), rs6680585-A (5A), rs61773390-T (4B.1), rs228697-G (4B.2) and
#       rs150812083-G (a rare branch of 4B.1); the reference haplotype is the main four-unit lineage, 4A. The
#       4B.1 and 4B.2 markers are among the strongest chronotype variants at PER3 (they were chosen because they
#       mark those signals); the question is whether repeat length adds anything to them. With the 5A marker in
#       the model, the rs1773135 term is 5B vs 4A and the rs6680585 term 5A vs 5B; 5A vs 4A is their sum, and
#       five-unit haplotypes as a whole vs 4A weight 5A and 5B by their frequencies (with standard errors from
#       the covariance of the joint effects); a 2-df Wald test asks whether 5A and 5B both equal 4A;
#   (b) stepwise conditional and joint analysis (COJO; P < 5e-8) of the 500-kb window, by this R implementation and
#       by GCTA 1.94.1 when it is available; the tag of repeat length conditional on the signals selected;
#   (c) the single-SNP approximate conditional test of the student guide (on rs61773390 only), for comparison;
#   (d) each with three LD reference panels: 1000 Genomes EUR (503), EUR without FIN (404), CEU+GBR (190);
#   (e) the R implementation checked against GCTA (joint model of the lineage markers).
# Traits: chronotype, morning person, and the accelerometer measures of sleep timing (midpoint, L5, M10).
# Outputs (outputs/tables): cojo_reference_panels.tsv, lineage_model.tsv, cojo_selected.tsv,
#   cojo_tag_conditional.tsv, single_snp_conditional.tsv, gcta_crosscheck.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))
source(file.path(ROOT, "scripts", "cojo.R"))

vcf <- file.path(DATA, "1000g", "1kgp2504_chr1_7.58-8.08Mb.vcf.gz")
smp <- read_tsv(file.path(DATA, "1000g", "1kgp_samples.tsv"))
panels <- list(
  EUR = smp$sample[smp$unrelated_2504 == "yes" & smp$superpopulation == "EUR"],
  EUR_noFIN = smp$sample[smp$unrelated_2504 == "yes" & smp$superpopulation == "EUR" & smp$population != "FIN"],
  CEU_GBR = smp$sample[smp$unrelated_2504 == "yes" & smp$population %in% c("CEU", "GBR")])
gcta <- Sys.which("gcta64"); if (gcta == "") gcta <- file.path(ROOT, "tools", "gcta64")
have_gcta <- file.exists(gcta)
gdir <- file.path(WORK, "gcta"); dir.create(gdir, showWarnings = FALSE)

# dosage matrix (ALT allele) of the biallelic SNVs and short indels of a panel
read_dosage <- function(samples) {
  s <- tempfile(); writeLines(samples, s); tmp <- tempfile()
  system2("bcftools", c("view", "-S", s, "-m2", "-M2", "-Ou", vcf, "|", "bcftools", "query",
                        "-f", shQuote("%POS\t%REF\t%ALT[\t%GT]\n")), stdout = tmp)
  d <- read.delim(tmp, header = FALSE, colClasses = "character")
  info <- data.frame(pos = as.integer(d[[1]]), ref = d[[2]], alt = d[[3]], stringsAsFactors = FALSE)
  g <- as.matrix(d[, -(1:3)])
  dos <- matrix(suppressWarnings(as.integer(substr(g, 1, 1)) + as.integer(substr(g, 3, 3))), nrow(g))
  keep <- nchar(info$ref) <= 50 & nchar(info$alt) <= 50 & rowSums(is.na(dos)) == 0
  list(info = info[keep, ], dos = dos[keep, , drop = FALSE])
}
db <- read.delim(text = bcf_query(file.path(DATA, "reference", "dbsnp157_per3_500kb.vcf.gz"), "%POS\t%REF\t%ALT\t%ID\n"),
                 header = FALSE, col.names = c("pos", "ref", "alt", "rsid"), stringsAsFactors = FALSE)
alts <- strsplit(db$alt, ",")
db <- data.frame(key = paste(rep(db$pos, lengths(alts)), rep(db$ref, lengths(alts)), unlist(alts)),
                 rsid = rep(db$rsid, lengths(alts)), stringsAsFactors = FALSE)
gw <- function(trait) read.delim(file.path(WORK, paste0("gwas_harmonised_", trait, ".tsv")), stringsAsFactors = FALSE)
palin <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")

# match a GWAS to a reference panel (by rsID; same alleles; frequencies within 0.2, palindromic SNPs only if
# MAF < 0.4 and frequencies within 0.1; at least two copies in the panel; INFO >= 0.8); LD matrix of the effect alleles
build <- function(trait, ref) {
  g <- gw(trait)
  rs <- db$rsid[match(paste(ref$info$pos, ref$info$ref, ref$info$alt), db$key)]
  m <- match(rs, g$rsid); k <- which(!is.na(m))
  g <- g[m[k], ]; info <- ref$info[k, ]; dos <- ref$dos[k, , drop = FALSE]
  same <- g$a1 == info$alt & g$a0 == info$ref; swap <- g$a1 == info$ref & g$a0 == info$alt
  dos[swap, ] <- 2 - dos[swap, ]
  f_ref <- rowMeans(dos) / 2
  pal <- palin(g$a1, g$a0)
  ok <- (same | swap) & abs(f_ref - g$freq1) < 0.2 &
        !(pal & (pmin(g$freq1, 1 - g$freq1) > 0.4 | abs(f_ref - g$freq1) > 0.1)) &
        pmin(f_ref, 1 - f_ref) * 2 * ncol(dos) >= 2 & g$freq1 > 0.001 & g$freq1 < 0.999 &
        (is.na(g$info) | g$info >= 0.8) & apply(dos, 1, sd) > 0
  g <- g[ok, ]; dos <- dos[ok, , drop = FALSE]; info <- info[ok, ]
  dup <- duplicated(g$rsid); g <- g[!dup, ]; dos <- dos[!dup, , drop = FALSE]; info <- info[!dup, ]
  s <- data.frame(snp = g$rsid, a1 = g$a1, a0 = g$a0, freq = g$freq1, b = g$beta, se = g$se, p = g$p, n = g$n,
                  pos38 = info$pos, freq_ref = rowMeans(dos) / 2, stringsAsFactors = FALSE)
  list(s = s, R = cor(t(dos)), dos = dos, matched = length(k))
}

write_bed <- function(prefix, b, samples) {
  dos <- b$dos; n <- ncol(dos); nb <- ceiling(n / 4)
  code <- matrix(0L, nrow(dos), nb * 4)
  code[, seq_len(n)] <- ifelse(dos == 2, 0L, ifelse(dos == 1, 2L, 3L))     # A1 = GWAS effect allele
  con <- file(paste0(prefix, ".bed"), "wb"); writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (i in seq_len(nrow(dos))) { q <- matrix(code[i, ], nrow = 4); writeBin(as.raw(q[1, ] + 4L * q[2, ] + 16L * q[3, ] + 64L * q[4, ]), con) }
  close(con)
  write.table(data.frame(1, b$s$snp, 0, b$s$pos38, b$s$a1, b$s$a0), paste0(prefix, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(data.frame(samples, samples, 0, 0, 0, -9), paste0(prefix, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  write.table(data.frame(SNP = b$s$snp, A1 = b$s$a1, A2 = b$s$a0, freq = b$s$freq, b = b$s$b, se = b$s$se, p = b$s$p, N = b$s$n),
              paste0(prefix, ".ma"), quote = FALSE, row.names = FALSE, sep = "\t")
}
run_gcta <- function(args) system2(gcta, args, stdout = FALSE, stderr = FALSE)

markers <- data.frame(snp = c("rs1773135", "rs6680585", "rs61773390", "rs228697", "rs150812083"),
                      allele = c("A", "A", "T", "G", "G"),
                      contrast = c("5B vs 4A", "5A vs 5B", "4B.1 vs 4A", "4B.2 vs 4A", "rare branch of 4B.1 (vs rest of 4B.1)"),
                      stringsAsFactors = FALSE)
traits_main <- c("chronotype", "morning", "acc_midpoint", "acc_l5", "acc_m10")

refs <- lapply(panels, read_dosage)
out <- list(panel = list(), lin = list(), sel = list(), tag = list(), ss = list(), chk = list())
for (pn in names(panels)) for (tr in traits_main) {
  key <- paste(pn, tr)
  b <- build(tr, refs[[pn]]); cj <- cojo_prepare(b$s, b$R)
  out$panel[[key]] <- data.frame(panel = pn, people = length(panels[[pn]]), trait = tr, variants_in_panel = nrow(refs[[pn]]$info),
                                 matched_to_gwas = b$matched, used = nrow(b$s), Vp = signif(cj$Vp, 4), N = round(cj$N))
  # (a) lineage model
  idx <- match(markers$snp, b$s$snp)
  if (anyNA(idx)) stop("lineage marker missing: ", tr, " ", pn, " ", paste(markers$snp[is.na(idx)], collapse = ","))
  jt <- cojo_joint(cj, idx); sg <- ifelse(b$s$a1[idx] == markers$allele, 1, -1)
  bj <- sg * jt$bJ; V <- attr(jt, "vcov") * outer(sg, sg)          # effects and covariance, marker alleles
  f <- ifelse(sg == 1, b$s$freq[idx], 1 - b$s$freq[idx])
  lin_m <- data.frame(panel = pn, trait = tr, marker = markers$snp, allele = markers$allele, contrast = markers$contrast,
                      type = "marker", w_5A = NA, freq = f, b = sg * jt$b, se = jt$se, p = jt$p, bJ = bj, bJ_se = jt$bJ_se, pJ = jt$pJ)
  # the five-unit lineages against 4A (markers 1 and 2: rs1773135 = five units, rs6680585 = 5A)
  w <- f[2] / f[1]                                                   # share of 5A among five-unit haplotypes
  ctr <- list(c(1, 1, 0, 0, 0), c(1, w, 0, 0, 0))
  est <- sapply(ctr, function(cc) sum(cc * bj)); ses <- sapply(ctr, function(cc) sqrt(drop(cc %*% V %*% cc)))
  chi2 <- drop(bj[1:2] %*% solve(V[1:2, 1:2], bj[1:2]))
  lin_d <- data.frame(panel = pn, trait = tr, marker = c("rs1773135+rs6680585", "rs1773135+w*rs6680585", "rs1773135,rs6680585"),
                      allele = c("A+A", "A+A", "A,A"),
                      contrast = c("5A vs 4A", "five units (5A and 5B by frequency) vs 4A", "5A and 5B vs 4A (2-df test)"),
                      type = c("derived", "derived", "test"), w_5A = c(NA, round(w, 4), NA), freq = c(f[2], f[1], NA),
                      b = NA, se = NA, p = NA, bJ = c(est, NA), bJ_se = c(ses, NA),
                      pJ = c(2 * pnorm(-abs(est / ses)), pchisq(chi2, 2, lower.tail = FALSE)))
  out$lin[[key]] <- rbind(lin_m, lin_d)
  # (b) stepwise selection: R implementation, and GCTA
  t_i <- match("rs1773135", b$s$snp); st <- if (b$s$a1[t_i] == "A") 1 else -1
  sels <- list(R = cojo_select(cj, p_thresh = 5e-8))
  pre <- file.path(gdir, paste0(pn, "_", tr))
  if (have_gcta) {
    write_bed(pre, b, panels[[pn]])
    unlink(paste0(pre, "_slct.jma.cojo"))
    run_gcta(c("--bfile", pre, "--cojo-file", paste0(pre, ".ma"), "--cojo-slct", "--cojo-p", "5e-8", "--out", paste0(pre, "_slct")))
    f <- paste0(pre, "_slct.jma.cojo")
    sels$GCTA <- if (file.exists(f)) match(read.delim(f, stringsAsFactors = FALSE)$SNP, b$s$snp) else integer(0)
  }
  for (meth in names(sels)) {
    sel <- sels[[meth]]
    if (length(sel)) {
      j <- cojo_joint(cj, sel)
      out$sel[[paste(key, meth)]] <- data.frame(panel = pn, trait = tr, method = meth, j, pos38 = b$s$pos38[sel], a1 = b$s$a1[sel],
                                                freq = b$s$freq[sel])
      cc <- cojo_cond(cj, sel)[t_i, ]
    } else cc <- data.frame(bC = b$s$b[t_i], bC_se = b$s$se[t_i], pC = b$s$p[t_i])
    out$tag[[paste(key, meth)]] <- data.frame(panel = pn, trait = tr, method = meth, signals = length(sel),
                                              conditioned_on = paste(b$s$snp[sel], collapse = ","),
                                              b = st * b$s$b[t_i], se = b$s$se[t_i], p = b$s$p[t_i],
                                              bC = st * cc$bC, bC_se = cc$bC_se, pC = cc$pC)
  }
  # (c) single-SNP conditional on rs61773390
  j2 <- match("rs61773390", b$s$snp); r <- b$R[t_i, j2]
  z1 <- b$s$b[t_i] / b$s$se[t_i]; z2 <- b$s$b[j2] / b$s$se[j2]; zc <- (z1 - r * z2) / sqrt(1 - r^2)
  out$ss[[key]] <- data.frame(panel = pn, trait = tr, r_five_with_4B1 = r * st * ifelse(b$s$a1[j2] == "T", 1, -1),
                              z_tag = st * z1, z_tag_given_rs61773390 = st * zc, p = 2 * pnorm(-abs(zc)))
  # (e) GCTA check of the lineage model
  if (have_gcta) {
    writeLines(markers$snp, paste0(pre, "_markers.txt"))
    run_gcta(c("--bfile", pre, "--cojo-file", paste0(pre, ".ma"), "--extract", paste0(pre, "_markers.txt"), "--cojo-joint",
               "--out", paste0(pre, "_joint")))
    gj <- read.delim(paste0(pre, "_joint.jma.cojo"), stringsAsFactors = FALSE)
    rj <- cojo_joint(cj, match(gj$SNP, b$s$snp))
    out$chk[[key]] <- data.frame(panel = pn, trait = tr, snp = gj$SNP, gcta_bJ = gj$bJ, gcta_se = gj$bJ_se, gcta_pJ = gj$pJ,
                                 r_bJ = rj$bJ, r_se = rj$bJ_se, r_pJ = rj$pJ,
                                 max_rel_diff_b = max(abs(gj$bJ - rj$bJ) / abs(gj$bJ)))
  }
}
tabs <- lapply(out, function(x) if (length(x)) do.call(rbind, x) else NULL)
write_tsv(tabs$panel, "cojo_reference_panels.tsv"); write_tsv(tabs$lin, "lineage_model.tsv")
write_tsv(tabs$sel, "cojo_selected.tsv"); write_tsv(tabs$tag, "cojo_tag_conditional.tsv")
write_tsv(tabs$ss, "single_snp_conditional.tsv")
if (!is.null(tabs$chk)) write_tsv(tabs$chk, "gcta_crosscheck.tsv")

cat("09: lineage model (EUR panel):\n")
x <- tabs$lin[tabs$lin$panel == "EUR", ]
print(data.frame(trait = x$trait, contrast = x$contrast, bJ = signif(x$bJ, 3), se = signif(x$bJ_se, 3), pJ = signif(x$pJ, 2)), row.names = FALSE)
cat("09: five units (5A and 5B by frequency) vs 4A in the lineage model, by panel:\n")
x <- tabs$lin[tabs$lin$contrast == "five units (5A and 5B by frequency) vs 4A", ]
print(data.frame(panel = x$panel, trait = x$trait, bJ = signif(x$bJ, 3), se = signif(x$bJ_se, 3), pJ = signif(x$pJ, 2)), row.names = FALSE)
cat("09: signals selected (P < 5e-8):\n")
print(tabs$sel[, c("panel", "trait", "method", "snp", "a1", "freq", "pJ")], row.names = FALSE)
cat("09: the tag, conditional on the selected signals:\n")
print(tabs$tag[, c("panel", "trait", "method", "signals", "b", "p", "bC", "bC_se", "pC")], row.names = FALSE)
cat("09: single-SNP conditional on rs61773390:\n"); print(tabs$ss, row.names = FALSE)
if (!is.null(tabs$chk)) cat(sprintf("09: GCTA vs R, joint lineage model: largest relative difference in bJ = %.2g over %d estimates\n",
                                    max(tabs$chk$max_rel_diff_b), nrow(tabs$chk)))
