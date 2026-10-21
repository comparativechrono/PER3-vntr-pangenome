# 15. The figures of the paper (PDF and 300-dpi PNG in outputs/figures), from the tables of steps 02-14.
#   Figure 1  The repeat and its four lineages in 462 long-read haplotypes; frequencies in 26 populations
#   Figure 2  Linkage disequilibrium between repeat length and SNVs; the tag pair
#   Figure 3  Effects of repeat length and of the lineages on sleep and circadian traits in UK Biobank
#   Figure 4  The 4B lineage, its two sub-lineages and the PER3 chronotype signal
#   Figure 5  The repeat in other primates
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

FONT <- tryCatch(if (any(grepl("Liberation Sans", system("fc-list", intern = TRUE, ignore.stderr = TRUE)))) "Liberation Sans" else "sans",
                 error = function(e) "sans", warning = function(w) "sans")
SUB_COL <- c("4B.1" = "#1baf7a", "4B.2" = "#6250d6", "five" = "#eb6834", other = "#c3c2b7")   # validated all-pairs
save_fig <- function(name, width_mm, height_mm, draw) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  if (capabilities("cairo")) {
    cairo_pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h, family = FONT); draw(); invisible(dev.off())
    png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 300, type = "cairo", family = FONT); draw(); invisible(dev.off())
  } else {
    pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h); draw(); invisible(dev.off())
    png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 300); draw(); invisible(dev.off())
  }
}
lab <- function(l) mtext(l, side = 3, line = 0.5, adj = 0, font = 2, cex = 1, at = grconvertX(0.012, "nfc", "user"), xpd = NA)
axis_style <- function(side, at = NULL, labels = TRUE, ...) axis(side, at = at, labels = labels, col = AXIS, col.ticks = AXIS, col.axis = INK2, lwd = 0.8, tcl = -0.25, mgp = c(2, 0.4, 0), cex.axis = 0.75, ...)
grid_h <- function(at) abline(h = at, col = GRID, lwd = 0.6)
grid_v <- function(at) abline(v = at, col = GRID, lwd = 0.6)

hp <- read_tsv(file.path(TAB, "hprc_haplotype_alleles.tsv"))
cat_units <- read_tsv(file.path(TAB, "unit_catalogue.tsv"))

# ---------------------------------------------------------------------------------------------- Figure 1
fig1 <- function() {
  layout(matrix(c(1, 1, 2, 3, 4, 4), 3, byrow = TRUE), heights = c(0.95, 1.35, 1.5), widths = c(1.15, 1))
  # A: gene model and the repeat
  par(mar = c(1.2, 1, 2, 1), xpd = FALSE)
  ex <- read.delim(file.path(DATA, "reference", "PER3_exons.bed"), header = FALSE)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  lab("A")
  g0 <- min(ex$V2); g1 <- max(ex$V3); sx <- function(p) 0.04 + 0.92 * (p - g0) / (g1 - g0)
  segments(sx(g0), 0.8, sx(g1), 0.8, col = INK2, lwd = 1)
  rect(sx(ex$V2), 0.74, sx(ex$V3), 0.86, col = ifelse(grepl("19", ex$V4), "#0b0b0b", "#898781"), border = NA)
  text(0.04, 0.95, "PER3 (ENST00000377532), chr1:7.78-7.85 Mb (GRCh38), 22 exons", adj = 0, cex = 0.75, col = INK2, font = 3)
  e19 <- ex[grepl("19", ex$V4), ]
  text(sx(e19$V2), 0.65, "exon 19", cex = 0.7, col = INK)
  # zoom: exon 19 and its five units
  z0 <- 0.18; z1 <- 0.82; zx <- function(p) z0 + (z1 - z0) * (p - 7829834) / (7830161 - 7829834)
  segments(sx(e19$V2), 0.72, z0, 0.42, col = AXIS, lwd = 0.6); segments(sx(e19$V3), 0.72, z1, 0.42, col = AXIS, lwd = 0.6)
  rect(z0, 0.2, z1, 0.4, col = "#f0efec", border = AXIS, lwd = 0.6)
  st <- 7829888 + 54 * (0:4)
  rect(zx(st), 0.22, zx(st + 53), 0.38, col = "white", border = INK, lwd = 0.7)
  text(zx(st + 27), 0.30, paste0("R", 1:5), cex = 0.72)
  text(z0 - 0.01, 0.30, "exon 19", adj = 1, cex = 0.7, col = INK2)
  text((z0 + z1) / 2, 0.06, "five 54-bp units in GRCh38 (chr1:7,829,888-7,830,157) = 5 x 18 amino acids, residues 981-1070 of the 1,210-aa protein",
       cex = 0.68, col = INK2)
  # B: the lineages, unit by unit
  par(mar = c(1.5, 1, 2.2, 0.5))
  plot(NA, xlim = c(0, 10), ylim = c(0.3, 4.9), axes = FALSE, xlab = "", ylab = "")
  lab("B")
  L <- list("5B" = c("R1", "R2", "R3", "R4", "R5"), "5A" = c("R1", "R2", "R3*", "R4", "R5"),
            "4B" = c("R1", "R3", "R4*", "R5"), "4A" = c("R1", "R3", "R4", "R5"))
  n <- table(hp$lineage); st_ <- tapply(hp$ser_thr, hp$lineage, function(x) x[1])
  for (i in seq_along(L)) {
    y <- i; l <- names(L)[i]; u <- L[[l]]
    rect(0, y - 0.32, 0.22, y + 0.32, col = LINEAGE_COL[l], border = NA)
    text(0.35, y, l, adj = 0, font = 2, cex = 0.85)
    for (j in seq_along(u)) {
      x0 <- 1.1 + (j - 1) * 1.15
      rect(x0, y - 0.3, x0 + 1.05, y + 0.3, col = "white", border = INK, lwd = 0.7)
      text(x0 + 0.525, y, sub("\\*", "", u[j]), cex = 0.72)
      if (grepl("\\*", u[j])) points(x0 + 0.95, y + 0.3, pch = 25, bg = INK, col = INK, cex = 0.55)
    }
    text(7.05, y, sprintf("%d haplotypes", n[[l]]), adj = 0, cex = 0.72, col = INK2)
    text(7.05, y - 0.27, sprintf("%d Ser/Thr", st_[[l]]), adj = 0, cex = 0.62, col = MUTED)
  }
  text(0, 4.75, "Lineages of the HPRC haplotypes (5 of 462 are other variants)", adj = 0, cex = 0.72, col = INK2)
  text(1.1, 0.42, "triangle: R4 c.8T>C, p.Met1037Thr (rs2640909) in 4B; R3 c.9A>G, synonymous (rs12023156) in 5A",
       adj = 0, cex = 0.6, col = MUTED)
  # C: lineage composition by group, HPRC (long reads) and 1000G (tags)
  lg_h <- read_tsv(file.path(TAB, "hprc_lineages_by_group.tsv"))
  lg_k <- read_tsv(file.path(TAB, "kg_lineages_by_group.tsv"))
  par(mar = c(2.6, 4.6, 2.2, 0.8))
  grp <- c(GROUPS, "All")
  ys <- rev(seq_along(grp)) * 2.6
  plot(NA, xlim = c(0, 1), ylim = c(0.8, max(ys) + 1.6), axes = FALSE, xlab = "", ylab = "")
  lab("C")
  grid_v(seq(0, 1, 0.25))
  for (i in seq_along(grp)) {
    g <- grp[i]
    hv <- unlist(lg_h[lg_h$group == g, c("4A", "4B", "5A", "5B", "other")]); hv <- hv / sum(hv)
    kv <- unlist(lg_k[lg_k$group == g, c("4A", "4B", "5A", "5B")])
    for (k in 1:2) {
      v <- if (k == 1) hv else kv; y <- ys[i] + (k == 1) * 0.95
      x <- c(0, cumsum(v)); cols <- LINEAGE_COL[names(v)]
      rect(x[-length(x)], y - 0.38, x[-1] - 0.004, y + 0.38, col = cols, border = NA)
      text(-0.02, y, c("long reads", "tags")[k], adj = 1, cex = 0.55, col = MUTED, xpd = NA)
    }
    text(-0.2, ys[i] + 0.47, g, adj = 1, cex = 0.75, font = 2, xpd = NA)
  }
  axis_style(1, at = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"))
  mtext("Share of haplotypes", side = 1, line = 1.4, cex = 0.65, col = INK2)
  legend("top", inset = c(0, -0.07), horiz = TRUE, legend = c("4A", "4B", "5A", "5B", "other"), fill = LINEAGE_COL,
         border = NA, bty = "n", cex = 0.62, x.intersp = 0.4, text.width = 0.12, xpd = NA)
  # D: five-unit frequency in 26 populations of 1000 Genomes
  fp <- read_tsv(file.path(TAB, "kg_freq_by_population.tsv"))
  fp <- fp[order(match(fp$superpopulation, GROUPS), fp$freq_five), ]
  par(mar = c(3.6, 3.2, 2.2, 0.8))
  x <- seq_len(nrow(fp)) + 1.2 * (match(fp$superpopulation, GROUPS) - 1)
  plot(NA, xlim = range(x) + c(-0.8, 0.8), ylim = c(0, 0.5), axes = FALSE, xlab = "", ylab = "")
  lab("D")
  grid_h(seq(0, 0.5, 0.1))
  segments(x, fp$ci_lower, x, fp$ci_upper, col = INK2, lwd = 0.9)
  points(x, fp$freq_five, pch = 21, bg = INK, col = "white", cex = 0.9, lwd = 0.8)
  axis_style(2, at = seq(0, 0.5, 0.1), las = 1)
  mtext("Five-unit allele frequency", side = 2, line = 2.1, cex = 0.65, col = INK2)
  text(x, -0.035, fp$population, srt = 90, adj = 1, cex = 0.55, xpd = NA, col = INK2)
  for (g in GROUPS) { xx <- x[fp$superpopulation == g]; text(mean(xx), 0.49, g, cex = 0.7, font = 2) }
  hf <- read_tsv(file.path(TAB, "hprc_freq_by_group.tsv"))
  text(max(x) + 0.6, 0.02, "1000 Genomes, 2,504 unrelated people; repeat length from rs1773135; bars: 95% CI", adj = 1, cex = 0.55, col = MUTED)
}
save_fig("Figure1", 180, 170, fig1)

# ---------------------------------------------------------------------------------------------- Figure 2
fig2 <- function() {
  ld <- read_tsv(file.path(TAB, "hprc_ld_vntr_snvs.tsv"))
  layout(matrix(1:2, 1), widths = c(2.1, 1))
  par(mar = c(3.4, 3.4, 2, 0.8))
  xr <- c(7.58, 8.08)
  plot(NA, xlim = xr, ylim = c(-0.12, 1.05), axes = FALSE, xlab = "", ylab = "")
  lab("A")
  grid_h(seq(0, 1, 0.2))
  o <- order(ld$r2)
  points(ld$pos[o] / 1e6, ld$r2[o], pch = 16, cex = 0.45, col = adjustcolor(INK2, 0.45))
  tg <- ld[ld$rsid %in% c("rs228701", "rs1773135"), ]
  points(tg$pos / 1e6, tg$r2, pch = 21, bg = "#eb6834", col = "white", cex = 1.1, lwd = 0.8)
  text(tg$pos[1] / 1e6 - 0.006, tg$r2[1] + 0.035, "rs228701 and rs1773135 (r2 = 0.99)", adj = 1, cex = 0.62, col = INK)
  # PER3 gene track
  ex <- read.delim(file.path(DATA, "reference", "PER3_exons.bed"), header = FALSE)
  rect(min(ex$V2) / 1e6, -0.1, max(ex$V3) / 1e6, -0.06, col = INK2, border = NA)
  text(max(ex$V3) / 1e6 + 0.005, -0.08, "PER3", adj = 0, cex = 0.6, font = 3, col = INK2)
  abline(v = 7.83, col = "#eb6834", lwd = 0.5, lty = 1)
  axis_style(1, at = seq(7.6, 8.05, 0.1)); axis_style(2, at = seq(0, 1, 0.2), las = 1)
  mtext("Position on chromosome 1 (Mb, GRCh38)", side = 1, line = 1.7, cex = 0.65, col = INK2)
  mtext(expression(italic(r)^2 * " with repeat length"), side = 2, line = 2.1, cex = 0.65, col = INK2)
  # B: r2 of the tag with the repeat, by group (HPRC), and between the two tags (1000G)
  t5 <- read_tsv(file.path(TAB, "hprc_tags_top.tsv")); t5 <- t5[t5$rsid == "rs1773135", ]
  kt <- read_tsv(file.path(TAB, "kg_tag_ld.tsv"))
  par(mar = c(3.4, 4.2, 2, 0.6))
  g <- c(GROUPS, "All"); v <- c(unlist(t5[paste0("r2_", GROUPS)]), t5$r2)
  y <- rev(seq_along(g))
  plot(NA, xlim = c(0.9, 1.0), ylim = c(0.5, length(g) + 0.5), axes = FALSE, xlab = "", ylab = "")
  lab("B")
  grid_v(seq(0.9, 1, 0.025))
  segments(0.9, y, v, y, col = AXIS, lwd = 0.8)
  points(v, y, pch = 21, bg = INK, col = "white", cex = 1.1)
  text(v - 0.003, y + 0.3, sprintf("%.3f", v), cex = 0.58, col = INK2, adj = 1)
  axis_style(1, at = seq(0.9, 1, 0.025), labels = c("0.9", "", "0.95", "", "1"))
  axis(2, at = y, labels = g, las = 1, tick = FALSE, cex.axis = 0.72, col.axis = INK, mgp = c(2, 0.3, 0))
  mtext(expression(italic(r)^2 * ", tag and repeat (long reads)"), side = 1, line = 1.7, cex = 0.62, col = INK2)
}
save_fig("Figure2", 180, 80, fig2)

# ---------------------------------------------------------------------------------------------- Figure 3
fig3 <- function() {
  tr <- read_tsv(file.path(TAB, "gwas_tag_results.tsv")); tr <- tr[tr$rsid == "rs1773135", ]
  lm_ <- read_tsv(file.path(TAB, "lineage_model.tsv")); lm_ <- lm_[lm_$panel == "EUR", ]
  sdt <- setNames(tr$sd_trait, tr$trait)
  fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(signif(p, 2), format = "fg", digits = 2, flag = "#"))
  layout(matrix(1:4, 1), widths = c(1.75, 1.05, 0.9, 0.9))
  # A: the five-unit allele, 13 traits, SD units
  par(mar = c(3.4, 12.4, 2.2, 2.6))
  y <- rev(seq_len(nrow(tr)))
  xl <- c(-0.032, 0.032)
  plot(NA, xlim = xl, ylim = c(0.5, nrow(tr) + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  lab("A")
  grid_v(seq(-0.03, 0.03, 0.01)); abline(v = 0, col = INK2, lwd = 0.7)
  segments(tr$ci_lower_sd, y, tr$ci_upper_sd, y, col = INK2, lwd = 1)
  points(tr$beta_sd, y, pch = 21, bg = ifelse(tr$p < 0.05 / 13, "#eb6834", INK), col = "white", cex = 1.1)
  axis(2, at = y, labels = tr$label, las = 1, tick = FALSE, cex.axis = 0.62, col.axis = INK, mgp = c(2, 0.2, 0))
  text(xl[2] + 0.001, y, fp(tr$p), cex = 0.52, col = MUTED, adj = 0, xpd = NA)
  text(xl[2] + 0.001, nrow(tr) + 0.75, "P", cex = 0.55, col = MUTED, adj = 0, xpd = NA, font = 3)
  axis_style(1, at = seq(-0.03, 0.03, 0.01), labels = c("-0.03", "", "", "0", "", "", "0.03"), gap.axis = -1)
  mtext("Five-unit allele (SD per allele, 95% CI)", side = 1, line = 1.7, cex = 0.58, col = INK2)
  # B-D: lineage model, three traits, every contrast against 4A (the rare branch of 4B.1 is in Table 3)
  con <- c("five units (5A and 5B by frequency) vs 4A", "5A vs 4A", "5B vs 4A", "4B.1 vs 4A", "4B.2 vs 4A")
  short <- c("five units vs 4A", "5A vs 4A", "5B vs 4A", "4B.1 vs 4A", "4B.2 vs 4A")
  pc <- c("#eb6834", LINEAGE_COL[["5A"]], LINEAGE_COL[["5B"]], "#1baf7a", "#6250d6")
  for (k in 1:3) {
    t <- c("chronotype", "morning", "acc_midpoint")[k]
    x <- lm_[lm_$trait == t, ]; x <- x[match(con, x$contrast), ]
    s <- if (t == "morning") 1 else sdt[[t]]
    b <- x$bJ / s; lo <- (x$bJ - 1.96 * x$bJ_se) / s; hi <- (x$bJ + 1.96 * x$bJ_se) / s
    par(mar = c(3.4, if (k == 1) 6.4 else 0.8, 2.2, 0.8))
    yy <- rev(seq_along(con))
    rng <- c(-1, 1) * max(abs(c(lo, hi))) * 1.08
    plot(NA, xlim = rng, ylim = c(0.5, length(con) + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    lab(LETTERS[k + 1])
    grid_v(pretty(rng, 4)); abline(v = 0, col = INK2, lwd = 0.7)
    segments(lo, yy, hi, yy, col = INK2, lwd = 1)
    points(b, yy, pch = 21, bg = pc, col = "white", cex = 1.15)
    if (k == 1) axis(2, at = yy, labels = short, las = 1, tick = FALSE, cex.axis = 0.62, col.axis = INK, mgp = c(2, 0.2, 0))
    axis_style(1, at = pretty(rng, 4), gap.axis = -1)
    mtext(c("Chronotype (SD)", "Morning person (log OR)", "Sleep midpoint, accel. (SD)")[k], side = 1, line = 1.7, cex = 0.58, col = INK2)
    text(rng[2] * 0.97, yy + 0.3, paste0("P = ", fp(x$pJ)), cex = 0.5, col = MUTED, adj = 1)
  }
}
save_fig("Figure3", 180, 82, fig3)

# ---------------------------------------------------------------------------------------------- Figure 4
fig4 <- function() {
  layout(matrix(c(1, 2), 1), widths = c(1, 1.45))
  # A: haplotype matrix, 1000G EUR, grouped by lineage and sub-lineage
  vcf <- file.path(DATA, "1000g", "1kgp_per3_gene.vcf.gz")
  smp <- read_tsv(file.path(DATA, "1000g", "1kgp_samples.tsv"))
  eur <- smp$sample[smp$unrelated_2504 == "yes" & smp$superpopulation == "EUR"]
  m <- haplotype_matrix(vcf, "chr1:7780000-7850000", eur)
  g1 <- function(pos, ref, alt) m$h[which(m$info$pos == pos & m$info$ref == ref & m$info$alt == alt), ]
  five <- g1(7830204, "A", "G") == 0
  sub <- ifelse(five, ifelse(g1(7830004, "A", "G") == 1, "5A", "5B"),
                ifelse(g1(7832049, "T", "A") == 1, ifelse(g1(7824465, "G", "T") == 1, "4B.1", ifelse(g1(7827519, "C", "G") == 1, "4B.2", "4B.0")), "4A"))
  cols <- data.frame(
    lab = c("rs2859387", "rs228663", "rs2640909", "rs10462020", "rs61773390", "rs150812083", "rs228697", "rs228654", "rs148515378", "rs1773135", "rs6680585", "rs12023156"),
    pos = c(7827188, 7832049, 7830057, 7820623, 7824465, 7809893, 7827519, 7837168, 7814731, 7830204, 7805003, 7830004),
    ref = c("G", "T", "T", "T", "G", "C", "C", "G", "C", "A", "G", "A"),
    alt = c("A", "A", "C", "G", "T", "G", "G", "A", "A", "G", "A", "G"),
    show = c("A", "A", "C", "G", "T", "G", "G", "A", "A", "A", "A", "G"),
    note = c("4A", "4B", "4B  M1037T", "4B.1  V647G", "4B.1", "4B.1 (rare)", "4B.2  P864A", "4B.2", "4B.2 (pQTL)", "five units", "5A", "5A"),
    stringsAsFactors = FALSE)
  M <- sapply(seq_len(nrow(cols)), function(j) { x <- g1(cols$pos[j], cols$ref[j], cols$alt[j]); if (cols$show[j] == cols$alt[j]) x else 1 - x })
  ord_l <- c("4A", "4B.1", "4B.2", "4B.0", "5A", "5B")
  o <- order(match(sub, ord_l), M[, 6], M[, 9])
  M <- M[o, ]; sub <- sub[o]
  par(mar = c(5.2, 3.2, 6.4, 0.5))
  nh <- nrow(M)
  plot(NA, xlim = c(-1.6, ncol(M)), ylim = c(nh, 0), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  lab("A")
  lc <- c(LINEAGE_COL[c("4A")], "4B.1" = LINEAGE_COL[["4B"]], "4B.2" = LINEAGE_COL[["4B"]], "4B.0" = LINEAGE_COL[["4B"]], LINEAGE_COL[c("5A", "5B")])
  runs <- function(x) { r <- rle(x); e <- cumsum(r$lengths); data.frame(v = r$values, from = e - r$lengths, to = e) }
  rs_ <- runs(sub); rect(-1.5, rs_$from, -0.9, rs_$to, col = lc[rs_$v], border = NA)
  for (i in seq_len(nrow(rs_))) {
    text(-1.6, (rs_$from[i] + rs_$to[i]) / 2, rs_$v[i], adj = 1, cex = 0.6, xpd = NA)
    if (rs_$from[i] > 0) segments(-1.5, rs_$from[i], ncol(M), rs_$from[i], col = INK2, lwd = 0.5)
  }
  for (j in seq_len(ncol(M))) {
    r <- runs(M[, j]); r <- r[r$v == 1, ]
    if (nrow(r)) rect(j - 0.92, r$from, j - 0.08, r$to, col = INK, border = NA)
  }
  rect(0, 0, ncol(M), nh, border = AXIS, lwd = 0.5)
  text(seq_len(ncol(M)) - 0.5, nh + 10, paste0(cols$lab, "-", cols$show), srt = 90, adj = 1, cex = 0.55, xpd = NA)
  text(seq_len(ncol(M)) - 0.5, -10, cols$note, srt = 90, adj = 0, cex = 0.55, col = INK2, xpd = NA)
  mtext("1,006 European haplotypes (1000 Genomes)", side = 2, line = 1.6, cex = 0.6, col = INK2)
  # B: regional association plot, chronotype
  g <- read.delim(file.path(WORK, "gwas_harmonised_chronotype.tsv"), stringsAsFactors = FALSE)
  ld <- read_tsv(file.path(TAB, "hprc_ld_vntr_snvs.tsv"))
  # lineage that each variant marks: r2 in 1000G EUR with 4B.1, 4B.2 and repeat length
  db <- read.delim(text = bcf_query(file.path(DATA, "reference", "dbsnp157_per3_500kb.vcf.gz"), "%POS\t%ID\n", region = "chr1:7740000-7920000"),
                   header = FALSE, col.names = c("pos", "rsid"))
  mm <- haplotype_matrix(vcf, "chr1:7780000-7850000", eur)
  ok <- rowSums(is.na(mm$h)) == 0 & !grepl(",", mm$info$alt)
  hh <- mm$h[ok, ]; inf <- mm$info[ok, ]
  r2s <- function(y) as.vector(suppressWarnings(cor(t(hh), y)))^2
  r41 <- r2s(as.integer(sub[order(o)] == "4B.1")); r42 <- r2s(as.integer(sub[order(o)] == "4B.2")); r5 <- r2s(as.integer(five))
  inf$rsid <- db$rsid[match(inf$pos, db$pos)]
  cls <- ifelse(r41 >= 0.8, "4B.1", ifelse(r42 >= 0.8, "4B.2", ifelse(r5 >= 0.8, "five", "other")))
  g$cls <- cls[match(g$rsid, inf$rsid)]; g$cls[is.na(g$cls)] <- "other"
  g$pos38 <- g$bp37 - 60060
  g <- g[g$pos38 > 7780000 & g$pos38 < 7850000, ]
  par(mar = c(6.4, 3.6, 6.4, 0.8))
  plot(NA, xlim = c(7.78, 7.85), ylim = c(-5, 26), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  lab("B")
  grid_h(seq(0, 25, 5))
  abline(h = -log10(5e-8), col = MUTED, lwd = 0.6, lty = 2)
  pchs <- c("4B.1" = 24, "4B.2" = 25, five = 23, other = 21)
  for (cl in c("other", "five", "4B.2", "4B.1")) {
    k <- g$cls == cl
    points(g$pos38[k] / 1e6, -log10(g$p[k]), pch = pchs[cl], bg = SUB_COL[cl], col = "white", lwd = 0.5, cex = ifelse(cl == "other", 0.6, 0.95))
  }
  lbl <- data.frame(s = c("rs61773390", "rs228697", "rs150812083", "rs1773135"), dx = c(0.0012, 0.0012, 0.0009, 0.0012),
                    dy = c(1.4, 1.6, -1.4, 1.6), adj = c(0, 0, 0, 0))
  for (i in seq_len(nrow(lbl))) {
    k <- which(g$rsid == lbl$s[i])
    if (length(k)) text(g$pos38[k] / 1e6 + lbl$dx[i], -log10(g$p[k]) + lbl$dy[i], lbl$s[i], adj = lbl$adj[i], cex = 0.58, col = INK)
  }
  axis_style(1, at = seq(7.78, 7.85, 0.01)); axis_style(2, at = seq(0, 25, 5), las = 1)
  mtext("Position on chromosome 1 (Mb, GRCh38)", side = 1, line = 1.7, cex = 0.62, col = INK2)
  mtext(expression(-log[10] * italic(P) * ", chronotype (UK Biobank, n = 449,734)"), side = 2, line = 2.2, cex = 0.62, col = INK2)
  ex <- read.delim(file.path(DATA, "reference", "PER3_exons.bed"), header = FALSE)
  segments(min(ex$V2) / 1e6, -3, max(ex$V3) / 1e6, -3, col = INK2, lwd = 0.8)
  rect(ex$V2 / 1e6, -3.7, ex$V3 / 1e6, -2.3, col = INK2, border = NA)
  text(min(ex$V2) / 1e6 - 0.0008, -3, "PER3", adj = 1, cex = 0.6, font = 3, col = INK2)
  legend(x = 7.78, y = 35.5, legend = c("marks 4B.1", "marks 4B.2", "marks repeat length", "other"), pch = pchs,
         pt.bg = SUB_COL, col = "white", bty = "n", cex = 0.6, pt.cex = 1, ncol = 2, xpd = NA, x.intersp = 0.6,
         text.width = 0.022)
  text(7.78, 29.2, expression("variants with " * italic(r)^2 >= 0.8 * " with the sub-lineage or with repeat length (1000 Genomes EUR)"),
       adj = 0, cex = 0.52, col = MUTED, xpd = NA)
}
save_fig("Figure4", 180, 120, fig4)

# ---------------------------------------------------------------------------------------------- Figure 5
fig5 <- function() {
  d <- read_tsv(file.path(TAB, "primate_units_dna.tsv"))
  grp_order <- c("Hominidae (great apes)", "Hylobatidae (gibbons)", "Cercopithecidae (Old World monkeys)", "Platyrrhini (New World monkeys)")
  grp_lab <- c("Great apes and human", "Gibbons", "Old World monkeys", "New World monkeys")
  d <- d[d$group %in% grp_order, ]
  sp <- unique(d[, c("species", "common_name", "group")])
  first <- c("Homo sapiens", "Pan troglodytes", "Pan paniscus", "Gorilla gorilla gorilla", "Pongo abelii", "Pongo pygmaeus")
  sp <- sp[order(match(sp$group, grp_order), ifelse(sp$species %in% first, match(sp$species, first), 99), sp$species), ]
  gi <- match(sp$group, grp_order)
  pos <- seq_len(nrow(sp)) + 1.4 * (gi - 1) + 1.1 * gi          # room for a heading above each group
  y <- max(pos) + 1 - pos
  par(mar = c(3.6, 15.5, 0.8, 4))
  plot(NA, xlim = c(0.5, 7.5), ylim = c(0.3, max(y) + 0.6), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  grid_v(1:7)
  for (i in seq_len(nrow(sp))) {
    x <- d[d$species == sp$species[i], ]
    filled <- x$long_read_assembly == "yes" | x$source == "GRCh38"
    off <- if (nrow(x) > 1) seq(-0.17, 0.17, length.out = nrow(x)) else 0
    points(x$units, y[i] + off, pch = 21, bg = ifelse(filled, INK, "white"), col = INK, cex = 0.95, lwd = 0.9)
  }
  labs <- ifelse(sp$common_name == "", sp$species, paste0(sp$common_name, " (", sp$species, ")"))
  labs <- sub("^Rhesus", "rhesus", sub("Colobus\\)", "colobus)", labs))
  axis(2, at = y, labels = labs, las = 1, tick = FALSE, cex.axis = 0.58, col.axis = INK, mgp = c(2, 0.2, 0))
  for (g in seq_along(grp_order)) { k <- which(gi == g); text(grconvertX(0.01, "nfc", "user"), max(y[k]) + 1.05, grp_lab[g], adj = 0, cex = 0.66, font = 2, col = INK, xpd = NA) }
  axis_style(1, at = 1:7)
  mtext("Units (one dot per haplotype)", side = 1, line = 1.6, cex = 0.6, col = INK2)
  hy <- y[sp$species == "Homo sapiens"]
  text(5.55, hy, "(4 or 5 in the HPRC)", adj = 0, cex = 0.5, col = MUTED, xpd = NA)
}
save_fig("Figure5", 120, 165, fig5)
cat("15: figures written to", FIG, "\n")
