# Conditional and joint analysis of GWAS summary statistics (Yang et al. 2012, Nat Genet 44:369), base R.
# Follows the model of GCTA-COJO: X'X is approximated from the allele frequencies and sample sizes of the GWAS and
# the LD (genotype correlations) of a reference sample. Checked against GCTA 1.94.1 in 09_conditional.R.

# s: data frame with columns snp, freq (of the effect allele), b, se, p, n; R: LD correlation matrix (same order,
# same effect allele). Returns a list with the matrices needed by the functions below.
cojo_prepare <- function(s, R) {
  stopifnot(nrow(s) == nrow(R))
  vx <- 2 * s$freq * (1 - s$freq)
  D <- vx * s$n
  yy_j <- D * s$se^2 * (s$n - 1) + D * s$b^2              # y'y implied by each SNP
  Vp <- median(yy_j / (s$n - 1))                             # phenotypic variance
  Neff <- (Vp - vx * s$b^2) / (vx * s$se^2) + 1             # effective N of each SNP
  N <- median(Neff)
  Dn <- vx * Neff
  nmin <- outer(Neff, Neff, pmin)
  B <- nmin * sqrt(outer(vx, vx)) * R
  diag(B) <- Dn
  list(s = s, R = R, vx = vx, D = Dn, B = B, Vp = Vp, N = N, yy = Vp * (N - 1))
}

# joint model of the SNPs with indices idx
cojo_joint <- function(cj, idx) {
  B <- cj$B[idx, idx, drop = FALSE]; D <- cj$D[idx]; b <- cj$s$b[idx]
  Bi <- solve(B)
  bJ <- as.vector(Bi %*% (D * b))
  sigma2 <- (cj$yy - sum(bJ * D * b)) / (cj$N - length(idx))
  seJ <- sqrt(diag(Bi) * sigma2)
  out <- data.frame(snp = cj$s$snp[idx], b = cj$s$b[idx], se = cj$s$se[idx], p = cj$s$p[idx],
                    bJ = bJ, bJ_se = seJ, pJ = 2 * pnorm(-abs(bJ / seJ)), stringsAsFactors = FALSE)
  attr(out, "vcov") <- Bi * sigma2                         # covariance matrix of the joint effects
  out
}

# conditional effect of every SNP given the SNPs with indices cond
cojo_cond <- function(cj, cond) {
  B_ss <- cj$B[cond, cond, drop = FALSE]; Bi <- solve(B_ss)
  D <- cj$D; b <- cj$s$b
  bJ <- as.vector(Bi %*% (D[cond] * b[cond]))
  sigma2 <- (cj$yy - sum(bJ * D[cond] * b[cond])) / (cj$N - length(cond) - 1)
  B_js <- cj$B[, cond, drop = FALSE]
  num <- D * b - as.vector(B_js %*% (Bi %*% (D[cond] * b[cond])))
  den <- D - rowSums((B_js %*% Bi) * B_js)
  bC <- num / den
  seC <- sqrt(sigma2 / den)
  out <- data.frame(snp = cj$s$snp, b = b, se = cj$s$se, p = cj$s$p, bC = bC, bC_se = seC,
                    pC = 2 * pnorm(-abs(bC / seC)), r2_with_conditioning = 1 - den / D, stringsAsFactors = FALSE)
  out[cond, c("bC", "bC_se", "pC")] <- NA
  out
}

# stepwise selection of independent signals, as GCTA --cojo-slct (p threshold, collinearity cut-off)
cojo_select <- function(cj, p_thresh = 5e-8, collinear = 0.9, max_snps = 50) {
  sel <- integer(0); dropped <- integer(0)
  repeat {
    if (length(sel) == 0) {
      cand <- data.frame(p = cj$s$p, r2 = 0)
    } else {
      cc <- cojo_cond(cj, sel); cand <- data.frame(p = cc$pC, r2 = cc$r2_with_conditioning)
    }
    cand$p[c(sel, dropped)] <- NA
    ok <- !is.na(cand$p) & cand$p < p_thresh & cand$r2 < collinear
    if (!any(ok) || length(sel) >= max_snps) break
    j <- which(ok)[which.min(cand$p[ok])]
    sel <- c(sel, j)
    # backward step: drop SNPs that are no longer significant in the joint model (and do not add them again)
    repeat {
      jt <- cojo_joint(cj, sel)
      if (all(jt$pJ < p_thresh) || length(sel) == 1) break
      w <- which.max(jt$pJ); dropped <- c(dropped, sel[w]); sel <- sel[-w]
    }
  }
  sel
}
