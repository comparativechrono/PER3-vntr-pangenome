# Shared R helpers for the PER3 VNTR analysis (base R only).
args_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(f), ".."))
}
ROOT <- Sys.getenv("PER3_ROOT"); if (ROOT == "") ROOT <- args_root()
DATA <- file.path(ROOT, "data")
TAB <- file.path(ROOT, "outputs", "tables")
FIG <- file.path(ROOT, "outputs", "figures")
WORK <- file.path(ROOT, "outputs", "work")
for (d in c(TAB, FIG, WORK)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

read_tsv <- function(path, ...) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
write_tsv <- function(x, name) write.table(x, file.path(TAB, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Wilson score interval for a proportion
wilson <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  cbind(lower = centre - half, upper = centre + half)
}

# Exact test of Hardy-Weinberg equilibrium (Wigginton, Cutler & Abecasis 2005)
hwe_exact <- function(n_hom1, n_het, n_hom2) {
  n_hom_r <- min(n_hom1, n_hom2); n_hom_c <- max(n_hom1, n_hom2)
  rare <- 2 * n_hom_r + n_het; N <- n_hom1 + n_het + n_hom2
  probs <- numeric(rare + 1)
  mid <- floor(rare * (2 * N - rare) / (2 * N)); if ((mid %% 2) != (rare %% 2)) mid <- mid + 1
  probs[mid + 1] <- 1; curr_hom_r <- (rare - mid) / 2; curr_hom_c <- N - mid - curr_hom_r
  h <- mid; hr <- curr_hom_r; hc <- curr_hom_c
  while (h >= 2) {
    probs[h - 1] <- probs[h + 1] * h * (h - 1) / (4 * (hr + 1) * (hc + 1)); h <- h - 2; hr <- hr + 1; hc <- hc + 1
  }
  h <- mid; hr <- curr_hom_r; hc <- curr_hom_c
  while (h <= rare - 2) {
    probs[h + 3] <- probs[h + 1] * 4 * hr * hc / ((h + 2) * (h + 1)); h <- h + 2; hr <- hr - 1; hc <- hc - 1
  }
  probs <- probs / sum(probs)
  p <- sum(probs[probs <= probs[n_het + 1] * (1 + 1e-8)])
  min(1, p)
}

fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(signif(p, 2), format = "fg", digits = 2, flag = "#")))

# bcftools query into a data frame
bcf_query <- function(vcf, fmt, region = NULL, samples = NULL, extra = character()) {
  cmd <- c("query", "-f", shQuote(fmt), if (!is.null(region)) c("-r", region),
           if (!is.null(samples)) c("-s", paste(samples, collapse = ",")), extra, vcf)
  x <- system2("bcftools", cmd, stdout = TRUE)
  x
}

# phased haplotype matrix (variants x haplotypes, 0/1, NA for missing) from a VCF region
haplotype_matrix <- function(vcf, region, samples, biallelic_snv_only = FALSE) {
  tmp <- tempfile(fileext = ".txt")
  s <- tempfile(fileext = ".txt"); writeLines(samples, s)
  extra <- c("-S", s)
  if (biallelic_snv_only) extra <- c(extra, "-i", shQuote('TYPE="snp" && N_ALT=1'))
  system2("bcftools", c("query", "-f", shQuote("%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"), "-r", region, extra, vcf),
          stdout = tmp)
  d <- read.delim(tmp, header = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  info <- data.frame(chrom = d[[1]], pos = as.integer(d[[2]]), id = d[[3]], ref = d[[4]], alt = d[[5]],
                     stringsAsFactors = FALSE)
  g <- as.matrix(d[, -(1:5), drop = FALSE])
  a1 <- substr(g, 1, 1); a2 <- substr(g, 3, 3)
  h <- matrix(NA_integer_, nrow(g), 2 * ncol(g))
  h[, seq(1, 2 * ncol(g), 2)] <- suppressWarnings(as.integer(a1))
  h[, seq(2, 2 * ncol(g), 2)] <- suppressWarnings(as.integer(a2))
  colnames(h) <- paste0(rep(samples, each = 2), "_hap", 1:2)
  list(info = info, h = h)
}

r2 <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA); cor(x[ok], y[ok])^2 }
rr <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA); cor(x[ok], y[ok]) }
dprime <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  pA <- mean(x); pB <- mean(y); D <- mean(x * y) - pA * pB
  Dmax <- if (D > 0) min(pA * (1 - pB), (1 - pA) * pB) else min(pA * pB, (1 - pA) * (1 - pB))
  if (Dmax == 0) NA else D / Dmax
}

GROUPS <- c("AFR", "AMR", "EAS", "EUR", "SAS")
GROUP_NAMES <- c(AFR = "African", AMR = "American", EAS = "East Asian", EUR = "European", SAS = "South Asian")
# Lineage colours: four-unit lineages cool (blue, aqua), five-unit lineages warm (yellow, red); validated as an
# adjacent-pair categorical palette in this order (CVD dE >= 9.1, normal-vision dE >= 20.8); 'other' is neutral grey.
LINEAGE_COL <- c("4A" = "#2a78d6", "4B" = "#1baf7a", "5A" = "#eda100", "5B" = "#e34948", other = "#898781")
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
SEQ_BLUE <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
