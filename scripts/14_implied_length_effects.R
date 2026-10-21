# 14. What association of repeat length (five vs four units) with chronotype would the lineage effects alone
#     produce in each population? A candidate-gene study that compares four- and five-unit alleles compares
#     mixtures of lineages, and the mixtures differ between populations. With the joint lineage effects estimated
#     in UK Biobank (step 09, EUR LD panel; 4A = 0) and the lineage composition of each 1000 Genomes population
#     (step 06 and the 4B sub-lineages), the expected per-allele difference between five- and four-unit haplotypes
#     is mean(effect | five units) - mean(effect | four units). This assumes that the European estimates apply
#     elsewhere; it is an illustration, not a test.
# Output (outputs/tables): implied_length_effects.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

lm_ <- read_tsv(file.path(TAB, "lineage_model.tsv")); lm_ <- lm_[lm_$panel == "EUR" & lm_$trait == "chronotype", ]
b <- setNames(lm_$bJ, lm_$marker)
eff <- c("4A" = 0, "4B.0" = 0, "4B.1" = b[["rs61773390"]], "4B.1r" = b[["rs61773390"]] + b[["rs150812083"]],
         "4B.2" = b[["rs228697"]], "5B" = b[["rs1773135"]], "5A" = b[["rs1773135"]] + b[["rs6680585"]])
sd_y <- read_tsv(file.path(TAB, "gwas_tag_results.tsv")); sd_y <- sd_y$sd_trait[sd_y$trait == "chronotype"][1]

vcf <- file.path(DATA, "1000g", "1kgp_per3_gene.vcf.gz")
smp <- read_tsv(file.path(DATA, "1000g", "1kgp_samples.tsv"))
s <- system2("bcftools", c("query", "-l", vcf), stdout = TRUE); smp <- smp[match(s, smp$sample), ]
m <- haplotype_matrix(vcf, "chr1:7800000-7835000", s)
g1 <- function(pos, ref, alt) m$h[which(m$info$pos == pos & m$info$ref == ref & m$info$alt == alt), ]
five <- g1(7830204, "A", "G") == 0
lin <- ifelse(five, ifelse(g1(7830004, "A", "G") == 1, "5A", "5B"),
              ifelse(g1(7832049, "T", "A") == 1,
                     ifelse(g1(7824465, "G", "T") == 1, ifelse(g1(7809893, "C", "G") == 1, "4B.1r", "4B.1"),
                            ifelse(g1(7827519, "C", "G") == 1, "4B.2", "4B.0")), "4A"))
unrel <- rep(smp$unrelated_2504 == "yes", each = 2)
pop <- rep(smp$population, each = 2); grp <- rep(smp$superpopulation, each = 2)
calc <- function(k) {
  e <- eff[lin[k]]; f <- five[k]
  comp <- table(factor(lin[k], names(eff))) / sum(k)
  d <- mean(e[f]) - mean(e[!f])
  data.frame(haplotypes = sum(k), freq_five = round(mean(f), 3),
             share_4B1_of_four = round(sum(lin[k] %in% c("4B.1", "4B.1r")) / sum(!f), 3),
             share_4B2_of_four = round(sum(lin[k] == "4B.2") / sum(!f), 3),
             implied_beta_five = signif(d, 3), implied_beta_five_sd = signif(d / sd_y, 3),
             implied_55_vs_44_sd = signif(2 * d / sd_y, 3))
}
rows <- rbind(
  do.call(rbind, lapply(c(GROUPS), function(g) data.frame(level = "group", name = g, calc(unrel & grp == g)))),
  do.call(rbind, lapply(sort(unique(pop[unrel])), function(p) data.frame(level = "population", name = p, calc(unrel & pop == p)))))
write_tsv(rows, "implied_length_effects.tsv")
cat("14: lineage effects used (chronotype score units, 4A = 0):\n"); print(round(eff, 4))
cat(sprintf("14: chronotype SD (from the summary statistics): %.3f\n", sd_y))
print(rows[rows$level == "group", ], row.names = FALSE)
r <- rows[rows$level == "population", ]
cat(sprintf("14: implied per-allele difference over 26 populations: %.4f to %.4f SD (5/5 vs 4/4: %.4f to %.4f SD)\n",
            min(r$implied_beta_five_sd), max(r$implied_beta_five_sd), min(r$implied_55_vs_44_sd), max(r$implied_55_vs_44_sd)))
