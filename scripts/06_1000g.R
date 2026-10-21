# 06. The repeat in 3,202 short-read genomes (1000 Genomes high coverage, phased), through its tag SNPs.
#   - repeat length from rs1773135 (A = five units; checked against rs228701)
#   - lineages: four-unit haplotypes are 4B if they carry rs228663-A, else 4A; five-unit haplotypes are 5A if they
#     carry rs12023156-G, else 5B (the best tags of each lineage in the HPRC, step 05)
#   - checks against the long-read alleles of the people in both data sets
#   - what the short-read structural-variant calls make of the repeat
# Outputs (outputs/tables): kg_freq_by_population.tsv, kg_freq_by_group.tsv, kg_lineages_by_group.tsv,
#   kg_lineages_by_population.tsv, kg_validation.tsv, kg_validation_discordant.tsv, kg_sv_calls.tsv, kg_tag_ld.tsv,
#   kg_sv_carriers.tsv, kg_sv_by_group.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

vcf <- file.path(DATA, "1000g", "1kgp_per3_gene.vcf.gz")
smp <- read_tsv(file.path(DATA, "1000g", "1kgp_samples.tsv"))
samples <- system2("bcftools", c("query", "-l", vcf), stdout = TRUE)
smp <- smp[match(samples, smp$sample), ]
stopifnot(identical(smp$sample, samples))
m <- haplotype_matrix(vcf, "chr1:7829700-7832100", samples)
get <- function(pos, ref, alt) {
  i <- which(m$info$pos == pos & m$info$ref == ref & m$info$alt == alt)
  stopifnot(length(i) == 1); m$h[i, ]
}
rs228701 <- get(7829744, "C", "G")      # C = five units
rs1773135 <- get(7830204, "A", "G")     # A = five units
rs228663 <- get(7832049, "T", "A")      # A = 4B
rs12023156 <- get(7830004, "A", "G")    # G = 5A
rs2640909 <- get(7830057, "T", "C")     # C = Thr1037, inside unit R4 of 4B
hap_sample <- rep(samples, each = 2)
hap_group <- rep(smp$superpopulation, each = 2)
hap_pop <- rep(smp$population, each = 2)
unrel <- rep(smp$unrelated_2504 == "yes", each = 2)
stopifnot(!anyNA(rs228701), !anyNA(rs1773135), !anyNA(rs228663), !anyNA(rs12023156))

five <- as.integer(rs1773135 == 0)
lineage <- ifelse(five == 1, ifelse(rs12023156 == 1, "5A", "5B"), ifelse(rs228663 == 1, "4B", "4A"))

# the two tags agree?
tag_ld <- do.call(rbind, lapply(c("All", GROUPS), function(g) {
  k <- unrel & (g == "All" | hap_group == g)
  data.frame(group = g, haplotypes = sum(k), r2_rs228701_rs1773135 = round(r2(rs228701[k], rs1773135[k]), 4),
             discordant_haplotypes = sum((rs228701[k] == 0) != (rs1773135[k] == 0)),
             r2_rs2640909_rs228663 = round(r2(rs2640909[k], rs228663[k]), 3))
}))
write_tsv(tag_ld, "kg_tag_ld.tsv")

# frequencies in the 2,504 unrelated people
fr <- function(k) { ci <- wilson(sum(five[k]), sum(k)); c(haplotypes = sum(k), five = sum(five[k]), freq_five = round(mean(five[k]), 3), ci_lower = round(ci[1], 3), ci_upper = round(ci[2], 3)) }
byg <- do.call(rbind, lapply(c(GROUPS, "All"), function(g) data.frame(group = g, t(fr(unrel & (g == "All" | hap_group == g))))))
write_tsv(byg, "kg_freq_by_group.tsv")
pops <- unique(smp[smp$unrelated_2504 == "yes", c("population", "population_name", "superpopulation")])
byp <- do.call(rbind, lapply(seq_len(nrow(pops)), function(i) data.frame(pops[i, ], t(fr(unrel & hap_pop == pops$population[i])))))
byp <- byp[order(byp$superpopulation, byp$freq_five), ]
write_tsv(byp, "kg_freq_by_population.tsv")

lin_tab <- function(k) { t <- table(factor(lineage[k], c("4A", "4B", "5A", "5B"))); c(haplotypes = sum(k), round(t / sum(k), 3)) }
lg <- do.call(rbind, lapply(c(GROUPS, "All"), function(g) data.frame(group = g, t(lin_tab(unrel & (g == "All" | hap_group == g))), check.names = FALSE)))
write_tsv(lg, "kg_lineages_by_group.tsv")
lp <- do.call(rbind, lapply(seq_len(nrow(pops)), function(i) data.frame(pops[i, ], t(lin_tab(unrel & hap_pop == pops$population[i])), check.names = FALSE)))
write_tsv(lp[order(lp$superpopulation, lp$population), ], "kg_lineages_by_population.tsv")

# validation against the long reads (people in both data sets; genotypes compared without phase)
hp <- read_tsv(file.path(TAB, "hprc_haplotype_alleles.tsv"))
both <- intersect(unique(hp$sample), samples)
pair <- function(x) paste(sort(x), collapse = "/")
val <- do.call(rbind, lapply(both, function(s) {
  i <- which(hap_sample == s); j <- which(hp$sample == s)
  data.frame(sample = s, group = hp$superpopulation[j[1]],
             units_short = pair(4 + five[i]), units_long = pair(hp$units[j]),
             lineage_short = pair(lineage[i]), lineage_long = pair(hp$lineage[j]))
}))
val$units_agree <- val$units_short == val$units_long
val$lineage_agree <- val$lineage_short == val$lineage_long
val$lineage_agree_excl_other <- ifelse(grepl("other", val$lineage_long), NA, val$lineage_agree)
write_tsv(val[!val$units_agree | !val$lineage_agree, ], "kg_validation_discordant.tsv")
summ <- data.frame(people = nrow(val), units_concordant = sum(val$units_agree),
                   lineage_concordant = sum(val$lineage_agree),
                   people_without_other_lineage = sum(!is.na(val$lineage_agree_excl_other)),
                   lineage_concordant_without_other = sum(val$lineage_agree_excl_other, na.rm = TRUE))
write_tsv(summ, "kg_validation.tsv")

# the structural-variant calls of the short-read data: 54-bp deletions at the repeat
dlen <- nchar(m$info$ref) - nchar(m$info$alt)
sv <- m$info$pos >= 7829834 & m$info$pos <= 7830161 & dlen >= 50
svc <- data.frame(m$info[sv, c("pos", "ref", "alt")], deleted_bp = dlen[sv],
                  carriers_haplotypes = rowSums(m$h[sv, , drop = FALSE], na.rm = TRUE), haplotypes = ncol(m$h))
svc$ref <- paste0(substr(svc$ref, 1, 12), "... (", nchar(svc$ref), " bp)")
any_del <- colSums(m$h[sv, , drop = FALSE], na.rm = TRUE) > 0
svc <- rbind(svc, data.frame(pos = NA, ref = "any of these deletions", alt = "", deleted_bp = NA,
                             carriers_haplotypes = sum(any_del), haplotypes = ncol(m$h)))
svc$four_unit_by_tag <- sum(five == 0)
write_tsv(svc, "kg_sv_calls.tsv")
# each haplotype with a deletion call: where it is, what the tags say, and the structure that the panel's variants
# give on that haplotype (rebuilt in step 01, units read in step 02)
sva <- read_tsv(file.path(TAB, "kg_sv_alleles.tsv"))
k <- which(any_del)
car <- data.frame(haplotype = paste0(hap_sample[k], "_hap", ifelse(k %% 2 == 1, 1, 2)), sample = hap_sample[k],
                  population = hap_pop[k], superpopulation = hap_group[k], unrelated = ifelse(unrel[k], "yes", "no"),
                  deletion_at = sapply(k, function(j) paste(m$info$pos[sv][which(m$h[sv, j] == 1)], collapse = ",")),
                  rs228701 = ifelse(rs228701[k] == 0, "C", "G"), rs1773135 = ifelse(rs1773135[k] == 0, "A", "G"),
                  units_by_tag = 4 + five[k], stringsAsFactors = FALSE)
car$units_short_read <- sva$units[match(car$haplotype, sva$haplotype)]
car$structure_short_read <- sva$structure[match(car$haplotype, sva$haplotype)]
car$long_read_alleles <- sapply(car$sample, function(s) if (s %in% hp$sample) paste(sort(hp$structure[hp$sample == s]), collapse = " / ") else "")
write_tsv(car, "kg_sv_carriers.tsv")
svg <- do.call(rbind, lapply(c("All", GROUPS), function(g) {
  kk <- if (g == "All") rep(TRUE, length(five)) else hap_group == g
  data.frame(group = g, haplotypes = sum(kk), with_deletion_call = sum(any_del & kk),
             on_five_unit_tag = sum(any_del & kk & five == 1), on_four_unit_tag = sum(any_del & kk & five == 0),
             four_unit_by_tag = sum(kk & five == 0), freq_deletion_on_five_unit_tag = round(sum(any_del & kk & five == 1) / sum(kk), 4))
}))
write_tsv(svg, "kg_sv_by_group.tsv")

cat("06: the two tags in 1000G (2,504 unrelated):\n"); print(tag_ld, row.names = FALSE)
cat("06: five-unit frequency by group (1000G, via rs1773135):\n"); print(byg, row.names = FALSE)
cat(sprintf("06: range over 26 populations: %s %.3f to %s %.3f\n", byp$population[which.min(byp$freq_five)], min(byp$freq_five),
            byp$population[which.max(byp$freq_five)], max(byp$freq_five)))
cat("06: lineages by group (1000G):\n"); print(lg, row.names = FALSE)
cat("06: validation against long reads:\n"); print(summ, row.names = FALSE); print(val[!val$units_agree | !val$lineage_agree, ], row.names = FALSE)
cat("06: SV calls:\n"); print(svc, row.names = FALSE)
cat("06: haplotypes with a deletion call, by group and tag:\n"); print(svg, row.names = FALSE)
cat("06: structures that the calls give:\n"); print(table(car$structure_short_read, car$units_by_tag))
