# 05. Linkage disequilibrium between the VNTR and every biallelic SNV within 250 kb, in the 462 phased long-read
#     haplotypes of the HPRC (decomposed 'wave' VCF of the same graph); and the SNVs that best tag each lineage.
# Outputs (outputs/tables): hprc_ld_vntr_snvs.tsv (all SNVs), hprc_tags_top.tsv, hprc_ld_extent.tsv,
#                           hprc_lineage_tags.tsv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))

a <- read_tsv(file.path(TAB, "hprc_haplotype_alleles.tsv"))
vcf <- file.path(DATA, "hprc", "hprc_per3_500kb.vcf.gz")
samples <- system2("bcftools", c("query", "-l", vcf), stdout = TRUE)
m <- haplotype_matrix(vcf, "chr1:7580000-8080000", samples, biallelic_snv_only = TRUE)
a <- a[match(colnames(m$h), a$haplotype), ]
stopifnot(identical(a$haplotype, colnames(m$h)))
vntr <- as.integer(a$units == 5)
grp <- a$superpopulation

# rsIDs (dbSNP build 157)
db <- read.delim(text = bcf_query(file.path(DATA, "reference", "dbsnp157_per3_500kb.vcf.gz"), "%POS\t%REF\t%ALT\t%ID\n"),
                 header = FALSE, col.names = c("pos", "ref", "alt", "rsid"), stringsAsFactors = FALSE)
alts <- strsplit(db$alt, ","); db_long <- data.frame(pos = rep(db$pos, lengths(alts)), ref = rep(db$ref, lengths(alts)),
                                                      alt = unlist(alts), rsid = rep(db$rsid, lengths(alts)))
key_db <- paste(db_long$pos, db_long$ref, db_long$alt)
rs <- db_long$rsid[match(paste(m$info$pos, m$info$ref, m$info$alt), key_db)]

h <- m$h
complete <- rowSums(is.na(h)) == 0
poly <- rowSums(h, na.rm = TRUE) > 0 & rowSums(h, na.rm = TRUE) < rowSums(!is.na(h))
keep <- complete & poly
cat(sprintf("05: %d biallelic SNVs; %d complete and polymorphic in the 462 haplotypes\n", nrow(h), sum(keep)))
h <- h[keep, , drop = FALSE]; info <- m$info[keep, ]; rs <- rs[keep]

r_all <- as.vector(cor(t(h), vntr))
res <- data.frame(pos = info$pos, rsid = rs, ref = info$ref, alt = info$alt, alt_freq = round(rowMeans(h), 4),
                  r = round(r_all, 4), r2 = round(r_all^2, 4))
res$dprime <- round(apply(h, 1, dprime, y = vntr), 3)
for (g in GROUPS) {
  k <- grp == g
  rg <- suppressWarnings(as.vector(cor(t(h[, k]), vntr[k])))
  res[[paste0("r2_", g)]] <- round(rg^2, 4)
}
res$distance_kb <- round((res$pos - mean(c(7829888, 7830157))) / 1000, 2)
res$five_unit_allele <- ifelse(res$r > 0, res$alt, res$ref)
write_tsv(res, "hprc_ld_vntr_snvs.tsv")
top <- res[order(-res$r2), ][1:25, ]
write_tsv(top, "hprc_tags_top.tsv")

# extent of LD: the farthest SNV (either side) with r2 above each threshold
ext <- do.call(rbind, lapply(c(0.95, 0.8, 0.5, 0.2), function(t) {
  x <- res[res$r2 > t, ]
  data.frame(r2_above = t, snvs = nrow(x), from = min(x$pos), to = max(x$pos),
             span_kb = round((max(x$pos) - min(x$pos)) / 1000, 1),
             farthest_kb = round(max(abs(x$distance_kb)), 1))
}))
write_tsv(ext, "hprc_ld_extent.tsv")

# best SNV tags of each lineage (indicator of the lineage vs all other haplotypes)
lt <- do.call(rbind, lapply(c("4A", "4B", "5A", "5B"), function(l) {
  y <- as.integer(a$lineage == l)
  rr_ <- as.vector(cor(t(h), y))
  o <- order(-rr_^2)[1:5]
  data.frame(lineage = l, haplotypes = sum(y), rank = 1:5, pos = info$pos[o], rsid = rs[o], ref = info$ref[o],
             alt = info$alt[o], r2 = round(rr_[o]^2, 3),
             lineage_allele = ifelse(rr_[o] > 0, info$alt[o], info$ref[o]))
}))
write_tsv(lt, "hprc_lineage_tags.tsv")

cat("05: the best tags of repeat length (HPRC):\n")
print(head(top[, c("pos", "rsid", "ref", "alt", "five_unit_allele", "alt_freq", "r2", "dprime", paste0("r2_", GROUPS))], 8), row.names = FALSE)
cat("05: extent of LD with the repeat:\n"); print(ext, row.names = FALSE)
cat("05: best tag of each lineage:\n"); print(lt[lt$rank == 1, ], row.names = FALSE)
