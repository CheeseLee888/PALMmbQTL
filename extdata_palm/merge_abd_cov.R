rm(list = ls())

abd <- read.table("input/pheno_stool_norrf_count.txt", header = TRUE, row.names = 1)
cov <- read.table("input/covar_stool.txt", header = TRUE, row.names = 1)
# cov$SeqDepth <- NULL
# covar <- covar[,c(1:7,14)]

abd$IID <- rownames(abd)
rownames(abd) <- NULL
abd <- abd[, c("IID", setdiff(names(abd), "IID"))]

cov$IID <- rownames(cov)
rownames(cov) <- NULL
cov <- cov[, c("IID", setdiff(names(cov), "IID"))]

merged <- merge(cov, abd, by = "IID", sort = FALSE)


covarCols <- c(
  "PC01","PC02","PC03","PC04","PC05",
  "SAMPLE_COLLECTION_AGE_MONTHS_DEV","SEX","CSECT","EVERBREASTFED",
  "EARLY_CATS_DOGS_DEV","ANTIBIOTICSEXPOSURE_DEV",
  "BIRTH_SEASON_DEV_Summer","BIRTH_SEASON_DEV_Autumn","BIRTH_SEASON_DEV_Winter"
)

# 强制转 numeric（防止 factor/character/logical）
# for (cc in covarCols) if (cc %in% names(merged)) merged[[cc]] <- as.numeric(merged[[cc]])

# merged2 <- merged[rep(seq_len(nrow(merged)), each = 10), ]
# write.table(merged, "input/nasal_bialleic_cov_pheno_merged.txt", sep = "\t", quote = FALSE, row.names = FALSE)


write.table(merged, "input/stool_bialleic_cov_pheno_merged.txt", sep = "\t", quote = FALSE, row.names = FALSE)

