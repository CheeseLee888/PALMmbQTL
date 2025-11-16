rm(list = ls())

# args <- commandArgs(trailingOnly = TRUE)
# genus <- as.numeric(args[1]) # 2:g_Escherichia_Shigella, 58: g_Actinomyces, 56: g_Acinetobacter

genus <- "g_Blautia"
library(GMMAT)
library(phyloseq)
library(parallel)
# library(Matrix)
# devtools::install_github("variani/lme4qtl")

base_dir <- "input"
Y <- read.table(file.path(base_dir, "pheno_stool_norrf_count.txt"), header = TRUE, row.names = 1)
Y.g <- Y[, which(startsWith(colnames(Y), "g_"))]

covar <- read.table(file.path(base_dir, "covar_stool.txt"), header = T, row.names = 1)
covar <- covar[,c(1:7,14)]

# rownames(Y.g) <- rownames(covar) <- paste0("s", 1:nrow(Y.g))

# pheno <- cbind(rownames(Y.g), Y.g[,genus], covar)
pheno <- data.frame(
  id    = rownames(Y.g),
  pheno = Y.g[, genus],
  covar,
  check.names = FALSE
)
rownames(pheno) <- NULL
colnames(pheno) <- c("id", "genus", "pc1", "pc2", "pc3", "pc4", "pc5", "age", "sex", "seqdepth")


grm <- diag(1, nrow = nrow(pheno))
colnames(grm) <- rownames(grm) <- pheno$id
# model0_gmmat <- glmmkin(genus ~ offset(log(seqdepth)) + pc1 + pc2 + pc3 + pc4 + pc5 + age + sex, data = pheno, kins = grm,
#                         id = "id", family = poisson(link = "log"))
model0_gmmat <- glmmkin(genus ~ offset(log(seqdepth)) + pc1 + pc2 + pc3 + pc4 + pc5 + age + sex, data = pheno,
                        kins = grm,
                        id = "id", family = poisson(link = "log"))
print("done")
# gdsfile <- file.path(base_dir, "stool_bialleic_merged_data.gds")
# glmm.score(model0_gmmat, infile = gdsfile, center = T, ncores = detectCores() - 1, 
#            outfile = paste0("output/", colnames(Y.g)[genus], ".testoutfile_c.txt"))

plinkfile <- file.path(base_dir, "stool_bialleic_merged_data")
glmm.score(model0_gmmat, infile = plinkfile, center = T, 
           outfile = paste0("output/", colnames(Y.g)[genus], ".plink.txt"))

# nohup Rscript run_PALM_mbQTL_score_stool_nohup.R 2 > mbqtl_2.log 2>&1 &