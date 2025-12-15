#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(PALM)
  library(snpStats)
})

option_list <- list(
    make_option("--inFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--correct",
        type = "character", default = "NULL",
        help = ""
    ),
    make_option("--NULLmodelFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--PALMOutputFile",
        type = "character", default = "",
        help = ""
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

load(opt$NULLmodelFile)  # load modglmm

# read genotype data and make it a data.frame
bed <- paste0(opt$inFile, ".bed")
bim <- paste0(opt$inFile, ".bim")
fam <- paste0(opt$inFile, ".fam")
for (f in c(bed, bim, fam)) if (!file.exists(f)) stop("Missing PLINK file: ", f)

plink <- snpStats::read.plink(bed, bim, fam)

G <- plink$genotypes

iid <- rownames(G)
if (is.null(iid)) stop("No rownames (IID) found in genotype matrix from read.plink().")

# Convert to numeric 0/1/2/NA matrix, then to data.frame
geno_mat <- as(G, "numeric")   # returns matrix with 0/1/2 and NA
rownames(geno_mat) <- iid
colnames(geno_mat) <- colnames(G)

geno <- as.matrix(as.data.frame(geno_mat, check.names = FALSE))


res <- palm.get.summary(
  null.obj = modglmm,
  covariate.interest = geno,
  correct = opt$correct
)

write.table(
  res,
  file = opt$PALMOutputFile,
  sep = "\t",
  row.names = TRUE
)



