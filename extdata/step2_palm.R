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
geno <- as(G, "numeric")   # returns matrix with 0/1/2 and NA
rownames(geno) <- iid
colnames(geno) <- colnames(G)


geno<-geno[,1:3]

res <- palm.get.summary(
  null.obj = modglmm,
  covariate.interest = geno,
  correct = opt$correct
)

# ---------- split by pheno and write {pheno}_step2_palm.txt ----------
res <- as.data.frame(res, check.names = FALSE)

# Automatically infer the study prefix (usually "Study")
prefix <- sub("\\.est\\..*$", "", grep("\\.est\\.", colnames(res), value = TRUE)[1])
if (is.na(prefix) || prefix == "") prefix <- "Study"

est_pat    <- paste0("^", prefix, "\\.est\\.")
stderr_pat <- paste0("^", prefix, "\\.stderr\\.")

est_cols    <- grep(est_pat, colnames(res), value = TRUE)
stderr_cols <- grep(stderr_pat, colnames(res), value = TRUE)

if (length(est_cols) == 0 || length(stderr_cols) == 0) {
  stop("Cannot find est/stderr columns in res. Example colnames(res): ",
       paste(head(colnames(res), 5), collapse = ", "))
}

snp_est    <- sub(est_pat,    "", est_cols)
snp_stderr <- sub(stderr_pat, "", stderr_cols)
common_snp <- intersect(snp_est, snp_stderr)
if (length(common_snp) == 0) stop("No matched SNPs between est and stderr columns.")

# keep SNP order as in est columns
common_snp <- snp_est[snp_est %in% common_snp]
est_map    <- setNames(est_cols,    snp_est)
stderr_map <- setNames(stderr_cols, snp_stderr)

if (is.null(rownames(res)) || any(rownames(res) == "")) {
  stop("res has no rownames (phenotype names). Please ensure rownames(res)=pheno names.")
}

# opt$PALMOutputFile can be a "directory" or "prefix"; here we treat it as a directory for clarity
out_dir <- opt$PALMOutputFile
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (pheno in rownames(res)) {
  out <- data.frame(
    SNP    = common_snp,
    est    = as.numeric(res[pheno, est_map[common_snp], drop = TRUE]),
    stderr = as.numeric(res[pheno, stderr_map[common_snp], drop = TRUE]),
    check.names = FALSE
  )

  out_file <- file.path(out_dir, paste0(pheno, "_step2_palm.txt"))
  write.table(out, file = out_file, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}

cat("Wrote per-pheno files to: ", out_dir, "\n")



