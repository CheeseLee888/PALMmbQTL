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
    ),
    make_option("--chrom",
        type = "character", default = "",
        help = ""
    ),
    make_option("--cluster",
        type = "character", default = NULL,
        help = "")
)

opt <- parse_args(OptionParser(option_list = option_list))

load(opt$NULLmodelFile)  # load modglmm

# read genotype data and make it a data.frame
bed <- paste0(opt$inFile, ".bed")
bim <- paste0(opt$inFile, ".bim")
fam <- paste0(opt$inFile, ".fam")
for (f in c(bed, bim, fam)) if (!file.exists(f)) stop("Missing PLINK file: ", f)

# read fam to get cluster info (make sure IDs align with abdFile's)
cat("Reading PLINK .fam file for cluster info. \n")
fam_data <- read.table(fam, stringsAsFactors = FALSE)
colnames(fam_data) <- c("FID","IID","PID","MID","SEX","PHENO")

cluster <- fam_data$FID
names(cluster) <- fam_data$IID
# ----------------------------

plink <- snpStats::read.plink(bed, bim, fam)

G <- plink$genotypes

iid <- rownames(G)
if (is.null(iid)) stop("No rownames (IID) found in genotype matrix from read.plink().")

# Convert to numeric 0/1/2/NA matrix, then to data.frame
geno <- as(G, "numeric")   # returns matrix with 0/1/2 and NA
rownames(geno) <- iid
colnames(geno) <- colnames(G)

# Subset for quick testing (every 10th SNP)
# geno <- geno[, seq(1, ncol(geno), by = 10), drop = FALSE]

# --------------------------
# Subset by chromosome if specified
# --------------------------
if (!is.null(opt$chrom) && nzchar(opt$chrom)) {
    message("Subsetting genotype data for chromosome: ", opt$chrom)
    chrom <- opt$chrom
    chrom <- sub("^chr", "", chrom, ignore.case = TRUE)

    chr_vec <- sub("^chr([0-9]+).*", "\\1", colnames(geno))  # abstract chr from SNP IDs
    keep <- which(chr_vec == chrom)

    if (length(keep) == 0L) stop("No SNPs found for --chrom=", opt$chrom)
    geno <- geno[, keep, drop = FALSE]
}

# normalize correct from optparse (character) to R NULL
if (is.null(opt$correct) || !nzchar(opt$correct) || toupper(opt$correct) == "NULL") {
  opt$correct <- NULL
}
# normalize cluster from optparse to R NULL
if (is.null(opt$cluster) || !nzchar(opt$cluster) || toupper(opt$cluster) == "NULL") {
  opt$cluster <- NULL
}

if (is.null(opt$cluster)) {
  cat("No cluster provided; running palm.get.summary without cluster.\n")
  res <- palm.get.summary(
    null.obj = modglmm,
    covariate.interest = geno,
    correct = opt$correct
  )
}else{
  cat("Cluster provided; running palm.get.summary with cluster.\n")
  res <- palm.get.summary(
    null.obj = modglmm,
    covariate.interest = geno,
    correct = opt$correct,
    cluster = cluster
  )
}

# ---------- split by pheno and write {pheno}_step2_palm.txt ----------
# res <- as.data.frame(res, check.names = FALSE)
# res is list returned by palm.get.summary()

stopifnot(length(res) >= 1)
study_names <- names(res)
if (is.null(study_names) || any(study_names == "")) study_names <- paste0("Study", seq_along(res))
names(res) <- study_names

res_df_list <- lapply(study_names, function(d) {
  est_df <- as.data.frame(res[[d]]$est, check.names = FALSE)
  se_df  <- as.data.frame(res[[d]]$stderr, check.names = FALSE)

  colnames(est_df) <- paste0(d, ".est.", colnames(est_df))
  colnames(se_df)  <- paste0(d, ".stderr.", colnames(se_df))

  n_df <- data.frame(tmp = res[[d]]$n)
  colnames(n_df) <- paste0(d, ".n")

  cbind(est_df, se_df, n_df)
})

res <- res_df_list[[1]]
if (length(res_df_list) > 1) {
  for (i in 2:length(res_df_list)) {
    # feature rows align (usually rownames are feature IDs); if not aligned, merge by rownames
    res <- cbind(res, res_df_list[[i]])
  }
}

# inherit rownames（feature IDs）
rownames(res) <- rownames(res_df_list[[1]])


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
    ## compute p-value (Wald Z test)
    out$pval <- 2 * pnorm(-abs(out$est / out$stderr))

    # suffix: add _chr{chrom} only if --chrom is specified
    chr_suffix <- ""
    if (!is.null(opt$chrom) && nzchar(opt$chrom)) {
    chr_clean  <- sub("^chr", "", opt$chrom, ignore.case = TRUE)
    chr_suffix <- paste0("_chr", chr_clean)
    }

    out_file <- file.path(
    out_dir,
    paste0("palm1_step2", chr_suffix, "_", pheno, ".txt")
    )

    write.table(out, file = out_file, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}

cat("Wrote per-pheno files to: ", out_dir, "\n")



