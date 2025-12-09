#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(PALM)
})

option_list <- list(
  make_option("--abdFile", type = "character",
              help = "TSV file with sample_id + taxa columns"),
  make_option("--covFile", type = "character",
              help = "TSV covariate file with sample_id column"),
  make_option("--sampleIDCol", type = "character", default = "sample_id",
              help = "Sample ID column name shared by abdFile and covFile"),
  make_option("--covariateInterest", type = "character",
              help = "Column name in covFile used as exposure (e.g. sex)"),
  make_option("--covariateAdjust", type = "character", default = NA,
              help = "Comma-separated adjust covariates, e.g. 'age,PC1,PC2,batch'"),
  make_option("--outPrefix", type = "character", default = "output/PALM_result",
              help = "Output prefix (without extension)"),
  make_option("--correct", type = "character", default = "median",
              help = "Correction method (e.g. 'median' or 'NULL')")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$abdFile) || is.null(opt$covFile) ||
    is.null(opt$covariateInterest)) {
  stop("Must provide --abdFile, --covFile, and --covariateInterest.")
}

## 1) Read data -------------------------------------------------------------

abd <- read.table(opt$abdFile, header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.table(opt$covFile, header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)

sid <- opt$sampleIDCol
if (!sid %in% colnames(abd))
  stop("sample ID column ", sid, " not found in abdFile")
if (!sid %in% colnames(cov))
  stop("sample ID column ", sid, " not found in covFile")

## Align sample order
common <- intersect(abd[[sid]], cov[[sid]])
if (length(common) == 0) {
  stop("No overlapping sample IDs between abdFile and covFile.")
}

abd <- abd[match(common, abd[[sid]]), , drop = FALSE]
cov <- cov[match(common, cov[[sid]]), , drop = FALSE]

## 2) Construct rel.abd (matrix) -------------------------------------------

taxa_cols <- setdiff(colnames(abd), sid)
if (length(taxa_cols) == 0) {
  stop("No taxa columns found in abdFile (only sample_id?).")
}

rel.abd <- as.matrix(abd[, taxa_cols, drop = FALSE])
rownames(rel.abd) <- common

## 3) covariate.interest (matrix, force numeric) ---------------------------

if (!opt$covariateInterest %in% colnames(cov)) {
  stop("covariateInterest ", opt$covariateInterest,
       " not found in covFile.")
}
ci <- cov[[opt$covariateInterest]]

# If not numeric: e.g., sex = "M"/"F", automatically convert to 0/1
if (!is.numeric(ci)) {
  u <- sort(unique(ci[!is.na(ci)]))
  if (length(u) == 2L) {
    # Binary, convert to 0/1, second level is 1 (similar to model.matrix convention)
    ci_num <- as.numeric(ci == u[2L])
    message("Converted covariate.interest ", opt$covariateInterest,
            " from character to 0/1 numeric (", u[1L], "=0, ", u[2L], "=1).")
    ci <- ci_num
  } else {
    stop("covariate.interest ", opt$covariateInterest,
         " is non-numeric with ", length(u),
         " levels; please provide a numeric column or a binary variable.")
  }
}

cov_int <- matrix(ci, ncol = 1,
                  dimnames = list(common, opt$covariateInterest))


## 4) covariate.adjust (optional) ------------------------------------------

cov_adj <- NULL
if (!is.na(opt$covariateAdjust)) {
  adj_cols <- strsplit(opt$covariateAdjust, ",")[[1]]
  adj_cols <- trimws(adj_cols)
  missing <- setdiff(adj_cols, colnames(cov))
  if (length(missing)) {
    stop("Adjust covariates not found in covFile: ",
         paste(missing, collapse = ", "))
  }
  cov_adj <- as.data.frame(cov[, adj_cols, drop = FALSE])
  rownames(cov_adj) <- common

  # Convert character columns to factors for PALM/model.matrix to create dummy variables
  cov_adj[] <- lapply(cov_adj, function(x) {
    if (is.character(x)) {
      factor(x)
    } else {
      x
    }
  })
}


## 5) Single study directly call PALM (not wrapped as list) ----------------

cat("Running PALM::palm() (single-study mode) ...\n")
palm.res <- PALM::palm(
  rel.abd            = rel.abd,
  covariate.interest = cov_int,
  covariate.adjust   = cov_adj,
  prev.filter        = 0,
  correct            = opt$correct
)

## For single study, the return value is a data.frame,
## not a list divided by covariate names, so treat it as a table:
res.df <- palm.res

## 6) Output results -------------------------------------------------------------

out_rds <- paste0(opt$outPrefix, ".rds")
out_tsv <- paste0(opt$outPrefix, ".tsv")
dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)

saveRDS(res.df, out_rds)
write.table(res.df, file = out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("PALM result written to:\n  ", out_rds, "\n  ", out_tsv, "\n")
