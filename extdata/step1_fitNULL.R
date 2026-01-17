#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = F)

## load R libraries

library(SAIGEQTL)

require(optparse) # install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--grmFile",
    type = "character", default = "",
    help = "Path to save the GRM RDS file"
  ),
  make_option("--abdFile",
    type = "character", default = "",
    help = "Required. Path to the abundance file."
  ),
  make_option("--covFile",
    type = "character", default = "",
    help = "Required. Path to the covariate file."
  ), 
  make_option("--covarColList",
    type = "character", default = "",
    help = "List of covariates (comma separated)"
  ),
  make_option("--offsetCol",
    type = "character", default = "",
    help = "offset column"
  ),
  make_option("--sampleIDColinabdFile",
    type = "character", default = "IID",
    help = "Column name of sample IDs in the abundance file, e.g. IID"
  ),
  make_option("--sampleIDColincovFile",
    type = "character", default = "IID",
    help = "Column name of sample IDs in the covariate file, e.g. IID"
  ),
  make_option("--outputPrefix",
    type = "character", default = "~/",
    help = "Required. Path and prefix of the output files [default='~/']"
  ),
  make_option("--isCovariateOffset",
    type = "logical", default = TRUE,
    help = "Optional. Whether to estimate fixed effect coeffciets. [default, 'TRUE']"
  ),
  make_option("--useGRMtoFitNULL", type = "logical", default = TRUE, help = "")
)


## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

## covariates: if covarColList is ALL, infer from covFile header
if (opt$covarColList!="all") {
  cat("Using user-specified covariates from --covarColList\n")
  covars <- strsplit(opt$covarColList, ",")[[1]]
} else {
  cat("Using all covariates from covFile header\n")
  cov <- read_table_with_id(opt$covFile)
  cov_header <- colnames(cov)

  ## Remove ID and offset columns
  drop_cols <- c(opt$sampleIDColincovFile, opt$offsetCol)

  covars <- setdiff(cov_header, drop_cols)

  if (!length(covars)) {
    stop("No covariate columns found in covFile after removing ID and offset columns. Please set '--covarColList= '(empty) and try again.\n")
  }
  print(covars)
}

# convertoNumeric <- function(x, stringOutput) {
#   y <- tryCatch(expr = as.numeric(x), warning = function(w) {
#     return(NULL)
#   })
#   if (is.null(y)) {
#     stop(stringOutput, " is not numeric\n")
#   } else {
#     cat(stringOutput, " is ", y, "\n")
#   }
#   return(y)
# }

set.seed(1)
fitNULLGLMM_multiV(
  grmFile = opt$grmFile,
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  covarColList = covars,
  offsetCol = opt$offsetCol,
  sampleIDColinabdFile = opt$sampleIDColinabdFile,
  sampleIDColincovFile = opt$sampleIDColincovFile,
  outputPrefix = opt$outputPrefix,
  isCovariateOffset = opt$isCovariateOffset,
  useGRMtoFitNULL = opt$useGRMtoFitNULL
)
