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
  make_option("--isRemoveZerosinPheno",
    type = "logical", default = FALSE,
    help = "Optional. Whether to remove zeros in the phenotype"
  ),
  make_option("--traitType", type = "character", default = "count", help = "Required. binary or quantitative [default=count]"),
  make_option("--invNormalize",
    type = "logical", default = FALSE,
    help = "Optional. Only for quantitative. Whether to perform the inverse normalization for the phenotype [default='FALSE']"
  ),
  make_option("--covarColList",
    type = "character", default = "",
    help = "List of covariates (comma separated)"
  ),
  make_option("--qCovarColList",
    type = "character", default = "",
    help = "List of categorical covariates (comma separated). All categorical covariates must also be in covarColList"
  ),
  make_option("--offsetCol",
    type = "character", default = "",
    help = "offset column"
  ),
  make_option("--varWeightsCol",
    type = "character", default = NULL,
    help = "variance weight column"
  ),
  make_option("--longlCol",
    type = "character", default = "",
    help = ""
  ),
  make_option("--sampleIDColinabdFile",
    type = "character", default = "IID",
    help = "Column name of sample IDs in the abundance file, e.g. IID"
  ),
  make_option("--sampleIDColincovFile",
    type = "character", default = "IID",
    help = "Column name of sample IDs in the covariate file, e.g. IID"
  ),
  make_option("--cellIDColinphenoFile",
    type = "character", default = "",
    help = "Column name of cell IDs in the phenotype file, e.g. barcode"
  ),
  make_option("--tol",
    type = "numeric", default = 0.02,
    help = "Optional. Tolerance for fitting the null GLMM to converge [default=0.02]."
  ),
  make_option("--maxiter",
    type = "integer", default = 20,
    help = "Optional. Maximum number of iterations used to fit the null GLMM [default=20]."
  ),
  make_option("--nThreads",
    type = "integer", default = 1,
    help = "Optional. Number of threads (CPUs) to use [default=1]."
  ),
  make_option("--memoryChunk",
    type = "numeric", default = 2,
    help = "Optional. Size (Gb) for each memory chunk [default=2]"
  ),
  make_option("--LOCO",
    type = "logical", default = FALSE,
    help = "Whether to apply the leave-one-chromosome-out (LOCO) approach when fitting the null model using the full GRM [default=FALSE]."
  ),
  make_option("--outputPrefix",
    type = "character", default = "~/",
    help = "Required. Path and prefix of the output files [default='~/']"
  ),
  make_option("--sparseGRMFile",
    type = "character", default = "",
    help = "Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"
  ),
  make_option("--sparseGRMSampleIDFile",
    type = "character", default = "",
    help = "Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM [default=NULL]"
  ),
  make_option("--isCovariateTransform",
    type = "logical", default = FALSE,
    help = "Optional. Whether use qr transformation on covariates [default='FALSE']."
  ),
  make_option("--useSparseGRMtoFitNULL",
    type = "logical", default = FALSE,
    help = "Optional. Whether to use sparse GRM to fit the null model [default='FALSE']."
  ),
  make_option("--sexCol",
    type = "character", default = "",
    help = "Optional. Column name for sex in the phenotype file, e.g Sex"
  ),
  make_option("--FemaleCode",
    type = "character", default = "1",
    help = "Optional. Values in the column for sex in the phenotype file are used for females [default, '1']"
  ),
  make_option("--FemaleOnly",
    type = "logical", default = FALSE,
    help = "Optional. Whether to run Step 1 for females only [default=FALSE]. if TRUE, --sexCol and --FemaleCode need to be specified"
  ),
  make_option("--MaleCode",
    type = "character", default = "0",
    help = "Optional. Values in the column for sex in the phenotype file are used for males [default, '0']"
  ),
  make_option("--MaleOnly",
    type = "logical", default = FALSE,
    help = "Optional. Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol and --MaleCode need to be specified"
  ),
  make_option("--SampleIDIncludeFile",
    type = "character", default = "",
    help = "Path to the file that contains one column for IDs of samples who will be include for null model fitting."
  ),
  make_option("--isCovariateOffset",
    type = "logical", default = TRUE,
    help = "Optional. Whether to estimate fixed effect coeffciets. [default, 'TRUE']"
  ),
  make_option("--useGRMtoFitNULL", type = "logical", default = TRUE, help = ""),
  make_option("--isShrinkModelOutput",
    type = "logical", default = FALSE,
    help = "Optional. Whether to remove unnecessary objects for step2 from the model output. [default, 'FALSE']"
  )
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

qcovars <- strsplit(opt$qCovarColList, ",")[[1]]
convertoNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) {
    return(NULL)
  })
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  } else {
    cat(stringOutput, " is ", y, "\n")
  }
  return(y)
}


# BLASctl_installed <- require(RhpcBLASctl)
# if (BLASctl_installed) {
#   # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
#   original_num_threads <- blas_get_num_procs()
#   blas_set_num_threads(1)
# }


# set seed
set.seed(1)
fitNULLGLMM_multiV(
  grmFile = opt$grmFile,
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
  traitType = opt$traitType,
  invNormalize = opt$invNormalize,
  covarColList = covars,
  qCovarCol = qcovars,
  offsetCol = opt$offsetCol,
  varWeightsCol = opt$varWeightsCol,
  longlCol = opt$longlCol,
  sampleIDColinabdFile = opt$sampleIDColinabdFile,
  sampleIDColincovFile = opt$sampleIDColincovFile,
  cellIDColinphenoFile = opt$cellIDColinphenoFile,
  tol = opt$tol,
  maxiter = opt$maxiter,
  nThreads = opt$nThreads,
  memoryChunk = opt$memoryChunk,
  LOCO = opt$LOCO,
  outputPrefix = opt$outputPrefix,
  sparseGRMFile = opt$sparseGRMFile,
  sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
  isCovariateTransform = opt$isCovariateTransform,
  useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
  sexCol = opt$sexCol,
  FemaleCode = opt$FemaleCode,
  FemaleOnly = opt$FemaleOnly,
  MaleCode = opt$MaleCode,
  MaleOnly = opt$MaleOnly,
  SampleIDIncludeFile = opt$SampleIDIncludeFile,
  isCovariateOffset = opt$isCovariateOffset,
  useGRMtoFitNULL = opt$useGRMtoFitNULL,
  isShrinkModelOutput = opt$isShrinkModelOutput
)



# if (BLASctl_installed) {
#   # Restore originally configured BLAS thread count
#   blas_set_num_threads(original_num_threads)
# }
