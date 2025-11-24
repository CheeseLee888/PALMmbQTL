#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = FALSE)

## load R libraries
library(SAIGEQTL)
require(optparse)

print(sessionInfo())

## command line options
option_list <- list(
  make_option(
    "--abdFile",
    type    = "character", default = "",
    help    = "Required. Input abundance file (samples in columns, features in rows)."
  ),
  make_option(
    "--covFile",
    type    = "character", default = "",
    help    = "Required. Input covariate file (samples in rows, covariates in columns)."
  ),
  make_option(
    "--sampleIDColinabdFile",
    type    = "character", default = "",
    help    = "Optional. Column name for the sample IDs in the abundance file (e.g. \"IID\"). If empty, the first column is used as ID."
  ),
  make_option(
    "--sampleIDColincovFile",
    type    = "character", default = "",
    help    = "Optional. Column name for the sample IDs in the covariate file (e.g. \"IID\"). If empty, the first column is used as ID."
  ),
  make_option(
    "--mergeOutFile",
    type    = "character", default = "",
    help    = "Required. Output merged file name."
  )
)

## parse options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

## sanity checks for required files / arguments -------------------------------

if (opt$abdFile == "" || !file.exists(opt$abdFile)) {
  stop("ERROR: --abdFile must be provided and must exist.")
}

if (opt$covFile == "" || !file.exists(opt$covFile)) {
  stop("ERROR: --covFile must be provided and must exist.")
}

if (opt$mergeOutFile == "") {
  stop("ERROR: --mergeOutFile must be provided.")
}

## handle ID column arguments: "" -> NULL for our read_table_with_id ----------

abd_id_col <- if (nzchar(opt$sampleIDColinabdFile)) opt$sampleIDColinabdFile else NULL
cov_id_col <- if (nzchar(opt$sampleIDColincovFile)) opt$sampleIDColincovFile else NULL

## read tables with ID handling -----------------------------------------------

# read_table_with_id(path, id_col = NULL) is assumed to:
# - if id_col is NULL: use the first column as ID and rename it to IID
# - if id_col is not NULL: require that column, rename it to IID
abd <- read_table_with_id(opt$abdFile, id_col = abd_id_col)
cov <- read_table_with_id(opt$covFile, id_col = cov_id_col)

## merge abundance and covariate files ----------------------------------------

# merge_abd_cov(abd, cov) is assumed to:
# - require both tables to contain 'IID'
# - merge by 'IID' and keep IID as the first column
merged <- merge_abd_cov(abd, cov)

## save merged file -----------------------------------------------------------

write.table(
  merged,
  file      = opt$mergeOutFile,
  sep       = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote     = FALSE
)
cat(sprintf("Merged file saved to '%s'.\n", opt$mergeOutFile))