#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--abdFile",
    type = "character", default = "",
    help = "Path to the microbiome abundance data file"
  ),
  make_option("--outputFile",
    type = "character", default = "",
    help = "Path to the output file for feature information"
  ),
  make_option("--genoFile",
    type = "character", default = "",
    help = "Path to the genotype file (without extension) for family information"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate abdFile before reading
if (is.null(opt$abdFile) || !nzchar(opt$abdFile)) {
  stop("The --abdFile argument is missing or empty. Please provide a valid file path.")
}
if (!file.exists(opt$abdFile)) {
  stop("The specified abdFile does not exist: ", opt$abdFile)
}

# Summarize microbiome abundance data
abd_data <- read.table(opt$abdFile, header = TRUE, sep = "\t")
feature_info <- data.frame(
  FeatureID = colnames(abd_data)[-1],
  Prevalence = colSums(!is.na(abd_data[, -1])) / nrow(abd_data),
  AvgProportion = colMeans(abd_data[, -1], na.rm = TRUE)
)

# Write feature information to the output file
if (is.null(opt$outputFile) || !nzchar(opt$outputFile)) {
  stop("The --outputFile argument is missing or empty. Please provide a valid file path.")
}
write.table(feature_info, file = opt$outputFile, sep = "\t", quote = FALSE, row.names = FALSE)

# Read family information from the PLINK .fam file derived from --genoFile
fam_file <- paste0(opt$genoFile, ".fam")
if (!file.exists(fam_file)) {
  stop("The corresponding .fam file does not exist: ", fam_file)
}

fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

# Total samples and families
total_samples <- nrow(fam_data)
total_families <- length(unique(fam_data$FID))

# Append this information to the feature info file
cat("Total samples: ", total_samples, "\n", file = opt$outputFile, append = TRUE)
cat("Total families: ", total_families, "\n", file = opt$outputFile, append = TRUE)

cat("Feature information written to: ", opt$outputFile, "\n")
cat("Total samples and families information appended to: ", opt$outputFile, "\n")