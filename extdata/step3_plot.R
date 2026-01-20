#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(qqman)
})

# ============================
# Parse command-line options
# ============================
option_list <- list(
  make_option(
    "--indir",
    type = "character",
    default = "",
    help = "Input directory containing result files"
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "",
    help = "Output directory for plots and top hits"
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

# ----------------------------
# Validate options
# ----------------------------
if (opt$indir == "" || !dir.exists(opt$indir)) {
  stop("Invalid --indir: directory does not exist or is empty.")
}
if (opt$outdir == "") {
  stop("Invalid --outdir: must be provided.")
}

indir  <- opt$indir
outdir <- opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================
# Discover input files
# ============================
files <- list.files(indir, full.names = TRUE)
if (length(files) == 0) {
  stop("No files found in --indir.")
}

cat("Total files found:", length(files), "\n")

# ============================
# Process each file
# ============================
required_cols <- c("SNP", "CHR", "POS", "pval")

for (infile in files) {

  prefix <- sub("\\.[^.]+$", "", basename(infile))
  cat("\n[Processing]", basename(infile), "\n")

  # Try reading file
  dt <- tryCatch(
    fread(infile, sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE),
    error = function(e) {
      message("  ! Skip (read error): ", e$message)
      return(NULL)
    }
  )
  if (is.null(dt)) next

  # Check required columns
  if (!all(required_cols %in% colnames(dt))) {
    message("  ! Skip (missing required columns)")
    next
  }

  man <- data.frame(
    SNP = as.character(dt$SNP),
    CHR = dt$CHR,
    BP  = dt$POS,
    P   = dt$pval,
    stringsAsFactors = FALSE
  )

  if (nrow(man) == 0) {
    message("  ! Skip (no valid variants after cleaning)")
    next
  }

  cat("  Variants used:", nrow(man), "\n")

  # ----------------------------
  # Manhattan plot
  # ----------------------------
  jpeg(
    file.path(outdir, paste0("manhattan_", prefix, ".jpg")),
    width = 11, height = 5.5, units = "in", res = 300
  )
  qqman::manhattan(
    man,
    chr = "CHR",
    bp  = "BP",
    p   = "P",
    snp = "SNP",
    logp = TRUE,
    suggestiveline = -log10(1e-5),
    genomewideline = -log10(5e-8),
    cex = 0.6,
    main = "Manhattan Plot"
  )
  dev.off()

  # ----------------------------
  # QQ plot
  # ----------------------------
  jpeg(
    file.path(outdir, paste0("qq_", prefix, ".jpg")),
    width = 6, height = 6, units = "in", res = 300
  )
  qqman::qq(man$P, main = "QQ Plot")
  dev.off()

  # ----------------------------
  # Top hits
  # ----------------------------
  top <- man[order(man$P), ][seq_len(min(10, nrow(man))), ]
  fwrite(
    top,
    file.path(outdir, paste0("top10_", prefix, ".txt")),
    sep = "\t"
  )

  message("  Done: ", prefix)
}

message("\nAll done. Results written to: ", outdir)
