#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
#'
#' @param inFile character. Path to the input file containing the genotype data in PLINK binary format (bed/bim/fam) or VCF format.
#' @param NULLmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param PALMOutputFile character. Prefix of the output files containing assoc test results
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0.05.
#' @return PALMOutputFile
#' @export
SPAGMMATtest <- function(inFile = "",
                         phenoCol = "",
                         NULLmodelFile = "",
                         PALMOutputFile = "",
                         minMAF = 0.05) {

  load(NULLmodelFile)
  if (!exists("null_list")) stop("NULLmodelFile does not contain object 'null_list'.")
  if (length(null_list) == 0) stop("null_list is empty in NULLmodelFile (no phenotype models available).")

  ## If phenoCol is empty, run all phenotypes
  if (is.null(phenoCol) || toupper(phenoCol) %in% c("", "NULL", "ALL")) {
    pheno_keys <- names(null_list)
    if (is.null(pheno_keys) || length(pheno_keys) == 0) stop("null_list has no names (pheno keys).")
    cat(sprintf("No phenotype specified — running all %d phenotypes.\n", length(pheno_keys)))
  } else {
    if (!phenoCol %in% names(null_list)) {
      stop(sprintf("Phenotype '%s' not found in NULLmodelFile.", phenoCol))
    }
    pheno_keys <- phenoCol
    cat(sprintf("Phenotype '%s' found in NULLmodelFile.\n", phenoCol))
  }

  for (pheno_key in pheno_keys) {
    pheno_name <- null_list[[pheno_key]]$pheno_name
    modglmm    <- null_list[[pheno_key]]$modglmm

    # construct outfile name: use provided prefix if given, else default prefix
    if (identical(PALMOutputFile, "") || is.null(PALMOutputFile)) {
      # out_file <- paste0("PALM_output_", pheno_name, ".txt")
      stop("Please provide PALMOutputFile prefix.")
    } else {
      out_file <- paste0(PALMOutputFile, "_", pheno_name, ".txt")
    }

    cat(sprintf("Starting glmmscore test for phenotype '%s' ...\n", pheno_name))
    GMMAT::glmm.score(
      obj = modglmm,
      infile = inFile,
      center = TRUE,
      outfile = out_file,
      MAF.range = c(minMAF, 0.5)
    )
    cat(sprintf("glmmscore test for '%s' finished. Output: %s\n", pheno_name, out_file))

    # --------------------------
    # Post-process: keep only SNP SCORE VAR PVAL
    # --------------------------
    res <- data.table::fread(out_file, data.table = FALSE)

    # normalize column names (only what we need)
    if (!"SNP" %in% names(res)) {
      stop("Missing SNP column in file: ", out_file)
    }
    if (!"SCORE" %in% names(res) || !"VAR" %in% names(res)) {
      stop("Missing SCORE or VAR column in file: ", out_file)
    }

    # compute
    res$est    <- res$SCORE / res$VAR
    res$stderr <- sqrt(1 / res$VAR)
    # res$stat  <- res$est / res$stderr
    # res$pval   <- 2 * pnorm(-abs(res$stat))
    res$pval <- 1 - pchisq((res$est / res$stderr)^2, df = 1)

    # keep only required columns
    # res_out <- res[, c("SNP", "est", "stderr", "pval")]
    res_out <- res[, c("SNP", "CHR", "POS", "est", "stderr", "pval")]
    # res_out <- res

    data.table::fwrite(
      res_out,
      out_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
    cat(sprintf("Post-processed output saved (SNP,SCORE,VAR,PVAL): %s\n", out_file))
  }

}
