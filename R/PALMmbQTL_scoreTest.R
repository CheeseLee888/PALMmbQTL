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

  # checkArgsListBool(
  # )
  # checkArgsListNumeric(
  # )


  load(NULLmodelFile)

  ## extract all pheno_name
  all_pheno_names <- vapply(
    null_list,
    function(x) x$pheno_name,
    FUN.VALUE = character(1L)
  )

  ## If phenoCol is empty, run all phenotypes
  if (identical(phenoCol, "") || is.null(phenoCol)) {
    cat(sprintf("No phenotype specified — running all %d phenotypes.\n", length(all_pheno_names)))
    idxs <- seq_along(all_pheno_names)
  } else {
    ## use which to find the corresponding index
    idx <- which(all_pheno_names == phenoCol)

    if (length(idx) == 0) {
      stop(sprintf("Phenotype '%s' not found in NULLmodelFile.", phenoCol))
    } else {
      cat(sprintf("Phenotype '%s' found in NULLmodelFile at index %d.\n", phenoCol, idx))
    }
    idxs <- idx
  }

  for (i in idxs) {
    pheno_name <- all_pheno_names[i]
    modglmm <- null_list[[i]]$modglmm

    # construct outfile name: use provided prefix if given, else default prefix
    if (identical(PALMOutputFile, "") || is.null(PALMOutputFile)) {
      # out_file <- paste0("PALM_output_", pheno_name, ".txt")
      stop("Please provide PALMOutputFile prefix.")
    } else {
      out_file <- paste0(PALMOutputFile, "_", pheno_name, ".txt")
    }

    cat(sprintf("Starting glmmscore test for phenotype '%s' (index %d)...\n", pheno_name, i))
    GMMAT::glmm.score(
      obj = modglmm,
      infile = inFile,
      center = T,
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
    res$stat  <- res$est / res$stderr
    res$pval   <- 2 * pnorm(-abs(res$stat))

    # keep only required columns
    res_out <- res[, c("SNP", "est", "stderr", "stat", "pval")]

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
