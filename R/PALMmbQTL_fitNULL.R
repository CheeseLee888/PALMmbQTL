#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param grmFile character. Path to the GRM file in RDS format. The RDS file contains a list with two elements: sample.id and K. sample.id is a vector of sample IDs. K is the GRM matrix with row and column names as sample IDs. By default, "".
#' @param abdFile character. Path to the abundance file.
#' @param covFile character. Path to the covariate file.
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age").
#' @param offsetCol character. Offset column name in the merged data. If not specified, SeqDepth will be calculated from the abundance file and used as the offset.
#' @param sampleIDColinabdFile character. Column name for the sample IDs in the abundance file e.g. "IID".
#' @param sampleIDColincovFile character. Column name for the sample IDs in the covariate file e.g. "IID".
#' @param outputPrefix character. Path to the output files with prefix.
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, TRUE.
#' @param useGRMtoFitNULL logical. Whether to use the GRM to fit the NULL model. If FALSE, an identity matrix will be used. By default, TRUE.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV <- function(grmFile = "",
                               abdFile = "",
                               covFile = "",
                               covarColList = NULL,
                               offsetCol = "",
                               sampleIDColinabdFile = "IID",
                               sampleIDColincovFile = "IID",
                               outputPrefix = "",
                               isCovariateOffset = TRUE,
                               useGRMtoFitNULL = TRUE) {
  ## set up output files
  modelOut <- paste0(outputPrefix, ".rda")
  # file.create(modelOut, showWarnings = TRUE)

  ## sanity checks for required files / arguments -------------------------------
  if (abdFile == "" || !file.exists(abdFile)) {
    stop("ERROR: abdFile must be provided and must exist.")
  }

  if (covFile == "" || !file.exists(covFile)) {
    stop("ERROR: covFile must be provided and must exist.")
  }

  ## read tables with ID handling -----------------------------------------------
  abd <- read_table_with_id(abdFile)
  cov <- read_table_with_id(covFile)
  cat("Abundance and covariate files have been read\n")

  ## read phenotype list ------------------------------------------------------
  pheno_names <- setdiff(colnames(abd), sampleIDColinabdFile)
  pheno_list <- abd[, pheno_names, drop = FALSE]
  # Simple sanity check
  if (any(pheno_list < 0, na.rm = TRUE)) {
    warning("Abundance table contains negative values; 
            please make sure abdFile is raw counts when using SeqDepth.")
  }

  ## merge abundance and covariate files ----------------------------------------

  # merge_abd_cov(abd, cov) is assumed to:
  # - require both tables to contain 'IID'
  # - merge by 'IID' and keep IID as the first column
  merged <- merge_abd_cov(abd, cov)
  data <- merged
  cat("Abundance and covariate files have been merged\n")
  cat(colnames(data), "\n")

  ## ------------------------------------------------------------------------
  ## Force the use of SeqDepth as offset：
  ## - If the user does not specify offsetCol: Generate SeqDepth by summing rows in abd
  ## - If the user specifies offsetCol: Check if the column exists and is valid
  ## ------------------------------------------------------------------------

  if (offsetCol == "") {
    cat("No offset column is specified. Calculate SeqDepth from abundance file...\n")

    seqdepth <- rowSums(pheno_list, na.rm = TRUE)
    names(seqdepth) <- abd[[sampleIDColinabdFile]]
    
    data$SeqDepth <- seqdepth[match(data$IID, names(seqdepth))]

    # Check for NA and non-positive values
    if (any(is.na(data$SeqDepth))) {
      warning("Some samples in merged data do not have SeqDepth (no matching abd).")
    }
    if (any(data$SeqDepth <= 0, na.rm = TRUE)) {
      stop("SeqDepth contains non-positive values. Please check abdFile.")
    }

    # Set offsetCol to SeqDepth for use in the formula later
    offsetCol <- "SeqDepth"
    
  } else {
    cat("Offset column", offsetCol, "is specified.\n")

    # Ensure the column exists
    if (!offsetCol %in% colnames(data)) {
      stop("Specified offsetCol '", offsetCol, "' not found in merged data.")
    }

    # Also check if values are > 0
    if (any(data[[offsetCol]] <= 0, na.rm = TRUE)) {
      stop("Offset column '", offsetCol, "' contains non-positive values.")
    }
  }


  if (length(covarColList) > 0) {
    cat("Covariates: ", covarColList, "\n")
    if (!all(covarColList %in% colnames(data))) {
      stop("ERROR! all covariates in covarColList must be in the merged data\n")
    }
  }


  cat("Start fitting the NULL GLMM\n")

  pheno_id <- data[["IID"]]
  # GRM
  if (useGRMtoFitNULL) {
    cat("GRM will be used to fit the NULL model\n")
    if (!file.exists(grmFile)) {
      stop("ERROR! grmFile ", grmFile, " does not exist\n")
    }
    grm_obj <- readRDS(grmFile)
    grm_K   <- grm_obj$K
    grm_id  <- grm_obj$sample.id

    cat("grm id:\n")
    cat(grm_id[1:5], "\n")
    cat("pheno id:\n")
    cat(pheno_id[1:5], "\n")

    if (!setequal(grm_id, pheno_id)) {
      stop("ERROR! the sample IDs in the GRM file are not the same as those in the phenotype file\n")
    }else {
      cat("All sample IDs in the GRM file are the same as those in the phenotype file\n")
    }
    if (!all(grm_id == pheno_id)) {
      cat("GRM ID and phenotype ID are not in the same order; reorder GRM ...\n")
      idx      <- match(pheno_id, grm_id)
      grm_K    <- grm_K[idx, idx, drop = FALSE]
    }

  } else {
    cat("Identity matrix will be used to fit the NULL model\n")
    grm_K <- diag(nrow(data))
  }
  rownames(grm_K) <- colnames(grm_K) <- pheno_id

  # Core step1 for PALM-mbQTL
  null_list <- list()
  failed_pheno <- character()
  cat("Total phenotypes to be processed: ", length(pheno_names), "\n")

  for (i in seq_along(pheno_names)) {
    pheno_name <- pheno_names[i]
    cat("[", i, "/", length(pheno_names), "] Fitting NULL GLMM for pheno: ", pheno_name, "\n")

    # construct the formula
    if (length(covarColList) > 0) {
      formula <- paste0(pheno_name, "~", paste0(covarColList,
        collapse = "+"
      ))
    } else {
      formula <- paste0(pheno_name, "~ 1")
    }
    
    if(isCovariateOffset & offsetCol != ""){
      formula <- paste0(formula, "+offset(log(", offsetCol, "))")
    }
    cat("formula is ", formula, "\n")

    # run glmmkin (skip failed phenotypes)
    modglmm <- tryCatch(
      {
        callr::r(
          func = function(formula, data, grm_K) {
            GMMAT::glmmkin(
              formula,
              data   = data,
              kins   = grm_K,
              id     = "IID",
              family = poisson(link = "log")
            )
          },
          args = list(formula = as.formula(formula), data = data, grm_K = grm_K),
          show = FALSE
        )
      },
      error = function(e) {
        failed_pheno <<- c(failed_pheno, pheno_name)
        NULL
      }
    )

    # if failed, skip this phenotype
    if (is.null(modglmm)) {
      next
    }

    cat("pheno: ", pheno_name, ", glmmkin succeed.\n")

    null_list[[pheno_name]] <- list(
      pheno_name  = pheno_name,
      modglmm     = modglmm
    )
  }

  if (length(failed_pheno) > 0) {
    cat("Total failed phenotypes:", length(failed_pheno), "\n")
    writeLines(
      failed_pheno,
      con = paste0(outPrefix, "_failed_pheno.txt")
    )
    cat("List of failed phenotypes has been saved to ",
      paste0(outputPrefix, "_failed_pheno.txt"), "\n")
  }

  save(null_list, file = modelOut)
  cat("NULL model has been saved to ", modelOut, "\n")

}
