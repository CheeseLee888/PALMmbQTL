#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param grmFile character. Path to the GRM file in RDS format. The RDS file contains a list with two elements: sample.id and K. sample.id is a vector of sample IDs. K is the GRM matrix with row and column names as sample IDs. By default, "".
#' @param abdFile character. Path to the abundance file.
#' @param covFile character. Path to the covariate file.
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age").
#' @param offsetCol character. Offset column name in the merged data. If not specified, SeqDepth will be calculated from the abundance file and used as the offset.
#' @param outputPrefix character. Path to the output files with prefix.
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, TRUE.
#' @param useGRMtoFitNULL logical. Whether to use the GRM to fit the NULL model. If FALSE, an identity matrix will be used. By default, TRUE.
#' @param batch_idx integer. Index of the batch to process (starting from 1). By default, 1. If you already have batch files from 1 to N-1, you can set batch_idx to N to continue processing the remaining phenotypes.
#' @param batch_size integer. Number of phenotypes to save in each batch file. By default, 5.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV <- function(grmFile = "",
                               abdFile = "",
                               covFile = "",
                               covarColList = character(0),
                               offsetCol = "",
                               outputPrefix = "",
                               isCovariateOffset = TRUE,
                               useGRMtoFitNULL = TRUE,
                               batch_idx = 1,
                               batch_size = 5) {
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
  pheno_names <- setdiff(colnames(abd), "IID")
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
  # cat(colnames(data), "\n")

  ## ------------------------------------------------------------------------
  ## Force the use of SeqDepth as offset：
  ## - If the user does not specify offsetCol: Generate SeqDepth by summing rows in abd
  ## - If the user specifies offsetCol: Check if the column exists and is valid
  ## ------------------------------------------------------------------------

  if (offsetCol == "") {
    cat("No offset column is specified. Calculate SeqDepth from abundance file...\n")

    seqdepth <- rowSums(pheno_list, na.rm = TRUE)
    names(seqdepth) <- abd[[ "IID" ]]
    
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

    if (!all(grm_id == pheno_id)) {
      stop("GRM ID and phenotype ID are not in the same order!\n")
    }

  } else {
    cat("Identity matrix will be used to fit the NULL model\n")
    grm_K <- diag(nrow(data))
  }
  rownames(grm_K) <- colnames(grm_K) <- pheno_id

  # Core step1 for PALM-mbQTL
  null_list <- list()
  cat("Total phenotypes to be processed: ", length(pheno_names), "\n")

  ## --- add: batch saving settings ---
  # batch_idx is user-controlled, must be >= 1
  if (batch_idx < 1L) {
    stop("batch_idx must be >= 1")
  }

  start_i <- (batch_idx - 1L) * batch_size + 1L
  end_i   <- length(pheno_names)

  if (start_i > end_i) {
    stop("start_i exceeds number of phenotypes: start_i = ",
        start_i, ", total = ", end_i)
  }

  cat("Start processing phenotypes from index ", start_i,
      " to ", end_i, "\n", sep = "")

  batch_dir  <- paste0(outputPrefix, "_batches")
  failed_pheno_file <- paste0(outputPrefix, "_failed_pheno.txt")
  if (batch_idx == 1 && dir.exists(batch_dir)) {
    cat("batch_idx == 1: cleaning existing batch directory:", batch_dir, "\n")
    unlink(list.files(batch_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
  }
  if (batch_idx == 1 && file.exists(failed_pheno_file)) {
    cat("batch_idx == 1: removing failed pheno file:", failed_pheno_file, "\n")
    unlink(failed_pheno_file, force = TRUE)
  }
  if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)
  ## ---------------------------------

  for (i in seq.int(start_i, end_i)) {
    pheno_name <- pheno_names[i]
    cat("[", i, "/", length(pheno_names), "] Fitting NULL GLMM for pheno: ", pheno_name, "\n")

    # construct the formula
    if (length(covarColList) > 0) {
      formula <- paste0(pheno_name, "~", paste0(covarColList, collapse = "+"))
    } else {
      formula <- paste0(pheno_name, "~ 1")
    }

    if (isCovariateOffset && offsetCol != "") {
      formula <- paste0(formula, "+offset(log(", offsetCol, "))")
    }
    if (i == start_i || i == end_i) {
      cat("formula is ", formula, "\n")
    }

    # run glmmkin (skip failed phenotypes)
    if (useGRMtoFitNULL) {
      modglmm <- tryCatch(
        {
          GMMAT::glmmkin(
            formula,
            data   = data,
            kins   = grm_K,
            id     = "IID",
            family = poisson(link = "log")
          )
        },
        error = function(e) {
          cat("pheno: ", pheno_name, ", glmmkin failed.\n")
          write(
            pheno_name,
            file   = failed_pheno_file,
            append = TRUE
          )
          NULL
        }
      )
    }else{
      cat("No kins matrix provided to glmmkin\n")
      modglmm <- tryCatch(
        {
          GMMAT::glmmkin(
            formula,
            data   = data,
            id     = "IID",
            family = poisson(link = "log")
          )
        },
        error = function(e) {
          cat("pheno: ", pheno_name, ", glmmkin failed.\n")
          write(
            pheno_name,
            file   = failed_pheno_file,
            append = TRUE
          )
          NULL
        }
      )
    }

    if (!is.null(modglmm)) {
      cat("pheno: ", pheno_name, ", glmmkin succeed.\n")
      null_list[[pheno_name]] <- list(
        pheno_name = pheno_name,
        modglmm    = modglmm
      )
    }

    ## --- add: save every <batch_size> phenotypes OR at the end, then clear null_list ---
    if (i %% batch_size == 0L || i == length(pheno_names)) {
      batch_file <- file.path(batch_dir, sprintf("batch_%04d.rds", batch_idx))
      saveRDS(null_list, file = batch_file)
      cat("Saved batch ", batch_idx, " to ", batch_file,
          " (n=", length(null_list), ")\n", sep = "")
      # cat("Contains pheno indexes: ",
      #     paste0(seq.int(i - length(null_list) + 1, i), collapse = ", "), "\n")

      # clear in-memory list to reduce memory usage
      null_list <- list()
      invisible(gc())
      batch_idx <- batch_idx + 1
    }
    ## ------------------------------------------------------------------------
  }


  ## --- add: merge all batch files into one null_list, then save to modelOut ---
  batch_files <- list.files(batch_dir, pattern = "^batch_[0-9]{4}\\.rds$", full.names = TRUE)
  batch_files <- sort(batch_files)
  if (length(batch_files) == 0) stop("No batch files found in ", batch_dir)

  null_list <- list()
  for (bf in batch_files) {
    tmp <- readRDS(bf)
    if (length(tmp) > 0) {
      # warn if duplicates (shouldn't happen normally)
      dup <- intersect(names(null_list), names(tmp))
      if (length(dup) > 0) warning("Duplicate phenotypes when merging: ", paste(dup, collapse = ", "))
      null_list <- c(null_list, tmp)
    }
  }
  cat("Merged all batch files.\n")

  save(null_list, file = modelOut)
  cat("NULL model has been saved to ", modelOut,
      " (total successful phenotypes: ", length(null_list), ")\n", sep = "")

  ## --- clean up batch files ---
  unlink(batch_dir, recursive = TRUE, force = TRUE)
  cat("Temporary batch files have been removed from ", batch_dir, "\n", sep = "")
}
