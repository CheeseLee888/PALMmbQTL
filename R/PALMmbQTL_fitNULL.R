#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param grmFile character. Path to the GRM file in RDS format. The RDS file contains a list with two elements: sample.id and K. sample.id is a vector of sample IDs. K is the GRM matrix with row and column names as sample IDs. By default, "".
#' @param abdFile character. Path to the abundance file.
#' @param covFile character. Path to the covariate file.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary", "quantitative", "count" or "count_nb". By default, "count".
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the null model. All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex").
#' @param eCovarCol vector of characters. Covariates of environmental factors/cell context to be used in the null model. All covariates listed in eCovarCol must be also in covarColList,  e,g c("cellType").
#' @param sampleIDColinabdFile character. Column name for the sample IDs in the abundance file e.g. "IID".
#' @param sampleIDColincovFile character. Column name for the sample IDs in the covariate file e.g. "IID".
#' @param cellIDColinphenoFile character. Column name for the cell IDs in the phenotype file e.g. "barcode".
#' @param tol numeric. The tolerance for fitting the null model to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param nThreads integer. Number of threads to be used. By default, 1
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. By default, FALSE
#' @param outputPrefix character. Path to the output files with prefix.
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. By default, ""
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. By default, ""
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null model. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1'
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, ''
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, FALSE.
#' @param isShrinkModelOutput logical. remove unnecessary objects for step2 from the model output. By default, FALSE.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV <- function(grmFile = "",
                               abdFile = "",
                               covFile = "",
                               phenoCol = "",
                               isRemoveZerosinPheno = FALSE,
                               traitType = "count",
                               invNormalize = FALSE,
                               covarColList = NULL,
                               qCovarCol = NULL,
                               eCovarCol = NULL,
                               offsetCol = "",
                               varWeightsCol = NULL,
                               longlCol = "",
                               sampleIDColinabdFile = "IID",
                               sampleIDColincovFile = "IID",
                               cellIDColinphenoFile = "",
                               tol = 0.02,
                               maxiter = 20,
                               nThreads = 1,
                               skipModelFitting = FALSE,
                               memoryChunk = 2,
                               LOCO = FALSE,
                               outputPrefix = "",
                               sparseGRMFile = "",
                               sparseGRMSampleIDFile = "",
                               isCovariateTransform = FALSE,
                               useSparseGRMtoFitNULL = FALSE,
                               sexCol = "",
                               FemaleCode = 1,
                               FemaleOnly = FALSE,
                               MaleCode = 0,
                               MaleOnly = FALSE,
                               SampleIDIncludeFile = "",
                               isCovariateOffset = FALSE,
                               useGRMtoFitNULL = TRUE,
                               isShrinkModelOutput = FALSE) {
  ## set up output files
  modelOut <- paste0(outputPrefix, ".rda")

  if (skipModelFitting) {
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
  } else {
    file.create(modelOut, showWarnings = TRUE)
  }

  # if (!useGRMtoFitNULL) {
  #   useSparseGRMtoFitNULL <- FALSE
  #   cat("No GRM will be used to fit the NULL model and nThreads is set to 1\n")
  # }

  if (nThreads > 1) {
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads will be used ", "\n")
  }
  # set_g_omp_num_threads(nThreads)

  if (FemaleOnly & MaleOnly) {
    stop("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job\n")
  }

  if (FemaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_FemaleOnly")
    cat(
      "Female-specific model will be fitted. Samples coded as ",
      FemaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  } else if (MaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_MaleOnly")
    cat(
      "Male-specific model will be fitted. Samples coded as ",
      MaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  }


  if (longlCol == "") {
    checkColList <- c(phenoCol, covarColList, "IID")
  } else {
    checkColList <- c(phenoCol, covarColList, "IID", longlCol)
  }

  if (isCovariateOffset & length(offsetCol) != "") {
    cat("Use offset term: ", offsetCol, "\n")
    checkColList <- c(checkColList, offsetCol)
  }else{
    cat("No offset term is used\n")
  }

  ## sanity checks for required files / arguments -------------------------------
  if (abdFile == "" || !file.exists(abdFile)) {
    stop("ERROR: abdFile must be provided and must exist.")
  }

  if (covFile == "" || !file.exists(covFile)) {
    stop("ERROR: covFile must be provided and must exist.")
  }

  ## read tables with ID handling -----------------------------------------------
  abd <- read_table_with_id(abdFile, id_col = sampleIDColinabdFile)
  cov <- read_table_with_id(covFile, id_col = sampleIDColincovFile)
  cat("Abundance and covariate files have been read\n")

  ## merge abundance and covariate files ----------------------------------------

  # merge_abd_cov(abd, cov) is assumed to:
  # - require both tables to contain 'IID'
  # - merge by 'IID' and keep IID as the first column
  merged <- merge_abd_cov(abd, cov)

        # data <- data.table::fread(phenoFile,
        #   header = T,
        #   stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        # )

  # # select required columns
  # data <- merged[, checkColList, drop = FALSE]
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

    # 1) Only calculate sequencing depth from count columns in abd (exclude sample ID column)
    abd_counts <- abd[, setdiff(colnames(abd), sampleIDColinabdFile), drop = FALSE]

    # Simple sanity check
    if (any(abd_counts < 0, na.rm = TRUE)) {
      warning("Abundance table contains negative values; 
              please make sure abdFile is raw counts when using SeqDepth.")
    }

    seqdepth <- rowSums(abd_counts, na.rm = TRUE)
    names(seqdepth) <- abd[[sampleIDColinabdFile]]
    
    data$SeqDepth <- seqdepth[match(data$IID, names(seqdepth))]

    # 2) Check for NA and non-positive values
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
  

  if (isRemoveZerosinPheno) {
    data <- data[which(data[, which(colnames(data) == phenoCol)] > 0), ]
    cat("Removing all zeros in the phenotype\n")
    if (nrow(data) == 0) {
      stop("ERROR: no samples are left after removing zeros in the phenotype\n")
    }
  }


  if (length(qCovarCol) > 0) {
    cat(qCovarCol, "are categorical covariates\n")
    if (!all(qCovarCol %in% covarColList)) {
      stop("ERROR! all covariates in qCovarCol must be in covarColList\n")
    } else {
      for (q in qCovarCol) {
        data[, q] <- as.factor(data[, q])
      }
    }
  }

  if (length(eCovarCol) > 0) {
    cat(eCovarCol, "are environmental covariates\n")
    if (!all(eCovarCol %in% covarColList)) {
      stop("ERROR! all covariates in eCovarCol must be in covarColList\n")
    }
  }


  if (length(covarColList) > 0) {
    cat(covarColList, "are sample-level covariates\n")
    if (!all(covarColList %in% covarColList)) {
      stop("ERROR! all covariates in covarColList must be in covarColList\n")
    }
  }

  if (FemaleOnly | MaleOnly) {
    if (!sexCol %in% colnames(data)) {
      stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
    } else {
      if (FemaleOnly) {
        data <- data[which(data[, which(colnames(data) ==
          sexCol)] == FemaleCode), ]
        if (nrow(data) == 0) {
          stop(
            "ERROR! no samples in the phenotype are coded as ",
            FemaleCode, " in the column ", sexCol,
            "\n"
          )
        }
      } else if (MaleOnly) {
        data <- data[which(data[, which(colnames(data) ==
          sexCol)] == MaleCode), ]
        if (nrow(data) == 0) {
          stop(
            "ERROR! no samples in the phenotype are coded as ",
            MaleCode, " in the column ", sexCol, "\n"
          )
        }
      }
    }
  }

  # construct the formula
  if (length(covarColList) > 0) {
    formula <- paste0(phenoCol, "~", paste0(covarColList,
      collapse = "+"
    ))
    hasCovariate <- TRUE
  } else {
    formula <- paste0(phenoCol, "~ 1")
    hasCovariate <- FALSE
  }
  
  if(isCovariateOffset & offsetCol != ""){
    formula <- paste0(formula, "+offset(log(", offsetCol, "))")
  }
  
  cat("formula is ", formula, "\n")
  

  # if (!is.null(sampleListwithGeno)) {
  #   dataMerge <- merge(mmat_nomissing, sampleListwithGeno,
  #     by.x = "IID", by.y = "IIDgeno"
  #   )
  #   dataMerge_sort <- dataMerge[with(dataMerge, order(IndexGeno)), ]
  # } else {
  #   dataMerge_sort <- mmat_nomissing
  #   dataMerge_sort$IIDgeno <- dataMerge_sort$IID
  # }
  


  if (!hasCovariate) {
    print("No covariate is includes so isCovariateOffset = FALSE")
    isCovariateOffset <- FALSE
  }

  # if (isCovariateTransform & hasCovariate) {
  #   cat("qr transformation has been performed on covariates\n")
  #   out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
  #   formulaNewList <- c(phenoCol, " ~ ", out.transform$Param.transform$X_name[1])
  #   if (length(out.transform$Param.transform$X_name) > 1) {
  #     for (i in c(2:length(out.transform$Param.transform$X_name))) {
  #       formulaNewList <- c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
  #     }
  #   }
  #   formulaNewList <- paste0(formulaNewList, collapse = "")
  #   formulaNewList <- paste0(formulaNewList, "-1")
  #   formula.new <- as.formula(paste0(formulaNewList, collapse = ""))
  #   data.new <- as.data.frame(cbind(out.transform$Y, out.transform$X1))
  #   colnames(data.new) <- c(phenoCol, out.transform$Param.transform$X_name)
  #   cat("colnames(data.new) is ", colnames(data.new), "\n")
  #   cat(
  #     "out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr),
  #     "\n"
  #   )

  #   if (length(offsetCol) > 0) {
  #     data.new <- cbind(data.new, dataMerge_sort[, which(colnames(dataMerge_sort) == offsetCol)])
  #     colnames(data.new)[ncol(data.new)] <- offsetCol
  #   }
  # } else {
  #   formula.new <- formula.null
  #   data.new <- dataMerge_sort
  #   out.transform <- NULL
  # }


  # if (useSparseGRMtoFitNULL) {
  #   if (!isSparseGRMIdentity) {
  #     getsubGRM_orig(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID)
  #   } else {
  #     sparseGRM <- Matrix:::sparseMatrix(i = as.vector(1:nrow(data)), j = as.vector(1:nrow(data)), x = rep(1, nrow(data)), symmetric = TRUE)
  #     rownames(sparseGRM) <- colnames(sparseGRM) <- data[[sampleIDColinphenoFile]]
  #   }
  #   gc()
  # }


  # print(dataMerge_sort$IID[1:200])
  # print(any(duplicated(dataMerge_sort$IID)))


  # if (longlCol != "") {
  #   covarianceIdxMat <- set_covarianceidx_Mat()
  # } else {
  #   covarianceIdxMat <- NULL
  # }


  # if (traitType == "binary") {
  #   stop("ERROR: This traitType is not supported in the current version.\n")
  # } else if (traitType == "quantitative") {
  #   stop("ERROR: This traitType is not supported in the current version.\n")
  # } else if (traitType == "count") {
  #   cat(phenoCol, " is a count trait\n")
  #   miny <- min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
  #   if (miny < 0) {
  #     stop("ERROR! phenotype value needs to be non-negative \n")
  #   }
  # } else if (traitType == "count_nb") {
  #   stop("ERROR: This traitType is not supported in the current version.\n")
  # }


  if (!skipModelFitting) {
    cat("Start fitting the NULL GLMM\n")
    t_begin <- proc.time()
    print(t_begin)

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
    if (traitType != "count_nb") {
      system.time(modglmm <- GMMAT::glmmkin(
        formula, 
        data = data,
        kins = grm_K,
        id = "IID", 
        family = poisson(link = "log")
        ))
      cat("glmmkin succeed!\n")
    } else {
      stop("ERROR: This traitType is not supported in the current version.\n")
    }

    
    # if (length(eCovarCol) > 0) {
    #   cat(eCovarCol, "are environmental covariates\n")
    #   modglmm$eMat <- data.new[, which(colnames(data.new) %in% eCovarCol), drop = F]
    #   for (em in 1:ncol(modglmm$eMat)) {
    #     modglmm$eMat[, em] <- (modglmm$eMat[, em] - mean(modglmm$eMat[, em])) / (sd(modglmm$eMat[, em]))
    #   }
    # }
    
    # if (length(covarColList) > 0) {
    #   cat(covarColList, "are sample-level covariates\n")
    
    #   covarColList <- c(covarColList, sampleCovarCol_q_names)
    #   modglmm$sampleXMat <- modglmm$X[, which(colnames(modglmm$X) %in% covarColList), drop = F]
    #   modglmm$sampleXMat <- cbind(modglmm$X[, 1], modglmm$sampleXMat)
    #   uniqsampleind <- which(!duplicated(modglmm$sampleID))
    #   modglmm$sampleXMat <- modglmm$sampleXMat[uniqsampleind, ]
    # }


    t_end <- proc.time()
    print(t_end)
    cat("t_end - t_begin, fitting the NULL model took\n")
    print(t_end - t_begin)

  } else {
    cat("Skip fitting the NULL GLMM\n")
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
  }

  save(modglmm, file = modelOut)

}
