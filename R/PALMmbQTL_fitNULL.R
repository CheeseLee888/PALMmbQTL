#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param phenoFile character. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarCol. All categorical covariates must also be included in covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary", "quantitative", "count" or "count_nb". By default, "count".
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the null model. All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex").
#' @param eCovarCol vector of characters. Covariates of environmental factors/cell context to be used in the null model. All covariates listed in eCovarCol must be also in covarColList,  e,g c("cellType").
#' @param sampleIDColinphenoFile character. Column name for the sample IDs in the phenotype file e.g. "IID".
#' @param cellIDColinphenoFile character. Column name for the cell IDs in the phenotype file e.g. "barcode".
#' @param tol numeric. The tolerance for fitting the null model to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500.
#' @param nThreads integer. Number of threads to be used. By default, 1
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkersForVarRatio integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Initial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is 0,0, the second element will be 0.5 for binary traits and the initial tau vector for quantitative traits is 1,0
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. By default, FALSE
#' @param traceCVcutoff numeric. The threshold for coefficient of variation (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variation (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkersForVarRatio will be increased by 10. By default, 0.001
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix. By default, ""
#' @param IsOverwriteVarianceRatioFile logical. Whether to overwrite the variance ratio file if the file exists. By default, FALSE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. By default, ""
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. By default, ""
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples that are included in the sparse GRM. By default, 2000
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM.
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null model. By default, FALSE
#' @param useSparseGRMforVarRatio logical. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE
#' @param minCovariateCount integer. If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues. By default, -1 (no covariates will be excluded)
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @param includeNonautoMarkersforVarRatio logical. Whether to allow for non-autosomal markers for variance ratio. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1'
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, ''
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, FALSE.
#' @param isStoreSigma logical. Whether to store sigma matrix. By default, FALSE. If number of individuals is greater than 10,000, this option may use large memory
#' @param isShrinkModelOutput logical. remove unnecessary objects for step2 from the model output. By default, FALSE.
#' @param isExportResiduals logical. export a residual vector. By default, FALSE.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV <- function(grmFile = "",
                               phenoFile = "",
                               phenoCol = "",
                               isRemoveZerosinPheno = FALSE,
                               traitType = "count",
                               invNormalize = FALSE,
                               covarColList = NULL,
                               qCovarCol = NULL,
                               eCovarCol = NULL,
                               sampleCovarCol = NULL,
                               offsetCol = NULL,
                               varWeightsCol = NULL,
                               longlCol = "",
                               sampleIDColinphenoFile = "IID",
                               cellIDColinphenoFile = "",
                               tol = 0.02,
                               maxiter = 20,
                               tolPCG = 1e-5,
                               maxiterPCG = 500,
                               nThreads = 1,
                               SPAcutoff = 2,
                               numMarkersForVarRatio = 30,
                               skipModelFitting = FALSE,
                               memoryChunk = 2,
                               tauInit = c(0, 0),
                               LOCO = FALSE,
                               isLowMemLOCO = FALSE,
                               traceCVcutoff = 0.0025,
                               ratioCVcutoff = 0.001,
                               outputPrefix = "",
                               outputPrefix_varRatio = "",
                               IsOverwriteVarianceRatioFile = FALSE,
                               sparseGRMFile = "",
                               sparseGRMSampleIDFile = "",
                               numRandomMarkerforSparseKin = 1000,
                               relatednessCutoff = 0.125,
                               isCateVarianceRatio = FALSE,
                               cateVarRatioIndexVec = NULL,
                               cateVarRatioMinMACVecExclude = c(10, 20.5),
                               cateVarRatioMaxMACVecInclude = c(20.5),
                               isCovariateTransform = TRUE,
                               isDiagofKinSetAsOne = FALSE,
                               minCovariateCount = -1,
                               minMAFforGRM = 0.01,
                               maxMissingRateforGRM = 0.15,
                               useSparseGRMtoFitNULL = FALSE,
                               useSparseGRMforVarRatio = FALSE,
                               includeNonautoMarkersforVarRatio = FALSE,
                               sexCol = "",
                               FemaleCode = 1,
                               FemaleOnly = FALSE,
                               MaleCode = 0,
                               MaleOnly = FALSE,
                               SampleIDIncludeFile = "",
                               isCovariateOffset = FALSE,
                               skipVarianceRatioEstimation = FALSE,
                               nrun = 30,
                               VmatFilelist = "",
                               VmatSampleFilelist = "",
                               VcellmatFilelist = "",
                               VcellmatSampleFilelist = "",
                               useGRMtoFitNULL = TRUE,
                               isStoreSigma = FALSE,
                               isShrinkModelOutput = FALSE,
                               isExportResiduals = FALSE) {
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


  if (!file.exists(phenoFile)) {
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  } else {
    if (longlCol == "") {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile)
    } else {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile, longlCol)
    }

    # if (cellIDColinphenoFile != "") {
    #   cat(cellIDColinphenoFile, "is the cell ID column\n")
    #   checkColList <- c(checkColList, cellIDColinphenoFile)
    # }

    if (length(offsetCol) > 0) {
      cat(offsetCol, "is the offset term\n")
      checkColList <- c(checkColList, offsetCol)
    }

    # if (length(varWeightsCol) > 0) {
    #   cat(varWeightsCol, " is the weights for variance\n")
    #   checkColList <- c(checkColList, varWeightsCol)
    # }


    ## check whether the phenotype file is large
    cmd <- paste0("du ", phenoFile, "| awk '{print $1}' > ", outputPrefix, "_", phenoCol, "_size_temp")
    system(cmd)
    datasize <- data.table::fread(paste0(outputPrefix, "_", phenoCol, "_size_temp"), header = F, data.table = F)
    isphenoFileLarge <- FALSE
    if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
      if (datasize[1, 1] > 200000) {
        isphenoFileLarge <- TRUE
      }
    } else {
      if (datasize[1, 1] > 500000) {
        isphenoFileLarge <- TRUE
      }
    }

    if (isphenoFileLarge) {
      catcmd <- ifelse(grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile), "gunzip -c ", "cat ")
      cmd <- paste0(catcmd, phenoFile, " | head -n 1 | sed 's/[\\t ]/\\n/g' | awk '{print $1\"\\t\"NR}' > ", outputPrefix, "_", phenoCol, "_lineNum_temp")
      system(cmd)

      checkColListDataFrame <- data.frame(colna = checkColList)
      phenoFilephenoCol_lineNum <- data.table::fread(paste0(outputPrefix, "_", phenoCol, "_lineNum_temp"), header = F, data.table = F)

      phenoFilephenoCol_lineNum_checkColList <- merge(checkColListDataFrame, phenoFilephenoCol_lineNum, by.x = 1, by.y = 1)

      write.table(phenoFilephenoCol_lineNum_checkColList[, 2], paste0(outputPrefix, "_", phenoCol, "_colnames_subset_temp"), quote = F, col.names = F, row.names = F)

      cmdb <- paste0(catcmd, phenoFile, " | cut -f $(tr '\\n' ',' < ", outputPrefix, "_", phenoCol, "_colnames_subset_temp | sed 's/,$//') > ", outputPrefix, "_", phenoCol, "_subcols_temp")
      system(cmdb)

      phenoFiletemp <- paste0(outputPrefix, "_", phenoCol, "_subcols_temp")

      data <- data.table::fread(phenoFiletemp,
        header = T,
        stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F
      )

      file.remove(paste0(outputPrefix, "_", phenoCol, "_colnames_subset_temp"))
      file.remove(paste0(outputPrefix, "_", phenoCol, "_lineNum_temp"))
      file.remove(paste0(outputPrefix, "_", phenoCol, "_subcols_temp"))
    } else { # !isphenoFileLarge

      if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
        data <- data.table::fread(
          cmd = paste0(
            "gunzip -c ",
            phenoFile
          ), header = T, stringsAsFactors = FALSE,
          colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        )
      } else {
        data <- data.table::fread(phenoFile,
          header = T,
          stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        )
      }
    }

    file.remove(paste0(outputPrefix, "_", phenoCol, "_size_temp"))


    if (isRemoveZerosinPheno) {
      data <- data[which(data[, which(colnames(data) == phenoCol)] > 0), ]
      cat("Removing all zeros in the phenotype\n")
      if (nrow(data) == 0) {
        stop("ERROR: no samples are left after removing zeros in the phenotype\n")
      }
    }



    if (SampleIDIncludeFile != "") {
      if (!file.exists(SampleIDIncludeFile)) {
        stop("ERROR! SampleIDIncludeFile ", SampleIDIncludeFile, " does not exsit\n")
      } else {
        sampleIDInclude <- data.table::fread(SampleIDIncludeFile, header = F, stringsAsFactors = FALSE, colClasses = c("character"), data.table = F)
        sampleIDInclude <- as.vector(sampleIDInclude[!duplicated(sampleIDInclude), ])
        cat(length(sampleIDInclude), " non-duplicated sample IDs were found in SampleIDIncludeFile\n")
        data <- data[which(as.vector(data[, which(colnames(data) == sampleIDColinphenoFile)]) %in% sampleIDInclude), , drop = F]
        cat(nrow(data), " samples in sampleIDInclude have non-missing phenotypes and covariates\n")
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


    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
      if (!all(sampleCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in sampleCovarCol must be in covarColList\n")
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
    
    # formula.null <- as.formula(formula)
    # mmat <- model.matrix(formula.null, data, na.action = NULL)
    # mmat <- cbind(mmat, data[, which(colnames(data) == phenoCol), drop = F])
    # colnames(mmat)[ncol(mmat)] <- phenoCol

    # if (length(sampleCovarCol) > 0) {
    #   cat(sampleCovarCol, "are sample-level covariates\n")
    #   # check which sample-level covariates are categorical and record the names after factorizing in the data frame
    #   if (length(qCovarCol) > 0) {
    #     if (any(sampleCovarCol %in% qCovarCol)) {
    #       sampleCovarCol_q <- sampleCovarCol[which(sampleCovarCol %in% qCovarCol)]
    #       formula_sq <- paste0("~", paste0(sampleCovarCol_q, collapse = "+"))
    #       formula_sq.null <- as.formula(formula_sq)
    #       mmat_sq <- model.matrix(formula_sq.null, data, na.action = NULL)
    #       sampleCovarCol_q_names <- colnames(mmat_sq)[-1]
    #       rm(mmat_sq)
    #     } else {
    #       sampleCovarCol_q_names <- NULL
    #     }
    #   } else {
    #     sampleCovarCol_q_names <- NULL
    #   }
    # }

    # coln <- 1
    # if (length(offsetCol) > 0) {
    #   mmat <- cbind(mmat, data[, which(colnames(data) == offsetCol), drop = F])
    #   colnames(mmat)[ncol(mmat)] <- offsetCol
    #   coln <- coln + 1
    # }

    # if (length(varWeightsCol) > 0) {
    #   mmat <- cbind(mmat, data[, which(colnames(data) == varWeightsCol), drop = F])
    #   colnames(mmat)[ncol(mmat)] <- varWeightsCol

    #   coln <- coln + 1
    # }

    # if (length(covarColList) > 0) {
    #   if (length(qCovarCol) > 0) {
    #     covarColList <- colnames(mmat)[2:(ncol(mmat) - coln)]
    #     formula <- paste0(phenoCol, "~", paste0(covarColList, collapse = "+"))
    #     formula.null <- as.formula(formula)
    #   }
    # }

    # mmat$IID <- data[, which(sampleIDColinphenoFile == colnames(data))]
    # if (cellIDColinphenoFile != "") {
    #   mmat$barcode <- data[, which(cellIDColinphenoFile == colnames(data))]
    # }
    # if (longlCol != "") {
    #   mmat$longlVar <- data[, which(longlCol == colnames(data))]
    # }

    # mmat_nomissing <- mmat[complete.cases(mmat), ]
    # mmat_nomissing$IndexPheno <- seq(1, nrow(mmat_nomissing),
    #   by = 1
    # )
    # cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")

    # if (length(varWeightsCol) > 0) {
    #   varWeights <- mmat_nomissing[, which(colnames(mmat_nomissing) == varWeightsCol)]
    # } else {
    #   varWeights <- NULL
    # }
    # if (sparseGRMSampleIDFile != "") {
    #   sampleListwithGenov0 <- data.table::fread(sparseGRMSampleIDFile,
    #     header = F, , colClasses = c("character"), data.table = F
    #   )
    #   colnames(sampleListwithGenov0) <- c("IIDgeno")
    #   cat(length(sampleListwithGenov0$IIDgeno), " samples are in the sparse GRM\n")
    #   mmat_nomissing <- mmat_nomissing[which(mmat_nomissing$IID %in% sampleListwithGenov0$IIDgeno), ]
    #   cat(nrow(mmat_nomissing), " samples who have non-missing phenotypes are also in the sparse GRM\n")
    # }


    # if (longlCol == "") {
    #   if (any(duplicated(mmat_nomissing$IID))) {
    #     cat("Duplicated sample IDs are detected in the phenotype file. Assuming repeated measurements\n")
    #   }
    # } else {
    #   cat("Longitudinal variable ", longlCol, " is specified\n")
    #   if (!any(duplicated(mmat_nomissing$IID))) {
    #     stop("No duplicated sample IDs are detected in the phenotype file\n")
    #   }
    # }


    # if (!is.null(sampleListwithGeno)) {
    #   dataMerge <- merge(mmat_nomissing, sampleListwithGeno,
    #     by.x = "IID", by.y = "IIDgeno"
    #   )
    #   dataMerge_sort <- dataMerge[with(dataMerge, order(IndexGeno)), ]
    # } else {
    #   dataMerge_sort <- mmat_nomissing
    #   dataMerge_sort$IIDgeno <- dataMerge_sort$IID
    # }

    # print("Test")
    # print(head(dataMerge_sort))

    # rm(mmat)
    # rm(mmat_nomissing)
    # gc()
    # isSparseGRMIdentity <- FALSE
    # if (useGRMtoFitNULL) {
    #   indicatorGenoSamplesWithPheno <- (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)

    #   if (length(unique(dataMerge_sort$IIDgeno)) < length(unique(sampleListwithGeno$IIDgeno))) {
    #     cat(
    #       length(unique(sampleListwithGeno$IIDgeno)) - length(unique(dataMerge_sort$IIDgeno)),
    #       " samples in geno file do not have phenotypes\n"
    #     )
    #   }
    #   cat(length(unique(dataMerge_sort$IIDgeno)), " samples will be used for analysis\n")
    # } else {
    #   indicatorGenoSamplesWithPheno <- rep(TRUE, nrow(dataMerge_sort))
    # }

    # if (any(duplicated(dataMerge_sort$IID))) {
    #   cat(nrow(dataMerge_sort), " observations will be used for analysis\n")
    #   set_I_mat_inR(dataMerge_sort$IID)
    #   if (longlCol != "") {
    #     set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
    #   }
    # } else {
    #   if (!useGRMtoFitNULL) {
    #     cat("No duplicated IDs are observed in the phenotype file, so the identity matrix will be used as a sparse GRM will be used to fit the null model\n")
    #     isSparseGRMIdentity <- TRUE
    #     useSparseGRMtoFitNULL <- TRUE
    #     useGRMtoFitNULL <- TRUE
    #   }
    # }
    # set_useGRMtoFitNULL(useGRMtoFitNULL)
  }


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


  # if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
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

    pheno_id <- data[[sampleIDColinphenoFile]]
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
      cat(grm_id[1:5])
      cat("pheno id:\n")
      cat(pheno_id[1:5])

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
        id = sampleIDColinphenoFile, 
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
    
    # if (length(sampleCovarCol) > 0) {
    #   cat(sampleCovarCol, "are sample-level covariates\n")
    
    #   sampleCovarCol <- c(sampleCovarCol, sampleCovarCol_q_names)
    #   modglmm$sampleXMat <- modglmm$X[, which(colnames(modglmm$X) %in% sampleCovarCol), drop = F]
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

  fastSave(modglmm, file = modelOut)

}
