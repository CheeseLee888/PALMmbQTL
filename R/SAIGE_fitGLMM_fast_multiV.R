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
fitNULLGLMM_multiV <- function(plinkFile = "",
                               bedFile = "",
                               bimFile = "",
                               famFile = "",
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

  if (plinkFile != "") {
    bimFile <- paste0(plinkFile, ".bim")
    bedFile <- paste0(plinkFile, ".bed")
    famFile <- paste0(plinkFile, ".fam")
  }
  setgenoNULL()

  if (!useGRMtoFitNULL) {
    useSparseGRMtoFitNULL <- FALSE
    useSparseGRMforVarRatio <- FALSE
    cat("No GRM will be used to fit the NULL model and nThreads is set to 1\n")
  }


  if (useSparseGRMtoFitNULL & bedFile == "") {
    cat("Sparse GRM is used to fit the null model and plink file is not specified, so variance ratios won't be estimated\n")
    skipVarianceRatioEstimation <- TRUE
  }

  if (!skipVarianceRatioEstimation) {
    SPAGMMATOut <- paste0(outputPrefix, "_", numMarkersForVarRatio, "markers.SAIGE.results.txt")

    if (outputPrefix_varRatio == "") {
      outputPrefix_varRatio <- outputPrefix
    }
    varRatioFile <- paste0(outputPrefix_varRatio, ".varianceRatio.txt")

    if (!file.exists(varRatioFile)) {
      file.create(varRatioFile, showWarnings = TRUE)
    } else {
      if (!IsOverwriteVarianceRatioFile) {
        stop(
          "WARNING: The variance ratio file ", varRatioFile,
          " already exists. The new variance ratios will be output to ",
          varRatioFile, ". In order to avoid overwriting the file, please remove the ",
          varRatioFile, " or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s). Otherwise, specify --IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten with new variance ratio(s)\n"
        )
      } else {
        cat("The variance ratio file ", varRatioFile, " already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten\n")
      }
    }
  } else {
    cat("Variance ratio estimation will be skipped\n.")
    useSparseGRMforVarRatio <- FALSE
  }


  if (useSparseGRMtoFitNULL) {
    cat("Leave-one-chromosome-out is not applied\n")
  }



  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!file.exists(sparseGRMFile)) {
      stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
    }
    if (!file.exists(sparseGRMSampleIDFile)) {
      stop(
        "sparseGRMSampleIDFile ", sparseGRMSampleIDFile,
        " does not exist!"
      )
    }
  }


  if (nThreads > 1) {
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads will be used ", "\n")
  }
  set_g_omp_num_threads(nThreads)

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


  sampleListwithGeno <- NULL
  if ((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | !skipVarianceRatioEstimation) {
    if (!file.exists(bedFile)) {
      stop("ERROR! bed file does not exsit\n")
    }
    if (!file.exists(bimFile)) {
      stop("ERROR! bim file does not exsit\n")
    }


    if (!file.exists(famFile)) {
      stop("ERROR! fam file does not exsit\n")
    } else {
      sampleListwithGenov0 <- data.table::fread(famFile, header = F, colClasses = list(character = 1:4), data.table = FALSE)
      colnames(sampleListwithGenov0) <- c(
        "FIDgeno", "IIDgeno",
        "father", "mother", "sex", "phe"
      )
      sampleListwithGeno <- NULL
      sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
      sampleListwithGeno <- data.frame(sampleListwithGeno)
      sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
        by = 1
      )
      cat(nrow(sampleListwithGeno), " samples have genotypes\n")
    }
  } else {
    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
      sampleListwithGenov0 <- data.table::fread(sparseGRMSampleIDFile,
        header = F, , colClasses = c("character"), data.table = F
      )
      colnames(sampleListwithGenov0) <- c("IIDgeno")
      sampleListwithGeno <- NULL
      sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
      sampleListwithGeno <- data.frame(sampleListwithGeno)
      sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
        by = 1
      )
      cat(nrow(sampleListwithGeno), " samples are in the sparse GRM\n")
    }
  }
  if (!file.exists(phenoFile)) {
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  } else {
    if (longlCol == "") {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile)
    } else {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile, longlCol)
    }

    if (cellIDColinphenoFile != "") {
      cat(cellIDColinphenoFile, "is the cell ID column\n")
      checkColList <- c(checkColList, cellIDColinphenoFile)
    }

    if (length(offsetCol) > 0) {
      cat(offsetCol, "is the offset term\n")
      checkColList <- c(checkColList, offsetCol)
    }

    if (length(varWeightsCol) > 0) {
      cat(varWeightsCol, " is the weights for variance\n")
      checkColList <- c(checkColList, varWeightsCol)
    }


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


    print("HERERE2")

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
    
    formula.null <- as.formula(formula)
    mmat <- model.matrix(formula.null, data, na.action = NULL)
    mmat <- cbind(mmat, data[, which(colnames(data) == phenoCol), drop = F])
    colnames(mmat)[ncol(mmat)] <- phenoCol

    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
      # check which sample-level covariates are categorical and record the names after factorizing in the data frame
      if (length(qCovarCol) > 0) {
        if (any(sampleCovarCol %in% qCovarCol)) {
          sampleCovarCol_q <- sampleCovarCol[which(sampleCovarCol %in% qCovarCol)]
          formula_sq <- paste0("~", paste0(sampleCovarCol_q, collapse = "+"))
          formula_sq.null <- as.formula(formula_sq)
          mmat_sq <- model.matrix(formula_sq.null, data, na.action = NULL)
          sampleCovarCol_q_names <- colnames(mmat_sq)[-1]
          rm(mmat_sq)
        } else {
          sampleCovarCol_q_names <- NULL
        }
      } else {
        sampleCovarCol_q_names <- NULL
      }
    }



    coln <- 1
    if (length(offsetCol) > 0) {
      mmat <- cbind(mmat, data[, which(colnames(data) == offsetCol), drop = F])
      colnames(mmat)[ncol(mmat)] <- offsetCol
      coln <- coln + 1
    }

    if (length(varWeightsCol) > 0) {
      mmat <- cbind(mmat, data[, which(colnames(data) == varWeightsCol), drop = F])
      colnames(mmat)[ncol(mmat)] <- varWeightsCol

      coln <- coln + 1
    }



    if (length(covarColList) > 0) {
      if (length(qCovarCol) > 0) {
        covarColList <- colnames(mmat)[2:(ncol(mmat) - coln)]
        formula <- paste0(phenoCol, "~", paste0(covarColList, collapse = "+"))
        formula.null <- as.formula(formula)
      }
    }

    mmat$IID <- data[, which(sampleIDColinphenoFile == colnames(data))]
    if (cellIDColinphenoFile != "") {
      mmat$barcode <- data[, which(cellIDColinphenoFile == colnames(data))]
    }
    if (longlCol != "") {
      mmat$longlVar <- data[, which(longlCol == colnames(data))]
    }

    mmat_nomissing <- mmat[complete.cases(mmat), ]
    mmat_nomissing$IndexPheno <- seq(1, nrow(mmat_nomissing),
      by = 1
    )
    cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")

    if (length(varWeightsCol) > 0) {
      varWeights <- mmat_nomissing[, which(colnames(mmat_nomissing) == varWeightsCol)]
    } else {
      varWeights <- NULL
    }
    if (sparseGRMSampleIDFile != "") {
      sampleListwithGenov0 <- data.table::fread(sparseGRMSampleIDFile,
        header = F, , colClasses = c("character"), data.table = F
      )
      colnames(sampleListwithGenov0) <- c("IIDgeno")
      cat(length(sampleListwithGenov0$IIDgeno), " samples are in the sparse GRM\n")
      mmat_nomissing <- mmat_nomissing[which(mmat_nomissing$IID %in% sampleListwithGenov0$IIDgeno), ]
      cat(nrow(mmat_nomissing), " samples who have non-missing phenotypes are also in the sparse GRM\n")
    }


    if (longlCol == "") {
      if (any(duplicated(mmat_nomissing$IID))) {
        cat("Duplicated sample IDs are detected in the phenotype file. Assuming repeated measurements\n")
      }
    } else {
      cat("Longitudinal variable ", longlCol, " is specified\n")
      if (!any(duplicated(mmat_nomissing$IID))) {
        stop("No duplicated sample IDs are detected in the phenotype file\n")
      }
    }


    if (!is.null(sampleListwithGeno)) {
      dataMerge <- merge(mmat_nomissing, sampleListwithGeno,
        by.x = "IID", by.y = "IIDgeno"
      )
      dataMerge_sort <- dataMerge[with(dataMerge, order(IndexGeno)), ]
    } else {
      dataMerge_sort <- mmat_nomissing
      dataMerge_sort$IIDgeno <- dataMerge_sort$IID
    }

    print("Test")
    print(head(dataMerge_sort))

    rm(mmat)
    rm(mmat_nomissing)
    gc()
    isSparseGRMIdentity <- FALSE
    if (useGRMtoFitNULL) {
      indicatorGenoSamplesWithPheno <- (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)

      if (length(unique(dataMerge_sort$IIDgeno)) < length(unique(sampleListwithGeno$IIDgeno))) {
        cat(
          length(unique(sampleListwithGeno$IIDgeno)) - length(unique(dataMerge_sort$IIDgeno)),
          " samples in geno file do not have phenotypes\n"
        )
      }
      cat(length(unique(dataMerge_sort$IIDgeno)), " samples will be used for analysis\n")
    } else {
      indicatorGenoSamplesWithPheno <- rep(TRUE, nrow(dataMerge_sort))
    }

    if (any(duplicated(dataMerge_sort$IID))) {
      cat(nrow(dataMerge_sort), " observations will be used for analysis\n")
      set_I_mat_inR(dataMerge_sort$IID)
      if (longlCol != "") {
        set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
      }
    } else {
      if (!useGRMtoFitNULL) {
        cat("No duplicated IDs are observed in the phenotype file, so the identity matrix will be used as a sparse GRM will be used to fit the null model\n")
        isSparseGRMIdentity <- TRUE
        useSparseGRMtoFitNULL <- TRUE
        useGRMtoFitNULL <- TRUE
      }
    }
    set_useGRMtoFitNULL(useGRMtoFitNULL)
  }



  print("Test3")
  print(head(dataMerge_sort))

  if (traitType == "quantitative" & invNormalize) {
    stop("ERROR: This traitType is not supported in the current version.\n")
  }
  print("Test4")
  print(head(dataMerge_sort))
  if (traitType == "binary" & (length(covarColList) > 0)) {
    stop("ERROR: This traitType is not supported in the current version.\n")
  }
  if (!hasCovariate) {
    print("No covariate is includes so isCovariateOffset = FALSE")
    isCovariateOffset <- FALSE
  }

  if (isCovariateTransform & hasCovariate) {
    cat("qr transformation has been performed on covariates\n")
    out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
    formulaNewList <- c(phenoCol, " ~ ", out.transform$Param.transform$X_name[1])
    if (length(out.transform$Param.transform$X_name) > 1) {
      for (i in c(2:length(out.transform$Param.transform$X_name))) {
        formulaNewList <- c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
      }
    }
    formulaNewList <- paste0(formulaNewList, collapse = "")
    formulaNewList <- paste0(formulaNewList, "-1")
    formula.new <- as.formula(paste0(formulaNewList, collapse = ""))
    data.new <- as.data.frame(cbind(out.transform$Y, out.transform$X1))
    colnames(data.new) <- c(phenoCol, out.transform$Param.transform$X_name)
    cat("colnames(data.new) is ", colnames(data.new), "\n")
    cat(
      "out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr),
      "\n"
    )

    if (length(offsetCol) > 0) {
      data.new <- cbind(data.new, dataMerge_sort[, which(colnames(dataMerge_sort) == offsetCol)])
      colnames(data.new)[ncol(data.new)] <- offsetCol
    }
  } else {
    formula.new <- formula.null
    data.new <- dataMerge_sort
    out.transform <- NULL
  }

  if (traitType == "binary") {
    stop("ERROR: This traitType is not supported in the current version.\n")
  } else if (traitType == "quantitative") {
    stop("ERROR: This traitType is not supported in the current version.\n")
  } else if (traitType == "count") {
    if (length(offsetCol) == 0) {
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
    }
  } else if (traitType == "count_nb") {
    stop("ERROR: This traitType is not supported in the current version.\n")
  }
  mmat <- model.matrix(formula.new, data = data.new, na.action = NULL)

  if (isCovariateOffset) {
    # covoffset <- mmat[, -1, drop = F] %*% modwitcov$coefficients[-1]
    print("isCovariateOffset=TRUE, so fixed effects coefficnets won't be estimated.")
    formula.new.withCov <- formula.new
    formula_nocov <- paste0(phenoCol, "~ 1")
    formula.new <- as.formula(formula_nocov)
    hasCovariate <- FALSE
  } else {
    covoffset <- rep(0, nrow(data.new))
  }

  # data.new$covoffset <- covoffset


  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!isSparseGRMIdentity) {
      getsubGRM_orig(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID)
    } else {
      sparseGRM <- Matrix:::sparseMatrix(i = as.vector(1:nrow(data)), j = as.vector(1:nrow(data)), x = rep(1, nrow(data)), symmetric = TRUE)
      rownames(sparseGRM) <- colnames(sparseGRM) <- data[[sampleIDColinphenoFile]]
    }
    gc()
  }

  # allow for multiple variance components
  set_Vmat_vec_orig(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID)


  print(dataMerge_sort$IID[1:200])
  print(any(duplicated(dataMerge_sort$IID)))

  if (any(duplicated(dataMerge_sort$IID))) {
    print("HERE")
    if (longlCol == "") {
      print("HERE1")
      print(useGRMtoFitNULL)
      if (useGRMtoFitNULL) {
        print("HERE2")
      }
    }
  }





  if (longlCol != "") {
    covarianceIdxMat <- set_covarianceidx_Mat()
  } else {
    covarianceIdxMat <- NULL
  }

  if (!skipVarianceRatioEstimation) {
    isVarianceRatioinGeno <- TRUE
    if (isCateVarianceRatio) {
      minMAC_varRatio <- min(cateVarRatioMinMACVecExclude)
      maxMAC_varRatio <- max(cateVarRatioMaxMACVecInclude)
      cat("Categorical variance ratios will be estimated. Please make sure there are at least 200 markers in each MAC category.\n")
    } else {
      minMAC_varRatio <- 20
      maxMAC_varRatio <- -1 # will randomly select markers from the plink file and leave them out when constructing GRM
    }
    setminMAC_VarianceRatio(minMAC_varRatio, maxMAC_varRatio, isVarianceRatioinGeno)
  }

  # set up parameters
  if (minMAFforGRM > 0) {
    cat(
      "Markers in the Plink file with MAF < ", minMAFforGRM,
      " will be removed before constructing GRM\n"
    )
  }
  if (maxMissingRateforGRM > 0) {
    cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
  }

  setminMAFforGRM(minMAFforGRM)
  setmaxMissingRateforGRM(maxMissingRateforGRM)



  if (traitType == "binary") {
    stop("ERROR: This traitType is not supported in the current version.\n")
  } else if (traitType == "quantitative") {
    stop("ERROR: This traitType is not supported in the current version.\n")
  } else if (traitType == "count") {
    cat(phenoCol, " is a count trait\n")
    miny <- min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
    if (miny < 0) {
      stop("ERROR! phenotype value needs to be non-negative \n")
    }

    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      }
      Xorig <- NULL
    } else {
      # gc()
      # if (length(offsetCol) == 0) {
      # } else {
      #   offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
      # }
    }
  } else if (traitType == "count_nb") {
    stop("ERROR: This traitType is not supported in the current version.\n")
  }


  obj.noK <- NULL


  # print("isStoreSigma")
  # print(isStoreSigma)

  if (!skipModelFitting) {
    cat("Start fitting the NULL GLMM\n")
    t_begin <- proc.time()
    print(t_begin)



    set_isSparseGRM(useSparseGRMtoFitNULL)
    set_useGRMtoFitNULL(useGRMtoFitNULL)

    # Core step1 for PALM-mbQTL
    if (traitType != "count_nb") {
      system.time(modglmm <- GMMAT::glmmkin(
        formula, 
        data = data,
        kins = sparseGRM,
        id = sampleIDColinphenoFile, 
        family = poisson(link = "log")
        ))
      cat("glmmkin succeed!\n")
    } else {
      stop("ERROR: This traitType is not supported in the current version.\n")
    }

    
    if (length(eCovarCol) > 0) {
      cat(eCovarCol, "are environmental covariates\n")
      modglmm$eMat <- data.new[, which(colnames(data.new) %in% eCovarCol), drop = F]
      for (em in 1:ncol(modglmm$eMat)) {
        modglmm$eMat[, em] <- (modglmm$eMat[, em] - mean(modglmm$eMat[, em])) / (sd(modglmm$eMat[, em]))
      }
    }
    
    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
    
      sampleCovarCol <- c(sampleCovarCol, sampleCovarCol_q_names)
      modglmm$sampleXMat <- modglmm$X[, which(colnames(modglmm$X) %in% sampleCovarCol), drop = F]
      modglmm$sampleXMat <- cbind(modglmm$X[, 1], modglmm$sampleXMat)
      uniqsampleind <- which(!duplicated(modglmm$sampleID))
      modglmm$sampleXMat <- modglmm$sampleXMat[uniqsampleind, ]
    }


    t_end <- proc.time()
    print(t_end)
    cat("t_end - t_begin, fitting the NULL model took\n")
    print(t_end - t_begin)


    if (bedFile != "") {
      subSampleInGeno <- dataMerge_sort$IndexGeno
      if (is.null(dataMerge_sort$IndexGeno)) {
        subSampleInGeno <- dataMerge_sort$IndexPheno
      }
    
      print(subSampleInGeno[1:1000])
      print(head(dataMerge_sort))
      print("HEREHRE")
    
      subSampleInGeno_unique <- subSampleInGeno[!duplicated(subSampleInGeno)]
    
      setgeno(bedFile, bimFile, famFile, subSampleInGeno_unique, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    }
  } else {
    cat("Skip fitting the NULL GLMM\n")
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
    load(modelOut)

    # need check
    subSampleInGeno <- dataMerge_sort$IndexGeno
    if (is.null(dataMerge_sort$IndexGeno)) {
      subSampleInGeno <- dataMerge_sort$IndexPheno
    }

    print(subSampleInGeno[1:1000])
    print(head(dataMerge_sort))
    print("HEREHRE")

    subSampleInGeno_unique <- subSampleInGeno[!duplicated(subSampleInGeno)]

    setgeno(bedFile, bimFile, famFile, subSampleInGeno_unique, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)



    if (any(duplicated(modglmm$sampleID))) {
      set_I_mat_inR(modglmm$sampleID)
    }
    set_dup_sample_index(as.numeric(factor(modglmm$sampleID, levels = unique(modglmm$sampleID))))
  }

  if (!skipVarianceRatioEstimation) {
    cat("Start estimating variance ratios\n")
    extractVarianceRatio_multiV(
      obj.glmm.null = modglmm,
      obj.glm.null = fit0, maxiterPCG = maxiterPCG,
      tolPCG = tolPCG, numMarkers = numMarkersForVarRatio, varRatioOutFile = varRatioFile,
      ratioCVcutoff = ratioCVcutoff, testOut = SPAGMMATOut,
      bedFile = bedFile, bimFile = bimFile, famFile = famFile, chromosomeStartIndexVec = chromosomeStartIndexVec,
      chromosomeEndIndexVec = chromosomeEndIndexVec,
      isCateVarianceRatio = isCateVarianceRatio, cateVarRatioIndexVec = cateVarRatioIndexVec,
      useSparseGRMforVarRatio = useSparseGRMforVarRatio, sparseGRMFile = sparseGRMFile,
      sparseGRMSampleIDFile = sparseGRMSampleIDFile,
      numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
      relatednessCutoff = relatednessCutoff, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL,
      nThreads = nThreads, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      minMAFforGRM = minMAFforGRM, isDiagofKinSetAsOne = isDiagofKinSetAsOne,
      includeNonautoMarkersforVarRatio = includeNonautoMarkersforVarRatio, isStoreSigma = isStoreSigma, useGRMtoFitNULL = useGRMtoFitNULL
    )
  } else {
    cat("Skip estimating variance ratios\n")
  }
  closeGenoFile_plink()

  fastSave(modglmm, file = modelOut)

  if (isExportResiduals) {
    b <- as.numeric(factor(modglmm$sampleID, levels = unique(modglmm$sampleID)))
    I_mat <- 1.0 * Matrix::sparseMatrix(i = seq_along(b), j = b, x = rep(1, length(b)))
    res_sample <- as.vector(t(I_mat) %*% modglmm$residuals)
    data.table::fwrite(
      data.frame(sampleID = unique(modglmm$sampleID), residuals = res_sample),
      paste0(outputPrefix, ".sample.residuals.txt"),
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      na = "NA"
    )
    data.table::fwrite(
      data.frame(barcode = modglmm$barcode, residuals = modglmm$residuals),
      paste0(outputPrefix, ".residuals.txt"),
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      na = "NA"
    )
  }
}



extractVarianceRatio_multiV <- function(obj.glmm.null,
                                        obj.glm.null,
                                        maxiterPCG = 500,
                                        tolPCG = 0.01,
                                        numMarkers,
                                        varRatioOutFile,
                                        ratioCVcutoff,
                                        testOut,
                                        bedFile,
                                        bimFile,
                                        famFile,
                                        chromosomeStartIndexVec,
                                        chromosomeEndIndexVec,
                                        isCateVarianceRatio,
                                        cateVarRatioIndexVec,
                                        useSparseGRMforVarRatio,
                                        sparseGRMFile,
                                        sparseGRMSampleIDFile,
                                        numRandomMarkerforSparseKin,
                                        relatednessCutoff,
                                        useSparseGRMtoFitNULL,
                                        nThreads,
                                        cateVarRatioMinMACVecExclude,
                                        cateVarRatioMaxMACVecInclude,
                                        minMAFforGRM,
                                        isDiagofKinSetAsOne,
                                        includeNonautoMarkersforVarRatio,
                                        isStoreSigma = FALSE,
                                        useGRMtoFitNULL = TRUE) {
  obj.noK <- obj.glmm.null$obj.noK
  if (file.exists(testOut)) {
    file.remove(testOut)
  }
  bimPlink <- data.table::fread(bimFile, header = F, data.table = FALSE)
  if (sum(sapply(bimPlink[, 1], is.numeric)) != nrow(bimPlink)) {
    stop("ERROR: chromosome column in plink bim file is no numeric!\n")
  }

  print(family)
  eta <- obj.glmm.null$linear.predictors
  mu <- obj.glmm.null$fitted.values
  mu.eta <- family$mu.eta(eta)

  var_weights <- obj.glmm.null$varWeights
  sqrtW <- mu.eta / sqrt(family$variance(mu))


  print("mu[1:20]")
  print(mu[1:20])
  W <- sqrtW^2 ## (mu*(1-mu) for binary)
  W <- W * var_weights
  print("mu[1:20]")
  print(W[1:20])

  X <- obj.glmm.null$X


  set_isSparseGRM(useSparseGRMtoFitNULL)
  set_useGRMtoFitNULL(useGRMtoFitNULL)


  Sigma_iX_noLOCO <- getSigma_X_multiV(W, tauVecNew, X, maxiterPCG, tolPCG, LOCO = FALSE)


  y <- obj.glmm.null$y

  if (any(duplicated(obj.glmm.null$sampleID))) {
    dupSampleIndex <- as.numeric(factor(obj.glmm.null$sampleID, levels = unique(obj.glmm.null$sampleID)))
  }


  ## randomize the marker orders to be tested
  if (FALSE) {
    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
      sparseSigma <- getSparseSigma(
        bedFile = bedFile, bimFile = bimFile, famFile = famFile,
        outputPrefix = varRatioOutFile,
        sparseGRMFile = sparseGRMFile,
        sparseGRMSampleIDFile = sparseGRMSampleIDFile,
        numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
        relatednessCutoff = relatednessCutoff,
        minMAFforGRM = minMAFforGRM,
        nThreads = nThreads,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        obj.glmm.null = obj.glmm.null,
        W = W, tauVecNew = tauVecNew
      )
      if (length(tauVecNew) > 2) {
        sparseSigma <- sparseSigma + getProdTauKmat(tauVecNew[3:length(tauVecNew)])
      }
    }
  }


  mMarkers <- gettotalMarker()
  listOfMarkersForVarRatio <- list()
  MACvector <- getMACVec()
  isVarianceRatioinGeno <- getIsVarRatioGeno()

  if (isVarianceRatioinGeno) {
    MACvector_forVarRatio <- getMACVec_forVarRatio()
    Indexvector_forVarRatio <- getIndexVec_forVarRatio()
    cat("length(MACvector): ", length(MACvector), "\n")
    cat("length(MACvector_forVarRatio): ", length(MACvector_forVarRatio), "\n")

    if (length(MACvector_forVarRatio) > 0) {
      MACdata <- data.frame(MACvector = MACvector_forVarRatio, geno_ind = rep(1, length(MACvector_forVarRatio)), indexInGeno = seq(1, length(MACvector_forVarRatio)))
    } else {
      stop("No markers were found for variance ratio estimation. Please make sure there are at least 200 markers in each MAC category\n")
    }
  } else {
    MACdata <- data.frame(MACvector = MACvector, geno_ind = rep(0, length(MACvector)), indexInGeno = seq(1, length(MACvector)))
  }

  if (!isCateVarianceRatio) {
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex <- 1:nrow(MACdata)
    listOfMarkersForVarRatio[[1]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec <- c(1)
  } else {
    cat("Categorical variance ratios will be estimated.\n")

    if (is.null(cateVarRatioIndexVec)) {
      cateVarRatioIndexVec <- rep(1, length(cateVarRatioMinMACVecExclude))
    }
    numCate <- length(cateVarRatioIndexVec)
    for (i in 1:(numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[i] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[i])
      listOfMarkersForVarRatio[[i]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    }

    if (length(cateVarRatioMaxMACVecInclude) == (numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate])
    } else {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }

    listOfMarkersForVarRatio[[numCate]] <- sample(MACindex, size = length(MACindex), replace = FALSE)

    for (k in 1:length(cateVarRatioIndexVec)) {
      if (k <= length(cateVarRatioIndexVec) - 1) {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], "\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      } else {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }
  }



  b <- as.numeric(factor(obj.glmm.null$sampleID, levels = unique(obj.glmm.null$sampleID)))
  I_mat <- Matrix::sparseMatrix(i = seq_along(b), j = b, x = rep(1, length(b)))

  freqVec <- getAlleleFreqVec()
  Nnomissing <- length(mu)
  varRatioTable <- NULL




  Vsample0 <- as.vector(crossprod(obj.noK$V, I_mat))
  Xsample0 <- obj.glmm.null$sampleXMat
  XVsample0 <- t(Xsample0 * Vsample0)
  XVXsample0 <- crossprod(Xsample0, t(XVsample0))
  XVXsample_inv0 <- solve(XVXsample0)
  XXVXsample_inv0 <- Xsample0 %*% XVXsample_inv0
  XVX_inv_XVsample0 <- XXVXsample_inv0 * Vsample0

  for (k in 1:length(listOfMarkersForVarRatio)) {
    if (cateVarRatioIndexVec[k] == 1) {
      numMarkers0 <- numMarkers
      varRatio_sparseGRM_vec <- NULL
      varRatio_NULL_vec <- NULL
      varRatio_NULL_sample_vec <- NULL
      varRatio_NULL_noXadj_vec <- NULL


      if (!is.null(obj.glmm.null$eMat)) {
        varRatio_NULL_eg_mat <- NULL
        varRatio_NULL_eg_vec <- NULL
        varRatio_sparse_eg_mat <- NULL
        varRatio_sparse_eg_vec <- NULL
      }

      indexInMarkerList <- 1
      numTestedMarker <- 0
      ratioCV <- ratioCVcutoff + 0.1

      while (ratioCV > ratioCVcutoff) {
        while (numTestedMarker < numMarkers0) {
          macdata_i <- listOfMarkersForVarRatio[[k]][indexInMarkerList]
          i <- (MACdata$indexInGeno)[macdata_i]
          genoInd <- (MACdata$geno_ind)[macdata_i]
          cat(i, "th marker in geno ", genoInd, "\n")
          cat("MAC: ", (MACdata$MACvector)[macdata_i], "\n")
          if (genoInd == 0) {
            G0 <- Get_OneSNP_Geno(i - 1)
          } else if (genoInd == 1) {
            G0 <- Get_OneSNP_Geno_forVarRatio(i - 1)
          }

          if (sum(G0) / (2 * length(G0)) > 0.5) {
            G0 <- 2 - G0
          }
          G0sample <- G0
          print("length(G0)   aaaaa")
          print(length(G0))
          cat("G0", G0[1:10], "\n")
          print(dim(I_mat))
          print(length(G0sample))
          G0 <- as.numeric(I_mat %*% G0sample)
          cat("G0", G0[1:10], "\n")

          CHR <- bimPlink[Indexvector_forVarRatio[i] + 1, 1]
          cat("CHR ", CHR, "\n")
          print(bimPlink[Indexvector_forVarRatio[i] + 1, ])
          if (sum(G0) / (2 * length(G0)) > 0.5) {
            G0 <- 2 - G0
          }
          print("length(G0)")
          print(length(G0))
          NAset <- which(G0 == 0)
          AC <- sum(G0)
          print("length(NAset)")
          print(length(NAset))
          indexInMarkerList <- indexInMarkerList + 1
          if ((CHR >= 1 & CHR <= 22 & AC > 0 & AC < length(G0)) | includeNonautoMarkersforVarRatio) {
            AF <- AC / (2 * Nnomissing)
            if (CHR >= 1 & CHR <= 22) {
              autoMarker <- TRUE
            } else {
              autoMarker <- FALSE
            }

            G <- G0 - obj.noK$XXVX_inv %*% (obj.noK$XV %*% G0) # G1 is X adjusted


            set_isSparseGRM(useSparseGRMtoFitNULL)
            set_useGRMtoFitNULL(useGRMtoFitNULL)


            Sigma_iG <- getSigma_G_multiV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO = FALSE)
            Sigma_iX <- Sigma_iX_noLOCO

            var1 <- crossprod(G, Sigma_iG) - crossprod(G, Sigma_iX) %*% solve(crossprod(X, Sigma_iX)) %*% crossprod(X, Sigma_iG)
            cat("AC ", AC, "\n")
            S <- innerProduct(G, obj.glmm.null$residuals * var_weights)
            cat("S is ", S, "\n")
            cat("var1 is ", var1, "\n")
            p_exact <- pchisq(S^2 / var1, df = 1, lower.tail = F)
            cat("p_exact ", p_exact, "\n")



            if (!is.null(obj.glmm.null$eMat)) {
              var1GE_vec <- NULL
              var2sparseGE_vec <- NULL
              getildeMat <- NULL
              getilde_sample0_Mat <- NULL
              for (ne in 1:ncol(obj.glmm.null$eMat)) {
                evec <- obj.glmm.null$eMat[, ne]


                print("evec[1:100]")
                print(evec[1:100])
                GE <- G0 * evec
                print("length(GE)")
                print(length(GE))
                GE_tilde <- GE - obj.noK$XXVX_inv %*% (obj.noK$XV %*% GE)
                getildeMat <- cbind(getildeMat, GE_tilde)

                Sigma_iGE <- getSigma_G_multiV(W, tauVecNew, GE_tilde, maxiterPCG, tolPCG, LOCO = FALSE)
                var1GE <- crossprod(GE_tilde, Sigma_iGE) - crossprod(GE_tilde, Sigma_iX) %*% solve(crossprod(X, Sigma_iX)) %*% crossprod(X, Sigma_iGE)
                var1GE_vec <- c(var1GE_vec, var1GE)
                S_GE <- innerProduct(GE_tilde, obj.glmm.null$residuals * var_weights)
                p_exact_GE <- pchisq(S_GE^2 / var1GE, df = 1, lower.tail = F)
                cat("p_exact_GE ", p_exact_GE, "\n")
                cat("S_GE ", S_GE, "\n")
                cat("var1GE ", var1GE, "\n")

                if (useSparseGRMforVarRatio) {
                  set_isSparseGRM(useSparseGRMforVarRatio)
                  Sigma_iGE_sparse <- getSigma_G_noV(W, tauVecNew, GE_tilde, maxiterPCG, tolPCG, LOCO = FALSE)
                  var2_a_GE <- crossprod(GE_tilde, Sigma_iGE_sparse)
                  var2sparseGRM_GE <- var2_a_GE[1, 1]
                  var2sparseGE_vec <- c(var2sparseGE_vec, var2sparseGRM_GE)
                } else {
                  if (any(duplicated(obj.glmm.null$sampleID))) {
                    Sigma_iGE_sparse <- getSigma_G_V(W, tauVal, tauVecNew[1], GE_tilde, maxiterPCG, tolPCG)
                    var2_a_GE <- crossprod(GE_tilde, Sigma_iGE_sparse)
                    var2sparseGRM_GE <- var2_a_GE[1, 1]
                    var2sparseGE_vec <- c(var2sparseGE_vec, var2sparseGRM_GE)
                  } else {
                    var2sparseGE_vec <- c(var2sparseGE_vec, var1GE)
                  }
                }
              }
            }


            G_noXadj <- as.vector(G0sample - mean(G0sample))
            G0_sample_tilde <- G0sample - XXVXsample_inv0 %*% XVsample0 %*% G0sample

            if (any(duplicated(obj.glmm.null$sampleID))) {
              if (isStoreSigma) {
                Sigma_iG <- (obj.glmm.null$spSigma) %*% G0_sample_tilde
                var2_a <- crossprod(G0_sample_tilde, Sigma_iG)
              } else {
                G0_sample_tilde_I <- as.vector(I_mat %*% G0_sample_tilde)
                Sigma_iG <- getSigma_G_multiV(W, tauVecNew, G0_sample_tilde_I, maxiterPCG, tolPCG, LOCO = FALSE)
                var2_a <- crossprod(G0_sample_tilde_I, Sigma_iG)
              }
              var2sparseGRM <- var2_a[1, 1]
              cat("var2sparseGRM Here ", var2sparseGRM, "\n")
              varRatio_sparseGRM_vec <- c(varRatio_sparseGRM_vec, var1 / var2sparseGRM)
            } else {
              varRatio_sparseGRM_vec <- c(varRatio_sparseGRM_vec, 1)
            }


            if (obj.glmm.null$traitType == "binary") {
              stop("ERROR: This traitType is not supported in the current version.\n")
            } else if (obj.glmm.null$traitType == "quantitative") {
              stop("ERROR: This traitType is not supported in the current version.\n")
            } else if (obj.glmm.null$traitType == "count") {
              tmpW <- mu * var_weights
              tmpWI <- as.vector(crossprod(tmpW, I_mat))
              var2null <- innerProduct(tmpW, G * G)
              var2null_sample <- innerProduct(tmpWI, G0_sample_tilde * G0_sample_tilde)
              var2null_noXadj <- innerProduct(tmpWI, G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  var22nullGE <- innerProduct(tmpW, GE_tilde * GE_tilde)

                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "count_nb") {
              stop("ERROR: This traitType is not supported in the current version.\n")
            }

            cat("mu\n")
            print(mu[1:100])
            cat("AC ", AC, "\n")
            cat("var1 ", var1, "\n")
            cat("var2null ", var2null, "\n")
            cat("var2null_noXadj ", var2null_noXadj, "\n")
            varRatio_NULL_vec <- c(varRatio_NULL_vec, var1 / var2null)
            varRatio_NULL_sample_vec <- c(varRatio_NULL_sample_vec, var1 / var2null_sample)
            varRatio_NULL_noXadj_vec <- c(varRatio_NULL_noXadj_vec, var1 / var2null_noXadj)
            if (!is.null(obj.glmm.null$eMat)) {
              varRatio_NULL_eg_mat <- rbind(varRatio_NULL_eg_mat, var1GE_vec / var2nullGE_vec)
              varRatio_sparse_eg_mat <- rbind(varRatio_sparse_eg_mat, var1GE_vec / var2sparseGE_vec)
            }
            numTestedMarker <- numTestedMarker + 1
          } else {
            indexInMarkerList <- indexInMarkerList + 1
          }

          if (indexInMarkerList - 1 == length(listOfMarkersForVarRatio[[k]])) {
            numTestedMarker <- numMarkers0
          }
        } # end of while(numTestedMarker < numMarkers)


        print("varRatio_NULL_vec")
        print(varRatio_NULL_vec)
        print("varRatio_NULL_noXadj_vec")
        print(varRatio_NULL_noXadj_vec)

        ratioCV <- calCV(varRatio_NULL_noXadj_vec)

        if (ratioCV > ratioCVcutoff) {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
          numMarkers0 <- numMarkers0 + 10
          cat("try ", numMarkers0, " markers\n")
        } else {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
        }

        if (indexInMarkerList - 1 == length(listOfMarkersForVarRatio[[k]])) {
          ratioCV <- ratioCVcutoff
          cat("no more markers are available in the MAC category ", k, "\n")
          print(indexInMarkerList - 1)
        }
      } # end of while(ratioCV > ratioCVcutoff)

      if (length(varRatio_sparseGRM_vec) > 0) {
        cat("varRatio_sparseGRM_vec\n")
        print(varRatio_sparseGRM_vec)

        varRatio_sparse <- mean(varRatio_sparseGRM_vec)
        cat("varRatio_sparse", varRatio_sparse, "\n")
        varRatioTable <- rbind(varRatioTable, c(varRatio_sparse, "sparse", k))
      }
      varRatio_null <- mean(varRatio_NULL_vec)
      varRatio_null_sample <- mean(varRatio_NULL_sample_vec)
      cat("varRatio_null", varRatio_null, "\n")

      varRatio_null_noXadj <- mean(varRatio_NULL_noXadj_vec)
      cat("varRatio_null_noXadj", varRatio_null_noXadj, "\n")

      if (!is.null(obj.glmm.null$eMat)) {
        varRatio_NULL_eg_vec <- as.vector(colMeans(varRatio_NULL_eg_mat))
        varRatio_sparse_eg_vec <- as.vector(colMeans(varRatio_sparse_eg_mat))
        for (ne in 1:ncol(obj.glmm.null$eMat)) {
          varRatioTable <- rbind(varRatioTable, c(varRatio_NULL_eg_vec[ne], "null", 0))
          varRatioTable <- rbind(varRatioTable, c(varRatio_sparse_eg_vec[ne], "sparse", 0))
        }
        print(varRatio_NULL_eg_vec)
        print(varRatio_NULL_eg_mat)
      }

      varRatioTable <- rbind(varRatioTable, c(varRatio_null, "null", k))
      varRatioTable <- rbind(varRatioTable, c(varRatio_null_noXadj, "null_noXadj", k))
      varRatioTable <- rbind(varRatioTable, c(varRatio_null_sample, "null_sample", k))


    } else { # if(cateVarRatioVec[k] == 1)
      varRatioTable <- rbind(varRatioTable, c(1, "null", k))
      varRatioTable <- rbind(varRatioTable, c(1, "null_noXadj", k))
      varRatioTable <- rbind(varRatioTable, c(1, "null_sample", k))
      if (length(varRatio_sparseGRM_vec) > 0) {
        varRatioTable <- rbind(varRatioTable, c(1, "sparse", k))
      }
    }
  }
  write.table(varRatioTable, varRatioOutFile, quote = F, col.names = F, row.names = F)
  data <- read.table(varRatioOutFile, header = F)
  print(data)
}
