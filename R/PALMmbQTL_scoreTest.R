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

  isGroupTest <- FALSE
  cat("single-variant association test will be performed\n")

  # if (GMMATmodel_varianceRatio_multiTraits_File != "") {
  #   if (NULLmodelFile != "" | varianceRatioFile != "") {
  #     stop("GMMATmodel_varianceRatio_multiTraits_File is specified while varianceRatioFile and/or NULLmodelFile are also specified. Please check\n")
  #   } else {
  #     modelvrfile_data <- data.table::fread(GMMATmodel_varianceRatio_multiTraits_File, header = F, data.table = F)
  #     if (ncol(modelvrfile_data) != 3) {
  #       stop("GMMATmodel_varianceRatio_multiTraits_File needs to have 3 columns: phenotype name, model file, and variance ratio file\n")
  #     } else {
  #       GMMATmodelFile_vec <- modelvrfile_data[, 2]
  #       NULLmodelFile <- paste(GMMATmodelFile_vec, collapse = ",")
  #       varianceRatioFile_vec <- modelvrfile_data[, 3]
  #       varianceRatioFile <- paste(varianceRatioFile_vec, collapse = ",")
  #       phenotype_name_vec <- as.character(modelvrfile_data[, 1])
  #     }
  #   }
  # } else {
  #   GMMATmodelFile_vec <- unlist(strsplit(NULLmodelFile, split = ","))
  #   phenotype_name_vec <- as.character(seq(1, length(GMMATmodelFile_vec)))
  # }

  load(NULLmodelFile)

  ## extract all pheno_name
  all_pheno_names <- vapply(
    null_list,
    function(x) x$pheno_name,
    FUN.VALUE = character(1L)
  )

  ## use which to find the corresponding index
  idx <- which(all_pheno_names == phenoCol)

  if (length(idx) == 0) {
    stop(sprintf("Phenotype '%s' not found in NULLmodelFile.", phenoCol))
  }else{
    cat(sprintf("Phenotype '%s' found in NULLmodelFile at index %d.\n", phenoCol, idx))
  }

  ## extract the corresponding model based on idx
  modglmm <- null_list[[idx]]$modglmm


  if (!isGroupTest) {
    cat("Starting glmmscore test...\n")
    GMMAT::glmm.score(
      obj = modglmm,
      infile = inFile,
      center = T, 
      outfile = PALMOutputFile,
      MAF.range = c(minMAF, 0.5)
      )
    cat("glmmscore test for", phenoCol, "finished.\n")
  } else {
    stop("Group-based test is not yet implemented in PALM-mbQTL.\n")
  }
}
