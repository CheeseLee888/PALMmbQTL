#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
#'
#' @param inFile character. Path to the input file containing the genotype data in PLINK binary format (bed/bim/fam) or VCF format.
#' @param NULLmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param PALMOutputFile character. Prefix of the output files containing assoc test results
#' @return PALMOutputFile
#' @export
SPAGMMATtest <- function(inFile = "",
                         NULLmodelFile = "",
                         PALMOutputFile = "") {

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
  
  if (!isGroupTest) {
    cat("Starting glmmscore test...\n")
    GMMAT::glmm.score(
      obj = modglmm,
      infile = inFile,
      center = T, 
      outfile = PALMOutputFile
      )
    cat("glmmscore test finished.\n")
  } else {
    stop("Group-based test is not yet implemented in PALM-mbQTL.\n")
  }
}
