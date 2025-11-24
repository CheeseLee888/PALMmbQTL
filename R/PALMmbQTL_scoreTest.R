#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
#'
#' @param bgenFile character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported
#' @param bgenFileIndex character. Path to the .bgi file (index of the bgen file)
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the bgen file. The file does not contain header lines.
#' @param vcfFile character. Path to vcf file
#' @param vcfFileIndex character. Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .csi file using 'tabix --csi -p vcf file.vcf.gz'
#' @param vcfField character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".
#' @param savFile character. Path to sav file
#' @param savFileIndex character. Path to index for sav file .s1r
#' @param bedFile character. Path to bed file (PLINK)
#' @param bimFile character. Path to bim file (PLINK)
#' @param famFile character. Path to fam file (PLINK)
#' @param AlleleOrder character. alt-first or ref-first for bgen or PLINK files. By default, alt-first
#' @param idstoIncludeFile character. Path to a file containing variant ids to be included from the dosage file. The file does not have a header and each line is for a marker ID. Variant ids are in the format chr:pos_ref/alt
#' @param rangestoIncludeFile character. Path to a file containing genome regions to be included from the dosage file. The file contains three columns for chromosome, start, and end respectively with no header. Note for vcf and sav files, only the first line in the file will be used.
#' @param chrom character. If LOCO is specified, chrom is required. chrom is also required for VCF/BCF/SAV input. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, 1 does not match chr1.
#' @param is_imputed_data logical. Whether the dosages/genotypes imputed are imputed. If TRUE, the program will output the imputed info score. By default, FALSE.
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 0.5. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param minInfo numeric. Minimum imputation info of markers to test. By default, 0.
#' @param maxMissing numeric. Maximum missing rate for markers to be tested. By default, 0.15
#' @param impute_method character. Imputation method for missing dosages. best_guess, mean or minor. best_guess: missing dosages imputed as best guessed genotyes round(2*allele frequency). mean: missing dosages are imputed as mean (2*allele frequency). minor: missing dosages are imputed as minor allele homozygotes. By default, minor
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. If TRUE, --chrom is required. By default, TRUE
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param GMMATmodel_varianceRatio_multiTraits_File character. Path to the input file containing 3 columns: phenotype name, model file, and variance ratio file. Each line is for one phenotype. This file is used when multiple phenotypes are analyzed simutaneously
#' @param SAIGEOutputFile character. Prefix of the output files containing assoc test results
#' @param markers_per_chunk character. Number of markers to be tested and output in each chunk in the single-variant assoc tests. By default, 10000
#' @param groups_per_chunk character. Number of groups/sets to be read in and tested in each chunk in the set-based assoc tests. By default, 100
#' @param is_output_moreDetails logical. Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns homN_Allele2_cases, hetN_Allelelogical2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls will be output. By default, FALSE
#' @param is_overwrite_output logical. Whether to overwrite the output file if it exists. If FALSE, the program will continue the unfinished analysis instead of starting over from the beginining. By default, TRUE
#' @param maxMAF_in_groupTest. vector of numeric. Max MAF for markers tested in group test seperated by comma. e.g. c(0.0001,0.001,0.01). By default, c(0.01)
#' @param maxMAC_in_groupTest. vector of numeric. Max MAC for markers tested in group test seperated by comma. This vector will be combined with maxMAF_in_groupTest. e.g. c(1) to only test singletons. By default, c(0) and no Max MAC cutoffs are applied.
#' @param minGroupMAC_in_BurdenTest numeric. Only applied when only Burden tests are performed (r.corr=1). Minimum minor allele count in the Burden test for the psueodo marker. By default, 5
#' @param  annotation_in_groupTest. vector of character. annotations of markers to be tested in the set-based tests. using ; to combine multiple annotations in the same test. e.g. c("lof","missense;lof","missense;lof;synonymous")  will test lof variants only, missense+lof variants, and missense+lof+synonymous variants. By default:  c("lof","missense;lof","missense;lof;synonymous")
#' @param groupFile character. Path to the file containing the group information for gene-based tests. Each gene/set has 2 or 3 lines in the group file. The first element is the gene/set name. The second element in the first line is to indicate whether this line contains variant IDs (var), annotations (anno), or weights (weight). The line for weights is optional. If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use weights.beta to change the parameters for the Beta distribution. The variant ids must be in the format chr:pos_ref/alt. Elements are seperated by tab or space.
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file that was used in Step 1
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM. By default, 0
#' @param MACCutoff_to_CollapseUltraRare numeric. MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10.
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param weights.beta vector of numeric with two elements. parameters for the beta distribution to weight genetic markers in gene-based tests. By default, "c(1,25)".
#' @param r.corr numeric. bewteen 0 and 1. parameters for gene-based tests. If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O tests will be performed and results for Burden tests and SKAT tests will be output too.  By default, 0.
#' @param markers_per_chunk_in_groupTest numeric. Number of markers in each chunk when calculating the variance covariance matrix in the set/group-based tests. By default, 100.
#' @param condition character. For conditional analysis. Variant ids are in the format chr:pos_ref/alt and seperated by by comma. e.g."chr3:101651171:C:T,chr3:101651186:G:A".
#' @param weights_for_condition. matrix of numeric. weights for conditioning markers for gene- or region-based tests. The nrow equals to the number of conditioning markers and the ncol equals to the number of sets of weights specified in the group file, e.g. c(1,2,3). If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use weights.beta to change the parameters for the Beta distribution.
#' @param SPAcutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param dosage_zerod_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 0.2
#' @param dosage_zerod_MAC_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 10
#' @param is_single_in_groupTest logical.  Whether to output single-variant assoc test results when perform group tests. Note, single-variant assoc test results will always be output when SKAT and SKAT-O tests are conducted with r.corr=0. This parameter should only be used when only Burden tests are condcuted with r.corr=1. By default, TRUE
#' @param is_equal_weight_in_groupTest logical. Whether equal weights are used in group Test. If TRUE, equal weights will be included in the group test. By default, FALSE
#' @param is_output_markerList_in_groupTest logical. Whether to output the marker lists included in the set-based tests for each mask. By default, FALSE
#' @param is_Firth_beta logical. Whether to estimate effect sizes using approx Firth, only for binary traits. By default, FALSE
#' @param pCutoffforFirth numeric. p-value cutoff to use approx Firth to estiamte the effect sizes. Only for binary traits. The effect sizes of markers with p-value <= pCutoffforFirth will be estimated using approx Firth. By default, 0.01.
#' @param IsOutputlogPforSingle logical. Whether to output log(Pvalue) for single-variant assoc tests. By default, FALSE. If TRUE, the log(Pvalue) instead of original P values will be output (Not activated)
#' @param X_PARregion character. ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes. In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be multiplied by 2 (Not activated).
#' @param is_rewrite_XnonPAR_forMales logical. Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Note, only use is_rewrite_XnonPAR_forMales=TRUE when the specified VCF or Bgen file only has variants on chromosome X. When is_rewrite_XnonPAR_forMales=TRUE, the program does not check the chromosome value by assuming all variants are on chromosome X (Not activated)
#' @param sampleFile_male character. Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter
#' @return SAIGEOutputFile
#' @export
SPAGMMATtest <- function(inFile = "",
                         AlleleOrder = "alt-first", # new
                         idstoIncludeFile = "",
                         rangestoIncludeFile = "",
                         chrom = "", # for vcf file
                         max_missing = 0.15, # new
                         impute_method = "best_guess", # "mean", "minor", "best_guess"     #new
                         min_MAC = 0.5,
                         min_MAF = 0,
                         min_Info = 0,
                         is_imputed_data = FALSE, # new
                         GMMATmodelFile = "",
                         LOCO = TRUE,
                         varianceRatioFile = "",
                         GMMATmodel_varianceRatio_multiTraits_File = "",
                         cateVarRatioMinMACVecExclude = c(10.5, 20.5),
                         cateVarRatioMaxMACVecInclude = c(20.5),
                         SPAcutoff = 2,
                         SAIGEOutputFile = "",
                         markers_per_chunk = 10000,
                         groups_per_chunk = 100,
                         markers_per_chunk_in_groupTest = 100, # new
                         condition = "",
                         sparseGRMFile = "",
                         sparseGRMSampleIDFile = "",
                         VmatFilelist = "",
                         VmatSampleFilelist = "",
                         relatednessCutoff = 0,
                         groupFile = "",
                         weights.beta = c("1,25"),
                         weights_for_condition = NULL,
                         r.corr = 0,
                         dosage_zerod_cutoff = 0.2,
                         dosage_zerod_MAC_cutoff = 10,
                         is_output_moreDetails = FALSE, # new
                         MACCutoff_to_CollapseUltraRare = 10,
                         annotation_in_groupTest = c("lof", "missense;lof", "missense;lof;synonymous"), # new
                         maxMAF_in_groupTest = c(0.01),
                         minMAF_in_groupTest_Exclude = NULL,
                         maxMAC_in_groupTest = c(0),
                         minMAC_in_groupTest_Exclude = NULL,
                         minGroupMAC_in_BurdenTest = 5,
                         is_Firth_beta = FALSE,
                         pCutoffforFirth = 0.01,
                         is_overwrite_output = TRUE,
                         is_single_in_groupTest = TRUE,
                         is_SKATO = FALSE,
                         is_equal_weight_in_groupTest = FALSE,
                         is_output_markerList_in_groupTest = FALSE,
                         pval_cutoff_for_fastTest = 0.05,
                         is_fastTest = FALSE,
                         pval_cutoff_for_gxe = 0.001,
                         is_noadjCov = TRUE,
                         is_sparseGRM = TRUE,
                         max_MAC_use_ER = 4,
                         is_EmpSPA = FALSE) {
  if (!(impute_method %in% c("best_guess", "mean", "minor"))) {
    stop("impute_method should be 'best_guess', 'mean' or 'minor'.")
  }

  checkArgsListBool(
    is_imputed_data = is_imputed_data,
    LOCO = LOCO,
    is_output_moreDetails = is_output_moreDetails,
    is_overwrite_output = is_overwrite_output
  )
  cat("dosage_zerod_cutoff ", dosage_zerod_cutoff, "\n")
  checkArgsListNumeric(
    start = 1,
    end = 250000000,
    max_missing = max_missing,
    min_MAC = min_MAC,
    min_MAF = min_MAF,
    min_Info = min_Info,
    SPAcutoff = SPAcutoff,
    dosage_zerod_cutoff = dosage_zerod_cutoff,
    dosage_zerod_MAC_cutoff = dosage_zerod_MAC_cutoff,
    markers_per_chunk = markers_per_chunk,
    groups_per_chunk = groups_per_chunk,
    minGroupMAC_in_BurdenTest = minGroupMAC_in_BurdenTest,
    max_MAC_use_ER = max_MAC_use_ER
  )




  ## check and create the output file
  OutputFile <- SAIGEOutputFile
  OutputFileIndex <- NULL
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }


  if (!is.null(minMAF_in_groupTest_Exclude)) {
    min_MAF <- min(min_MAF, min(minMAF_in_groupTest_Exclude))
  }
  if (!is.null(minMAC_in_groupTest_Exclude)) {
    min_MAC <- min(min_MAF, min(minMAC_in_groupTest_Exclude))
  }



  if (groupFile == "") {
    isGroupTest <- FALSE
    cat("single-variant association test will be performed\n")
  } else {
    isGroupTest <- TRUE
    Check_File_Exist(groupFile, "groupFile")
    cat("group-based test will be performed\n")

    cat("maxMAF_in_groupTest ", maxMAF_in_groupTest, "\n")
    cat("minMAF_in_groupTest_Exclude ", minMAF_in_groupTest_Exclude, "\n")


    checkArgsList_for_Region(
      MACCutoff_to_CollapseUltraRare,
      # DosageCutoff_for_UltraRarePresence,
      maxMAF_in_groupTest = maxMAF_in_groupTest,
      minMAF_in_groupTest_Exclude = minMAF_in_groupTest_Exclude,
      maxMAC_in_groupTest = maxMAC_in_groupTest,
      minMAC_in_groupTest_Exclude = minMAC_in_groupTest_Exclude,
      markers_per_chunk_in_groupTest = markers_per_chunk_in_groupTest
    )




    IsOutputlogPforSingle <- FALSE # to check
  }

  # if (GMMATmodel_varianceRatio_multiTraits_File != "") {
  #   if (GMMATmodelFile != "" | varianceRatioFile != "") {
  #     stop("GMMATmodel_varianceRatio_multiTraits_File is specified while varianceRatioFile and/or GMMATmodelFile are also specified. Please check\n")
  #   } else {
  #     modelvrfile_data <- data.table::fread(GMMATmodel_varianceRatio_multiTraits_File, header = F, data.table = F)
  #     if (ncol(modelvrfile_data) != 3) {
  #       stop("GMMATmodel_varianceRatio_multiTraits_File needs to have 3 columns: phenotype name, model file, and variance ratio file\n")
  #     } else {
  #       GMMATmodelFile_vec <- modelvrfile_data[, 2]
  #       GMMATmodelFile <- paste(GMMATmodelFile_vec, collapse = ",")
  #       varianceRatioFile_vec <- modelvrfile_data[, 3]
  #       varianceRatioFile <- paste(varianceRatioFile_vec, collapse = ",")
  #       phenotype_name_vec <- as.character(modelvrfile_data[, 1])
  #     }
  #   }
  # } else {
  #   GMMATmodelFile_vec <- unlist(strsplit(GMMATmodelFile, split = ","))
  #   phenotype_name_vec <- as.character(seq(1, length(GMMATmodelFile_vec)))
  # }

  load(GMMATmodelFile)

  if (!LOCO) {
    print("LOCO = FASLE and leave-one-chromosome-out is not applied")
  }
  
  isSparseGRM <- is_sparseGRM
  cat("isSparseGRM ", isSparseGRM, "\n")


  if (!isGroupTest) {
    OutputFile <- SAIGEOutputFile
    
    cat("Starting glmmscore test...\n")
    GMMAT::glmm.score(
      obj = modglmm,
      infile = inFile,
      center = T, 
      outfile = OutputFile
      )
  } else {
    stop("Group-based test is not yet implemented in PALM-mbQTL.\n")
  }
}
