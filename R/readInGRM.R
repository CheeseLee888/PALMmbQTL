## Automatically generate GRM from:
## - existing prefix.gds (SNP or SeqArray), or
## - PLINK bed/bim/fam
## Only the final GRM is saved as RDS; temporary SNP GDS files are removed.

generateGRM <- function(genoFile = "", grmFile = "") {

  if (!nzchar(genoFile)) stop("Please provide 'genoFile' (PLINK prefix or .gds).")
  if (!nzchar(grmFile))  stop("Please provide 'grmFile' (output RDS path).")

  ## 1. Prepare a SNP GDS (temporary if we need to convert) -----------------------

  snp_gds_use <- NULL   # path to the SNP GDS we will actually open
  tmp_created <- FALSE  # whether snp_gds_use is a temporary file we should delete
  
  # genoFile is the full path
  # genoFile = dir + "/" + prefix + "." + ext
  # Get directory of genoFile
  dir <- dirname(genoFile)
  # Get file prefix
  prefix <- tools::file_path_sans_ext(basename(genoFile))
  # Determine genotype file type based on extension
  ext <- tools::file_ext(basename(genoFile))

  if (ext == "") { # use bim/bed/fam
    bed_fn   <- paste0(genoFile, ".bed")
    bim_fn   <- paste0(genoFile, ".bim")
    fam_fn   <- paste0(genoFile, ".fam")
    if (!all(file.exists(bed_fn, bim_fn, fam_fn))) {
      stop("PLINK files not found for prefix '", genoFile, "': need .bed, .bim and .fam.")
    }
    message("Detected bed/bim/fam file.")
    
    ## Create a temporary SNP GDS from PLINK
    tmp_snp_gds <- paste0(genoFile, "_tmp_snp.gds")
    
    snpgdsBED2GDS(
      bed.fn      = bed_fn,
      fam.fn      = fam_fn,
      bim.fn      = bim_fn,
      out.gdsfn   = tmp_snp_gds,
      snpfirstdim = TRUE
    )
    
    message("PLINK -> temporary SNP GDS conversion finished.")
    snp_gds_use <- tmp_snp_gds
    tmp_created <- TRUE

  } else if (tolower(ext) == "gds"){ # use gds
    message("Detected gds file.")
    
    ## Try to open as SNP GDS first (without modifying it)
    is_snp_gds <- FALSE
    gds_try <- try(snpgdsOpen(genoFile), silent = TRUE)
    
    if (!inherits(gds_try, "try-error")) {
      ## Successfully opened as SNP GDS
      is_snp_gds <- TRUE
      snpgdsClose(gds_try)
    }
    
    if (is_snp_gds) {
      message("File is SNP GDS (SNP_ARRAY). Using it directly.")
      snp_gds_use <- genoFile
      tmp_created <- FALSE
    } else {
      message("File is not SNP GDS. Assuming SeqArray GDS and converting to a temporary SNP GDS ...")
      
      ## Convert SeqArray GDS -> temporary SNP GDS
      tmp_snp_gds <- file.path(dir, paste0(prefix, "_tmp_snp.gds"))
      
      seqfile_try <- try(seqOpen(genoFile), silent = TRUE)
      if (inherits(seqfile_try, "try-error")) {
        stop("Failed to open ", genoFile, " as SeqArray GDS. ",
            "It is neither SNP GDS nor SeqArray GDS.")
      }
      seqfile <- seqfile_try
      
      # If needed, add seqSetFilter(seqfile, ...) here for filtering/subsetting
      
      seqGDS2SNP(
        seqfile,
        out.gdsfn           = tmp_snp_gds,
        compress.geno       = "ZIP_RA",
        compress.annotation = "ZIP_RA",
        optimize            = TRUE
      )
      
      seqClose(seqfile)
      
      message("SeqArray -> temporary SNP GDS conversion finished.")
      snp_gds_use <- tmp_snp_gds
      tmp_created <- TRUE
    }
  } else {
    stop("Unsupported file extension for genoFile. "
        "\n  - Supported: no extension (PLINK bed/bim/fam) or .gds")
  }

  ## 2. Compute GRM from SNP GDS -----------------------------------------------

  message("Computing GRM from SNP GDS ... ")

  geno <- snpgdsOpen(snp_gds_use)

  grm_res <- snpgdsGRM(
    geno,
    method     = "GCTA",
    num.thread = 4
  )

  K         <- grm_res$grm       # GRM matrix
  sample.id <- grm_res$sample.id # sample order

  snpgdsClose(geno)
  ## If we created a temporary SNP GDS, remove it now
  if (tmp_created && file.exists(snp_gds_use)) {
    message("Removing temporary SNP GDS ... ")
    file.remove(snp_gds_use)
  }

  ## 3. Save GRM as RDS --------------------------------------------
  saveRDS(list(K = K, sample.id = sample.id), file = grmFile)

  message("GRM saved to:")
  message("  - ", grmFile)
}
