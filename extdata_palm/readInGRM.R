## Automatically generate GRM from:
## - existing prefix.gds (SNP or SeqArray), or
## - PLINK bed/bim/fam
## Only the final GRM is saved as RDS; temporary SNP GDS files are removed.

library(SeqArray)
library(SNPRelate)

base_dir <- "input"
prefix   <- file.path(base_dir, "stool_bialleic_merged_data")

## Possible input files
gds_fn   <- paste0(prefix, ".gds")   # original GDS (kept untouched)
bed_fn   <- paste0(prefix, ".bed")
bim_fn   <- paste0(prefix, ".bim")
fam_fn   <- paste0(prefix, ".fam")

## 1. Prepare a SNP GDS (temporary if we need to convert) -----------------------

snp_gds_use <- NULL   # path to the SNP GDS we will actually open
tmp_created <- FALSE  # whether snp_gds_use is a temporary file we should delete

if (file.exists(gds_fn)) {
  message("Detected .gds file: ", gds_fn)
  
  ## Try to open as SNP GDS first (without modifying it)
  is_snp_gds <- FALSE
  gds_try <- try(snpgdsOpen(gds_fn), silent = TRUE)
  
  if (!inherits(gds_try, "try-error")) {
    ## Successfully opened as SNP GDS
    is_snp_gds <- TRUE
    snpgdsClose(gds_try)
  }
  
  if (is_snp_gds) {
    message("File is SNP GDS (SNP_ARRAY). Using it directly: ", gds_fn)
    snp_gds_use <- gds_fn
    tmp_created <- FALSE
  } else {
    message("File is not SNP GDS. Assuming SeqArray GDS and converting to a temporary SNP GDS...")
    
    ## Convert SeqArray GDS -> temporary SNP GDS
    tmp_snp_gds <- paste0(prefix, "_tmp_snp.gds")
    
    seqfile_try <- try(seqOpen(gds_fn), silent = TRUE)
    if (inherits(seqfile_try, "try-error")) {
      stop("Failed to open ", gds_fn, " as SeqArray GDS. ",
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
    
    message("SeqArray -> temporary SNP GDS conversion finished: ", tmp_snp_gds)
    snp_gds_use <- tmp_snp_gds
    tmp_created <- TRUE
  }
  
} else if (all(file.exists(bed_fn, bim_fn, fam_fn))) {
  message("No .gds found. Detected PLINK bed/bim/fam: ", bed_fn)
  
  ## Create a temporary SNP GDS from PLINK
  tmp_snp_gds <- paste0(prefix, "_tmp_snp.gds")
  
  snpgdsBED2GDS(
    bed.fn      = bed_fn,
    fam.fn      = fam_fn,
    bim.fn      = bim_fn,
    out.gdsfn   = tmp_snp_gds,
    snpfirstdim = TRUE
  )
  
  message("PLINK -> temporary SNP GDS conversion finished: ", tmp_snp_gds)
  snp_gds_use <- tmp_snp_gds
  tmp_created <- TRUE
  
} else {
  stop("No usable genotype input found:\n",
       "  - .gds: ", gds_fn, "\n",
       "  - or PLINK: ", bed_fn, " / .bim / .fam\n")
}

## 2. Compute GRM from SNP GDS -----------------------------------------------

message("Computing GRM from SNP GDS: ", snp_gds_use)

genofile <- snpgdsOpen(snp_gds_use)

grm_res <- snpgdsGRM(
  genofile,
  method     = "GCTA",
  num.thread = 4
)

K         <- grm_res$grm       # GRM matrix
sample.id <- grm_res$sample.id # sample order

snpgdsClose(genofile)

## If we created a temporary SNP GDS, remove it now
if (tmp_created && file.exists(snp_gds_use)) {
  message("Removing temporary SNP GDS: ", snp_gds_use)
  file.remove(snp_gds_use)
}

## 3. Save GRM as RDS using prefix --------------------------------------------

K_rds_path   <- paste0(prefix, "_GRM_K.rds")
sid_rds_path <- paste0(prefix, "_GRM_sampleid.rds")

saveRDS(K,         file = K_rds_path)
saveRDS(sample.id, file = sid_rds_path)

message("GRM saved to:")
message("  - ", K_rds_path)
message("  - ", sid_rds_path)
