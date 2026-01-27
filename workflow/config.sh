
################################# parameter settings below #################################
##################################### required below #######################################
# PALM method; 1 for 'PALM' or 2 for 'PALM-mbQTL'
PALMmethod=2 

# input and output folders
inputFolder=input
outputFolder=output

# input files
genoFile=${inputFolder}/geno
abdFile=${inputFolder}/abd.txt
covFile=${inputFolder}/cov.txt

# chromosome number to analyze (for step2 only); set to 'NULL' for all chromosomes
chrom=NULL

##################################### optional below #######################################
# Default paths for other files, can be modified as needed
grmFile=${inputFolder}/grm.rds
outputSnpFile=${outputFolder}/info_snp.txt
outputFeatureFile=${outputFolder}/info_feature.txt
palm1_step1_prefix=${outputFolder}/palm1_step1_allpheno
palm2_step1_prefix=${outputFolder}/palm2_step1_allpheno
if [[ "${chrom}" == "NULL" ]]; then
  palm1_step2_prefix=${outputFolder}/palm1_step2_allchr
  palm2_step2_prefix=${outputFolder}/palm2_step2_allchr
else
  palm1_step2_prefix=${outputFolder}/palm1_step2_chr${chrom}
  palm2_step2_prefix=${outputFolder}/palm2_step2_chr${chrom}
fi




###################### export variables for scripts (do not modify) ########################
export PALMmethod
export inputFolder
export outputFolder
export genoFile
export abdFile
export covFile
export chrom
export grmFile
export outputSnpFile
export outputFeatureFile
export palm1_step1_prefix
export palm1_step2_prefix
export palm2_step1_prefix
export palm2_step2_prefix