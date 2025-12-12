################################# parameter settings below #################################
######### required below #########
inputFolder=input
outputFolder=output

genoFile=./${inputFolder}/geno
abdFile=./${inputFolder}/abd.txt
covFile=./${inputFolder}/cov.txt
sampleIDColinabdFile=IID
sampleIDColincovFile=IID
phenoCol=g_1

# offsetCol=SeqDepth
covarColList=AGE,SEX

######### optional below #########
mergeOutFile=./${inputFolder}/merged.txt
grmFile=./${inputFolder}/grm.rds
step1_prefix=./${outputFolder}/allpheno_step1
step2_prefix=./${outputFolder}/${phenoCol}_step2



################################# workflow below (do not modify) #################################
# mkdir -p "${outputFolder}"
# # step0: generate GRM from genotype data
# echo "Generating GRM from genotype data..."
# pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_generateGRM.R \
#     --genoFile=${genoFile} \
#     --grmFile=${grmFile}

# # step1: fit null model for all phenotypes
# echo "Fitting null model for all phenotypes..."
# pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_fitNULL.R \
#     --abdFile=${abdFile} \
#     --covFile=${covFile} \
#     --sampleIDColincovFile=${sampleIDColincovFile} \
#     --grmFile=${grmFile} \
#     --covarColList=${covarColList} \
#     --outputPrefix=${step1_prefix} \
#     --useGRMtoFitNULL=TRUE \
#     --sampleIDColinabdFile=${sampleIDColinabdFile} \
#     --sampleIDColincovFile=${sampleIDColincovFile} \
#     --traitType=count


# step2: score test for phenoCol
echo "Performing score test for ${phenoCol}..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_scoreTest.R \
    --inFile=${genoFile} \
    --phenoCol=${phenoCol} \
    --chrom=${chrom} \
    --PALMOutputFile=${step2_prefix}.txt \
    --NULLmodelFile=${step1_prefix}.rda \
    --minMAF=0
