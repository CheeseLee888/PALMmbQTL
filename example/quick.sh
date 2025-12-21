################################# parameter settings below #################################
######### PALMmethod option ##########
PALMmethod=1 # 1 for 'PALM' or 2 for 'PALM-mbQTL'

######### required below #########
inputFolder=input
outputFolder=output

genoFile=./${inputFolder}/geno
abdFile=./${inputFolder}/abd.txt
covFile=./${inputFolder}/cov.txt
sampleIDColinabdFile=IID
sampleIDColincovFile=IID
phenoCol=g_1
chrom=1
# offsetCol=SeqDepth
covarColList=AGE,SEX

######### optional below #########
grmFile=./${inputFolder}/grm.rds
step1_prefix=./${outputFolder}/allpheno_step1
step2_prefix=./${outputFolder}/${phenoCol}_step2${chrom:+_chr${chrom}}
step1_prefix_palm=./${outputFolder}/allpheno_step1_palm
# step2_prefix_palm





################################# workflow below (do not modify) #################################
if [[ "${PALMmethod}" != 1 && "${PALMmethod}" != 2 ]]; then
    echo "PALMmethod must be specified as 1 (PALM) or 2 (PALM-mbQTL)."
    exit 1
fi

if [[ "${PALMmethod}" == 1 ]]; then
    echo "Running PALM method..."
else
    echo "Running PALM-mbQTL method..."
fi

mkdir -p "${outputFolder}"

# step0: generate GRM from genotype data
if [[ "${PALMmethod}" == 2 ]]; then
    echo "Generating GRM from genotype data..."
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_generateGRM.R \
        --genoFile=${genoFile} \
        --grmFile=${grmFile}
fi

# step1: fit null model for all phenotypes
echo "Fitting null model for all phenotypes..."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${step1_prefix_palm}
else
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_fitNULL.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --grmFile=${grmFile} \
        --covarColList=${covarColList} \
        --outputPrefix=${step1_prefix} \
        --useGRMtoFitNULL=TRUE \
        --sampleIDColinabdFile=${sampleIDColinabdFile} \
        --sampleIDColincovFile=${sampleIDColincovFile}
fi

# step2: score test for phenoCol
echo "Performing score test..."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_palm.R \
        --inFile=${genoFile} \
        --correct=NULL \
        --NULLmodelFile=${step1_prefix_palm}.rda \
        --PALMOutputFile=${outputFolder} \
        --chrom=${chrom}
else
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_scoreTest.R \
        --inFile=${genoFile} \
        --phenoCol=${phenoCol} \
        --chrom=${chrom} \
        --PALMOutputFile=${step2_prefix}.txt \
        --NULLmodelFile=${step1_prefix}.rda \
        --minMAF=0
fi