################################# parameter settings below #################################
######### required below #########
PALMmethod=2 # 1 for 'PALM' or 2 for 'PALM-mbQTL'

inputFolder=input
outputFolder=output

genoFile=./${inputFolder}/geno
chrom=1

######### optional below #########
palm1_step1_prefix=./${outputFolder}/palm1_step1_allpheno
palm1_step2_prefix=./${outputFolder}/palm1_step2${chrom:+_chr${chrom}}
palm2_step1_prefix=./${outputFolder}/palm2_step1_allpheno
palm2_step2_prefix=./${outputFolder}/palm2_step2${chrom:+_chr${chrom}}



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

# step2: score test for phenoCol
echo "Performing score test..."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_palm.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm1_step1_prefix}.rda \
        --PALMOutputFile=${palm1_step2_prefix} \
        --chrom=${chrom} \
        --correct=NULL \
        --cluster=NULL
else
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_scoreTest.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm2_step1_prefix}.rda \
        --PALMOutputFile=${palm2_step2_prefix} \
        --chrom=${chrom} \
        --phenoCol=${phenoCol} \
        --minMAF=0
fi