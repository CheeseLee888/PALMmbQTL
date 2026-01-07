################################# parameter settings below #################################
######### required below #########
PALMmethod=1 # 1 for 'PALM' or 2 for 'PALM-mbQTL'

inputFolder=input
outputFolder=output

genoFile=./${inputFolder}/geno
abdFile=./${inputFolder}/abd.txt
covFile=./${inputFolder}/cov.txt
chrom=1

######### optional below #########
grmFile=./${inputFolder}/grm.rds
palm1_step1_prefix=./${outputFolder}/palm1_step1_allpheno
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

mkdir -p "${outputFolder}"

# step0: generate GRM from genotype data (only for PALM-mbQTL)
if [[ "${PALMmethod}" == 2 ]]; then
    if [[ -f "${grmFile}" ]]; then
        echo "GRM already exists at: ${grmFile}"
        echo "Skip generating GRM and reuse the existing GRM."
    else
        echo "Generating GRM from genotype data..."
        pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_generateGRM.R \
            --genoFile=${genoFile} \
            --grmFile=${grmFile}
    fi
fi

# step1: fit null model for all phenotypes
echo "Fitting null model for all phenotypes..."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${palm1_step1_prefix}
else
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_fitNULL.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --grmFile=${grmFile} \
        --outputPrefix=${palm2_step1_prefix} \
        --covarColList=all \
        --useGRMtoFitNULL=TRUE
fi

# step2: score test for phenoCol
echo "Performing score test..."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_palm.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm1_step1_prefix}.rda \
        --PALMOutputFile=${outputFolder} \
        --chrom=${chrom} \
        --correct=NULL \
        --cluster=
else
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_scoreTest.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm2_step1_prefix}.rda \
        --PALMOutputFile=${palm2_step2_prefix} \
        --chrom=${chrom} \
        --phenoCol=${phenoCol} \
        --minMAF=0
fi

# Step 3: Generate feature information
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step3_info.R \
    --abdFile=${abdFile} \
    --genoFile=${genoFile} \
    --outputFile=${outputFolder}/feature_info.txt