################################# parameter settings below #################################
######### required below #########
PALMmethod=2 # 1 for 'PALM' or 2 for 'PALM-mbQTL'

inputFolder=input
outputFolder=output

genoFile=./${inputFolder}/geno
abdFile=./${inputFolder}/abd.txt
covFile=./${inputFolder}/cov.txt

######### optional below #########
grmFile=./${inputFolder}/grm.rds
palm1_step1_prefix=./${outputFolder}/palm1_step1_allpheno
palm2_step1_prefix=./${outputFolder}/palm2_step1_allpheno



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

# step0: check input files
echo "Checking input files..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_checkInput.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --genoFile=${genoFile}


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