################################# parameter settings below #################################
######### required below #########
PALMmethod=1 # 1 for 'PALM' or 2 for 'PALM-mbQTL'

inputFolder=input
outputFolder=output

genoFile=geno
abdFile=abd.txt
covFile=cov.txt
sampleIDColinabdFile=IID
sampleIDColincovFile=IID
chrom=1

######### optional below #########
grmFile=grm.rds
palm1_step1_prefix=palm1_step1_allpheno
palm2_step1_prefix=palm2_step1_allpheno
palm2_step2_prefix=palm2_step2${chrom:+_chr${chrom}}





###################### change to absolute paths in docker (do not modify) ########################
workDir="$(pwd)"
inputFolder="${workDir}/${inputFolder}"
outputFolder="${workDir}/${outputFolder}"
genoFile="${inputFolder}/${genoFile}"
abdFile="${inputFolder}/${abdFile}"
covFile="${inputFolder}/${covFile}"
grmFile="${inputFolder}/${grmFile}"
palm1_step1_prefix="${outputFolder}/${palm1_step1_prefix}"
palm2_step1_prefix="${outputFolder}/${palm2_step1_prefix}"
palm2_step2_prefix="${outputFolder}/${palm2_step2_prefix}"

echo "[INFO] Using input folder: ${inputFolder}"
echo "[INFO] Using output folder: ${outputFolder}"
echo "[INFO] Using genotype file: ${genoFile}"
echo "[INFO] Using abundance file: ${abdFile}"
echo "[INFO] Using covariate file: ${covFile}"
echo "[INFO] Using GRM file: ${grmFile}"





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
        step0_generateGRM.R \
            --genoFile=${genoFile} \
            --grmFile=${grmFile}
    fi
fi

# step1: fit null model for all phenotypes
echo "Fitting null model for all phenotypes..."
if [[ "${PALMmethod}" == 1 ]]; then
    step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${palm1_step1_prefix}
else
    step1_fitNULL.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --grmFile=${grmFile} \
        --covarColList=all \
        --outputPrefix=${palm2_step1_prefix} \
        --useGRMtoFitNULL=TRUE \
        --sampleIDColinabdFile=${sampleIDColinabdFile} \
        --sampleIDColincovFile=${sampleIDColincovFile}
fi

# step2: score test for phenoCol
echo "Performing score test..."
if [[ "${PALMmethod}" == 1 ]]; then
    step2_palm.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm1_step1_prefix}.rda \
        --PALMOutputFile=${outputFolder} \
        --chrom=${chrom} \
        --correct=NULL
else
    step2_scoreTest.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm2_step1_prefix}.rda \
        --PALMOutputFile=${palm2_step2_prefix} \
        --chrom=${chrom} \
        --phenoCol=${phenoCol} \
        --minMAF=0
fi

# Step 3: Generate feature information
step3_info.R \
    --abdFile=${abdFile} \
    --genoFile=${genoFile} \
    --outputFile=${outputFolder}/feature_info.txt