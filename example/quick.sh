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





################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"
# step0: generate GRM from genotype data
echo "Generating GRM from genotype data..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_generateGRM.R \
    --genoFile=${genoFile} \
    --grmFile=${grmFile}

# step1 & step2: fit null model and score test for one or more phenotypes
# Support multiple phenotypes separated by commas in `phenoCol`, or phenoCol=all to use all phenos from abdFile
if [ -z "${phenoCol}" ]; then
    echo "phenoCol is empty. Please set phenoCol to one or more phenotype column names, or phenoCol=all."
    exit 1
fi

# Check if phenoCol is "all" - if so, extract all phenotype column names from abdFile header
if [ "${phenoCol}" = "all" ]; then
    echo "phenoCol is 'all'; extracting all phenotype columns from abdFile header..."
    # Read first line (header) from abdFile and split by tab
    header=$(head -n 1 "${abdFile}")
    # Convert to array by splitting on tab
    IFS=$'\t' read -ra HEADER_COLS <<< "${header}"
    
    # Build PHENOS array, excluding IID and metadata columns
    PHENOS=()
    for col in "${HEADER_COLS[@]}"; do
        # Skip IID, SeqDepth, and other known metadata columns
        if [[ "${col}" != "${sampleIDColinabdFile}" && "${col}" != "" ]]; then
            PHENOS+=("${col}")
        fi
    done
    
    if [ ${#PHENOS[@]} -eq 0 ]; then
        echo "No phenotype columns found in abdFile (after excluding IID, etc.)"
        exit 1
    fi
    
    echo "Found ${#PHENOS[@]} phenotype columns: ${PHENOS[*]}"
else
    # remove spaces around commas and split into array
    phenoColClean=$(echo "${phenoCol}" | sed 's/[[:space:]]//g')
    IFS=',' read -ra PHENOS <<< "${phenoColClean}"
fi
for pheno in "${PHENOS[@]}"; do
    if [ -z "${pheno}" ]; then
        continue
    fi
    # step1: fit null model
    echo "Fitting null model for ${pheno}..."
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_fitNULL.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --sampleIDColincovFile=${sampleIDColincovFile} \
        --grmFile=${grmFile} \
        --phenoCol=${pheno} \
        --covarColList=${covarColList} \
        --outputPrefix=./${outputFolder}/${pheno}_step1 \
        --useGRMtoFitNULL=TRUE \
        --sampleIDColinabdFile=${sampleIDColinabdFile} \
        --sampleIDColincovFile=${sampleIDColincovFile} \
        --traitType=count

    # step2: score test
    step1prefix=./${outputFolder}/${pheno}_step1
    step2prefix=./${outputFolder}/${pheno}_step2
    echo "Performing score test for ${pheno}..."
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_scoreTest.R \
        --inFile=${genoFile} \
        --PALMOutputFile=${step2prefix}.txt \
        --NULLmodelFile=${step1prefix}.rda \
        --minMAF=0
done
