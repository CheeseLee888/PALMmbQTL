################################# parameter settings below #################################
######### required below #########
inputFolder=input
outputFolder=output

genoFile=./${inputFolder}/geno
abdFile=./${inputFolder}/abd.txt
covFile=./${inputFolder}/cov.txt
sampleIDColinabdFile=IID
sampleIDColincovFile=IID

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

# step1: fit null model for all phenotypes
echo "Fitting null model for all phenotypes..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_fitNULL.R \
--abdFile=${abdFile} \
--covFile=${covFile} \
--sampleIDColincovFile=${sampleIDColincovFile} \
--grmFile=${grmFile} \
--covarColList=${covarColList} \
--outputPrefix=./${outputFolder}/${pheno}_step1 \
--useGRMtoFitNULL=TRUE \
--sampleIDColinabdFile=${sampleIDColinabdFile} \
--sampleIDColincovFile=${sampleIDColincovFile} \
--traitType=count


for pheno in "${PHENOS[@]}"; do

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
