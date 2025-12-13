################################# parameter settings below #################################
######### model option ##########
model=palm # 'palm' or 'palmmbqtl'

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
step2_prefix_palm=./${outputFolder}/${phenoCol}_step2_palm${chrom:+_chr${chrom}}




################################# workflow below (do not modify) #################################
if(${model}!='palm' and ${model}!='palmmbqtl'){
    stop("Model must be specified as palm or palmmbqtl.")
}

mkdir -p "${outputFolder}"

# step0: generate GRM from genotype data
if(${model}=='palm'){
    echo "Generating GRM from genotype data..."
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_generateGRM.R \
        --genoFile=${genoFile} \
        --grmFile=${grmFile}
}

# step1: fit null model for all phenotypes
echo "Fitting null model for all phenotypes..."
if(${model}=='palm'){
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_palm.R \
        --rel.abd "${abdFile}" \
        --covariate.adjust "${covFile}" \
        --depth "${offsetCol}" \
        --outputPrefix=${step1_prefix_palm}
}else{
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step1_fitNULL.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --grmFile=${grmFile} \
        --covarColList=${covarColList} \
        --outputPrefix=${step1_prefix} \
        --useGRMtoFitNULL=TRUE \
        --sampleIDColinabdFile=${sampleIDColinabdFile} \
        --sampleIDColincovFile=${sampleIDColincovFile} \
}


# step2: score test for phenoCol
echo "Performing score test for ${phenoCol}..."
if(${model}=='palm'){
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_palm.R \
        --null.obj=${step1_prefix_palm}.rda \
        --covariate.interest=${genoFile} \
        --correct=NULL \
        --PALMOutputFile=${step2_prefix_palm}.txt
}else{
    pixi run --manifest-path=../pixi.toml Rscript ../extdata/step2_scoreTest.R \
        --inFile=${genoFile} \
        --phenoCol=${phenoCol} \
        --chrom=${chrom} \
        --PALMOutputFile=${step2_prefix}.txt \
        --NULLmodelFile=${step1_prefix}.rda \
        --minMAF=0
}