# Navigate to SAIGEQTL directory first

# choose covariates 1-7 (aligned with Kirin's analysis)

pheno=g_Blautia
# pheno2=g_Bifidobacterium
pixi run --manifest-path=../pixi.toml Rscript step1_fitNULLGLMM_qtl.R \
    --phenoFile=./input/stool_bialleic_cov_pheno_merged.txt \
    --phenoCol=${pheno} \
    --covarColList=PC01,PC02,PC03,PC04,PC05,SAMPLE_COLLECTION_AGE_MONTHS_DEV,SEX \
    --sampleCovarColList=PC01,PC02,PC03,PC04,PC05,SAMPLE_COLLECTION_AGE_MONTHS_DEV,SEX \
    --sampleIDColinphenoFile=IID \
    --offsetCol=SeqDepth \
    --traitType=count \
    --outputPrefix=./output/${pheno}_step1 \
    --skipVarianceRatioEstimation=TRUE \
    --useGRMtoFitNULL=FALSE \
    --isCovariateOffset=TRUE \
    --isCovariateTransform=FALSE \
    --plinkFile=./input/stool_bialleic_merged_data \
    --LOCO=FALSE \
    --isShrinkModelOutput=FALSE \
    --useGRMtoFitNULL=TRUE \
    --grmFile=./input/stool_bialleic_merged_grm.rds
