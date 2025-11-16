pheno=g_Blautia

step1prefix=./output/${pheno}_step1
step2prefix=./output/${pheno}_allchr_step2

# pixi run --manifest-path=../pixi.toml Rscript step2_tests_qtl.R \
#     --bedFile=./input/stool_bialleic_merged_data.bed \
#     --bimFile=./input/stool_bialleic_merged_data.bim \
#     --famFile=./input/stool_bialleic_merged_data.fam \
#     --SAIGEOutputFile=${step2prefix} \
#     --LOCO=FALSE \
#     --GMMATmodelFile=${step1prefix}.rda \
#     --varianceRatioFile=${step1prefix}.varianceRatio.txt \
#     --is_overwrite_output=TRUE

pixi run --manifest-path=../pixi.toml Rscript step2_tests_qtl.R \
    --inFile=./input/stool_bialleic_merged_data.gds \
    --SAIGEOutputFile=${step2prefix}.txt \
    --LOCO=FALSE \
    --GMMATmodelFile=${step1prefix}.rda