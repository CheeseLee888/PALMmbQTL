#!/usr/bin/env bash
set -euo pipefail
################################# parameter settings below #################################
inputFolder=input
outputFolder=output
genoFile=${inputFolder}/geno
abdFile=${inputFolder}/abd.txt
outputSnpFile=${outputFolder}/info_snp.txt
outputFeatureFile=${outputFolder}/info_feature.txt

################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

echo "Start: Generate SNP information."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputSnpFile}
echo "Finish: Generate SNP information."

echo "Start: Generate feature information."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_featureInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputFeatureFile}
echo "Finish: Generate feature information."