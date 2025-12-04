#!/usr/bin/env bash
set -euo pipefail

# 
inputFolder=input
outputFolder=output

abdFile=./${inputFolder}/abd.tsv
covFile=./${inputFolder}/cov.tsv
sampleIDCol=sample_id

# 
covariateInterest=sex
covariateAdjust="age,PC1,PC2,batch"

outPrefix=./${outputFolder}/PALM_${covariateInterest}

mkdir -p "${outputFolder}"

echo "[INFO] Running PALM on example data ..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/run_PALM.R \
  --abdFile "${abdFile}" \
  --covFile "${covFile}" \
  --sampleIDCol "${sampleIDCol}" \
  --covariateInterest "${covariateInterest}" \
  --covariateAdjust "${covariateAdjust}" \
  --outPrefix "${outPrefix}"

echo "[INFO] Done. Results in ${outPrefix}.rds and ${outPrefix}.tsv"
