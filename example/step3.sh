#!/usr/bin/env bash
set -euo pipefail
################################# parameter settings below #################################
inFile=./output/palm2_step2_chr1_g1.txt
outdir=./plot



################################# workflow below (do not modify) #################################
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step3_plot.R \
    --inFile=${inFile} \
    --outdir=${outdir}