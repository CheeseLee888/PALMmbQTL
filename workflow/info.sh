#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"


################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

echo "Start: Generate SNP information."
pixi run --manifest-path=pixi.toml Rscript extdata/step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputSnpFile}
echo "Finish: Generate SNP information."

echo "Start: Generate feature information."
pixi run --manifest-path=pixi.toml Rscript extdata/step0_featureInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputFeatureFile}
echo "Finish: Generate feature information."

echo "Start: Generate sample information."
pixi run --manifest-path=pixi.toml Rscript extdata/step0_sampleInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputSampleFile}
echo "Finish: Generate sample information."