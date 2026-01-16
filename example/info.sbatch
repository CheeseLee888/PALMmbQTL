#!/bin/bash
#SBATCH --job-name=palm-info
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/scratch/group/ztang2/pli297/logs/%x_%A_%a.out
#SBATCH --error=/mnt/scratch/group/ztang2/pli297/logs/%x_%A_%a.err
set -euo pipefail

################################ parameter settings below #################################
inputFolder=input
outputFolder=output
genoFile=geno
abdFile=abd.txt
outputSnpFile=info_snp.txt
outputFeatureFile=info_feature.txt



###################### change to absolute paths in Platform R ########################
WORK=/mnt/scratch/group/ztang2/pli297
SIF=/mnt/scratch/group/ztang2/pli297/palmmbqtl.sif

inputFolder="${WORK}/example/${inputFolder}"
outputFolder="${WORK}/example/${outputFolder}"
genoFile="${inputFolder}/${genoFile}"
abdFile="${inputFolder}/${abdFile}"
outputSnpFile="${outputFolder}/${outputSnpFile}"
outputFeatureFile="${outputFolder}/${outputFeatureFile}"

################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"
echo "Generating SNP information from genotype data..."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputSnpFile}

echo "Generating feature information from phenotype data..."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_featureInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputFeatureFile}