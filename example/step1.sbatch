#!/bin/bash
#SBATCH --job-name=palm-step1
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/scratch/group/ztang2/pli297/logs/%x_%A_%a.out
#SBATCH --error=/mnt/scratch/group/ztang2/pli297/logs/%x_%A_%a.err
set -euo pipefail

################################ parameter settings below #################################
######### required below #########
PALMmethod=1 # 1 for 'PALM' or 2 for 'PALM-mbQTL'

inputFolder=input
outputFolder=output

genoFile=geno
abdFile=abd.txt
covFile=cov.txt

######### optional below #########
grmFile=grm.rds
palm1_step1_prefix=palm1_step1_allpheno
palm2_step1_prefix=palm2_step1_allpheno



###################### change to absolute paths in Platform R ########################
WORK=/mnt/scratch/group/ztang2/pli297
SIF=/mnt/scratch/group/ztang2/pli297/palmmbqtl.sif

inputFolder="${WORK}/example/${inputFolder}"
outputFolder="${WORK}/example/${outputFolder}"
genoFile="${inputFolder}/${genoFile}"
abdFile="${inputFolder}/${abdFile}"
covFile="${inputFolder}/${covFile}"
grmFile="${inputFolder}/${grmFile}"
palm1_step1_prefix="${outputFolder}/${palm1_step1_prefix}"
palm2_step1_prefix="${outputFolder}/${palm2_step1_prefix}"

################################# workflow below (do not modify) #################################
if [[ "${PALMmethod}" != 1 && "${PALMmethod}" != 2 ]]; then
    echo "PALMmethod must be specified as 1 (PALM) or 2 (PALM-mbQTL)."
    exit 1
fi

if [[ "${PALMmethod}" == 1 ]]; then
    echo "Running PALM method..."
else
    echo "Running PALM-mbQTL method..."
fi

mkdir -p "${outputFolder}"

# step0: check input files
echo "Checking input files..."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_checkInput.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --genoFile=${genoFile}


# step0: generate GRM from genotype data (only for PALM-mbQTL)
if [[ "${PALMmethod}" == 2 ]]; then
    if [[ -f "${grmFile}" ]]; then
        echo "GRM already exists at: ${grmFile}"
        echo "Skip generating GRM and reuse the existing GRM."
    else
        echo "Generating GRM from genotype data..."
        apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
            step0_generateGRM.R \
            --genoFile=${genoFile} \
            --grmFile=${grmFile}
    fi
fi

# step1: fit null model for all phenotypes
echo "Fitting null model for all phenotypes..."
if [[ "${PALMmethod}" == 1 ]]; then
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${palm1_step1_prefix}
else
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step1_fitNULL.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --grmFile=${grmFile} \
        --outputPrefix=${palm2_step1_prefix} \
        --covarColList=all \
        --useGRMtoFitNULL=TRUE
fi