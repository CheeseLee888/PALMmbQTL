#!/bin/bash
#SBATCH --job-name=palm-step2
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
# --- parallel setting ---
chrom=${SLURM_ARRAY_TASK_ID}

######### optional below #########
palm1_step1_prefix=palm1_step1_allpheno
palm1_step2_prefix=palm1_step2${chrom:+_chr${chrom}}
palm2_step1_prefix=palm2_step1_allpheno
palm2_step2_prefix=palm2_step2${chrom:+_chr${chrom}}



###################### change to absolute paths in Platform R ########################
WORK=/mnt/scratch/group/ztang2/pli297
SIF=/mnt/scratch/group/ztang2/pli297/palmmbqtl.sif

inputFolder="${WORK}/example/${inputFolder}"
outputFolder="${WORK}/example/${outputFolder}"
genoFile="${inputFolder}/${genoFile}"
palm1_step1_prefix="${outputFolder}/${palm1_step1_prefix}"
palm1_step2_prefix="${outputFolder}/${palm1_step2_prefix}"
palm2_step1_prefix="${outputFolder}/${palm2_step1_prefix}"
palm2_step2_prefix="${outputFolder}/${palm2_step2_prefix}"

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

# step2: score test for phenoCol
echo "Performing score test..."
if [[ "${PALMmethod}" == 1 ]]; then
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step2_palm.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm1_step1_prefix}.rda \
        --PALMOutputFile=${palm1_step2_prefix} \
        --chrom=${chrom} \
        --correct=NULL \
        --cluster=NULL
else
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step2_scoreTest.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm2_step1_prefix}.rda \
        --PALMOutputFile=${palm2_step2_prefix} \
        --chrom=${chrom} \
        --phenoCol=${phenoCol} \
        --minMAF=0
fi