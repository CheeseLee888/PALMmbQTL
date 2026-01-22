#!/bin/bash
#SBATCH --job-name=palm-step3
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/scratch/group/ztang2/pli297/logs/%x_%A_%a.out
#SBATCH --error=/mnt/scratch/group/ztang2/pli297/logs/%x_%A_%a.err

set -euo pipefail

################################ parameter settings below #################################
WORK=/mnt/scratch/group/ztang2/pli297
SIF=/mnt/scratch/group/ztang2/pli297/palmmbqtl.sif

# input file list (one file per line)
FILELIST=${WORK}/output/wisc_palm2/plot_files.txt

# output folder for plots
OUTDIR=${WORK}/output/wisc_palm2_plot

mkdir -p ${WORK}/logs
mkdir -p ${OUTDIR}

################################ workflow below (do not modify) #################################
# pick the infile for this array task
INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FILELIST}")

if [[ -z "${INFILE}" ]]; then
  echo "ERROR: empty INFILE for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID} processing: ${INFILE}"

# run plotting (single-file mode)
apptainer exec --cleanenv --bind ${WORK}:${WORK} ${SIF} \
    step3_plot.R \
    --inFile="${INFILE}" \
    --outdir="${OUTDIR}"

echo "Done: ${INFILE}"
