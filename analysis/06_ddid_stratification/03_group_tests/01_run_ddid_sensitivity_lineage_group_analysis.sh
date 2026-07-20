#!/bin/bash
#SBATCH --job-name=scDRS-group-downstream
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 0-18:00
#SBATCH --mem=90G
#SBATCH --array=0-3         # 4 tasks: 2 groups (ddid/noddid) × 2 types (EN/IN)
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_group_downstream_%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_group_downstream_%A_%a.err

echo "JobID: ${SLURM_JOB_ID}, ArrayID: ${SLURM_ARRAY_TASK_ID}"
date
source "$HOME/envs/scdrs_env/bin/activate"

DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data

# Array layout:
# 0 = ddid   + EN
# 1 = ddid   + IN
# 2 = noddid + EN
# 3 = noddid + IN
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    GROUP="ddid";   TYPE="EN"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    GROUP="ddid";   TYPE="IN"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 2 ]; then
    GROUP="noddid"; TYPE="EN"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 3 ]; then
    GROUP="noddid"; TYPE="IN"
else
    echo "ERROR: unexpected SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Resolved parameters: GROUP = ${GROUP}, TYPE = ${TYPE}"

# TRAIT="ASC_Pfdr_${GROUP}_new_DMN"
TRAIT="ASC_Pfdr_${GROUP}_perproband_new_DMN"

H5AD_FILE=${DATA_PATH}/SnMultiome_Wang2025/obj_rna_raw_${TYPE}_slingshot.h5ad
SCORE_FOLDER="/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_ddid_noddid/score_file/${GROUP}_genes"
RESULT_FOLDER="/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_ddid_noddid/downstream_analysis/${GROUP}_genes/${TYPE}"

mkdir -p "${RESULT_FOLDER}"

echo "H5AD_FILE:     ${H5AD_FILE}"
echo "SCORE_FOLDER:  ${SCORE_FOLDER}"
echo "RESULT_FOLDER: ${RESULT_FOLDER}"

# Check score file exists
SCORE_FILE="${SCORE_FOLDER}/${TRAIT}.full_score.gz"
if [ ! -f "${SCORE_FILE}" ]; then
    echo "ERROR: score file not found: ${SCORE_FILE}"
    exit 1
fi

echo "Running scDRS downstream of ${TRAIT} for TYPE=${TYPE}, GROUP=${GROUP}..."

scdrs perform-downstream \
    --h5ad-file "${H5AD_FILE}" \
    --score-file "${SCORE_FILE}" \
    --out-folder "${RESULT_FOLDER}" \
    --group-analysis pseudotime_bin20,pseudotime_bin50,pseudotime_bin100,pseudotime_bin200,mclust,Group,mclust_Group \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count False

if [ $? -ne 0 ]; then
    echo "ERROR: perform-downstream failed for GROUP=${GROUP}, TYPE=${TYPE}"
    exit 1
fi

echo "Finished downstream analysis: GROUP = ${GROUP}, TYPE = ${TYPE}"
date