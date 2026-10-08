#!/bin/bash
#SBATCH --job-name=scDRS-group-downstream
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 0-18:00
#SBATCH --mem=90G
#SBATCH --array=0-7         # 8 tasks: 4 gene-set groups × 2 lineage types (EN/IN)
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_group_downstream_%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_group_downstream_%A_%a.err

echo "JobID: ${SLURM_JOB_ID}, ArrayID: ${SLURM_ARRAY_TASK_ID}"
date
source "$HOME/envs/scdrs_env/bin/activate"

DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data

# Array layout: 4 gene-set groups × 2 lineage types
# 0 = p_hat DDID (ASDDH)              + EN
# 1 = p_hat DDID (ASDDH)              + IN
# 2 = p_hat NoDDID (ASDDL+ASDDA)      + EN
# 3 = p_hat NoDDID (ASDDL+ASDDA)      + IN
# 4 = per-proband NoDDID (ASDDA, FDR<0.001) + EN
# 5 = per-proband NoDDID (ASDDA, FDR<0.001) + IN
# 6 = per-proband NoDDID (ASDDA, FDR<0.01)  + EN
# 7 = per-proband NoDDID (ASDDA, FDR<0.01)  + IN

# Gene-set group (0-3) and lineage type (EN/IN)
GS_INDEX=$(( SLURM_ARRAY_TASK_ID / 2 ))
if (( SLURM_ARRAY_TASK_ID % 2 == 0 )); then
    TYPE="EN"
else
    TYPE="IN"
fi

case "${GS_INDEX}" in
    0) GROUP="ddid";   SUFFIX="new_DMN" ;;
    1) GROUP="noddid"; SUFFIX="new_DMN" ;;
    2) GROUP="noddid"; SUFFIX="perproband_new_DMN" ;;
    3) GROUP="noddid"; SUFFIX="perproband_new_DMN_fdr01" ;;
    *)
        echo "ERROR: unexpected SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
        exit 1
        ;;
esac

echo "Resolved parameters: GROUP=${GROUP}, SUFFIX=${SUFFIX}, TYPE=${TYPE}"

TRAIT="ASC_Pfdr_${GROUP}_${SUFFIX}"

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