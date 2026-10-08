#!/bin/bash
#SBATCH --job-name=scDRS_group
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 0-18:00
#SBATCH --mem=90G
#SBATCH --array=0-3
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_group_%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_group_%A_%a.err

echo "JobID: ${SLURM_JOB_ID}, ArrayID: ${SLURM_ARRAY_TASK_ID}"
date
source "$HOME/envs/scdrs_env/bin/activate"

DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data
H5AD_FILE=${DATA_PATH}/SnMultiome_Wang2025/obj_rna_raw.h5ad
COV_FILE=${DATA_PATH}/SnMultiome_Wang2025/donorID_sex_ngene.cov

# Array 0 = p_hat DDID (ASDDH, 121 genes)
# Array 1 = p_hat NoDDID (ASDDL + ASDDA, 128 genes)
# Array 2 = per-proband NoDDID (ASDDA alone, 44 genes, FDR<0.001)
# Array 3 = per-proband NoDDID (ASDDA alone, 105 genes, FDR<0.01)
case "${SLURM_ARRAY_TASK_ID}" in
    0)
        GROUP="ddid"
        SUFFIX="new_DMN"
        ;;
    1)
        GROUP="noddid"
        SUFFIX="new_DMN"
        ;;
    2)
        GROUP="noddid"
        SUFFIX="perproband_new_DMN"
        ;;
    3)
        GROUP="noddid"
        SUFFIX="perproband_new_DMN_fdr01"
        ;;
    *)
        echo "ERROR: unexpected SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
        exit 1
        ;;
esac

P_FILE=${DATA_PATH}/RVAS_result/full_results_wcounts_pfdr_${GROUP}_${SUFFIX}.tsv
TRAIT_SCORE="ASC_TADA_Pfdr_${GROUP}_${SUFFIX}"
TRAIT="ASC_Pfdr_${GROUP}_${SUFFIX}"

echo "Running for GROUP=${GROUP}, SUFFIX=${SUFFIX}"

# Count number of genes in the file (excluding header)
GENE_NUMBER=$(tail -n +2 ${P_FILE} | wc -l)
echo "Number of genes in ${P_FILE}: ${GENE_NUMBER}"

GS_DIR=${DATA_PATH}/scDRS/gs_file
GS_FILE=${GS_DIR}/${TRAIT_SCORE}.${GENE_NUMBER}.gs
OUT_FOLDER="/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_ddid_noddid/score_file/${GROUP}_genes"

mkdir -p "${GS_DIR}"
mkdir -p "${OUT_FOLDER}"

# Generate GS file if it doesn't exist
GS_FILE=${GS_DIR}/${TRAIT_SCORE}.${GENE_NUMBER}.gs

if [ ! -f "${GS_FILE}" ]; then
    echo "Generating GS file: ${GS_FILE}"
    scdrs munge-gs \
        --out-file "${GS_FILE}" \
        --pval-file "${P_FILE}" \
        --weight zscore \
        --n-max "${GENE_NUMBER}"
else
    echo "GS file already exists: ${GS_FILE}"
fi

# Compute scDRS Score
echo "Running scDRS compute-score for ${GROUP} with ${GENE_NUMBER} genes..."
scdrs compute-score \
    --h5ad-file ${H5AD_FILE} \
    --h5ad-species human \
    --gs-file ${GS_FILE} \
    --gs-species human \
    --out-folder ${OUT_FOLDER} \
    --cov-file ${COV_FILE} \
    --flag-filter-data True \
    --flag-raw-count False \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

echo "Finished scDRS Score computation for ${GROUP}"

# Perform Downstream Analysis
echo "Running scDRS downstream for ${GROUP}..."
RESULT_FOLDER="/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_ddid_noddid/downstream_analysis/${GROUP}_genes"
mkdir -p ${RESULT_FOLDER}

scdrs perform-downstream \
    --h5ad-file ${H5AD_FILE} \
    --score-file ${OUT_FOLDER}/${TRAIT}.full_score.gz \
    --out-folder ${RESULT_FOLDER} \
    --group-analysis type,subclass,class,region_summary,region_combine,Group \
    --flag-filter-data True \
    --flag-raw-count False

scdrs perform-downstream \
    --h5ad-file ${H5AD_FILE} \
    --score-file ${OUT_FOLDER}/${TRAIT}.full_score.gz \
    --out-folder ${RESULT_FOLDER} \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count False

echo "Finished downstream analysis for ${GROUP}"
date