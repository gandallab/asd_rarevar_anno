#!/bin/bash
#SBATCH --job-name=scDRS-downstream-mclust-group-region
#SBATCH -n 32               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-18:00          # Runtime in D-HH:MM
#SBATCH --mem=90G          # Memory pool for all cores
#SBATCH --array=0-1         # 2 tasks: 2 lineages (EN, IN), 253 genes only
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_downstream_mclust_group_region_%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_downstream_mclust_group_region_%A_%a.err

# Reuses the scDRS score files already computed by analysis/01_scdrs/02_compute/
# 01_run_scdrs_rvas.sh -- only perform-downstream is (re-)run here, against the
# mclust_Group_region group key built by the "mclust_Group_region" section of
# analysis/01_scdrs/01_input_prep/01_2_prepare_downstream_metadata.ipynb (added directly
# onto obj_rna_raw_{EN,IN}_slingshot.h5ad -- no separate file). 

echo "JobID: ${SLURM_JOB_ID}, ArrayID: ${SLURM_ARRAY_TASK_ID}"
date
source "$HOME/envs/scdrs_env/bin/activate"
DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data

TRAIT="ASC_Pfdr"
GENE_NUMBER=253
TYPES=(EN IN)

TYPE=${TYPES[$SLURM_ARRAY_TASK_ID]}

echo "Resolved parameters: TYPE = ${TYPE}, GENE_NUMBER = ${GENE_NUMBER}"
H5AD_FILE=${DATA_PATH}/SnMultiome_Wang2025/obj_rna_raw_${TYPE}_slingshot.h5ad

SCORE_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/score_file/${GENE_NUMBER}_genes
RESULT_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/downstream_analysis/${GENE_NUMBER}_genes/${TYPE}

mkdir -p "${RESULT_FOLDER}"

echo "Running scDRS downstream of ${TRAIT} with ${GENE_NUMBER} genes, group=mclust_Group_region, TYPE=${TYPE}..."
echo "H5AD_FILE: ${H5AD_FILE}"
echo "SCORE_FOLDER: ${SCORE_FOLDER}"
echo "RESULT_FOLDER: ${RESULT_FOLDER}"

scdrs perform-downstream \
    --h5ad-file "${H5AD_FILE}" \
    --score-file "${SCORE_FOLDER}"/"${TRAIT}".full_score.gz \
    --out-folder "${RESULT_FOLDER}" \
    --group-analysis mclust_Group_region \
    --flag-filter-data True \
    --flag-raw-count False

echo "Finished downstream analysis: TYPE = ${TYPE}, GENE_NUMBER = ${GENE_NUMBER}"
date
