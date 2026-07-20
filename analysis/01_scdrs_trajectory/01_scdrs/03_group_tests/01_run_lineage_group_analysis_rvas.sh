#!/bin/bash
#SBATCH --job-name=scDRS-downstream
#SBATCH -n 32               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-18:00          # Runtime in D-HH:MM
#SBATCH --mem=90G          # Memory pool for all cores
#SBATCH --array=0-7         # 8 tasks: 4 gene numbers × 2 types
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_downstream_%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_downstream_%A_%a.err 

echo "JobID: ${SLURM_JOB_ID}, ArrayID: ${SLURM_ARRAY_TASK_ID}"
date
source "$HOME/envs/scdrs_env/bin/activate"
DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data

TRAIT="ASC_Pfdr"
GENE_NUMBERS=(253 416 696 951) 
TYPES=(EN IN)

gene_idx=$((SLURM_ARRAY_TASK_ID % 4)) # 4 gene sets
type_idx=$((SLURM_ARRAY_TASK_ID / 4)) # 4 gene sets
GENE_NUMBER=${GENE_NUMBERS[$gene_idx]}
TYPE=${TYPES[$type_idx]}

echo "Resolved parameters: TYPE = ${TYPE}, GENE_NUMBER = ${GENE_NUMBER}"
H5AD_FILE=${DATA_PATH}/SnMultiome_Wang2025/obj_rna_raw_${TYPE}_slingshot.h5ad

SCORE_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/score_file/${GENE_NUMBER}_genes
RESULT_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/downstream_analysis/${GENE_NUMBER}_genes/${TYPE}

mkdir -p "${RESULT_FOLDER}"

echo "Running scDRS downstream of ${TRAIT} with ${GENE_NUMBER} genes for ${TYPE}..."
echo "H5AD_FILE: ${H5AD_FILE}"
echo "SCORE_FOLDER: ${SCORE_FOLDER}"
echo "RESULT_FOLDER: ${RESULT_FOLDER}"

scdrs perform-downstream \
    --h5ad-file "${H5AD_FILE}" \
    --score-file "${SCORE_FOLDER}"/"${TRAIT}".full_score.gz \
    --out-folder "${RESULT_FOLDER}" \
    --group-analysis pseudotime_bin20,pseudotime_bin50,pseudotime_bin100,pseudotime_bin200,mclust,Group \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count False
echo "Finished downstream analysis: TYPE = ${TYPE}, GENE_NUMBER = ${GENE_NUMBER}"
date

