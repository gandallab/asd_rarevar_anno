#!/bin/bash
#SBATCH --job-name=scDRS-GWAS
#SBATCH -n 32               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-18:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=150G          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-1          
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_GWAS_%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_GWAS_%A_%a.err 

echo "JobID: ${SLURM_JOB_ID}, ArrayID: ${SLURM_ARRAY_TASK_ID}"
date
source $HOME/envs/scdrs_env/bin/activate

TRAIT="ASD_Matoba_2020"
DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data
Z_FILE=${DATA_PATH}/GWAS_result/MAGMA/step2/ASD_Matoba_2020.genes.tsv
H5AD_FILE=${DATA_PATH}/SnMultiome_Wang2025/obj_rna_raw.h5ad
COV_FILE=${DATA_PATH}/SnMultiome_Wang2025/donorID_sex_ngene.cov 

# fdr < 0.05 gives genes < 100, so follow n-min 100 to ensure enough genes for scDRS
GENE_NUMBERS=(100 200)
GENE_NUMBER=${GENE_NUMBERS[$SLURM_ARRAY_TASK_ID]}

echo "Finished gs file generation with top ${GENE_NUMBER} genes"

GS_DIR=${DATA_PATH}/scDRS/gs_file
mkdir -p "${GS_DIR}"

# Generate GS file if it doesn't exist
GS_FILE=${GS_DIR}/${TRAIT}.${GENE_NUMBER}.gs

if [ ! -f "${GS_FILE}" ]; then
    echo "Generating GS file: ${GS_FILE}"
    scdrs munge-gs \
        --out-file "${GS_FILE}" \
        --zscore-file "${Z_FILE}" \
        --n-max "${GENE_NUMBER}"  \
        --weight zscore
else
    echo "GS file already exists: ${GS_FILE}"
fi

# Compute scDRS Score
echo "Running scDRS compute-score with ${GENE_NUMBER} genes..."

OUT_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_GWAS_v3/score_file/${GENE_NUMBER}_genes
mkdir -p ${OUT_FOLDER}

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

# echo "Finished scDRS Score computation for GENE_NUMBER = ${GENE_NUMBER}"

# Perform Downstream Analysis
echo "Running scDRS downstream with ${GENE_NUMBER} genes..."
SCORE_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_GWAS_v3/score_file/${GENE_NUMBER}_genes
RESULT_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_GWAS_v3/downstream_analysis/${GENE_NUMBER}_genes

mkdir -p ${RESULT_FOLDER}

scdrs perform-downstream \
    --h5ad-file $H5AD_FILE \
    --score-file ${SCORE_FOLDER}/${TRAIT}.full_score.gz \
    --out-folder $RESULT_FOLDER \
    --group-analysis type,subclass,class,region_summary,region_combine,Group \
    --flag-filter-data True \
    --flag-raw-count False

scdrs perform-downstream \
    --h5ad-file $H5AD_FILE \
    --score-file ${SCORE_FOLDER}/${TRAIT}.full_score.gz \
    --out-folder $RESULT_FOLDER \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count False
echo "Finished downstream analysis GENE_NUMBER = ${GENE_NUMBER}"
date
