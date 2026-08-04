#!/bin/bash
#SBATCH --job-name=scDRS-downstream-subclass-donor-region
#SBATCH -n 32               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-18:00          # Runtime in D-HH:MM
#SBATCH --mem=90G          # Memory pool for all cores
#SBATCH -o /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_downstream_subclass_donor_region_%j.out
#SBATCH -e /mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/logs/scDRS_downstream_subclass_donor_region_%j.err

# Per-donor group-level Z-score: group key = (EN/IN broad subclass) x donor x region,
# pooling ALL types and ALL age Groups within each donor -- built by the
# "subclass_donor_region" section of analysis/01_scdrs/01_input_prep/01_2_prepare_downstream_metadata.ipynb, on obj_rna_raw_subclass_donor_region_paired.h5ad (restricted to the 11-donor paired cohort). 

echo "JobID: ${SLURM_JOB_ID}"
date
source "$HOME/envs/scdrs_env/bin/activate"
DATA_PATH=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data

TRAIT="ASC_Pfdr"
GENE_NUMBER=253

H5AD_FILE=${DATA_PATH}/SnMultiome_Wang2025/obj_rna_raw_subclass_donor_region_paired.h5ad
SCORE_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/score_file/${GENE_NUMBER}_genes
RESULT_FOLDER=/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/downstream_analysis/${GENE_NUMBER}_genes/subclass_donor_region_paired

mkdir -p "${RESULT_FOLDER}"

echo "Running scDRS downstream of ${TRAIT} with ${GENE_NUMBER} genes, group=subclass_donor_region (paired-donor cohort)..."
echo "H5AD_FILE: ${H5AD_FILE}"
echo "SCORE_FOLDER: ${SCORE_FOLDER}"
echo "RESULT_FOLDER: ${RESULT_FOLDER}"

scdrs perform-downstream \
    --h5ad-file "${H5AD_FILE}" \
    --score-file "${SCORE_FOLDER}"/"${TRAIT}".full_score.gz \
    --out-folder "${RESULT_FOLDER}" \
    --group-analysis subclass_donor_region \
    --flag-filter-data True \
    --flag-raw-count False

echo "Finished downstream analysis: GENE_NUMBER = ${GENE_NUMBER} (subclass_donor_region, paired-donor cohort)"
date
