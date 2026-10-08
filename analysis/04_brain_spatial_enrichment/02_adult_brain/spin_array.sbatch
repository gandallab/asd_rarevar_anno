#!/bin/bash
#SBATCH --job-name=spin_tests
#SBATCH --output=logs/spin_tests_%A_%a.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-3

# =============================================================================
# spin_array.sbatch  --  spin tests for the 4 gene sets, as a SLURM array
#
# One array task per gene set (ASD, SYN, GR, MORPH). Each task emits that gene
# set's 6 tests, using the shared spins + gradient nulls that spin_prep.sbatch
# already cached:
#     4 correlations  : map vs AP, C1, C2, C3      (<prefix>spin_results_*.csv)
#     occipital focal : mean over the 4 occipital regions, as a group (union)
#                                                    (<prefix>occipital_focal_*.csv)
#     per-region      : spin-corrected enrichment for every cortical region
#                       (separate mode over all regions) (<prefix>cortex_perregion_*.csv)
# All tasks write into one flat outputs/analysis/04_brain_spatial_enrichment/, kept apart
# by --prefix naming the notebook that consumes them: 04_01_ for ASD, 04_02_ for the topics.
# 4 gene sets x 6 tests = 24 permutation tests. The 4 correlations come out of a
# single run (they share the gradient nulls), so the natural array unit is the
# gene set; occipital (union) and per-region (separate) are two run_spin_test
# calls per task. Depends on spin_prep.sbatch -- submit via submit_spin_tests.sh.
# =============================================================================

set -euo pipefail
# inside a SLURM job $0 is the spooled copy, so cd to the launch dir (code/) explicitly
cd "$SLURM_SUBMIT_DIR"

# Load spin env
#module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate spin
export PYTHONUNBUFFERED=1     # stream progress to the log live instead of block-buffering

# gene set for this array index, and its regional-enrichment map
GENE_SETS=(ASD SYN GR MORPH)
PREFIXES=(04_01_ASD_ 04_02_SYN_ 04_02_GR_ 04_02_MORPH_)
GS=${GENE_SETS[$SLURM_ARRAY_TASK_ID]}
PREFIX=${PREFIXES[$SLURM_ARRAY_TASK_ID]}
# project root is three levels up from code/analysis/04_brain_spatial_enrichment
OUT=$(cd ../../.. && pwd)/outputs/analysis/04_brain_spatial_enrichment
MAP=$OUT/${PREFIX}regional_enrichment.csv
mkdir -p "$OUT"
echo "[task $SLURM_ARRAY_TASK_ID] gene set=$GS  map=$MAP"

# (1) 4 correlations (AP/C1/C2/C3) + occipital focal test (the 4 occipital
#     regions tested together as one ROI, i.e. --roi-mode union)
python run_spin_test.py \
    --n-spin 5000 --seed 0 \
    --asd-map "$MAP" --asd-cols z_vs_bg \
    --roi-mode union --roi-name occipital \
    --prefix "$PREFIX" --outdir "$OUT" --cachedir "$OUT/spin_cache"

# (2) per-region spin-corrected enrichment: --roi-mode separate over ALL regions,
#     one spin-corrected p (+ FDR across regions) per cortical region
python run_spin_test.py \
    --n-spin 5000 --seed 0 \
    --asd-map "$MAP" --asd-cols z_vs_bg \
    --roi-mode separate --roi-regions all --roi-name cortex \
    --prefix "$PREFIX" --outdir "$OUT" --cachedir "$OUT/spin_cache"

echo "[task $SLURM_ARRAY_TASK_ID] $GS done -> $OUT/${PREFIX}*"
