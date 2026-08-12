#!/bin/bash
#SBATCH --job-name=spin_prep
#SBATCH --output=logs/spin_prep_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4

# =============================================================================
# spin_prep.sbatch  --  one-time shared prep for the spin-test array
#
# Generates the artifacts that are the SAME for every gene set, so the array
# tasks don't each recompute them (the 5000-spin, 41k surface array is ~3.3 GB):
#   * spins        (outputs/analysis/03_brain_spatial_enrichment/spin_cache/spins_41k_n5000_seed0.npy)
#   * gradient nulls for AP/C1/C2/C3  (nulls_41k_n5000_seed0_4comp.npy)
# We trigger this by running the full pipeline once on the ASD map; the array
# tasks (spin_array.sbatch) then load these from cache. Submit via
# submit_spin_tests.sh, which chains the array after this job with afterok.
# =============================================================================

set -euo pipefail

# run from the directory sbatch was launched in (submit_spin_tests.sh cd's to code/
# first, so this is the dir holding run_spin_test.py, whose paths are ../data ../outputs).
# NB: inside a SLURM job $0 is the spooled copy, so use SLURM_SUBMIT_DIR, not dirname $0.
cd "$SLURM_SUBMIT_DIR"
OUT=$(cd ../../.. && pwd)/outputs/analysis/03_brain_spatial_enrichment
mkdir -p logs "$OUT/spin_cache"

# activate the env with neuromaps + statsmodels ('spin' env)
#module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate spin
export PYTHONUNBUFFERED=1     # stream progress to the log live instead of block-buffering

# Populate the shared spins + gradient (AP/C1/C2/C3) nulls in ../outputs/spin_cache.
# Running on ASD also writes ASD's correlation + occipital results (harmless; the
# ASD array task reloads them from cache).
python run_spin_test.py \
    --n-spin 5000 --seed 0 \
    --asd-map "$OUT/03_01_ASD_regional_enrichment.csv" \
    --asd-cols z_vs_bg \
    --roi-mode union --roi-name occipital \
    --prefix 03_01_ASD_ --outdir "$OUT" \
    --cachedir "$OUT/spin_cache"

echo "prep done: shared spins + gradient nulls cached in $OUT/spin_cache"
