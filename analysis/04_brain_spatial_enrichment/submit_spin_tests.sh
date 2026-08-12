#!/bin/bash
# =============================================================================
# submit_spin_tests.sh  --  submit the shared prep, then the 4-gene-set array
#
# The array (spin_array.sbatch) reuses the spins + gradient nulls that
# spin_prep.sbatch caches, so it must not start until prep succeeds. We chain
# them with SLURM's afterok dependency.
#
#   Usage:  ./submit_spin_tests.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

PREP=$(sbatch --parsable spin_prep.sbatch)
echo "submitted prep  : job $PREP"

ARRAY=$(sbatch --parsable --dependency=afterok:$PREP spin_array.sbatch)
echo "submitted array : job $ARRAY (starts after $PREP succeeds)"
echo
echo "results -> outputs/analysis/04_brain_spatial_enrichment/ as 04_01_ASD_* and 04_02_{SYN,GR,MORPH}_*"
