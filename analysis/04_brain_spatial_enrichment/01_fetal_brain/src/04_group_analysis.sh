#!/bin/bash
# =============================================================================
# Fig2-occipital-only: stage 4 (src/00-08)
# scDRS group-level association over the 12 occipital panel groups.
#
# For each group,
# scDRS takes the 95th percentile (q95) of the spots' norm_score and compares it
# with the q95 of each of the 1,000 control gene sets in the same spots.
# assoc_mcz is the z of the observed q95 against the control q95s, and assoc_mcp
# is the one-sided upper Monte-Carlo p, which cannot fall below 1/1001. A
# one-sided test can show enrichment but not depletion: a group with a large
# negative assoc_mcz is simply not enriched. The filter flags match stage 2, so
# both stages work on the same spots.
#
# The norm_scores and the control sets both come from this object, so assoc_mcz
# is standardised within it: the same spots scored in a different object would
# give a different z. The Benjamini-Hochberg correction is applied afterwards,
# in stage 5, over these 12 groups.
# =============================================================================
set -euo pipefail

export NUMBA_CACHE_DIR=$PWD/.numba_cache
NT=4
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
export NUMEXPR_NUM_THREADS=$NT VECLIB_MAXIMUM_THREADS=$NT NUMBA_NUM_THREADS=$NT

mkdir -p results/occipital/downstream logs

scdrs perform-downstream \
    --h5ad-file export/visium_A1_groups.h5ad \
    --score-file "results/occipital/score_file/253_genes/@.full_score.gz" \
    --out-folder results/occipital/downstream \
    --group-analysis panel_group \
    --flag-filter-data True \
    --flag-raw-count True \
    2>&1 | tee logs/04_group_analysis.log

# scDRS writes tqdm progress bars to the same stream as its banner and its
# substrate counts, so the raw log is mostly carriage returns. Keep a cleaned
# copy: later stages read the scored spot and gene counts back out of it rather
# than restating them, and those lines are invisible under the progress output.
python - "logs/04_group_analysis.log" "logs/04_group_analysis.clean.log" <<'PYEOF'
import sys
raw = open(sys.argv[1], errors="replace").read()
lines = [ln.split("\r")[-1] for ln in raw.split("\n")]
keep = [ln for ln in lines if "it/s]" not in ln and "?it/s" not in ln]
open(sys.argv[2], "w").write("\n".join(keep))
PYEOF

echo "stage 4 complete"
