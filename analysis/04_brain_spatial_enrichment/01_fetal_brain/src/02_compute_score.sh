#!/bin/bash
# =============================================================================
# Fig2-occipital-only: stage 2 (src/00-08)
# scDRS scoring of the occipital section on its own.
#
# Every flag matches the reference driver
# (src/reference/run_scdrs_rvas.sh.reference), including
# --flag-return-ctrl-raw-score False, which only omits the control sets' raw
# scores from the output; the normalised control scores that every downstream
# test needs are kept by the next flag.
#
# --n-ctrl 1000 draws 1,000 control gene sets. Each control gene comes from the
# same mean-variance bin as the gene-set gene it stands in for (ctrl_match_opt
# mean_var; 20 mean x 20 variance bins, hardcoded in the CLI). The controls'
# norm_scores are kept in the output (--flag-return-ctrl-norm-score True) because
# the group tests of stage 4 and the null comparison of stage 8 are built on
# them.
#
# Other settings are fixed by the package rather than by this script: weight_opt
# vs, which weights each gene by the inverse of its technical SD and is separate
# from the per-gene TADA weights carried in the .gs file; adj_prop None; and
# random_seed 0 (the CLI never passes a seed).
# =============================================================================
set -euo pipefail

export NUMBA_CACHE_DIR=$PWD/.numba_cache
NT=4
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
export NUMEXPR_NUM_THREADS=$NT VECLIB_MAXIMUM_THREADS=$NT NUMBA_NUM_THREADS=$NT

OUT_FOLDER=results/occipital/score_file/253_genes
mkdir -p "${OUT_FOLDER}" logs

scdrs compute-score \
    --h5ad-file export/visium_A1_raw.h5ad \
    --h5ad-species human \
    --gs-file genesets/ASC_TADA_Pfdr.253.gs \
    --gs-species human \
    --out-folder "${OUT_FOLDER}" \
    --cov-file export/visium_A1_ngene.cov \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    2>&1 | tee logs/02_compute_score.log

# scDRS writes tqdm progress bars to the same stream as its banner and its
# substrate counts, so the raw log is mostly carriage returns. Keep a cleaned
# copy: later stages read the scored spot and gene counts back out of it rather
# than restating them, and those lines are invisible under the progress output.
python - "logs/02_compute_score.log" "logs/02_compute_score.clean.log" <<'PYEOF'
import sys
raw = open(sys.argv[1], errors="replace").read()
lines = [ln.split("\r")[-1] for ln in raw.split("\n")]
keep = [ln for ln in lines if "it/s]" not in ln and "?it/s" not in ln]
open(sys.argv[2], "w").write("\n".join(keep))
PYEOF

echo "stage 2 complete"
