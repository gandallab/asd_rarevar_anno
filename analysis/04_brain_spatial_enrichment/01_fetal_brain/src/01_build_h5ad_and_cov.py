#!/usr/bin/env python
"""
Fig2-occipital-only: stage 1 (src/00-08)
Build the raw-count AnnData for occipital section A1 (FB080-O1) and its scDRS
covariate file.

scDRS standardises norm_score against control gene sets matched within the
dataset it is given, so a score is a statement about the object it was computed
in. Scores from this run should not be set beside scores from any other scoring
run, here or elsewhere in the manuscript.

Covariates are `const` (the intercept) and n_genes, the number of genes detected
per spot; scDRS regresses them out of expression before scoring. With a single
section from a single donor, section, donor and sex are all constant and
therefore unidentifiable, so none of them enters.

Inputs   export/A1_counts.mtx, export/A1_barcodes.tsv, export/A1_features.tsv
Outputs  export/visium_A1_raw.h5ad, export/visium_A1_ngene.cov,
         results/occipital/01_substrate.json
"""
import os, json

NT = "4"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = NT
os.environ["NUMBA_NUM_THREADS"] = NT
os.environ["NUMBA_CACHE_DIR"] = os.getcwd() + "/.numba_cache"

import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from scipy.io import mmread

os.makedirs("results/occipital", exist_ok=True)

# Counts are written as float32 (values unchanged). scDRS 1.0.2 normalises in
# place through a deprecated scanpy routine that cannot handle an integer matrix
# under scanpy 1.11.x.
counts = sp.csr_matrix(mmread("export/A1_counts.mtx").T).astype(np.float32)
gene_ids = pd.read_csv("export/A1_features.tsv", header=None)[0].astype(str)
spot_barcodes = pd.read_csv("export/A1_barcodes.tsv", header=None)[0].astype(str)
adata = ad.AnnData(counts, obs=pd.DataFrame(index=spot_barcodes),
                   var=pd.DataFrame(index=gene_ids))
adata.obs.index.name = None
adata.var.index.name = None
# Barcodes carry the section prefix (A1_) so that spots are identified the same
# way throughout. Visium barcodes come from a shared
# whitelist, so the same barcode can occur in both sections.
adata.obs_names = [f"A1_{b}" for b in adata.obs_names]

assert adata.n_obs == 3591, (
    f"section A1 should have 3,591 spots, found {adata.n_obs}")
assert adata.n_vars == 18085, (
    f"expected the 18,085-gene Visium feature space, found {adata.n_vars}")
assert adata.obs_names.is_unique, "spot keys must be unique within the section"
assert adata.X.min() >= 0 and not np.isnan(adata.X.sum()), (
    "counts must be non-negative and contain no NaN")
frac_int = float(np.mean(adata.X.data == np.round(adata.X.data)))
assert frac_int == 1.0, (
    "scDRS is run with --flag-raw-count True and needs raw integer counts; "
    f"only {frac_int:.4f} of the non-zero entries are integral")

# Genes detected per spot, counted on the full matrix before scDRS applies its
# own spot and gene filters.
n_genes = np.asarray((adata.X > 0).sum(1)).ravel().astype(int)
covariates = pd.DataFrame({"const": 1, "n_genes": n_genes}, index=adata.obs_names)
covariates.to_csv("export/visium_A1_ngene.cov", sep="\t", index=True,
                  index_label="index")
adata.write_h5ad("export/visium_A1_raw.h5ad")

summary = {
    "section": "A1 / FB080-O1 (occipital)",
    "n_spots": int(adata.n_obs),
    "n_genes_input": int(adata.n_vars),
    "n_genes_per_spot_median": int(np.median(n_genes)),
    "n_spots_below_scdrs_min_genes_250": int((n_genes < 250).sum()),
    "covariates": list(covariates.columns),
    "section_covariate": "dropped - degenerate with a single section",
    "comparability": "norm_score is standardised within this object and is not comparable across scoring runs",
}
with open("results/occipital/01_substrate.json", "w") as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
