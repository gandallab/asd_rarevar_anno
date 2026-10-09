#!/usr/bin/env python
"""
Fig2-occipital-only: stage 6 of the pipeline (src/00-08)
Record the software that produced the analysis.

Versions are read from the LIVE environments rather than transcribed, so the
table cannot drift from the environment that actually ran the pipeline. The R
row set is read from the capture written by stage 0's environment; the Python
rows come from importlib.metadata in the interpreter running this script.

The UMAP stack (umap-learn, pynndescent, llvmlite) is listed because it sets the
embedding drawn in panels c-f. UMAP is stochastic and sensitive to all three, so
that embedding should be expected to shift under a different combination even
though the scores behind the colours will not.

Inputs   handoff/versions_{python,r}.json, written by capture_versions() below
Outputs  SOFTWARE_VERSIONS.csv, SOFTWARE_VERSIONS.md
"""
import json, os, platform, re, sys
from importlib import metadata

import pandas as pd

PY_PACKAGES = ["scdrs", "scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib",
               "statsmodels", "numba", "h5py", "pyarrow", "scikit-learn", "fire",
               "legacy-api-wrap", "natsort", "tqdm", "umap-learn", "pynndescent",
               "llvmlite"]


def capture_versions():
    """Read the Python side of the environment and the scDRS banner from the
    scoring log. Earlier versions of this pipeline expected the capture to be
    supplied from outside the deposit, which meant a clean checkout could not
    produce the version table at all."""
    captured = {"python": sys.version.split()[0], "platform": platform.platform(),
                "machine": platform.machine()}
    for package in PY_PACKAGES:
        try:
            captured[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            captured[package] = None
    log = open("logs/02_compute_score.clean.log", errors="replace").read()
    banner = re.search(r"\* Version (\S+)", log)
    captured["scdrs_cli_masthead"] = banner.group(1) if banner else None
    os.makedirs("handoff", exist_ok=True)
    json.dump(captured, open("handoff/versions_python.json", "w"), indent=1)
    return captured


V_PY = capture_versions()
V_R = json.load(open("handoff/versions_r.json"))

# ---- software version table ---------------------------------------------------------
version_rows = []
STAGE_R = "Stage 0 - Seurat export (R)"
STAGE_PY = "Stages 1-5 and 7 - scDRS scoring, group analysis, figures (Python)"
version_rows.append({"stage": STAGE_R, "component": "R", "version": V_R["R"]})
version_rows.append({"stage": STAGE_R, "component": "platform", "version": V_R["platform"]})
for k in ["Seurat", "SeuratObject", "sctransform", "Matrix", "irlba", "uwot",
          "RcppAnnoy", "Rcpp", "arrow", "dplyr", "png", "future", "rlang", "ggplot2"]:
    version_rows.append({"stage": STAGE_R, "component": k,
                         "version": V_R.get(k) if V_R.get(k) else "NOT INSTALLED"})
version_rows.append({"stage": STAGE_PY, "component": "python", "version": V_PY["python"]})
version_rows.append({"stage": STAGE_PY, "component": "platform", "version": V_PY["platform"]})
version_rows.append({"stage": STAGE_PY, "component": "scdrs (CLI masthead)",
                     "version": V_PY["scdrs_cli_masthead"]})
# These three set the UMAP of panels c-f (stage 3). They are listed only when the
# environment capture recorded them.
UMAP_PACKAGES = ["umap-learn", "pynndescent", "llvmlite"]
for k in ["scdrs", "scanpy", "anndata", "numpy", "scipy", "pandas", "statsmodels",
          "matplotlib", "numba", "h5py", "pyarrow", "scikit-learn", "fire",
          "legacy-api-wrap", "natsort", "tqdm"] + [p for p in UMAP_PACKAGES if V_PY.get(p)]:
    version_rows.append({"stage": STAGE_PY, "component": k,
                         "version": V_PY.get(k) if V_PY.get(k) else "NOT INSTALLED"})
version_table = pd.DataFrame(version_rows)
version_table.to_csv("SOFTWARE_VERSIONS.csv", index=False)

version_md = ["# Software versions", "",
              ("Read from the live environments at the time of the run "
               "(`importlib.metadata.version` for Python, `packageVersion` for R), not "
               "transcribed from script headers."), ""]
for stage in [STAGE_R, STAGE_PY]:
    version_md += [f"## {stage}", "", "| component | version |", "|---|---|"]
    version_md += [f"| {r.component} | {r.version} |"
                   for r in version_table[version_table.stage == stage].itertuples()]
    version_md += [""]
version_md += [("`glmGamPoi` is not installed. Seurat's `SCTransform` uses it when "
                "present, so this run used the native `nb_offset` fit of sctransform "
                "instead. The choice of backend changes the variance estimates, and it "
                "cannot be recovered from the script text."), ""]
unrecorded = [f"`{p}`" for p in UMAP_PACKAGES if not V_PY.get(p)]
if unrecorded:
    names = (", ".join(unrecorded[:-1]) + " and " + unrecorded[-1]
             if len(unrecorded) > 1 else unrecorded[0])
    version_md += [(f"The occipital-only UMAP (stage 3) depends on {names}, for which no "
                    f"version was captured in this record, so that embedding cannot be "
                    f"rebuilt from this table alone."), ""]
open("SOFTWARE_VERSIONS.md", "w").write("\n".join(version_md))

