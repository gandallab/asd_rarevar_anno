# Software versions

Read from the live environments at the time of the run (`importlib.metadata.version` for Python, `packageVersion` for R), not transcribed from script headers.

## Stage 0 - Seurat export (R)

| component | version |
|---|---|
| R | 4.4.3 |
| platform | aarch64-apple-darwin20.0.0 |
| Seurat | 5.5.1 |
| SeuratObject | 5.4.0 |
| sctransform | 0.4.3 |
| Matrix | 1.7.6 |
| irlba | 2.3.7 |
| uwot | 0.2.5 |
| RcppAnnoy | 0.0.23 |
| Rcpp | 1.1.2 |
| arrow | 25.0.0 |
| dplyr | 1.2.1 |
| png | 0.1.9 |
| future | 1.75.0 |
| rlang | 1.3.0 |
| ggplot2 | 4.0.3 |

## Stages 1-5 and 7 - scDRS scoring, group analysis, figures (Python)

| component | version |
|---|---|
| python | 3.11.16 |
| platform | macOS-27.0.1-arm64-arm-64bit |
| scdrs (CLI masthead) | 1.0.2 |
| scdrs | 1.0.2 |
| scanpy | 1.11.5 |
| anndata | 0.12.19 |
| numpy | 1.26.4 |
| scipy | 1.12.0 |
| pandas | 2.3.3 |
| statsmodels | 0.15.0 |
| matplotlib | 3.11.1 |
| numba | 0.67.0 |
| h5py | 3.16.0 |
| pyarrow | 25.0.0 |
| scikit-learn | 1.9.0 |
| fire | 0.7.1 |
| legacy-api-wrap | 1.5 |
| natsort | 8.4.0 |
| tqdm | 4.70.0 |
| umap-learn | 0.5.12 |
| pynndescent | 0.6.0 |
| llvmlite | 0.49.0 |

`glmGamPoi` is not installed. Seurat's `SCTransform` uses it when present, so this run used the native `nb_offset` fit of sctransform instead. The choice of backend changes the variance estimates, and it cannot be recovered from the script text.
