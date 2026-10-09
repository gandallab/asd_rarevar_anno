# ASD-253 rare-variant enrichment across layers of second-trimester human visual cortex (10x Visium, GW20)

This module scores the ASD-253 rare-variant gene set in a 10x Visium section of
human occipital cortex at gestational week 20 and tests, by layer and by
cortical area, whether the set is enriched. The section spans primary and
secondary visual cortex, so every result here is a V1/V2 result. The gene set
(`ASC_TADA_Pfdr.253.gs`) has 253 genes, each weighted by a z-score derived from
the ASC TADA table, and scoring uses the scDRS command-line pipeline that
produced the other scDRS results in the manuscript.

The output is `Fig2-occipital-only` (panels a–g). Panels a and b map the
authors' laminar annotation and the per-spot score onto the tissue; c–f are a
UMAP of the same spots coloured by layer, by score, and with V1 and V2 spots
highlighted; panel g is the layer panel, plotting the manuscript's Enrichment
Z-score (`assoc_mcz`) for each of the 12 layer-by-area groups.

A prefrontal section from the same deposit was scored alongside this one in an
earlier version of the analysis. It has been dropped. The two sections come from
different donors, so no prefrontal-versus-occipital contrast is separable from a
donor difference, and nothing in this module reads it. The V1-versus-V2 contrast
that remains lies within one section and one donor.

Dropping it does not change any published number, and that was checked rather
than assumed. Running stages 0 to 5 of this module from the deposited occipital
Seurat object alone — the prefrontal object never opened — reproduces the group
table exactly: all 14 columns of `results/occipital/layer_group_stats.csv` agree
across all 12 groups, including `assoc_mcz`, `assoc_mcp`, `assoc_fdr`,
`hetero_mcz`, `obs_q95` and the bootstrap bounds. Every one of the 3,589
per-spot `norm_score` values is identical, and the scored substrate is the same
13,908 genes and 242 gene-set genes. The earlier two-section run pooled the
sections for scoring, so this was not a foregone conclusion; `norm_score` is
standardised within the scored object, and pooling changed it.

Start with `METHODS_occipital_only.md`. `STATISTIC_GLOSSARY_asd253.md` sorts out
the three z-valued statistics that the manuscript calls "enrichment Z" in one
form or another, which are not interchangeable. `SOFTWARE_VERSIONS.md` (also as
`.csv`) lists package versions read from the live environments.

## Source data

The section comes from Qian X, Coleman K, Jiang S, et al., "Spatial
transcriptomics reveals human cortical layer and area specification", *Nature*
644, 153–163 (2025), https://doi.org/10.1038/s41586-025-09010-1. The deposited
Seurat object is at Zenodo, https://doi.org/10.5281/zenodo.14422018. Section A1
is GW20 occipital cortex, 3,591 spots.

The layer and area labels are the authors'. The cluster-to-label mapping in
`src/00_prepare_from_seurat.R` (`AUTHOR_MAP`) is reproduced from the
`RenameIdents` call in `Fig4/Visium_A1_Fig4_h_i_j.R` of the authors' analysis
repository,
https://github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-Layer-and-Area-Specification,
and matches that script at commit `0f5fa04451d2ccac63bb01e26184cfc65162cfff`.
The deposited object does not carry the annotation itself. The authors dropped
their cluster 14 (128 spots), and the groups here exclude it.

## Stages

Run them in numerical order, with one exception: stage 7 reads a table written
by stage 8, so stage 8 runs first. The choices that determine the statistics are
explained in the script that makes them and in the Methods.

| stage | script | does |
|---|---|---|
| 0 | `00_prepare_from_seurat.R` | exports counts, coordinates and the authors' layer labels from the Seurat object |
| 1 | `01_build_h5ad_and_cov.py` | raw-count AnnData and the scDRS covariate file |
| 2 | `02_compute_score.sh` | `scdrs compute-score` |
| 3 | `03_layers_and_embedding.py` | attaches the authors' labels, defines the 12 groups, computes the display embedding |
| 4 | `04_group_analysis.sh` | `scdrs perform-downstream --group-analysis` |
| 5 | `05_figure.py` | group table and the seven-panel figure |
| 6 | `06_versions.py` | captures package versions from the live environments |
| 7 | `07_methods.py` | checks the pattern claims against the group table, then writes the Methods |
| 8 | `08_null_vs_bootstrap.py` | supplementary diagnostic: control-gene-set null against a spatial block bootstrap |

Stages 6 to 8 write documentation and a diagnostic; they add no analysis.

## Inputs

None of the data is committed here. A run needs the following.

| what | where it comes from |
|---|---|
| `Visium_A1_brain_011124.rds` | Zenodo 10.5281/zenodo.14422018 (Qian et al.). Stage 0 reads it; point `VISIUM_RDS_DIR` at the directory it is unpacked into. |
| `A1_author_labels.csv` | Written by stage 0 from `AUTHOR_MAP` (see Source data). |
| `ASC_TADA_Pfdr.253.gs` | Built with `scdrs munge-gs --weight zscore` from the ASC TADA table; not redistributed here. |
| `handoff/versions_r.json` | R package versions. Stage 6 captures the Python side itself and reads this for the R side. |
| `export/`, `results/`, `logs/` | Written by the earlier stages. |
| `ASC_Pfdr.full_score.gz` (tens of MB) | Written by stage 2 and read by stages 3, 4 and 8. It holds the per-spot scores and all 1,000 control sets, and is too large to commit. |

## Analysis choices

The gene set is weighted. The collaborator's driver
(`src/reference/run_scdrs_rvas.sh.reference`) calls `munge-gs --weight zscore`,
and the supplied `.gs` file carries 253 distinct weights (6.2039–10.0000), which
the CLI passes to `score_cell` as `gene_weight`. Earlier drafts of the
manuscript Methods said the set "was run unweighted (binary)". The weights are
kept, and the Methods are being amended to match what was run. Rescoring with
uniform weights changes the pattern very little (Spearman correlation of
`assoc_mcz` 0.993, mean |Δ| 0.09) and flips one marginal call, V2|SP; see
`results/occipital/weighted_vs_uniform.csv`. Whether the atlas panels used the
same weighted file still needs confirming before the Visium and atlas numbers
are pooled.

The layer panel plots `assoc_mcz`, the manuscript's "Enrichment Z-score", and
not the group mean. For each group, scDRS compares the 95th percentile of the
per-spot scores with the 95th percentile of each control gene set in the same
spots, and `assoc_mcz` is the z of the observed value against that control
distribution. Its P value, `assoc_mcp`, is one-sided upper: it asks only whether
the group exceeds its controls. A proliferative zone with a large negative
`assoc_mcz` is therefore not enriched, and since the test has no lower tail it
says nothing about depletion. An earlier version of the panel plotted the group
mean with a spatial block-bootstrap confidence interval. That was replaced
because the plotted point and the significance call came from different
statistics, and because the control gene-set null SD of the group mean exceeds
the bootstrap SE in all 12 groups (2.0–14.4×, median 5.6×), so an interval
excluding zero was not equivalent to significance.

`L2/3` is a collapsed band. The authors split V1's upper plate into L2 and L3,
while V2 carries a single L2/3; the three are collapsed into one band so the two
areas share a vocabulary. It is called `upper plate` in the tables and code and
shown as L2/3 in the figure, and the V1-specific separation of L2 from L3 cannot
be recovered downstream.

The covariates are `const` and `n_genes`. With one section and one donor, both a
section indicator and a donor term would be constant and therefore degenerate;
sex is likewise unidentifiable.

## Limitations

- One donor and one section. Nothing here speaks to between-donor variability,
  and the V1-versus-V2 difference is a within-section observation rather than a
  tested areal effect.
- The scDRS version used for the manuscript's other panels cannot be recovered:
  the collaborator's driver activates a virtualenv by path
  (`src/reference/run_scdrs_rvas.sh.reference`, line 13) and records no version.
  These runs used scDRS 1.0.2. The driver does set `--n-ctrl 1000` explicitly
  (line 57), and that is the value used here.
- `--flag-filter-data True` removes two spots below `min_genes=250`, reduces the
  gene universe from 18,085 to 13,908 at `min_cells=50`, and leaves 242 of the
  253 gene-set genes scored. Both thresholds are applied within this section.
- With `--n-ctrl 1000` the smallest attainable P value is 1/1001, so groups at
  that floor are reported as *P* < 0.001 and their FDR values are
  resolution-limited.

## Reproducibility and software

`glmGamPoi` was not installed, so `SCTransform` in stage 0 used sctransform's
own model fitting. This cannot be recovered from the script text.

`umap-learn`, `pynndescent` and `llvmlite` determine the embedding in panels
c–f, and are recorded in `SOFTWARE_VERSIONS.csv` along with the rest of the
stack. UMAP is stochastic and sensitive to the versions of all three, so that
embedding should be expected to shift under a different combination even though
the scores behind the colours will not. Bit-identical reproduction was
established on macOS arm64 only. Scoring is invariant to thread count (4 versus
16) and to `--flag-return-ctrl-raw-score`.

`environment_python.yml` and `environment_r.txt` pin what was used, and
`SOFTWARE_VERSIONS.md` records every component as read from the live
environments at run time. numpy and scipy are held below 2 and 1.13
respectively because scDRS 1.0.2 relies on interfaces removed in later releases.
Count matrices are stored as float32, because scDRS 1.0.2 calls a deprecated
scanpy normalisation helper that divides in place and fails on an integer matrix
under scanpy 1.11.5; the values are unchanged.

`ruff.toml` records which lint rules this module ignores and why.

## Reference scripts

`src/reference/` holds the collaborator's own scripts, unmodified, for
conformance checking: the scDRS driver used for the manuscript's other panels,
and the script that draws the laminar panel and defines its conventions. This
pipeline does not call them.

The collaborator's panel script lists the layers from pia to ventricle and
inverts the y-axis to put the pia at the top. The figure here lists them from
ventricle to pia on an ordinary axis, which gives the same picture, but the
offsets that separate the areas within a layer then stack them in the opposite
order. The offsets from one script cannot be reused in the other without
reversing their sign.
