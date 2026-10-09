#!/usr/bin/env python
"""
Fig2-occipital-only: stage 3 (src/00-08)
Attach the authors' layer annotation to the scored occipital spots, define the
panel groups, and compute a single-section embedding.

Layer labels are the authors' own and are used directly. Nothing is transferred
between sections or inferred from clustering, so the layer assignment carries no
label-transfer uncertainty.

A single section has no batch to correct, so there is no integration step: the
embedding is a plain normalise / HVG / PCA / UMAP. It is for display only, and
no statistic in the figure depends on it.

Panel groups (12, held in one composite column):
  occipital|<layer>        VZ, iSVZ, oSVZ, IZ: germinal zones, not split by area
  V1|<layer>, V2|<layer>   SP, L5/6, L4, upper plate: cortical plate, split by area
Every other spot is labelled "excluded": the opposing cortical wall, the authors'
dropped cluster 14, and any spot removed by the scDRS cell filter.

Inputs   results/occipital/score_file/253_genes/ASC_Pfdr.full_score.gz,
         export/A1_author_labels.csv, export/A1_coords_scales.parquet,
         export/visium_A1_raw.h5ad
Outputs  export/visium_A1_groups.h5ad, results/occipital/per_spot_scores.parquet,
         results/occipital/spot_layer_assignment.csv,
         results/occipital/umap_A1.parquet, results/occipital/03_assignment.json
"""
import os, json

NT = "4"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = NT
os.environ["NUMBA_NUM_THREADS"] = NT
os.environ["NUMBA_CACHE_DIR"] = os.getcwd() + "/.numba_cache"

import numpy as np, pandas as pd
import anndata as ad, scanpy as sc

SCORE_FILE = "results/occipital/score_file/253_genes/ASC_Pfdr.full_score.gz"
# The authors annotate V1's upper plate as separate L2 and L3 bands, V2's as a
# single L2/3. The three labels are merged so that both areas share one layer
# vocabulary. The merged band is keyed "upper plate" in the tables and shown as
# L2/3 in the figure.
LAYER_COLLAPSE = {"L2": "upper plate", "L3": "upper plate", "L2/3": "upper plate"}
GERMINAL, PLATE = ["VZ", "iSVZ", "oSVZ", "IZ"], ["SP", "L5/6", "L4", "upper plate"]

# scDRS output. Spots removed by its cell filter are absent.
full_score = pd.read_csv(SCORE_FILE, sep="\t", index_col=0)
per_spot = pd.DataFrame({
    "key": full_score.index.astype(str),
    "norm_score": full_score["norm_score"].to_numpy(np.float64),
    "raw_score": full_score["raw_score"].to_numpy(np.float64),
    "mc_pval": full_score["pval"].to_numpy(np.float64),
})
per_spot["barcode"] = [key.split("_", 1)[1] for key in per_spot.key]
del full_score

# The authors' labels carry layer, area and wall in one string (L4-V1, oSVZ-2,
# iSVZ-L). Removing the area suffix (-V1, -V2) and any cluster number or wall
# suffix (oSVZ-1 to oSVZ-3; -L, -L1, -L2) leaves the layer.
author_labels = pd.read_csv("export/A1_author_labels.csv")
author_labels["layer_a"] = (
    author_labels.author_label.str.replace(r"-V[12]$", "", regex=True)
                              .str.replace(r"-(L?\d|L)$", "", regex=True))
# A trailing -L (oSVZ-L1, oSVZ-L2, iSVZ-L) marks the opposing cortical wall.
author_labels["wall"] = np.where(
    author_labels.author_label.str.contains(r"-L\d?$|^iSVZ-L$", regex=True),
    "opposing", "main")
# Only cortical-plate labels carry an area; germinal labels have none.
author_labels["area"] = np.where(
    author_labels.author_label.str.endswith("-V1"), "V1",
    np.where(author_labels.author_label.str.endswith("-V2"), "V2", None))
author_labels["layer"] = author_labels.layer_a.map(
    lambda layer: LAYER_COLLAPSE.get(layer, layer))
# Cluster 14 was dropped by the authors and carries neither layer nor area.
author_labels.loc[author_labels.author_label == "14", ["layer", "area"]] = None

coords = (pd.read_parquet("export/A1_coords_scales.parquet")
            .rename(columns={"row_raw": "imagerow", "col_raw": "imagecol"}))
spots = coords.merge(
    author_labels[["barcode", "author_label", "layer", "wall", "area"]],
    on="barcode", validate="1:1")
spots["key"] = "A1_" + spots.barcode
assert len(spots) == 3591 and spots.key.is_unique, (
    "all 3,591 section-A1 spots need one coordinate row and one author label, "
    f"with unique keys; got {len(spots)} spots")


# Germinal zones are pooled across V1 and V2, because the annotation does not
# split them. Cortical-plate bands are split by area.
def panel_group(spot):
    if spot.wall != "main" or pd.isna(spot.layer):
        return "excluded"
    if spot.layer in GERMINAL:
        return f"occipital|{spot.layer}"
    if spot.layer in PLATE and spot.area in ("V1", "V2"):
        return f"{spot.area}|{spot.layer}"
    return "excluded"


spots["panel_group"] = spots.apply(panel_group, axis=1)
spots = spots.merge(per_spot.drop(columns="key"), on="barcode", how="left",
                    validate="1:1")
# A spot missing from the scDRS output was removed by its cell filter. It keeps
# its row here but belongs to no group.
spots["scored"] = spots.norm_score.notna()
spots.loc[~spots.scored, "panel_group"] = "excluded"
spots.to_csv("results/occipital/spot_layer_assignment.csv", index=False)
spots[spots.scored].drop(columns=["imagerow", "imagecol"]).to_parquet(
    "results/occipital/per_spot_scores.parquet", index=False)

# Single-section embedding, for display only. It is built from exactly the spots
# and genes scDRS scored (the same min_genes and min_cells filters), so a spot's
# position and its score refer to the same count matrix.
emb = ad.read_h5ad("export/visium_A1_raw.h5ad")
sc.pp.filter_cells(emb, min_genes=250)
sc.pp.filter_genes(emb, min_cells=50)
assert emb.n_obs == int(spots.scored.sum()), (
    "the embedding must contain exactly the spots scDRS scored: the min_genes=250 "
    f"filter left {emb.n_obs}, scDRS scored {int(spots.scored.sum())}")
sc.pp.normalize_total(emb, target_sum=1e4)
sc.pp.log1p(emb)
sc.pp.highly_variable_genes(emb, n_top_genes=2000)
emb = emb[:, emb.var.highly_variable].copy()
sc.pp.scale(emb, max_value=10)
sc.tl.pca(emb, n_comps=30, svd_solver="arpack", random_state=0)
sc.pp.neighbors(emb, n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.umap(emb, random_state=0)
umap_table = pd.DataFrame(emb.obsm["X_umap"], columns=["u1", "u2"],
                          index=emb.obs_names)
umap_table["barcode"] = [key.split("_", 1)[1] for key in umap_table.index]
umap_table.reset_index(drop=True).to_parquet("results/occipital/umap_A1.parquet",
                                             index=False)

# perform-downstream reads each spot's group from this AnnData. Excluded and
# unscored spots carry the label "excluded"; scDRS reports a row for that label
# too, and stage 5 drops it.
grouped = ad.read_h5ad("export/visium_A1_raw.h5ad")
group_labels = spots.set_index("key").loc[grouped.obs_names, "panel_group"]
grouped.obs["panel_group"] = pd.Categorical(group_labels.values)
grouped.write_h5ad("export/visium_A1_groups.h5ad")

group_sizes = spots.panel_group.value_counts().to_dict()
summary = {
    "n_spots_total": len(spots),
    "n_spots_scored": int(spots.scored.sum()),
    "spots_dropped_by_scdrs_cell_filter": spots.loc[~spots.scored, "key"].tolist(),
    "n_groups_in_panel": int(
        spots.loc[spots.panel_group != "excluded", "panel_group"].nunique()),
    "n_spots_excluded_from_panel": int((spots.panel_group == "excluded").sum()),
    "group_sizes": group_sizes,
    "layer_labels": "authors' own annotation, direct - no label transfer",
    "embedding": "scanpy normalise/HVG(2000)/scale/PCA(30)/neighbors(15)/UMAP, no integration",
    "n_spots_in_embedding": len(umap_table),
}
with open("results/occipital/03_assignment.json", "w") as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
