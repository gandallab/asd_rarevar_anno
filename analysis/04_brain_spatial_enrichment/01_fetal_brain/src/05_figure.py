#!/usr/bin/env python
"""
Fig2-occipital-only: stage 5 (src/00-08)
Group table and the seven-panel manuscript figure: six spot panels (a-f) in two
rows of three, a colourbar lane, then the layer panel (g) with the shared layer
legend beneath it.

Panel g plots assoc_mcz for each group, the manuscript's "Enrichment Z-score",
with markers filled where assoc_fdr < 0.05. All markers share one size (s = 26);
group sizes are in the n_cell column of the group table. Zero, the expectation
under the matched control gene sets, is drawn as a dashed line. There is no null
band because assoc_mcz is already standardised against that null. Colours,
display names (L2/3, V1+V2) and marker encoding follow
src/reference/yundan_panel_i_assoc_mcz.py.reference, so this panel is drawn on the
same convention as the other scDRS panels in the manuscript. (In the regional
analyses "enrichment Z" is a different, resampling-based statistic; see
STATISTIC_GLOSSARY_asd253.md.)

assoc_mcz is plotted rather than the group mean because scDRS tests each group on
the 95th percentile (q95) of norm_score. A plotted mean and a filled or hollow
call would describe two different statistics; assoc_mcz is the standardised form
of the tested q95, so a marker's position and its call agree. The group mean is
still tabulated in layer_group_stats.csv (mean_norm_score, with a 95% spatial
block-bootstrap interval in ci_lo and ci_hi, and n_blocks) but is not drawn,
because the bootstrap interval is the wrong yardstick for the test. It resamples
spots with the gene set held fixed, whereas the scDRS null holds the spots fixed
and varies the gene set. The control gene-set null SD of the group mean exceeds
the bootstrap SE in all 12 groups (2.0-14.4x, median 5.6x), so an interval that
excludes zero does not imply significance. Stage 8 plots this comparison.

Panel a labels the four cortical-plate bands in situ, once per area. The germinal
bands are shared by V1 and V2 and carry no areal label. Every panel is section
A1: the figure contains no prefrontal data, and the areal series are V1, V2 and
their union in the shared germinal zones.

assoc_mcp is one-sided (upper tail), so it can show enrichment but not
depletion. A proliferative zone with a large negative assoc_mcz is not enriched;
that is not evidence that it is depleted.

Inputs   results/occipital/downstream/ASC_Pfdr.scdrs_group.panel_group,
         results/occipital/{per_spot_scores.parquet,spot_layer_assignment.csv,
         umap_A1.parquet}
Outputs  results/occipital/layer_group_stats.csv,
         outputs/figures/Fig2/Fig2-occipital-only.{png,pdf,svg},
         results/occipital/05_figure_qc.json
"""
import os, json

os.environ["NUMBA_CACHE_DIR"] = os.getcwd() + "/.numba_cache"

import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

LAYER_ORDER = ["VZ", "iSVZ", "oSVZ", "IZ", "SP", "L5/6", "L4", "upper plate"]
PROLIFERATIVE = ["VZ", "iSVZ", "oSVZ"]          # the proliferative germinal zones
# Germinal groups are one series, shared by V1 and V2.
SERIES = {"occipital": "occipital (not split)", "V1": "V1", "V2": "V2"}
N_BOOT, TILE_SPACINGS, BOOT_SEED = 2000, 5, 77
FDR_LEVEL = 0.05

spots = pd.read_csv("results/occipital/spot_layer_assignment.csv")
per_spot = pd.read_parquet("results/occipital/per_spot_scores.parquet")
umap_table = (pd.read_parquet("results/occipital/umap_A1.parquet")
              .merge(per_spot[["barcode", "norm_score", "panel_group"]],
                     on="barcode", validate="1:1")
              .merge(spots[["barcode", "layer", "area"]], on="barcode", validate="1:1"))
scdrs_groups = pd.read_csv(
    "results/occipital/downstream/ASC_Pfdr.scdrs_group.panel_group",
    sep="\t", index_col=0)
scdrs_groups = scdrs_groups[scdrs_groups.index != "excluded"].copy()
assert len(scdrs_groups) == 12, (
    "the group table should hold the 12 occipital panel groups once the "
    f"'excluded' row is removed, found {len(scdrs_groups)}")

# Spots with both a layer label and a score. These include the opposing wall,
# which is drawn in the spatial panels but not tested.
spots_shown = spots[spots.layer.notna() & spots.scored].copy()
assert len(umap_table) == int(spots.scored.sum()) == 3589, (
    "the embedding must hold exactly the 3,589 spots scDRS scored; found "
    f"{len(umap_table)} in the embedding and {int(spots.scored.sum())} scored")

# Spatial block bootstrap of the group mean (kept in the table, not plotted). The
# section is tiled into square blocks ~5 nearest-neighbour spot spacings wide;
# within a group, the spots of one tile form a block. Blocks are resampled with
# replacement (as many as the group has, 2,000 times, seed 77) and the group mean
# norm_score is recomputed each time. The 2.5th and 97.5th percentiles of those
# means give a 95% interval that allows for short-range spatial correlation,
# because neighbouring spots are resampled together.
def tile_id(xy, k=TILE_SPACINGS):
    nn_dist, _ = cKDTree(xy).query(xy, k=2)
    # Tile width: k times the median distance between nearest-neighbour spots.
    tile_size = k * float(np.median(nn_dist[:, 1]))
    return ((xy[:, 0] // tile_size).astype(np.int64) * 100_000
            + (xy[:, 1] // tile_size).astype(np.int64))


grouped_spots = spots[(spots.panel_group != "excluded") & spots.scored].copy()
grouped_spots["tile"] = tile_id(grouped_spots[["imagecol", "imagerow"]].to_numpy(float))
rng = np.random.default_rng(BOOT_SEED)
rows = []
for group in scdrs_groups.index:
    members = grouped_spots[grouped_spots.panel_group == group]
    scores = members.norm_score.to_numpy(float)
    blocks = [scores[members.tile.to_numpy() == tile]
              for tile in np.unique(members.tile)]
    draws = rng.integers(0, len(blocks), size=(N_BOOT, len(blocks)))
    boot_means = np.array([np.concatenate([blocks[j] for j in draw]).mean()
                           for draw in draws])
    # obs_q95 is recomputed from the per-spot scores; the scDRS group file reports
    # assoc_mcz and assoc_mcp but not the observed q95 itself.
    rows.append({"panel_group": group,
                 "n_cell": int(scdrs_groups.loc[group, "n_cell"]),
                 "n_blocks": len(blocks),
                 "mean_norm_score": float(scores.mean()),
                 "ci_lo": float(np.quantile(boot_means, .025)),
                 "ci_hi": float(np.quantile(boot_means, .975)),
                 "obs_q95": float(np.quantile(scores, .95)),
                 "assoc_mcz": float(scdrs_groups.loc[group, "assoc_mcz"]),
                 "assoc_mcp": float(scdrs_groups.loc[group, "assoc_mcp"]),
                 "hetero_mcp": float(scdrs_groups.loc[group, "hetero_mcp"]),
                 "hetero_mcz": float(scdrs_groups.loc[group, "hetero_mcz"])})
group_stats = pd.DataFrame(rows)
group_stats["series"] = [SERIES[group.split("|")[0]]
                         for group in group_stats.panel_group]
group_stats["layer"] = pd.Categorical(
    [group.split("|")[1] for group in group_stats.panel_group],
    LAYER_ORDER, ordered=True)
# Benjamini-Hochberg over the 12 groups, applied to scDRS's one-sided Monte-Carlo p.
group_stats["assoc_fdr"] = multipletests(group_stats.assoc_mcp.values,
                                         method="fdr_bh")[1]
group_stats = group_stats.sort_values(["layer", "series"]).reset_index(drop=True)
group_stats = group_stats[["panel_group", "series", "layer", "n_cell", "n_blocks",
                           "mean_norm_score", "ci_lo", "ci_hi", "obs_q95",
                           "assoc_mcz", "assoc_mcp", "assoc_fdr",
                           "hetero_mcz", "hetero_mcp"]]
group_stats.to_csv("results/occipital/layer_group_stats.csv", index=False)

# One viridis colour per layer, from the ventricular zone to the upper plate.
LAYER_COLOUR = dict(zip(
    LAYER_ORDER,
    plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(LAYER_ORDER)))))
# Colours and display names follow the collaborator's panel script, so one colour
# stands for one area across the manuscript's scDRS panels.
SERIES_COLOUR = {"V1": "#6B3A2A", "V2": "#D4908F", "occipital (not split)": "#888888"}
# Vertical offsets within a layer row: V1 above V2, the shared germinal series centred.
SERIES_OFFSET = {"V1": 0.16, "V2": -0.16, "occipital (not split)": 0.0}
LABEL_MAP = {"upper plate": "L2/3"}
SERIES_MAP = {"occipital (not split)": "V1+V2"}
DISPLAY_LAYERS = [LABEL_MAP.get(layer, layer) for layer in LAYER_ORDER]
# Diverging scale centred on zero, the control-set expectation. The limits are the
# 1st and 99th percentiles of norm_score over all scored spots, so that a few
# extreme spots do not compress the scale.
colour_norm = TwoSlopeNorm(vcenter=0.0, **dict(zip(
    ["vmin", "vmax"], np.percentile(umap_table.norm_score, [1, 99]))))
mpl.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "savefig.dpi": 300})


def bare(a):
    a.set_xticks([])
    a.set_yticks([])
    for s in a.spines.values():
        s.set_visible(False)


def letter(a, s, dx):
    a.text(dx, 1.06, s, transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom")


fig = plt.figure(figsize=(15.6, 6.6))
grid = GridSpec(2, 9, figure=fig, height_ratios=[1.0, 1.0],
                width_ratios=[1, 1, 1, 1, 1, 1, 0.34, 1, 1], hspace=0.40, wspace=0.26)

# The spatial panels show a square window of VIEW_SPAN low-resolution image
# pixels, centred on the tissue.
VIEW_SPAN, SPOT_MARKER_SIZE = 300.0, 2.55
cx = (spots_shown.imagecol.min() + spots_shown.imagecol.max()) / 2
cy = (spots_shown.imagerow.min() + spots_shown.imagerow.max()) / 2


def square(a):
    a.set_xlim(cx - VIEW_SPAN / 2, cx + VIEW_SPAN / 2)
    a.set_ylim(cy - VIEW_SPAN / 2, cy + VIEW_SPAN / 2)
    a.set_aspect("equal", adjustable="box")


ax = {}
ax["a"] = fig.add_subplot(grid[0, 0:2])
# y is imagerow with the axis not inverted: imagerow increases toward the pia,
# which therefore lies at the top of the panel.
for layer in LAYER_ORDER:
    in_layer = spots_shown[spots_shown.layer == layer]
    ax["a"].scatter(in_layer.imagecol, in_layer.imagerow, c=[LAYER_COLOUR[layer]],
                    s=SPOT_MARKER_SIZE, linewidths=0)
ax["a"].set_title("Occipital layer annotation, authors' labels\n"
                  "A1 / FB080-O1  (pia at top)", loc="left", fontsize=7)
for area in ("V2", "V1"):
    in_area = spots_shown[spots_shown.area == area]
    if len(in_area):
        ax["a"].annotate(area, (in_area.imagecol.median(), ax["a"].get_ylim()[1] - 12),
                         ha="center", fontsize=6.5)

# Laminar labels sit on the tissue, once per area, for the four cortical-plate
# bands only; the authors' annotation does not split the germinal zones by area.
# Each label is placed at the spot of its band whose x position lies at the
# quantile given in STAGGER. The quantiles step toward the middle of the section
# with depth, so no two labels share an x position and they follow the oblique
# bands in radial order.
PLATE_LAYERS = ["SP", "L5/6", "L4", "upper plate"]
STAGGER = {"V2": {"upper plate": 0.08, "L4": 0.20, "L5/6": 0.32, "SP": 0.44},
           "V1": {"upper plate": 0.92, "L4": 0.80, "L5/6": 0.68, "SP": 0.56}}
plate_label_artists = []
for area in ("V2", "V1"):
    for layer in PLATE_LAYERS:
        band = spots_shown[(spots_shown.layer == layer) & (spots_shown.area == area)]
        if len(band) < 25:
            continue
        cxs, cys = band.imagecol.to_numpy(), band.imagerow.to_numpy()
        j = int(np.argmin(np.abs(cxs - np.quantile(cxs, STAGGER[area][layer]))))
        plate_label_artists.append(ax["a"].annotate(
            LABEL_MAP.get(layer, layer), (cxs[j], cys[j]), fontsize=4.0, ha="center",
            va="center", color="0.08",
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.80)))

ax["b"] = fig.add_subplot(grid[0, 2:4])
score_points = ax["b"].scatter(spots_shown.imagecol, spots_shown.imagerow,
                               c=spots_shown.norm_score, cmap="RdBu_r",
                               norm=colour_norm, s=SPOT_MARKER_SIZE, linewidths=0)
ax["b"].set_title("Occipital ASD-253 per-spot scDRS Z-score", loc="left", fontsize=7)

# All six spot panels draw the same 3,461 layer-annotated spots. The embedding was
# computed on all 3,589 scored spots, but the 128 spots of the authors' dropped
# cluster 14 have no layer; they are omitted so the layer and score panels show
# the same tissue. The counts are saved in 14_figure_qc.json.
umap_shown = umap_table[umap_table.layer.notna()].copy()

ax["c"] = fig.add_subplot(grid[0, 4:6])
for layer in LAYER_ORDER:
    in_layer = umap_shown[umap_shown.layer == layer]
    ax["c"].scatter(in_layer.u1, in_layer.u2, c=[LAYER_COLOUR[layer]], s=1.6,
                    linewidths=0)
ax["c"].set_title("Occipital-only UMAP\nby layer (no integration step)",
                  loc="left", fontsize=7)

ax["d"] = fig.add_subplot(grid[1, 0:2])
ax["d"].scatter(umap_shown.u1, umap_shown.u2, c=umap_shown.norm_score, cmap="RdBu_r",
                norm=colour_norm, s=1.6, linewidths=0)
ax["d"].set_title("same UMAP\nby ASD-253 Z-score", loc="left", fontsize=7)

# Panels e and f highlight the V1 and V2 spots on the same UMAP. The grey backdrop
# is every other spot, including the germinal-zone spots: the authors' annotation
# shares those between V1 and V2, so they have no area label and are grey in both
# panels. The backdrop is drawn as the complement of the highlighted set, so no
# spot is plotted twice.
for k, area in [("e", "V1"), ("f", "V2")]:
    ax[k] = fig.add_subplot(grid[1, (2 if k == "e" else 4):(4 if k == "e" else 6)])
    highlighted = umap_shown.area == area
    ax[k].scatter(umap_shown.u1[~highlighted], umap_shown.u2[~highlighted],
                  c="0.89", s=1.2, linewidths=0)
    ax[k].scatter(umap_shown.u1[highlighted], umap_shown.u2[highlighted],
                  c=SERIES_COLOUR[area], s=1.6, linewidths=0)
    ax[k].set_title(f"same UMAP\n{area} spots  ({int(highlighted.sum()):,})",
                    loc="left", fontsize=7, color=SERIES_COLOUR[area])

for k in "ab":
    square(ax[k])
    bare(ax[k])
for k in "cdef":
    bare(ax[k])

# Panel g: assoc_mcz by layer.
ax["g"] = fig.add_subplot(grid[0, 7:9])
# Layer rows run from the ventricular zone (bottom) to the pia (top).
layer_row = {layer: i for i, layer in enumerate(LAYER_ORDER)}
for ser, col in SERIES_COLOUR.items():
    series_rows = group_stats[group_stats.series == ser]
    if not len(series_rows):
        continue
    y = series_rows.layer.map(layer_row).to_numpy(float) + SERIES_OFFSET[ser]
    mcz = series_rows.assoc_mcz.to_numpy()
    # One marker size for every group. Group size (n_cell in layer_group_stats.csv,
    # 119-698 spots) is tabulated rather than drawn: assoc_mcz is a z against each
    # group's own control-set null, whose width depends on the number of spots, so
    # its magnitude is not comparable between groups of unequal size.
    marker_size = np.full(len(series_rows), 26.0)
    significant = (series_rows.assoc_fdr < FDR_LEVEL).to_numpy()
    legend_label = SERIES_MAP.get(ser, ser)
    ax["g"].scatter(mcz[significant], y[significant], s=marker_size[significant],
                    facecolor=col, edgecolor="0.15", linewidths=0.5, zorder=3,
                    label=legend_label)
    ax["g"].scatter(mcz[~significant], y[~significant], s=marker_size[~significant],
                    facecolor="white", edgecolor=col, linewidths=0.9, zorder=3,
                    label=None if significant.any() else legend_label)
ax["g"].axvline(0.0, color="0.55", lw=0.7, ls="--", zorder=0)
# Dotted line between the germinal zones (below) and the cortical plate (above).
ax["g"].axhline(3.5, color="0.75", lw=0.7, ls=":", zorder=0)
ax["g"].set_yticks(range(len(LAYER_ORDER)))
ax["g"].set_yticklabels(DISPLAY_LAYERS, fontsize=6.6)
ax["g"].yaxis.tick_right()
ax["g"].yaxis.set_label_position("right")
ax["g"].set_ylabel("layer vocabulary, V1-specific L2/L3 split collapsed\n"
                   "ventricle \u2192 pia", fontsize=6.6)
ax["g"].set_xlabel("Enrichment Z-score (scDRS assoc_mcz)", fontsize=6.6)
for s in ("top", "right"):
    ax["g"].spines[s].set_visible(False)

n_sig = int((group_stats.assoc_fdr < FDR_LEVEL).sum())
proliferative_rows = group_stats[group_stats.layer.isin(PROLIFERATIVE)]
n_prolif_sig = int(((proliferative_rows.assoc_fdr < FDR_LEVEL)
                    & (proliferative_rows.assoc_mcz > 0)).sum())
assert n_prolif_sig == 0, (
    "the panel title states that no proliferative zone is enriched, but:\n"
    f"{proliferative_rows[proliferative_rows.assoc_fdr < FDR_LEVEL]}")
ax["g"].set_title(f"ASD-253 by layer: {n_sig} of {len(group_stats)} groups enriched\n"
                  f"(FDR<{FDR_LEVEL:g}); no proliferative zone enriched",
                  loc="left", fontsize=7)
xlo = group_stats.assoc_mcz.min() - 0.45
xhi = group_stats.assoc_mcz.max() + 1.15          # right margin reserved for the legend
ax["g"].set_xlim(xlo, xhi)
ax["g"].set_ylim(-0.6, len(LAYER_ORDER) + 0.9)      # headroom reserved for the legend
handles, labels = ax["g"].get_legend_handles_labels()
handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                      markeredgecolor="0.35", markersize=5))
labels.append(f"hollow: FDR\u2265{FDR_LEVEL:g}")
ax["g"].legend(handles, labels, frameon=False, fontsize=6.0, loc="upper left", ncol=2,
               handletextpad=0.25, columnspacing=0.9, labelspacing=0.35)

fig.canvas.draw()
fig_width_in = fig.get_size_inches()[0]
for k, panel_ax in ax.items():
    letter(panel_ax, k, dx=-0.14 / (panel_ax.get_position().width * fig_width_in))
# Shared layer legend, beneath panel g.
legend_slot = grid[1, 7:9].get_position(fig)
fig.legend(handles=[Patch(facecolor=LAYER_COLOUR[layer],
                          label=LABEL_MAP.get(layer, layer))
                    for layer in LAYER_ORDER[::-1]],
           frameon=False, fontsize=6.2, loc="upper left",
           bbox_to_anchor=(legend_slot.x0, legend_slot.y1),
           bbox_transform=fig.transFigure,
           ncol=2, columnspacing=1.1,
           handlelength=0.95, handleheight=0.95, labelspacing=0.44, handletextpad=0.42,
           title="layer", title_fontsize=6.4)

fig.canvas.draw()
# Colourbar for panels b and d, in the lane between the spot panels and panel g.
lane = grid[:, 6].get_position(fig)
pos_b, pos_d = ax["b"].get_position(), ax["d"].get_position()
cax = fig.add_axes([lane.x0 + 0.55 * lane.width, pos_d.y0, 0.010, pos_b.y1 - pos_d.y0])
cb = fig.colorbar(score_points, cax=cax)
cax.yaxis.set_ticks_position("left")
cax.yaxis.set_label_position("left")
cb.set_label("ASD-253 scDRS norm_score (per-spot z vs matched-control null)", fontsize=6.0)
cb.ax.tick_params(labelsize=5.5)
cb.outline.set_visible(False)
cb.solids.set_rasterized(False)

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../.."))
FIG_DIR = os.path.join(PROJ, "outputs/figures/Fig2")
os.makedirs(FIG_DIR, exist_ok=True)
fig.savefig(os.path.join(FIG_DIR, "Fig2-occipital-only.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG_DIR, "Fig2-occipital-only.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIG_DIR, "Fig2-occipital-only.svg"), bbox_inches="tight")

# Checks on the rendered figure.
rr = fig.canvas.get_renderer()
rasterized_flags = {k: sorted({bool(c.get_rasterized()) for c in a.collections})
                    for k, a in ax.items()}
assert all(v == [False] for v in rasterized_flags.values()), (
    f"spot markers must stay vector objects in the PDF and SVG: {rasterized_flags}")
from matplotlib.collections import PathCollection
markers_per_panel = {k: sum(len(c.get_offsets()) for c in a.collections
                            if isinstance(c, PathCollection)) for k, a in ax.items()}
assert markers_per_panel["a"] == markers_per_panel["b"] == len(spots_shown), (
    "panels a and b must each draw every scored, layer-annotated spot: "
    f"{markers_per_panel} against {len(spots_shown)} spots")
assert markers_per_panel["c"] == markers_per_panel["d"] == len(umap_shown), (
    "panels c and d must each draw every scored, layer-annotated spot: "
    f"{markers_per_panel} against {len(umap_shown)} spots")
assert markers_per_panel["e"] == markers_per_panel["f"] == len(umap_shown), (
    "panels e and f must draw the highlighted spots plus the grey complement, "
    f"which together make up the spots of c and d: {markers_per_panel} "
    f"against {len(umap_shown)} spots")
assert markers_per_panel["g"] == len(group_stats), (
    f"panel g must carry one marker per group: {markers_per_panel['g']} drawn "
    f"for {len(group_stats)} groups")
assert len(spots_shown) == len(umap_shown), (
    "the spatial and UMAP panels must show the same spots: "
    f"{len(spots_shown)} against {len(umap_shown)}")
geo = {k: (round(ax[k].get_window_extent(rr).width, 2),
           round(ax[k].get_window_extent(rr).height, 2)) for k in "ab"}
assert len(set(geo.values())) == 1, (
    f"panels a and b must be the same size so the maps compare directly: {geo}")
clash = []
for k, a in ax.items():
    text_artists = [t for t in ([a.title] + list(a.texts)) if t.get_text().strip()]
    text_boxes = [(t, t.get_window_extent(rr)) for t in text_artists]
    clash += [(k, x.get_text()[:28], y.get_text()[:28])
              for i, (x, bx) in enumerate(text_boxes)
              for y, by in text_boxes[i + 1:] if bx.overlaps(by)]
assert not clash, f"titles and annotations overlap within a panel: {clash}"
colourbar_box = cax.get_tightbbox(rr)
for panel, axis in ax.items():
    assert not colourbar_box.overlaps(axis.get_tightbbox(rr)), (
        f"the colourbar must not overlap any panel; it overlaps panel {panel}")
legend_box = ax["g"].get_legend().get_window_extent(rr)
(lx0, ly0), (lx1, ly1) = ax["g"].transData.inverted().transform(
    [[legend_box.x0, legend_box.y0], [legend_box.x1, legend_box.y1]])
# The legend must not cover any plotted marker.
covered = []
for ser in SERIES_COLOUR:
    series_rows = group_stats[group_stats.series == ser]
    if not len(series_rows):
        continue
    y_marker = series_rows.layer.map(layer_row).to_numpy(float) + SERIES_OFFSET[ser]
    in_legend_rows = (y_marker >= ly0) & (y_marker <= ly1)
    mcz = series_rows.assoc_mcz.to_numpy()
    covered += [f"{group} marker" for group, ok in
                zip(series_rows.panel_group,
                    in_legend_rows & (mcz >= lx0) & (mcz <= lx1)) if ok]
assert not covered, f"legend covers plotted assoc_mcz markers: {covered}"

# The laminar labels of panel a are placed from the data and must not overlap.
label_boxes = [(t.get_text(), t.get_window_extent(rr)) for t in plate_label_artists]
label_clashes = [(x, y) for i, (x, bx) in enumerate(label_boxes)
                 for y, by in label_boxes[i + 1:] if bx.overlaps(by)]
assert not label_clashes, f"laminar labels in panel a overlap: {label_clashes}"
assert len(plate_label_artists) == 8, (
    "expected one laminar label for each of the four cortical-plate bands in V1 "
    "and V2 (a band needs at least 25 spots to be labelled), "
    f"got {len(plate_label_artists)}")

qc = {
    "n_spots_scored": int(spots.scored.sum()),
    "n_spots_displayed_spatial": len(spots_shown),
    "n_spots_in_umap": len(umap_table),
    "n_v1_highlighted": int((umap_shown.area == "V1").sum()),
    "n_v2_highlighted": int((umap_shown.area == "V2").sum()),
    "n_germinal_backdrop_both_panels": int(umap_shown.area.isna().sum()),
    "markers_per_panel": {k: int(v) for k, v in markers_per_panel.items()},
    "total_vector_markers": int(sum(markers_per_panel.values())),
    "all_panels_vector": all(v == [False] for v in rasterized_flags.values()),
    "text_collisions": len(clash),
    "panel_ab_pixel_size": geo["a"],
    "colour_scale_vmin_vmax": [float(colour_norm.vmin), float(colour_norm.vmax)],
    "n_groups": len(group_stats),
    "n_groups_enriched_fdr05": n_sig,
    "n_proliferative_enriched": n_prolif_sig,
    "bootstrap": {"B": N_BOOT, "tile_k_hex_spacings": TILE_SPACINGS, "seed": BOOT_SEED,
                  "n_blocks_range": [int(group_stats.n_blocks.min()),
                                     int(group_stats.n_blocks.max())]},
}
with open("results/occipital/05_figure_qc.json", "w") as fh:
    json.dump(qc, fh, indent=1)
print(json.dumps(qc, indent=1))
print()
print(group_stats.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
