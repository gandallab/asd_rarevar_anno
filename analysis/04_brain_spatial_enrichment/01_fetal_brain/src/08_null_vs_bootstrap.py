#!/usr/bin/env python
"""
Fig2-occipital-only: stage 8 (src/00-08)
Supplementary diagnostic: why a group whose mean scDRS Z-score is below 1 can
still be enriched, and why the bootstrap interval is the wrong yardstick for the
enrichment call.

For each of the 12 groups the figure shows the observed mean norm_score, its 95%
spatial block-bootstrap interval (from stage 5), and the 2.5th to 97.5th
percentile range of the group mean across the 1,000 control gene sets. Under the
null, a single spot's norm_score has mean 0 and SD of about 1 against its own
control sets. The null for a group mean is much narrower, because averaging over
many spots shrinks the control-set spread, so a mean well below 1 can still fall
outside it.

The bootstrap resamples spots with the gene set held fixed. The scDRS null holds
the spots fixed and varies the gene set, which is the larger source of
variability: across the 12 groups the null SD of the group mean is 2.0-14.4 times
the bootstrap SE (median 5.6). An interval that excludes zero is therefore not
the criterion for significance. Marker fill follows the scDRS group call, which
tests the 95th percentile of norm_score (assoc_fdr), not the mean, and the two
tests need not agree.

Outputs  Fig-supp-null-vs-bootstrap.{png,pdf},
         results/occipital/mean_vs_q95_test_comparison.csv  (read by stage 7)
"""
import os, json

os.environ["NUMBA_CACHE_DIR"] = os.getcwd() + "/.numba_cache"

import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def apply_figure_style(sizes=(8, 7, 6)):
    """Title / label / tick point sizes, vector text, no autolayout."""
    title_pt, label_pt, tick_pt = sizes
    mpl.rcParams.update({"figure.dpi": 100, "savefig.dpi": 300,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "figure.autolayout": False,
                         "axes.titlesize": title_pt, "axes.labelsize": label_pt,
                         "xtick.labelsize": tick_pt, "ytick.labelsize": tick_pt,
                         "legend.fontsize": tick_pt, "axes.spines.top": False,
                         "axes.spines.right": False})


def panel_letter(a, s, dx=-0.09, dy=1.02):
    a.text(dx, dy, s, transform=a.transAxes, fontsize=9, fontweight="bold",
           va="bottom", ha="left")


apply_figure_style(sizes=(8, 7, 6))

# Full scDRS output: each scored spot's observed norm_score and its norm_score
# under each of the 1,000 control gene sets.
full_score = pd.read_csv(
    "results/occipital/score_file/253_genes/ASC_Pfdr.full_score.gz",
    sep="\t", index_col=0)
ctrl_cols = [col_name for col_name in full_score.columns
             if col_name.startswith("ctrl_norm_score")]
ctrl_norm = full_score[ctrl_cols].to_numpy(np.float64)
obs_score = full_score["norm_score"].to_numpy(np.float64)

# Group of each scored spot (None for spots outside the 12 groups).
grouped_spots = pd.read_csv("results/occipital/spot_layer_assignment.csv")
grouped_spots = grouped_spots[(grouped_spots.panel_group != "excluded")
                              & grouped_spots.scored]
group_of_spot = pd.Series(
    [dict(zip("A1_" + grouped_spots.barcode, grouped_spots.panel_group)).get(k)
     for k in full_score.index],
    index=full_score.index)
group_stats = pd.read_csv("results/occipital/layer_group_stats.csv")

rows = []
for group in group_stats.panel_group:
    in_group = (group_of_spot == group).to_numpy()
    # Null of the group mean: for each control gene set, the mean norm_score over
    # the group's spots.
    null_means = ctrl_norm[in_group].mean(0)
    obs_mean = obs_score[in_group].mean()
    stats = group_stats[group_stats.panel_group == group].iloc[0]
    rows.append({"panel_group": group, "series": stats.series, "layer": stats.layer,
                 "n": int(in_group.sum()), "obs_mean": obs_mean,
                 "null_lo": float(np.quantile(null_means, .025)),
                 "null_hi": float(np.quantile(null_means, .975)),
                 "null_sd": float(null_means.std()),
                 # One-sided Monte-Carlo p for the group mean, with the +1
                 # correction so that it cannot reach zero. This is not the q95
                 # test that decides marker fill.
                 "mean_mcp": ((1 + int((null_means >= obs_mean).sum()))
                              / (1 + len(null_means))),
                 "ci_lo": stats.ci_lo, "ci_hi": stats.ci_hi,
                 # Standard error implied by the width of the 95% interval.
                 "boot_se": (stats.ci_hi - stats.ci_lo) / (2 * 1.96),
                 "assoc_fdr": stats.assoc_fdr})
null_vs_boot = pd.DataFrame(rows)
null_vs_boot["null_sd_over_boot_se"] = null_vs_boot.null_sd / null_vs_boot.boot_se
# Two-sided check against the central 95% of the control-set means; the scDRS
# group test itself is one-sided.
null_vs_boot["outside_null"] = ((null_vs_boot.obs_mean > null_vs_boot.null_hi)
                                | (null_vs_boot.obs_mean < null_vs_boot.null_lo))
null_vs_boot.to_csv("results/occipital/mean_vs_q95_test_comparison.csv", index=False)

LAYER_ORDER = ["VZ", "iSVZ", "oSVZ", "IZ", "SP", "L5/6", "L4", "upper plate"]
null_vs_boot["layer"] = pd.Categorical(null_vs_boot.layer, LAYER_ORDER, ordered=True)
null_vs_boot = null_vs_boot.sort_values(["layer", "series"]).reset_index(drop=True)
SERIES_COLOUR = {"V1": "#3C8CC3", "V2": "#4F9E4F", "occipital (not split)": "#6E7B8B"}

# One row per group. Grey bar: central 95% of the control-set group means. Coloured
# line: 95% block-bootstrap interval of the observed mean. Marker: observed mean,
# filled when the scDRS group call has FDR < 0.05. The ratio of null SD to
# bootstrap SE is saved in the table (null_sd_over_boot_se) rather than plotted.
fig, ax = plt.subplots(figsize=(7.4, 4.6))

y = np.arange(len(null_vs_boot))
for i, row in null_vs_boot.iterrows():
    colour = SERIES_COLOUR[row.series]
    ax.hlines(i, row.null_lo, row.null_hi, color="0.82", lw=7, zorder=1,
              capstyle="butt")
    ax.hlines(i, row.ci_lo, row.ci_hi, color=colour, lw=2.2, zorder=3)
    ax.plot(row.obs_mean, i, "o", ms=5.2,
            mfc=colour if row.assoc_fdr < 0.05 else "white",
            mec=colour, mew=1.1, zorder=4)
ax.axvline(0, color="0.45", lw=0.8, ls="--", zorder=0)
ax.set_yticks(y)
ax.set_yticklabels([f"{row.series.replace(' (not split)','')} {row.layer}"
                    for _, row in null_vs_boot.iterrows()])
ax.set_xlabel("group mean ASD-253 scDRS Z-score (norm_score)")
ax.set_title("A mean below 1 can still clear the gene-set null, and the\n"
             "bootstrap interval is not the yardstick the p-value uses", loc="left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# Extra room on the right for the two annotations.
RIGHT_MARGIN = 2.55
ax.set_xlim(min(null_vs_boot.null_lo.min(), null_vs_boot.obs_mean.min()) - 0.30,
            max(null_vs_boot.null_hi.max(), null_vs_boot.obs_mean.max()) + RIGHT_MARGIN)
ax.set_ylim(-0.8, len(null_vs_boot) - 0.2)

# The two annotated groups. The wording of both annotations (outside or inside
# the null band) is fixed text; only the means and p-values are read from the
# table, so the wording needs re-checking if the inputs change.
iz = int(null_vs_boot.index[null_vs_boot.panel_group == "occipital|IZ"][0])
v2sp = int(null_vs_boot.index[null_vs_boot.panel_group == "V2|SP"][0])
# The annotation text starts to the right of every drawn element. The leader lines
# do cross other rows on their way to their markers; that is intended.
DATA_RIGHT = float(np.max(np.c_[null_vs_boot.null_hi, null_vs_boot.ci_hi,
                                null_vs_boot.obs_mean]))
xtext = DATA_RIGHT + 0.18
ax.annotate(f"mean {null_vs_boot.loc[iz,'obs_mean']:.2f} is OUTSIDE the null\n"
            f"band (p={null_vs_boot.loc[iz,'mean_mcp']:.3f}) though it is < 1",
            xy=(null_vs_boot.loc[iz, "obs_mean"], iz), xytext=(xtext, iz - 0.55),
            fontsize=6, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="0.35",
                            shrinkA=0, shrinkB=3))
ax.annotate(f"mean {null_vs_boot.loc[v2sp,'obs_mean']:.2f} is INSIDE the null\n"
            f"band (p={null_vs_boot.loc[v2sp,'mean_mcp']:.3f}) yet drawn filled",
            xy=(null_vs_boot.loc[v2sp, "obs_mean"], v2sp), xytext=(xtext, v2sp + 0.75),
            fontsize=6, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="0.35",
                            shrinkA=0, shrinkB=3))

# Legend below the axes, clear of the data.
ax.legend(handles=[
    Line2D([0], [0], color="0.82", lw=7, label="control gene-set null, 2.5\u201397.5%"),
    Line2D([0], [0], color="0.45", lw=2.2, label="95% block-bootstrap CI"),
    Line2D([0], [0], marker="o", color="none", mfc="0.45", mec="0.45", ms=5,
           label="group q95 FDR < 0.05"),
    Line2D([0], [0], marker="o", color="none", mfc="white", mec="0.45", ms=5,
           label="group q95 FDR \u2265 0.05")],
    frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
    fontsize=6, handletextpad=0.6, columnspacing=1.6)

panel_letter(ax, "a")

fig.savefig("Fig-supp-null-vs-bootstrap.png", dpi=300, bbox_inches="tight")
fig.savefig("Fig-supp-null-vs-bootstrap.pdf", bbox_inches="tight")

# Layout checks on the saved figure.
r_ = fig.canvas.get_renderer()
leg_texts = {id(t) for t in ax.get_legend().get_texts()}
texts = [(t, t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible() and id(t) not in leg_texts]
ovl = [(a.get_text()[:26], b.get_text()[:26])
       for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
assert xtext > DATA_RIGHT, (
    "annotation text must start to the right of every plotted interval and "
    f"marker: text at {xtext:.3f}, data reach {DATA_RIGHT:.3f}")
assert ax.get_xlim()[1] > xtext + 2.0, "right margin too narrow for the annotation text"
ann_left = [(t.get_text()[:18], round(xtext, 3)) for t, _ in texts
            if t.get_text().startswith("mean ")]
assert len(ann_left) == 2, (
    f"both annotations (IZ and V2|SP) must be present: {ann_left}")
print(json.dumps({
    "n_groups": len(null_vs_boot),
    "median_null_sd_over_boot_se": round(
        float(null_vs_boot.null_sd_over_boot_se.median()), 2),
    "range_ratio": [round(float(null_vs_boot.null_sd_over_boot_se.min()), 2),
                    round(float(null_vs_boot.null_sd_over_boot_se.max()), 2)],
    "groups_mean_inside_null_but_filled": sorted(
        null_vs_boot.loc[(~null_vs_boot.outside_null)
                         & (null_vs_boot.assoc_fdr < 0.05), "panel_group"]),
    "groups_mean_outside_null": sorted(
        null_vs_boot.loc[null_vs_boot.outside_null, "panel_group"]),
    "text_overlaps": ovl,
    "annotation_text_x": round(xtext, 3),
    "rightmost_drawn_data_x": round(DATA_RIGHT, 3),
}, indent=1))
