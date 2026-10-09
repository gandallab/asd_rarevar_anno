#!/usr/bin/env python
"""
smith-asc-figure-2-spatial, stage 7 of the pipeline (src/00-08): documentation.
Computes the pattern summaries from the occipital group table and writes the Methods
Methods document for the figure.

Counts, group statistics, the uniform-weight comparison and the null-versus-bootstrap
ratios in the text are read at run time from the results tables and the scDRS logs.
The statements about the shape of the laminar profile (where it peaks, where it falls,
which groups are not enriched) are checked against the group table before they are
written. Typed into the template and not tested against data: the group definitions and
panel inventory, the remark that the radial series is monotonic in the group mean and in
obs_q95, the laminar profile
run (V2|SP), and the statement that the two spots lost to the scDRS cell filter are the
spots lost to the scDRS cell filter.

This script is not part of the analysis and cannot be run from the public repository.
It reads results/occipital/*, results/panel_i_group_stats.csv, results/03_assignment.json,
the scDRS stdout in logs/ and the gene set in genesets/, and takes versions from
handoff/versions_python.json or, failing that, SOFTWARE_VERSIONS.csv. Of these, only
results/occipital/weighted_vs_uniform.csv and SOFTWARE_VERSIONS.csv are shipped.

Outputs  METHODS_occipital_only.md, results/occipital/07_pattern_checks.json
"""
import json, os, re
import numpy as np, pandas as pd

if os.path.exists("handoff/versions_python.json"):
    V_PY = json.load(open("handoff/versions_python.json"))
else:                                 # no live capture available; read the shipped table
    sv = pd.read_csv("SOFTWARE_VERSIONS.csv")
    V_PY = dict(zip(sv.component, sv.version))

substrate = json.load(open("results/occipital/01_substrate.json"))
assignment = json.load(open("results/occipital/03_assignment.json"))
figure_qc = json.load(open("results/occipital/05_figure_qc.json"))
group_stats = pd.read_csv("results/occipital/layer_group_stats.csv")
PLATE = ["SP", "L5/6", "L4", "upper plate"]  # data keys; displayed as L2/3
RADIAL = ["VZ", "iSVZ", "oSVZ", "IZ"]
DISPLAY_NAME = {"upper plate": "L2/3"}


def as_text(values):
    """Comma-separated values, printed exactly as the numbers are stored."""
    return ", ".join(repr(v) for v in values)


# ---- pattern summaries ------------------------------------------------------
# Computed on assoc_mcz, the statistic the layer panel plots. The group mean of
# norm_score gives a different laminar shape (the radial series rises monotonically in
# the mean but not in assoc_mcz), so the text has to describe the plotted statistic.
STAT = "assoc_mcz"
checks = {"statistic": STAT}
for ser in ["V1", "V2"]:
    rows = group_stats[group_stats.series == ser].set_index("layer").reindex(PLATE)
    mcz = rows[STAT].to_numpy(float)
    checks[f"{ser}_plate_mcz"] = [round(float(x), 3) for x in mcz]
    checks[f"{ser}_monotonic_rising"] = bool(np.all(np.diff(mcz) > 0))
    checks[f"{ser}_argmax_layer"] = PLATE[int(np.argmax(mcz))]
    checks[f"{ser}_falls_at_outermost"] = bool(mcz[-1] < mcz[-2])
radial = (group_stats[group_stats.series == "occipital (not split)"]
          .set_index("layer").reindex(RADIAL))
radial_mcz = radial[STAT].to_numpy(float)
checks["radial_mcz"] = [round(float(x), 3) for x in radial_mcz]
checks["radial_monotonic_rising"] = bool(np.all(np.diff(radial_mcz) > 0))
checks["radial_negative_through_oSVZ"] = bool(np.all(radial_mcz[:3] < 0) and radial_mcz[3] > 0)
v1_mcz = (group_stats[group_stats.series == "V1"].set_index("layer")
          .reindex(PLATE)[STAT].to_numpy(float))
v2_mcz = (group_stats[group_stats.series == "V2"].set_index("layer")
          .reindex(PLATE)[STAT].to_numpy(float))
checks["V1_exceeds_V2_in_all_plate_layers"] = bool(np.all(v1_mcz > v2_mcz))
checks["V1_minus_V2_by_plate_layer"] = dict(
    zip(PLATE, [round(float(x), 3) for x in v1_mcz - v2_mcz]))
checks["V1_minus_V2_smallest_layer"] = PLATE[int(np.argmin(v1_mcz - v2_mcz))]
checks["n_enriched_fdr05"] = int((group_stats.assoc_fdr < 0.05).sum())
checks["groups_not_enriched"] = sorted(
    group_stats.loc[group_stats.assoc_fdr >= 0.05, "panel_group"])
checks["n_cell_range"] = [int(group_stats.n_cell.min()), int(group_stats.n_cell.max())]
# A group can have a positive q95 and still not be enriched: the test is against matched
# control gene sets scored in the same spots, not against zero.
checks["groups_positive_q95_but_not_enriched"] = sorted(
    group_stats.loc[(group_stats.obs_q95 > 0) & (group_stats.assoc_fdr >= 0.05), "panel_group"])
# Monte-Carlo resolution: with 1,000 control sets the smallest attainable p is 1/1001,
# so groups sitting on that floor must be reported as "<", not "=".
checks["n_groups_at_mc_floor"] = int((group_stats.assoc_mcp <= 1 / 1001 + 1e-12).sum())
fdr_at_mc_floor = float(
    group_stats.loc[group_stats.assoc_mcp <= 1 / 1001 + 1e-12, "assoc_fdr"].max())

# Ratio of the control gene-set null sd of the group mean to the block-bootstrap se, read
# from the saved comparison table (stage 8). Section 4 quotes its range and median.
null_vs_bootstrap = pd.read_csv("results/occipital/mean_vs_q95_test_comparison.csv")
null_sd_over_boot_se = null_vs_bootstrap.null_sd_over_boot_se
checks["n_groups"] = len(null_vs_bootstrap)
checks["n_groups_null_sd_gt_boot_se"] = int((null_sd_over_boot_se > 1).sum())
checks["null_sd_over_boot_se_min"] = float(null_sd_over_boot_se.min())
checks["null_sd_over_boot_se_median"] = float(null_sd_over_boot_se.median())
checks["null_sd_over_boot_se_max"] = float(null_sd_over_boot_se.max())

# Sections 4 and 5 state these patterns outright; check them against the table first.
assert not (checks["V1_monotonic_rising"] or checks["V2_monotonic_rising"]), (
    "section 5 says the cortical-plate profile is not monotonic in either area")
assert checks["V1_falls_at_outermost"] and checks["V2_falls_at_outermost"], (
    "section 5 says the profile falls at the outermost band in both areas")
assert checks["V1_argmax_layer"] == checks["V2_argmax_layer"] == "L4", (
    f"section 5 says enrichment peaks at L4 in both areas; peaks are "
    f"{checks['V1_argmax_layer']} (V1) and {checks['V2_argmax_layer']} (V2)")
assert checks["radial_negative_through_oSVZ"] and not checks["radial_monotonic_rising"], (
    "section 5 says VZ, iSVZ and oSVZ are negative, IZ is positive, and the radial "
    "series is not monotonic in assoc_mcz")
assert checks["groups_not_enriched"] == ["occipital|VZ", "occipital|iSVZ", "occipital|oSVZ"], (
    f"section 5 says the three groups that are not enriched are the proliferative zones; "
    f"found {checks['groups_not_enriched']}")
assert (group_stats.loc[group_stats.assoc_fdr >= 0.05, STAT] < 0).all(), (
    "section 5 says every group that is not enriched has a negative assoc_mcz")
assert checks["n_groups_null_sd_gt_boot_se"] == checks["n_groups"], (
    "section 4 says the control gene-set null is wider than the bootstrap in every group")

# ---- the substrate the scDRS command line actually scored ---------------------
# Read back from its own stdout rather than assumed: n_gene is the gene universe
# surviving filter_genes, and the Trait line's n_gene is how many of the 253
# gene-set genes were left among them. Both are smaller than the deposited
# feature space because a gene must be detected in 50 spots of this section.
scoring_log = open("logs/02_compute_score.clean.log", errors="replace").read()
N_GENES = int(re.search(r"n_cell=\d+, n_gene=(\d+)", scoring_log).group(1))
N_GS_SCORED = int(re.search(r"Trait=\S+, n_gene=(\d+)", scoring_log).group(1))
checks.update({"n_genes_scored": N_GENES, "n_geneset_genes_scored": N_GS_SCORED})

# ---- gene-set weights, read from the .gs file the CLI was given -----------------
# Earlier drafts of the Methods called the set unweighted; the file shows otherwise.
GROUP_TABLE = "results/occipital/layer_group_stats.csv"
gs_field = open("genesets/ASC_TADA_Pfdr.253.gs").read().strip().split("\n")[1].split("\t")[1]
gs_weights = [float(x.split(":")[1]) for x in gs_field.split(",")]
N_GS_GENES = len(gs_weights)
N_DISTINCT_WEIGHTS, WEIGHT_MIN, WEIGHT_MAX = len(set(gs_weights)), min(gs_weights), max(gs_weights)

# Uniform-weight sensitivity: same flags and spots, with the .gs file replaced by
# GENE:1. The rescoring was run separately and its table is shipped; section 2 quotes
# three numbers from it.
uniform_comparison = pd.read_csv("results/occipital/weighted_vs_uniform.csv")
from scipy.stats import spearmanr
UNIFORM_RHO = float(spearmanr(uniform_comparison.assoc_mcz_weighted,
                              uniform_comparison.assoc_mcz_uniform)[0])
UNIFORM_MEAN_ABS_DELTA = float(uniform_comparison.delta_assoc_mcz.abs().mean())
flipped_call = uniform_comparison[uniform_comparison["enriched_weighted_fdr_lt_0.05"]
                                  != uniform_comparison["enriched_uniform_fdr_lt_0.05"]]
assert len(flipped_call) == 1, (
    f"section 2 says the uniform-weight rescoring changes exactly one enrichment call; "
    f"changed: {flipped_call.group.tolist()}")
UNIFORM_FLIP_GROUP = str(flipped_call.group.iloc[0])
UNIFORM_FLIP_FROM = float(flipped_call.assoc_fdr_weighted.iloc[0])
UNIFORM_FLIP_TO = float(flipped_call.assoc_fdr_uniform.iloc[0])
checks.update({"n_distinct_gene_weights": N_DISTINCT_WEIGHTS,
               "gene_weight_range": [WEIGHT_MIN, WEIGHT_MAX],
               "uniform_vs_weighted_spearman": round(UNIFORM_RHO, 4),
               "uniform_vs_weighted_mean_abs_delta": round(UNIFORM_MEAN_ABS_DELTA, 4),
               "uniform_vs_weighted_call_flips": [UNIFORM_FLIP_GROUP]})
json.dump(checks, open("results/occipital/07_pattern_checks.json", "w"), indent=1)

# Packages that set the occipital-only UMAP. The text names those whose versions were not
# recorded; the recorded ones are in SOFTWARE_VERSIONS.csv.
UMAP_PACKAGES = ["umap-learn", "pynndescent", "llvmlite"]
umap_unrecorded = [f"`{p}`" for p in UMAP_PACKAGES if not V_PY.get(p)]
umap_note = ""
if umap_unrecorded:
    names = (", ".join(umap_unrecorded[:-1]) + " and " + umap_unrecorded[-1]
             if len(umap_unrecorded) > 1 else umap_unrecorded[0])
    umap_note = (f"\nThe embedding in panels c–f depends on {names}, for which no version "
                 f"was recorded or pinned in the environment files.")

n_enriched = checks["n_enriched_fdr05"]
n_positive_q95_not_enriched = len(checks["groups_positive_q95_but_not_enriched"])
single = n_positive_q95_not_enriched == 1
bootstrap = figure_qc["bootstrap"]
v1_minus_v2_text = ", ".join(f"{DISPLAY_NAME.get(k, k)} {v:+.3f}"
                             for k, v in checks["V1_minus_V2_by_plate_layer"].items())
smallest_gap_layer = checks["V1_minus_V2_smallest_layer"]
smallest_gap_layer = DISPLAY_NAME.get(smallest_gap_layer, smallest_gap_layer)
top_group = group_stats.loc[group_stats[STAT].idxmax(), "panel_group"]

methods_text = f"""# Methods — ASD-253 laminar enrichment, GW20 occipital cortex

Everything the figure rests on is described here. A prefrontal section from the
same deposit was scored alongside this one in an earlier version of the analysis
and has been dropped: the two sections come from different donors, so no
prefrontal-versus-occipital contrast is separable from a donor difference. No
panel, table or statistic below involves it.

## 1. Substrate and covariates

Section {substrate['section']}, from a GW20 donor, was scored with `scdrs
compute-score`. `norm_score` is standardised against control gene sets within the
scored object, so it is a statement about this section and carries no meaning
outside it: scores, and the colour scale built on them, should not be compared
with any other scoring run.

The covariates are `const` and `n_genes`. With a single section a section
indicator would be constant and therefore degenerate, so none is included; donor
is likewise constant and cannot be adjusted for.

Layer labels are the authors' own and are used directly, with no transfer step
and so no label-transfer uncertainty. There is no integration step either: with
one section there is no batch to correct. The embedding drawn in panels c-f
(recipe: {assignment['embedding']}) was computed on the same filtered matrix that
scDRS scored. It is used for display only, and no statistic in this figure
depends on it.

## 2. Scoring

The scDRS flags are identical to those of the reference driver
(`src/reference/run_scdrs_rvas.sh.reference`), including
`--flag-return-ctrl-raw-score False`: the per-spot control scores are needed only
for the diagnostic of stage 7, which reads them from the full score file.

The gene set is weighted. `munge-gs --weight zscore` converts each gene's TADA
statistic into a weight. The supplied `.gs` file carries {N_DISTINCT_WEIGHTS} distinct values
from {WEIGHT_MIN:.4f} to {WEIGHT_MAX:.4f}, and the CLI passes them to `score_cell` as
`gene_weight`, so they enter the score. Earlier drafts of the manuscript Methods
described the set as unweighted (binary); the Methods are being amended to match what
was run. Rescoring this section with uniform weights (identical flags, the `.gs` file
replaced by `GENE:1`) leaves the laminar pattern intact (Spearman correlation of
`assoc_mcz` {UNIFORM_RHO:.3f}, mean |Δ| {UNIFORM_MEAN_ABS_DELTA:.2f}) and changes exactly one
call: {UNIFORM_FLIP_GROUP}, `assoc_fdr` {UNIFORM_FLIP_FROM:.3f} → {UNIFORM_FLIP_TO:.3f}. The
comparison over {len(uniform_comparison)} groups is in
`results/occipital/weighted_vs_uniform.csv`. The call for {UNIFORM_FLIP_GROUP} depends on the
weighting and on the choice of statistic, so it should be read as marginal rather
than established.

Of {substrate['n_spots']:,} spots, {figure_qc['n_spots_scored']:,} survive
`--flag-filter-data True`; the two removed fall below `filter_cells(min_genes=250)`.
The gene filter (`filter_genes(min_cells=50)`) reduces {substrate['n_genes_input']:,}
genes to {N_GENES:,}, and leaves {N_GS_SCORED} of the {N_GS_GENES} gene-set genes
scored. Both thresholds are applied within this section, so a gene must be detected in
50 of its spots to be counted at all.

## 3. Groups and multiple testing

The {len(group_stats)} groups are the four radial occipital compartments, which are not split by
area (VZ, iSVZ, oSVZ, IZ), and the four cortical-plate bands, which are split by area
(V1 and V2 × SP, L5/6, L4, upper plate, displayed as L2/3). The opposing cortical wall,
the authors' dropped cluster 14 and the two spots removed by the filter are excluded
from the statistics. Each of the six spot panels (a–f) draws
{figure_qc['n_spots_displayed_spatial']:,} spots, and
{assignment['n_spots_excluded_from_panel']:,} spots are excluded from the groups.

Benjamini–Hochberg correction is applied over these {len(group_stats)} groups.

`assoc_mcz` is a within-object statistic and does not transfer between scoring runs.
scDRS standardises against control gene sets drawn from the scored object, so the gene
universe, the mean- and variance-matched control bins built over it, the covariate
model and the set of gene-set genes surviving the gene filter all belong to this run.
Rescoring the same spots in a different object — pooled with another section, or with a
different gene filter — moves the z and can move an enrichment call across the 0.05
line, as the uniform-weight comparison in section 2 shows for {UNIFORM_FLIP_GROUP}.
Compare enrichment calls with other analyses, not magnitudes.

## 4. The layer panel plots the manuscript's Enrichment Z-score

Panel g plots `assoc_mcz` for each group, labelled "Enrichment Z-score" as in the rest
of the manuscript. Markers are filled where `assoc_fdr` < 0.05 and hollow otherwise.
All markers are one size, and spot counts are in the `n_cell` column of
`layer_group_stats.csv`. Zero is the matched-control expectation. Colours, display
names and marker conventions follow the collaborator's panel script
(`src/reference/yundan_panel_i_assoc_mcz.py.reference`), so this panel is drawn on the
same convention as the other scDRS panels in the manuscript.

`assoc_mcz` is plotted because it is the standardised form of the group q95, the
statistic that `assoc_mcp` tests and on which the filled or hollow call rests, so a
group's position and its significance share one statistic. The group mean of
`norm_score` is not the tested statistic, and the two can disagree:
{n_positive_q95_not_enriched} {'group' if single else 'groups'}
({', '.join(checks['groups_positive_q95_but_not_enriched'])})
{'has' if single else 'have'} a positive `obs_q95` yet {'is' if single else 'are'} not
enriched, because the matched-control sets reach higher still in those spots. A
spot-block bootstrap interval for the group mean, with the gene set held fixed, is not
a substitute for the test either. The control gene-set null of the group mean is wider
than the bootstrap standard error in {checks['n_groups_null_sd_gt_boot_se']} of
{checks['n_groups']} groups, by a factor of {checks['null_sd_over_boot_se_min']:.1f} to
{checks['null_sd_over_boot_se_max']:.1f} (median {checks['null_sd_over_boot_se_median']:.1f}), so
an interval that excludes zero is not equivalent to significance. The group mean, its
95% bootstrap interval (square tiles of about {bootstrap['tile_k_hex_spacings']}
nearest-neighbour hex spacings, {bootstrap['n_blocks_range'][0]}–{bootstrap['n_blocks_range'][1]}
blocks per group, {bootstrap['B']:,} resamples, seed {bootstrap['seed']}) and `obs_q95`
remain in the group table and are not plotted. The diagnostic comparison is
`Fig-supp-null-vs-bootstrap.{{png,pdf}}` with
`results/occipital/mean_vs_q95_test_comparison.csv`.

Draft legend for panel g. Enrichment of the ASD-253 gene set by layer in one
second-trimester (GW20) occipital section from one donor. Each marker is one group of
spots ({len(group_stats)} groups). Position is `assoc_mcz`, the scDRS Monte-Carlo z of the
group's 95th-percentile `norm_score` against 1,000 matched control gene sets, and the
dashed line at zero is the control expectation. Filled markers have `assoc_fdr` < 0.05
(one-sided upper test, Benjamini–Hochberg over the {len(group_stats)} groups); hollow markers
do not. Negative values mean that a group is not enriched; no test for depletion was
run.

Panels a–f map the spots of section A1: a shows the authors' laminar annotation and b
the per-spot `norm_score`; c shows the occipital-only UMAP coloured by layer and d the
same UMAP by `norm_score`; e and f show the same UMAP with V1 and V2 spots highlighted
({int(figure_qc['n_v1_highlighted']):,} and {int(figure_qc['n_v2_highlighted']):,} spots).
Panel g is the layer panel. Panels b and d share one diverging colour scale, centred on
zero (the matched-control expectation) and spanning the 1st to 99th percentile of
`norm_score` over the spots in the embedding. The grey backdrop in e and f is the
complement of the highlighted set. It includes the germinal spots, which carry no areal
label because the authors' annotation does not split them by area, so they are grey in
both panels. The six spot panels draw the same
{int(figure_qc['n_spots_displayed_spatial']):,} layer-annotated scored spots.

Only the four cortical-plate bands are split by area in the authors' annotation (V1
alone separates L2 from L3, which this vocabulary collapses), whereas the germinal
bands are shared between the two areas and carry no areal label. The collapsed band is
labelled L2/3 in the figure, after the source annotation's own name for it in V2, but
the data keep the earlier name: the group key in `{GROUP_TABLE}`, in the scDRS group
output and in the code is `upper plate`, and only the tick and legend labels are
remapped (`LABEL_MAP` in the figure script). The collaborator's panel renderer
hard-codes that remapping, so the key was left unchanged. The two names denote one set
of spots.

This figure contains no prefrontal data: every panel is section A1, and the areal
series are V1, V2 and their union in the shared germinal zones. The prefrontal
comparison, and the donor confound it would carry, is out of scope.

`assoc_mcp` is one-sided (upper tail), as in the main analysis. A proliferative zone
with a large negative `assoc_mcz` is not enriched. The same 1,000 control sets would
support a lower-tail test, but none was run, so this is not evidence of depletion.

## 5. Results

Of the {len(group_stats)} groups, {n_enriched} reach `assoc_fdr` < 0.05. The three that do not are
the proliferative germinal zones ({', '.join(checks['groups_not_enriched'])}), all with
negative `assoc_mcz`.

The values below are `assoc_mcz`, the statistic the panel plots, taken from the full
group table. The laminar profile is not monotonic and should not be described as rising
continuously outward.

- V1 cortical plate, SP → L2/3: {as_text(checks['V1_plate_mcz'])}. The profile peaks at
  {checks['V1_argmax_layer']} and falls at the outermost band.
- V2 cortical plate, SP → L2/3: {as_text(checks['V2_plate_mcz'])}. The profile peaks at
  {checks['V2_argmax_layer']} and falls at the outermost band.
- Radial series, VZ → IZ: {as_text(checks['radial_mcz'])}. VZ, iSVZ and oSVZ are negative
  and only IZ is positive. The series is not monotonic in `assoc_mcz` (it falls from VZ
  to iSVZ), although it is monotonic in the group mean and in `obs_q95`. The three
  statistics describe different things and should not be substituted for one another.

Enrichment is therefore confined to the cortical plate and the intermediate zone,
peaks at L4 in both areas, and falls toward the pia. The largest value is
{group_stats[STAT].max():.3f} in {top_group}.

V1 exceeds V2 in {'all four' if checks['V1_exceeds_V2_in_all_plate_layers'] else 'some'}
cortical-plate bands (V1 − V2: {v1_minus_v2_text}).
The smallest gap is in {smallest_gap_layer}, where the two areas are effectively
indistinguishable, so the consistent direction
across bands is the observation, not any per-band magnitude. In addition, `assoc_mcz`
is a z against each group's own control-set null, whose width depends on the number of
spots in the group ({checks['n_cell_range'][0]}–{checks['n_cell_range'][1]} here), so it
cannot be compared across groups of unequal size in the way an effect size could. No
formal V1-versus-V2 test is reported, and none of these differences has been tested.

The V1-V2 contrast that remains lies within one section and one donor, so it is not
confounded by donor. It is not free of
position effects: V1 and V2 occupy different parts of the section (V1 on the right and
V2 on the left in the displayed orientation), so the contrast also carries any
difference in section quality or cutting plane between those positions. It is within
donor but not within structure.

Of the {len(group_stats)} groups, {checks['n_groups_at_mc_floor']} sit at the Monte-Carlo floor
`assoc_mcp` = 1/1001 with `--n-ctrl 1000`. They should be reported as *p* < 0.001 and
not as *p* = 0.001, and their FDR values ({fdr_at_mc_floor:.3f}) are likewise limited by the
resolution of the Monte-Carlo test and are not exact.

## 6. Limitations specific to this variant

The data are one donor, one section and one timepoint (GW20), so nothing here is
replicated, and the bootstrap does not supply replication. The V1/V2 boundary is the
authors' annotation of that single section. Cluster 14, which the authors dropped, is
not drawn and is excluded from every group statistic. No sensitivity analysis that
includes those spots was run, so the effect of excluding them on the V1 and V2 group
estimates is untested. The opposing cortical wall is identified from the suffix of the
authors' labels (-L, -L1, -L2) alone.

Spots are not independent observations. Neighbouring spots are spatially correlated, and
the scDRS null resamples genes and not spots, so the Monte-Carlo p-values do not account
for that correlation. The group table also reports `hetero_mcp`, a test for
heterogeneity of the score within each group; where it is small, the group statistic
summarises a score that varies within the group.

The laminar gradient may largely reflect the cell-type composition of the layers
(neurons in the cortical plate and intermediate zone, progenitors in the germinal
zones). It shows where the aggregate signal is concentrated, not which layers are
vulnerable, and no cell-type deconvolution was done. The embedding comes from a
different pipeline from any other embedding in the manuscript, so the UMAP panels should not
be compared across figures, and the scores are not comparable either (section 1).

Software: scDRS {V_PY['scdrs']}, scanpy {V_PY['scanpy']}, numpy {V_PY['numpy']},
Python {V_PY['python']}. No R stage runs in this variant; the archived export matrices
it reads were produced by the R stage documented in the main methods.{umap_note}
"""

open("METHODS_occipital_only.md", "w").write(methods_text)
print(f"addendum: {len(methods_text):,} chars")
print(json.dumps(checks, indent=1))
