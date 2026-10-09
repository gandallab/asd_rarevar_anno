# Methods — ASD-253 laminar enrichment, GW20 occipital cortex

Everything the figure rests on is described here. A prefrontal section from the
same deposit was scored alongside this one in an earlier version of the analysis
and has been dropped: the two sections come from different donors, so no
prefrontal-versus-occipital contrast is separable from a donor difference. No
panel, table or statistic below involves it.

## 1. Substrate and covariates

Section A1 / FB080-O1 (occipital), from a GW20 donor, was scored with `scdrs
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
(recipe: scanpy normalise/HVG(2000)/scale/PCA(30)/neighbors(15)/UMAP, no integration) was computed on the same filtered matrix that
scDRS scored. It is used for display only, and no statistic in this figure
depends on it.

## 2. Scoring

The scDRS flags are identical to those of the reference driver
(`src/reference/run_scdrs_rvas.sh.reference`), including
`--flag-return-ctrl-raw-score False`: the per-spot control scores are needed only
for the diagnostic of stage 7, which reads them from the full score file.

The gene set is weighted. `munge-gs --weight zscore` converts each gene's TADA
statistic into a weight. The supplied `.gs` file carries 253 distinct values
from 6.2039 to 10.0000, and the CLI passes them to `score_cell` as
`gene_weight`, so they enter the score. Earlier drafts of the manuscript Methods
described the set as unweighted (binary); the Methods are being amended to match what
was run. Rescoring this section with uniform weights (identical flags, the `.gs` file
replaced by `GENE:1`) leaves the laminar pattern intact (Spearman correlation of
`assoc_mcz` 0.993, mean |Δ| 0.09) and changes exactly one
call: V2|SP, `assoc_fdr` 0.035 → 0.059. The
comparison over 12 groups is in
`results/occipital/weighted_vs_uniform.csv`. The call for V2|SP depends on the
weighting and on the choice of statistic, so it should be read as marginal rather
than established.

Of 3,591 spots, 3,589 survive
`--flag-filter-data True`; the two removed fall below `filter_cells(min_genes=250)`.
The gene filter (`filter_genes(min_cells=50)`) reduces 18,085
genes to 13,908, and leaves 242 of the 253 gene-set genes
scored. Both thresholds are applied within this section, so a gene must be detected in
50 of its spots to be counted at all.

## 3. Groups and multiple testing

The 12 groups are the four radial occipital compartments, which are not split by
area (VZ, iSVZ, oSVZ, IZ), and the four cortical-plate bands, which are split by area
(V1 and V2 × SP, L5/6, L4, upper plate, displayed as L2/3). The opposing cortical wall,
the authors' dropped cluster 14 and the two spots removed by the filter are excluded
from the statistics. Each of the six spot panels (a–f) draws
3,461 spots, and
645 spots are excluded from the groups.

Benjamini–Hochberg correction is applied over these 12 groups.

`assoc_mcz` is a within-object statistic and does not transfer between scoring runs.
scDRS standardises against control gene sets drawn from the scored object, so the gene
universe, the mean- and variance-matched control bins built over it, the covariate
model and the set of gene-set genes surviving the gene filter all belong to this run.
Rescoring the same spots in a different object — pooled with another section, or with a
different gene filter — moves the z and can move an enrichment call across the 0.05
line, as the uniform-weight comparison in section 2 shows for V2|SP.
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
1 group
(occipital|oSVZ)
has a positive `obs_q95` yet is not
enriched, because the matched-control sets reach higher still in those spots. A
spot-block bootstrap interval for the group mean, with the gene set held fixed, is not
a substitute for the test either. The control gene-set null of the group mean is wider
than the bootstrap standard error in 12 of
12 groups, by a factor of 2.0 to
14.4 (median 5.6), so
an interval that excludes zero is not equivalent to significance. The group mean, its
95% bootstrap interval (square tiles of about 5
nearest-neighbour hex spacings, 14–47
blocks per group, 2,000 resamples, seed 77) and `obs_q95`
remain in the group table and are not plotted. The diagnostic comparison is
`Fig-supp-null-vs-bootstrap.{png,pdf}` with
`results/occipital/mean_vs_q95_test_comparison.csv`.

Draft legend for panel g. Enrichment of the ASD-253 gene set by layer in one
second-trimester (GW20) occipital section from one donor. Each marker is one group of
spots (12 groups). Position is `assoc_mcz`, the scDRS Monte-Carlo z of the
group's 95th-percentile `norm_score` against 1,000 matched control gene sets, and the
dashed line at zero is the control expectation. Filled markers have `assoc_fdr` < 0.05
(one-sided upper test, Benjamini–Hochberg over the 12 groups); hollow markers
do not. Negative values mean that a group is not enriched; no test for depletion was
run.

Panels a–f map the spots of section A1: a shows the authors' laminar annotation and b
the per-spot `norm_score`; c shows the occipital-only UMAP coloured by layer and d the
same UMAP by `norm_score`; e and f show the same UMAP with V1 and V2 spots highlighted
(780 and 606 spots).
Panel g is the layer panel. Panels b and d share one diverging colour scale, centred on
zero (the matched-control expectation) and spanning the 1st to 99th percentile of
`norm_score` over the spots in the embedding. The grey backdrop in e and f is the
complement of the highlighted set. It includes the germinal spots, which carry no areal
label because the authors' annotation does not split them by area, so they are grey in
both panels. The six spot panels draw the same
3,461 layer-annotated scored spots.

Only the four cortical-plate bands are split by area in the authors' annotation (V1
alone separates L2 from L3, which this vocabulary collapses), whereas the germinal
bands are shared between the two areas and carry no areal label. The collapsed band is
labelled L2/3 in the figure, after the source annotation's own name for it in V2, but
the data keep the earlier name: the group key in `results/occipital/layer_group_stats.csv`, in the scDRS group
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

Of the 12 groups, 9 reach `assoc_fdr` < 0.05. The three that do not are
the proliferative germinal zones (occipital|VZ, occipital|iSVZ, occipital|oSVZ), all with
negative `assoc_mcz`.

The values below are `assoc_mcz`, the statistic the panel plots, taken from the full
group table. The laminar profile is not monotonic and should not be described as rising
continuously outward.

- V1 cortical plate, SP → L2/3: 2.057, 3.709, 4.811, 4.065. The profile peaks at
  L4 and falls at the outermost band.
- V2 cortical plate, SP → L2/3: 2.046, 3.349, 3.908, 3.1. The profile peaks at
  L4 and falls at the outermost band.
- Radial series, VZ → IZ: -3.174, -3.364, -2.719, 2.27. VZ, iSVZ and oSVZ are negative
  and only IZ is positive. The series is not monotonic in `assoc_mcz` (it falls from VZ
  to iSVZ), although it is monotonic in the group mean and in `obs_q95`. The three
  statistics describe different things and should not be substituted for one another.

Enrichment is therefore confined to the cortical plate and the intermediate zone,
peaks at L4 in both areas, and falls toward the pia. The largest value is
4.811 in V1|L4.

V1 exceeds V2 in all four
cortical-plate bands (V1 − V2: SP +0.011, L5/6 +0.360, L4 +0.903, L2/3 +0.966).
The smallest gap is in SP, where the two areas are effectively
indistinguishable, so the consistent direction
across bands is the observation, not any per-band magnitude. In addition, `assoc_mcz`
is a z against each group's own control-set null, whose width depends on the number of
spots in the group (119–698 here), so it
cannot be compared across groups of unequal size in the way an effect size could. No
formal V1-versus-V2 test is reported, and none of these differences has been tested.

The V1-V2 contrast that remains lies within one section and one donor, so it is not
confounded by donor. It is not free of
position effects: V1 and V2 occupy different parts of the section (V1 on the right and
V2 on the left in the displayed orientation), so the contrast also carries any
difference in section quality or cutting plane between those positions. It is within
donor but not within structure.

Of the 12 groups, 6 sit at the Monte-Carlo floor
`assoc_mcp` = 1/1001 with `--n-ctrl 1000`. They should be reported as *p* < 0.001 and
not as *p* = 0.001, and their FDR values (0.002) are likewise limited by the
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

Software: scDRS 1.0.2, scanpy 1.11.5, numpy 1.26.4,
Python 3.11.16. No R stage runs in this variant; the archived export matrices
it reads were produced by the R stage documented in the main methods.
