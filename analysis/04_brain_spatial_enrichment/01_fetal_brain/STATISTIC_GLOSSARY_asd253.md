# Which z is which: ASD-253 manuscript

Three z-valued statistics appear in the manuscript under similar names, and they
are not interchangeable. The manuscript's own figures already separate two of
them within a single figure. The definitions below are taken from the
manuscript Methods and from the column names of the deposited tables.

## 1. Mean scDRS z: a descriptive group summary

This is the arithmetic mean of scDRS `norm_score` over the cells or spots in a
group. `norm_score` is itself a per-cell z of the focal gene set against that
cell's matched control sets. Under the null, each cell's score has mean 0 and SD
about 1 relative to its own control sets. That is a property of the null and not
a constraint on the observed scores: across the spots of this section the
observed SD is 2.05.

The manuscript uses the group mean descriptively in panel b of the subplate
figure, where the axis is labelled verbatim `mean scDRS z`. The Visium layer
panel no longer plots it: `Fig2-occipital-only` plots `assoc_mcz` (§2), whereas
earlier versions of that panel plotted the mean.

The mean is also a test statistic, but only for differences between donor
groups and never for the association of the gene set itself. Each such test
takes the donor as the unit of analysis and uses a permutation, rank or
mixed-model null instead of the control-gene-set null:

- a within-donor anterior-posterior contrast (donor-level deltas of
  mean `norm_score`, paired two-sided Wilcoxon, with a control-set version);
- a between-sex difference of those deltas (exact permutation over all sex-label
  assignments, corroborated by a donor-level random-slope mixed model,
  `norm_score ~ V1 x sex + (1+V1|donor)`);
- a developmental-window contrast of donor means, fetal minus postnatal
  (Mann-Whitney and Kruskal-Wallis, FDR-corrected).

The analyses behind those three contrasts are not part of this deposit and
cannot be reproduced from it, so their results are left out here; see the
manuscript. Effect sizes for them are quoted in "scDRS units", that is, in units
of mean `norm_score`.

The gene-set association itself is always the q95 statistic of §2. A figure that
plots the mean while taking its significance call from `assoc_fdr` mixes the two
statistics. The occipital layer panel did this at one stage and was changed to
plot `assoc_mcz`, so that the plotted point and the call come from the same
statistic.

The null for a group mean is much tighter than the per-cell null SD of about 1
(0.175–0.791 across the 12 occipital groups, and tighter the more spots a group
has). A group mean therefore cannot be read against a ±1.96 rule of thumb.

## 2. `assoc_mcz`: the manuscript's scDRS group association statistic

The Methods define the group-level association as scDRS's native statistic,
"the 95th percentile (q95) of within-group per-cell normalized scores versus
the q95 of each control set", which yields a Monte-Carlo z (`assoc_mcz`) and P
value, BH-corrected across groups within each atlas.

It is reported in every `scdrs3_*` table, whose columns are `group`, `n_cell`,
`obs_q95_score`, `assoc_mcz`, `assoc_mcp`, `n_fdr01`, `prop_fdr01` and
`assoc_fdr`, and in every z quoted in the scDRS results text (the
Cajal-Retzius values, for example). This is the scDRS z of the paper. It is
computed on the q95 and not on the mean. `assoc_mcp` is one-sided upper, so a
group whose q95 lies below its control sets receives a negative `assoc_mcz` and
no enrichment call; the P value has no lower tail and cannot show depletion. The
Methods also state that the statistic is not comparable across atlases, because
it "scales with group size".

## 3. "Enrichment Z": a size-matched resampling statistic

This is not an scDRS statistic and is not computed from single-cell scores. In
the regional atlases (Methods §2.4), each gene's expression was z-scored across
regions, the ASD-gene mean per region formed the observed profile, and
significance was assessed against "a null of 10,000 size-matched gene sets
resampled from brain-expressed genes only" (above the 10th percentile of mean
expression), with BH-FDR correction.

It is used in the 34-parcel regional analysis ("per-region enrichment Z") and in
the module-level analysis, whose table `asd_fig_module_z.csv` carries `z_full`
and `fdr_full` (all-gene background) next to `z_neuronal` and `fdr_neuronal`
(brain-expressed background). Both backgrounds are reported because testing
against all measured genes "inflates significance for anything neuronal". The
subplate figure labels these axes `ASD-253 enrichment z` and `enrichment z`.

The term is overloaded. The collaborator's panel script for the Visium laminar
panel (`src/reference/yundan_panel_i_assoc_mcz.py.reference`) labels its x-axis
"Enrichment Z-score" but plots `assoc_mcz`, the statistic of §2. "Enrichment z"
therefore means the size-matched resampling z in the regional and module
analyses and `assoc_mcz` in the Visium laminar panels. Both figures produced
here follow the panel script and state `assoc_mcz` on the axis.

## Summary

| term in the paper | statistic | computed on | tests what |
|---|---|---|---|
| mean scDRS z | mean of `norm_score` | per-cell scDRS scores | differences between donor groups (donor-level permutation / rank / mixed model); never the gene-set association |
| z / `assoc_mcz` | Monte-Carlo z of q95 | per-cell scDRS scores | group vs matched control gene sets |
| enrichment z (regional/module) | resampling z | regional or module expression | gene set vs size-matched brain-expressed sets |
| "Enrichment Z-score" (Visium laminar panel) | the same `assoc_mcz` as the `z / assoc_mcz` row, not the resampling z of the row directly above | per-spot scDRS scores | group vs matched control gene sets |

## Gene-set weighting: the Visium figures are weighted

Earlier drafts of the Methods described the gene set as run "unweighted
(binary)". The reference driver for the scDRS runs calls
`scdrs munge-gs --weight zscore`, and the supplied `ASC_TADA_Pfdr.253.gs`
carries 253 distinct per-gene weights spanning 6.2039-10.0000 (TADA z-scores
clipped at 10). The scDRS CLI passes those to `score_cell` as `gene_weight`, so a
run using that file is weighted, and both Visium figures here used it.

The weights are kept, and the Methods are being amended to match what was run.
A uniform-weight rescoring of the occipital section changes the laminar pattern
very little (Spearman of `assoc_mcz` 0.993, mean |Δ| 0.09) and flips one
marginal call, V2|SP; the 12-group comparison is shipped as
`results/occipital/weighted_vs_uniform.csv`. Whether the atlas panels were
scored with the same weighted file still needs confirming before Visium and
atlas numbers are pooled.

`weight_opt='vs'`, which the Methods also mention, is something else:
variance standardisation of gene expression inside `score_cell`, not a per-gene
prior weight. Both figures use `weight_opt='vs'` (the scDRS default), as in the
manuscript.
