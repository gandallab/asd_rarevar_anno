"""Shared machinery for the AHBA spatial-enrichment notebooks (01, 02, 05).

Holds the three things those notebooks would otherwise each redefine: the regional
gene-set enrichment test, the Desikan-Killiany surface renderer, and a handful of
empirical-null helpers. Keeping one copy means the maps in Fig 1, Fig 2 and the
module supplement are guaranteed to be drawn on the same surface, at the same
density, with the same medial-wall treatment.

Import from a notebook in this directory:
    import spatial_helpers as sh

Paths are resolved from this file's location, so notebooks work from any working directory.
"""
import os

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import zscore, spearmanr
from statsmodels.stats.multitest import multipletests
from nilearn import datasets, plotting
from nilearn.surface import load_surf_mesh
from scipy.spatial import cKDTree

# --- 0) Constants --- #

# Surface density for every brain render. "fsaverage" is the native 164k mesh: small
# regions (bankssts is ~136 vertices at fsaverage5 vs ~2,200 at 164k) and the
# significance outlines come out publication-crisp, at roughly 7 s per panel.
# SET THIS TO "fsaverage5" (10k verts/hemi) FOR FAST DRAFTS -- the map is constant
# within each DK region, so the coarse mesh is visually near-identical and ~16x faster.
DENSITY = "fsaverage"

MEDIAL_WALL = "#cfc4b0"                     # DK vertices with no AHBA data
ASD_COLOR = "#4d4d4d"
TOPIC_COLORS = {"SYN": "#4a2377", "GR": "#ea801c", "MORPH": "#0d7d87"}
GENE_SET_COLORS = {"ASD": ASD_COLOR, **TOPIC_COLORS}
OCCIPITAL = ["pericalcarine", "lateraloccipital", "lingual", "cuneus"]

# Project root, resolved from this file (code/analysis/04_brain_spatial_enrichment/).
PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))

# Everything AHBA-derived lives together: the region x gene matrix, the atlas metadata,
# the Dear et al. gradient scores and the two DK .annot label files. Sets the default for
# every data_dir / annot_dir argument below, so a future move is a one-line change.
AHBA_DIR = f"{PROJ}/data/AHBA"
ANALYSIS_DIR = f"{PROJ}/outputs/analysis/04_brain_spatial_enrichment"
SPIN_CACHE = f"{ANALYSIS_DIR}/spin_cache"
ANNOT_LH, ANNOT_RH = f"{AHBA_DIR}/lh.aparc.annot", f"{AHBA_DIR}/rh.aparc.annot"


# --- 1) AHBA expression and the regional enrichment test --- #

def load_ahba(data_dir=AHBA_DIR):
    """AHBA Desikan-Killiany region x gene matrix and region metadata (built by nb 00)."""
    expr = pd.read_csv(f"{data_dir}/ahba_region_by_gene_5donor.csv", index_col=0)
    info = pd.read_csv(f"{data_dir}/ahba_atlas_info.csv").set_index("id")
    return expr, info


def brain_expressed_background(expr, quantile=0.10):
    """Genes above the 10th percentile of mean cross-regional AHBA expression.

    The gene sets tested here are already known to be neuronal, so a whole-transcriptome
    background would confound "expressed in brain" with "regionally patterned". This is
    the background used throughout the analysis.
    """
    gene_mean = expr.mean(axis=0)
    return gene_mean[gene_mean >= gene_mean.quantile(quantile)].index.tolist()


def zscore_expression(expr, background):
    """Per-gene z-score across regions, restricted to the background (regions x genes)."""
    return expr[background].apply(zscore, axis=0)


def regional_enrichment(Z, info, gene_set, n_perm=10000, seed=0):
    """Per-region enrichment of `gene_set` vs `n_perm` size-matched random background sets.

    The observed per-region score is the mean z-scored expression over the set's genes;
    the null is the same statistic for random gene sets of identical size drawn from the
    brain-expressed background. `z_vs_bg` = (observed - null mean) / null SD is the
    background-corrected enrichment that every downstream map and spin test uses. The
    empirical p-value is one-sided (enrichment) and BH-FDR corrected across region labels.
    Left/right homologues are collapsed by DK label.
    """
    in_bg = [g for g in gene_set if g in Z.columns]
    Zvals = Z.to_numpy()
    observed = Z[in_bg].mean(axis=1).to_numpy()

    rng = np.random.default_rng(seed)
    null = np.empty((n_perm, Z.shape[0]))
    for b in range(n_perm):
        null[b] = Zvals[:, rng.choice(Z.shape[1], len(in_bg), replace=False)].mean(axis=1)

    label = np.array([info.loc[i, "label"] for i in Z.index])
    structure = np.array([info.loc[i, "structure"] for i in Z.index])
    rows = []
    for L in pd.unique(label):
        m = label == L
        obs_L, null_L = observed[m].mean(), null[:, m].mean(axis=1)
        rows.append({"region": L, "structure": structure[m][0], "obs": obs_L,
                     "z_vs_bg": (obs_L - null_L.mean()) / null_L.std(),
                     "p_emp": (np.sum(null_L >= obs_L) + 1) / (n_perm + 1)})
    enr = pd.DataFrame(rows)
    enr["fdr"] = multipletests(enr["p_emp"], method="fdr_bh")[1]
    return enr.sort_values("z_vs_bg", ascending=False).reset_index(drop=True)


def cortical_vmax(enr):
    """Symmetric colour limit from the cortical regions only.

    Subcortical |z| is larger, so letting it set vmax washes out the cortical surface.
    """
    return float(np.nanmax(np.abs(enr.loc[enr.structure == "cortex", "z_vs_bg"])))


# --- 2) Desikan-Killiany geometry --- #

def dk_region_order(data_dir=AHBA_DIR):
    """The 34 left-hemisphere DK regions in the order run_spin_test.py uses.

    Taken from the Dear et al. gradient-score table so that any map built here lines up
    with the cached spin nulls without re-deriving the parcel order.
    """
    scores = pd.read_csv(f"{data_dir}/ahba_dme_scores_in_dk.csv", index_col=0)
    return [s.replace("lh_", "") for s in scores.index]


def read_dk_annot(annot_path):
    """DK vertex labels and region names, hemisphere prefix stripped so LH and RH match."""
    labels, _, names = nib.freesurfer.read_annot(annot_path)
    names = [n.decode() if isinstance(n, bytes) else n for n in names]
    return labels, [n[3:] if n[:3] in ("lh_", "rh_") else n for n in names]


_RESAMPLE_INDEX = {}


def _index_164k_to(density, hemi):
    """Vertex index mapping native 164k -> a coarser fsaverage density, cached.

    fsaverage spheres are nested, so every coarse vertex coincides exactly with a 164k
    vertex and takes its label: the resample is lossless, not interpolated.
    """
    key = (density, hemi)
    if key not in _RESAMPLE_INDEX:
        unit = lambda x: x / np.linalg.norm(x, axis=1, keepdims=True)
        src = datasets.fetch_surf_fsaverage("fsaverage")[f"sphere_{hemi}"]
        tgt = datasets.fetch_surf_fsaverage(density)[f"sphere_{hemi}"]
        _RESAMPLE_INDEX[key] = cKDTree(unit(load_surf_mesh(src).coordinates)) \
            .query(unit(load_surf_mesh(tgt).coordinates))[1]
    return _RESAMPLE_INDEX[key]


def ap_position(data_dir=AHBA_DIR):
    """Mean pial Y coordinate per DK region: the anterior-posterior axis, in mm.

    Reconstructed from the surface (rather than read from a table) so the scatter axis
    matches the A-P map that run_spin_test.py actually spin-tests.
    """
    fs5 = datasets.fetch_surf_fsaverage("fsaverage5")
    labels, names = read_dk_annot(f"{data_dir}/lh.aparc.annot")
    labels = labels[_index_164k_to("fsaverage5", "left")]
    pial_y = load_surf_mesh(fs5["pial_left"]).coordinates[:, 1]
    return {names[c]: pial_y[labels == c].mean() for c in np.unique(labels) if c >= 0}


# --- 3) Surface maps --- #

def diverging_cmap(high_color=None):
    """Colormap for enrichment Z, with out-of-range values routed to the medial wall.

    `high_color=None` gives RdBu_r (used for the all-ASD map); passing a GO-topic colour
    gives black -> white -> topic, so zero stays white and the enriched end carries topic
    identity while the depleted end is a shared hue comparable across topics.
    """
    cmap = (plt.get_cmap("RdBu_r").copy() if high_color is None else
            LinearSegmentedColormap.from_list("topic", ["#000000", "#ffffff", high_color]).copy())
    cmap.set_over(MEDIAL_WALL)
    return cmap


def plot_dk_surface(z_by_region, cmap, vmax, outline=(), outline_color="black",
                    density=DENSITY, hemi="right", views=("lateral", "medial"),
                    figsize=(4.4, 2.4), title=None, title_color="black",
                    clabel="enrichment Z", bare=False, annot_dir=AHBA_DIR):
    """Paint per-region values on an inflated fsaverage surface, one axis per view.

    Regions named in `outline` are traced (used to mark spin-significant regions).
    No-data vertices are pushed to an over-range sentinel so the colormap's `set_over`
    paints them MEDIAL_WALL rather than leaving a hole in the scale. `bare=True` drops
    the title, colorbar and margins, giving a clean panel for figure assembly.
    """
    surf = datasets.fetch_surf_fsaverage(density)
    labels, names = read_dk_annot(f"{annot_dir}/{'lh' if hemi == 'left' else 'rh'}.aparc.annot")
    if density != "fsaverage":
        labels = labels[_index_164k_to(density, hemi)]

    vertex_z = np.full(labels.shape, np.nan)
    for code, nm in enumerate(names):
        if nm in z_by_region:
            vertex_z[labels == code] = z_by_region[nm]
    sentinel = np.where(np.isnan(vertex_z), vmax * 10, vertex_z)

    fig, axes = plt.subplots(1, len(views), figsize=figsize, subplot_kw={"projection": "3d"})
    axes = np.atleast_1d(axes)
    for view, ax in zip(views, axes):
        plotting.plot_surf_stat_map(surf[f"infl_{hemi}"], sentinel, hemi=hemi, view=view,
                                    cmap=cmap, vmax=vmax, bg_map=surf[f"sulc_{hemi}"],
                                    bg_on_data=True, axes=ax, colorbar=False)
        codes = [c for c, nm in enumerate(names) if nm in outline]
        if codes:
            plotting.plot_surf_contours(surf[f"infl_{hemi}"], labels, levels=codes,
                                        colors=[outline_color] * len(codes), axes=ax, figure=fig)
    if bare:
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=-0.10)
        return fig
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(-vmax, vmax)); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.03, shrink=0.55, aspect=22)
    cb.set_label(clabel, fontsize=6); cb.ax.tick_params(labelsize=5)
    if title:
        fig.suptitle(title, y=0.70, fontsize=8, color=title_color)
    return fig


def colorbar_panel(cmap, vmax, label="Z", figsize=(0.9, 0.8), bar_cm=(0.15, 0.30)):
    """Standalone vertical colorbar sized in cm, for assembling main-text figures."""
    cm2in = 1 / 2.54
    fig = plt.figure(figsize=figsize)
    cax = fig.add_axes([0.30, 0.28, bar_cm[0] * cm2in / figsize[0], bar_cm[1] * cm2in / figsize[1]])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(-vmax, vmax)); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="vertical")
    cb.outline.set_linewidth(0.3)
    cb.set_ticks([-round(vmax), 0, round(vmax)])
    cb.ax.tick_params(labelsize=4, length=1, width=0.3, pad=1)
    cb.set_label(label, fontsize=5, labelpad=2, rotation=0)
    return fig


# --- 4) Empirical-null helpers --- #

def masked_spearman(x, y):
    """Spearman over parcels where both maps are finite.

    Spun maps carry NaN wherever a rotated parcel lands on the medial wall; dropping
    those rather than imputing them is how neuromaps' own comparisons handle it.
    """
    m = np.isfinite(x) & np.isfinite(y)
    return spearmanr(x[m], y[m])[0] if m.sum() >= 3 else np.nan


def emp_p(observed, null, two_sided=True):
    """Empirical p with the observed value included in the null (never returns 0)."""
    null = np.asarray(null)[np.isfinite(null)]
    if two_sided:
        return (np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)
    return (np.sum(null >= observed) + 1) / (len(null) + 1)


def pstr(p):
    """p-value for a figure label: 3 decimals, or m x 10^-n mathtext below 1e-3."""
    if p >= 1e-3:
        return f"{p:.3f}"
    mantissa, exponent = f"{p:.0e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def plot_spin_null(ax, null, observed, color, title=None, bins=50, xlim=None, lw=1.25):
    """Spin-null distribution with the observed statistic marked."""
    null = np.asarray(null)[np.isfinite(null)]
    ax.hist(null, bins=bins, color="0.90", edgecolor="black", linewidth=0.05, density=True)
    ax.axvline(0, color="0.6", lw=0.5, ls=":")
    ax.axvline(observed, color=color, lw=lw)
    if title:
        ax.set_title(title, fontsize=5.5, color="black")
    ax.set_yticks([])
    if xlim:
        ax.set_xlim(*xlim)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    return ax


# --- 5) Occipital sensitivity analyses --- #

def occipital_gene_scores(Z, info, regions=None):
    """Per-gene occipital score: each gene's z-expression averaged over the occipital regions.

    Region-balanced within each occipital label, so a label sampled in both hemispheres
    does not count twice. A gene set's occipital enrichment is then just the mean of this
    series over the set -- algebraically identical to averaging the per-region set score
    over the four labels, but O(n) per call, which is what makes the 5,000-draw matched
    resampling run in seconds.
    """
    regions = OCCIPITAL if regions is None else regions
    label = np.array([info.loc[i, "label"] for i in Z.index])
    weights = np.zeros(Z.shape[0])
    for L in regions:
        rows = np.where(label == L)[0]
        if len(rows):
            weights[rows] = 1.0 / (len(regions) * len(rows))
    return pd.Series(weights @ Z.to_numpy(), index=list(Z.columns))


def matched_null_occipital(gene_set, gene_occ, background, covariates, n_perm=5000, seed=0):
    """Occipital enrichment vs a null matched on expression, CDS length and constraint.

    The brain-expressed background does not control for the fact that these genes are
    long, constrained and highly expressed -- properties that covary with cortical
    expression topography. Here the null sets are matched on the JOINT distribution of
    those three covariates (tertile bins on each, crossed into up to 27 strata, drawing
    the tested set's count from each stratum). An observed value still exceeding the
    matched null means an effect beyond those gene properties. The unmatched null is
    returned alongside so the two can be shown together.

    `covariates` is a DataFrame indexed by gene with columns expr, cds, LOEUF.
    """
    binned = covariates[["expr", "cds", "LOEUF"]].apply(
        lambda c: pd.qcut(c, 3, labels=False, duplicates="drop"))
    strata = binned["expr"].astype(str) + "_" + binned["cds"].astype(str) + "_" + binned["LOEUF"].astype(str)
    by_stratum = {s: strata.index[strata == s].tolist() for s in strata.unique()}

    in_set = [g for g in gene_set if g in gene_occ.index]
    observed = gene_occ.loc[in_set].mean()
    counts = strata.loc[[g for g in in_set if g in covariates.index]].value_counts()
    rng = np.random.default_rng(seed)

    def draw_matched():
        out = []
        for s, k in counts.items():
            out += list(rng.choice(by_stratum[s], min(k, len(by_stratum[s])), replace=False))
        return out

    score = lambda genes: gene_occ.loc[[g for g in genes if g in gene_occ.index]].mean()
    null_m = np.array([score(draw_matched()) for _ in range(n_perm)])
    null_u = np.array([score(list(rng.choice(background, len(in_set), replace=False)))
                       for _ in range(n_perm)])
    return dict(n=len(in_set), obs=observed, null_u=null_u, null_m=null_m,
                p_u=(np.sum(null_u >= observed) + 1) / (n_perm + 1),
                p_m=(np.sum(null_m >= observed) + 1) / (n_perm + 1))


def dropocc_ap_test(enrichment_map, ap, ap_null, dk_regions, occipital=None):
    """A-P spin correlation with and without the four occipital regions.

    The A-P gradient and the occipital focal test are the same signal seen two ways;
    dropping the occipital regions asks whether the gradient is broader than occipital.
    Returns rho and spin p for the full 34 regions and for the 30 remaining.
    """
    occipital = OCCIPITAL if occipital is None else occipital
    occ_idx = [dk_regions.index(r) for r in occipital]
    keep = np.array([i for i in range(len(dk_regions)) if i not in occ_idx])

    def test(idx):
        obs = spearmanr(enrichment_map[idx], ap[idx])[0]
        null = np.array([spearmanr(enrichment_map[idx], ap_null[idx, s])[0]
                         for s in range(ap_null.shape[1])])
        return obs, emp_p(obs, null)

    rho_full, p_full = test(np.arange(len(dk_regions)))
    rho_drop, p_drop = test(keep)
    return dict(rho_full=rho_full, p_full=p_full, rho_dropocc=rho_drop, p_dropocc=p_drop)


def drop_top_k(gene_set, gene_occ, max_k=30):
    """Occipital enrichment of the remainder after removing the top k occipital genes.

    Leave-one-out only rules out a single dominant gene. Ranking the set by per-gene
    occipital score and removing the top k asks whether a *handful* carries the effect:
    staying above chance after ~10-20 removals means the signal is distributed.
    Returns (ks, enrichment, n_genes).
    """
    values = np.sort(gene_occ.loc[[g for g in gene_set if g in gene_occ.index]].values)[::-1]
    ks = np.arange(0, min(max_k, len(values) - 5) + 1)
    return ks, np.array([values[k:].mean() for k in ks]), len(values)


def map_nonuniformity(Z, info, gene_set, n_perm=20000, seed=0):
    """Is the map spatially structured, or flat? Spatial variance vs random gene sets.

    For the observed set and for `n_perm` random size-matched sets, the per-region
    enrichment Z (background-corrected, the same quantity as z_vs_bg) is computed across
    the 34 cortical regions and its variance taken -- one number for how spatially spread
    out the map is. This is a gene-set resampling test, not a spin test: it asks whether
    the set has *any* regional structure, not whether that structure aligns with an axis.
    Returns (observed regional Z, observed variance, null variances).
    """
    label = np.array([info.loc[i, "label"] for i in Z.index])
    structure = np.array([info.loc[i, "structure"] for i in Z.index])
    ctx_labels = [L for L in pd.unique(label) if structure[label == L][0] == "cortex"]
    Zvals = Z.to_numpy()
    col_of = {g: i for i, g in enumerate(Z.columns)}

    def collapsed(gene_idx):
        per_region = Zvals[:, gene_idx].mean(1)
        return np.array([per_region[label == L].mean() for L in ctx_labels])

    idx = np.array([col_of[g] for g in gene_set if g in col_of])
    rng = np.random.default_rng(seed)
    nulls = np.array([collapsed(rng.choice(Z.shape[1], len(idx), replace=False))
                      for _ in range(n_perm)])
    mu, sd = nulls.mean(0), nulls.std(0)
    obs_z = (collapsed(idx) - mu) / sd
    return obs_z, obs_z.var(), ((nulls - mu) / sd).var(axis=1)


def cached_ap_spin_null(n_spin=1000, density="41k", seed=0, cache_dir=SPIN_CACHE,
                        data_dir=AHBA_DIR):
    """Spun anterior-posterior axis (34 regions x n_spin), generated once and cached.

    Used by the drop-occipital A-P test in notebooks 01 and 02. Same Cornblath parcel-spin
    null as run_spin_test.py, so the drop-occipital p-values are on the same footing as
    the primary ones. Returns (ap, ap_null) with `ap` the observed A-P axis.
    """
    from run_spin_test import resample_parc_to_density
    from neuromaps.datasets import fetch_atlas
    from neuromaps.images import annot_to_gifti
    from neuromaps.nulls import cornblath
    from neuromaps.nulls.spins import get_parcel_centroids, gen_spinsamples, vertices_to_parcels

    dk = dk_region_order(data_dir)
    fs = fetch_atlas("fsaverage", density)
    parc = resample_parc_to_density(
        annot_to_gifti((f"{data_dir}/lh.aparc.annot", f"{data_dir}/rh.aparc.annot")), density)
    y_vert = np.concatenate([np.asarray(nib.load(str(fs["pial"][h])).darrays[0].data)[:, 1]
                             for h in (0, 1)])
    ap = np.asarray(vertices_to_parcels(y_vert, parc)).ravel()[:len(dk)]

    path = f"{cache_dir}/apnull_{density}_n{n_spin}_seed{seed}.npy"
    if os.path.exists(path):
        return ap, np.load(path)
    coords, hemi = get_parcel_centroids(fs["sphere"], method="surface")
    spins = gen_spinsamples(coords, hemi, n_rotate=n_spin, seed=seed, verbose=0)
    ap68 = pd.Series(ap, index=range(1, len(dk) + 1)).reindex(range(1, 69)).to_numpy()
    ap_null = cornblath(data=ap68, atlas="fsaverage", density=density, parcellation=parc,
                        n_perm=n_spin, spins=spins)[:34]
    os.makedirs(cache_dir, exist_ok=True)
    np.save(path, ap_null)
    return ap, ap_null
