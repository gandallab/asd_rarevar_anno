#!/usr/bin/env python
"""
Spin test: ASD gene-set enrichment map vs the A-P axis and the AHBA transcriptomic
gradients C1/C2/C3, using Dear et al.'s neuromaps `cornblath` parcel-spin
(maps_null_test.py) on the Desikan-Killiany atlas.

The expensive steps -- generating the surface spins, null-projecting the components,
and spinning the ASD map for the focal test -- are all cached. The spins and component
nulls depend only on the atlas/gradients/A-P axis (not the ASD map), so they are reused
across any maps you test; the focal ASD-map null is keyed by the map's values, so it is
reused across ROIs and reruns but regenerated if the map changes. Re-running with the
same --density/--n-spin/--seed loads the caches instead of regenerating (--force skips
them). At n_spin=5000 this means iterating on --roi-regions is essentially free.

In addition to the gradient/A-P correlations, the ASD map itself is spun (same
cornblath spins) to test **focal ROI enrichment**: the observed ASD score over an ROI
of DK region(s) vs its distribution under map rotation. --roi-mode union tests the
mean over a set of regions (default = the 4 occipital regions, reproducing the old
occipital test); --roi-mode separate runs one spin test per region (spin-corrected
per-region enrichment). This keeps a single industry-standard
(Alexander-Bloch/Vasa/cornblath) null throughout.

Run from this directory (which holds maps_null_test.py). Default paths are resolved from
this file's location: data from `<project>/data/AHBA/`, results and caches under
`<project>/outputs/analysis/04_brain_spatial_enrichment/`, e.g.:
    python run_spin_test.py --n-spin 1000
    python run_spin_test.py --n-spin 5000 --density 41k      # 41k needs wb_command
    python run_spin_test.py --n-spin 1000 --force            # ignore cache, regenerate

Requires: neuromaps, statsmodels, nibabel. No Connectome Workbench (wb_command) needed:
the DK annot is resampled to the spin density by sphere nearest-neighbour (see below).
"""
import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd
import nibabel as nib


def resample_parc_to_density(parc, target_density, src_density="164k"):
    """Nearest-neighbour resample DK label giftis from src_density to target_density
    using fsaverage sphere coordinates.

    Reproduces Dear et al.'s fsaverage_to_fsaverage(method="nearest") without Connectome
    Workbench. fsaverage spheres are nested, so every target-density vertex coincides
    with a source-density vertex and takes its label exactly -- identical to
    `wb_command -label-resample ... NEAREST`.
    """
    from neuromaps.datasets import fetch_atlas
    from scipy.spatial import cKDTree
    src = fetch_atlas("fsaverage", src_density)["sphere"]
    tgt = fetch_atlas("fsaverage", target_density)["sphere"]
    unit = lambda x: x / np.linalg.norm(x, axis=1, keepdims=True)
    out = []
    for h in (0, 1):
        s = unit(np.asarray(nib.load(str(src[h])).darrays[0].data))
        t = unit(np.asarray(nib.load(str(tgt[h])).darrays[0].data))
        idx = cKDTree(s).query(t)[1]
        labels = np.asarray(parc[h].darrays[0].data)[idx].astype(np.int32)
        g = nib.gifti.GiftiImage()
        g.add_gifti_data_array(nib.gifti.GiftiDataArray(labels, intent="NIFTI_INTENT_LABEL"))
        out.append(g)
    return tuple(out)

# make `import maps_null_test` work regardless of the cwd the job is submitted from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# project root, so the defaults below hold wherever the job is launched from
PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
AHBA = f"{PROJ}/data/AHBA"
ANALYSIS = f"{PROJ}/outputs/analysis/04_brain_spatial_enrichment"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-spin", type=int, default=1000,
                   help="number of spin permutations (Dear et al. used 500; 1000+ gives tighter p)")
    p.add_argument("--density", default="41k",
                   help="fsaverage density for the spin surface. 41k = Dear et al.'s density "
                        "(the default). The DK annot is native 164k and is resampled to this "
                        "density by nearest-neighbour on the sphere (no wb_command needed); "
                        "pass 164k to skip resampling.")
    p.add_argument("--seed", type=int, default=0, help="RNG seed for the spins (reproducibility)")
    p.add_argument("--scores", default=f"{AHBA}/ahba_dme_scores_in_dk.csv",
                   help="Dear et al. gradient scores on DK (index = lh_<region>; cols C1,C2,C3)")
    p.add_argument("--asd-map", default=f"{ANALYSIS}/04_01_ASD_regional_enrichment.csv",
                   help="ASD regional enrichment table (must contain a 'region' column + --asd-col)")
    p.add_argument("--asd-cols", nargs="+", default=["z_vs_bg", "obs"],
                   help="one or more columns of --asd-map to spin-test as maps, in one job "
                        "(z_vs_bg = background-corrected enrichment z; obs = raw set-score)")
    p.add_argument("--annot-lh", default=f"{AHBA}/lh.aparc.annot")
    p.add_argument("--annot-rh", default=f"{AHBA}/rh.aparc.annot")
    p.add_argument("--roi-regions", nargs="+",
                   default=["pericalcarine", "lateraloccipital", "lingual", "cuneus"],
                   help="DK region label(s) for the focal enrichment test: a single region or a "
                        "list. Pass 'all' for every cortical region (useful with --roi-mode "
                        "separate to scan them). Default = the 4 occipital regions.")
    p.add_argument("--roi-mode", choices=["union", "separate"], default="union",
                   help="union: one test on the mean over all --roi-regions (default; the "
                        "occipital-cluster test with the default regions). separate: one spin test "
                        "per region -> spin-corrected per-region enrichment (FDR across regions).")
    p.add_argument("--roi-name", default="occipital",
                   help="label for the ROI in the output rows and output filename")
    p.add_argument("--method", default="spearmanr", help="correlation metric (spearmanr | pearsonr)")
    p.add_argument("--outdir", default=ANALYSIS)
    p.add_argument("--prefix", default="",
                   help="prefix for the output filenames, identifying the notebook and gene set "
                        "that the tested map came from (e.g. 04_01_ASD_, 04_02_GR_). Results from every "
                        "gene set share one flat --outdir, so this is what keeps them apart.")
    p.add_argument("--cachedir", default=f"{ANALYSIS}/spin_cache",
                   help="where spins/nulls are cached and reused across runs")
    p.add_argument("--force", action="store_true",
                   help="regenerate spins and nulls even if a matching cache exists")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.cachedir, exist_ok=True)

    import maps_null_test as mnt
    from neuromaps.images import annot_to_gifti
    from neuromaps.datasets import fetch_atlas
    from neuromaps.nulls import cornblath
    from neuromaps.nulls.spins import get_parcel_centroids, gen_spinsamples, vertices_to_parcels

    # --- surfaces + DK parcellation on the spin surface ------------------------ #
    fs = fetch_atlas("fsaverage", args.density)
    parc = annot_to_gifti((args.annot_lh, args.annot_rh))          # native 164k
    if args.density != "164k":
        # Resample the DK labels 164k -> spin density by nearest-neighbour on the sphere.
        # This reproduces Dear et al.'s fsaverage_to_fsaverage(method="nearest") WITHOUT
        # Connectome Workbench: fsaverage spheres are nested, so each target vertex
        # coincides with a source vertex and inherits its label exactly (identical to
        # wb_command -label-resample ... NEAREST).
        parc = resample_parc_to_density(parc, args.density)

    # --- components to spin: Dear's C1/C2/C3 + the anterior-posterior (Y) axis -- #
    scores_raw = pd.read_csv(args.scores, index_col=0)
    dk_regions = [s.replace("lh_", "") for s in scores_raw.index]
    scores = scores_raw.reset_index(drop=True)
    scores.index = range(1, len(scores) + 1)                        # integer 1..34
    y_vert = np.concatenate([nib.load(str(fs["pial"][h])).darrays[0].data[:, 1] for h in (0, 1)])
    ap = np.asarray(vertices_to_parcels(y_vert, parc)).ravel()      # mean Y per parcel (68); LH first
    scores["AP"] = ap[:len(dk_regions)]
    n_comp = scores.shape[1]
    print(f"components ({n_comp}): {list(scores.columns)} | regions: {len(scores)}")

    # --- ASD map(s), aligned to the same regions/order ------------------------ #
    # one column per requested source column (e.g. z_vs_bg = background-corrected,
    # obs = raw set-score), so all are spin-tested against every component in one job.
    enr_df = pd.read_csv(args.asd_map).set_index("region")
    missing = [c for c in args.asd_cols if c not in enr_df.columns]
    if missing:
        sys.exit(f"ERROR: columns not in {args.asd_map}: {missing}. Available: {list(enr_df.columns)}")
    asd = pd.DataFrame(
        {f"ASD_{col}": [enr_df[col].get(r, np.nan) for r in dk_regions] for col in args.asd_cols},
        index=range(1, len(dk_regions) + 1))
    for col in asd.columns:
        n_match = int(asd[col].notna().sum())
        print(f"{col}: regions matched {n_match}/{len(dk_regions)}")
        if n_match != len(dk_regions):
            sys.exit(f"ERROR: {col} label mismatch -- check region naming vs the gradient regions.")

    tag = f"{args.density}_n{args.n_spin}_seed{args.seed}"

    # --- spins (cached; depend only on density/n/seed) ------------------------- #
    spins_path = os.path.join(args.cachedir, f"spins_{tag}.npy")
    if os.path.exists(spins_path) and not args.force:
        spins = np.load(spins_path)
        print(f"[cache] loaded spins {spins.shape} <- {spins_path}")
    else:
        print(f"generating {args.n_spin} spins @ {args.density} (seed {args.seed}) ...")
        coords, hemiid = get_parcel_centroids(fs["sphere"], method="surface")
        spins = gen_spinsamples(coords, hemiid, n_rotate=args.n_spin, seed=args.seed, verbose=1)
        np.save(spins_path, spins)
        print(f"[cache] saved spins {spins.shape} -> {spins_path}")

    # --- component nulls (cached; depend on scores + spins, not the ASD map) --- #
    nulls_path = os.path.join(args.cachedir, f"nulls_{tag}_{n_comp}comp.npy")
    if os.path.exists(nulls_path) and not args.force:
        nulls = np.load(nulls_path)
        print(f"[cache] loaded nulls {nulls.shape} <- {nulls_path}")
    else:
        print("projecting component nulls (cornblath) ...")
        base = os.path.dirname(nulls_path)
        name = os.path.basename(nulls_path)[:-len(".npy")]
        mnt.generate_nulls_from_components(scores, spins, atlas="dk", parcellation_img=parc,
                                           density=args.density, n=args.n_spin,
                                           save_dir=base + os.sep, save_name=name)
        nulls = np.load(nulls_path)

    # --- correlate the ASD map against each component ------------------------- #
    # keep the FDR column for transparency; report whichever suits your plan
    # (e.g. A-P as the pre-specified primary test, C1/C2/C3 as secondary).
    res = mnt.correlate_maps_with_null_scores(nulls, scores, asd, method=args.method)
    out_csv = os.path.join(args.outdir, f"{args.prefix}spin_results_{tag}.csv")
    res.to_csv(out_csv, index=False)
    print("\n" + res.to_string(index=False))
    print(f"\nsaved -> {out_csv}")

    # --- focal ROI enrichment: spin the ASD map itself ------------------------- #
    # A location test (not a gradient correlation): is the ASD score over an ROI of DK
    # region(s) higher than under map rotation? Each ASD map is spun with the SAME
    # cornblath spins, and the observed ROI statistic is compared to its rotated
    # distribution (one-sided empirical p; enrichment = higher than chance). This keeps
    # the focal result under the single Alexander-Bloch/Vasa/cornblath null.
    #   union    : one test on the MEAN over all --roi-regions (default; the occipital-
    #              cluster test with the default 4 regions).
    #   separate : one test PER region -> spin-corrected per-region enrichment; FDR is
    #              applied across the tested regions, within each map.
    roi_regions = list(dk_regions) if args.roi_regions == ["all"] else args.roi_regions
    missing_roi = [r for r in roi_regions if r not in dk_regions]
    if missing_roi:
        sys.exit(f"ERROR: --roi-regions not in the {len(dk_regions)} DK regions: {missing_roi}")
    roi_idx = [dk_regions.index(r) for r in roi_regions]

    focal_rows = []
    for col in asd.columns:
        data68 = asd[col].reindex(range(1, 69)).to_numpy()          # right hemi -> NaN
        # Cache the spun ASD map (68 x n_spin). It depends only on the map values + spins,
        # NOT on the ROI or --roi-mode, so iterating on --roi-regions reuses it. Keyed by a
        # hash of the map values (auto-invalidates if the map content changes) plus the tag.
        maphash = hashlib.md5(np.ascontiguousarray(data68, dtype=np.float64).tobytes()).hexdigest()[:10]
        fnull_path = os.path.join(args.cachedir, f"focalnull_{tag}_{col}_{maphash}.npy")
        if os.path.exists(fnull_path) and not args.force:
            null_map = np.load(fnull_path)
            print(f"[cache] loaded focal null {null_map.shape} ({col}) <- {fnull_path}")
        else:
            null_map = cornblath(data=data68, atlas="fsaverage", density=args.density,
                                 parcellation=parc, n_perm=args.n_spin, spins=spins)   # 68 x n
            np.save(fnull_path, null_map)
            print(f"[cache] saved focal null {null_map.shape} ({col}) -> {fnull_path}")
        lh_null = null_map[:34]                                      # LH parcels x spins
        vals = asd[col].to_numpy()
        # NaN-aware empirical p: a spin whose ROI lands entirely on the medial wall
        # returns NaN for that parcel (no data). Ignore those parcels/spins rather than
        # counting them as "not exceeding" (which would bias p downward). Mirrors how
        # Dear's compare_images drops NaNs in the correlation tests.
        def emp_p(obs, null_stat):
            valid = np.isfinite(null_stat)
            return (np.sum(null_stat[valid] >= obs) + 1) / (valid.sum() + 1)
        if args.roi_mode == "union":
            obs = np.nanmean(vals[roi_idx])
            null_stat = np.nanmean(lh_null[roi_idx], axis=0)        # ROI mean per spin (NaN parcels ignored)
            p = emp_p(obs, null_stat)
            n_valid = int(np.isfinite(null_stat).sum())
            focal_rows.append({"C": f"{args.roi_name}_focal", "map": col, "enrich_z": obs,
                               "p": p, "q": np.nan})
            print(f"{args.roi_name} focal ({col}): obs={obs:.3f} vs null mean "
                  f"{np.nanmean(null_stat):.3f} ({n_valid}/{args.n_spin} valid spins), spin p={p:.4f}")
        else:  # separate: one test per region
            # max-statistic FWE: the null of the single largest region statistic per spin.
            # A region is FWE-significant if its observed stat beats that max-null. This is
            # the permutation-native family-wise correction across the tested regions -- it
            # controls the chance of ANY false positive, no distributional assumptions.
            max_null = np.nanmax(lh_null[roi_idx], axis=0)         # max over tested regions per spin
            mvalid = np.isfinite(max_null)
            for r, i in zip(roi_regions, roi_idx):
                obs = vals[i]
                p = emp_p(obs, lh_null[i])                          # marginal (uncorrected) spin p
                p_fwe = (np.sum(max_null[mvalid] >= obs) + 1) / (mvalid.sum() + 1)
                focal_rows.append({"C": r, "map": col, "enrich_z": obs, "p": p,
                                   "p_fwe": p_fwe, "q": np.nan})

    focal_df = pd.DataFrame(focal_rows)
    if args.roi_mode == "separate":
        # per-region family: FDR (Benjamini-Hochberg) across the tested regions, per map
        from statsmodels.stats.multitest import multipletests
        focal_df["q"] = focal_df.groupby("map")["p"].transform(
            lambda pv: multipletests(pv, method="fdr_bh")[1])
        print(f"\n{args.roi_mode} per-region focal tests "
              f"({len(roi_regions)} regions x {asd.shape[1]} maps):")
        print(focal_df.to_string(index=False))

    label = f"{args.roi_name}_focal" if args.roi_mode == "union" else f"{args.roi_name}_perregion"
    focal_csv = os.path.join(args.outdir, f"{args.prefix}{label}_{tag}.csv")
    focal_df.to_csv(focal_csv, index=False)
    print(f"saved -> {focal_csv}")


if __name__ == "__main__":
    main()
