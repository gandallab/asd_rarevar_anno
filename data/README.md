# `data/`

Source inputs for the ASD rare-variant functional annotation analyses. Everything
versioned here is external and published — supplementary tables, reference atlases and
curated gene lists. Anything derived by our own code lives in `outputs/analysis/`, not
here.

Script numbers below (e.g. `04_01`) refer to the file-name prefixes of the notebooks in
`analysis/` modules.

---

## Not included in this repository

**Rare-variant statistics and de novo mutation counts are not distributed with this
repo** — they are not ours to publish. Obtain them from the original sources and place
them at the paths below to reproduce the analyses.

| Path | Contents | Source |
|---|---|---|
| `full_results_wcounts_2025-10-08.txt` | TADA gene-level results: BF, FDR, LOEUF, and per-gene proband/sibling DNM counts (PTV, Mis2, Mis1, Del, Dup) | Autism Sequencing Consortium |
| `DDID/counts_asd_ddid.xlsx` | Per-gene DNM counts, ASD probands **with** developmental delay / intellectual disability | Autism Sequencing Consortium |
| `DDID/counts_asd_noddid.xlsx` | Per-gene DNM counts, ASD probands **without** DD/ID | Autism Sequencing Consortium |
| `DDID/Kaplanis_DDD_ASD_counts_by_gene_new_mis_cats_2025-04-09.txt` | Per-gene DNM counts from the Deciphering Developmental Disorders study, sex-stratified | Kaplanis et al., DDD |

Used by: `00_01` (formats the TADA table into `00_01_ASD-rarevar-stats.RDS`, which most
downstream analyses read instead of the raw file), `04_01`, `04_02`, `05_00`, `05_01`,
`05_03`, `07_00`, `07_01`.

A few large files are also excluded for size rather than licensing — see `.gitignore`:
`BrainSpan/expression_matrix.csv`, `Gandal2022/TableS3.xlsx`,
`networks/Gandal2022_table_s5.xlsx`, `networks/Nano2025_table_s1-33.xlsx`, and
`AHBA/abagen_data/`.

The single-cell data for the scDRS/trajectory/AUCell/DD-ID modules is also not posted
here — see **Single-cell data (lab HPC)** below.

---

## Contents

### `AHBA/`
Allen Human Brain Atlas microarray, processed to a Desikan-Killiany region × gene matrix.

- `ahba_region_by_gene_5donor.csv` — 83 regions × ~15,500 genes, 5 donors
- `ahba_atlas_info.csv` — DK region metadata (id, label, structure, hemisphere)
- `ahba_dme_scores_in_dk.csv` — Dear et al. C1/C2/C3 transcriptomic gradient scores on DK;
  also fixes the 34-region parcel order the spin tests use
- `lh.aparc.annot`, `rh.aparc.annot` — FreeSurfer DK labels on fsaverage (164k)
- `abagen_data/` — raw donor archives (~6 GB, **not versioned**); re-downloaded by `04_00`

Used by: `04_00` (builds the matrix), `04_01`, `04_02`, `04_05`, `dk_atlas_schematic`,
`run_spin_test.py`. The three notebooks reach it through `spatial_helpers.AHBA_DIR`.

### `BrainSpan/`
Developmental bulk RNA-seq (`expression_matrix.csv` plus row/column metadata and the
BrainSpan `readme.txt`). Used by `03_04` for the prenatal expression trajectories.

### `Gandal2022/`
Cortical transcriptomic dysregulation in ASD.

- `TableS4.xlsx` — pairwise attenuation of regional identity (ARI), both the permutation
  and bootstrap metrics
- `TableS3.xlsx` — per-region ASD-vs-control differential expression (57 MB, **not
  versioned**)

Used by: `04_02` §6 (convergence of risk-gene topography with transcriptomic attenuation).

### `SynGO/`
`syngo_genes.xlsx`, `syngo_annotations.xlsx`, `syngo_ontologies.xlsx` — the SynGO synaptic
gene ontology release. Used by `05_00`, `05_01`, `05_02`, `05_03`, `06_00`, `06_02`,
`06_03`, `07_01`, and `figures/fig03/Fig3.qmd`.

### `Tsyporin2026/`
- `Supplementary_Table_5_GMs_shared.xlsx` — sensorimotor (S) and association (A) shared
  gene modules defining the areal-patterning axis
- `MORPH_patterning_genes.xlsx` — patterning genes highlighted in the text

Used by `04_03` (module enrichment for S/A areal-identity genes).

### `Wang2025/`
`Wang2025_type_markers.csv`, `Wang2025_background_genes.csv`, `Wang2025_table_s3.xlsx` —
cell-type marker sets. Used by `03_01_AUCell.R`.

### `gene_annotations/`
Published gene→category membership lists. Distinct from `networks/`: these assign genes to
a functional category, they do not encode gene-gene relationships.

- `TFs_Ensembl_v_1.01.txt` — transcription factors, AnimalTFDB v1.01 (Ensembl IDs)
- `kepi_a_2139067_sm6382(1).xlsx` — epigenetic regulators, Engelen et al. 2023
  (doi:10.1080/15548627.2022.2139067), Supplementary Table S2

Used by `05_00` (mutation-rate-matched enrichment vs a brain-expressed background).

### `networks/`
Published co-expression (WGCNA) and gene regulatory network module definitions — the
networks the whole project annotates. `Gandal2022_table_s5.xlsx`,
`Gandal2022_table_s6.xlsx`, `Wamsley2024_table_s5.xlsx`, `Wang2025_table_s13.xlsx`,
`Wen2024_table_s6.xlsx`. Used by `00_00` to build `00_00_gene-networks-all.RDS`.

### `proteomics/`
`synaptic_proteomics_enrichments.xlsx` — synaptic proteomics enrichment results. Used by
`figures/fig03/Fig3.qmd`.

### Loose files
- `gene_cds_length.csv` — union CDS length per gene, used as a covariate in the
  gene-property-matched nulls in `04_01` and `04_02`

---

## Single-cell data (lab HPC)

The single-cell inputs for the scDRS/trajectory/AUCell/DD-ID modules are not stored in
this repository. They live on the republica HPC:

- Project data root: `/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/`
- Project result root: `/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/`

Input data locations:

- Single-cell atlas (snMultiome Wang 2025): `.../data/SnMultiome_Wang2025/`
- Rare variant results: `.../data/RVAS_result/`
- GWAS and MAGMA preparation: `.../data/GWAS_result/`
- SynGO reference: `.../data/SynGO/`

Frequently referenced files:

- `SnMultiome_Wang2025/obj_rna_raw.h5ad`
- `SnMultiome_Wang2025/obj_rna_raw_EN_slingshot.h5ad` / `obj_rna_raw_IN_slingshot.h5ad`
- `SnMultiome_Wang2025/obj_rna_raw_EN_slingshot_paired.h5ad` / `obj_rna_raw_IN_slingshot_paired.h5ad` — EN/IN lineage objects restricted to the fixed 11-donor paired (PFC + V1) cohort; written by the "mclust_Group_region, paired-donor cohort only" section of `analysis/01_scdrs_trajectory/01_scdrs/01_input_prep/01_2_prepare_downstream_metadata.ipynb`; consumed by `analysis/01_scdrs_trajectory/01_scdrs/03_group_tests/04_run_mclust_group_region_analysis_paired.sh`
- `SnMultiome_Wang2025/obj_rna_raw_subclass_donor_region_paired.h5ad` — paired-cohort EN/IN cells pooled across all types/ages, keyed by broad subclass x donor x region; written by the "subclass_donor_region" section of the same notebook; consumed by `analysis/01_scdrs_trajectory/01_scdrs/03_group_tests/05_run_subclass_donor_region_analysis.sh`
- `SnMultiome_Wang2025/donorID_sex_ngene.cov`
- `SnMultiome_Wang2025/metadata.csv`
