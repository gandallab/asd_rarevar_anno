# Data

No data files are stored in this repository. This directory documents where source
data lives on the republica HPC and what each path is used for.

## Data roots

- Project data root for scDRS/trajectory/aucell/ddid analysis: `/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/`
- Project result root for scDRS/trajectory/aucell/ddid analysis: `/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/`

## Input data locations

- Single-cell atlas (snMultiome Wang 2025): `.../data/SnMultiome_Wang2025/`
- Rare variant results: `.../data/RVAS_result/`
- GWAS and MAGMA preparation: `.../data/GWAS_result/`
- SynGO reference: `.../data/SynGO/`

## Frequently referenced files

- `SnMultiome_Wang2025/obj_rna_raw.h5ad`
- `SnMultiome_Wang2025/obj_rna_raw_EN_slingshot.h5ad` / `obj_rna_raw_IN_slingshot.h5ad`
- `SnMultiome_Wang2025/obj_rna_raw_EN_slingshot_paired.h5ad` / `obj_rna_raw_IN_slingshot_paired.h5ad` — EN/IN lineage objects restricted to the fixed 11-donor paired (PFC + V1) cohort; written by the "mclust_Group_region, paired-donor cohort only" section of `analysis/01_scdrs/01_input_prep/01_2_prepare_downstream_metadata.ipynb`; consumed by `analysis/01_scdrs/03_group_tests/04_run_mclust_group_region_analysis_paired.sh`
- `SnMultiome_Wang2025/obj_rna_raw_subclass_donor_region_paired.h5ad` — paired-cohort EN/IN cells pooled across all types/ages, keyed by broad subclass x donor x region; written by the "subclass_donor_region" section of the same notebook; consumed by `analysis/01_scdrs/03_group_tests/05_run_subclass_donor_region_analysis.sh`
- `SnMultiome_Wang2025/donorID_sex_ngene.cov`
- `SnMultiome_Wang2025/metadata.csv`
- `SynGO/gene_topic_HotNet_category.xlsx`
- `RVAS_result/2026.01.07_ASD-Zscore-FDRtransform.RDS`
- `RVAS_result/counts_asd_ddid.xlsx` / `counts_asd_noddid.xlsx`
- `RVAS_result/Kaplanis_DDD_ASD_counts_by_gene_new_mis_cats_2025-04-09.txt`
- `RVAS_result/full_results_wcounts_2025-10-08.txt` — full results with variant counts; used by `figures/supp_fig/05_plot_regional_enrichment.R`
- `RVAS_result/full_results_fdr001.txt` — genes at FDR < 0.001; used by `analysis/04_ddid_stratification/01_input_prep/01_classify_ddid_gene_groups.R` and `figures/supp_fig/05_plot_regional_enrichment.R` (gene universe for regional factor assignment)
- `RVAS_result/full_results_fdr01.txt` — genes at FDR < 0.01; used by same script
- `RVAS_result/enriched_genes_per_factor.txt` — long-format gene x NMF-factor relative loading + robust-expression flag (one row per gene/factor membership from Aivazidis et al); used by `figures/supp_fig/05_plot_regional_enrichment.R` to assign genes to regional factors (GZ/CP/TH)

## External reference tools

- MAGMA gene annotation (GRCh37): `/mnt/isilon/gandal_lab/liaoyd/tools/magma/NCBI37.3.gene.loc`
- MAGMA gene annotation (GRCh38): `/mnt/isilon/gandal_lab/liaoyd/tools/magma/NCBI38.gene.loc`
- Allen Brain Atlas taxonomy: `/mnt/isilon/gandal_lab/liaoyd/data/allen_brain_atlas/WHB_taxonomy/science.add7046_table_s3.xlsx`

