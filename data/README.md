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
- `SnMultiome_Wang2025/obj_rna_raw_EN_slingshot.h5ad`
- `SnMultiome_Wang2025/obj_rna_raw_IN_slingshot.h5ad`
- `SnMultiome_Wang2025/donorID_sex_ngene.cov`
- `SnMultiome_Wang2025/metadata.csv`
- `SynGO/gene_topic_HotNet_category.xlsx`
- `RVAS_result/2026.01.07_ASD-Zscore-FDRtransform.RDS`
- `RVAS_result/counts_asd_ddid.xlsx` / `counts_asd_noddid.xlsx`
- `RVAS_result/Kaplanis_DDD_ASD_counts_by_gene_new_mis_cats_2025-04-09.txt`
- `RVAS_result/full_results_wcounts_2025-10-08.txt` — full results with variant counts; used by `figures/supp_fig/03_plot_regional_enrichment.R`
- `RVAS_result/full_results_fdr001.txt` — genes at FDR < 0.001; used by `analysis/04_ddid_stratification/01_input_prep/01_classify_ddid_gene_groups.R`
- `RVAS_result/full_results_fdr01.txt` — genes at FDR < 0.01; used by same script

## External reference tools

- MAGMA gene annotation (GRCh37): `/mnt/isilon/gandal_lab/liaoyd/tools/magma/NCBI37.3.gene.loc`
- MAGMA gene annotation (GRCh38): `/mnt/isilon/gandal_lab/liaoyd/tools/magma/NCBI38.gene.loc`
- Allen Brain Atlas taxonomy: `/mnt/isilon/gandal_lab/liaoyd/data/allen_brain_atlas/WHB_taxonomy/science.add7046_table_s3.xlsx`

