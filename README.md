# ASD Rare Variant Annotation Project

Research analysis project for ASD rare variant annotation and downstream single-cell
interpretation.

## System requirements

**Operating system:** tested on macOS Sequoia 15.3.1 (Apple silicon). **Hardware:** the R notebooks and figure scripts run on a standard desktop. scDRS, AUCell and the spin tests are run on a computing cluster (SLURM).

### R dependencies

R 4.4.0 with the following packages (versions tested):

| Category | Packages |
|----------|----------|
| **Tidyverse / data wrangling** | tidyverse 2.0.0, data.table 1.16.4, readxl 1.4.5, openxlsx 4.2.8.1, arrow 23.0.1.1, janitor 2.2.1, R.utils 2.13.0, scales 1.4.0 |
| **Visualization** | ggplot2 4.0.1, cowplot 1.2.0, patchwork 1.3.2, gridExtra 2.3, gridGraphics 0.5-1, ggplotify 0.1.3, ragg 1.3.3 |
| **ggplot2 extensions** | ggrepel 0.9.6, ggtext 0.1.2, ggh4x 0.3.1, ggforce 0.5.0, ggside 0.4.1, ggExtra 0.11.0, ggbeeswarm 0.7.2, ggpubr 0.6.2, ggdendro 0.2.0, ggnewscale 0.5.2, ggraph 2.2.2 |
| **Color / palettes** | RColorBrewer 1.1-3, viridis 0.6.5, paletteer 1.6.0, pals 1.9, colorspace 2.1-1, circlize 0.4.16, ggsci 4.1.0 |
| **Heatmaps** | ComplexHeatmap 2.20.0, pheatmap 1.0.13 |
| **Network / graph** | igraph 2.2.1, tidygraph 1.3.1, STRINGdb 2.20.0 |
| **Single-cell / genomics** | Seurat 5.3.1, SeuratDisk 0.0.0.9021, Signac 1.16.0, AUCell 1.26.0, anndata 0.8.0, reticulate 1.44.0, edgeR 4.2.2, slingshot 2.7.0 |
| **Genome annotation** | AnnotationHub 3.16.1, TxDb.Hsapiens.UCSC.hg19.knownGene 3.2.2, TxDb.Hsapiens.UCSC.hg38.knownGene 3.21.0, AnnotationDbi 1.66.0, GenomicFeatures 1.56.0, org.Hs.eg.db 3.19.1, GO.db 3.19.1, ensembldb 2.28.1, rtracklayer 1.64.0 |
| **Functional enrichment** | clusterProfiler 4.12.6, gprofiler2 0.2.4 |
| **Statistical modeling** | lme4 1.1-38, lmerTest 3.1-3, emmeans 2.0.0, mgcv 1.9-1, glmmSeq 0.5.7, metafor 4.6-0, mclust 6.1.2, mixtools 2.0.0.1 |
| **Text / topic modeling** | tidytext 0.4.3, stm 1.3.8 |
| **Other** | dendextend 1.19.1, ape 5.8-1, png 0.1-8, future 1.67.0, knitr 1.51, Matrix 1.7-1, eulerr 8.0.0, tidymodels 1.4.1, vip 0.4.1, ggfortify 0.4.19, tidyomics 1.4.0, writexl 1.5.4 |
| **Also loaded (versions not yet recorded)** | here, WriteXLS, scCustomize, scater, scran |

### Python dependencies

Python 3.9.14 (scDRS / trajectory modules):

| Package | Version |
|---------|---------|
| scdrs | 1.0.4 |
| scanpy | 1.10.3 |
| anndata | 0.10.9 |
| numpy | 1.26.4 |
| scipy | 1.13.1 |
| pandas | 2.3.3 |
| statsmodels | 0.14.5 |
| matplotlib | 3.9.4 |

Python 3.10 (spatial enrichment module, `analysis/04_brain_spatial_enrichment/`), via two pinned conda environments in that directory:

- [`ahba.yml`](analysis/04_brain_spatial_enrichment/ahba.yml) — AHBA matrix build (abagen 0.1.3; numpy 1.26, pandas 1.5, nilearn 0.10)
- [`spin.yml`](analysis/04_brain_spatial_enrichment/spin.yml) — spin tests and spatial enrichment (neuromaps 0.0.7; numpy 2.2, pandas 2.3, nilearn 0.14)

Create each with `conda env create -f <file>.yml`; see the header of each file for kernel registration.

### External tools

| Tool | Version | Purpose |
|------|---------|---------|
| MAGMA | v1.10 | GWAS gene-level analysis for scDRS common-variant input |
| MapMyCells | RRID: SCR_024672 | Allen Institute cell-type mapping (web service) |
| SynGO | v1.2 | Synaptic gene ontology database |

## Installation guide

Install R packages and Python dependencies listed above. R packages can be installed from CRAN/Bioconductor. Typical install time is ~30 minutes on a standard desktop.

## Data

[`data/`](data/) contains the published source inputs the analyses read: supplementary
tables from the network studies, reference atlases (AHBA, BrainSpan), SynGO, and curated
gene lists. See [`data/README.md`](data/README.md) for the full inventory and sources.

### External datasets not included in this repository

The following publicly available datasets are required for full reproduction but are not distributed here due to size or licensing. Obtain them from the original sources:

| Dataset | Source |
|---------|--------|
| Rare variant association results (gene-level counts, TADA) | [Satterstrom et al. 2026](https://www.medrxiv.org/content/10.64898/2026.08.24.26360398v1), Tables S7, S9, S17 |
| Developing snMultiome human brain atlas (Wang et al. 2025) | https://doi.org/10.5061/dryad.2280gb612 |
| uBrain mid-fetal cortical microarray (Ball et al. 2024) | https://zenodo.org/records/10622337 |
| Prenatal Visium spatial transcriptome (Aivazidis et al. 2025) | https://zenodo.org/records/14422018 |
| Pre-synaptic mass spectrometry (Dumrongprechachan et al. 2022) | https://www.ebi.ac.uk/pride/archive?keyword=PXD030864 |
| ASD GWAS summary statistics (Matoba et al. 2020) | https://bitbucket.org/steinlabunc/spark_asd_sumstats |
| BrainSpan developmental transcriptome (`expression_matrix.csv`) | https://www.brainspan.org/static/download.html |
| Allen Human Brain Atlas microarray (~4 GB; downloaded automatically by abagen in `04_00`) | https://human.brain-map.org/static/download |

## Demo

Each `analysis/` module contains numbered scripts that run on `data/` and the external datasets listed above. Expected outputs are provided in `outputs/analysis/` (statistical results) and `outputs/figures/` (rendered figures). To verify the setup, render any figure notebook (e.g. `quarto render figures/fig03/Fig3.qmd`), which reads from pre-computed outputs and completes in under one minute. Full re-execution of compute-intensive modules (scDRS, AUCell) requires the single-cell atlas on the HPC; see [`data/README.md`](data/README.md) for data locations. The HPC scripts for these modules (`analysis/01_scdrs_trajectory/`, `analysis/03_module_enrichments/03_01_AUCell.R`, `analysis/05_SynGO/05_04_aucell/`, `analysis/07_ddid_stratification/01_input_prep/`–`03_group_tests/`, and the `.R` plotting scripts in `figures/fig01/`, `fig03/`, `fig05/` and `supp_fig/`) use absolute cluster paths that must be edited to point at your local copies of the data.

## Repository layout

```
.
├── analysis/        # Executable analysis modules
│   ├── _shared/     # Shared R helpers used across modules
│   ├── 00_data_format/           # Format network + ASD rare-variant input data
│   ├── 01_scdrs_trajectory/      # scDRS disease-relevance scoring + lineage trajectory
│   ├── 02_module_identification/ # Stepwise conditional module identification
│   ├── 03_module_enrichments/    # GO/topic enrichments, module AUCell, BrainSpan development
│   ├── 04_brain_spatial_enrichment/ # AHBA cortical spin tests, developmental spatial enrichment
│   ├── 05_SynGO/                 # SynGO/HotNet enrichment, partition rate ratios, HotNet AUCell
│   ├── 06_MORPH_annotation/      # MORPH subclustering and TF-regulon/MEF2C module analyses
│   └── 07_ddid_stratification/   # DD/ID comorbidity-stratified analyses
├── figures/         # Manuscript figure reproduction by figure number
│   ├── fig01/ … fig05/
│   └── supp_fig/
├── functions/       # Shared plotting/network-analysis helper functions
├── outputs/
│   ├── analysis/    # Reusable statistical outputs, mirrored per analysis module
│   ├── figures/     # Rendered figure PDFs/PNGs
│   └── tables/      # Supplementary tables for submission
├── data/            # Published source inputs (supp. tables, atlases, gene lists)
├── assets/          # Quarto CSS and static assets
├── index.qmd        # Quarto site home page
└── _quarto.yml      # Quarto site configuration
```

## Analysis modules

Each module follows a stepwise numbered layout. Common subdirectory patterns:
`01_input_prep/`, `02_compute/`, `03_group_tests/`, `04_export/`.

File-name prefixes follow `<module>_<step>_<name>` (e.g.
`02_module_identification/02_00_*`), so the first number always matches the module
directory. `outputs/analysis/` subdirectories and their file prefixes mirror the same
numbering.

### [`analysis/_shared/`](analysis/_shared/)

Shared utilities used across modules:

- [`cluster_labels.R`](analysis/_shared/cluster_labels.R) — EN/IN lineage cluster label maps; used by AUCell group-test scripts, [`figures/fig01`](figures/fig01/), and [`figures/fig05`](figures/fig05/)
- [`syngo_helpers.R`](analysis/_shared/syngo_helpers.R) — SynGO loading, HotNet gene-set construction, AUCell ordering/color helpers; used by AUCell compute scripts and HotNet AUCell figure scripts
- [`scdrs_aucell_helpers.R`](analysis/_shared/scdrs_aucell_helpers.R) — scDRS × AUCell grouping, timepoint, control-column, and lineage-mapping helpers; used by AUCell scDRS correlation scripts

### [`analysis/00_data_format/`](analysis/00_data_format/)

Format raw inputs for downstream analysis.

- [`00_00_format-network-data.qmd`](analysis/00_data_format/00_00_format-network-data.qmd) — format gene co-expression network data
- [`00_01_format-ASD-rarevar-data.qmd`](analysis/00_data_format/00_01_format-ASD-rarevar-data.qmd) — format ASD rare-variant statistics

### [`analysis/01_scdrs_trajectory/`](analysis/01_scdrs_trajectory/)

scDRS disease-relevance scoring, lineage trajectory inference, and group-level lineage
association testing.

**[`01_scdrs/`](analysis/01_scdrs_trajectory/01_scdrs/)**
- [`01_input_prep/`](analysis/01_scdrs_trajectory/01_scdrs/01_input_prep/) — prepare single-cell data (`01_0`), scDRS covariates (`01_1`), downstream group-key metadata incl. `mclust_Group_region` and `subclass_donor_region` (`01_2`), RVAS gene set (`02`), and MAGMA gene Z-scores GRCh37/38 (`03_0`–`03_2`)
- [`02_compute/`](analysis/01_scdrs_trajectory/01_scdrs/02_compute/) — run scDRS for RVAS ([`01_run_scdrs_rvas.sh`](analysis/01_scdrs_trajectory/01_scdrs/02_compute/01_run_scdrs_rvas.sh)) and GWAS ([`02_run_scdrs_gwas.sh`](analysis/01_scdrs_trajectory/01_scdrs/02_compute/02_run_scdrs_gwas.sh))
- [`03_group_tests/`](analysis/01_scdrs_trajectory/01_scdrs/03_group_tests/) — lineage-level group association tests for RVAS and GWAS (`01`, `02`), mclust x Group x region tests pooled and paired-donor cohort (`03`, `04`), and per-donor subclass x region tests on the paired cohort (`05`)

**[`02_trajectory/`](analysis/01_scdrs_trajectory/02_trajectory/)**
- [`01_lineage_inference/`](analysis/01_scdrs_trajectory/02_trajectory/01_lineage_inference/) — generate EN and IN lineage objects and primary trajectory outputs
- [`02_postprocess/`](analysis/01_scdrs_trajectory/02_trajectory/02_postprocess/) — add metadata and prepare trajectory plotting inputs for [`figures/fig01`](figures/fig01/)
- [`03_secondary_analysis/`](analysis/01_scdrs_trajectory/02_trajectory/03_secondary_analysis/) — MapMyCell and MGE/SST subtype follow-up analyses

### [`analysis/02_module_identification/`](analysis/02_module_identification/)

Stepwise conditional analysis to identify ASD rare-variant modules.

- [`02_00_MAIN-linear-stepwise-forward-selection.qmd`](analysis/02_module_identification/02_00_MAIN-linear-stepwise-forward-selection.qmd) — main stepwise conditional analysis
- [`02_01_conditional-LR-sensitivity.qmd`](analysis/02_module_identification/02_01_conditional-LR-sensitivity.qmd) — binomial (ASD gene binary) sensitivity analysis

### [`analysis/03_module_enrichments/`](analysis/03_module_enrichments/)

GO/topic-model enrichments and developmental trajectories for identified modules.

- [`03_00_GO-enrichments-and-topic-modeling.qmd`](analysis/03_module_enrichments/03_00_GO-enrichments-and-topic-modeling.qmd) — module GO and cell-type enrichments, topic modeling
- [`03_01_AUCell.R`](analysis/03_module_enrichments/03_01_AUCell.R) — score gene sets via AUCell
- [`03_02_assign-genes-to-GO-topics.qmd`](analysis/03_module_enrichments/03_02_assign-genes-to-GO-topics.qmd) — assign genes to GO topics
- [`03_03_GR-subclustering.qmd`](analysis/03_module_enrichments/03_03_GR-subclustering.qmd) — GR subclustering
- [`03_04_BrainSpan-development.qmd`](analysis/03_module_enrichments/03_04_BrainSpan-development.qmd) — BrainSpan developmental trajectories

### [`analysis/04_brain_spatial_enrichment/`](analysis/04_brain_spatial_enrichment/)

Cortical spatial enrichment of ASD risk genes and gene programs against the Allen Human
Brain Atlas (AHBA), with spin-test significance and mubrain mid-fetal cortical region x tissue layer enrichments.

- [`04_00_build-ahba-matrix.ipynb`](analysis/04_brain_spatial_enrichment/04_00_build-ahba-matrix.ipynb) — build the AHBA region × gene expression matrix
- [`04_01_asd-spatial-enrichment.ipynb`](analysis/04_brain_spatial_enrichment/04_01_asd-spatial-enrichment.ipynb) — ASD risk-gene cortical enrichment + spin tests
- [`04_02_topic-spatial-enrichment.ipynb`](analysis/04_brain_spatial_enrichment/04_02_topic-spatial-enrichment.ipynb) — GO-topic/program spatial enrichment, Gandal 2022 dysregulation comparison
- [`04_03_tsyporin-SA-axis.qmd`](analysis/04_brain_spatial_enrichment/04_03_tsyporin-SA-axis.qmd) — sensorimotor–association axis overlap (Tsyporin patterning genes)
- [`04_04_mubrain-developmental-enrichment.qmd`](analysis/04_brain_spatial_enrichment/04_04_mubrain-developmental-enrichment.qmd) — developmental (mubrain) spatial enrichment
- [`04_05_module-spatial-maps.ipynb`](analysis/04_brain_spatial_enrichment/04_05_module-spatial-maps.ipynb) — per-module cortical maps and pairwise spatial correlations
- [`run_spin_test.py`](analysis/04_brain_spatial_enrichment/run_spin_test.py), [`maps_null_test.py`](analysis/04_brain_spatial_enrichment/maps_null_test.py), [`spatial_helpers.py`](analysis/04_brain_spatial_enrichment/spatial_helpers.py) — spin-test/null-model computation helpers
- [`spin_prep.sbatch`](analysis/04_brain_spatial_enrichment/spin_prep.sbatch), [`spin_array.sbatch`](analysis/04_brain_spatial_enrichment/spin_array.sbatch), [`submit_spin_tests.sh`](analysis/04_brain_spatial_enrichment/submit_spin_tests.sh), [`spin.yml`](analysis/04_brain_spatial_enrichment/spin.yml), [`ahba.yml`](analysis/04_brain_spatial_enrichment/ahba.yml) — HPC job scripts and environments
- [`dk_atlas_schematic.ipynb`](analysis/04_brain_spatial_enrichment/dk_atlas_schematic.ipynb) — Desikan-Killiany atlas schematic for figures

### [`analysis/05_SynGO/`](analysis/05_SynGO/)

SynGO/HotNet gene-set enrichment, partition rate-ratio tests, and AUCell correlation
with scDRS.

- [`05_00_gene-set-enrichments-vs-brain-background.qmd`](analysis/05_SynGO/05_00_gene-set-enrichments-vs-brain-background.qmd) — test gene sets for DNM enrichment vs. brain background
- [`05_01_SynGO-MAIN.R`](analysis/05_SynGO/05_01_SynGO-MAIN.R) — main SynGO enrichment analysis
- [`05_02_SynGO-HotNet-tree-prep.qmd`](analysis/05_SynGO/05_02_SynGO-HotNet-tree-prep.qmd) — SynGO HotNet trees (prep + BP plot)
- [`05_03_independent-partition-rate-ratio.qmd`](analysis/05_SynGO/05_03_independent-partition-rate-ratio.qmd) — rate-ratio calculations on independent SynGO partitions
- [`05_04_aucell/01_compute/`](analysis/05_SynGO/05_04_aucell/01_compute/) — score HotNet and topic gene sets via AUCell
- [`05_04_aucell/02_group_tests/`](analysis/05_SynGO/05_04_aucell/02_group_tests/) — correlate HotNet and topic AUCell scores with scDRS scores

### [`analysis/06_MORPH_annotation/`](analysis/06_MORPH_annotation/)

MORPH subclustering and TF-regulon/MEF2C-focused module analyses.

- [`06_00_module-HotNet-enrichments.qmd`](analysis/06_MORPH_annotation/06_00_module-HotNet-enrichments.qmd) — module HotNet and SynGO enrichment
- [`06_01_MORPH-subclustering.qmd`](analysis/06_MORPH_annotation/06_01_MORPH-subclustering.qmd) — MORPH subclustering
- [`06_02_ASD-TF-regulon-examples.qmd`](analysis/06_MORPH_annotation/06_02_ASD-TF-regulon-examples.qmd) — ASD TF regulon examples (incl. MEF2C)
- [`06_03_module-plots.qmd`](analysis/06_MORPH_annotation/06_03_module-plots.qmd) — module network plots

### [`analysis/07_ddid_stratification/`](analysis/07_ddid_stratification/)

DD/ID comorbidity-stratified enrichment and burden analyses.

- [`01_input_prep/`](analysis/07_ddid_stratification/01_input_prep/) — classify genes into DDID/NoDDID groups by per-proband rate and p_hat; save stratified gene lists for downstream scDRS
- [`02_compute/`](analysis/07_ddid_stratification/02_compute/) — run DD/ID sensitivity scDRS
- [`03_group_tests/`](analysis/07_ddid_stratification/03_group_tests/) — DD/ID sensitivity lineage group analysis
- [`04_gene_set_DDID_stratification/`](analysis/07_ddid_stratification/04_gene_set_DDID_stratification/) — module-level (`07_00`) and SynGO-level (`07_01`) O/E DD/ID enrichment analyses

Figure-specific plotting scripts are in [`figures/fig05/`](figures/fig05/).

### [`functions/`](functions/)

Shared plotting and network-analysis helpers used across `analysis/03_module_enrichments`
and `analysis/05_SynGO`:

- [`ORA.R`](functions/ORA.R) — overrepresentation analysis on two gene lists
- [`build_go_network.R`](functions/build_go_network.R) — gene–gene GO cosine-similarity network construction
- [`hypergeometric-overlap.R`](functions/hypergeometric-overlap.R) — hypergeometric gene-list overlap test
- [`plot_GO_similarity_graphs.R`](functions/plot_GO_similarity_graphs.R) — GO cosine-similarity "ball-and-stick" network plots
- [`syngo-tree-graph.R`](functions/syngo-tree-graph.R) — pruned ggraph tree from a SynGO CMU domain subset

[`setup.R`](setup.R) (repo root) is the library/theme/path setup script sourced at the top
of the R notebooks.

## Figures

Each `figures/figXX/` directory reproduces one manuscript figure from prepared files
in `outputs/` and documented external data sources. Panels are not recomputed from
raw data here.

| Directory | Contents |
|-----------|----------|
| [`fig01/`](figures/fig01/) | scDRS UMAP (`01`), lineage z-score (`02`), trajectory (`03`), and paired PFC/V1 per-donor subclass Z-score dumbbell plots (`04`) |
| [`fig02/`](figures/fig02/) | [`Fig2.qmd`](figures/fig02/Fig2.qmd) |
| [`fig03/`](figures/fig03/) | AUCell proportion heatmap (`01`), [`Fig3.qmd`](figures/fig03/Fig3.qmd) |
| [`fig04/`](figures/fig04/) | [`Fig4.qmd`](figures/fig04/Fig4.qmd) |
| [`fig05/`](figures/fig05/) | DD/ID gene-set distribution and lineage comparison dotplot |
| [`module_slides/`](figures/module_slides/) | [`module-slides.qmd`](figures/module_slides/module-slides.qmd) — per-module summary slides (stats, GO/SynGO enrichments, spatial maps, network plots) |
| [`supp_fig/`](figures/supp_fig/) | mclust type composition markers (`00`), scDRS RVAS vs GWAS group comparison (`01`), paired PFC/V1 mclust-group Z-diff (`02`), HotNet UMAP/pseudotime (`03`), topic × scDRS correlation (`04_1`), HotNet × scDRS correlation (`04_2`), regional enrichment burden analysis (`05`) |

## Outputs

- [`outputs/analysis/`](outputs/analysis/) — reusable statistical outputs mirrored per analysis module (module stats, GO enrichments, AUCell scores, SynGO/MORPH results, DD/ID enrichment tables)
- [`outputs/figures/`](outputs/figures/) — rendered figure PDFs/PNGs, including per-module summary slides ([`module_slides/`](outputs/figures/module_slides/)) and the final supplementary figures ([`supplement/`](outputs/figures/supplement/), `Supp-Fig-01.png`–`Supp-Fig-30.png`, numbered as in the supplement text)
- [`outputs/tables/`](outputs/tables/) — submission-ready supplementary tables (`TableS1`–`TableS18`, numbered as in the supplement text)
