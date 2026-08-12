# ASD Rare Variant Annotation Project

Research analysis project for ASD rare variant annotation and downstream single-cell
interpretation.

## Repository layout

```
.
├── analysis/        # Executable analysis modules
│   ├── _shared/     # Shared R helpers used across modules
│   ├── 00_data_format/           # Format network + ASD rare-variant input data
│   ├── 01_scdrs_trajectory/      # scDRS disease-relevance scoring + lineage trajectory
│   ├── 02_module_identification/ # Stepwise conditional module identification
│   ├── 03_module_enrichments/    # GO/topic enrichments, AUCell, BrainSpan development
│   ├── 04_brain_spatial_enrichment/ # AHBA cortical spin tests, developmental spatial enrichment
│   ├── 05_SynGO/                 # SynGO/HotNet enrichment, partition rate ratios, AUCell
│   ├── 06_MORPH_annotation/      # MORPH subclustering and TF-regulon/MEF2C module analyses
│   └── 07_ddid_stratification/   # DD/ID comorbidity-stratified analyses
├── figures/         # Manuscript figure reproduction by figure number
│   ├── fig01/ … fig05/
│   └── supp_fig/
├── functions/       # Shared plotting/network-analysis helper functions
├── outputs/
│   ├── analysis/    # Reusable statistical outputs, mirrored per analysis module
│   ├── figures/     # Rendered figure PDFs/PNGs/SVGs
│   └── tables/      # Submission-ready supplementary tables
├── data/            # Path notes (no data files stored in this repo)
├── cmu/             # Collaborator (CMU) working files for DNM/SynGO enrichments — reference only, not part of the reproducible pipeline
├── assets/          # Quarto CSS and static assets
├── index.qmd        # Quarto site home page
└── _quarto.yml      # Quarto site configuration
```

## Analysis modules

Each module follows a stepwise numbered layout. Common subdirectory patterns:
`01_input_prep/`, `02_compute/`, `03_group_tests/`, `04_export/`.

File names within modules `02`–`07` keep the numbering of the development repo
(`rlsmith1/ASD-rarevar-annot`, where `01_scdrs_trajectory` does not exist and modules
run `01`–`06`), so file prefixes are offset by one from the directory numbers here.
Files are copied verbatim from that repo; `outputs/analysis/` subdirectories mirror
the `analysis/` module numbering.

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

- [`01_00_MAIN-linear-stepwise-forward-selection.qmd`](analysis/02_module_identification/01_00_MAIN-linear-stepwise-forward-selection.qmd) — main stepwise conditional analysis
- [`01_01_conditional-LR-sensitivity.qmd`](analysis/02_module_identification/01_01_conditional-LR-sensitivity.qmd) — binomial (ASD gene binary) sensitivity analysis

### [`analysis/03_module_enrichments/`](analysis/03_module_enrichments/)

GO/topic-model enrichments and developmental trajectories for identified modules.

- [`02_00_GO-enrichments-and-topic-modeling.qmd`](analysis/03_module_enrichments/02_00_GO-enrichments-and-topic-modeling.qmd) — module GO and cell-type enrichments, topic modeling
- [`02_01_AUCell.R`](analysis/03_module_enrichments/02_01_AUCell.R) — score gene sets via AUCell
- [`02_02_assign-genes-to-GO-topics.qmd`](analysis/03_module_enrichments/02_02_assign-genes-to-GO-topics.qmd) — assign genes to GO topics
- [`02_03_GR-subclustering.qmd`](analysis/03_module_enrichments/02_03_GR-subclustering.qmd) — GR subclustering
- [`02_04_BrainSpan-development.qmd`](analysis/03_module_enrichments/02_04_BrainSpan-development.qmd) — BrainSpan developmental trajectories
- [`MIKE_AHBA_ASD_analysis.ipynb`](analysis/03_module_enrichments/MIKE_AHBA_ASD_analysis.ipynb) — Allen Human Brain Atlas follow-up analysis

### [`analysis/04_brain_spatial_enrichment/`](analysis/04_brain_spatial_enrichment/)

Cortical spatial enrichment of ASD risk genes and gene programs against the Allen Human
Brain Atlas (AHBA), with spin-test significance (Alexander-Bloch permutations) and
developmental follow-ups.

- [`03_00_build-ahba-matrix.ipynb`](analysis/04_brain_spatial_enrichment/03_00_build-ahba-matrix.ipynb) — build the AHBA region × gene expression matrix
- [`03_01_asd-spatial-enrichment.ipynb`](analysis/04_brain_spatial_enrichment/03_01_asd-spatial-enrichment.ipynb) — ASD risk-gene cortical enrichment + spin tests
- [`03_02_topic-spatial-enrichment.ipynb`](analysis/04_brain_spatial_enrichment/03_02_topic-spatial-enrichment.ipynb) — GO-topic/program spatial enrichment, Gandal 2022 dysregulation comparison
- [`03_03_tsyporin-SA-axis.qmd`](analysis/04_brain_spatial_enrichment/03_03_tsyporin-SA-axis.qmd) — sensorimotor–association axis overlap (Tsyporin patterning genes)
- [`03_04_mubrain-developmental-enrichment.qmd`](analysis/04_brain_spatial_enrichment/03_04_mubrain-developmental-enrichment.qmd) — developmental (mubrain) spatial enrichment
- [`03_05_module-spatial-maps.ipynb`](analysis/04_brain_spatial_enrichment/03_05_module-spatial-maps.ipynb) — per-module cortical maps and pairwise spatial correlations
- [`run_spin_test.py`](analysis/04_brain_spatial_enrichment/run_spin_test.py), [`maps_null_test.py`](analysis/04_brain_spatial_enrichment/maps_null_test.py), [`spatial_helpers.py`](analysis/04_brain_spatial_enrichment/spatial_helpers.py) — spin-test/null-model computation helpers
- [`spin_prep.sbatch`](analysis/04_brain_spatial_enrichment/spin_prep.sbatch), [`spin_array.sbatch`](analysis/04_brain_spatial_enrichment/spin_array.sbatch), [`submit_spin_tests.sh`](analysis/04_brain_spatial_enrichment/submit_spin_tests.sh), [`spin.yml`](analysis/04_brain_spatial_enrichment/spin.yml), [`ahba.yml`](analysis/04_brain_spatial_enrichment/ahba.yml) — HPC job scripts and environments
- [`dk_atlas_schematic.ipynb`](analysis/04_brain_spatial_enrichment/dk_atlas_schematic.ipynb) — Desikan-Killiany atlas schematic for figures

### [`analysis/05_SynGO/`](analysis/05_SynGO/)

SynGO/HotNet gene-set enrichment, partition rate-ratio tests, and AUCell correlation
with scDRS.

- [`04_00_gene-set-enrichments-vs-brain-background.qmd`](analysis/05_SynGO/04_00_gene-set-enrichments-vs-brain-background.qmd) — test gene sets for DNM enrichment vs. brain background
- [`04_01_SynGO-MAIN.R`](analysis/05_SynGO/04_01_SynGO-MAIN.R) — main SynGO enrichment analysis
- [`04_02_SynGO-HotNet-tree-prep.qmd`](analysis/05_SynGO/04_02_SynGO-HotNet-tree-prep.qmd) — SynGO HotNet trees (prep + BP plot)
- [`04_03_independent-partition-rate-ratio.qmd`](analysis/05_SynGO/04_03_independent-partition-rate-ratio.qmd) — rate-ratio calculations on independent SynGO partitions
- [`03_03_aucell/01_compute/`](analysis/05_SynGO/03_03_aucell/01_compute/) — score HotNet and topic gene sets via AUCell
- [`03_03_aucell/02_group_tests/`](analysis/05_SynGO/03_03_aucell/02_group_tests/) — correlate HotNet and topic AUCell scores with scDRS scores

### [`analysis/06_MORPH_annotation/`](analysis/06_MORPH_annotation/)

MORPH subclustering and TF-regulon/MEF2C-focused module analyses.

- [`05_00_module-HotNet-enrichments.qmd`](analysis/06_MORPH_annotation/05_00_module-HotNet-enrichments.qmd) — module HotNet and SynGO enrichment
- [`05_01_MORPH-subclustering.qmd`](analysis/06_MORPH_annotation/05_01_MORPH-subclustering.qmd) — MORPH subclustering
- [`05_02_ASD-TF-regulon-examples.qmd`](analysis/06_MORPH_annotation/05_02_ASD-TF-regulon-examples.qmd) — ASD TF regulon examples (incl. MEF2C)
- [`05_03_module-plots.qmd`](analysis/06_MORPH_annotation/05_03_module-plots.qmd) — module network plots

### [`analysis/07_ddid_stratification/`](analysis/07_ddid_stratification/)

DD/ID comorbidity-stratified enrichment and burden analyses.

- [`01_input_prep/`](analysis/07_ddid_stratification/01_input_prep/) — classify genes into DDID/NoDDID groups by per-proband rate and p_hat; save stratified gene lists for downstream scDRS
- [`02_compute/`](analysis/07_ddid_stratification/02_compute/) — run DD/ID sensitivity scDRS
- [`03_group_tests/`](analysis/07_ddid_stratification/03_group_tests/) — DD/ID sensitivity lineage group analysis
- [`04_gene_set_DDID_stratification/`](analysis/07_ddid_stratification/04_gene_set_DDID_stratification/) — module-level (`06_00`) and SynGO-level (`06_01`) O/E DD/ID enrichment analyses

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
of the R notebooks (as `code/setup.R` in the development repo).

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
- [`outputs/figures/`](outputs/figures/) — rendered figure PDFs/PNGs/SVGs, including per-module summary slides ([`module_slides/`](outputs/figures/module_slides/))
- [`outputs/tables/`](outputs/tables/) — submission-ready supplementary tables

## Data

No data files are stored in this repository. See [`data/README.md`](data/README.md).
