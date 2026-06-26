# ASD Rare Variant Annotation Project

Research analysis project for ASD rare variant annotation and downstream single-cell
interpretation. 

## Repository layout

```
.
├── analysis/        # Executable analysis modules
│   ├── _shared/     # Shared R helpers used across modules
│   ├── 01_scdrs/    # scDRS disease-relevance scoring
│   ├── 02_trajectory/  # Lineage trajectory inference and pseudotime
│   ├── 03_aucell/   # AUCell single-cell gene-set enrichment
│   ├── 04_ddid_stratification/  # DD/ID comorbidity-stratified analyses
│   └── legacy/      # Superseded notebooks (reference only)
├── figures/         # Manuscript figure reproduction by figure number
│   ├── fig01/ … fig05/
│   └── supp_fig/
├── output/          
│   └── manuscript/  # Submission-ready deliverables
├── data/            # Path notes 
├── assets/          # Quarto CSS and static assets
├── index.qmd        # Quarto site home page
└── _quarto.yml      # Quarto site configuration
```

## Analysis modules

Each module follows a stepwise numbered layout. Common subdirectory patterns:
`01_input_prep/`, `02_compute/`, `03_group_tests/`, `04_export/`.

### `analysis/_shared/`

Shared utilities used across modules:

- `cluster_labels.R` — EN/IN lineage cluster label maps; used by AUCell group-test scripts and `figures/fig05`
- `syngo_helpers.R` — SynGO loading, HotNet gene-set construction, AUCell ordering/color helpers; used by AUCell compute scripts and HotNet AUCell figure scripts
- `scdrs_aucell_helpers.R` — scDRS × AUCell grouping, timepoint, control-column, and lineage-mapping helpers; used by AUCell scDRS correlation scripts

### `analysis/01_scdrs/`

scDRS disease-relevance scoring and group-level lineage association testing.

- `01_input_prep/` — prepare single-cell data (`01_0`), scDRS covariates (`01_1`), RVAS gene set (`02`), and MAGMA gene Z-scores GRCh37/38 (`03_0`–`03_2`)
- `02_compute/` — run scDRS for RVAS (`01_run_scdrs_rvas.sh`) and GWAS (`02_run_scdrs_gwas.sh`)
- `03_group_tests/` — lineage-level group association tests for RVAS and GWAS

### `analysis/02_trajectory/`

Trajectory inference, pseudotime binning, and lineage-level analyses.

- `01_lineage_inference/` — generate EN and IN lineage objects and primary trajectory outputs
- `02_postprocess/` — add metadata and prepare trajectory plotting inputs for `figures/fig01`
- `03_secondary_analysis/` — MapMyCell and MGE/SST subtype follow-up analyses


### `analysis/03_aucell/`

Single-cell AUCell scoring and downstream enrichment analyses.

- `01_compute/` — score HotNet and topic gene sets via AUCell
- `02_group_tests/` — correlate HotNet and topic AUCell scores with scDRS scores

### `analysis/04_ddid_stratification/`

DD/ID comorbidity-stratified enrichment and burden analyses.

- `01_input_prep/` — classify genes into DDID/NoDDID groups by per-proband rate and p_hat; save stratified gene lists for downstream scDRS
- `02_compute/` — run DD/ID sensitivity scDRS
- `03_group_tests/` — DD/ID sensitivity lineage group analysis

Figure-specific plotting scripts are in `figures/fig05/`.


## Figures

Each `figures/figXX/` directory reproduces one manuscript figure from prepared files
in `output/` and documented external data sources. Panels are not recomputed from
raw data here.

| Directory | Contents |
|-----------|----------|
| `fig01/` | scDRS UMAP, lineage z-score, and trajectory plots |
| `fig02/` | *(scripts in progress)* |
| `fig03/` | AUCell proportion heatmap |
| `fig04/` | Topic × scDRS correlation |
| `fig05/` | DD/ID gene-set distribution and lineage comparison dotplot |
| `supp_fig/` | HotNet–scDRS correlation (`01`), HotNet supplementary panels (`02`), regional enrichment burden analysis (`03`), scDRS RVAS vs GWAS group comparison (`04`) |

## Output

`output/manuscript/` holds submission-ready deliverables: final figure PDFs/PNGs,
exported tables, and submission-ready bundles.

## Data

No data files are stored in this repository. See [`data/README.md`](data/README.md).
