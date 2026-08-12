# AUCell analysis on co-expression modules #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# NOTE: This script needs be run on an HPC cluster! 
# Run using the following code:
# sbatch --job-name=AUCell --mem=128G --cpus-per-task=4 --time=8:00:00 \
# --output=/mnt/isilon/gandal_lab/smithr30/ASD-rarevar-annot/outputs/analysis/03_module_enrichments/03_01_AUCell/AUCell_run_%j.log \
# --wrap="bash -lc 'module load R/4.5.1 && Rscript /mnt/isilon/gandal_lab/smithr30/ASD-rarevar-annot/code/analysis/03_module_enrichments/03_01_AUCell.R'"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


library(AUCell)
library(Seurat)
library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

# Set directories
project_dir   <- "/mnt/isilon/gandal_lab/smithr30/ASD-rarevar-annot/"
out_dir       <- file.path(project_dir, "outputs/")
wang_data_dir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"

# All AUCell outputs (final CSVs + intermediate .rds caches) live here
aucdir        <- file.path(out_dir, "analysis", "03_module_enrichments", "03_01_AUCell/")
dir.create(aucdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Module data loading and gene set construction
# ============================================================

# Gene network assignments
df_network_data <- readRDS(file.path(out_dir, "analysis", "00_data_format", "00_00_gene-networks-all.RDS"))

# Significant ASD modules
df_stepwise_res <- readRDS(
    file.path(out_dir, "analysis", "02_module_identification", "02_00_stepwise-linear-res-step1p0.001-step2p0.001.RDS")
)$lm_results %>% 
    mutate(module_label = factor(module_label, levels = unique(.$module_label)))

sig_module_genes <- df_network_data %>%
    unnest(cols = c(data)) %>%
    filter(module %in% df_stepwise_res$module) %>%
    ungroup() %>%
    dplyr::select(module, gene_name) %>%
    distinct() %>%
    group_by(module) %>%
    summarise(genes = list(gene_name), .groups = "drop") %>%
    tibble::deframe()


# ============================================================
# Module and cell-type ordering (for output matrix row/col order)
# ============================================================
gs_order <- df_stepwise_res$module

cell_order <- c(
  "RG-vRG", "RG-tRG", "RG-oRG",
  "IPC-EN", "EN-Newborn", "EN-IT-Immature",
  "EN-L2_3-IT", "EN-L4-IT", "EN-L5-IT", "EN-L6-IT",
  "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP", "EN-L6-CT", "EN-L6b",
  "IPC-Glia",
  "IN-dLGE-Immature", "IN-CGE-Immature", "IN-MGE-Immature",
  "IN-CGE-VIP", "IN-CGE-SNCG", "IN-Mix-LAMP5",
  "IN-MGE-SST", "IN-MGE-PV",
  "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous",
  "OPC", "Oligodendrocyte-Immature", "Oligodendrocyte",
  "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown"
)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Compute AUCell module enrichment per cell, threshold, and write the      ║
# ║  high-proportion matrix used by Fig2 (+ per-module thresholds).           ║
# ║  Requires: obj (SCT counts, ~74 GB) ONLY if cell_rankings.rds is absent   ║
# ║            meta (lightweight cell metadata CSV) for Steps 5–6             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# Cell metadata: expected columns: ID (cell barcode), type (cell type), Ident (donor ID)
meta <- fread(paste0(wang_data_dir, "metadata.csv")) %>%
  rename(cell_id = ID, cell_type = type, donor_id = Ident)

# ── Step 1. Build expression matrix and cell rankings ────────────────────────
# The ~74 GB Seurat object is loaded ONLY when rankings must be built; with a
# cached cell_rankings.rds it is never read, which drastically lowers the memory need.
rankings_cache <- paste0(aucdir, "cell_rankings.rds")
if (file.exists(rankings_cache)) {
  message("Step 1: Loading cached cell rankings...")
  cell_rankings <- readRDS(rankings_cache)
} else {
  message("Step 1: Building cell rankings (this is slow)...")
  obj <- readRDS(paste0(wang_data_dir, "snMultiome_atlas_Seurat_object.rds"))
  DefaultAssay(obj) <- "SCT"
  counts <- GetAssayData(obj, assay = "SCT", layer = "counts")
  cell_rankings <- AUCell_buildRankings(counts, nCores = 4, plotStats = FALSE)
  saveRDS(cell_rankings, rankings_cache)
  message("Step 1: Done. cell_rankings saved.")
}


# ── Step 2. Calculate AUC scores ─────────────────────────────────────────────
auc_cache <- paste0(aucdir, "cell_AUC.rds")
if (file.exists(auc_cache)) {
  message("Step 2: Loading cached AUC scores...")
  cell_AUC <- readRDS(auc_cache)
} else {
  message("Step 2: Calculating AUC scores...")
  cell_AUC <- AUCell_calcAUC(
    sig_module_genes,
    cell_rankings,
    aucMaxRank = ceiling(0.05 * nrow(cell_rankings))
  )
  saveRDS(cell_AUC, auc_cache)
  message("Step 2: Done. cell_AUC saved.")
}


# ── Step 3. Explore thresholds ───────────────────────────────────────────────
assignment_cache <- paste0(aucdir, "cells_assignment.rds")
if (file.exists(assignment_cache)) {
  message("Step 3: Loading cached threshold assignments...")
  cells_assignment <- readRDS(assignment_cache)
} else {
  message("Step 3: Exploring thresholds...")
  cells_assignment <- AUCell_exploreThresholds(cell_AUC, plotHist = FALSE, assign = TRUE)
  saveRDS(cells_assignment, assignment_cache)
  message("Step 3: Done. cells_assignment saved.")
}

print(getThresholdSelected(cells_assignment))


# ── Step 4. Extract Global_k1 threshold automatically for all modules ─────────
message("Step 4: Extracting Global_k1 thresholds...")
thresholds <- sapply(gs_order, function(gs_name) {
  thr <- cells_assignment[[gs_name]]$aucThr$thresholds
  if ("Global_k1" %in% rownames(thr)) thr["Global_k1", "threshold"] else NA
})
message("Step 4: Done.")


# ── Step 5. Assign cells to High/Low groups based on Global_k1 ───────────────
message("Step 5: Assigning cells to High/Low groups...")
auc_matrix <- getAUC(cell_AUC)

group_df <- as.data.frame(t(auc_matrix)) %>%
  tibble::rownames_to_column("cell_id")

for (gs_name in gs_order) {
  group_df[[gs_name]] <- ifelse(group_df[[gs_name]] >= thresholds[gs_name], "High", "Low")
}

group_df <- group_df %>%
  left_join(meta %>% select(cell_id, cell_type), by = "cell_id")
message("Step 5: Done.")


# ── Step 6. High-group proportion matrix (module × cell type, % of cells "High")
message("Step 6: Building high-group proportion matrix...")
cell_type_levels <- unique(meta$cell_type)
high_matrix <- t(sapply(gs_order, function(gs_name) {
  groups <- group_df[[gs_name]]
  sapply(cell_type_levels, function(ct) {
    ct_idx <- meta$cell_type == ct
    sum(groups[ct_idx] == "High", na.rm = TRUE) / sum(ct_idx)
  })
}))
cell_order_valid <- cell_order[cell_order %in% colnames(high_matrix)]
high_matrix      <- high_matrix[, cell_order_valid]
high_matrix_pct  <- high_matrix * 100
message("Step 6: Done.")


# ── Step 7. Save outputs ─────────────────────────────────────────────────────
message("Step 7: Saving outputs...")

# High-group proportion matrix (module × cell type, % of cells "High") — Fig2 / supplement table
write.csv(high_matrix_pct, paste0(aucdir, "high_proportion_pct.csv"))

# Global_k1 AUC threshold per module (reproducibility / supplement)
thresholds_df <- data.frame(
  gene_set  = gs_order,
  threshold = thresholds[gs_order],
  thr_name  = "Global_k1"
)
write.csv(thresholds_df, paste0(aucdir, "thresholds.csv"), row.names = FALSE)

message("Done. Outputs written to ", aucdir)


