library(AUCell)
library(pheatmap)
library(Seurat)
library(readxl)
library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(ggpubr)

source(file.path("analysis", "_shared", "syngo_helpers.R"))

datapath  <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
aucdir    <- paste0(resultdir, "WangNature/AUCell/topic/")

# ============================================================
# Gene set construction
# ============================================================
gene_module <- read_excel("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SynGO/gene_topic_HotNet_category.xlsx")
gene_module %>%
  count(GO_topic, name = "n_genes") %>%
  arrange(desc(n_genes))
# GO_topic n_genes
# 1 CHR          155
# 2 SYN           54
# 3 MORPH         44

all_genesets <- gene_module %>%
  group_by(GO_topic) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  tibble::deframe()

# ============================================================
# Global ordering and color palettes
# ============================================================
gs_order = names(all_genesets)
cell_order <- get_aucell_cell_order()
celltype_colors <- get_aucell_celltype_colors()
lineage_colors <- get_aucell_lineage_colors()

# Cell metadata: read directly from obj's metadata CSV
# Expected columns: ID (cell barcode), type (cell type), Ident (donor ID)
meta <- fread(paste0(datapath, "metadata.csv")) %>%
  rename(cell_id = ID, cell_type = type, donor_id = Ident)

obj <- readRDS(paste0(datapath, "snMultiome_atlas_Seurat_object.rds"))

# ── Step 1. Build expression matrix and cell rankings ────────────────────────
DefaultAssay(obj) <- "SCT"
counts <- GetAssayData(obj, assay = "SCT", layer = "counts")

pdf(paste0(aucdir, "gene_cell_number_histogram.pdf"), width = 7, height = 5)
cell_rankings <- AUCell_buildRankings(counts, nCores = 4, plotStats = TRUE)
dev.off()

# ── Step 2. Calculate AUC scores ─────────────────────────────────────────────
cell_AUC <- AUCell_calcAUC(
  all_genesets,
  cell_rankings,
  aucMaxRank = ceiling(0.05 * nrow(cell_rankings))
)

# ── Step 3. Explore thresholds ───────────────────────────────────────────────
cells_assignment <- AUCell_exploreThresholds(cell_AUC, plotHist = FALSE, assign = TRUE)

pdf(paste0(aucdir, "AUCell_topic_all_thresholds.pdf"), width = 25, height = 6)
par(mfrow = c(1, length(all_genesets)))
for (gs_name in names(all_genesets)) {
  AUCell_exploreThresholds(cell_AUC[gs_name, ], plotHist = TRUE, assign = TRUE)
}
dev.off()

for (gs_name in names(all_genesets)) {
  png(paste0(aucdir, "AUCell_", gsub(" ", "_", gs_name), "_threshold.png"),
      width = 6, height = 5, units = "in", res = 300)
  par(mfrow = c(1, 1))
  AUCell_exploreThresholds(cell_AUC[gs_name, ], plotHist = TRUE, assign = TRUE)
  dev.off()
}

print(getThresholdSelected(cells_assignment))

thresholds <- sapply(gs_order, function(gs_name) {
  thr <- cells_assignment[[gs_name]]$aucThr$thresholds
  if ("Global_k1" %in% rownames(thr)) thr["Global_k1", "threshold"] else NA
})

# ── Step 5. Assign cells to groups based on thresholds ───────────────────────
auc_matrix <- getAUC(cell_AUC)

group_df <- as.data.frame(t(auc_matrix)) %>%
  tibble::rownames_to_column("cell_id")

for (gs_name in gs_order) {
  group_df[[gs_name]] <- ifelse(group_df[[gs_name]] >= thresholds[gs_name], "High", "Low")
}

group_df <- group_df %>%
  left_join(meta %>% select(cell_id, cell_type), by = "cell_id")

# ── Step 6. Build confusion matrix (auto-threshold, cell counts & proportions)
cellLabels       <- data.frame(cell_type = meta$cell_type, row.names = meta$cell_id)
cell_type_levels <- unique(meta$cell_type)

confMatrix <- t(sapply(cells_assignment, function(x) {
  tab <- table(cellLabels[x$assignment, ])
  counts_per_type <- tab[cell_type_levels]
  names(counts_per_type) <- cell_type_levels
  counts_per_type
}))
confMatrix[is.na(confMatrix)] <- 0

cell_type_counts <- table(meta$cell_type)[cell_type_levels]
confMatrix_prop  <- sweep(confMatrix, 2, as.numeric(cell_type_counts), "/")
confMatrix_prop  <- confMatrix_prop[, cell_order]


# ── Step 7. High-group proportion matrix (threshold-based) ───────────────────
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


# ── Step 8. Stacked bar data (cell type × group proportion) ──────────────────
long_df <- purrr::map_dfr(gs_order, function(gs_name) {
  group_df %>%
    group_by(cell_type, group = .data[[gs_name]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cell_type) %>%
    mutate(prop = n / sum(n), gene_set = gs_name,
           group = factor(group, levels = c("Low", "High")))
}) %>%
  mutate(cell_type = factor(cell_type, levels = cell_order),
         gene_set  = factor(gene_set,  levels = gs_order))


# ── Step 9. Prepare cell-level AUC data frame ────────────────────────────────
auc_df         <- as.data.frame(t(auc_matrix))
auc_df$cell_id <- rownames(auc_df)

auc_cell <- auc_df %>%
  left_join(meta %>% select(cell_id, cell_type, donor_id), by = "cell_id") %>%
  mutate(
    lineage = case_when(
      grepl("^RG",           cell_type) ~ "RG",
      grepl("^EN|^IPC-EN",   cell_type) ~ "EN",
      grepl("^IN|^IPC-Glia", cell_type) ~ "IN",
      TRUE ~ "Other"
    )
  )


# ── Step 10. Pseudotime trajectory data ──────────────────────────────────────
load_metadata <- function(lineage_name) {
  data.table::fread(file.path(
    resultdir, "WangNature/slingshot", lineage_name,
    paste0(lineage_name, "_lineage_slingshot_results.csv")
  )) %>%
    as.data.frame() %>%
    select(ID, type, average_pseudotime) %>%
    rename(cell_id = ID, cell_type = type) %>%
    filter(!is.na(average_pseudotime))
}

auc_cell_long <- bind_rows(
  auc_df %>% inner_join(load_metadata("EN"), by = "cell_id") %>% mutate(lineage = "EN lineage"),
  auc_df %>% inner_join(load_metadata("IN"), by = "cell_id") %>% mutate(lineage = "IN lineage")
) %>%
  tidyr::pivot_longer(cols = all_of(gs_order), names_to = "gene_set", values_to = "AUC") %>%
  mutate(
    gene_set  = factor(gene_set,  levels = gs_order),
    cell_type = factor(cell_type, levels = names(celltype_colors)),
    lineage   = factor(lineage,   levels = c("EN lineage", "IN lineage"))
  ) %>%
  filter(!is.na(cell_type))


# ── Step 11. Lineage Wilcoxon tests ──────────────────────────────────────────
auc_lineage_long <- auc_cell %>%
  mutate(lineage = factor(lineage, levels = c("RG", "EN", "IN", "Other"))) %>%
  tidyr::pivot_longer(cols = all_of(gs_order), names_to = "signature", values_to = "AUC") %>%
  mutate(signature = factor(signature, levels = gs_order)) %>%
  filter(!is.na(lineage))

auc_pb <- auc_cell %>%
  mutate(lineage = factor(lineage, levels = c("RG", "EN", "IN", "Other"))) %>%
  group_by(donor_id, lineage) %>%
  summarise(across(all_of(gs_order), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

comparisons <- combn(c("RG", "EN", "IN", "Other"), 2, simplify = FALSE)

auc_sig_tests <- purrr::map_dfr(gs_order, function(sig_name) {
  purrr::map_dfr(comparisons, function(comp) {
    g1   <- auc_pb %>% filter(lineage == comp[1]) %>% pull(sig_name)
    g2   <- auc_pb %>% filter(lineage == comp[2]) %>% pull(sig_name)
    test <- wilcox.test(g1, g2)
    data.frame(signature = sig_name, group1 = comp[1], group2 = comp[2],
               pvalue = test$p.value,
               median_g1 = median(g1, na.rm = TRUE), median_g2 = median(g2, na.rm = TRUE))
  })
}) %>%
  mutate(
    FDR     = p.adjust(pvalue, method = "BH"),
    p_label = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**",
                        FDR < 0.05  ~ "*",   TRUE        ~ "ns"),
    signature = factor(signature, levels = gs_order)
  )

# ── Step 12. Save all outputs ─────────────────────────────────────────────────

# AUCell objects
saveRDS(cell_AUC,         paste0(aucdir, "cell_AUC.rds"))
saveRDS(cells_assignment, paste0(aucdir, "cells_assignment.rds"))
# saveRDS(cell_rankings,  paste0(aucdir, "cell_rankings.rds"))  # large; optional

# AUC score matrix (cell × gene set)
write.csv(auc_df,
          paste0(aucdir, "AUC_scores.csv"))

# Threshold summary (Global_k1 for all topic)
thresholds_df <- data.frame(
  gene_set  = gs_order,
  threshold = thresholds[gs_order],
  thr_name  = "Global_k1"
)
write.csv(thresholds_df,
          paste0(aucdir, "thresholds.csv"), row.names = FALSE)

# Per-cell group assignments
write.csv(group_df,
          paste0(aucdir, "cell_group_assignment.csv"))

# Stacked bar long format
write.csv(long_df,
          paste0(aucdir, "celltype_proportion_stacked.csv"), row.names = FALSE)

# High-group proportion matrix
write.csv(high_matrix_pct,
          paste0(aucdir, "high_proportion_pct.csv"))

# Confusion matrix (auto-threshold)
write.csv(confMatrix,
          paste0(aucdir, "confMatrix_counts.csv"))
write.csv(confMatrix_prop,
          paste0(aucdir, "confMatrix_proportion.csv"))

# Cell-level data with lineage annotation
write.csv(auc_cell,
          paste0(aucdir, "auc_cell.csv"), row.names = FALSE)

# Pseudotime trajectory data
write.csv(auc_cell_long,
          paste0(aucdir, "AUC_trajectory_pseudotime.csv"), row.names = FALSE)

# Lineage Wilcoxon test results
write.csv(auc_sig_tests,
          paste0(aucdir, "AUC_lineage_wilcoxon.csv"), row.names = FALSE)

# Lineage long format (for violin)
write.csv(auc_lineage_long,
          paste0(aucdir, "auc_lineage_long.csv"), row.names = FALSE)
