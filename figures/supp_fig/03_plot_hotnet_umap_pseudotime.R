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
aucdir    <- paste0(resultdir, "WangNature/AUCell/")

# ============================================================
# Global ordering and color palettes
# ============================================================
gs_order <- get_hotnet_gs_order()
celltype_colors <- get_aucell_celltype_colors()
lineage_colors <- get_aucell_lineage_colors()

# ── Read saved data ───────────────────────────────────────────────────────────
meta <- fread(paste0(datapath, "metadata.csv")) %>% rename(cell_id = ID, cell_type = type, donor_id = Ident)
auc_df          <- fread(paste0(aucdir, "AUC_scores.csv")) %>% as.data.frame() %>% tibble::column_to_rownames("V1") 
auc_cell_long   <- fread(paste0(aucdir, "AUC_trajectory_pseudotime.csv"))
auc_lineage_long <- fread(paste0(aucdir, "auc_lineage_long.csv"))
auc_sig_tests   <- fread(paste0(aucdir, "AUC_lineage_wilcoxon.csv"))
umap_coords     <- fread(paste0(aucdir, "umap_coords.csv"))

# Re-apply factor levels after reading from CSV
auc_lineage_long$signature <- factor(auc_lineage_long$signature, levels = gs_order)
auc_lineage_long$lineage   <- factor(auc_lineage_long$lineage,   levels = c("RG", "EN", "IN", "Other"))
auc_cell_long$gene_set  <- factor(auc_cell_long$gene_set,  levels = gs_order)
auc_cell_long$cell_type <- factor(auc_cell_long$cell_type, levels = names(celltype_colors))
auc_cell_long$lineage   <- factor(auc_cell_long$lineage,   levels = c("EN lineage", "IN lineage"))
auc_sig_tests$signature <- factor(auc_sig_tests$signature, levels = gs_order)

# ── Plot Violin: AUC score by lineage with pseudobulk significance ─────────
auc_lineage_median <- auc_lineage_long %>%
  group_by(lineage, signature) %>%
  summarise(median_AUC = median(AUC, na.rm = TRUE), .groups = "drop") # cell level AUC value

auc_sig_tests2 <- auc_sig_tests %>%
  filter(p_label != "ns") %>%
  mutate(signature = as.character(signature)) %>%
  group_by(signature) %>%
  arrange(group1, group2, .by_group = TRUE) %>%
  mutate(comp_id = row_number()) %>%
  left_join(
    auc_lineage_long %>%
      mutate(signature = as.character(signature)) %>%
      group_by(signature) %>%
      summarise(y_base = max(AUC, na.rm = TRUE), .groups = "drop"),
    by = "signature"
  ) %>%
  mutate(y.position = y_base * (1 + 0.08 * comp_id)) %>%
  ungroup()

make_lineage_violin <- function(gs_subset, show_title = TRUE) {
  en_long_sub <- auc_lineage_long %>%
    filter(signature %in% gs_subset) %>%
    mutate(signature = factor(signature, levels = gs_subset))
  
  tests_sub <- auc_sig_tests2 %>%
    filter(signature %in% gs_subset) %>%
    mutate(signature = factor(signature, levels = gs_subset))
  
  median_sub <- auc_lineage_median %>%
    filter(signature %in% gs_subset) %>%
    mutate(signature = factor(signature, levels = gs_subset))
  
  p <- ggplot(en_long_sub, aes(x = lineage, y = AUC, fill = lineage)) +
    geom_violin(scale = "width", linewidth = 0.3) +
    geom_boxplot(width = 0.1, outlier.size = 0.1, fill = "white", linewidth = 0.3) +
    geom_text(data = median_sub,
              aes(x = lineage, y = median_AUC, label = sprintf("%.2f", median_AUC)),
              size = 1.7, vjust = -0.8, fontface = "bold", color = "black",
              inherit.aes = FALSE) +
    stat_pvalue_manual(data = tests_sub, label = "p_label",
                       xmin = "group1", xmax = "group2", y.position = "y.position",
                       size = 1.7, inherit.aes = FALSE) +
    facet_wrap(~ signature, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = lineage_colors) +
    theme_classic(base_size = 5) +
    theme(
      plot.title = element_text(size = 6),
      strip.background = element_blank(),
      strip.text = element_text(size = 5),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "AUC score")
  
  if (show_title) {
    p <- p + labs(title = "AUC score comparison across lineage (pseudobulk Wilcoxon)")
  }
  p
}

gs_pre  <- gs_order[grepl("^Presynaptic|^Active zone", gs_order)]
gs_post <- gs_order[grepl("^Postsynaptic",             gs_order)]

p_violin_pre_lineage  <- make_lineage_violin(gs_pre,  show_title = TRUE)
p_violin_post_lineage <- make_lineage_violin(gs_post, show_title = FALSE)

p_violin_lineage <- p_violin_pre_lineage / 
  (p_violin_post_lineage + patchwork::plot_spacer() +
     patchwork::plot_layout(widths = c(2, 1)))
ggsave(paste0(aucdir, "AUC_violin_lineage_sig.pdf"), plot= p_violin_lineage, width = 5, height = 5, dpi = 300)


# ── Plot Pseudotime trajectory ────────────────────────────────────────────
p_pseudotime <- ggplot(auc_cell_long, aes(x = average_pseudotime, y = AUC)) +
  geom_point(aes(color = cell_type), size = 0.1, alpha = 0.15) +
  geom_smooth(aes(group = 1), method = "gam", color = "black",
              linewidth = 0.5, se = TRUE, fill = "grey70", alpha = 0.2) +
  facet_grid(lineage ~ gene_set, scales = "free") +
  scale_color_manual(values = celltype_colors, name = "Cell type", na.translate = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(x = "Average pseudotime", y = "AUC score",
       title = "AUC score of SynGO HotNet gene set along pseudotime") +
  theme_classic(base_size = 5) +
  theme(
    plot.title = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    panel.grid.major.x = element_blank()
  )

# ggsave(paste0(aucdir, "AUC_trajectory_pseudotime.pdf"), plot = p_pseudotime, width = 12, height = 5, dpi = 300)
ggsave(paste0(aucdir, "AUC_trajectory_pseudotime.png"), plot = p_pseudotime, width = 12, height = 5, dpi = 300)


# ── Plot UMAP ─────────────────────────────────────────────────────────────
plot_df <- auc_df %>%
  left_join(umap_coords, by = "cell_id") %>%
  left_join(meta %>% select(cell_id, cell_type), by = "cell_id")

auc_min <- min(auc_df[, gs_order], na.rm = TRUE)
auc_max <- max(auc_df[, gs_order], na.rm = TRUE)

plot_list <- purrr::map(gs_order, function(name) {
  ggplot(plot_df, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = .data[[name]])) +
    geom_point(size = 0.1, alpha = 0.3) +
    scale_color_viridis_c(limits = c(auc_min, auc_max)) +
    labs(title = name, color = "AUC") +
    theme_void(base_size = 5) +
    theme(
      plot.title = element_text(size = 5, hjust = 0.5),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.3, "cm")
    )
})

p_umap = patchwork::wrap_plots(plot_list, ncol = 5, nrow = 1) +
  patchwork::plot_layout(guides = "collect")
ggsave(paste0(aucdir, "UMAP_all_genesets.png"), plot = p_umap, width = 10, height = 2, dpi = 600)
# ggsave(paste0(aucdir, "UMAP_all_genesets.pdf"), width = 9, height = 5, dpi = 300)

p_full <- (p_umap / p_pseudotime / p_violin_lineage) +
  patchwork::plot_layout(heights = c(1, 2.5, 2), guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )
ggsave(paste0(aucdir, "AUCell_combined.png"),plot = p_full, width = 7, height = 6.5, dpi = 600)
