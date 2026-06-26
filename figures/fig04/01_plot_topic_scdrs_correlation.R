library(data.table)
library(dplyr)
library(ggplot2)

source(file.path("analysis", "_shared", "scdrs_aucell_helpers.R"))

resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
aucdir <- paste0(resultdir, "WangNature/AUCell/")

celltype_groups <- get_scdrs_celltype_groups()
ct_order <- get_scdrs_ct_order(celltype_groups)
timepoint_labels <- get_scdrs_timepoint_labels()
timepoint_order <- get_scdrs_timepoint_order()

cor_overall_ct <- fread(paste0(aucdir, "scDRS_topic_AUCell_cor_overall_celltype.csv")) %>%
  as.data.frame()
gs_order <- unique(as.character(cor_overall_ct$gene_set))

cor_overall_ct <- cor_overall_ct %>%
  mutate(
    lineage_label = factor(lineage_label, levels = ct_order),
    gene_set = factor(gene_set, levels = gs_order),
    sig_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

cor_results <- fread(paste0(aucdir, "scDRS_topic_AUCell_cor_celltype.csv")) %>%
  as.data.frame() %>%
  mutate(
    lineage_label = factor(lineage_label, levels = ct_order),
    gene_set = factor(gene_set, levels = gs_order),
    Group = factor(Group, levels = timepoint_order)
  )

group_last_ct <- sapply(celltype_groups, function(types) types[length(types)])
separator_xpos <- sapply(group_last_ct[-length(group_last_ct)], function(ct) {
  which(ct_order == ct) + 0.5
})

p_overall_ct <- ggplot() +
  geom_tile(data = expand.grid(lineage_label = ct_order, gene_set = gs_order, stringsAsFactors = FALSE),
            aes(x = lineage_label, y = gene_set),
            fill = "white", color = "grey85", linewidth = 0.2) +
  geom_tile(data = cor_overall_ct %>% filter(n_cells > 150),
            aes(x = lineage_label, y = gene_set, fill = r_disease),
            color = "grey85", linewidth = 0.2) +
  geom_text(data = cor_overall_ct %>% filter(sig_label != ""),
            aes(x = lineage_label, y = gene_set, label = sig_label),
            size = 3, color = "white", vjust = 0.8) +
  geom_vline(xintercept = separator_xpos, color = "grey40", linewidth = 0.4) +
  scale_x_discrete(limits = ct_order) +
  scale_y_discrete(limits = rev(gs_order)) +
  scale_fill_gradientn(colors = c("white", "#fcbba1", "#b2182b", "#67000d"), name = "Spearman r") +
  theme_classic(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6)
  ) +
  labs(title = "Overall: cell type x gene set")

ggsave(paste0(aucdir, "scDRS_topic_AUCell_cor_overall_ct.png"),
       plot = p_overall_ct, width = 7, height = 4, dpi = 300)
ggsave(paste0(aucdir, "scDRS_topic_AUCell_cor_overall_ct.pdf"),
       plot = p_overall_ct, width = 5, height = 2, unit = "in", dpi = 300)

all_tiles_ct <- expand.grid(
  lineage_label = ct_order,
  Group = timepoint_order,
  gene_set = gs_order,
  stringsAsFactors = FALSE
) %>%
  mutate(
    lineage_label = factor(lineage_label, levels = ct_order),
    Group = factor(Group, levels = timepoint_order),
    gene_set = factor(gene_set, levels = gs_order)
  )

group_sizes <- sapply(celltype_groups, length)
separator_ypos_ct <- cumsum(rev(group_sizes)[-length(group_sizes)]) + 0.5

p_heatmap_ct <- ggplot() +
  geom_tile(data = all_tiles_ct,
            aes(x = Group, y = lineage_label),
            fill = "white", color = "grey85", linewidth = 0.2) +
  geom_tile(data = cor_results %>% filter(!is.na(r_plot) & n_cells > 150),
            aes(x = Group, y = lineage_label, fill = r_plot),
            color = "grey85", linewidth = 0.2) +
  geom_hline(yintercept = separator_ypos_ct, color = "grey40", linewidth = 0.4) +
  scale_fill_gradientn(
    colors = c("white", "#fcbba1", "#b2182b", "#67000d"),
    na.value = "white",
    name = "Spearman r\n(MC FDR < 0.05)"
  ) +
  scale_y_discrete(limits = rev(ct_order)) +
  scale_x_discrete(labels = timepoint_labels) +
  facet_wrap(~ gene_set, nrow = 1) +
  theme_classic(base_size = 6) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6)
  ) +
  labs(title = "scDRS x AUCell correlation: cell type level\n(MC-corrected, FDR < 0.05; white = not significant)")

ggsave(paste0(aucdir, "scDRS_topic_AUCell_cor_celltype_heatmap.png"),
       plot = p_heatmap_ct, width = 8, height = 5, dpi = 300)
ggsave(paste0(aucdir, "scDRS_topic_AUCell_cor_celltype_heatmap.pdf"),
       plot = p_heatmap_ct, width = 4.5, height = 3, unit = "in", dpi = 300)
