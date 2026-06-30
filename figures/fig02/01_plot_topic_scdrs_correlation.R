library(data.table)
library(dplyr)
library(ggplot2)

source(file.path("analysis", "_shared", "scdrs_aucell_helpers.R"))

resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
aucdir <- paste0(resultdir, "WangNature/AUCell/")

celltype_groups <- get_scdrs_celltype_groups()
ct_order <- get_scdrs_ct_order(celltype_groups)

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
