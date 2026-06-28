library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

source(file.path("analysis", "_shared", "scdrs_aucell_helpers.R"))
source(file.path("analysis", "_shared", "syngo_helpers.R"))

resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
aucdir <- paste0(resultdir, "WangNature/AUCell/")

gs_order <- get_hotnet_gs_order()
celltype_groups <- get_scdrs_celltype_groups()
ct_order <- get_scdrs_ct_order(celltype_groups)
timepoint_colors <- get_scdrs_timepoint_colors()
timepoint_labels <- get_scdrs_timepoint_labels()
timepoint_order <- get_scdrs_timepoint_order()

cor_overall_ct <- fread(paste0(aucdir, "scDRS_AUCell_cor_overall_celltype.csv")) %>%
  as.data.frame() %>%
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

cor_results <- fread(paste0(aucdir, "scDRS_AUCell_cor_celltype.csv")) %>%
  as.data.frame() %>%
  mutate(
    lineage_label = factor(lineage_label, levels = ct_order),
    gene_set = factor(gene_set, levels = gs_order),
    Group = factor(Group, levels = timepoint_order)
  )

quintile_en_it_l4 <- fread(paste0(aucdir, "scDRS_AUCell_quintile_en_it_l4.csv")) %>%
  as.data.frame() %>%
  mutate(gene_set = factor(gene_set, levels = gs_order), Group = factor(Group, levels = timepoint_order))
quintile_en_it_l5 <- fread(paste0(aucdir, "scDRS_AUCell_quintile_en_it_l5.csv")) %>%
  as.data.frame() %>%
  mutate(gene_set = factor(gene_set, levels = gs_order), Group = factor(Group, levels = timepoint_order))
quintile_en_newborn <- fread(paste0(aucdir, "scDRS_AUCell_quintile_en_newborn.csv")) %>%
  as.data.frame() %>%
  mutate(gene_set = factor(gene_set, levels = gs_order), Group = factor(Group, levels = timepoint_order))
quintile_ipc_en <- fread(paste0(aucdir, "scDRS_AUCell_quintile_ipc_en.csv")) %>%
  as.data.frame() %>%
  mutate(gene_set = factor(gene_set, levels = gs_order), Group = factor(Group, levels = timepoint_order))

# lineage correlation plot              
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
  scale_fill_gradientn(colors = c("#2166ac", "#F7F7F7", "#b2182b"), name = "Spearman r") +
  theme_classic(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6)
  ) +
  labs(title = "Overall: cell type x gene set") +
  coord_flip()

ggsave(paste0(aucdir, "scDRS_AUCell_cor_overall_ct.pdf"),
       plot = p_overall_ct, width = 4, height = 2.5, unit = "in", dpi = 300)

# lineage*group correlation plot                        
sig_gs_pre <- c("Active zone", "Presynaptic membrane")
sig_gs_post <- c("Postsynaptic organization")

cor_pre <- cor_results %>%
  filter(lineage_label %in% c("EN-L4-IT", "EN-L4-IT-V1", "EN-L5-IT"),
         gene_set %in% sig_gs_pre) %>%
  mutate(
    sig_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    lineage_label = factor(lineage_label, levels = c("EN-L4-IT", "EN-L4-IT-V1", "EN-L5-IT")),
    gene_set = factor(gene_set, levels = sig_gs_pre)
  )

cor_post <- cor_results %>%
  filter(lineage_label %in% c("IPC-EN", "EN-Newborn"),
         gene_set == "Postsynaptic organization") %>%
  mutate(
    sig_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    lineage_label = factor(lineage_label, levels = c("IPC-EN", "EN-Newborn"))
  ) %>%
  filter(sig_label != "") %>%
  mutate(Group = factor(Group, levels = timepoint_order))

cor_pre <- cor_pre %>%
  mutate(Group = factor(Group, levels = timepoint_order))

fill_scale <- scale_fill_gradientn(
  colors = c("#F7F7F7", "#FCBBA1", "#CB181D", "#67000D"),
  na.value = "white",
  name = "Spearman r"
)

p_pre <- ggplot(cor_pre, aes(x = Group, y = lineage_label, fill = r_disease)) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  geom_text(aes(label = sig_label), size = 2.5, color = "white", vjust = 0.8) +
  fill_scale +
  scale_y_discrete(limits = rev(c("EN-L4-IT", "EN-L4-IT-V1", "EN-L5-IT"))) +
  scale_x_discrete(labels = timepoint_labels) +
  facet_wrap(~ gene_set, nrow = 1, scales = "free_x") +
  theme_classic(base_size = 6) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  )

p_post <- ggplot(cor_post, aes(x = Group, y = lineage_label, fill = r_disease)) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  geom_text(aes(label = sig_label), size = 2.5, color = "white", vjust = 0.8) +
  fill_scale +
  scale_y_discrete(limits = rev(c("IPC-EN", "EN-Newborn"))) +
  scale_x_discrete(labels = timepoint_labels) +
  facet_wrap(~ gene_set, nrow = 1, scales = "free_x") +
  theme_classic(base_size = 6) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_blank()
  )

p_sig_time <- p_pre / (p_post + patchwork::plot_spacer())
p_ab <- (p_overall_ct + p_sig_time) + patchwork::plot_layout(widths = c(1, 1))
p_abc <- p_ab / p_quintile_l4_sig + patchwork::plot_layout(heights = c(1.5, 1))
ggsave(paste0(aucdir, "scDRS_AUCell_combined_abc.png"),
       plot = p_abc, width = 6, height = 4.7, unit = "in", dpi = 300)

# quintile plot              
plot_quintile <- function(data, title, timepoint_colors, timepoint_labels,
                          show_x = TRUE, show_strip = TRUE, legend_position = "none") {
  ggplot(data, aes(x = auc_quintile, y = mean_scdrs, color = Group, group = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 0.2) +
    geom_point(size = 0.5) +
    geom_errorbar(aes(ymin = mean_scdrs - se_scdrs, ymax = mean_scdrs + se_scdrs),
                  width = 0.1, linewidth = 0.2) +
    facet_wrap(~ gene_set, nrow = 1, scales = "free_y") +
    scale_color_manual(values = timepoint_colors, labels = timepoint_labels) +
    scale_x_continuous(breaks = 0:4, labels = c("Q1", "Q2", "Q3", "Q4", "Q5")) +
    theme_classic(base_size = 6) +
    theme(
      strip.background = element_blank(),
      strip.text = if (show_strip) element_text(size = 6) else element_blank(),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 6),
      axis.title.x = if (show_x) element_text(size = 6) else element_blank(),
      axis.text.x = if (show_x) element_text(size = 6) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank()
    ) +
    labs(x = "AUC score quintile", y = "Average scDRS", title = title)
}

quintile_en_it_l4_sig <- quintile_en_it_l4 %>%
  filter(gene_set %in% sig_gs_pre) %>%
  mutate(gene_set = factor(gene_set, levels = sig_gs_pre))
quintile_en_it_l5_sig <- quintile_en_it_l5 %>%
  filter(gene_set %in% sig_gs_pre) %>%
  mutate(gene_set = factor(gene_set, levels = sig_gs_pre))
quintile_en_newborn_sig <- quintile_en_newborn %>%
  filter(gene_set %in% sig_gs_post) %>%
  mutate(gene_set = factor(gene_set, levels = sig_gs_post))
quintile_ipc_en_sig <- quintile_ipc_en %>%
  filter(gene_set %in% sig_gs_post) %>%
  mutate(gene_set = factor(gene_set, levels = sig_gs_post))

p_quintile_l4_sig <- plot_quintile(quintile_en_it_l4_sig, "EN-L4-IT", timepoint_colors, timepoint_labels, FALSE, TRUE, "none")
p_quintile_l5_sig <- plot_quintile(quintile_en_it_l5_sig, "EN-L5-IT", timepoint_colors, timepoint_labels, FALSE, TRUE, "none")
p_quintile_ipc_sig <- plot_quintile(quintile_ipc_en_sig, "IPC-EN", timepoint_colors, timepoint_labels, TRUE, TRUE, "none")
p_quintile_newborn_sig <- plot_quintile(quintile_en_newborn_sig, "EN-Newborn", timepoint_colors, timepoint_labels, TRUE, TRUE, "right")
      
p_quintile_combined <- p_overall_ct +
  (p_quintile_l4_sig / p_quintile_l5_sig / (p_quintile_ipc_sig + p_quintile_newborn_sig)) +
  patchwork::plot_layout(widths = c(1, 1.5))

ggsave(paste0(aucdir, "scDRS_AUCell_quintile_combined.pdf"),
       plot = p_quintile_combined, width = 6.5, height = 5, unit = "in", dpi = 300)
