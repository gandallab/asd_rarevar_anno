library(cowplot)
library(ggdendro)
library(ape)
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
cell_order <- get_aucell_cell_order()
gs_levels <- get_hotnet_gs_levels()

# Define collapsed cell type groups
neuron_groups <- c(
  list("RG" = cell_order[grepl("^RG", cell_order)]),
  setNames(as.list(cell_order[grepl("^EN|^IPC-EN",  cell_order)]),
           cell_order[grepl("^EN|^IPC-EN",  cell_order)]),
  setNames(as.list(cell_order[grepl("^IN|^IPC-Glia", cell_order)]),
           cell_order[grepl("^IN|^IPC-Glia", cell_order)]),
  list("Other" = cell_order[grepl(
    "^Astrocyte|^OPC|^Oligodendrocyte|^Cajal|^Microglia|^Vascular|^Unknown",
    cell_order)])
)

group_map <- purrr::imap_dfr(neuron_groups, function(types, grp) {
  data.frame(cell_type = types, neuron_group = grp)
})

group_df_raw <- read.csv(paste0(aucdir, "cell_group_assignment.csv"),
                         row.names = 1, check.names = FALSE) %>%
  left_join(group_map, by = "cell_type")

col_order <- c(
  "RG",
  cell_order[grepl("^EN|^IPC-EN",   cell_order)],
  cell_order[grepl("^IN|^IPC-Glia", cell_order)],
  "Other"
)

high_grouped <- purrr::map_dfr(names(gs_levels), function(gs_name) {
  group_df_raw %>%
    filter(!is.na(neuron_group)) %>%
    group_by(neuron_group) %>%
    summarise(
      pct_high = sum(.data[[gs_name]] == "High", na.rm = TRUE) / n() * 100,
      .groups  = "drop"
    ) %>%
    mutate(gene_set = gs_name)
}) %>%
  tidyr::pivot_wider(names_from = neuron_group, values_from = pct_high) %>%
  tibble::column_to_rownames("gene_set") %>%
  .[gs_order, col_order]

display_labels_grouped <- matrix(
  ifelse(round(high_grouped) == 0, "", sprintf("%.0f", as.matrix(high_grouped))),
  nrow = nrow(high_grouped), dimnames = dimnames(high_grouped)
)

# plot
tree_str <- "(RG,(((IPC-EN,EN-Newborn),((EN-IT-Immature,(EN-L2_3-IT,EN-L4-IT,EN-L5-IT,EN-L6-IT)),(EN-Non-IT-Immature,(EN-L5-ET,EN-L5_6-NP,EN-L6-CT,EN-L6b)))),(IPC-Glia,(IN-dLGE-Immature,(IN-CGE-Immature,(IN-CGE-VIP,IN-CGE-SNCG,IN-Mix-LAMP5)),(IN-MGE-Immature,(IN-MGE-SST,IN-MGE-PV))))),Other);"
tree <- ape::read.tree(text = tree_str)

p_dend <- ggplotify::as.ggplot(function() {
  par(mar = c(0, 0, 0, 0))
  ape::plot.phylo(
    tree,
    direction   = "downwards",
    show.tip.label = FALSE,
    edge.width  = 0.5
  )
})

leaf_order <- tree$tip.label 
heatmap_long <- as.data.frame(as.table(as.matrix(high_grouped[gs_order, leaf_order]))) %>%
  rename(gene_set = Var1, cell_type = Var2, pct = Freq) %>%
  mutate(
    gene_set  = factor(gene_set,  levels = rev(gs_order)),
    cell_type = factor(cell_type, levels = leaf_order)
  )

p_heatmap <- ggplot(heatmap_long, aes(x = cell_type, y = gene_set, fill = pct)) +
  geom_tile(color = "grey80", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = c("#FFFFFF","#C2B5D1","#866CA4","#4A2377"), #colorRampPalette(c("white", "#4A2377"))(4)
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    name   = "Proportion"
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y      = element_text(size = 6),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    axis.line        = element_blank(),
    legend.key.size  = unit(0.3, "cm"),
    legend.position = "right", 
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6),
    plot.margin      = margin(0, 0, 0, 0)
  )

p_combined_heatmap <- cowplot::plot_grid(
  p_dend,
  p_heatmap,
  ncol        = 1,
  rel_heights = c(0.5, 1.5),
  align       = "v",
  axis        = "lr"
)
# p_combined_heatmap
ggsave(paste0(aucdir, "High_proportion_heatmap_grouped_bio_dend.pdf"),
       plot = p_combined_heatmap, width = 5, height = 2.3, units = "in", dpi = 300)
