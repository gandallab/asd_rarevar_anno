library(data.table)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(Seurat)

source(file.path("analysis", "_shared", "cluster_labels.R"))

datapath <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
outdir    <- file.path(resultdir, "WangNature/slingshot/mclust_type_composition")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

en_labels <- get_en_labels()
in_labels <- get_in_labels()

# ── Predominant-label rule: a single type >70% of the cluster gets that label outright.
# Otherwise, the combined label lists every type individually comprising >20% of the
# cluster (ordered by proportion) 
TOP1_SINGLE <- 0.70
MIN_SHARE   <- 0.20

# ── Load one lineage's per-cell (mclust, type) table and build the composition +
# predominant-label summary, cross-checked against the existing assigned label.
# `cluster_col` is the RAW numeric mclust id column ("mclust25" for EN, "mclust22" for
# IN) -- NOT the "mclust" column also present in these CSVs, which is already prefixed
# ("ENlineage_13"/"INlineage_9") and won't match cluster_labels.R's plain-numeric keys
# ("1", "2", ...). ──
build_composition <- function(lineage_name, csv_path, cluster_col, label_map) {
  dt <- fread(csv_path)
  dt[, mclust := as.character(get(cluster_col))]
  dt[, type   := as.character(type)]
  
  composition <- dt[, .N, by = .(mclust, type)]
  setnames(composition, "N", "n")
  composition[, prop := n / sum(n), by = mclust]
  setorder(composition, mclust, -prop)
  
  summary_dt <- composition[, {
    ord <- order(-prop)
    types_sorted <- type[ord]
    props_sorted <- prop[ord]
    top1_type <- types_sorted[1]
    top1_prop <- props_sorted[1]
    top2_type <- if (.N >= 2) types_sorted[2] else NA_character_
    top2_prop <- if (.N >= 2) props_sorted[2] else 0
    top3_type <- if (.N >= 3) types_sorted[3] else NA_character_
    top3_prop <- if (.N >= 3) props_sorted[3] else 0
    
    if (top1_prop > TOP1_SINGLE) {
      quant_label <- top1_type
    } else {
      included <- types_sorted[props_sorted > MIN_SHARE]
      quant_label <- if (length(included) >= 1) paste(included, collapse = "_")
      else paste0(top1_type, " (mixed)")
    }
    
    .(top1_type = top1_type, top1_prop = top1_prop,
      top2_type = top2_type, top2_prop = top2_prop,
      top3_type = top3_type, top3_prop = top3_prop,
      quant_label = quant_label)
  }, by = mclust]
  
  # mclust cluster id here is the raw numeric mclust label ("1", "2", ...), not the
  # "ENlineage_"/"INlineage_"-prefixed id used downstream in scDRS's group key -- match
  # cluster_labels.R's assigned label on that raw id for comparison.
  summary_dt[, assigned_label := label_map[mclust]]
  summary_dt[, lineage := lineage_name]
  summary_dt <- summary_dt[order(as.integer(mclust))]
  
  composition[, lineage := lineage_name]
  list(composition = composition, summary = summary_dt)
}

en  <- build_composition("EN", file.path(resultdir, "WangNature/slingshot/EN/EN_lineage_mclust.csv"), "mclust25", en_labels)
in_ <- build_composition("IN", file.path(resultdir, "WangNature/slingshot/IN/IN_lineage_mclust.csv"), "mclust22", in_labels)

composition_all <- rbindlist(list(en$composition, in_$composition), use.names = TRUE)
summary_all     <- rbindlist(list(en$summary, in_$summary), use.names = TRUE)

fwrite(composition_all, file.path(outdir, "mclust_type_composition.csv"))
fwrite(summary_all,     file.path(outdir, "mclust_type_predominant_label.csv"))

# ── One shared, named colour palette across the UNION of EN's and IN's atlas types for the stacked bar and heatmap plots, so the same type is always the same colour
all_types <- sort(unique(composition_all$type))
shared_pal <- setNames(DiscretePalette_scCustomize(num_colors = length(all_types), palette = "polychrome"), all_types)

# Visualize: stacked bar (proportion per mclust cluster, fill = atlas type) and a
# companion heatmap, one of each per lineage. 
plot_composition <- function(composition, summary_dt, lineage_name, pal) {
  # NA assigned_label sorts after every real label (not before), tie-broken by cluster id.
  sort_key <- ifelse(is.na(summary_dt$assigned_label), paste0("zzz_", summary_dt$mclust), summary_dt$assigned_label)
  cluster_order <- summary_dt[order(sort_key, as.integer(mclust)), mclust]
  axis_labels <- setNames(
    paste0(ifelse(is.na(summary_dt$assigned_label), "?", summary_dt$assigned_label), " (", summary_dt$mclust, ")"),
    summary_dt$mclust
  )[cluster_order]
  
  is_flagged <- setNames(summary_dt$top1_prop <= TOP1_SINGLE, summary_dt$mclust)[cluster_order]
  label_colour <- unname(ifelse(is_flagged, "#B2182B", "black"))
  label_face   <- unname(ifelse(is_flagged, "bold", "plain"))
  
  composition <- composition[mclust %in% cluster_order]
  composition[, mclust := factor(mclust, levels = cluster_order)]
  
  p_bar <- ggplot(composition, aes(x = mclust, y = prop, fill = type)) +
    geom_col() +
    scale_x_discrete(labels = axis_labels) +
    scale_fill_manual(values = pal) +
    labs(x = "mclust cluster (assigned label)", y = "Proportion of cells",
         fill = "Atlas type", title = paste0(lineage_name, ": mclust vs atlas type composition")) +
    theme_minimal() +
    theme(axis.title = element_text(size = 8,face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6,
                                     colour = label_colour, face = label_face),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.position = "top")

ggsave(file.path(outdir, paste0(lineage_name, "_mclust_type_stacked_bar.pdf")),
       p_bar, width = max(8, 0.35 * length(cluster_order)), height = 6)

home_cluster <- composition[composition[, .I[which.max(n)], by = type]$V1]
home_cluster[, home_label := axis_labels[as.character(mclust)]]
type_sort_key <- ifelse(is.na(home_cluster$home_label), paste0("zzz_", home_cluster$type), home_cluster$home_label)
type_order <- home_cluster[order(type_sort_key, type), type]
composition[, type := factor(type, levels = type_order)]

p_heat <- ggplot(composition, aes(x = mclust, y = type, fill = prop)) +
  geom_tile() +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1)) +
  labs(x = "mclust cluster (assigned label)", y = "Atlas type", fill = "Proportion",
       title = paste0(lineage_name, ": mclust vs atlas type composition (heatmap)")) +
  theme_minimal() +
  theme(axis.title = element_text(size = 8,face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6,
                                   colour = label_colour, face = label_face),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))

ggsave(file.path(outdir, paste0(lineage_name, "_mclust_type_heatmap.pdf")),
       p_heat, width = max(8, 0.35 * length(cluster_order)), height = 6)

invisible(list(bar = p_bar, heat = p_heat, n_clusters = length(cluster_order), cluster_order = cluster_order))
}

en_plots <- plot_composition(en$composition, en$summary, "EN", shared_pal)
in_plots <- plot_composition(in_$composition, in_$summary, "IN", shared_pal)

supp_fig <- (en_plots$bar / in_plots$bar) +
  plot_annotation(
    tag_levels = "A",
    title = "mclust subcluster composition by atlas cell type",
    caption = paste(
      "Cluster labels in red/bold: no single atlas `type` exceeds 70% of cells in that",
      "mclust subcluster; these were assigned a combined two- (or three-) type label",
      "(see Methods)."
    )
  )

ggsave(file.path(outdir, "SuppFig_mclust_type_composition.pdf"), supp_fig,
       width = 6,
       height = 5.5)

# Save the two bar plots for 10_plot_lineage_marker_dotplots.R to pick up and combine
# with its marker DotPlots into one composition+marker supplementary figure.
saveRDS(list(EN = en_plots$bar, IN = in_plots$bar), file.path(outdir, "composition_bar_plots.rds"))

# Canonical marker-gene (from Wang et al) DotPlot per mclust cluster, EN and IN lineages
marker_ref <- get_marker_reference()
EN_TYPES <- c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature",
              "EN-L2_3-IT", "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature",
              "EN-L5-ET", "EN-L5_6-NP", "EN-L6-CT", "EN-L6b")
IN_TYPES <- c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-Glia", "IN-dLGE-Immature",
              "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-Mix-LAMP5",
              "IN-MGE-Immature", "IN-MGE-SST", "IN-MGE-PV")

EN_MARKER_REF <- marker_ref[EN_TYPES]
IN_MARKER_REF <- marker_ref[IN_TYPES]
EN_MARKERS <- unique(unlist(EN_MARKER_REF, use.names = FALSE))
IN_MARKERS <- unique(unlist(IN_MARKER_REF, use.names = FALSE))

dat_en_marker <- readRDS(file.path(datapath, "dat_EN_lineage_mclust_slingshot.rds"))
DefaultAssay(dat_en_marker) <- "SCT"
dat_en_marker$mclust25[dat_en_marker$mclust25 == 20] <- 22
dat_in_marker <- readRDS(file.path(datapath, "dat_IN_lineage_mclust_slingshot.rds"))
DefaultAssay(dat_in_marker) <- "SCT"
dat_in_marker$mclust22[dat_in_marker$mclust22 == "21"] <- "15"

plot_marker_dotplot <- function(obj, cluster_col, label_map, features, lineage_name, marker_ref_subset, cluster_order) {
  cluster_id <- as.character(obj@meta.data[[cluster_col]])
  assigned_label <- label_map[cluster_id]
  obj@meta.data$mclust_label <- paste0(ifelse(is.na(assigned_label), "?", assigned_label), " (", cluster_id, ")")
  present <- unique(cluster_id)
  cluster_order <- intersect(cluster_order, present)
  level_order <- paste0(ifelse(is.na(label_map[cluster_order]), "?", label_map[cluster_order]), " (", cluster_order, ")")
  obj@meta.data$mclust_label <- factor(obj@meta.data$mclust_label, levels = level_order)
  Idents(obj) <- "mclust_label"
  
  features <- features[features %in% rownames(obj)]
  
  gene_type_label <- get_marker_gene_types(features, marker_ref_subset)
  wrap_type_label <- function(label, per_line = 2) {
    parts <- strsplit(label, " / ", fixed = TRUE)[[1]]
    if (length(parts) <= per_line) return(paste(parts, collapse = " / "))
    chunks <- split(parts, ceiling(seq_along(parts) / per_line))
    paste(vapply(chunks, paste, character(1), collapse = " / "), collapse = "\n")
  }
  gene_type_label_wrapped <- vapply(gene_type_label, wrap_type_label, character(1))
  gene_axis_labels <- setNames(paste0(features, "\n(", gene_type_label_wrapped, ")"), features)
  
  p <- DotPlot_scCustom(obj, features = features, dot.scale = 4) +
    scale_x_discrete(labels = gene_axis_labels) +
    labs(title = paste0(lineage_name, ": marker gene expression by mclust cluster"),
         x = "Marker gene (reference cell type)", y = "mclust cluster (assigned label)") +
    theme(axis.title = element_text(size = 8,face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.position = "bottom")
  
  ggsave(file.path(outdir, paste0(lineage_name, "_mclust_marker_dotplot.pdf")),
         p, width = 12, height = 8, unit = "in")
  
  invisible(list(plot = p, n_clusters = length(cluster_order), cluster_order = cluster_order))
}

en_dot <- plot_marker_dotplot(dat_en_marker, "mclust25", en_labels, EN_MARKERS, "EN", EN_MARKER_REF, en_plots$cluster_order)
in_dot <- plot_marker_dotplot(dat_in_marker, "mclust22", in_labels, IN_MARKERS, "IN", IN_MARKER_REF, in_plots$cluster_order)

composition = readRDS(file.path(outdir, "composition_bar_plots.rds"))
no_x_axis <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

predominant_label <- read.csv(file.path(outdir, "mclust_type_predominant_label.csv"))
cluster_axis_theme <- function(lineage_name, cluster_order) {
  sub <- predominant_label[predominant_label$lineage == lineage_name, ]
  flagged <- setNames(sub$top1_prop <= 0.7, as.character(sub$mclust))
  is_flagged <- unname(ifelse(is.na(flagged[cluster_order]), FALSE, flagged[cluster_order]))
  theme(axis.text.x = element_text(colour = ifelse(is_flagged, "#B2182B", "black"),
                                   face = ifelse(is_flagged, "bold", "plain")))
}

combined_fig <- (composition$EN + no_x_axis + composition$IN + no_x_axis) /
  ((en_dot$plot + scale_x_discrete() + coord_flip() + cluster_axis_theme("EN", en_dot$cluster_order)) +
     (in_dot$plot + scale_x_discrete() + coord_flip() + cluster_axis_theme("IN", in_dot$cluster_order))) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(outdir, "SuppFig_mclust_and_markers.pdf"), combined_fig,
       width = 7.5, height = 10, units = "in", dpi = 300)