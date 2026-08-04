library(data.table)
library(dplyr)

source(file.path("analysis", "_shared", "cluster_labels.R"))

# ══════════════════════════════════════════════════════════════════════════════
# PFC vs V1 comparison via scDRS's native GROUP-LEVEL enrichment Z-score
# (mclust cluster x Group x region), run twice: once pooling all donors (03_run_mclust_group_region_analysis.sh), 
# once restricted to the 11-donor paired cohort (h5ad filtered before scDRS's group-analysis
# runs (04_run_mclust_group_region_analysis_paired.sh).
#
# scDRS's --group-analysis mclust_Group_region pools ALL cells in a (cluster, Group,
# region) group and gives one Z-score (assoc_mcz) per group.
# ══════════════════════════════════════════════════════════════════════════════

resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
gene_number <- 253
trait <- "ASC_Pfdr"
min_cells_region <- 150
age_order <- c("First_trimester", "Second_trimester", "Infancy", "Adolescence")

en_labels <- get_en_labels()
in_labels <- get_in_labels()

group_dir <- file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis", paste0(gene_number, "_genes"))

# EN-L4-IT and EN-L4-IT-V1 merged into cluster id "27" upstream (01_2_prepare_downstream_metadata.ipynb)
label_overrides <- c("27" = "EN-L4-IT & EN-L4-IT-V1")

rg_labels <- c("vRG", "vRG_oRG", "oRG_tRG", "vRG_oRG_tRG", "oRG and tRG", "vRG and oRG",
               "RG-vRG", "RG-oRG", "RG-tRG")

load_group_z <- function(lineage_name, folder_name, lineage_prefix, label_map) {
  path <- file.path(group_dir, folder_name, paste0(trait, ".scdrs_group.mclust_Group_region"))
  dt <- fread(path, sep = "\t")
  dt <- dt[grepl(paste0("^", lineage_prefix, "lineage_"), group)]
  dt[, cluster   := sub(paste0("^", lineage_prefix, "lineage_(\\d+)_.*$"), "\\1", group)]
  remainder       <- sub(paste0("^", lineage_prefix, "lineage_\\d+_(.*)$"), "\\1", dt$group)
  dt[, region    := sub("^.*_(PFC|V1)$", "\\1", remainder)]
  dt[, age_group := sub("_(PFC|V1)$", "", remainder)]
  dt[, label     := fifelse(cluster %in% names(label_overrides), label_overrides[cluster], label_map[cluster])]
  dt <- dt[!is.na(label) & region %in% c("PFC", "V1")]
  dt[, .(lineage = lineage_name, cluster, label, Group = age_group, region, n_cell, assoc_mcz)]
}

# ── "EN-L4-IT" replaced by the merged "EN-L4-IT & EN-L4-IT-V1" label (see label_overrides
# above); labs_present in compute_family_boundary() are display_labels ("EN-IT-Immature
# 1", "EN-IT-Immature 2", ...) for any label with >1 cluster, so it strips the trailing
# " N" rank suffix before matching against these plain (unexpanded) family lists, while
# still using each entry's actual position in labs_present for the boundary. ──
canonical_mclust_groups <- list(
  EN = c("IPC-EN", "EN-Newborn_IPC-EN", "EN-Newborn", "EN-Newborn_EN-IT-Immature",
         "EN-IT-Immature", "EN-IT-Immature_EN-L2_3-IT", "EN-L2_3-IT", "EN-L2_3_4-IT", "EN-L4-IT & EN-L4-IT-V1",
         "EN-L4_5-IT", "EN-L5-IT", "EN-L6-IT",
         "EN-Non-IT-Immature", "EN-L5_6-NP_EN-L5-ET_EN-Non-IT-Immature",
         "EN-L6-CT_EN-Non-IT-Immature", "EN-L6-CT_EN-Non-IT-Immature_EN-L6b", "EN-L6b"),
  IN = c("IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-Immature_IN-dLGE-Immature", "IN-CGE-VIP",
         "IN-CGE-VIP_IN-CGE-Immature", "IN-CGE-SNCG", "IN-CGE-LAMP5",
         "IN-MGE-Immature", "IN-MGE-Immature_IN-MGE-SST-upper", "IN-MGE-SST-upper", "IN-MGE-SST-deep", "IN-MGE-PV", "Tri-IPC")
)
canonical_mclust_order <- unlist(canonical_mclust_groups, use.names = FALSE)

IT_family_mclust    <- c("EN-IT-Immature", "EN-IT-Immature_EN-L2_3-IT", "EN-L2_3-IT", "EN-L2_3_4-IT",
                         "EN-L4-IT & EN-L4-IT-V1", "EN-L4_5-IT", "EN-L5-IT", "EN-L6-IT")
NonIT_family_mclust <- c("EN-Non-IT-Immature", "EN-L5_6-NP_EN-L5-ET_EN-Non-IT-Immature",
                         "EN-L6-CT_EN-Non-IT-Immature", "EN-L6-CT_EN-Non-IT-Immature_EN-L6b", "EN-L6b")
CGE_family_mclust   <- c("IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-Immature_IN-dLGE-Immature",
                         "IN-CGE-VIP", "IN-CGE-VIP_IN-CGE-Immature", "IN-CGE-SNCG", "IN-CGE-LAMP5")
MGE_family_mclust   <- c("IN-MGE-Immature", "IN-MGE-Immature_IN-MGE-SST-upper", "IN-MGE-SST-upper", "IN-MGE-SST-deep", "IN-MGE-PV")

compute_family_boundary <- function(labs_present, family_before, family_after) {
  base_labels <- sub(" [0-9]+$", "", labs_present)
  idx_before <- which(base_labels %in% family_before)
  idx_after  <- which(base_labels %in% family_after)
  if (length(idx_before) == 0 || length(idx_after) == 0) return(NA_real_)
  (max(idx_before) + min(idx_after)) / 2
}

prenatal_groups  <- c("First_trimester", "Second_trimester")
postnatal_groups <- c("Infancy", "Adolescence")

compute_prepost_boundary <- function(groups_present_chrono) {
  bottom_to_top <- rev(groups_present_chrono)
  y_pos <- setNames(seq_along(bottom_to_top), bottom_to_top)
  last_pre  <- tail(intersect(groups_present_chrono, prenatal_groups), 1)
  first_post <- head(intersect(groups_present_chrono, postnatal_groups), 1)
  if (length(last_pre) == 0 || length(first_post) == 0) return(NA_real_)
  (y_pos[[last_pre]] + y_pos[[first_post]]) / 2
}

library(ggplot2)
library(scales)
library(RColorBrewer)
heatmap_colors <- brewer.pal(11, "RdBu")

plot_z_diff_dotplot <- function(dt, title, row_order = canonical_mclust_order) {
  dt <- dt[!is.na(z_diff_V1_minus_PFC)]
  # Expand each canonical label into however many display_labels it actually has in this
  # data (e.g. "EN-Newborn" -> "EN-Newborn 1".."EN-Newborn 4"), at that label's position
  # in row_order -- a label with just one cluster expands to itself unchanged.
  display_map <- unique(dt[, .(label, display_label)])
  present <- unlist(lapply(row_order, function(lab) sort(display_map[label == lab, display_label])))
  leftover <- setdiff(unique(dt$display_label), present)
  label_order <- c(present, leftover)
  dt[, label := factor(display_label, levels = label_order)]
  dt[, Group := factor(Group, levels = rev(age_order))]

  zmax <- max(abs(dt$z_diff_V1_minus_PFC), na.rm = TRUE)
  color_breaks <- c(-zmax, 0, zmax)
  color_labels <- as.character(round(color_breaks, 2))
  color_labels[1]                    <- paste0(color_labels[1], "\n(PFC higher)")
  color_labels[length(color_labels)] <- paste0(color_labels[length(color_labels)], "\n(V1 higher)")

  groups_present_chrono <- age_order[age_order %in% unique(dt$Group)]
  hline_y <- compute_prepost_boundary(groups_present_chrono)

  en_present <- label_order[label_order %in% unique(dt[lineage == "EN", label])]
  in_present <- label_order[label_order %in% unique(dt[lineage == "IN", label])]
  n_en <- length(en_present)
  en_in_boundary_x <- if (n_en > 0 && length(in_present) > 0) n_en + 0.5 else NA_real_
  en_vline_x <- compute_family_boundary(en_present, IT_family_mclust, NonIT_family_mclust)
  in_vline_x_raw <- compute_family_boundary(in_present, CGE_family_mclust, MGE_family_mclust)
  in_vline_x <- if (!is.na(in_vline_x_raw)) in_vline_x_raw + n_en else NA_real_

  p <- ggplot(dt, aes(x = label, y = Group))
  p <- p + geom_hline(yintercept = hline_y, color = "black", linewidth = 0.5, linetype = "dashed")
  p <- p + geom_vline(xintercept = en_in_boundary_x, color = "black", linewidth = 0.8, linetype = "solid")
  p <- p + geom_vline(xintercept = en_vline_x, color = "black", linewidth = 0.5, linetype = "dashed")
  p <- p + geom_vline(xintercept = in_vline_x, color = "black", linewidth = 0.5, linetype = "dashed")

  p +
    geom_point(aes(size = n_cells_total, fill = z_diff_V1_minus_PFC), shape = 21, color = "black", stroke = 0.5) +
    scale_fill_gradientn(
      colours = rev(heatmap_colors),
      limits = c(-zmax, zmax),
      name = "Z diff\n(V1 - PFC)", breaks = color_breaks,
      labels = color_labels,
      guide = guide_colorbar(barheight = unit(2, "cm"), ticks = TRUE)
    ) +
    scale_size_continuous(range = c(1, 6), name = "n cells\n(PFC + V1)") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, size = 6),
      axis.text.y  = element_text(size = 6),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title   = element_text(size = 8, face = "bold"),
      plot.caption = element_text(size = 8, hjust = 0),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 6)
    )
}

# ── Full pipeline (load -> filter -> widen -> display_label -> save -> plot), reused for
# both the all-donor and paired-donor-cohort variants -- folder_suffix picks which scDRS
# output folder to read ("" = all donors, "_paired" = the 11-donor cohort), file_suffix/
# title_suffix keep the two variants' outputs from colliding. ──
run_zdiff_analysis <- function(folder_suffix, file_suffix, title_suffix) {
  dt_z <- rbindlist(list(
    load_group_z("EN", paste0("EN", folder_suffix), "EN", en_labels),
    load_group_z("IN", paste0("IN", folder_suffix), "IN", in_labels)
  ), use.names = TRUE)

  dt_z <- dt_z[!label %in% rg_labels]
  dt_z <- dt_z[Group != "Third_trimester"]

  keep_keys <- dt_z[, .(n_regions = uniqueN(region), min_n = min(n_cell)), by = .(lineage, cluster, label, Group)
  ][n_regions == 2 & min_n >= min_cells_region, .(lineage, cluster, label, Group)]
  dt_z <- dt_z[keep_keys, on = .(lineage, cluster, label, Group), nomatch = 0]

  wide <- dcast(dt_z, lineage + cluster + label + Group ~ region, value.var = c("assoc_mcz", "n_cell"))
  setnames(wide, c("assoc_mcz_PFC", "assoc_mcz_V1"), c("z_PFC", "z_V1"))
  wide[, z_diff_V1_minus_PFC := z_V1 - z_PFC]
  wide[, n_cells_total := n_cell_PFC + n_cell_V1]

  # display_label: label as-is when it covers one cluster, "label N" when it covers
  # several (N = rank by cluster id, consistent across Group) -- same convention as
  # fig05's prep_ordered() (row_label = paste0(label, " ", rank)).
  cluster_map <- unique(wide[, .(label, cluster)])
  cluster_map[, cluster_num := suppressWarnings(as.integer(cluster))]
  setorder(cluster_map, label, cluster_num)
  cluster_map[, n_cluster := .N, by = label]
  cluster_map[, rank := seq_len(.N), by = label]
  cluster_map[, display_label := fifelse(n_cluster > 1, paste0(label, " ", rank), label)]
  wide <- merge(wide, cluster_map[, .(label, cluster, display_label)], by = c("label", "cluster"))

  wide[, Group := factor(Group, levels = age_order)]
  setorder(wide, lineage, label, cluster, Group)

  fwrite(wide, file.path(resultdir, "WangNature/slingshot",
                         paste0("paired_PFC_V1_mclust_group_Zdiff", file_suffix, "_min", min_cells_region, "cells.csv")))

  p_dot <- plot_z_diff_dotplot(copy(wide), paste0("Group-level scDRS Z-score, V1 vs PFC (mclust", title_suffix, "), n_cells >= ", min_cells_region))
  ggsave(file.path(resultdir, "WangNature/slingshot",
                   paste0("paired_PFC_V1_mclust_group_Zdiff", file_suffix, "_dotplot_min", min_cells_region, "cells.pdf")),
         p_dot, width = 9, height = 3.5, units = "in")

  wide
}

results_all_donors   <- run_zdiff_analysis(folder_suffix = "",        file_suffix = "",                title_suffix = "")
results_paired_donors <- run_zdiff_analysis(folder_suffix = "_paired", file_suffix = "_paired_donors",  title_suffix = ", paired-donor cohort")
