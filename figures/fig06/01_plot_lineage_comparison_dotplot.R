library(data.table)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
library(patchwork)
library(ggplotify)
library(mgcv)
library(cowplot)
library(scales)

source(file.path("analysis", "_shared", "cluster_labels.R"))

# ── Paths ──────────────────────────────────────────────────────────────────────
datapath  <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"

trait_ddid       <- "ASC_Pfdr_ddid_new_DMN" # ASD_DH
trait_noddid     <- "ASC_Pfdr_noddid_new_DMN" # ASD_DL+ASD_DA
trait_noddid_new <- "ASC_Pfdr_noddid_perproband_new_DMN" # ASD_DA

en_labels <- get_en_labels()
in_labels <- get_in_labels()

# ── Load scDRS group-level enrichment results (cluster level) ──────────────────
group_dir <- file.path(
  resultdir, "WangNature/scDRS/run_ddid_noddid/downstream_analysis")

lineage_scDRS_ddid <- fread(file.path(
  group_dir, "ddid_genes",
  paste0(trait_ddid, ".scdrs_group.mclust.csv")
))
lineage_scDRS_noddid <- fread(file.path(
  group_dir, "noddid_genes",
  paste0(trait_noddid, ".scdrs_group.mclust.csv")
))
lineage_scDRS_noddid_new <- fread(file.path(
  group_dir, "noddid_genes",
  paste0(trait_noddid_new, ".scdrs_group.mclust.csv")
))

# ── Helper: compute cluster-level enrichment comparison ───────────────────────
compute_cluster_comparison <- function(lineage_scDRS_ddid, lineage_scDRS_noddid,
                                       cell_type_prefix) {
  merge(
    lineage_scDRS_ddid %>%
      filter(grepl(cell_type_prefix, group)) %>%
      dplyr::select(cluster, assoc_mcz, assoc_mcp, assoc_mcp_fdr) %>%
      dplyr::rename(z_ddid = assoc_mcz, p_ddid = assoc_mcp, fdr_ddid = assoc_mcp_fdr),
    lineage_scDRS_noddid %>%
      filter(grepl(cell_type_prefix, group)) %>%
      dplyr::select(cluster, assoc_mcz, assoc_mcp, assoc_mcp_fdr) %>%
      dplyr::rename(z_noddid = assoc_mcz, p_noddid = assoc_mcp, fdr_noddid = assoc_mcp_fdr),
    by = "cluster"
  ) %>%
    mutate(cluster = as.character(cluster))
}

# ── Helper: prepare combined EN+IN result_data with joint FDR correction ───────
prepare_result_data <- function(cluster_comparison_EN, cluster_comparison_IN,
                                node_order_EN, node_order_IN,
                                sig_clusters_EN, sig_clusters_IN,
                                resultdir, file_suffix = "simple") {
  
  order_EN <- node_order_EN[node_order_EN %in% sig_clusters_EN]
  order_IN <- node_order_IN[node_order_IN %in% sig_clusters_IN]
  
  prep_ordered <- function(comparison, order, labels) {
    data.frame(cluster = order) %>%
      left_join(comparison %>% mutate(cluster = as.character(cluster)),
                by = "cluster") %>%
      mutate(label = labels[cluster]) %>%
      group_by(label) %>%
      mutate(
        n_cluster = n(),
        rank      = row_number(),
        row_label = ifelse(n_cluster > 1, paste0(label, " ", rank), label)
      ) %>%
      ungroup() %>%
      mutate(row_label = factor(row_label, levels = row_label),
             x_pos     = as.numeric(row_label))
  }
  
  df_EN <- prep_ordered(cluster_comparison_EN, order_EN, en_labels)
  df_IN <- prep_ordered(cluster_comparison_IN, order_IN, in_labels)
  
  # Joint FDR: EN + IN together, separately for ddid and noddid
  n_EN <- nrow(df_EN)
  n_IN <- nrow(df_IN)
  
  fdr_ddid_joint   <- p.adjust(c(df_EN$p_ddid,   df_IN$p_ddid),   method = "fdr")
  fdr_noddid_joint <- p.adjust(c(df_EN$p_noddid, df_IN$p_noddid), method = "fdr")
  
  df_EN <- df_EN %>%
    mutate(fdr_ddid   = fdr_ddid_joint[seq_len(n_EN)],
           fdr_noddid = fdr_noddid_joint[seq_len(n_EN)])
  df_IN <- df_IN %>%
    mutate(fdr_ddid   = fdr_ddid_joint[n_EN + seq_len(n_IN)],
           fdr_noddid = fdr_noddid_joint[n_EN + seq_len(n_IN)])
  
  add_enrichment <- function(df) {
    df %>%
      mutate(enrichment_type = case_when(
        fdr_ddid  < 0.05 & fdr_noddid >= 0.05 ~ "ddid only",
        fdr_noddid < 0.05 & fdr_ddid  >= 0.05 ~ "noddid only",
        fdr_ddid  < 0.05 & fdr_noddid < 0.05  ~ "both",
        TRUE                                    ~ "neither"
      )) %>%
      dplyr::select(cluster, label, row_label, x_pos,
                    z_ddid, z_noddid, p_ddid, p_noddid,
                    fdr_ddid, fdr_noddid, enrichment_type)
  }
  
  result_EN <- add_enrichment(df_EN)
  result_IN <- add_enrichment(df_IN)
  
  for (ct in c("EN", "IN")) {
    rd  <- if (ct == "EN") result_EN else result_IN
    csv <- file.path(resultdir, "WangNature/slingshot", ct,
                     paste0("scDRS_ddid_vs_noddid_cluster_heatmap_",
                            trait_ddid, "_", ct, "_", file_suffix, "_result_data.csv"))
    write.table(rd, file = csv, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Saved result data: ", csv)
  }
  
  invisible(list(EN = result_EN, IN = result_IN))
}

# ── Helper: combined EN+IN cluster-level heatmap ─────────────────
# rows    = "both": ddid + noddid rows; "noddid": noddid only
# show_x  = FALSE: suppress x-axis text and ticks (use when panel below shares the same columns)
make_combined_heatmap <- function(res_EN_data, res_IN_data,
                                  theme_fig, z_limits = NULL,
                                  rows = c("both", "noddid"),
                                  show_x = TRUE,
                                  show_legend = TRUE) {
  rows <- match.arg(rows)
  
  df_en <- res_EN_data %>%
    mutate(cell_group = "EN",
           x_pos      = as.numeric(factor(row_label, levels = row_label)))
  n_en <- nrow(df_en)
  
  df_in <- res_IN_data %>%
    mutate(cell_group = "IN",
           x_pos      = as.numeric(factor(row_label, levels = row_label)) + n_en)
  
  make_long_rows <- function(df, score_rows) {
    rows_list <- lapply(score_rows, function(sr) {
      z_col   <- if (sr == "ddid") "z_ddid"   else "z_noddid"
      fdr_col <- if (sr == "ddid") "fdr_ddid"  else "fdr_noddid"
      p_col   <- if (sr == "ddid") "p_ddid"    else "p_noddid"
      df %>% mutate(
        score_type = sr,
        z_val      = .data[[z_col]],
        sig_label  = case_when(
          .data[[fdr_col]] < 0.001 ~ "***",
          .data[[fdr_col]] < 0.01  ~ "**",
          .data[[fdr_col]] < 0.05  ~ "*",
          .data[[p_col]]   < 0.05  ~ "+",
          TRUE                     ~ ""
        )
      )
    })
    bind_rows(rows_list) %>%
      mutate(score_type = factor(score_type,
                                 levels = if (length(score_rows) > 1)
                                   c("noddid", "ddid")
                                 else score_rows))
  }
  
  score_rows  <- if (rows == "both") c("noddid", "ddid") else "noddid"
  df_combined <- bind_rows(make_long_rows(df_en, score_rows),
                           make_long_rows(df_in, score_rows))
  
  x_ref    <- bind_rows(df_en, df_in) %>% arrange(x_pos)
  x_breaks <- x_ref$x_pos
  x_labels <- as.character(x_ref$row_label)
  
  p <- ggplot(df_combined, aes(x = x_pos, y = score_type)) +
    geom_tile(aes(fill = z_val), color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 3, color = "white") +
    scale_fill_gradientn(colours = c("#F7F7F7", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"),
                         name = "Z score",
                         limits = z_limits) + 
    scale_x_continuous(breaks = x_breaks, labels = x_labels,
                       expand = c(0.01, 0.01)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() + theme_fig +
    theme(
      axis.text.y     = element_text(size = 6),
      panel.grid      = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.title    = element_text(size = 6),
      legend.text     = element_text(size = 6)
    )
  
  if (show_x) {
    p <- p + theme(
      axis.text.x  = element_text(angle = 0, hjust = 1, size = 6),
      axis.ticks.x = element_line()
    )
  } else {
    p <- p + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# ── Cluster-level node orders and sig_clusters ────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
node_order_EN <- c(
  "9", "8", "22", "3", "1", "12", "16",
  "4", "21", "13", "15", "23", "19", "24",
  "5", "17", "6", "10", "18", "7", "2", "14", "11"
)

node_order_IN <- c(
  "4", "7", "8", "13", "19", "17", "20", "14",
  "12", "6", "3", "5", "10", "2", "22", "18",
  "1", "15", "16", "9", "11"
)

whole_set <- fread(
  "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/run_v8/downstream_analysis/253_genes/ASC_Pfdr.scdrs_group.mclust.csv",
  sep = "\t"
)

sig_clusters_EN <- whole_set %>%
  filter(grepl("^EN|^IPC", label), assoc_mcp_fdr < 0.05) %>%
  pull(cluster) %>% as.character()

sig_clusters_IN <- whole_set %>%
  filter(grepl("^IN", label), assoc_mcp_fdr < 0.05) %>%
  pull(cluster) %>% as.character()

# ── Compute cluster comparisons ────────────────────────────────────────────────
cluster_comparison_EN <- compute_cluster_comparison(
  lineage_scDRS_ddid, lineage_scDRS_noddid, "^EN")
cluster_comparison_IN <- compute_cluster_comparison(
  lineage_scDRS_ddid, lineage_scDRS_noddid, "^IN")

# noddid_new (per-proband) comparison
cluster_comparison_EN_new <- compute_cluster_comparison(
  lineage_scDRS_ddid, lineage_scDRS_noddid_new, "^EN")
cluster_comparison_IN_new <- compute_cluster_comparison(
  lineage_scDRS_ddid, lineage_scDRS_noddid_new, "^IN")

# ── Prepare cluster result data (joint FDR: EN+IN together) ───────────────────
res <- prepare_result_data(
  cluster_comparison_EN = cluster_comparison_EN,
  cluster_comparison_IN = cluster_comparison_IN,
  node_order_EN = node_order_EN, node_order_IN = node_order_IN,
  sig_clusters_EN = sig_clusters_EN, sig_clusters_IN = sig_clusters_IN,
  resultdir = resultdir, file_suffix = "simple"
)
res_EN <- res$EN
res_IN <- res$IN

res_new <- prepare_result_data(
  cluster_comparison_EN = cluster_comparison_EN_new,
  cluster_comparison_IN = cluster_comparison_IN_new,
  node_order_EN = node_order_EN, node_order_IN = node_order_IN,
  sig_clusters_EN = sig_clusters_EN, sig_clusters_IN = sig_clusters_IN,
  resultdir = resultdir, file_suffix = "perproband_simple"
)
res_EN_new <- res_new$EN
res_IN_new <- res_new$IN

# ══════════════════════════════════════════════════════════════════════════════
# ── Build all panels ──────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
theme_fig <- theme_classic() +
  theme(
    legend.position = "right",
    text            = element_text(size = 6),
    axis.text       = element_text(size = 6),
    axis.title      = element_text(size = 6),
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 6),
    plot.margin     = margin(5, 5, 5, 5)
  )

# ── Global z limits for cluster-level heatmaps ──────────────────
all_z_cluster <- c(
  res_EN$z_ddid,     res_EN$z_noddid,
  res_IN$z_ddid,     res_IN$z_noddid,
  res_EN_new$z_noddid, res_IN_new$z_noddid
)
global_z_limits <- c(min(all_z_cluster, na.rm = TRUE),
                     max(all_z_cluster, na.rm = TRUE))
message("Cluster-level z limits: ", global_z_limits[1], " - ", global_z_limits[2])

# ── cluster-level ddid vs noddid (EN + IN, both rows) ───────────────
# x-axis suppressed: cluster labels are shown in the age-group panels below
panel_a_ht <- make_combined_heatmap(
  res_EN_data = res_EN, res_IN_data = res_IN,
  theme_fig = theme_fig, z_limits = global_z_limits, rows = "both",
  show_x = FALSE, show_legend = FALSE
)

# ── cluster-level noddid_new (per-proband, noddid row only) ──────────
panel_b_ht <- make_combined_heatmap(
  res_EN_data = res_EN_new, res_IN_data = res_IN_new,
  theme_fig = theme_fig, z_limits = global_z_limits, rows = "noddid"
)

# ══════════════════════════════════════════════════════════════════════════════
# ── Build dotplot: cluster × age_group, DDID + NoDDID + NoDDID-pp ──
# ══════════════════════════════════════════════════════════════════════════════
library(scales)

trait_noddid_new <- "ASC_Pfdr_noddid_perproband_new_DMN"

age_group_order <- c(
  "First_trimester", "Second_trimester", "Third_trimester",
  "Infancy", "Adolescence"
)
age_group_labels_dot <- c(
  "First_trimester"  = "1st tri.",
  "Second_trimester" = "2nd tri.",
  "Third_trimester"  = "3rd tri.",
  "Infancy"          = "Infancy",
  "Adolescence"      = "Adolescence"
)

load_group_data_dot <- function(path, lineage_prefix, sig_clusters,
                                label_map, node_order, min_cells = 150) {
  fread(path) %>%
    mutate(
      cluster   = sub(paste0("^", lineage_prefix, "_(\\d+)_.*$"), "\\1", group),
      age_group = sub(paste0("^", lineage_prefix, "_\\d+_(.*)$"),  "\\1", group)
    ) %>%
    filter(cluster %in% sig_clusters,
           age_group %in% age_group_order,
           n_cell >= min_cells) %>%
    mutate(
      label     = label_map[cluster],
      age_group = factor(age_group, levels = age_group_order),
      cluster   = factor(cluster,   levels = node_order)
    ) %>%
    arrange(cluster) %>%
    group_by(label) %>%
    mutate(
      n_cluster  = n_distinct(cluster),
      clust_rank = as.integer(factor(cluster, levels = levels(cluster)))
    ) %>%
    ungroup() %>%
    mutate(line_label = ifelse(n_cluster > 1,
                               paste0(label, " ", clust_rank), label))
}

dat_dot_EN_ddid <- load_group_data_dot(
  file.path(group_dir, "ddid_genes/EN",
            paste0(trait_ddid, ".scdrs_group.mclust_Group")),
  "ENlineage", sig_clusters_EN, en_labels, node_order_EN
) %>% mutate(cell_group = "EN", trait = "DDID")

dat_dot_IN_ddid <- load_group_data_dot(
  file.path(group_dir, "ddid_genes/IN",
            paste0(trait_ddid, ".scdrs_group.mclust_Group")),
  "INlineage", sig_clusters_IN, in_labels, node_order_IN
) %>% mutate(cell_group = "IN", trait = "DDID")

dat_dot_EN_noddid <- load_group_data_dot(
  file.path(group_dir, "noddid_genes/EN",
            paste0(trait_noddid, ".scdrs_group.mclust_Group")),
  "ENlineage", sig_clusters_EN, en_labels, node_order_EN
) %>% mutate(cell_group = "EN", trait = "NoDDID")

dat_dot_IN_noddid <- load_group_data_dot(
  file.path(group_dir, "noddid_genes/IN",
            paste0(trait_noddid, ".scdrs_group.mclust_Group")),
  "INlineage", sig_clusters_IN, in_labels, node_order_IN
) %>% mutate(cell_group = "IN", trait = "NoDDID")

dat_dot_EN_noddid_new <- load_group_data_dot(
  file.path(group_dir, "noddid_genes/EN",
            paste0(trait_noddid_new, ".scdrs_group.mclust_Group")),
  "ENlineage", sig_clusters_EN, en_labels, node_order_EN
) %>% mutate(cell_group = "EN", trait = "NoDDID (per-proband)")

dat_dot_IN_noddid_new <- load_group_data_dot(
  file.path(group_dir, "noddid_genes/IN",
            paste0(trait_noddid_new, ".scdrs_group.mclust_Group")),
  "INlineage", sig_clusters_IN, in_labels, node_order_IN
) %>% mutate(cell_group = "IN", trait = "NoDDID (per-proband)")

df_dot <- bind_rows(
  bind_rows(dat_dot_EN_ddid,       dat_dot_IN_ddid)       %>% mutate(fdr_joint = p.adjust(assoc_mcp, method = "fdr")),
  bind_rows(dat_dot_EN_noddid,     dat_dot_IN_noddid)     %>% mutate(fdr_joint = p.adjust(assoc_mcp, method = "fdr")),
  bind_rows(dat_dot_EN_noddid_new, dat_dot_IN_noddid_new) %>% mutate(fdr_joint = p.adjust(assoc_mcp, method = "fdr"))
) %>%
  mutate(
    sig_label = case_when(
      fdr_joint < 0.001 ~ "***",
      fdr_joint < 0.01  ~ "**",
      fdr_joint < 0.05  ~ "*",
      assoc_mcp < 0.05  ~ "+",
      TRUE              ~ ""
    ),
    alpha_val = case_when(
      fdr_joint < 0.05  ~ 1,
      assoc_mcp < 0.05  ~ 0.5,
      TRUE              ~ 0.25
    ),
    trait     = factor(trait, levels = c("DDID", "NoDDID", "NoDDID (per-proband)"))
  )

# y_order for dotplot x axis
make_y_order_dot <- function(dat) {
  unique(dat %>% arrange(cluster) %>% distinct(cluster, line_label) %>% pull(line_label))
}
y_order_EN_dot <- make_y_order_dot(bind_rows(dat_dot_EN_ddid, dat_dot_EN_noddid))
y_order_IN_dot <- make_y_order_dot(bind_rows(dat_dot_IN_ddid, dat_dot_IN_noddid))
y_order_full_dot <- c(y_order_EN_dot, y_order_IN_dot)

n_EN_dot <- length(y_order_EN_dot)
n_IN_dot <- length(y_order_IN_dot)

# separator positions (x axis = y_order_full_dot index)
vline_EN_IN_dot     <- n_EN_dot + 0.5
last_newborn_idx    <- max(which(grepl("Newborn", y_order_EN_dot)))
vline_newborn_IT    <- last_newborn_idx + 0.5
last_mge_imm_idx    <- n_EN_dot + max(which(grepl("MGE-Immature", y_order_IN_dot)))
vline_IN_mge_cge    <- last_mge_imm_idx + 0.5
hline_pre_post      <- 2.5   # between Infancy (pos 2) and 3rd tri (pos 3) in rev order

panel_c_dot <- ggplot(df_dot, aes(x = line_label, y = age_group)) +
  annotate("rect",
           ymin = 2.5, ymax = 5.5, xmin = -Inf, xmax = Inf,
           fill = "#F5F5F5", alpha = 0.6) +
  geom_point(aes(size = n_cell, color = assoc_mcz, alpha = alpha_val), na.rm = TRUE) +
  scale_alpha_identity() +
  geom_text(aes(label = sig_label),
            size = 2.2, color = "white", fontface = "bold", na.rm = TRUE) +
  geom_hline(yintercept = hline_pre_post,
             color = "grey40", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = vline_EN_IN_dot,
             color = "grey40", linewidth = 0.8, linetype = "solid") +
  geom_vline(xintercept = vline_newborn_IT,
             color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = vline_IN_mge_cge,
             color = "grey50", linewidth = 0.4, linetype = "dashed") +
  scale_color_gradientn(
    colours = c("#2166AC", "#F7F7F7", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"),
    values  = scales::rescale(c(-2.5, 0, 1, 2, 3, 5)),
    limits  = c(-2.5, max(df_dot$assoc_mcz, na.rm = TRUE)),
    oob     = scales::squish,
    name    = "Z score"
  ) +
  scale_size_continuous(
    range  = c(2, 7.5),
    limits = c(0, 6000),
    breaks = c(100, 500, 2000, 6000),
    labels = c("100", "500", "2000", "6000"),
    name   = "n cells"
  ) +
  scale_x_discrete(limits = y_order_full_dot, expand = c(0.04, 0.04)) +
  scale_y_discrete(limits = rev(age_group_order),
                   labels = age_group_labels_dot,
                   expand = c(0.08, 0.08)) +
  annotate("text", x = n_EN_dot / 2,               y = 5.7,
           label = "EN", size = 2.5, fontface = "bold", color = "grey30") +
  annotate("text", x = n_EN_dot + n_IN_dot / 2,    y = 5.7,
           label = "IN", size = 2.5, fontface = "bold", color = "grey30") +
  facet_wrap(~ trait, ncol = 1) +
  labs(x = NULL, y = NULL) +
  theme_classic() + theme_fig +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, size = 5.5),
    axis.text.y      = element_text(size = 6),
    strip.text       = element_text(size = 6, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 6),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.8, "lines"),
    panel.spacing    = unit(0.5, "lines"),
    plot.margin      = margin(15, 5, 5, 5)
  )

# ── Export dotplot data ────────────────────────────────────────────────────────
write.table(
  df_dot %>% dplyr::select(trait, cell_group, cluster, line_label, age_group,
                           n_cell, n_ctrl, assoc_mcp, assoc_mcz,
                           fdr_joint, sig_label, alpha_val),
  file  = file.path(resultdir, "WangNature/slingshot",
                    paste0("scDRS_agegroup_dotplot_", trait_ddid, "_plot_data.tsv")),
  sep   = "\t", row.names = FALSE, quote = FALSE
)
message("Saved dotplot data.")

# ══════════════════════════════════════════════════════════════════════════════
# ── mutation rate scatter plot ───────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
count_perproband <- fread(
  "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_perproband_group.txt",
  sep = "\t"
)
count_DMN <- fread(
  "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_group.txt",
  sep = "\t"
)

p_hat_threshold <- sum(count_DMN$count_new_DMN_ddid) /
  (sum(count_DMN$count_new_DMN_ddid) + sum(count_DMN$count_new_DMN_noddid))

count_combined <- count_perproband %>%
  mutate(d_rate = rate_ddid - rate_noddid) %>%
  left_join(count_DMN %>% dplyr::select(Gene, Group_phat = Group), by = "Gene") %>%
  mutate(
    Group_combined = case_when(
      Group_perproband == "DDID"   & Group_phat == "DDID"   ~ "ASD-DH",
      Group_perproband == "NoDDID" & Group_phat == "NoDDID" ~ "ASD-DA",
      Group_perproband == "DDID"   & Group_phat == "NoDDID" ~ "ASD-DL"
    )
  )

panel_a_scatter <- ggplot(count_combined,
                          aes(x = d_rate, y = p_hat, color = Group_combined)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = p_hat_threshold, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c(
      "ASD-DH"                 = "#F79564",
      "ASD-DA"                   = "#68C3a5",
      "ASD-DL"                = "grey50"
    )
  ) +
  labs(
    x     = "Mutation rate difference (DDID - NoDDID)",
    y     = expression(hat(p)[DDID]),
    color = "Group"
  ) +
  theme_classic() + theme_fig

# ══════════════════════════════════════════════════════════════════════════════
# ── Assemble final figure ─────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
total_fig_width <- 6

# scatter (left) + placeholder (right)
row_a <- plot_grid(
  as.grob(panel_a_scatter), NULL,
  ncol       = 2,
  rel_widths = c(1, 1),
  labels     = c("a", ""), label_size = 8, label_fontface = "bold"
)

# cluster-level heatmap (ddid + noddid rows, no x-axis)
row_b <- plot_grid(
  as.grob(panel_a_ht),
  labels = "b", label_size = 8, label_fontface = "bold"
)

# cluster-level heatmap (noddid per-proband, with x-axis)
row_c <- plot_grid(
  as.grob(panel_b_ht),
  labels = "c", label_size = 8, label_fontface = "bold"
)

# dotplot (3 traits × 5 age groups)
row_d <- plot_grid(
  as.grob(panel_c_dot),
  labels = "d", label_size = 8, label_fontface = "bold"
)

# Relative heights:
#   row a: scatter plot
#   row b: 2 tile rows, no x-axis labels
#   row c: 1 tile row + rotated x-axis labels
#   row d: 3 traits × 5 age_group rows
h_a <- 5
h_b <- 2
h_c <- 1.5
h_d <- 15

final_plot <- plot_grid(
  row_a, row_b, row_c, row_d,
  ncol        = 1,
  rel_heights = c(h_a, h_b, h_c, h_d)
)

total_fig_height <- 7

ggsave(
  filename = file.path(resultdir, "WangNature/slingshot",
                       paste0("scDRS_ddid_vs_noddid_", trait_ddid,
                              "_heatmap_dotplot.pdf")),
  plot   = final_plot,
  width  = total_fig_width,
  height = total_fig_height,
  units  = "in"
)
message("Figure saved.")
